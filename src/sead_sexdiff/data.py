"""Data utilities for downloading resources from the SEA-AD S3 bucket.

This module provides a convenience function `download_data()` to fetch public
files from the SEA-AD single-cell profiling S3 bucket without requiring AWS
credentials. It supports filtering by suffix (e.g., ``.csv``, ``.html``), and
will always attempt to fetch the ``Changelog.html`` file within the specified
prefix.

Bucket browser: https://sea-ad-single-cell-profiling.s3.amazonaws.com/index.html

Example
-------
>>> from sead_sexdiff.data import download_data
>>> downloaded = download_data(prefix="PFC/RNAseq/", suffixes=[".csv", ".html"])  # doctest: +SKIP
>>> print(downloaded)  # doctest: +SKIP

Notes
-----
- The bucket is public and does not require credentials. We configure boto3 to
  use unsigned requests via botocore.
- Files are downloaded preserving their S3 key structure under the local
  destination directory (default: ``data/``).
"""

from __future__ import annotations

import gc
import os
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, List
import json

import boto3
from botocore import UNSIGNED
from botocore.client import Config
import pandas as pd
import scanpy as sc
import numpy as np
import random
import decoupler as dc
import anndata as ad
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


DEFAULT_BUCKET = "sea-ad-single-cell-profiling"
DATA_ROOT_ENV = "SEAD_DATA_ROOT"
DOWNLOAD_MANIFEST_NAME = "download_manifest.json"

# MTG per-donor object prefixes, one per omic
MTG_PREFIXES: dict[str, str] = {
    "rna": "MTG/RNAseq/donors_objects/",
    "atac": "MTG/ATACseq/donor_objects/",
}


def resolve_data_root(explicit: str | Path | None = None) -> Path:
    """Resolve the local data root.

    Parameters
    ----------
    explicit:
        Root given on the command line or by the caller. Takes precedence.

    Returns
    -------
    pathlib.Path
        ``explicit`` if given, else ``$SEAD_DATA_ROOT`` if set and non-empty,
        else ``data`` relative to the current working directory.
    """
    if explicit is not None and str(explicit) != "":
        return Path(explicit)
    env = os.environ.get(DATA_ROOT_ENV)
    if env:
        return Path(env)
    return Path("data")


class DownloadResult(list):
    """List of downloaded paths, with the rest of the download report attached.

    Attributes
    ----------
    skipped:
        Local paths already present at the remote size, not transferred.
    failed:
        S3 keys that were missing or failed to download.
    bytes_downloaded:
        Total bytes transferred.
    manifest_path:
        Path of the updated ``download_manifest.json``, or ``None``.
    """

    skipped: list[Path]
    failed: list[str]
    bytes_downloaded: int
    manifest_path: Path | None


def _ensure_parent_dir(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def _matches_suffixes(key: str, suffixes: Iterable[str] | None) -> bool:
    if not suffixes:
        return True
    for suf in suffixes:
        # Normalize: allow users to pass extensions with or without leading dot
        if not suf:
            # empty string matches everything (but ignore for safety)
            continue
        if suf.startswith("*"):
            # Simple wildcard pattern like *.csv — fallback to endswith after '*'
            suf = suf.lstrip("*")
        if not suf.startswith(".") and "/" not in suf:
            # If it's likely an extension without dot, add dot for convenience
            candidate = "." + suf
        else:
            candidate = suf
        if key.endswith(candidate) or key == suf:
            return True
    return False


def download_data(
    prefix: str = "PFC/RNAseq/",
    *,
    bucket: str = DEFAULT_BUCKET,
    suffixes: Iterable[str] | None = (".csv", ".html"),
    files: Iterable[str] | None = None,
    dest_root: str | Path | None = None,
    show_progress: bool = True,
) -> DownloadResult:
    """Download public SEA-AD data from S3 via boto3.

    Parameters
    ----------
    prefix:
        S3 prefix (folder-like path) to search within the bucket. For example
        "PFC/RNAseq/".
    bucket:
        S3 bucket name. Defaults to "sea-ad-single-cell-profiling".
    suffixes:
        File suffixes/extensions to download (e.g., [".csv", ".html"]). If
        ``None`` or empty, all objects under the prefix will be downloaded.
        You can also pass patterns like "*.csv" or bare extensions like
        "csv". Matching is performed by simple endswith checks.
    files:
        Specific file paths under the given ``prefix`` to download. These may
        be provided either as keys relative to ``prefix`` (e.g.,
        "public_datasets/SEA-AD-cell-annotation.2024-08-27.csv") or as full
        S3 keys starting with the same ``prefix``. When this option is
        provided and non-empty, it overrides ``suffixes`` filtering and only
        the specified files (plus the ``Changelog.html``, if present) are
        attempted.
    dest_root:
        Local root directory where files will be saved. The S3 key structure is
        preserved under this directory. ``None`` resolves through
        :func:`resolve_data_root` (``$SEAD_DATA_ROOT``, else ``data``).

    Returns
    -------
    DownloadResult
        A list of local file paths that were downloaded, carrying the skipped
        paths, failed keys, bytes transferred and manifest path as attributes.

    Notes
    -----
    - To download all files under a folder/prefix (e.g. the donor_objects
      folder), pass ``suffixes=None``. For example:
      ``download_data(prefix="PFC/RNAseq/donor_objects/", suffixes=None)``.
    - This function uses unsigned S3 requests, which is appropriate for public
      buckets. If the bucket ever becomes private, authentication would be
      required and this approach would need updating.
    - Files already present locally at the remote size are skipped. Every
      downloaded or present file is recorded (key, size, ETag, timestamp) in
      ``<dest_root>/<prefix>/download_manifest.json``.
    """

    # Normalize inputs
    if prefix and not prefix.endswith("/"):
        # Keep behavior consistent with folder-like prefixes
        prefix = prefix + "/"

    dest_root = resolve_data_root(dest_root)
    dest_root.mkdir(parents=True, exist_ok=True)

    # Create an anonymous S3 client (no credentials needed for public buckets)
    s3 = boto3.client("s3", config=Config(signature_version=UNSIGNED))

    downloaded = DownloadResult()

    # Internal progress tracker
    class _Progress:
        def __init__(self, total: int | None):
            self.total = total if (isinstance(total, int) and total > 0) else None
            self.seen = 0
            self.start = time.time()
            self._last_print = 0.0

        def cb(self, bytes_amount: int):
            # Callback invoked by boto3 on each chunk
            self.seen += int(bytes_amount)
            if not show_progress:
                return
            now = time.time()
            # Throttle updates to ~10 Hz
            if now - self._last_print < 0.1:
                return
            self._last_print = now
            self._print_line(end="\r")

        def _eta_str(self, rate: float) -> str:
            if rate <= 0:
                return "--:--"
            remaining = (self.total - self.seen) / rate if self.total else 0
            m, s = divmod(int(max(0, remaining)), 60)
            if m >= 60:
                h, m = divmod(m, 60)
                return f"{h:d}h{m:02d}m"
            return f"{m:02d}:{s:02d}"

        def _human(self, n: int) -> str:
            return _human_bytes(n)

        def _print_line(self, end: str = "\n"):
            elapsed = max(1e-6, time.time() - self.start)
            rate = self.seen / elapsed
            rate_str = self._human(int(rate)) + "/s"
            if self.total:
                pct = (self.seen / self.total) * 100
                eta = self._eta_str(rate)
                line = (
                    f"Downloading: {self._human(self.seen)} / {self._human(self.total)} "
                    f"({pct:5.1f}%) | ETA {eta} | {rate_str}"
                )
            else:
                line = f"Downloading: {self._human(self.seen)} | {rate_str}"
            sys.stdout.write(line + end)
            sys.stdout.flush()

        def finish(self):
            if show_progress:
                # Ensure a final full line with 100% if total known
                if self.total:
                    self.seen = max(self.seen, self.total)
                self._print_line(end="\n")

    # Build a queue of downloads: (key, local_path, remote_size, etag)
    download_queue: list[tuple[str, Path, int | None, str | None]] = []
    # Files present locally at the remote size: (key, local_path, size, etag)
    present: list[tuple[str, Path, int, str | None]] = []
    failed: list[str] = []

    seen: set[str] = set()

    def _consider(key: str, size: int | None, etag: str | None) -> None:
        if key in seen:
            return
        seen.add(key)
        local_path = dest_root / key
        try:
            if size is not None and local_path.exists() and local_path.stat().st_size == size:
                present.append((key, local_path, size, etag))
                return
        except OSError:
            pass
        download_queue.append((key, local_path, size, etag))

    def _head(key: str) -> tuple[int | None, str | None] | None:
        """Return (size, etag) for a key, or None when the object is absent."""
        try:
            head = s3.head_object(Bucket=bucket, Key=key)
        except Exception:
            return None
        size = int(head.get("ContentLength", -1))
        return (size if size >= 0 else None), head.get("ETag")

    # Always consider Changelog.html, silently skipping it when absent
    changelog_key = f"{prefix}Changelog.html"
    changelog = _head(changelog_key)
    if changelog is not None:
        _consider(changelog_key, *changelog)

    if files:
        # Only the requested files (plus the changelog) are considered
        for f in files:
            if not f:
                continue
            rel = str(f).lstrip("/")
            key = rel if rel.startswith(prefix) else f"{prefix}{rel}"
            head = _head(key)
            if head is None:
                print(f"[download] Not found in s3://{bucket}: {key}")
                failed.append(key)
                continue
            _consider(key, *head)
    else:
        # List all objects under the prefix and keep those matching the suffixes
        paginator = s3.get_paginator("list_objects_v2")
        for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
            for obj in page.get("Contents") or []:
                key = obj.get("Key")
                if not key or key.endswith("/") or key == changelog_key:
                    continue
                if not _matches_suffixes(key, suffixes):
                    continue
                size = obj.get("Size")
                _consider(key, size if isinstance(size, int) else None, obj.get("ETag"))

    # Execute queued downloads with progress
    total_bytes = sum(sz for _, _, sz, _ in download_queue if isinstance(sz, int)) or None
    prog = _Progress(total_bytes)
    fetched: list[tuple[str, Path, int, str | None]] = []
    for key, local_path, _sz, etag in download_queue:
        try:
            _ensure_parent_dir(local_path)
            s3.download_file(bucket, key, str(local_path), Callback=prog.cb)
        except Exception as exc:
            print(f"[download] Failed: {key} ({exc})")
            failed.append(key)
            continue
        downloaded.append(local_path)
        fetched.append((key, local_path, local_path.stat().st_size, etag))
    if download_queue:
        prog.finish()

    downloaded.skipped = [p for _, p, _, _ in present]
    downloaded.failed = failed
    downloaded.bytes_downloaded = sum(sz for _, _, sz, _ in fetched)
    downloaded.manifest_path = _update_download_manifest(
        dest_root / prefix, bucket, fetched=fetched, present=present
    )
    print(
        f"[download] {len(downloaded)} downloaded, {len(downloaded.skipped)} skipped, "
        f"{len(failed)} failed, {_human_bytes(downloaded.bytes_downloaded)} transferred"
    )
    return downloaded


def _human_bytes(n: int) -> str:
    size = float(n)
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size < 1024 or unit == "TB":
            return f"{size:.1f}{unit}"
        size /= 1024
    return f"{n}B"


def _update_download_manifest(
    prefix_dir: Path,
    bucket: str,
    *,
    fetched: list[tuple[str, Path, int, str | None]],
    present: list[tuple[str, Path, int, str | None]],
) -> Path | None:
    """Merge downloaded and already-present files into ``download_manifest.json``.

    Entries are keyed by S3 key. A freshly downloaded file gets the current
    timestamp; a skipped file keeps the timestamp of its previous entry, or
    its local modification time when it has none.
    """
    if not fetched and not present:
        return None
    path = prefix_dir / DOWNLOAD_MANIFEST_NAME
    entries: dict[str, dict] = {}
    if path.exists():
        try:
            entries = {e["key"]: e for e in json.loads(path.read_text())["files"]}
        except (ValueError, KeyError, TypeError):
            entries = {}
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    for key, _local, size, etag in fetched:
        entries[key] = {"key": key, "size": size, "etag": etag, "downloaded_at": now}
    for key, local, size, etag in present:
        prev = entries.get(key, {})
        if prev.get("size") == size and prev.get("downloaded_at"):
            stamp = prev["downloaded_at"]
        else:
            stamp = datetime.fromtimestamp(local.stat().st_mtime, timezone.utc).isoformat(
                timespec="seconds"
            )
        entries[key] = {"key": key, "size": size, "etag": etag or prev.get("etag"), "downloaded_at": stamp}
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {"bucket": bucket, "updated_at": now, "files": [entries[k] for k in sorted(entries)]}
    path.write_text(json.dumps(payload, indent=2) + "\n")
    return path


def subset_adata(file_path: str | Path, *, subclass_label: str = "Microglia-PVM"):
    """Read an AnnData H5AD file and return a subset by ``obs['Subclass']``.

    In addition to subsetting by ``Subclass == subclass_label``, this function
    applies two rudimentary QC filters to the per-donor object before returning:
    - Filter cells with too few detected genes: ``sc.pp.filter_cells(adata, min_genes=200)``
    - Filter genes detected in very few cells: ``sc.pp.filter_genes(adata, min_cells=3)``

    Parameters
    ----------
    file_path:
        Path to the ``.h5ad`` file to read.
    subclass_label:
        The label to match in ``adata.obs['Subclass']``. Defaults to
        "Microglia-PVM".

    Returns
    -------
    anndata.AnnData
        The subset ``AnnData`` object after the basic QC filters above. The
        original loaded object is deleted to free memory before returning.
    """
    # Import locally to avoid heavy import on module import
    import gc
    import anndata as ad
    import scanpy as sc

    file_path = Path(file_path)
    adata = ad.read_h5ad(str(file_path))
    if "Subclass" not in adata.obs:
        try:
            del adata
        finally:
            gc.collect()
        raise KeyError("Column 'Subclass' not found in adata.obs")

    subset = adata[adata.obs["Subclass"] == subclass_label].copy()

    # Apply basic QC on the per-donor subset
    # 1) filter cells by detected genes
    sc.pp.filter_cells(subset, min_genes=200)
    # 2) filter genes by number of cells
    sc.pp.filter_genes(subset, min_cells=3)
    # Drop the original to free memory
    try:
        del adata
    finally:
        gc.collect()
    return subset


def subset_and_concat_folder(
    folder: str | Path,
    *,
    subclass_label: str = "Microglia-PVM",
    pattern: str = "*.h5ad",
    recursive: bool = False,
):
    """Create a concatenated subset from all H5AD files in a folder.

    Iterates through files, subsetting by ``Subclass == subclass_label`` and
    additionally excluding samples where ``obs["Overall AD neuropathological Change"] == "Reference"``
    for each file. Concatenates the per-file subsets into a single AnnData while
    aggressively freeing intermediate objects to minimize memory usage.

    Parameters
    ----------
    folder:
        Directory containing input ``.h5ad`` files (e.g.
        ``data/PFC/RNAseq/donor_objects``).
    subclass_label:
        Label to match in ``obs['Subclass']`` during subsetting.
    pattern:
        Glob pattern of files to include (default: ``*.h5ad``).
    recursive:
        If True, search subdirectories recursively using ``rglob``.

    Returns
    -------
    anndata.AnnData
        Concatenation of per-file subsets after applying both filters:
        1) ``obs['Subclass'] == subclass_label`` and
        2) ``obs['Overall AD neuropathological Change'] != 'Reference'``.
    """
    import gc
    from typing import Iterable as _Iter
    import anndata as ad

    folder = Path(folder)
    if not folder.exists():
        raise FileNotFoundError(f"Folder not found: {folder}")

    files: _Iter[Path]
    files = folder.rglob(pattern) if recursive else folder.glob(pattern)
    files = sorted([p for p in files if p.is_file()])
    if not files:
        raise FileNotFoundError(f"No files matching {pattern!r} in {folder}")

    result: ad.AnnData | None = None
    for fp in files:
        subset = subset_adata(fp, subclass_label=subclass_label)
        # Apply additional subsetting condition
        if "Overall AD neuropathological Change" not in subset.obs:
            try:
                del subset
            finally:
                gc.collect()
            raise KeyError(
                "Column 'Overall AD neuropathological Change' not found in adata.obs"
            )
        subset = subset[
            subset.obs["Overall AD neuropathological Change"] != "Reference"
        ]
        if result is None:
            result = subset
        else:
            # Concatenate incrementally, then drop intermediates
            new_result = ad.concat([result, subset], join="outer", merge="unique")
            try:
                del result
                del subset
            finally:
                gc.collect()
            result = new_result

    assert result is not None  # for type checkers
    return result


def write_subset_from_folder(
    folder: str | Path = "data/PFC/RNAseq/donor_objects",
    *,
    subclass_label: str | None = None,
    out_dir: str | Path = "results/cleaned_files",
    out_name: str | None = None,
    compression: str | None = "lzf",
) -> Path | list[Path]:
    """High-level helper: subset donor files and write H5AD outputs.

    Parameters
    ----------
    folder:
        Input folder with donor ``.h5ad`` objects.
    subclass_label:
        Subclass value used for filtering. If ``None`` (default), create one
        subset file per unique value in ``adata.obs['Subclass']`` across the
        input donor files.
    out_dir:
        Output directory where the single concatenated file will be written.
    out_name:
        Output filename. Used only when ``subclass_label`` is provided. If not
        provided, a filename based on the subclass label will be generated
        automatically (e.g., ``Microglia-PVM_subset.h5ad``). Ignored when
        ``subclass_label`` is ``None``.
    compression:
        Compression passed to ``AnnData.write`` (e.g., "gzip").

    Returns
    -------
    pathlib.Path | list[pathlib.Path]
        Path to the written file when a single ``subclass_label`` is provided,
        or a list of paths when processing all subclasses.
    """
    import anndata as ad  # only for typing/IDE hints; actual write uses the object
    from collections import OrderedDict

    def _sanitize(name: str) -> str:
        # Simple filename sanitizer: keep alphanum, dash, underscore; replace others with '_'
        import re

        name = name.strip()
        name = re.sub(r"\s+", "_", name)
        name = re.sub(r"[^A-Za-z0-9._-]", "_", name)
        # collapse multiple underscores
        name = re.sub(r"_+", "_", name)
        return name

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # If a specific subclass is requested, process and write a single file
    if subclass_label is not None:
        label = str(subclass_label)
        auto_name = f"{_sanitize(label)}_subset.h5ad"
        out_path = out_dir / (out_name or auto_name)

        adata_concat = subset_and_concat_folder(folder, subclass_label=label)

        # Rename obs column before writing if present
        try:
            if "Overall AD neuropathological Change" in adata_concat.obs.columns:
                adata_concat.obs.rename(
                    columns={"Overall AD neuropathological Change": "ADNC"},
                    inplace=True,
                )
        except Exception:
            # Non-fatal; continue to write
            pass

        adata_concat.write(str(out_path), compression=compression)

        # Immediately create and write pseudobulk, then free memory
        try:
            pdata = pseudobulk_and_filter(adata_concat)
            stem = out_path.stem
            if stem.endswith("_subset"):
                stem = stem[: -len("_subset")]
            pseudobulk_path = out_dir / f"{stem}_pseudobulk.h5ad"
            pdata.write(str(pseudobulk_path), compression=compression)
        finally:
            # Ensure we free large objects aggressively
            try:
                del pdata
            except Exception:
                pass
            try:
                del adata_concat
            except Exception:
                pass
            import gc as _gc
            _gc.collect()
        return out_path

    # Otherwise, discover all subclasses and write one file per subclass
    # Use lightweight reader to avoid loading matrices and cache discovery
    folder = Path(folder)
    files = sorted([p for p in folder.glob("*.h5ad") if p.is_file()])
    if not files:
        raise FileNotFoundError(f"No files matching '*.h5ad' in {folder}")

    # Cache path (under output directory) for subclass discovery
    cache_path = out_dir / ".subclass_index.json"

    def _stat_for(p: Path) -> dict:
        try:
            st = p.stat()
            return {"size": int(st.st_size), "mtime": float(st.st_mtime)}
        except Exception:
            return {"size": None, "mtime": None}

    cached_index: dict[str, dict] | None = None
    if cache_path.exists():
        try:
            cached_index = json.loads(cache_path.read_text())
        except Exception:
            cached_index = None

    # Determine if cache is valid for current files
    cache_valid = False
    if cached_index and isinstance(cached_index, dict):
        files_map = cached_index.get("files") or {}
        try:
            cache_valid = all(
                (
                    str(fp) in files_map
                    and files_map[str(fp)].get("size") == _stat_for(fp)["size"]
                    and abs(files_map[str(fp)].get("mtime", 0.0) - _stat_for(fp)["mtime"]) < 1e-6
                )
                for fp in files
            )
        except Exception:
            cache_valid = False

    # Ordered unique subclasses for deterministic writing order
    subclasses_od: OrderedDict[str, None] = OrderedDict()
    per_file_subclasses: dict[str, list[str]] = {}

    if cache_valid:
        # Rebuild ordered dict by iterating files in current order
        for fp in files:
            vals = cached_index["files"][str(fp)].get("subclasses") or []
            vals = [str(v) for v in vals]
            per_file_subclasses[str(fp)] = vals
            for v in vals:
                subclasses_od.setdefault(v, None)
    else:
        for fp in files:
            try:
                meta = ad.read_h5ad(fp, backed='r')
                obs = meta.obs
                del meta
                vals: list[str] = []
                if obs is not None and "Subclass" in obs.columns:
                    for val in pd.unique(obs["Subclass"].astype(str)):
                        subclasses_od.setdefault(str(val), None)
                        vals.append(str(val))
                per_file_subclasses[str(fp)] = vals
            except Exception:
                # Skip problematic files, continue
                continue

        # Write cache
        try:
            cache_payload = {
                "files": {
                    str(fp): {
                        "size": _stat_for(fp)["size"],
                        "mtime": _stat_for(fp)["mtime"],
                        "subclasses": per_file_subclasses.get(str(fp), []),
                    }
                    for fp in files
                }
            }
            cache_path.write_text(json.dumps(cache_payload, indent=2))
        except Exception:
            # Cache failures are non-fatal
            pass

    if not subclasses_od:
        raise RuntimeError("No subclasses found across input files.")

    written: list[Path] = []
    for label in subclasses_od.keys():
        try:
            print(label)
            adata_concat = subset_and_concat_folder(folder, subclass_label=label)
            try:
                if "Overall AD neuropathological Change" in adata_concat.obs.columns:
                    adata_concat.obs.rename(
                        columns={"Overall AD neuropathological Change": "ADNC"},
                        inplace=True,
                    )
            except Exception:
                pass
            fname = f"{_sanitize(label)}_subset.h5ad"
            out_path = out_dir / fname
            adata_concat.write(str(out_path), compression=compression)
            written.append(out_path)
        except Exception:
            # Continue with next subclass if one fails
            continue

        # Immediately compute and write pseudobulk, then free memory
        try:
            pdata = pseudobulk_and_filter(adata_concat)
            stem = out_path.stem
            if stem.endswith("_subset"):
                stem = stem[: -len("_subset")]
            pseudobulk_path = out_dir / f"{stem}_pseudobulk.h5ad"
            pdata.write(str(pseudobulk_path), compression=compression)
        finally:
            # Ensure we free large objects aggressively
            try:
                del pdata
            except Exception:
                pass
            try:
                del adata_concat
            except Exception:
                pass
            import gc as _gc
            _gc.collect()

    if not written:
        raise RuntimeError("Failed to create any subclass subset files.")
    return written



def prep_anndata(adata_, *, min_cells: int = 3, min_genes: int = 200):
    """Prepare an AnnData object exported from R/Seurat for Scanpy workflows.

    This helper fixes common conversion issues when moving between R and Python
    and applies basic QC filtering:
    - Ensure ``obs`` column dtypes are properly interpreted (e.g. byte strings
      become UTF-8 strings) by reconstructing the object through a temporary
      ``DataFrame`` merge.
    - Filter out genes that are expressed in a small number of cells and cells
      with too few detected genes.

    Parameters
    ----------
    adata_ : anndata.AnnData
        Input AnnData object. Typically created by reading a ``.h5ad`` that was
        exported from R/Seurat or other pipelines where categorical/text columns
        may arrive as bytes or mixed dtypes.
    min_cells : int, default 3
        Keep genes that are expressed in at least ``min_cells`` cells. Passed to
        ``scanpy.pp.filter_genes`` as ``min_cells``.
    min_genes : int, default 200
        Keep cells that have at least ``min_genes`` detected genes. Passed to
        ``scanpy.pp.filter_cells`` as ``min_genes``.

    Returns
    -------
    anndata.AnnData
        A new AnnData object with fixed dtypes and filtered cells/genes.

    Notes
    -----
    - The function returns a fresh AnnData instance; the input object is not
      modified in place.
    - The dtype-fix step reconstructs ``X``, ``obs`` and ``var`` via a
      ``pandas.DataFrame`` to coerce potential byte-string columns to proper
      Python strings and to ensure consistent dtypes across columns.
    """

    def fix_dtypes(adata_):
        df = pd.DataFrame(adata_.X.A, index=adata_.obs_names, columns=adata_.var_names)
        df = df.join(adata_.obs)
        return sc.AnnData(df[adata_.var_names], obs=df.drop(columns=adata_.var_names))

    adata_ = fix_dtypes(adata_)
    sc.pp.filter_genes(adata_, min_cells=min_cells)
    sc.pp.filter_cells(adata_, min_genes=min_genes)
    return adata_


def pseudobulk_and_filter(
    adata,
    *,
    sample_col: str = "Donor ID",
    groups_col: str = "Supertype",
    mode: str = "sum",
    layer: str | None = "UMIs",
    min_cells: int = 10,
    min_counts: int = 1000,
    pseudobulk_kwargs: dict | None = None,
    filter_kwargs: dict | None = None,
):
    """Create pseudobulk profiles and filter low-quality samples using decoupler.

    This helper wraps two decoupler preprocessing utilities in sequence:
    1) ``dc.pp.pseudobulk`` to aggregate single-cell counts into sample/group
       pseudobulk profiles.
    2) ``dc.pp.filter_samples`` to filter out samples with too few cells or
       total counts.

    Parameters
    ----------
    adata : anndata.AnnData
        Input single-cell AnnData object.
    sample_col : str, default "Donor ID"
        Column in ``adata.obs`` identifying the sample (e.g., donor ID).
    groups_col : str, default "Supertype"
        Column in ``adata.obs`` defining the biological group/cluster to
        aggregate within (e.g., cell type).
    mode : str, default "sum"
        Aggregation mode for ``dc.pp.pseudobulk`` (e.g., "sum", "mean").
    layer : str | None, default "UMIs"
        Name of the layer to use for counts. If ``None``, use ``adata.X``.
    min_cells : int, default 10
        Minimum number of cells required per sample/group to keep in
        ``dc.pp.filter_samples``.
    min_counts : int, default 1000
        Minimum total counts required per sample/group to keep in
        ``dc.pp.filter_samples``.
    pseudobulk_kwargs : dict | None, optional
        Extra keyword arguments forwarded to ``dc.pp.pseudobulk``.
    filter_kwargs : dict | None, optional
        Extra keyword arguments forwarded to ``dc.pp.filter_samples``.

    Returns
    -------
    anndata.AnnData
        The pseudobulk AnnData after filtering.

    Notes
    -----
    - Requires the ``decoupler`` package to be installed and importable as
      ``decoupler``.
    """

    pb_kwargs = dict(
        adata=adata,
        sample_col=sample_col,
        groups_col=groups_col,
        mode=mode,
        layer=layer,
    )
    if pseudobulk_kwargs:
        pb_kwargs.update(pseudobulk_kwargs)

    pdata = dc.pp.pseudobulk(**pb_kwargs)

    fs_kwargs = dict(min_cells=min_cells, min_counts=min_counts)
    if filter_kwargs:
        fs_kwargs.update(filter_kwargs)

    dc.pp.filter_samples(pdata, **fs_kwargs)
    return pdata


def differential_expression(
    pseudobulk_file: str | Path = "results/cleaned_files/microglia_pseudobulk.h5ad",
    *,
    supertype: str | None = None,
    results_dir: str | Path = "results/DE",
    contrast: np.ndarray | list[float] | tuple[float, ...] = (
        0, 0, 0, 0, 0, 0, 0, 1
    ),
):
    """Run DESeq2 differential expression on pseudobulked data.

    Steps:
    - Load pseudobulk AnnData file.
    - For each category in ``obs.Supertype`` (or a specific one if provided),
      create a ``DeseqDataSet`` with design ``~Sex + ADNC + Sex:ADNC``.
    - Fit the model, compute statistics for a single contrast, and summarize.
    - Save a results table and a volcano plot per supertype.

    Parameters
    ----------
    pseudobulk_file : str | Path
        Path to the pseudobulk AnnData file (``.h5ad``).
    supertype : str | None
        If provided, limit the analysis to this value in ``obs.Supertype``.
        If ``None``, run for all supertypes present in the data.
    results_dir : str | Path, default "results/DE"
        Directory to write results into. A CSV and a PDF are produced per
        supertype as ``<Supertype>_deseq2_results.csv`` and
        ``<Supertype>_volcano.pdf``.
    contrast : array-like of float, default (0,0,0,0,0,0,0,1)
        Contrast vector passed to ``pydeseq2.DeseqStats``. Provide a sequence
        of numbers matching your design matrix. Defaults to the requested
        vector with the interaction term.

    Behavior on non full-rank design
    --------------------------------
    If ``pydeseq2`` emits a ``UserWarning`` stating that "The design matrix is
    not full rank, so the model cannot be fitted", this function skips that
    supertype and proceeds to the next one without writing outputs for it.

    Returns
    -------
    dict[str, dict[str, Path]]
        Mapping of ``supertype`` to a dict with keys ``results_csv`` and
        ``volcano_pdf``.
    """

    pseudobulk_file = Path(pseudobulk_file)
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    # Load pseudobulk AnnData
    pdata = ad.read_h5ad(str(pseudobulk_file))

    # Validate columns
    if "Supertype" not in pdata.obs.columns:
        raise KeyError("Column 'Supertype' not found in pseudobulk.obs")
    for col in ("Sex", "ADNC"):
        if col not in pdata.obs.columns:
            raise KeyError(f"Column '{col}' not found in pseudobulk.obs")

    # Determine supertypes to process
    if supertype is None:
        supertypes = pd.Index(sorted(pd.unique(pdata.obs["Supertype"].astype(str))))
    else:
        supertypes = pd.Index([supertype])

    outputs: dict[str, dict[str, Path]] = {}

    for st in supertypes:
        st_mask = pdata.obs["Supertype"].astype(str) == st
        micro_pdata = pdata[st_mask].copy()
        if micro_pdata.n_obs == 0:
            # Skip empty groups (in case of mismatched dtype)
            continue

        # Ensure categorical covariates are strings/categorical
        for col in ("Sex", "ADNC"):
            # Convert bytes to str if any
            if micro_pdata.obs[col].dtype == object:
                micro_pdata.obs[col] = micro_pdata.obs[col].apply(
                    lambda x: x.decode() if isinstance(x, (bytes, bytearray)) else x
                )
            micro_pdata.obs[col] = micro_pdata.obs[col].astype("category")

        # Create DESeq2 dataset and fit
        # If pydeseq2 emits a UserWarning that the design matrix is not full rank,
        # skip this supertype and continue with the next one.
        import warnings

        try:
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "error",
                    message=r"The design matrix is not full rank, so the model cannot be fitted.*",
                    category=UserWarning,
                )
                dds = DeseqDataSet(adata=micro_pdata, design="~Sex + ADNC + Sex:ADNC")
                # Some versions could defer the check to deseq2(); keep the same guard
                with warnings.catch_warnings():
                    warnings.filterwarnings(
                        "error",
                        message=r"The design matrix is not full rank, so the model cannot be fitted.*",
                        category=UserWarning,
                    )
                    dds.deseq2()
        except UserWarning as uw:  # design matrix not full rank
            print(
                f"[de] Skipping supertype '{st}' due to non-full-rank design: {uw}")
            continue

        # Use provided contrast vector
        contrast_vec = np.array(contrast, dtype=float)
        ds = DeseqStats(dds, contrast=contrast_vec)
        ds.summary()

        # Save results table
        safe_super = str(st).replace("/", "-")
        out_csv = results_dir / f"{safe_super}_deseq2_results.csv"
        ds.results_df.to_csv(out_csv, index=True)

        # Volcano plot per supertype with name including supertype
        volcano_pdf = results_dir / f"{safe_super}_volcano.pdf"
        dc.pl.volcano(ds.results_df, x="log2FoldChange", y="pvalue", save=str(volcano_pdf))

        outputs[str(st)] = {"results_csv": out_csv, "volcano_pdf": volcano_pdf}

    return outputs

# ---------------------------------------------------------------------------
# MTG astrocyte donor pseudobulks (RNA + ATAC)
# ---------------------------------------------------------------------------

ADNC_COLUMN = "Overall AD neuropathological Change"

_ASTRO_SUPERTYPES = ["Astro_1", "Astro_2", "Astro_3", "Astro_4", "Astro_5", "Astro_6-SEAAD"]
# Microglia / PVM supertypes, the same labels in both releases (the Micro-PVM_*
# labels bundle microglia and perivascular macrophages); the immune subclass
# also holds lymphocytes and monocytes, which are not selected. Labels and the
# decision: openspec/changes/archive/2026-10-22-mtg-microglia-positive-control/spike.md.
_MICRO_SUPERTYPES = [
    "Micro-PVM_1", "Micro-PVM_2", "Micro-PVM_2_1-SEAAD", "Micro-PVM_2_3-SEAAD", "Micro-PVM_3-SEAAD", "Micro-PVM_4-SEAAD",
]
_IMMUNE_OTHER_SUPERTYPES = ["Lymphocyte", "Monocyte"]
_SHARED_SUBCLASSES = [
    "L2/3 IT", "L4 IT", "L5 IT", "L6 IT", "L6 IT Car3", "L5 ET", "L5/6 NP", "L6 CT",
    "L6b", "Lamp5", "Lamp5 Lhx6", "Pax6", "Sncg", "Vip", "Sst", "Sst Chodl", "Pvalb",
    "Chandelier", "OPC", "Oligodendrocyte", "Endothelial",
]

# Release-wide taxonomy, keyed by the date in the file name
# (``..._final-nuclei.<date>.h5ad``): the omic and the complete ``Subclass``
# inventory of the release (labels confirmed by the spikes in
# ``openspec/changes/archive/2026-10-16-mtg-astrocyte-donor-pseudobulks/spike.md``
# and ``openspec/changes/archive/2026-10-22-mtg-microglia-positive-control/spike.md``).
# Every ``Subclass`` value of a donor object must be listed, otherwise the run
# fails. The population a step selects comes from ``CELL_TYPES`` and is
# joined to this inventory by :func:`taxonomy_for`.
TAXONOMY: dict[str, dict] = {
    "2026-06-22": {
        "omic": "rna",
        "subclasses": ["Astrocyte"] + _SHARED_SUBCLASSES + ["Immune", "VLMC & Perivascular"],
    },
    "2024-12-06": {
        "omic": "atac",
        "subclasses": ["Astrocyte"] + _SHARED_SUBCLASSES + ["Microglia-PVM", "VLMC"],
    },
}

# Cell-type registry: canonical name → the population of each release
# (``subclass``: the ``Subclass`` label(s); ``supertypes``: raw → canonical
# ``Supertype`` labels of the population; ``other_supertypes``: the
# subclass's recorded non-population supertypes, which are not selected, so
# any other label is an error), the results root (default directories of
# every MTG step, :func:`cell_type_dirs`) and the preregistered cohort cells.
# The name is a flag (``--cell-type``), the mapping is not; entries are added
# by a change proposal only. ``astrocyte`` is the default and its root is
# ``results/``, so no astrocyte path moves.
DEFAULT_CELL_TYPE = "astrocyte"
CELL_TYPES: dict[str, dict] = {
    "astrocyte": {
        "description": "astrocytes: Subclass 'Astrocyte' in both releases, supertypes Astro_1 to Astro_6-SEAAD",
        "results_root": "results",
        "releases": {
            "2026-06-22": {"subclass": ["Astrocyte"], "supertypes": {s: s for s in _ASTRO_SUPERTYPES}, "other_supertypes": []},
            "2024-12-06": {"subclass": ["Astrocyte"], "supertypes": {s: s for s in _ASTRO_SUPERTYPES}, "other_supertypes": []},
        },
        "preregistered_cohorts": {
            "cohort_30": {"n": 80, "Female": [11, 10, 28], "Male": [9, 10, 12]},
            "cohort_50": {"n": 74},
        },
    },
    "microglia": {
        "description": "microglia / PVM: RNA Subclass 'Immune', ATAC Subclass 'Microglia-PVM', restricted to the "
                       "Micro-PVM_* supertypes (lymphocytes and monocytes not selected)",
        "results_root": "results/microglia",
        "releases": {
            "2026-06-22": {
                "subclass": ["Immune"],
                "supertypes": {s: s for s in _MICRO_SUPERTYPES},
                "other_supertypes": list(_IMMUNE_OTHER_SUPERTYPES),
            },
            "2024-12-06": {
                "subclass": ["Microglia-PVM"],
                "supertypes": {s: s for s in _MICRO_SUPERTYPES},
                "other_supertypes": list(_IMMUNE_OTHER_SUPERTYPES),
            },
        },
        # The metadata-only spike's cells (spike.md), written before any microglia pseudobulk
        "preregistered_cohorts": {
            "cohort_30": {"n": 76, "Female": [11, 10, 24], "Male": [9, 11, 11]},
            "cohort_50": {"n": 63, "Female": [8, 10, 18], "Male": [8, 10, 9]},
        },
    },
}
# The cluster scripts (scripts/*.sbatch, scripts/sync_pseudobulk.sh) derive a
# cell type's pseudobulk directory from this rule without importing Python
for _name, _entry in CELL_TYPES.items():
    if _entry["results_root"] != ("results" if _name == DEFAULT_CELL_TYPE else f"results/{_name}"):
        raise AssertionError(f"CELL_TYPES[{_name!r}]['results_root'] breaks the results/<name> rule of the cluster scripts")
del _name, _entry
# Preregistered cohort cells per cell type (from the registry)
PREREGISTERED_COHORTS: dict[str, dict] = {name: entry["preregistered_cohorts"] for name, entry in CELL_TYPES.items()}

# Raw count matrix per omic (h5ad element path)
MATRIX_SOURCE: dict[str, str] = {"rna": "layers/UMIs", "atac": "X"}
RNA_MIN_GENES = 200
ADNC_STAGES: dict[str, str] = {"Not AD": "0", "Low": "0", "Intermediate": "1", "High": "2"}
COHORT_THRESHOLDS: dict[str, int] = {"cohort_30": 30, "cohort_50": 50}
EXPECTED_DONORS: dict[str, int] = {"rna": 84, "atac": 86}
STAGES = ["0", "1", "2"]
# Default output directories under a cell type's results root
CELL_TYPE_SUBDIRS: dict[str, str] = {
    "pseudobulk": "pseudobulk",
    "motco_input": "motco_input",
    "motco_analysis": "motco_analysis",
    "motco_rank_sensitivity": "motco_rank_sensitivity",
    "motco_diagnostics": "motco_diagnostics",
    "motco_latent_sensitivity": "motco_latent_sensitivity",
    "stage_de": "stage_de",
    "stage_de_limma": "stage_de_limma",
    "positive_control": "positive_control",
}


def cell_type_entry(cell_type: str) -> dict:
    """Return the registry entry of a cell type; raise listing the registered names otherwise."""
    if cell_type not in CELL_TYPES:
        raise ValueError(f"Unknown cell type {cell_type!r}; registered cell types: {list(CELL_TYPES)}")
    return CELL_TYPES[cell_type]


def cell_type_dirs(cell_type: str = DEFAULT_CELL_TYPE) -> dict[str, Path]:
    """Default directories of the MTG steps for a cell type (``<results root>/<subdir>``).

    For ``astrocyte`` the root is ``results`` and the directories are the
    ones every step used before the registry existed.
    """
    root = Path(cell_type_entry(cell_type)["results_root"])
    return {key: root / sub for key, sub in CELL_TYPE_SUBDIRS.items()}


def _resolve_dir(value: str | Path | None, cell_type: str, key: str) -> Path:
    """An explicit directory, or the cell type's default for ``key``."""
    return Path(value) if value is not None else cell_type_dirs(cell_type)[key]


def check_manifest_cell_type(manifest: dict, path: str | Path, cell_type: str) -> list[str]:
    """Check the cell type a manifest records against the requested one; return findings.

    A manifest written before the registry existed has no ``cell_type`` key;
    it is read as ``DEFAULT_CELL_TYPE`` (astrocyte) and that is returned as a
    finding.

    Raises
    ------
    ValueError
        The recorded cell type differs from ``cell_type``; the message names
        the file and both cell types.
    """
    cell_type_entry(cell_type)
    recorded = manifest.get("cell_type")
    if (DEFAULT_CELL_TYPE if recorded is None else recorded) != cell_type:
        what = f"cell type {recorded!r}" if recorded is not None else f"no cell_type (read as {DEFAULT_CELL_TYPE!r})"
        raise ValueError(
            f"{path} records {what} but {cell_type!r} was requested; outputs of one cell type are never read as "
            "another's; nothing was written"
        )
    if recorded is None:
        return [f"{path} records no cell_type (written before the cell-type registry); read as {DEFAULT_CELL_TYPE!r}"]
    return []

_RELEASE_RE = re.compile(r"_final-nuclei\.(\d{4}-\d{2}-\d{2})\.h5ad$")


def release_of(path: str | Path) -> str:
    """Return the release date parsed from a donor object file name."""
    m = _RELEASE_RE.search(Path(path).name)
    if m is None:
        raise ValueError(f"{Path(path).name}: cannot parse release date from file name")
    return m.group(1)


def donor_id_of(path: str | Path) -> str:
    """Return the donor ID encoded in a donor object file name."""
    return Path(path).name.split("_SEAAD_")[0]


def release_taxonomy(release: str, cell_type: str = DEFAULT_CELL_TYPE) -> dict:
    """The taxonomy mapping of one release for a registered cell type.

    Returns
    -------
    dict
        ``omic``, ``cell_type``, ``subclass`` (the population's ``Subclass``
        label(s)), ``other_subclasses`` (the rest of the release inventory,
        in inventory order), ``supertypes`` (raw → canonical population
        labels) and ``other_supertypes`` (the subclass's recorded
        non-population labels). For ``astrocyte`` the labels are those of the
        mapping used before the registry existed.
    """
    entry = cell_type_entry(cell_type)
    if release not in TAXONOMY:
        raise ValueError(f"no taxonomy mapping for release {release}")
    population = entry["releases"].get(release)
    if population is None:
        raise ValueError(f"cell type {cell_type!r} has no mapping for release {release} in CELL_TYPES")
    inventory = TAXONOMY[release]["subclasses"]
    unknown = [s for s in population["subclass"] if s not in inventory]
    if unknown:
        raise ValueError(f"CELL_TYPES[{cell_type!r}] names Subclass label(s) {unknown} absent from TAXONOMY[{release!r}]")
    return {
        "omic": TAXONOMY[release]["omic"],
        "cell_type": cell_type,
        "subclass": list(population["subclass"]),
        "other_subclasses": [s for s in inventory if s not in population["subclass"]],
        "supertypes": dict(population["supertypes"]),
        "other_supertypes": list(population["other_supertypes"]),
    }


def taxonomy_for(path: str | Path, cell_type: str = DEFAULT_CELL_TYPE) -> dict:
    """Return the taxonomy mapping of the release a donor object belongs to (:func:`release_taxonomy`)."""
    release = release_of(path)
    try:
        return release_taxonomy(release, cell_type)
    except ValueError as err:
        raise ValueError(f"{Path(path).name}: {err}") from None


def validate_taxonomy(
    path: str | Path,
    subclasses: Iterable[str],
    population_supertypes: Iterable[str],
    *,
    cell_type: str = DEFAULT_CELL_TYPE,
) -> None:
    """Fail on any ``Subclass`` label, or ``Supertype`` label of the population's subclass, absent from the mapping.

    Parameters
    ----------
    path:
        Donor object, used for the release and the error message.
    subclasses:
        All ``Subclass`` values of the object.
    population_supertypes:
        ``Supertype`` values of the nuclei in the population's subclass; each
        must be a population label or a recorded non-population label.
    cell_type:
        Registered cell type (default ``astrocyte``).
    """
    tax = taxonomy_for(path, cell_type)
    known = set(tax["subclass"]) | set(tax["other_subclasses"])
    for label in sorted(set(map(str, subclasses)) - known):
        raise ValueError(
            f"{Path(path).name}: Subclass label {label!r} is not in the taxonomy "
            f"mapping for release {release_of(path)}"
        )
    known = set(tax["supertypes"]) | set(tax["other_supertypes"])
    for label in sorted(set(map(str, population_supertypes)) - known):
        raise ValueError(
            f"{Path(path).name}: {cell_type} Supertype label {label!r} is not in the "
            f"taxonomy mapping for release {release_of(path)}"
        )


def _single_value(values: pd.Series, path: Path, column: str) -> str:
    uniq = pd.unique(values.astype(str))
    if len(uniq) != 1:
        raise ValueError(f"{path.name}: expected one value in obs[{column!r}], found {list(uniq)}")
    return str(uniq[0])


def _read_obs_columns(f, columns: list[str], path: Path) -> pd.DataFrame:
    from anndata.io import read_elem

    group = f["obs"]
    missing = [c for c in columns if c not in group]
    if missing:
        raise KeyError(f"{path.name}: obs lacks column(s) {missing}")
    return pd.DataFrame({c: read_elem(group[c]) for c in columns})


def _read_rows(elem, rows: np.ndarray):
    """Read selected rows of an h5ad matrix element as a CSR matrix."""
    import h5py
    import scipy.sparse as sp
    from anndata.io import sparse_dataset

    if isinstance(elem, h5py.Dataset):
        n_features = elem.shape[1]
        if rows.size == 0:
            return sp.csr_matrix((0, n_features), dtype=elem.dtype)
        return sp.csr_matrix(elem[rows])
    ds = sparse_dataset(elem)
    if rows.size == 0:
        return sp.csr_matrix((0, ds.shape[1]), dtype=ds.dtype)
    return sp.csr_matrix(ds[rows])


def _check_counts(m, path: Path, source: str) -> None:
    d = m.data
    if d.size and (not np.all(np.isfinite(d)) or d.min() < 0 or np.any(d != np.round(d))):
        raise ValueError(
            f"{path.name}: matrix {source!r} contains non-integer or negative values"
        )


def _check_features(names: np.ndarray, reference: np.ndarray, path: Path) -> None:
    """Require the same ordered feature set as the reference."""
    if len(names) == len(reference) and np.array_equal(names, reference):
        return
    if len(names) == len(reference):
        n_diff = int(np.sum(names != reference))
    else:
        n_diff = len(set(names) ^ set(reference)) or abs(len(names) - len(reference))
    raise ValueError(
        f"{path.name}: feature set differs from the reference ({n_diff} differing "
        f"features; {len(names)} vs {len(reference)} features)"
    )


def _features_hash(names: Iterable[str]) -> str:
    import hashlib

    return hashlib.sha256("\n".join(map(str, names)).encode()).hexdigest()


def select_nuclei(
    f,
    path: str | Path,
    omic: str,
    *,
    cell_type: str = DEFAULT_CELL_TYPE,
    min_genes: int = RNA_MIN_GENES,
    reference_features: np.ndarray | None = None,
    read_counts: bool = True,
) -> dict:
    """Select the nuclei of a registered cell type in one open donor object and read their raw counts.

    The one nucleus-selection code path of the MTG steps: ``pseudobulk_donor``
    sums what it returns and ``technical_donor`` summarizes the same nuclei,
    so "the same nuclei as the pseudobulk" is a property of the code, not a
    claim. Nuclei are those whose ``Subclass`` is the population's subclass
    label and whose ``Supertype`` is one of its population labels
    (``CELL_TYPES``; every astrocyte supertype for ``astrocyte``).

    Parameters
    ----------
    f:
        The donor object, opened with ``h5py.File(path, "r")``.
    path:
        Its path (release, taxonomy mapping and error messages).
    omic:
        ``"rna"`` (counts in ``layers/UMIs``) or ``"atac"`` (counts in ``X``).
    cell_type:
        Registered cell type (default ``astrocyte``).
    min_genes:
        RNA only: selected nuclei with fewer detected genes are dropped.
        ATAC applies no nucleus QC.
    reference_features:
        Ordered feature names the object must match. ``None`` skips the check.
    read_counts:
        Read the selected nuclei's rows of the count matrix. ``False`` is
        allowed for ATAC only, whose selection needs no counts (``counts`` is
        then ``None``); the RNA nucleus filter needs them.

    Returns
    -------
    dict
        ``donor_id``, ``omic``, ``release``, ``cell_type``, ``sex``, ``adnc``,
        ``excluded`` (``True`` for ``Reference`` donors), ``n_astro`` (nuclei
        of the population before QC; the name predates the registry and is
        kept so existing per-donor caches stay readable), ``var`` (feature
        annotations), ``rows`` (``obs`` positions of the selected nuclei,
        after QC), ``supertypes`` (their raw ``Supertype`` labels) and
        ``counts`` (their raw counts, nuclei × features int64 CSR; ``None``
        when not read). For an excluded donor ``rows`` and ``supertypes`` are
        empty and ``counts`` is ``None``.

    Notes
    -----
    Only ``obs`` columns, ``var`` and the selected rows of the count matrix
    are read, so the normalized RNA ``X`` is never loaded.
    """
    from anndata.io import read_elem

    path = Path(path)
    if omic not in MATRIX_SOURCE:
        raise ValueError(f"Unknown omic {omic!r}")
    tax = taxonomy_for(path, cell_type)
    if tax["omic"] != omic:
        raise ValueError(f"{path.name}: release {release_of(path)} is {tax['omic']}, not {omic}")
    if not read_counts and omic == "rna":
        raise ValueError("read_counts=False is for ATAC only: the RNA nucleus filter needs the counts")
    source = MATRIX_SOURCE[omic]

    obs = _read_obs_columns(
        f, ["Donor ID", "Sex", ADNC_COLUMN, "Subclass", "Supertype"], path
    )
    donor = _single_value(obs["Donor ID"], path, "Donor ID")
    sex = _single_value(obs["Sex"], path, "Sex")
    adnc = _single_value(obs[ADNC_COLUMN], path, ADNC_COLUMN)

    subclass = obs["Subclass"].astype(str).to_numpy()
    in_subclass = np.flatnonzero(np.isin(subclass, tax["subclass"]))
    supertype_all = obs["Supertype"].astype(str).to_numpy()[in_subclass]
    validate_taxonomy(path, subclass, supertype_all, cell_type=cell_type)
    population = np.isin(supertype_all, list(tax["supertypes"]))
    rows = in_subclass[population]
    supertype_raw = supertype_all[population]

    var = read_elem(f["var"])
    if reference_features is not None:
        _check_features(var.index.to_numpy(), reference_features, path)

    selection = {
        "donor_id": donor,
        "omic": omic,
        "release": release_of(path),
        "cell_type": cell_type,
        "sex": sex,
        "adnc": adnc,
        "excluded": adnc == "Reference",
        "n_astro": int(rows.size),
        "var": var,
        "rows": rows[:0],
        "supertypes": supertype_raw[:0],
        "counts": None,
    }
    if selection["excluded"]:
        return selection

    m = None
    if read_counts:
        if source not in f:
            raise KeyError(f"{path.name}: count matrix {source!r} not found")
        m = _read_rows(f[source], rows)
        _check_counts(m, path, source)
        m = m.astype(np.int64)
        if omic == "rna":
            n_genes = np.asarray((m > 0).sum(axis=1)).ravel()
            keep = n_genes >= min_genes
            m = m[keep]
            rows = rows[keep]
            supertype_raw = supertype_raw[keep]
    selection.update({"rows": rows, "supertypes": supertype_raw, "counts": m})
    return selection


def select_astrocyte_nuclei(f, path: str | Path, omic: str, **kwargs) -> dict:
    """:func:`select_nuclei` for ``cell_type="astrocyte"`` (the name the MTG steps used before the registry)."""
    return select_nuclei(f, path, omic, cell_type="astrocyte", **kwargs)


def pseudobulk_donor(
    path: str | Path,
    omic: str,
    *,
    cell_type: str = DEFAULT_CELL_TYPE,
    min_genes: int = RNA_MIN_GENES,
    reference_features: np.ndarray | None = None,
) -> dict:
    """Sum raw counts over the nuclei of a registered cell type in one donor object.

    Parameters
    ----------
    path:
        Donor ``.h5ad`` object (one donor per file).
    omic:
        ``"rna"`` (sums ``layers/UMIs``) or ``"atac"`` (sums ``X``).
    cell_type:
        Registered cell type (default ``astrocyte``).
    min_genes:
        RNA only: selected nuclei with fewer detected genes are dropped
        before summing. ATAC applies no nucleus QC.
    reference_features:
        Ordered feature names every donor of the omic must match. ``None``
        skips the check (first donor).

    Returns
    -------
    dict
        ``donor_id``, ``omic``, ``release``, ``cell_type``, ``sex``, ``adnc``,
        ``excluded`` (``True`` for ``Reference`` donors), ``n_astro``
        (population nuclei before QC), ``n_nuclei`` (summed),
        ``supertype_counts`` (canonical → count, all mapped supertypes),
        ``row`` (1 × features int64 CSR, ``None`` when excluded) and ``var``
        (feature annotations).

    Notes
    -----
    The nuclei are selected by :func:`select_nuclei`, the code path shared
    with the technical-covariates step; this function only sums them. Only
    ``obs`` columns, ``var`` and the selected rows of the count matrix are
    read, so the normalized RNA ``X`` is never loaded.
    """
    import h5py
    import scipy.sparse as sp

    path = Path(path)
    if omic not in MATRIX_SOURCE:
        raise ValueError(f"Unknown omic {omic!r}")
    tax = taxonomy_for(path, cell_type)
    if tax["omic"] != omic:
        raise ValueError(f"{path.name}: release {release_of(path)} is {tax['omic']}, not {omic}")

    with h5py.File(path, "r") as f:
        sel = select_nuclei(
            f, path, omic, cell_type=cell_type, min_genes=min_genes, reference_features=reference_features
        )

    record = {
        "donor_id": sel["donor_id"],
        "omic": omic,
        "release": sel["release"],
        "cell_type": cell_type,
        "sex": sel["sex"],
        "adnc": sel["adnc"],
        "excluded": sel["excluded"],
        "n_astro": sel["n_astro"],
        "n_nuclei": 0,
        "supertype_counts": {c: 0 for c in tax["supertypes"].values()},
        "row": None,
        "var": sel["var"],
    }
    if record["excluded"]:
        return record

    m = sel["counts"]
    record["n_nuclei"] = int(m.shape[0])
    record["row"] = sp.csr_matrix(np.asarray(m.sum(axis=0), dtype=np.int64))
    for raw, n in zip(*np.unique(sel["supertypes"], return_counts=True)):
        record["supertype_counts"][tax["supertypes"][str(raw)]] += int(n)
    del m, sel
    gc.collect()
    return record


# --- per-donor cache -------------------------------------------------------

def _save_cache(path: Path, record: dict, features_hash: str, source_key: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    row = record["row"]
    supertypes = list(record["supertype_counts"])
    tmp = path.with_suffix(".tmp.npz")
    np.savez_compressed(
        tmp,
        row_data=row.data if row is not None else np.zeros(0, np.int64),
        row_indices=row.indices if row is not None else np.zeros(0, np.int32),
        n_features=np.int64(row.shape[1] if row is not None else -1),
        donor_id=record["donor_id"],
        omic=record["omic"],
        release=record["release"],
        cell_type=record["cell_type"],
        sex=record["sex"],
        adnc=record["adnc"],
        excluded=record["excluded"],
        n_astro=record["n_astro"],
        n_nuclei=record["n_nuclei"],
        supertypes=np.array(supertypes, dtype=str),
        supertype_counts=np.array([record["supertype_counts"][s] for s in supertypes], np.int64),
        features_hash=features_hash,
        source_key=source_key,
    )
    os.replace(tmp, path)


def _load_cache(path: Path) -> dict:
    import scipy.sparse as sp

    with np.load(path, allow_pickle=False) as z:
        rec = {k: z[k].item() for k in (
            "donor_id", "omic", "release", "sex", "adnc", "excluded", "n_astro",
            "n_nuclei", "features_hash", "source_key", "n_features",
        )}
        # Caches written before the cell-type registry hold astrocytes
        rec["cell_type"] = z["cell_type"].item() if "cell_type" in z.files else DEFAULT_CELL_TYPE
        rec["supertype_counts"] = {
            str(s): int(n) for s, n in zip(z["supertypes"], z["supertype_counts"])
        }
        rec["row"] = None
        if not rec["excluded"]:
            data, idx = z["row_data"], z["row_indices"]
            rec["row"] = sp.csr_matrix(
                (data, idx, np.array([0, len(data)])), shape=(1, rec["n_features"])
            )
    return rec


# --- donor table and cohorts ----------------------------------------------

def recode_adnc(adnc: pd.Series, donor_ids: pd.Series) -> pd.Series:
    """Recode ADNC into the string stages ``0`` (Not AD, Low), ``1``, ``2``.

    Raises
    ------
    ValueError
        On any value outside ``ADNC_STAGES``, naming the donor and the value.
    """
    for donor, value in zip(donor_ids, adnc):
        if value not in ADNC_STAGES:
            raise ValueError(f"Donor {donor}: unexpected ADNC value {value!r}")
    return adnc.map(ADNC_STAGES)


def build_donor_table(meta: pd.DataFrame) -> tuple[pd.DataFrame, list[dict]]:
    """Build the donor table from per-donor, per-omic metadata.

    Parameters
    ----------
    meta:
        One row per donor × omic with columns ``donor_id``, ``omic``,
        ``sex``, ``adnc`` and ``n_nuclei``.

    Returns
    -------
    (pandas.DataFrame, list[dict])
        The table (``donor_id``, ``sex``, ``adnc_raw``, ``stage``, ``n_rna``,
        ``n_atac``, ``cohort_30``, ``cohort_50``), one row per non-Reference
        donor present in either omic, and the excluded Reference donors.

    Notes
    -----
    Sex and ADNC must agree between omics for donors present in both.
    Cohort flags require the threshold in **both** omics.
    """
    for donor, grp in meta.groupby("donor_id", sort=True):
        for col in ("sex", "adnc"):
            vals = dict(zip(grp["omic"], grp[col]))
            if len(set(vals.values())) > 1:
                raise ValueError(f"Donor {donor}: {col} disagrees between omics {vals}")

    first = meta.sort_values(["donor_id", "omic"]).drop_duplicates("donor_id")
    excluded = [
        {"donor_id": d, "omics": sorted(meta.loc[meta.donor_id == d, "omic"]), "reason": "ADNC Reference"}
        for d in first.loc[first.adnc == "Reference", "donor_id"]
    ]
    first = first[first.adnc != "Reference"]

    counts = meta.pivot_table(index="donor_id", columns="omic", values="n_nuclei", aggfunc="sum")
    table = pd.DataFrame({
        "donor_id": first["donor_id"].to_numpy(),
        "sex": first["sex"].to_numpy(),
        "adnc_raw": first["adnc"].to_numpy(),
    })
    table["stage"] = recode_adnc(table["adnc_raw"], table["donor_id"]).to_numpy()
    for omic in ("rna", "atac"):
        col = counts[omic] if omic in counts else pd.Series(dtype=float)
        table[f"n_{omic}"] = table["donor_id"].map(col).fillna(0).astype(int).to_numpy()
    for name, thr in COHORT_THRESHOLDS.items():
        table[name] = (table["n_rna"] >= thr) & (table["n_atac"] >= thr)
    return table.reset_index(drop=True), excluded


def check_sex_labels(table: pd.DataFrame, *, strict: bool = True) -> list[str]:
    """Require exactly two sex labels; return warnings when not strict."""
    labels = sorted(table["sex"].unique())
    if len(labels) == 2:
        return []
    msg = f"Expected exactly two sex labels, found {labels}"
    if strict or len(labels) > 2:
        raise ValueError(msg)
    return [msg]


def summarize_cohorts(table: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Count donors per sex × stage in each cohort.

    Returns
    -------
    (pandas.DataFrame, dict)
        Summary with columns ``cohort``, ``sex``, ``stage_0``, ``stage_1``,
        ``stage_2``, ``total`` (plus an ``all`` row per cohort), and per
        cohort ``n``, ``cells``, ``empty_cells`` and ``analysable`` (two sexes
        and no empty sex × stage cell).
    """
    sexes = sorted(table["sex"].unique())
    rows, info = [], {}
    for cohort in COHORT_THRESHOLDS:
        sub = table[table[cohort]]
        cells = {
            sex: [int(((sub.sex == sex) & (sub.stage == st)).sum()) for st in STAGES]
            for sex in sexes
        }
        for sex, n in cells.items():
            rows.append({"cohort": cohort, "sex": sex, **{f"stage_{st}": v for st, v in zip(STAGES, n)}, "total": sum(n)})
        tot = [sum(cells[s][i] for s in sexes) for i in range(len(STAGES))]
        rows.append({"cohort": cohort, "sex": "all", **{f"stage_{st}": v for st, v in zip(STAGES, tot)}, "total": len(sub)})
        empty = [f"{sex}/stage {st}" for sex in sexes for st, v in zip(STAGES, cells[sex]) if v == 0]
        info[cohort] = {
            "n": int(len(sub)),
            "cells": cells,
            "empty_cells": empty,
            "analysable": len(sexes) == 2 and not empty,
        }
    return pd.DataFrame(rows), info


def compare_preregistered(info: dict, cell_type: str = DEFAULT_CELL_TYPE) -> list[str]:
    """List deviations of the observed cohorts from the cell type's preregistered counts."""
    cell_type_entry(cell_type)
    findings = []
    if not PREREGISTERED_COHORTS[cell_type]:
        findings.append(f"no preregistered cohort cells for cell type {cell_type!r} in CELL_TYPES")
    for cohort, expected in PREREGISTERED_COHORTS[cell_type].items():
        obs = info.get(cohort, {})
        if obs.get("n") != expected["n"]:
            findings.append(f"{cohort}: n = {obs.get('n')}, preregistered {expected['n']}")
        for sex in ("Female", "Male"):
            if sex in expected and obs.get("cells", {}).get(sex) != expected[sex]:
                findings.append(
                    f"{cohort}: {sex} by stage = {obs.get('cells', {}).get(sex)}, "
                    f"preregistered {expected[sex]}"
                )
    return findings


# --- run -------------------------------------------------------------------

def _git_revision() -> dict:
    here = Path(__file__).resolve().parent
    try:
        rev = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=here, capture_output=True, text=True, check=True
        ).stdout.strip()
        dirty = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=no"],
            cwd=here, capture_output=True, text=True, check=True,
        ).stdout.strip() != ""
        return {"commit": rev, "dirty": dirty}
    except (OSError, subprocess.CalledProcessError):
        return {"commit": None, "dirty": None}


def _software_versions() -> dict:
    from importlib.metadata import PackageNotFoundError, version

    out = {"python": sys.version.split()[0]}
    for pkg in ("sead-sexdiff", "anndata", "numpy", "scipy", "pandas", "h5py", "boto3"):
        try:
            out[pkg] = version(pkg)
        except PackageNotFoundError:
            out[pkg] = None
    return out


def _select_donors(
    available: dict[str, dict[str, Path]],
    *,
    limit: int | None,
    donors: list[str] | None,
) -> dict[str, list[str]]:
    """Pick donors per omic: an explicit list, then the first ``limit``.

    With ``limit``, donors are taken in sorted order among those present in
    every requested omic, so a dry run covers the same donors in each.
    """
    universe = sorted(set().union(*available.values()))
    if donors:
        unknown = sorted(set(donors) - set(universe))
        if unknown:
            raise ValueError(f"Donor(s) not found under the data root: {unknown}")
        universe = [d for d in universe if d in set(donors)]
    if limit is not None:
        common = [d for d in universe if all(d in a for a in available.values())]
        universe = common[:limit]
    return {o: [d for d in universe if d in a] for o, a in available.items()}


def run_pseudobulk(
    *,
    omic: str = "both",
    cell_type: str = DEFAULT_CELL_TYPE,
    data_root: str | Path | None = None,
    out_dir: str | Path | None = None,
    limit: int | None = None,
    donors: list[str] | None = None,
    force: bool = False,
) -> dict:
    """Build donor-level pseudobulks of a registered cell type, the donor table and the manifest.

    Parameters
    ----------
    omic:
        ``"rna"``, ``"atac"`` or ``"both"``.
    cell_type:
        Registered cell type (``CELL_TYPES``; default ``astrocyte``): the
        nuclei summed and the default output directory.
    data_root:
        Data root (resolved with :func:`resolve_data_root`). Donor objects
        are read from ``<root>/<MTG_PREFIXES[omic]>``.
    out_dir:
        Output directory (default ``<results root>/pseudobulk`` of the cell
        type: ``results/pseudobulk`` for ``astrocyte``).
    limit:
        Process only the first ``limit`` donors (sorted, present in every
        requested omic). Marks the run partial.
    donors:
        Process only these donor IDs. Marks the run partial.
    force:
        Recompute cached donors and the reference feature set.

    Returns
    -------
    dict
        Paths of the written outputs plus ``cohort_info`` and ``findings``.

    Notes
    -----
    Donor objects are read one at a time and each result is cached under
    ``<out_dir>/per_donor/<omic>/<donor>.npz``; a rerun skips cached donors.
    A cached donor of another cell type is an error. Outputs:
    ``rna_counts.h5ad``, ``atac_counts.h5ad`` (donors × features, int64 CSR),
    ``donor_table.csv``, ``supertype_composition.csv``,
    ``cohort_summary.csv`` and ``manifest.json``.
    """
    import scipy.sparse as sp

    cell_type_entry(cell_type)
    root = resolve_data_root(data_root)
    out_dir = _resolve_dir(out_dir, cell_type, "pseudobulk")
    omics = ["rna", "atac"] if omic == "both" else [omic]
    print(f"[pseudobulk] Cell type: {cell_type}")
    print(f"[pseudobulk] Data root: {root.resolve()}")
    print(f"[pseudobulk] Output directory: {out_dir}")

    available: dict[str, dict[str, Path]] = {}
    for o in omics:
        folder = root / MTG_PREFIXES[o]
        files = sorted(folder.glob("*_final-nuclei.*.h5ad"))
        if not files:
            raise FileNotFoundError(f"No donor objects for {o} in {folder}")
        available[o] = {donor_id_of(p): p for p in files}
    selected = _select_donors(available, limit=limit, donors=donors)
    partial = bool(limit is not None or donors) or any(
        len(available[o]) != EXPECTED_DONORS[o] for o in omics
    ) or len(omics) < 2

    records: dict[str, list[dict]] = {}
    references: dict[str, pd.DataFrame] = {}
    for o in omics:
        cache_dir = out_dir / "per_donor" / o
        ref_path = cache_dir / "_var.h5ad"
        ref_var = None
        if ref_path.exists() and not force:
            ref_var = ad.read_h5ad(ref_path).var
        records[o] = []
        n = len(selected[o])
        for i, donor in enumerate(selected[o], 1):
            path = available[o][donor]
            cache = cache_dir / f"{donor}.npz"
            tag = f"[pseudobulk] {o} {i}/{n} {donor}:"
            if cache.exists() and not force:
                rec = _load_cache(cache)
                if rec["cell_type"] != cell_type:
                    raise ValueError(
                        f"{cache} holds cell type {rec['cell_type']!r} but {cell_type!r} was requested; use the cell "
                        "type's own output directory"
                    )
                print(f"{tag} cached ({rec['n_nuclei']} nuclei)")
            else:
                raw = pseudobulk_donor(
                    path, o, cell_type=cell_type,
                    reference_features=None if ref_var is None else ref_var.index.to_numpy(),
                )
                if raw["donor_id"] != donor:
                    raise ValueError(f"{path.name}: obs Donor ID {raw['donor_id']} differs from file name")
                if ref_var is None:
                    ref_var = raw["var"]
                    cache_dir.mkdir(parents=True, exist_ok=True)
                    ad.AnnData(var=ref_var).write_h5ad(ref_path)
                key = path.relative_to(root).as_posix()
                _save_cache(cache, raw, _features_hash(ref_var.index), key)
                rec = _load_cache(cache)
                if rec["excluded"]:
                    print(f"{tag} excluded (ADNC Reference)")
                else:
                    print(f"{tag} {rec['n_astro']} {cell_type} nuclei, {rec['n_nuclei']} summed")
                del raw
                gc.collect()
            records[o].append(rec)
        if ref_var is None:
            raise FileNotFoundError(f"No reference feature set for {o} in {ref_path}")
        ref_hash = _features_hash(ref_var.index)
        for rec in records[o]:
            if rec["features_hash"] != ref_hash:
                raise ValueError(
                    f"{cache_dir / (rec['donor_id'] + '.npz')}: cached feature set differs "
                    "from the reference; rerun with --force"
                )
        references[o] = ref_var

    meta = pd.DataFrame([
        {"donor_id": r["donor_id"], "omic": o, "sex": r["sex"], "adnc": r["adnc"], "n_nuclei": r["n_nuclei"]}
        for o in omics for r in records[o]
    ])
    table, excluded = build_donor_table(meta)
    warnings = check_sex_labels(table, strict=not partial)
    summary, cohort_info = summarize_cohorts(table)
    findings = compare_preregistered(cohort_info, cell_type) if not partial else []

    out_dir.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, Path] = {}
    table_idx = table.set_index("donor_id")
    composition = []
    for o in omics:
        kept = [r for r in records[o] if not r["excluded"]]
        X = sp.vstack([r["row"] for r in kept], format="csr").astype(np.int64)
        if X.nnz and X.data.min() < 0:
            raise ValueError(f"{o}: assembled matrix has negative values")
        obs = table_idx.loc[[r["donor_id"] for r in kept]].copy()
        obs.index.name = "donor_id"
        adata = ad.AnnData(X=X, obs=obs, var=references[o])
        path = out_dir / f"{o}_counts.h5ad"
        adata.write_h5ad(path, compression="gzip")
        outputs[f"{o}_counts"] = path
        for r in kept:
            for st, n_st in r["supertype_counts"].items():
                composition.append({"donor_id": r["donor_id"], "omic": o, "supertype": st, "n_nuclei": n_st})
        del adata, X
        gc.collect()

    comp = pd.DataFrame(composition, columns=["donor_id", "omic", "supertype", "n_nuclei"])
    totals = comp.groupby(["omic", "donor_id"])["n_nuclei"].sum()
    for (o, donor), tot in totals.items():
        if tot != table_idx.loc[donor, f"n_{o}"]:
            raise AssertionError(f"Composition total for {donor} ({o}) differs from donor table")

    outputs["donor_table"] = out_dir / "donor_table.csv"
    table.to_csv(outputs["donor_table"], index=False)
    outputs["supertype_composition"] = out_dir / "supertype_composition.csv"
    comp.to_csv(outputs["supertype_composition"], index=False)
    outputs["cohort_summary"] = out_dir / "cohort_summary.csv"
    summary.to_csv(outputs["cohort_summary"], index=False)

    def _release_files(o: str) -> list[dict]:
        dl_manifest = root / MTG_PREFIXES[o] / DOWNLOAD_MANIFEST_NAME
        entries = {}
        if dl_manifest.exists():
            entries = {e["key"]: e for e in json.loads(dl_manifest.read_text())["files"]}
        out = []
        for donor in selected[o]:
            p = available[o][donor]
            key = p.relative_to(root).as_posix()
            e = entries.get(key, {})
            out.append({
                "key": key,
                "size": p.stat().st_size,
                "etag": e.get("etag"),
                "in_download_manifest": bool(e) and e.get("size") == p.stat().st_size,
            })
        return out

    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "cell_type": cell_type,
        "partial": partial,
        "selection": {"omics": omics, "limit": limit, "donors": donors, "force": force},
        "git": _git_revision(),
        "software": _software_versions(),
        "data_root": str(root.resolve()),
        "releases": {
            o: {
                "prefix": MTG_PREFIXES[o],
                "release": sorted({r["release"] for r in records[o]}),
                "n_available": len(available[o]),
                "n_expected": EXPECTED_DONORS[o],
                "n_processed": len(records[o]),
                "n_features": int(references[o].shape[0]),
                "files": _release_files(o),
            }
            for o in omics
        },
        "taxonomy": {r: release_taxonomy(r, cell_type) for r in TAXONOMY},
        "qc": {
            "rna": {"min_genes": RNA_MIN_GENES, "applied_to": f"{cell_type} nuclei before summing"},
            "atac": {"min_genes": None, "note": "release final-nuclei QC only"},
            "nuclei_removed_by_qc": {
                o: int(sum(r["n_astro"] - r["n_nuclei"] for r in records[o] if not r["excluded"]))
                for o in omics
            },
        },
        "matrix_source": {o: MATRIX_SOURCE[o] for o in omics},
        "adnc_stages": ADNC_STAGES,
        "excluded_donors": excluded,
        "n_excluded_donors": len(excluded),
        "cohort_thresholds": COHORT_THRESHOLDS,
        "cohorts": cohort_info,
        "preregistered_cohorts": PREREGISTERED_COHORTS[cell_type],
        "findings": findings,
        "warnings": warnings,
        "outputs": {k: v.name for k, v in outputs.items()},
    }
    outputs["manifest"] = out_dir / "manifest.json"
    outputs["manifest"].write_text(json.dumps(manifest, indent=2) + "\n")

    print("[pseudobulk] Cohort summary (donors per sex x stage):")
    print(summary.to_string(index=False))
    for cohort, info in cohort_info.items():
        if not info["analysable"]:
            print(f"[pseudobulk] WARNING: {cohort} is not analysable (empty cells: {info['empty_cells']})")
    for msg in warnings + findings:
        print(f"[pseudobulk] WARNING: {msg}")
    if partial:
        print("[pseudobulk] Partial run: preregistered cohort comparison skipped")
    return {**outputs, "cohort_info": cohort_info, "findings": findings}


# ---------------------------------------------------------------------------
# MTG donor covariates (donor object metadata only)
# ---------------------------------------------------------------------------

# Source ``obs`` column → coded column of ``donor_covariates.csv``
COVARIATE_COLUMNS: dict[str, str] = {
    "Age at Death": "age_at_death",
    "PMI": "pmi",
    "RIN": "rin",
    "Years of education": "years_education",
    "Highest level of education": "education_level",
}
COVARIATE_NUMERIC = ["Age at Death", "PMI", "RIN", "Years of education"]
# Columns read per omic (``RIN`` exists only in the RNA objects)
COVARIATE_OMIC_COLUMNS: dict[str, list[str]] = {
    "rna": list(COVARIATE_COLUMNS),
    "atac": [c for c in COVARIATE_COLUMNS if c != "RIN"],
}
COVARIATE_CROSS_CHECK = ["Age at Death", "PMI", "Years of education", "Highest level of education"]
COVARIATE_RTOL = 1e-9

_CAPPED_RE = re.compile(r"^(\d+(?:\.\d+)?)\+$")
# String forms of a missing value as read from ``obs`` (NaN → ``"nan"``)
_MISSING_STRINGS = {"", "nan", "none", "<na>"}


def read_donor_covariates(path: str | Path, omic: str) -> dict[str, str]:
    """Read the donor-level covariate columns of one donor object.

    Parameters
    ----------
    path:
        Donor ``.h5ad`` object.
    omic:
        ``"rna"`` or ``"atac"``; selects the columns in ``COVARIATE_OMIC_COLUMNS``.

    Returns
    -------
    dict
        ``Donor ID`` and each covariate column → its single value, as a string.

    Notes
    -----
    Only the listed ``obs`` columns are read; ``X`` and the layers are never
    opened. A column holding more than one value fails with an error naming
    the file and the column.
    """
    import h5py

    path = Path(path)
    columns = ["Donor ID"] + COVARIATE_OMIC_COLUMNS[omic]
    with h5py.File(path, "r") as f:
        obs = _read_obs_columns(f, columns, path)
    return {c: _single_value(obs[c], path, c) for c in columns}


def code_numeric(value: str, *, allow_cap: bool = False) -> tuple[float, bool]:
    """Code one raw covariate value as a float.

    Parameters
    ----------
    value:
        Raw value as read from ``obs``.
    allow_cap:
        Accept the ``N+`` form (ages reported as ``90+``) and code it as ``N``.

    Returns
    -------
    (float, bool)
        The value (NaN when missing or non-numeric) and whether it was capped.
    """
    s = str(value).strip()
    if allow_cap:
        m = _CAPPED_RE.match(s)
        if m:
            return float(m.group(1)), True
    try:
        x = float(s)
    except ValueError:
        return float("nan"), False
    return (x if np.isfinite(x) else float("nan")), False


def _code_string(value: str) -> str | None:
    s = str(value).strip()
    return None if s.lower() in _MISSING_STRINGS else s


def code_covariates(raw: pd.DataFrame) -> tuple[pd.DataFrame, list[dict], int]:
    """Code raw RNA covariates into the ``donor_covariates.csv`` table.

    Parameters
    ----------
    raw:
        One row per donor with ``donor_id`` and the source columns of
        ``COVARIATE_COLUMNS`` as strings.

    Returns
    -------
    (pandas.DataFrame, list[dict], int)
        The table (``donor_id``, ``age_at_death``, ``age_capped``, ``pmi``,
        ``rin``, ``years_education``, ``education_level``) sorted by donor,
        the missing values (``donor_id``, ``column``, ``raw``) and the number
        of capped ages.

    Notes
    -----
    Nothing is imputed: a value that does not parse is left missing and listed.
    """
    rows, missing = [], []
    for rec in raw.sort_values("donor_id").to_dict("records"):
        row = {"donor_id": rec["donor_id"]}
        for src, name in COVARIATE_COLUMNS.items():
            if src in COVARIATE_NUMERIC:
                x, capped = code_numeric(rec[src], allow_cap=src == "Age at Death")
                row[name] = x
                if src == "Age at Death":
                    row["age_capped"] = capped
                if np.isnan(x):
                    missing.append({"donor_id": rec["donor_id"], "column": name, "raw": rec[src]})
            else:
                s = _code_string(rec[src])
                row[name] = s
                if s is None:
                    missing.append({"donor_id": rec["donor_id"], "column": name, "raw": rec[src]})
        rows.append(row)
    columns = ["donor_id", "age_at_death", "age_capped", "pmi", "rin", "years_education", "education_level"]
    table = pd.DataFrame(rows, columns=columns)
    return table, missing, int(table["age_capped"].sum())


def cross_check_covariates(rna_raw: pd.DataFrame, atac_raw: pd.DataFrame) -> list[dict]:
    """Compare the covariates shared by the RNA and ATAC objects of each donor.

    Parameters
    ----------
    rna_raw, atac_raw:
        Raw covariates per donor (``donor_id`` plus source columns, strings).

    Returns
    -------
    list[dict]
        One finding per disagreement: ``donor_id``, ``column``, ``rna``,
        ``atac`` (raw values). The RNA value is the one used.

    Notes
    -----
    Numeric columns are compared after coding (relative tolerance
    ``COVARIATE_RTOL``, capped flag included), since the ATAC objects store
    some of them as strings. Two missing values agree.
    """
    atac = atac_raw.set_index("donor_id")
    findings = []
    for rec in rna_raw.sort_values("donor_id").to_dict("records"):
        donor = rec["donor_id"]
        if donor not in atac.index:
            continue
        for col in COVARIATE_CROSS_CHECK:
            a_val = atac.loc[donor, col]
            if col in COVARIATE_NUMERIC:
                cap = col == "Age at Death"
                (x, cx), (y, cy) = code_numeric(rec[col], allow_cap=cap), code_numeric(a_val, allow_cap=cap)
                same = cx == cy and (
                    (np.isnan(x) and np.isnan(y))
                    or bool(np.isclose(x, y, rtol=COVARIATE_RTOL, atol=0.0))
                )
            else:
                same = _code_string(rec[col]) == _code_string(a_val)
            if not same:
                findings.append({"donor_id": donor, "column": col, "rna": rec[col], "atac": a_val})
    return findings


def run_covariates(
    *,
    omic: str = "both",
    cell_type: str = DEFAULT_CELL_TYPE,
    data_root: str | Path | None = None,
    out_dir: str | Path | None = None,
) -> dict:
    """Write the donor covariate table and its manifest.

    Parameters
    ----------
    omic:
        ``"both"`` (RNA values, cross-checked against ATAC) or ``"rna"``
        (no cross-check).
    cell_type:
        Registered cell type (default ``astrocyte``). The covariates are
        donor-level and do not depend on it; it sets the default output
        directory (beside the cell type's pseudobulk) and is recorded.
    data_root:
        Data root (resolved with :func:`resolve_data_root`).
    out_dir:
        Output directory (default the cell type's pseudobulk directory,
        ``results/pseudobulk`` for ``astrocyte``).

    Returns
    -------
    dict
        Paths ``covariates`` and ``manifest`` plus ``findings``.

    Notes
    -----
    Outputs: ``donor_covariates.csv`` (one row per donor with an RNA object)
    and ``covariates_manifest.json``.
    """
    if omic not in ("rna", "both"):
        raise ValueError(f"omic must be 'rna' or 'both', not {omic!r}")
    cell_type_entry(cell_type)
    root = resolve_data_root(data_root)
    out_dir = _resolve_dir(out_dir, cell_type, "pseudobulk")
    omics = ["rna", "atac"] if omic == "both" else ["rna"]
    print(f"[covariates] Cell type: {cell_type} (output directory only; the covariates are donor-level)")
    print(f"[covariates] Data root: {root.resolve()}")
    print(f"[covariates] Output directory: {out_dir}")

    raw: dict[str, pd.DataFrame] = {}
    files: dict[str, list[str]] = {}
    for o in omics:
        folder = root / MTG_PREFIXES[o]
        paths = sorted(folder.glob("*_final-nuclei.*.h5ad"))
        if not paths:
            raise FileNotFoundError(f"No donor objects for {o} in {folder}")
        recs = []
        for i, p in enumerate(paths, 1):
            rec = read_donor_covariates(p, o)
            if rec["Donor ID"] != donor_id_of(p):
                raise ValueError(f"{p.name}: obs Donor ID {rec['Donor ID']} differs from file name")
            rec["donor_id"] = rec.pop("Donor ID")
            recs.append(rec)
            print(f"[covariates] {o} {i}/{len(paths)} {rec['donor_id']}")
        raw[o] = pd.DataFrame(recs)
        files[o] = [p.relative_to(root).as_posix() for p in paths]

    table, missing, n_capped = code_covariates(raw["rna"])
    if "atac" in raw:
        findings = cross_check_covariates(raw["rna"], raw["atac"])
        rna_ids, atac_ids = set(raw["rna"].donor_id), set(raw["atac"].donor_id)
        cross_check = {
            "done": True,
            "columns": COVARIATE_CROSS_CHECK,
            "rtol": COVARIATE_RTOL,
            "n_compared": len(rna_ids & atac_ids),
            "rna_only_donors": sorted(rna_ids - atac_ids),
            "atac_only_donors": sorted(atac_ids - rna_ids),
        }
    else:
        findings = []
        cross_check = {"done": False, "note": "RNA-only run: no cross-omic check was done"}
    warnings = [
        f"{o}: {len(files[o])} donor objects read, expected {EXPECTED_DONORS[o]}"
        for o in omics if len(files[o]) != EXPECTED_DONORS[o]
    ]

    out_dir.mkdir(parents=True, exist_ok=True)
    outputs = {"covariates": out_dir / "donor_covariates.csv", "manifest": out_dir / "covariates_manifest.json"}
    table.to_csv(outputs["covariates"], index=False)
    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "cell_type": cell_type,
        "partial": bool(warnings) or len(omics) < 2,
        "git": _git_revision(),
        "software": _software_versions(),
        "data_root": str(root.resolve()),
        "value_source": "rna",
        "columns": {src: {"coded": name, "omics": [o for o in ("rna", "atac") if src in COVARIATE_OMIC_COLUMNS[o]]}
                    for src, name in COVARIATE_COLUMNS.items()},
        "coding": {
            "numeric": [COVARIATE_COLUMNS[c] for c in COVARIATE_NUMERIC],
            "string": ["education_level"],
            "age_cap": "'N+' (regex ^(\\d+(?:\\.\\d+)?)\\+$ after strip) is coded N with age_capped = true",
            "missing": "empty, nan, Unknown or any value float() cannot parse is left missing; nothing is imputed",
            "single_value": "each column must hold exactly one value per donor object",
        },
        "objects_read": {o: {"n": len(files[o]), "n_expected": EXPECTED_DONORS[o], "files": files[o]} for o in omics},
        "n_donors": int(len(table)),
        "n_age_capped": n_capped,
        "missing": missing,
        "cross_check": cross_check,
        "findings": findings,
        "warnings": warnings,
        "outputs": {k: v.name for k, v in outputs.items()},
    }
    outputs["manifest"].write_text(json.dumps(manifest, indent=2) + "\n")

    print(f"[covariates] {len(table)} donors, {n_capped} capped ages, {len(missing)} missing values")
    for m in missing:
        print(f"[covariates] Missing {m['column']} for {m['donor_id']} (raw {m['raw']!r})")
    for f in findings:
        print(f"[covariates] FINDING: {f['donor_id']} {f['column']}: RNA {f['rna']!r}, ATAC {f['atac']!r}")
    for w in warnings:
        print(f"[covariates] WARNING: {w}")
    if not cross_check["done"]:
        print(f"[covariates] {cross_check['note']}")
    return {**outputs, "findings": findings}


# ---------------------------------------------------------------------------
# MTG donor flags (SEA-AD's donor-level flag and staging labels, obs only)
# ---------------------------------------------------------------------------

# SEA-AD's donor-level obs columns → coded names. ``Severely Affected Donor``
# (DONOR_FLAG_SOURCE) is coded with DONOR_FLAG_VALUES; the staging labels are
# recorded as read and used by nothing (a later staging proposal may). The
# table is written beside, never into, donor_covariates.csv, whose hash every
# export and analysis manifest records, and it is an input of
# exclusion-sensitivity only. Flag and ask before changing.
DONOR_FLAG_COLUMNS: dict[str, str] = {
    "Severely Affected Donor": "severely_affected",
    "Cognitive Status": "cognitive_status",
    "Braak": "braak",
    "CERAD score": "cerad_score",
    "Thal": "thal",
    "LATE": "late",
}
# Read where present (the ATAC 2024-12-06 release holds it, the RNA 2026-06-22
# release does not); an absent column is recorded as a finding, not an error
DONOR_FLAG_ATAC_ONLY_COLUMNS: dict[str, str] = {"Continuous Pseudo-progression Score": "cps"}
DONOR_FLAG_SOURCE = "Severely Affected Donor"
DONOR_FLAG_VALUES: dict[str, bool] = {"Y": True, "N": False}
DONOR_FLAG_OUTPUTS: dict[str, str] = {"table": "donor_flags.csv", "manifest": "donor_flags_manifest.json"}
# The number of flagged donors SEA-AD reports (Gabitto et al. 2024: 11 of the
# 42 high-ADNC donors); a different count is recorded as a finding, never corrected
DONOR_FLAG_EXPECTED_FLAGGED = 11
# Reference donors are identified as the donor table identifies them (ADNC
# ``Reference``, which ADNC_STAGES does not map); they are listed, not tabled
DONOR_FLAG_REFERENCE_LABEL = "Reference"


def read_donor_flags(path: str | Path, omic: str) -> tuple[dict[str, str], list[str]]:
    """Read the donor-level flag and staging columns of one donor object.

    Parameters
    ----------
    path:
        Donor ``.h5ad`` object.
    omic:
        ``"rna"`` or ``"atac"`` (recorded in error messages only; the same
        columns are read from both).

    Returns
    -------
    (dict, list[str])
        ``Donor ID``, ``ADNC_COLUMN``, each column of ``DONOR_FLAG_COLUMNS`` and
        each column of ``DONOR_FLAG_ATAC_ONLY_COLUMNS`` the object holds → its
        single value as a string; and the ``DONOR_FLAG_ATAC_ONLY_COLUMNS`` the
        object lacks.

    Notes
    -----
    Only the listed ``obs`` columns are read; ``X`` and the layers are never
    opened. A column holding more than one value fails naming the file and
    the column; a missing column of ``DONOR_FLAG_COLUMNS`` (or ``Donor ID``,
    ``ADNC_COLUMN``) fails naming the file.
    """
    import h5py

    if omic not in MTG_PREFIXES:
        raise ValueError(f"omic must be one of {list(MTG_PREFIXES)}, not {omic!r}")
    path = Path(path)
    required = ["Donor ID", ADNC_COLUMN, *DONOR_FLAG_COLUMNS]
    with h5py.File(path, "r") as f:
        present = set(f["obs"].keys())
        optional = [c for c in DONOR_FLAG_ATAC_ONLY_COLUMNS if c in present]
        absent = [c for c in DONOR_FLAG_ATAC_ONLY_COLUMNS if c not in present]
        obs = _read_obs_columns(f, required + optional, path)
    return {c: _single_value(obs[c], path, c) for c in required + optional}, absent


def _flag_label(value: str) -> str | None:
    """A recorded label as read, with the missing forms (``nan``, empty) as ``None``."""
    return _code_string(value)


def run_donor_flags(
    *,
    omic: str = "both",
    cell_type: str = DEFAULT_CELL_TYPE,
    data_root: str | Path | None = None,
    out_dir: str | Path | None = None,
) -> dict:
    """Write SEA-AD's donor flag and staging labels (``donor_flags.csv``) and their manifest.

    Parameters
    ----------
    omic:
        ``"both"`` (RNA values, cross-checked against ATAC; the pseudo-progression
        score from ATAC) or ``"rna"`` (no cross-check; the score is then absent).
    cell_type:
        Registered cell type (default ``astrocyte``). The flags are
        donor-level and do not depend on it; it sets the default output
        directory (the cell type's pseudobulk directory) and is recorded.
    data_root:
        Data root (resolved with :func:`resolve_data_root`).
    out_dir:
        Output directory (default the cell type's pseudobulk directory).

    Returns
    -------
    dict
        Paths ``table`` and ``manifest`` plus ``findings``.

    Raises
    ------
    ValueError
        A non-Reference donor's ``Severely Affected Donor`` is not a key of
        ``DONOR_FLAG_VALUES`` (donor and value named), the flag or the
        Reference status disagrees between a donor's RNA and ATAC objects
        (donor and both values named), a column holds more than one value, or
        the ``Donor ID`` differs from the file name. Nothing is written.

    Notes
    -----
    Outputs: ``donor_flags.csv`` (one row per non-Reference donor: ``donor_id``,
    ``severely_affected``, ``cognitive_status``, ``braak``, ``cerad_score``,
    ``thal``, ``late``, ``cps``; labels as read, missing values empty) and
    ``donor_flags_manifest.json``. Only ``obs`` columns are read; no other
    file in the output directory is written. Reference donors (ADNC
    ``Reference``) are listed in the manifest with their labels. A label
    other than the flag that differs between the omics is a finding (both
    values recorded) and the RNA value is used; a donor with an ATAC object
    only takes its ATAC values (a finding).
    """
    if omic not in ("rna", "both"):
        raise ValueError(f"omic must be 'rna' or 'both', not {omic!r}")
    cell_type_entry(cell_type)
    root = resolve_data_root(data_root)
    out_dir = _resolve_dir(out_dir, cell_type, "pseudobulk")
    omics = ["rna", "atac"] if omic == "both" else ["rna"]
    print(f"[donor-flags] Cell type: {cell_type} (output directory only; the flags are donor-level)")
    print(f"[donor-flags] Data root: {root.resolve()}")
    print(f"[donor-flags] Output directory: {out_dir}")

    labels = {**DONOR_FLAG_COLUMNS, **DONOR_FLAG_ATAC_ONLY_COLUMNS}
    raw: dict[str, dict[str, dict[str, str]]] = {}
    files: dict[str, list[str]] = {}
    absent_columns: dict[str, dict[str, int]] = {}
    for o in omics:
        folder = root / MTG_PREFIXES[o]
        paths = sorted(folder.glob("*_final-nuclei.*.h5ad"))
        if not paths:
            raise FileNotFoundError(f"No donor objects for {o} in {folder}")
        raw[o], absent_columns[o] = {}, {}
        for i, p in enumerate(paths, 1):
            rec, absent = read_donor_flags(p, o)
            donor = rec.pop("Donor ID")
            if donor != donor_id_of(p):
                raise ValueError(f"{p.name}: obs Donor ID {donor} differs from file name")
            if donor in raw[o]:
                raise ValueError(f"{o}: donor {donor} has more than one object under {folder}")
            raw[o][donor] = rec
            for col in absent:
                absent_columns[o][col] = absent_columns[o].get(col, 0) + 1
            print(f"[donor-flags] {o} {i}/{len(paths)} {donor}")
        files[o] = [p.relative_to(root).as_posix() for p in paths]

    findings: list[str] = []
    for o in omics:
        for col, n in absent_columns[o].items():
            if n == len(files[o]):
                findings.append(f"{o}: obs lacks {col!r} in all {n} objects; {labels[col]} is not read from {o}")
            else:
                findings.append(f"{o}: obs lacks {col!r} in {n} of {len(files[o])} objects")
    if "atac" not in omics:
        findings.append(
            f"RNA-only run: no ATAC object was read, so no cross-omic check was done and "
            f"{', '.join(DONOR_FLAG_ATAC_ONLY_COLUMNS.values())} come(s) from the RNA objects only"
        )

    def is_reference(rec: dict[str, str]) -> bool:
        return str(rec[ADNC_COLUMN]).strip() == DONOR_FLAG_REFERENCE_LABEL

    rna, atac = raw["rna"], raw.get("atac", {})
    donors = sorted(set(rna) | set(atac))
    reference, disagreements, rows = [], [], []
    for donor in donors:
        r, a = rna.get(donor), atac.get(donor)
        if r is not None and a is not None and is_reference(r) != is_reference(a):
            raise ValueError(
                f"Donor {donor}: {ADNC_COLUMN} reads {r[ADNC_COLUMN]!r} in the RNA object and {a[ADNC_COLUMN]!r} in "
                "the ATAC object, so its Reference status differs between omics; nothing was written"
            )
        base = r if r is not None else a
        if is_reference(base):
            reference.append({
                "donor_id": donor,
                "omics": [o for o in omics if donor in raw[o]],
                "labels": {labels[c]: base.get(c) for c in labels if c in base},
            })
            continue
        for o, rec in (("RNA", r), ("ATAC", a)):
            if rec is not None and str(rec[DONOR_FLAG_SOURCE]).strip() not in DONOR_FLAG_VALUES:
                raise ValueError(
                    f"Donor {donor}: unexpected {DONOR_FLAG_SOURCE!r} value {rec[DONOR_FLAG_SOURCE]!r} in the {o} "
                    f"object (expected one of {list(DONOR_FLAG_VALUES)}); nothing was written"
                )
        if r is not None and a is not None:
            fr, fa = str(r[DONOR_FLAG_SOURCE]).strip(), str(a[DONOR_FLAG_SOURCE]).strip()
            if fr != fa:
                raise ValueError(
                    f"Donor {donor}: {DONOR_FLAG_SOURCE!r} is {fr!r} in the RNA object and {fa!r} in the ATAC object; "
                    "the flag must agree between omics; nothing was written"
                )
            for col in labels:
                if col == DONOR_FLAG_SOURCE or col not in r or col not in a:
                    continue
                if _flag_label(r[col]) != _flag_label(a[col]):
                    disagreements.append({"donor_id": donor, "column": col, "rna": r[col], "atac": a[col]})
        if r is None:
            findings.append(f"donor {donor} has an ATAC object only; its values are the ATAC ones")
        row = {"donor_id": donor, "severely_affected": DONOR_FLAG_VALUES[str(base[DONOR_FLAG_SOURCE]).strip()]}
        for col, name in labels.items():
            if col == DONOR_FLAG_SOURCE:
                continue
            # RNA value first; a column the RNA objects lack (the score) from ATAC
            source = r if r is not None and col in r else a
            row[name] = None if source is None or col not in source else _flag_label(source[col])
        rows.append(row)
    findings += [
        f"{d['donor_id']} {d['column']}: RNA {d['rna']!r}, ATAC {d['atac']!r} (RNA value used)" for d in disagreements
    ]
    columns = ["donor_id", "severely_affected"] + [n for c, n in labels.items() if c != DONOR_FLAG_SOURCE]
    table = pd.DataFrame(rows, columns=columns)
    n_flagged = int(table["severely_affected"].sum())
    if n_flagged != DONOR_FLAG_EXPECTED_FLAGGED:
        findings.append(
            f"{n_flagged} donors are flagged {DONOR_FLAG_SOURCE!r} = Y, not the {DONOR_FLAG_EXPECTED_FLAGGED} SEA-AD "
            "reports (Gabitto et al. 2024); recorded, not corrected"
        )
    missing = [
        {"donor_id": row["donor_id"], "column": name}
        for row in rows for name in columns[2:] if row[name] is None
    ]
    warnings = [
        f"{o}: {len(files[o])} donor objects read, expected {EXPECTED_DONORS[o]}"
        for o in omics if len(files[o]) != EXPECTED_DONORS[o]
    ]
    cross_check = (
        {
            "done": True,
            "columns": [c for c in labels if c != DONOR_FLAG_SOURCE],
            "flag": f"{DONOR_FLAG_SOURCE}: a disagreement is an error",
            "comparison": "labels compared as read (stripped; nan and empty are missing and agree)",
            "n_compared": len(set(rna) & set(atac)),
            "n_disagreements": len(disagreements),
            "rna_only_donors": sorted(set(rna) - set(atac)),
            "atac_only_donors": sorted(set(atac) - set(rna)),
        }
        if "atac" in omics else {"done": False, "note": "RNA-only run: no cross-omic check was done"}
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    outputs = {k: out_dir / v for k, v in DONOR_FLAG_OUTPUTS.items()}
    table.to_csv(outputs["table"], index=False)
    flagged = table.loc[table["severely_affected"], "donor_id"].tolist()
    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "cell_type": cell_type,
        "partial": bool(warnings) or len(omics) < 2,
        "git": _git_revision(),
        "software": _software_versions(),
        "data_root": str(root.resolve()),
        "value_source": "rna (a column the RNA objects lack from atac)",
        "columns": {
            src: {"coded": name, "omics": [o for o in omics if absent_columns[o].get(src, 0) < len(files[o])]}
            for src, name in labels.items()
        },
        "flag": {
            "column": DONOR_FLAG_SOURCE,
            "coded": DONOR_FLAG_COLUMNS[DONOR_FLAG_SOURCE],
            "values": dict(DONOR_FLAG_VALUES),
            "other_values": "an error naming the donor and the value",
            "expected_flagged": DONOR_FLAG_EXPECTED_FLAGGED,
            "expected_source": "Gabitto et al. 2024 (11 of the 42 high-ADNC donors)",
        },
        "coding": {
            "labels": "recorded as read (strings); empty, nan and none are missing and left empty; nothing is imputed",
            "single_value": "each column must hold exactly one value per donor object",
            "reference": f"donors whose {ADNC_COLUMN!r} is {DONOR_FLAG_REFERENCE_LABEL!r} are left out and listed",
        },
        "use": (
            "an input of exclusion-sensitivity only; the flag enters no model as a covariate or a variable, and no "
            "export, analysis or diagnostic step reads this table"
        ),
        "objects_read": {o: {"n": len(files[o]), "n_expected": EXPECTED_DONORS[o], "files": files[o]} for o in omics},
        "n_donors": int(len(table)),
        "n_flagged": n_flagged,
        "flagged_donors": flagged,
        "reference_donors": reference,
        "n_reference_donors": len(reference),
        "missing": missing,
        "cross_check": cross_check,
        "disagreements": disagreements,
        "findings": findings,
        "warnings": warnings,
        "outputs": {k: v.name for k, v in outputs.items()},
    }
    outputs["manifest"].write_text(json.dumps(manifest, indent=2) + "\n")

    print(f"[donor-flags] {len(table)} non-Reference donors, {n_flagged} flagged ({', '.join(flagged) or 'none'}); "
          f"{len(reference)} Reference donor(s) left out ({', '.join(d['donor_id'] for d in reference) or 'none'})")
    for msg in findings:
        print(f"[donor-flags] FINDING: {msg}")
    for w in warnings:
        print(f"[donor-flags] WARNING: {w}")
    if not cross_check["done"]:
        print(f"[donor-flags] {cross_check['note']}")
    return {**outputs, "findings": findings}


# ---------------------------------------------------------------------------
# MTG donor technical covariates (astrocyte nuclei of the raw donor objects)
# ---------------------------------------------------------------------------

# Predeclared gene sets as name patterns on the RNA gene symbols (``var``
# index), recorded as matched lists in the manifest. The ribosomal patterns
# take the cytoplasmic ribosomal protein genes (RPL*, RPS*, RPLP0-2, RPSA) and
# exclude pseudogenes (RPL13AP5), "-like" paralogs with a trailing L1 (RPL22L1)
# and kinases (RPS6KA1); the 13 mitochondrially encoded genes start with MT-;
# the hemoglobin set is HBB, HBA1, HBA2. A set that matches nothing is an error.
TECHNICAL_GENE_SETS: dict[str, tuple[str, ...]] = {
    "ribosomal": (r"^RP[LS]\d+[A-Z]?$", r"^RPLP[0-2]$", r"^RPSA$"),
    "mitochondrial": (r"^MT-",),
    "hemoglobin": (r"^HBB$", r"^HBA1$", r"^HBA2$"),
}
# ``obs`` columns read per omic for the selected astrocyte nuclei. RNA: four
# per-nucleus QC values and four library-level 10x metrics (constant within a
# library, so a donor spanning several libraries gets an astrocyte-weighted
# mean) plus ``load_name``, whose distinct values give ``n_libraries_rna``.
# ATAC: five library-level metrics (usually one library per donor object; an
# object spanning several gets the astrocyte-weighted mean and a constancy flag).
TECHNICAL_OBS_COLUMNS: dict[str, list[str]] = {
    "rna": [
        "Fraction mitochondrial UMIs",
        "Number of UMIs",
        "Genes detected",
        "Doublet score",
        "GEX_Reads_mapped_confidently_to_intronic_regions",
        "GEX_Fraction_of_transcriptomic_reads_in_cells",
        "GEX_Median_UMI_counts_per_cell",
        "GEX_sequencing_saturation",
        "load_name",
    ],
    "atac": [
        "ATAC_TSS_enrichment_score",
        "ATAC_Fraction_of_high_quality_fragments_overlapping_peaks",
        "ATAC_Fraction_of_high_quality_fragments_in_cells",
        "ATAC_Percent_duplicates",
        "ATAC_Median_high_quality_fragments_per_cell",
    ],
}
TECHNICAL_LIBRARY_COLUMN = "load_name"
TECHNICAL_OUTPUTS: dict[str, str] = {
    "table": "donor_technical.csv",
    "summary": "nucleus_technical_summary.csv",
    "manifest": "technical_manifest.json",
    "cache_dir": "per_donor_technical",
}
# Within-donor SDs are sample SDs; the supertype share (η²) is defined when at
# least two supertypes hold this many nuclei each
TECHNICAL_SD_DDOF = 1
TECHNICAL_MIN_SUPERTYPE_NUCLEI = 2


def match_gene_sets(
    names: Iterable[str], gene_sets: dict[str, tuple[str, ...]] = TECHNICAL_GENE_SETS
) -> dict[str, list[str]]:
    """Match the predeclared gene sets against feature names.

    Parameters
    ----------
    names:
        Feature names (the RNA ``var`` index: gene symbols), in matrix order.
    gene_sets:
        Set name → regular expressions; a name belongs to the set when any of
        them matches (``re.search``, so ``^``/``$`` anchor as written).

    Returns
    -------
    dict
        Set name → matched names in input order.

    Raises
    ------
    ValueError
        A set matches no name (naming the set and its patterns).
    """
    names = [str(n) for n in names]
    out: dict[str, list[str]] = {}
    for set_name, patterns in gene_sets.items():
        compiled = [re.compile(p) for p in patterns]
        matched = [n for n in names if any(c.search(n) for c in compiled)]
        if not matched:
            raise ValueError(
                f"gene set {set_name!r} matches no feature name (patterns {list(patterns)})"
            )
        out[set_name] = matched
    return out


def technical_column_name(column: str) -> str:
    """The snake_case output name of an ``obs`` QC column (``Fraction mitochondrial UMIs`` → ``fraction_mitochondrial_umis``)."""
    return re.sub(r"[^a-z0-9]+", "_", column.lower()).strip("_")


def _nan_mean_sd(values, *, ddof: int = TECHNICAL_SD_DDOF) -> tuple[float | None, float | None, int]:
    """Mean, sample SD and count of the finite values (``None`` when undefined)."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return None, None, 0
    sd = float(v.std(ddof=ddof)) if v.size > ddof else None
    return float(v.mean()), sd, int(v.size)


def _numeric_obs(series: pd.Series) -> np.ndarray:
    """An ``obs`` column as floats; categorical or string values that do not parse become NaN."""
    if isinstance(series.dtype, pd.CategoricalDtype) or series.dtype == object:
        return pd.to_numeric(series.astype(str), errors="coerce").to_numpy(dtype=float)
    return pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)


def supertype_share(
    values, groups, *, min_nuclei: int = TECHNICAL_MIN_SUPERTYPE_NUCLEI
) -> float | None:
    """Share of a per-nucleus value's variance explained by supertype within one donor (η²).

    Parameters
    ----------
    values:
        One value per nucleus (NaN values are dropped).
    groups:
        The nuclei's supertype labels.
    min_nuclei:
        The share is defined only when at least two supertypes hold this many
        nuclei each; it is then computed over all supertypes present.

    Returns
    -------
    float or None
        Between-supertype sum of squares over the total sum of squares, in
        ``[0, 1]``; ``None`` when undefined (one supertype, too few nuclei in
        every other supertype, or no variance).
    """
    v = np.asarray(values, dtype=float)
    g = np.asarray(groups).astype(str)
    ok = np.isfinite(v)
    v, g = v[ok], g[ok]
    labels, counts = np.unique(g, return_counts=True)
    if int((counts >= min_nuclei).sum()) < 2:
        return None
    grand = v.mean()
    ss_total = float(((v - grand) ** 2).sum())
    if ss_total <= 0:
        return None
    ss_between = sum(float(n) * (v[g == lab].mean() - grand) ** 2 for lab, n in zip(labels, counts))
    return float(ss_between / ss_total)


def between_donor_share(n, mean, sd, *, ddof: int = TECHNICAL_SD_DDOF) -> float | None:
    """Share of the pooled per-nucleus variance that lies between donors.

    Parameters
    ----------
    n, mean, sd:
        Per-donor nucleus count, mean and sample SD (``ddof``) of a
        per-nucleus value; a missing SD counts as no within-donor variance.

    Returns
    -------
    float or None
        ``Σ n_d (mean_d − grand)² / (that + Σ (n_d − ddof) sd_d²)``, the
        donor η² over all nuclei of all donors; ``None`` with fewer than two
        donors or no variance.
    """
    n = np.asarray(n, dtype=float)
    mean = np.asarray(mean, dtype=float)
    sd = np.asarray(sd, dtype=float)
    ok = np.isfinite(mean) & np.isfinite(n) & (n > 0)
    n, mean, sd = n[ok], mean[ok], np.where(np.isfinite(sd[ok]), sd[ok], 0.0)
    if n.size < 2:
        return None
    grand = float((n * mean).sum() / n.sum())
    between = float((n * (mean - grand) ** 2).sum())
    within = float((np.maximum(n - ddof, 0) * sd ** 2).sum())
    total = between + within
    return float(between / total) if total > 0 else None


def technical_donor(
    path: str | Path,
    omic: str,
    *,
    cell_type: str = DEFAULT_CELL_TYPE,
    gene_sets: dict[str, list[str]] | None = None,
    min_genes: int = RNA_MIN_GENES,
    reference_features: np.ndarray | None = None,
) -> dict:
    """Technical quality summaries of the nuclei of a registered cell type in one donor object.

    Parameters
    ----------
    path:
        Donor ``.h5ad`` object (one donor per file).
    omic:
        ``"rna"`` or ``"atac"``.
    gene_sets:
        RNA: set name → gene names (``match_gene_sets`` on the pseudobulk's
        feature set). ``None`` matches ``TECHNICAL_GENE_SETS`` on this
        object's ``var``.
    cell_type, min_genes, reference_features:
        Passed to :func:`select_nuclei`, so the nuclei are exactly the
        pseudobulk's of the same cell type.

    Returns
    -------
    dict
        JSON-serializable. Header (``donor_id``, ``omic``, ``release``,
        ``cell_type``, ``sex``, ``adnc``, ``excluded``, ``n_astro``
        (population nuclei before QC), ``n_nuclei``,
        ``columns``: source column → output name); for a non-excluded donor
        ``supertype_counts``, ``qc`` (output name → ``mean``, ``sd``, ``n``,
        ``n_missing``; ATAC adds ``n_unique`` and ``constant``), and for RNA
        ``gene_sets`` (set → ``n_genes``, ``mean``, ``sd``, ``n``,
        ``set_umis``, ``total_umis``, ``pseudobulk``, ``supertype_share``,
        ``n_supertypes``, ``per_supertype``), ``n_libraries``, ``libraries``
        (library → nuclei) and ``obs_consistency`` (max absolute difference
        between the ``UMIs`` row totals / nonzero counts and the ``Number of
        UMIs`` / ``Genes detected`` columns).

    Notes
    -----
    RNA: per selected nucleus, each set's UMI fraction is the set's UMIs over
    the nucleus's total UMIs in the raw ``UMIs`` layer (all genes); per donor,
    ``mean``/``sd`` are the nucleus-weighted mean and sample SD of those
    fractions, ``pseudobulk`` is the UMI-weighted fraction (sum of set UMIs
    over sum of total UMIs, the quantity ``rna_counts.h5ad`` gives), and
    ``supertype_share`` is :func:`supertype_share`. QC columns are summarized
    the same way (nucleus-weighted), so library-level metrics of a donor
    spanning several libraries become nucleus-weighted means. ATAC reads
    ``obs`` only; nothing reads the normalized RNA ``X``.
    """
    import h5py

    path = Path(path)
    if omic not in TECHNICAL_OBS_COLUMNS:
        raise ValueError(f"Unknown omic {omic!r}")
    columns = list(TECHNICAL_OBS_COLUMNS[omic])
    with h5py.File(path, "r") as f:
        sel = select_nuclei(
            f, path, omic, cell_type=cell_type, min_genes=min_genes, reference_features=reference_features,
            read_counts=omic == "rna",
        )
        qc = _read_obs_columns(f, columns, path)
    tax = taxonomy_for(path, cell_type)
    record = {
        "donor_id": sel["donor_id"],
        "omic": omic,
        "release": sel["release"],
        "cell_type": cell_type,
        "sex": sel["sex"],
        "adnc": sel["adnc"],
        "excluded": bool(sel["excluded"]),
        "n_astro": int(sel["n_astro"]),
        "n_nuclei": int(sel["rows"].size),
        "columns": {c: technical_column_name(c) for c in columns},
    }
    if record["excluded"]:
        return record

    rows = sel["rows"]
    qc = qc.iloc[rows].reset_index(drop=True)
    supertypes = np.array([tax["supertypes"][str(s)] for s in sel["supertypes"]], dtype=str)
    record["supertype_counts"] = {
        str(st): int(n) for st, n in zip(*np.unique(supertypes, return_counts=True))
    }
    qc_columns = columns
    if omic == "rna":
        m = sel["counts"]
        var_names = pd.Index(sel["var"].index.to_numpy().astype(str))
        if gene_sets is None:
            gene_sets = match_gene_sets(var_names)
        total = np.asarray(m.sum(axis=1)).ravel().astype(float)
        sets = {}
        for name, genes in gene_sets.items():
            pos = var_names.get_indexer(genes)
            if (pos < 0).any():
                missing = sorted(np.asarray(genes, dtype=str)[pos < 0].tolist())
                raise ValueError(f"{path.name}: gene set {name!r} names {missing} absent from var")
            set_umis = np.asarray(m[:, pos].sum(axis=1)).ravel().astype(float)
            frac = np.divide(set_umis, total, out=np.full_like(total, np.nan), where=total > 0)
            mean, sd, n = _nan_mean_sd(frac)
            per_supertype = {}
            for st in sorted(set(supertypes)):
                st_mean, _, st_n = _nan_mean_sd(frac[supertypes == st])
                per_supertype[st] = {"n": st_n, "mean": st_mean}
            sets[name] = {
                "n_genes": int(len(genes)),
                "mean": mean,
                "sd": sd,
                "n": n,
                "set_umis": int(set_umis.sum()),
                "total_umis": int(total.sum()),
                "pseudobulk": float(set_umis.sum() / total.sum()) if total.sum() > 0 else None,
                "supertype_share": supertype_share(frac, supertypes),
                "n_supertypes": int(len(per_supertype)),
                "per_supertype": per_supertype,
            }
        record["gene_sets"] = sets
        libraries = qc[TECHNICAL_LIBRARY_COLUMN].astype(str)
        record["n_libraries"] = int(libraries.nunique())
        record["libraries"] = {str(k): int(v) for k, v in libraries.value_counts().sort_index().items()}
        n_genes_row = np.asarray((m > 0).sum(axis=1)).ravel().astype(float)
        record["obs_consistency"] = {
            "number_of_umis_max_abs_diff": _json_number(np.nanmax(np.abs(_numeric_obs(qc["Number of UMIs"]) - total))),
            "genes_detected_max_abs_diff": _json_number(np.nanmax(np.abs(_numeric_obs(qc["Genes detected"]) - n_genes_row))),
        }
        qc_columns = [c for c in columns if c != TECHNICAL_LIBRARY_COLUMN]
        del m

    qc_out = {}
    for c in qc_columns:
        v = _numeric_obs(qc[c])
        mean, sd, n = _nan_mean_sd(v)
        entry = {"mean": mean, "sd": sd, "n": n, "n_missing": int(v.size - n)}
        if omic == "atac":
            entry["n_unique"] = int(np.unique(v[np.isfinite(v)]).size)
            entry["constant"] = entry["n_unique"] <= 1
        qc_out[technical_column_name(c)] = entry
    record["qc"] = qc_out
    del sel, qc
    gc.collect()
    return record


def _json_number(x) -> float | int | None:
    """A Python number for JSON (``None`` for NaN or infinity)."""
    if x is None:
        return None
    if isinstance(x, (int, np.integer)) and not isinstance(x, bool):
        return int(x)
    x = float(x)
    return x if np.isfinite(x) else None


def _cross_check_technical(
    donor: str, rec: dict, pb_fractions: dict, table: pd.DataFrame, pb_dir: Path,
    cell_type: str = DEFAULT_CELL_TYPE,
) -> dict:
    """Check one RNA technical record against the pseudobulk outputs (design D4).

    Raises
    ------
    ValueError
        A set's UMI-weighted fraction from the per-nucleus rows differs from
        the one computed from the donor's row of ``rna_counts.h5ad`` beyond
        ``COVARIATE_RTOL``; the selected nucleus count differs from the
        pseudobulk cache's ``n_nuclei`` or from ``donor_table.csv``; or the
        donor's presence in ``rna_counts.h5ad`` contradicts its exclusion.
        Each message names the donor and both values.
    """
    in_pb = donor in pb_fractions
    if rec["excluded"]:
        if in_pb:
            raise ValueError(
                f"{donor}: ADNC Reference in the donor object but present in {pb_dir / 'rna_counts.h5ad'}; "
                "the pseudobulk was not made from these donor objects"
            )
        return {"excluded": True}
    if not in_pb:
        raise ValueError(
            f"{donor}: not in {pb_dir / 'rna_counts.h5ad'}; the pseudobulk was not made from these donor "
            "objects (rerun `sead-sexdiff pseudobulk`)"
        )
    out = {"excluded": False, "n_nuclei": rec["n_nuclei"], "fractions": {}}
    for name, g in rec["gene_sets"].items():
        ref = pb_fractions[donor][name]
        tech = g["pseudobulk"]
        agree = tech is not None and bool(np.isclose(tech, ref, rtol=COVARIATE_RTOL, atol=0.0))
        out["fractions"][name] = {"technical": tech, "pseudobulk": ref, "agree": agree}
        if not agree:
            raise ValueError(
                f"{donor}: {name} UMI fraction {tech!r} from the per-nucleus UMIs rows differs from "
                f"{ref!r} computed from {pb_dir / 'rna_counts.h5ad'} (rtol {COVARIATE_RTOL}); the technical "
                "step did not read the nuclei the pseudobulk summed"
            )
    cache = pb_dir / "per_donor" / "rna" / f"{donor}.npz"
    out["pseudobulk_cache"] = str(cache) if cache.exists() else None
    out["pseudobulk_cache_n_nuclei"] = None
    if cache.exists():
        cached = _load_cache(cache)
        if cached["cell_type"] != cell_type:
            raise ValueError(
                f"{donor}: the pseudobulk cache {cache} holds cell type {cached['cell_type']!r} but {cell_type!r} "
                "was requested"
            )
        n_cache = int(cached["n_nuclei"])
        out["pseudobulk_cache_n_nuclei"] = n_cache
        if n_cache != rec["n_nuclei"]:
            raise ValueError(
                f"{donor}: {rec['n_nuclei']} {cell_type} nuclei selected but the pseudobulk cache {cache} "
                f"records n_nuclei = {n_cache}"
            )
    hit = table.loc[table["donor_id"] == donor, "n_rna"]
    out["donor_table_n_rna"] = int(hit.iloc[0]) if len(hit) else None
    if len(hit) and int(hit.iloc[0]) != rec["n_nuclei"]:
        raise ValueError(
            f"{donor}: {rec['n_nuclei']} {cell_type} nuclei selected but {pb_dir / 'donor_table.csv'} "
            f"records n_rna = {int(hit.iloc[0])}"
        )
    return out


def run_technical_covariates(
    omic: str = "both",
    *,
    cell_type: str = DEFAULT_CELL_TYPE,
    data_root: str | Path | None = None,
    pseudobulk_dir: str | Path | None = None,
    out_dir: str | Path | None = None,
    force: bool = False,
) -> dict:
    """Write the donor technical table of a registered cell type's nuclei and its manifest.

    Parameters
    ----------
    omic:
        ``"both"`` (RNA objects plus the ATAC library metrics) or ``"rna"``.
    cell_type:
        Registered cell type (default ``astrocyte``); the pseudobulk manifest
        must record the same one.
    data_root:
        Data root (resolved with :func:`resolve_data_root`) holding the MTG
        donor objects.
    pseudobulk_dir:
        Output directory of ``pseudobulk`` (``rna_counts.h5ad``,
        ``donor_table.csv``, ``manifest.json`` are read; the per-donor cache
        is read when present). Default: the cell type's pseudobulk directory.
    out_dir:
        Output directory (default the cell type's pseudobulk directory,
        ``results/pseudobulk`` for ``astrocyte``).
    force:
        Recompute donors already in the per-donor cache.

    Returns
    -------
    dict
        Paths ``table``, ``summary``, ``manifest`` plus ``findings`` and
        ``n_donors``.

    Raises
    ------
    FileNotFoundError
        A pseudobulk input or the donor objects are missing (the message
        names the ``pseudobulk`` command).
    ValueError
        A cross-check fails (``_cross_check_technical``), a cached record was
        made with different gene sets, QC columns, ``min_genes`` or cell type,
        the pseudobulk manifest records another cell type, or a donor
        object's ``Donor ID`` differs from its file name.

    Notes
    -----
    Per RNA donor object the cell type's nuclei are selected by
    :func:`select_nuclei` (the pseudobulk's code path), summarized
    by :func:`technical_donor` and cached under
    ``<out_dir>/per_donor_technical/rna/<donor>.json``; a rerun skips cached
    donors. Every donor's UMI-weighted set fractions are compared with the
    fractions computed from its row of ``rna_counts.h5ad`` and its nucleus
    count with the pseudobulk cache and ``donor_table.csv`` (design D4).
    Outputs: ``donor_technical.csv``, ``nucleus_technical_summary.csv`` and
    ``technical_manifest.json``; ``donor_covariates.csv`` and
    ``covariates_manifest.json`` are never opened for writing and their
    sha256 before and after the run are recorded.
    """
    import hashlib

    import scipy.sparse as sp

    log = "technical-covariates"
    if omic not in ("rna", "both"):
        raise ValueError(f"omic must be 'rna' or 'both', not {omic!r}")
    cell_type_entry(cell_type)
    root = resolve_data_root(data_root)
    pb = _resolve_dir(pseudobulk_dir, cell_type, "pseudobulk")
    out_dir = _resolve_dir(out_dir, cell_type, "pseudobulk")
    omics = ["rna", "atac"] if omic == "both" else ["rna"]
    for name in ("rna_counts.h5ad", "donor_table.csv", "manifest.json"):
        if not (pb / name).exists():
            raise FileNotFoundError(
                f"{pb / name} is missing; run `sead-sexdiff pseudobulk` (and sync) to produce it"
            )
    cov_files = {name: pb / name for name in ("donor_covariates.csv", "covariates_manifest.json")}
    cov_before = {name: _sha256_file(p) if p.exists() else None for name, p in cov_files.items()}
    pb_manifest = json.loads((pb / "manifest.json").read_text())
    cell_type_findings = check_manifest_cell_type(pb_manifest, pb / "manifest.json", cell_type)
    print(f"[{log}] Cell type: {cell_type}")
    print(f"[{log}] Data root: {root.resolve()}")
    print(f"[{log}] Pseudobulk directory: {pb}")
    print(f"[{log}] Output directory: {out_dir}")

    table = pd.read_csv(pb / "donor_table.csv", dtype={"donor_id": str, "stage": str, "sex": str})
    counts = ad.read_h5ad(pb / "rna_counts.h5ad")
    rna_names = counts.var_names.astype(str).to_numpy()
    matched = match_gene_sets(rna_names)
    gene_sets_hash = hashlib.sha256(json.dumps(matched, sort_keys=True).encode()).hexdigest()
    X = sp.csr_matrix(counts.X)
    pb_donors = counts.obs_names.astype(str).to_numpy()
    pb_total = np.asarray(X.sum(axis=1)).ravel().astype(float)
    index = pd.Index(rna_names)
    pb_set_umis = {
        name: np.asarray(X[:, index.get_indexer(genes)].sum(axis=1)).ravel().astype(float)
        for name, genes in matched.items()
    }
    pb_fractions = {
        d: {name: float(pb_set_umis[name][i] / pb_total[i]) if pb_total[i] > 0 else None for name in matched}
        for i, d in enumerate(pb_donors)
    }
    del counts, X
    for name, genes in matched.items():
        print(f"[{log}] gene set {name}: {len(genes)} genes")

    available: dict[str, dict[str, Path]] = {}
    for o in omics:
        folder = root / MTG_PREFIXES[o]
        files = sorted(folder.glob("*_final-nuclei.*.h5ad"))
        if not files:
            raise FileNotFoundError(f"No donor objects for {o} in {folder}")
        available[o] = {donor_id_of(p): p for p in files}

    expected_columns = {o: {c: technical_column_name(c) for c in TECHNICAL_OBS_COLUMNS[o]} for o in omics}
    cache_root = out_dir / TECHNICAL_OUTPUTS["cache_dir"]
    records: dict[str, dict[str, dict]] = {o: {} for o in omics}
    cross: dict[str, dict] = {}
    n_cached = {o: 0 for o in omics}
    for o in omics:
        n = len(available[o])
        for i, (donor, path) in enumerate(sorted(available[o].items()), 1):
            cache = cache_root / o / f"{donor}.json"
            tag = f"[{log}] {o} {i}/{n} {donor}:"
            if cache.exists() and not force:
                rec = json.loads(cache.read_text())
                stale = (
                    rec.get("gene_sets_sha256") != (gene_sets_hash if o == "rna" else None)
                    or rec.get("min_genes") != (RNA_MIN_GENES if o == "rna" else None)
                    or rec.get("columns") != expected_columns[o]
                )
                if stale:
                    raise ValueError(
                        f"{cache}: cached record was made with different gene sets, QC columns or min_genes; "
                        "rerun with --force"
                    )
                # Records written before the cell-type registry hold astrocytes
                cached_type = rec.get("cell_type", DEFAULT_CELL_TYPE)
                if cached_type != cell_type:
                    raise ValueError(
                        f"{cache} holds cell type {cached_type!r} but {cell_type!r} was requested; use the cell "
                        "type's own output directory"
                    )
                n_cached[o] += 1
                print(f"{tag} cached ({rec['n_nuclei']} nuclei)")
            else:
                rec = technical_donor(
                    path, o, cell_type=cell_type,
                    gene_sets=matched if o == "rna" else None,
                    reference_features=rna_names if o == "rna" else None,
                )
                if rec["donor_id"] != donor:
                    raise ValueError(f"{path.name}: obs Donor ID {rec['donor_id']} differs from file name")
                rec["gene_sets_sha256"] = gene_sets_hash if o == "rna" else None
                rec["min_genes"] = RNA_MIN_GENES if o == "rna" else None
                rec["source_key"] = path.relative_to(root).as_posix()
                cache.parent.mkdir(parents=True, exist_ok=True)
                tmp = cache.with_suffix(".tmp.json")
                tmp.write_text(json.dumps(rec, indent=1) + "\n")
                os.replace(tmp, cache)
                if rec["excluded"]:
                    print(f"{tag} excluded (ADNC Reference)")
                elif o == "rna":
                    frac = ", ".join(f"{s} {g['pseudobulk']:.4g}" for s, g in rec["gene_sets"].items())
                    print(f"{tag} {rec['n_astro']} {cell_type} nuclei, {rec['n_nuclei']} selected, {rec['n_libraries']} libraries; {frac}")
                else:
                    print(f"{tag} {rec['n_nuclei']} nuclei")
            records[o][donor] = rec
            if o == "rna":
                cross[donor] = _cross_check_technical(donor, rec, pb_fractions, table, pb, cell_type)
    for d in pb_donors:
        if d not in records["rna"]:
            raise ValueError(
                f"{pb / 'rna_counts.h5ad'} holds donor {d} but no RNA donor object was found under "
                f"{root / MTG_PREFIXES['rna']}"
            )

    # Table: one row per non-excluded RNA donor, sorted
    rna_recs = {d: r for d, r in records["rna"].items() if not r["excluded"]}
    atac_recs = {d: r for d, r in records.get("atac", {}).items() if not r["excluded"]}
    rna_qc = [expected_columns["rna"][c] for c in TECHNICAL_OBS_COLUMNS["rna"] if c != TECHNICAL_LIBRARY_COLUMN]
    atac_qc = [technical_column_name(c) for c in TECHNICAL_OBS_COLUMNS["atac"]]
    columns = ["donor_id", "n_nuclei_rna", "n_libraries_rna"]
    for s in TECHNICAL_GENE_SETS:
        columns += [f"{s}_mean", f"{s}_sd", f"{s}_pseudobulk", f"{s}_supertype_share"]
    for c in rna_qc:
        columns += [f"{c}_mean", f"{c}_sd"]
    columns += [f"{c}_mean" for c in atac_qc]
    rows, summary_rows, missing_atac, missing_values, non_constant = [], [], [], [], []
    for d in sorted(rna_recs):
        r = rna_recs[d]
        row = {"donor_id": d, "n_nuclei_rna": r["n_nuclei"], "n_libraries_rna": r["n_libraries"]}
        for s in TECHNICAL_GENE_SETS:
            g = r["gene_sets"][s]
            row[f"{s}_mean"], row[f"{s}_sd"] = g["mean"], g["sd"]
            row[f"{s}_pseudobulk"], row[f"{s}_supertype_share"] = g["pseudobulk"], g["supertype_share"]
            summary_rows.append({
                "donor_id": d, "gene_set": s, "n_nuclei": g["n"], "mean": g["mean"], "sd": g["sd"],
                "pseudobulk": g["pseudobulk"], "supertype_share": g["supertype_share"],
                "n_supertypes": g["n_supertypes"],
            })
        for c in rna_qc:
            q = r["qc"][c]
            row[f"{c}_mean"], row[f"{c}_sd"] = q["mean"], q["sd"]
            if q["n_missing"]:
                missing_values.append({"donor_id": d, "omic": "rna", "column": c, "n_missing": q["n_missing"]})
        a = atac_recs.get(d)
        if a is None:
            if "atac" in omics:
                missing_atac.append(d)
            for c in atac_qc:
                row[f"{c}_mean"] = None
        else:
            for c in atac_qc:
                q = a["qc"][c]
                row[f"{c}_mean"] = q["mean"]
                if q["n_missing"]:
                    missing_values.append({"donor_id": d, "omic": "atac", "column": c, "n_missing": q["n_missing"]})
                if not q["constant"]:
                    non_constant.append({"donor_id": d, "column": c, "n_unique": q["n_unique"]})
        rows.append(row)
    tech = pd.DataFrame(rows, columns=columns)
    summary = pd.DataFrame(
        summary_rows,
        columns=["donor_id", "gene_set", "n_nuclei", "mean", "sd", "pseudobulk", "supertype_share", "n_supertypes"],
    )
    atac_only = sorted(set(atac_recs) - set(rna_recs))
    excluded = [
        {"donor_id": d, "omic": o, "reason": "ADNC Reference"}
        for o in omics for d, r in sorted(records[o].items()) if r["excluded"]
    ]

    # Between-donor shares of the per-nucleus values, pooled over the table's donors
    between = {}
    for s in TECHNICAL_GENE_SETS:
        sub = summary[summary["gene_set"] == s]
        between[s] = between_donor_share(sub["n_nuclei"], sub["mean"], sub["sd"])
    for c in rna_qc:
        between[c] = between_donor_share(
            [rna_recs[d]["qc"][c]["n"] for d in sorted(rna_recs)],
            [rna_recs[d]["qc"][c]["mean"] for d in sorted(rna_recs)],
            [rna_recs[d]["qc"][c]["sd"] for d in sorted(rna_recs)],
        )
    undefined_share = {
        s: sorted(d for d in rna_recs if rna_recs[d]["gene_sets"][s]["supertype_share"] is None)
        for s in TECHNICAL_GENE_SETS
    }

    findings: list[str] = list(cell_type_findings)
    warnings: list[str] = []
    for o in omics:
        if len(available[o]) != EXPECTED_DONORS[o]:
            warnings.append(f"{o}: {len(available[o])} donor objects read, expected {EXPECTED_DONORS[o]}")
    if missing_atac:
        findings.append(f"{len(missing_atac)} donor(s) with an RNA object have no ATAC object: {missing_atac}")
    if atac_only:
        findings.append(f"{len(atac_only)} non-Reference donor(s) have an ATAC object but no RNA object (no technical row): {atac_only}")
    for s, donors in undefined_share.items():
        if donors:
            findings.append(
                f"{s}: supertype share undefined for {len(donors)} donor(s) (fewer than two supertypes with "
                f">= {TECHNICAL_MIN_SUPERTYPE_NUCLEI} nuclei, or no variance): {donors}"
            )
    if missing_values:
        findings.append(f"{len(missing_values)} donor × column cell(s) with missing QC values (listed in the manifest)")
    if non_constant:
        findings.append(f"{len(non_constant)} ATAC donor × column cell(s) are not constant within the object (listed in the manifest)")
    partial = bool(warnings) or len(omics) < 2 or bool(pb_manifest.get("partial", False))

    out_dir.mkdir(parents=True, exist_ok=True)
    outputs = {k: out_dir / TECHNICAL_OUTPUTS[k] for k in ("table", "summary", "manifest")}
    for key, frame in (("table", tech), ("summary", summary)):
        tmp = outputs[key].with_suffix(".tmp.csv")
        frame.to_csv(tmp, index=False, float_format=_csv_float)
        os.replace(tmp, outputs[key])
    cov_after = {name: _sha256_file(p) if p.exists() else None for name, p in cov_files.items()}
    if cov_after != cov_before:
        raise AssertionError(f"covariate files changed during the run: before {cov_before}, after {cov_after}")

    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "cell_type": cell_type,
        "partial": partial,
        "selection": {"omics": omics, "force": force},
        "git": _git_revision(),
        "software": _software_versions(),
        "data_root": str(root.resolve()),
        "pseudobulk": {
            "dir": str(pb.resolve()),
            "manifest_sha256": _sha256_file(pb / "manifest.json"),
            "manifest_created_at": pb_manifest.get("created_at"),
            "manifest_partial": pb_manifest.get("partial"),
            "rna_counts_sha256": _sha256_file(pb / "rna_counts.h5ad"),
            "donor_table_sha256": _sha256_file(pb / "donor_table.csv"),
        },
        "covariate_files": {
            name: {"sha256_before": cov_before[name], "sha256_after": cov_after[name], "written": False}
            for name in cov_files
        },
        "nucleus_selection": {
            "function": "select_nuclei (shared with pseudobulk_donor)",
            "cell_type": cell_type,
            "rna": {"matrix": MATRIX_SOURCE["rna"], "min_genes": RNA_MIN_GENES},
            "atac": {"matrix": None, "note": "obs columns only; the ATAC count matrix is not read"},
            "reference_excluded": True,
            "taxonomy": {r: release_taxonomy(r, cell_type) for r in TAXONOMY},
        },
        "gene_sets": {
            name: {"patterns": list(TECHNICAL_GENE_SETS[name]), "n_genes": len(genes), "genes": genes}
            for name, genes in matched.items()
        },
        "gene_sets_sha256": gene_sets_hash,
        "columns": {
            o: {c: technical_column_name(c) for c in TECHNICAL_OBS_COLUMNS[o]} for o in omics
        },
        "library_column": TECHNICAL_LIBRARY_COLUMN,
        "definitions": {
            "fraction": "per nucleus: UMIs of the set over the nucleus's total UMIs in layers/UMIs (all genes)",
            "<set>_mean, <set>_sd": f"nucleus-weighted mean and sample SD (ddof {TECHNICAL_SD_DDOF}) of the per-nucleus fractions",
            "<set>_pseudobulk": "sum of set UMIs over sum of total UMIs across the nuclei (the fraction rna_counts.h5ad gives)",
            "<set>_supertype_share": (
                "within-donor eta squared: between-supertype sum of squares of the per-nucleus fraction over its total "
                f"sum of squares, over the supertypes present; defined when at least two supertypes hold >= "
                f"{TECHNICAL_MIN_SUPERTYPE_NUCLEI} nuclei each, otherwise empty"
            ),
            "<qc>_mean, <qc>_sd": f"nucleus-weighted mean and sample SD of the obs column over the selected nuclei (library-level columns become {cell_type}-weighted means)",
            "n_libraries_rna": f"distinct values of obs['{TECHNICAL_LIBRARY_COLUMN}'] among the selected nuclei",
            "atac <qc>_mean": f"nucleus-weighted mean over the selected ATAC {cell_type} nuclei; 'constant' per donor and column in 'atac_constancy'",
            "between_donor_share": (
                "share of the pooled per-nucleus variance between donors, over the table's donors: "
                "sum_d n_d (mean_d - grand)^2 / (that + sum_d (n_d - ddof) sd_d^2)"
            ),
            "missing": "non-numeric or missing obs values are dropped from the means and SDs and listed; nothing is imputed",
        },
        "objects_read": {
            o: {
                "prefix": MTG_PREFIXES[o],
                "n": len(available[o]),
                "n_expected": EXPECTED_DONORS[o],
                "n_from_cache": n_cached[o],
                "files": [p.relative_to(root).as_posix() for _, p in sorted(available[o].items())],
            }
            for o in omics
        },
        "n_donors": int(len(tech)),
        "excluded_donors": excluded,
        "missing_atac_donors": missing_atac,
        "atac_only_donors": atac_only,
        "donors": {
            d: {
                "n_astro": rna_recs[d]["n_astro"],
                "n_nuclei": rna_recs[d]["n_nuclei"],
                "pseudobulk_cache_n_nuclei": cross[d].get("pseudobulk_cache_n_nuclei"),
                "donor_table_n_rna": cross[d].get("donor_table_n_rna"),
                "n_libraries": rna_recs[d]["n_libraries"],
                "libraries": rna_recs[d]["libraries"],
                "supertype_counts": rna_recs[d]["supertype_counts"],
                "obs_consistency": rna_recs[d]["obs_consistency"],
                "cross_check": cross[d]["fractions"],
                "n_nuclei_atac": atac_recs[d]["n_nuclei"] if d in atac_recs else None,
            }
            for d in sorted(rna_recs)
        },
        "cross_check": {
            "rtol": COVARIATE_RTOL,
            "n_donors": len(rna_recs),
            "fractions_agree": all(f["agree"] for d in rna_recs for f in cross[d]["fractions"].values()),
            "n_with_pseudobulk_cache": sum(1 for d in rna_recs if cross[d].get("pseudobulk_cache_n_nuclei") is not None),
            "nucleus_counts_agree": True,
            "note": "any disagreement stops the run, so a written manifest records agreement",
        },
        "n_libraries_distribution": {
            str(k): int(v) for k, v in sorted(pd.Series([rna_recs[d]["n_libraries"] for d in rna_recs]).value_counts().items())
        },
        "between_donor_share": between,
        "supertype_share_undefined": undefined_share,
        "atac_constancy": {
            d: {c: atac_recs[d]["qc"][c]["constant"] for c in atac_qc} for d in sorted(atac_recs) if d in rna_recs
        },
        "non_constant_atac": non_constant,
        "missing_values": missing_values,
        "cache": {
            "dir": str(cache_root),
            "note": "one JSON per donor and omic with the summaries and per-supertype counts and means (no per-nucleus arrays); skipped on rerun unless --force",
        },
        "findings": findings,
        "warnings": warnings,
        "outputs": {k: v.name for k, v in outputs.items()},
    }
    tmp = outputs["manifest"].with_suffix(".tmp.json")
    tmp.write_text(json.dumps(manifest, indent=2) + "\n")
    os.replace(tmp, outputs["manifest"])

    print(f"[{log}] {len(tech)} donors with a technical row; {len(excluded)} excluded; "
          f"{len(missing_atac)} missing ATAC; libraries per donor {manifest['n_libraries_distribution']}")
    for s in TECHNICAL_GENE_SETS:
        col = tech[f"{s}_pseudobulk"]
        share = tech[f"{s}_supertype_share"]
        print(
            f"[{log}] {s}: pseudobulk fraction median {col.median():.4g} (range {col.min():.4g}–{col.max():.4g}); "
            f"between-donor share {_fmt(between[s], '.3f')}; median supertype share {_fmt(share.median(), '.3f')} "
            f"({int(share.notna().sum())} donors defined)"
        )
    for msg in findings:
        print(f"[{log}] FINDING: {msg}")
    for msg in warnings:
        print(f"[{log}] WARNING: {msg}")
    return {**outputs, "findings": findings, "n_donors": int(len(tech))}


# ---------------------------------------------------------------------------
# MOTCO input export (per-cohort preprocessing of the pseudobulks)
# ---------------------------------------------------------------------------

EXPORT_EXCLUDED_CHROMOSOMES = ("chrX", "chrY", "chrM")
EXPORT_MIN_CPM = 1.0
EXPORT_PRIOR_COUNT = 2.0
EXPORT_TOP_N = 5000
EXPORT_COHORTS = ("cohort_30", "cohort_50")
STAGE_LABELS: dict[str, str] = {"0": "Not AD + Low", "1": "Intermediate", "2": "High"}
# edgeR calcNormFactors(method = "TMM") defaults
TMM_LOGRATIO_TRIM = 0.3
TMM_SUM_TRIM = 0.05
EXPORT_FLOAT_FORMAT = "%.10g"
# Selection rules of the export's step 6 (``preprocess_block``): the primary's
# ``variance`` and the two rules of the preprocessing-sensitivity registry.
# ``trend_ratio`` divides each feature's log2-CPM variance by a mean–variance
# trend: the median variance of TREND_N_BINS equal-count bins of mean log2-CPM
# (np.array_split of the features sorted by mean), each placed at its median
# mean, linearly interpolated (np.interp; flat beyond the end bins). ``mean``
# ranks by mean log2-CPM. Ties by name under every rule. Constants of the
# feature-selection sensitivity change; flag and ask before changing them.
SELECTION_RULES = ("variance", "trend_ratio", "mean")
TREND_N_BINS = 20
# Residualization (the ``residualized`` and ``revised_candidate`` variants):
# after log2-CPM and before the selection step each abundance-passing feature
# is replaced by its mean plus the residual of its OLS regression on an
# intercept and RESIDUALIZATION_COVARIATES. The per-omic log10 library sizes
# are the export's own (the step-3 library sizes: each donor's counts over the
# features kept at step 2); the gene-set pseudobulk fractions come from
# ``donor_technical.csv``, minus the fractions the technical-axis report
# classifies as cell state under its fixed rule (within-donor supertype share
# above the between-donor share), recorded here with the report's section as
# the reason. Sex, stage and the analysis covariates never enter the set.
RESIDUALIZATION_CANDIDATE_FRACTIONS = ("ribosomal_pseudobulk", "mitochondrial_pseudobulk", "hemoglobin_pseudobulk")
RESIDUALIZATION_CELL_STATE_EXCLUSIONS: dict[str, str] = {
    "hemoglobin_pseudobulk": (
        "reports/motco-astrocyte-technical-axis-2026-10-19.md, section 2, 'Classification under the spec's rule': "
        "median within-donor supertype share 0.00394 above the between-donor share 0.00285 (both negligible; "
        "median pseudobulk fraction 6.7e-07 of UMIs); not a residualization candidate"
    ),
}
RESIDUALIZATION_COVARIATES = ("log10_lib_rna", "log10_lib_atac") + tuple(
    f for f in RESIDUALIZATION_CANDIDATE_FRACTIONS if f not in RESIDUALIZATION_CELL_STATE_EXCLUSIONS
)
RESIDUALIZATION_FILE = "residualization.csv"
# Trend normalization (the ``trend_normalized`` and ``revised_candidate_trend``
# variants): after log2-CPM and before the residualization and the selection
# step, limma's fast cyclic loess (``normalizeCyclicLoess(method = "fast")``)
# over the abundance-passing features of a block. At each of LOESS_ITERATIONS
# iterations the reference is the per-feature mean over donors (A); each
# donor's deviation from it (M) is fitted by a lowess curve of M on A and the
# fit is subtracted. The lowess is statsmodels' (Cleveland's algorithm, as R's
# ``lowess``, which limma's ``loessFit`` calls without weights): frac = span,
# LOESS_ROBUSTNESS_ITERATIONS robustifying iterations (``loessFit``'s default
# ``iterations = 4`` passes ``iter = 3``), delta = LOESS_DELTA_FRACTION of the
# range of A (R's default). The span is limma's adaptive default,
# ``chooseLowessSpan(n, small.n, min.span, power)`` = min(min_span + (1 −
# min_span) · (small_n / n) ** power, 1) with n the block's number of
# abundance-passing features, so each block gets its own span. No label and
# no donor covariate enters; ``none`` is the primary's. The per-donor tilt
# accounting (slope of M on A, median M of the lowest and highest of
# NORMALIZATION_TILT_BINS equal-count bins of A, before and after) is
# label-free. Cross-checked against limma by
# ``scripts/check_cyclic_loess_limma.R``. Constants of the trend-normalization
# sensitivity change; flag and ask before changing them.
NORMALIZATION_RULES = ("none", "cyclic_loess")
LOESS_ITERATIONS = 3
LOESS_SPAN_RULE: dict[str, float] = {"small_n": 50, "min_span": 0.3, "power": 1 / 3}
LOESS_ROBUSTNESS_ITERATIONS = 3
LOESS_DELTA_FRACTION = 0.01
NORMALIZATION_TILT_BINS = 20
NORMALIZATION_TREND_FILE = "normalization_trend.csv"
NORMALIZATION_TILT_STATISTICS = ("slope", "low_bin", "high_bin")

# GENCODE v32 = the annotation of the 10x GRCh38-2020-A reference (Cell Ranger 6.0)
GENCODE_RELEASE = "32"
GENCODE_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/"
    "gencode.v32.basic.annotation.gtf.gz"
)
ANNOTATION_DIR = "annotation"
UNMAPPED_CHROMOSOME = "unmapped"

_PEAK_RE = re.compile(r"^(chr[0-9A-Za-z_]+):(\d+)-(\d+)$")
_GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')


def _sha256_file(path: Path) -> str:
    import hashlib

    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_gtf_gene_chromosomes(path: str | Path) -> dict[str, str]:
    """Map version-stripped Ensembl gene ids to chromosomes from a GTF.

    Parameters
    ----------
    path:
        GENCODE GTF (gzip or plain).

    Returns
    -------
    dict
        ``gene_id`` (without ``.N`` version) → ``seqname``, from ``gene`` lines.

    Notes
    -----
    GENCODE lists the pseudoautosomal genes twice, the chrY copy with a
    ``_PAR_Y`` suffix; those copies are skipped so PAR genes map to chrX.
    """
    import gzip

    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    mapping: dict[str, str] = {}
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 8)
            if len(fields) < 9 or fields[2] != "gene":
                continue
            m = _GENE_ID_RE.search(fields[8])
            if m is None:
                raise ValueError(f"{path.name}: gene line without gene_id: {line[:80]!r}")
            gid = m.group(1)
            if gid.endswith("_PAR_Y"):
                continue
            base = gid.split(".")[0]
            if mapping.get(base, fields[0]) != fields[0]:
                raise ValueError(f"{path.name}: gene {base} on {mapping[base]} and {fields[0]}")
            mapping[base] = fields[0]
    return mapping


def resolve_gene_annotation(data_root: str | Path | None = None) -> tuple[dict[str, str], dict]:
    """Return the GENCODE gene → chromosome map, downloading the GTF once.

    Parameters
    ----------
    data_root:
        Data root (resolved with :func:`resolve_data_root`). The GTF lives at
        ``<root>/annotation/gencode.v32.basic.annotation.gtf.gz``.

    Returns
    -------
    (dict, dict)
        The map from :func:`parse_gtf_gene_chromosomes` and its provenance
        (``source_url``, ``release``, ``file``, ``sha256``, ``n_genes``,
        ``downloaded``).
    """
    import urllib.request

    root = resolve_data_root(data_root)
    path = root / ANNOTATION_DIR / GENCODE_URL.rsplit("/", 1)[1]
    downloaded = False
    if not path.exists():
        print(f"[export] Downloading {GENCODE_URL}")
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(path.name + ".part")
        urllib.request.urlretrieve(GENCODE_URL, tmp)
        os.replace(tmp, path)
        downloaded = True
    mapping = parse_gtf_gene_chromosomes(path)
    info = {
        "source_url": GENCODE_URL,
        "release": f"GENCODE v{GENCODE_RELEASE} basic",
        "file": str(path.resolve()),
        "sha256": _sha256_file(path),
        "n_genes": len(mapping),
        "downloaded": downloaded,
    }
    return mapping, info


def rna_chromosomes(gene_ids: Iterable[str], mapping: dict[str, str]) -> np.ndarray:
    """Chromosome per gene (``"unmapped"`` when absent from the annotation)."""
    return np.array(
        [mapping.get(str(g).split(".")[0], UNMAPPED_CHROMOSOME) for g in gene_ids], dtype=object
    )


def atac_chromosomes(names: Iterable[str]) -> np.ndarray:
    """Chromosome per peak, parsed from ``chrN:start-end`` names.

    Raises
    ------
    ValueError
        On any name that does not parse, naming the first one.
    """
    names = [str(n) for n in names]
    out, bad = [], []
    for n in names:
        m = _PEAK_RE.match(n)
        if m is None:
            bad.append(n)
            out.append(None)
        else:
            out.append(m.group(1))
    if bad:
        raise ValueError(f"{len(bad)} peak name(s) do not parse as chrN:start-end, e.g. {bad[0]!r}")
    return np.array(out, dtype=object)


def _tmm_pair(
    obs: np.ndarray, ref: np.ndarray, n_o: float, n_r: float,
    *, logratio_trim: float, sum_trim: float,
) -> float:
    """edgeR ``.calcFactorTMM`` for one sample against the reference."""
    from scipy.stats import rankdata

    with np.errstate(divide="ignore", invalid="ignore"):
        log_r = np.log2((obs / n_o) / (ref / n_r))
        abs_e = (np.log2(obs / n_o) + np.log2(ref / n_r)) / 2
        v = (n_o - obs) / n_o / obs + (n_r - ref) / n_r / ref
    fin = np.isfinite(log_r) & np.isfinite(abs_e) & (abs_e > -1e10)
    log_r, abs_e, v = log_r[fin], abs_e[fin], v[fin]
    if log_r.size == 0:
        return 1.0
    if np.max(np.abs(log_r)) < 1e-6:
        return 1.0
    n = log_r.size
    lo_l = np.floor(n * logratio_trim) + 1
    hi_l = n + 1 - lo_l
    lo_s = np.floor(n * sum_trim) + 1
    hi_s = n + 1 - lo_s
    r_l, r_a = rankdata(log_r), rankdata(abs_e)
    keep = (r_l >= lo_l) & (r_l <= hi_l) & (r_a >= lo_s) & (r_a <= hi_s)
    with np.errstate(divide="ignore", invalid="ignore"):
        w = 1 / v[keep]
        num, den = log_r[keep] * w, w
        ok = np.isfinite(num) & np.isfinite(den)
        f = np.sum(num[ok]) / np.sum(den[ok])
    if not np.isfinite(f):
        f = 0.0
    return float(2.0 ** f)


def tmm_factors(
    counts: np.ndarray,
    *,
    lib_size: np.ndarray | None = None,
    logratio_trim: float = TMM_LOGRATIO_TRIM,
    sum_trim: float = TMM_SUM_TRIM,
) -> np.ndarray:
    """TMM normalization factors, following edgeR ``calcNormFactors(method = "TMM")``.

    Parameters
    ----------
    counts:
        Samples × features count matrix.
    lib_size:
        Library size per sample (default: row sums of ``counts``).
    logratio_trim, sum_trim:
        Trimming fractions for the M and A values (edgeR defaults 0.3, 0.05).

    Returns
    -------
    numpy.ndarray
        One factor per sample, rescaled to a geometric mean of 1.

    Notes
    -----
    The reference sample is the one whose upper-quartile factor is closest to
    the mean; M values are precision-weighted (delta method). Features with
    zero counts in every sample are dropped first, as in edgeR.
    """
    x = np.asarray(counts, dtype=float)
    lib = x.sum(axis=1) if lib_size is None else np.asarray(lib_size, dtype=float)
    x = x[:, (x > 0).any(axis=0)]
    n = x.shape[0]
    if x.shape[1] == 0 or n == 1:
        return np.ones(n)
    f75 = np.quantile(x, 0.75, axis=1) / lib
    if np.median(f75) < 1e-20:
        ref = int(np.argmax(np.sqrt(x).sum(axis=1)))
    else:
        ref = int(np.argmin(np.abs(f75 - f75.mean())))
    f = np.array([
        _tmm_pair(x[i], x[ref], lib[i], lib[ref], logratio_trim=logratio_trim, sum_trim=sum_trim)
        for i in range(n)
    ])
    return f / np.exp(np.mean(np.log(f)))


def log_cpm(counts: np.ndarray, lib_size: np.ndarray, *, prior_count: float = EXPORT_PRIOR_COUNT) -> np.ndarray:
    """log2-CPM with edgeR ``cpm(y, log = TRUE, prior.count)`` semantics.

    Parameters
    ----------
    counts:
        Samples × features count matrix.
    lib_size:
        Effective library size per sample (total count × TMM factor).
    prior_count:
        Average prior count; each sample gets ``prior_count · L / mean(L)``.

    Returns
    -------
    numpy.ndarray
        ``log2((count + p) / (L + 2p) · 1e6)``, samples × features.
    """
    x = np.asarray(counts, dtype=float)
    lib = np.asarray(lib_size, dtype=float)
    prior = prior_count * lib / lib.mean()
    adj = np.log(lib + 2 * prior)
    return (np.log(x + prior[:, None]) - adj[:, None] + np.log(1e6)) / np.log(2)


def smallest_cell(table: pd.DataFrame) -> tuple[int, str, dict]:
    """Size and name of the smallest sex × stage cell of a cohort's donors."""
    cells = {
        f"{sex}/stage {st}": int(((table.sex == sex) & (table.stage == st)).sum())
        for sex in sorted(table.sex.unique()) for st in STAGES
    }
    name = min(cells, key=lambda c: (cells[c], c))
    return cells[name], name, cells


def _excluded_features(
    names, chromosomes, *, excluded_chromosomes: Iterable[str], excluded_gene_sets: Iterable[str] = ()
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Step 2 of the export: features on excluded chromosomes, then features of excluded gene sets.

    Parameters
    ----------
    names, chromosomes:
        Feature names and chromosome per feature.
    excluded_chromosomes:
        Chromosomes whose features are removed (status ``sex_chromosome``).
    excluded_gene_sets:
        Names of ``TECHNICAL_GENE_SETS`` whose matching features are removed
        next (status ``excluded_gene_set``); a matched feature already on an
        excluded chromosome keeps ``sex_chromosome``. ATAC names match no set.

    Returns
    -------
    (on_chromosome, in_gene_set, matched)
        Boolean masks over the input features (``in_gene_set`` excludes the
        features already removed by chromosome) and, per excluded set, the
        matched names, the names removed at this step and the matched names
        already on an excluded chromosome.
    """
    names = np.asarray(names, dtype=str)
    chromosomes = np.asarray(chromosomes, dtype=object)
    on_chrom = np.isin(chromosomes, tuple(excluded_chromosomes))
    in_set = np.zeros(names.size, dtype=bool)
    matched: dict[str, dict] = {}
    for set_name in excluded_gene_sets:
        if set_name not in TECHNICAL_GENE_SETS:
            raise ValueError(f"unknown gene set {set_name!r}; TECHNICAL_GENE_SETS holds {list(TECHNICAL_GENE_SETS)}")
        compiled = [re.compile(p) for p in TECHNICAL_GENE_SETS[set_name]]
        hit = np.array([any(c.search(n) for c in compiled) for n in names], dtype=bool)
        matched[set_name] = {
            "n_matched": int(hit.sum()),
            "n_removed": int((hit & ~on_chrom).sum()),
            "matched": names[hit].tolist(),
            "removed": names[hit & ~on_chrom].tolist(),
            "on_excluded_chromosome": names[hit & on_chrom].tolist(),
        }
        in_set |= hit & ~on_chrom
    return on_chrom, in_set, matched


def trend_ratio(mean, variance, *, n_bins: int = TREND_N_BINS) -> tuple[np.ndarray, np.ndarray, list[dict]]:
    """Variance over a binned-median mean–variance trend (the ``trend_ratio`` selection rule).

    Parameters
    ----------
    mean, variance:
        Mean and population variance of log2-CPM per abundance-passing feature.
    n_bins:
        Number of equal-count bins of the mean (``TREND_N_BINS``).

    Returns
    -------
    (trend, ratio, bins)
        The trend at each feature's mean, ``variance / trend`` and the bins
        (``n``, ``lower``, ``upper``, ``median_mean``, ``median_variance``).

    Notes
    -----
    The features sorted by mean are cut with ``np.array_split`` into ``n_bins``
    bins (sizes differing by at most one); each bin's median variance is placed
    at its median mean and the trend is ``np.interp`` through these points,
    flat beyond the end bins. The ratio is undefined, and an error, when the
    bins' median means are not strictly increasing or a median variance is not
    positive.
    """
    mean = np.asarray(mean, dtype=float)
    variance = np.asarray(variance, dtype=float)
    if mean.size < n_bins:
        raise ValueError(f"trend_ratio needs at least {n_bins} features, got {mean.size}")
    order = np.argsort(mean, kind="stable")
    bins = []
    for chunk in np.array_split(order, n_bins):
        bins.append({
            "n": int(chunk.size),
            "lower": float(mean[chunk].min()),
            "upper": float(mean[chunk].max()),
            "median_mean": float(np.median(mean[chunk])),
            "median_variance": float(np.median(variance[chunk])),
        })
    xp = np.array([b["median_mean"] for b in bins])
    fp = np.array([b["median_variance"] for b in bins])
    if np.any(np.diff(xp) <= 0):
        raise ValueError("trend_ratio: the bins' median means are not strictly increasing (tied means); the trend is undefined")
    if np.any(fp <= 0):
        raise ValueError("trend_ratio: a bin's median variance is not positive; the trend is undefined")
    trend = np.interp(mean, xp, fp)
    return trend, variance / trend, bins


def residualize_values(values, design) -> tuple[np.ndarray, np.ndarray]:
    """Replace each column by its mean plus the residual of its OLS fit on an intercept and ``design``.

    Parameters
    ----------
    values:
        Donors × features (log2-CPM).
    design:
        Donors × covariates, aligned to the rows of ``values``; finite and of
        full column rank together with the intercept.

    Returns
    -------
    (residualized, r_squared)
        The residualized values (same means as the input, each column
        orthogonal to every design column and the intercept) and the per-feature
        R² of the fit (NaN for a constant column).
    """
    values = np.asarray(values, dtype=float)
    design = np.asarray(design, dtype=float)
    if design.ndim == 1:
        design = design[:, None]
    if design.shape[0] != values.shape[0]:
        raise ValueError(f"residualization design has {design.shape[0]} rows for {values.shape[0]} donors")
    if not np.all(np.isfinite(design)):
        raise ValueError("residualization design contains non-finite values")
    X = np.column_stack([np.ones(values.shape[0]), design])
    if np.linalg.matrix_rank(X) < X.shape[1]:
        raise ValueError("residualization design is rank deficient (a constant or collinear covariate)")
    beta = np.linalg.lstsq(X, values, rcond=None)[0]
    resid = values - X @ beta
    mean = values.mean(axis=0)
    ss_tot = ((values - mean) ** 2).sum(axis=0)
    ss_res = (resid ** 2).sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        r2 = np.where(ss_tot > 0, 1.0 - ss_res / ss_tot, np.nan)
    return resid + mean, r2


def loess_span(n_features: int, *, rule: dict | None = None) -> float:
    """limma's adaptive lowess span, ``chooseLowessSpan(n, small.n, min.span, power)``.

    Parameters
    ----------
    n_features:
        Number of features the curve is fitted over.
    rule:
        ``{"small_n", "min_span", "power"}``; ``LOESS_SPAN_RULE`` (limma's
        defaults 50, 0.3, 1/3) by default.

    Returns
    -------
    float
        ``min(min_span + (1 − min_span) · (small_n / n) ** power, 1)``.
    """
    rule = LOESS_SPAN_RULE if rule is None else rule
    if n_features < 1:
        raise ValueError(f"loess_span needs at least one feature, got {n_features}")
    span = rule["min_span"] + (1 - rule["min_span"]) * (rule["small_n"] / n_features) ** rule["power"]
    return float(min(span, 1.0))


def tilt_statistics(values, *, n_bins: int = NORMALIZATION_TILT_BINS) -> dict[str, np.ndarray]:
    """Per-donor abundance tilt of a donors × features block of log2-CPM values.

    Parameters
    ----------
    values:
        Donors × features.
    n_bins:
        Number of equal-count bins of A (``NORMALIZATION_TILT_BINS``).

    Returns
    -------
    dict
        Per donor (arrays in row order): ``slope``, the OLS slope (with an
        intercept) of M on A over all features; ``low_bin`` and ``high_bin``,
        the median M over the lowest and the highest of ``n_bins`` equal-count
        bins of A (``np.array_split`` of the features sorted by A, stable).
        A is the per-feature mean over donors and M the donor's value minus A,
        both of the values given. Label-free.
    """
    x = np.asarray(values, dtype=float)
    if x.shape[1] < n_bins:
        raise ValueError(f"tilt_statistics needs at least {n_bins} features, got {x.shape[1]}")
    a = x.mean(axis=0)
    m = x - a
    ac = a - a.mean()
    slope = (m * ac).sum(axis=1) / (ac ** 2).sum()
    bins = np.array_split(np.argsort(a, kind="stable"), n_bins)
    return {
        "slope": slope,
        "low_bin": np.median(m[:, bins[0]], axis=1),
        "high_bin": np.median(m[:, bins[-1]], axis=1),
    }


def cyclic_loess(
    values,
    *,
    iterations: int = LOESS_ITERATIONS,
    span: float | None = None,
) -> dict:
    """limma's fast cyclic loess (``normalizeCyclicLoess(method = "fast")``) of a donors × features block.

    Parameters
    ----------
    values:
        Donors × features (log2-CPM over the abundance-passing features).
    iterations:
        Cyclic iterations (``LOESS_ITERATIONS``); each recomputes the reference.
    span:
        Lowess span; ``None`` for limma's adaptive span from the number of
        features (``loess_span``).

    Returns
    -------
    dict
        ``values`` (the normalized block, same shape), ``span``,
        ``iterations``, ``robustness_iterations``, ``delta`` (per iteration,
        ``LOESS_DELTA_FRACTION`` of the range of A), ``n_features``,
        ``tilt_before`` / ``tilt_after`` (``tilt_statistics`` of the input and
        the normalized values) and ``max_abs_correction`` (per donor, the
        largest absolute change of a value).

    Notes
    -----
    At each iteration A is the per-feature mean over donors of the current
    values; for each donor the lowess fit of M = value − A on A
    (``statsmodels.nonparametric.smoothers_lowess.lowess`` with ``frac =
    span``, ``it = LOESS_ROBUSTNESS_ITERATIONS``, ``delta``) is subtracted. As
    in limma, A is fixed within an iteration, so the donor order does not
    matter. Reads nothing but the values.
    """
    from statsmodels.nonparametric.smoothers_lowess import lowess

    original = np.asarray(values, dtype=float)
    if original.ndim != 2:
        raise ValueError("cyclic_loess expects a donors × features matrix")
    if not np.all(np.isfinite(original)):
        raise ValueError("cyclic_loess: the block contains non-finite values")
    n_donors, n_feat = original.shape
    span = loess_span(n_feat) if span is None else float(span)
    x = original.copy()
    deltas = []
    for _ in range(iterations):
        a = x.mean(axis=0)
        delta = LOESS_DELTA_FRACTION * float(a.max() - a.min())
        deltas.append(delta)
        for i in range(n_donors):
            fit = lowess(
                x[i] - a, a, frac=span, it=LOESS_ROBUSTNESS_ITERATIONS, delta=delta, return_sorted=False
            )
            x[i] = x[i] - fit
    if not np.all(np.isfinite(x)):
        raise ValueError("cyclic_loess: the normalized block contains non-finite values")
    return {
        "values": x,
        "span": span,
        "iterations": int(iterations),
        "robustness_iterations": LOESS_ROBUSTNESS_ITERATIONS,
        "delta": deltas,
        "n_features": int(n_feat),
        "tilt_before": tilt_statistics(original),
        "tilt_after": tilt_statistics(x),
        "max_abs_correction": np.abs(x - original).max(axis=1),
    }


def preprocess_block(
    counts: np.ndarray,
    names: np.ndarray,
    chromosomes: np.ndarray,
    *,
    k: int,
    min_cpm: float = EXPORT_MIN_CPM,
    prior_count: float = EXPORT_PRIOR_COUNT,
    top_n: int = EXPORT_TOP_N,
    excluded_chromosomes: Iterable[str] = EXPORT_EXCLUDED_CHROMOSOMES,
    excluded_gene_sets: Iterable[str] = (),
    selection_rule: str = "variance",
    normalization: str = "none",
    residualize=None,
) -> dict:
    """Preprocess one omic block of one cohort (steps 2–6 of the export).

    Parameters
    ----------
    counts:
        Cohort donors × features raw counts (step 1 already applied).
    names, chromosomes:
        Feature names and chromosome per feature.
    k:
        Minimum number of donors with CPM ≥ ``min_cpm`` (smallest sex × stage cell).
    min_cpm, prior_count, top_n, excluded_chromosomes:
        Thresholds; defaults are the module constants.
    excluded_gene_sets:
        Names of ``TECHNICAL_GENE_SETS`` whose matching features are removed at
        step 2 after the chromosome removal (status ``excluded_gene_set``; a
        matched feature already on an excluded chromosome keeps
        ``sex_chromosome``). Empty for the primary export; ATAC names match no set.
    selection_rule:
        One of ``SELECTION_RULES``: ``variance`` (the primary's; dropped
        features carry status ``variance``), ``trend_ratio`` or ``mean``
        (dropped features carry ``selection``).
    normalization:
        One of ``NORMALIZATION_RULES``: ``none`` (the primary's; the values
        are left unchanged) or ``cyclic_loess``: after step 5 and before the
        residualization, the abundance-passing features' log2-CPM is
        normalized by limma's fast cyclic loess (``cyclic_loess``, adaptive
        span), and the residualization and step 6 see the normalized values.
    residualize:
        ``None`` (the primary's) or a donors × covariates design aligned to the
        rows of ``counts``: after step 5 (and the normalization) every
        abundance-passing feature is replaced by its mean plus the residual of
        its OLS fit on an intercept and the design (``residualize_values``),
        so units stay log2-CPM, and step 6 ranks the residualized values.

    Returns
    -------
    dict
        ``values`` (donors × kept features, log2-CPM, input feature order),
        ``kept`` (names), ``status``, ``variance`` (of the values step 6 ranks:
        the normalized and residualized values when those steps apply),
        ``mean`` and ``selection_score`` (equal to ``variance`` under
        ``variance``) per input feature, ``trend`` and ``ratio``
        (``trend_ratio`` only, else ``None``), ``trend_bins`` (``trend_ratio``
        only), ``r_squared`` (residualizing only, else ``None``),
        ``removed_by_gene_set``, ``selection_rule``, ``normalization`` (the
        rule), ``normalization_fit`` (``cyclic_loess`` output without the
        values, else ``None``), ``steps`` (feature counts after each step;
        ``after_gene_set`` only when sets were excluded),
        ``removed_by_chromosome``, ``lib_size_filter``, ``lib_size``,
        ``tmm_factors``, ``mean_log_cpm_range`` of the kept features,
        ``findings``, and the filtered ``filtered_counts`` /
        ``filtered_log_cpm`` (the step-5 log2-CPM, before normalization and
        residualization) / ``filtered_normalized`` (the normalized log2-CPM
        before residualization; ``None`` under ``none``) / ``filtered_names``.

    Notes
    -----
    Labels are never used: ``k`` is the only cohort-derived input. With the
    default keywords the computation is the primary export's, unchanged.
    """
    counts = np.asarray(counts, dtype=float)
    names = np.asarray(names, dtype=str)
    chromosomes = np.asarray(chromosomes, dtype=object)
    if selection_rule not in SELECTION_RULES:
        raise ValueError(f"unknown selection rule {selection_rule!r}; expected one of {SELECTION_RULES}")
    if normalization not in NORMALIZATION_RULES:
        raise ValueError(f"unknown normalization {normalization!r}; expected one of {NORMALIZATION_RULES}")
    excluded_gene_sets = tuple(excluded_gene_sets)
    n_feat = counts.shape[1]
    excluded = tuple(excluded_chromosomes)
    status = np.full(n_feat, "kept", dtype=object)
    variance = np.full(n_feat, np.nan)
    mean = np.full(n_feat, np.nan)
    score = np.full(n_feat, np.nan)
    trend = ratio = r_squared = trend_bins = None
    findings: list[str] = []

    # (2) excluded chromosomes, then excluded gene sets
    on_sex, in_set, matched = _excluded_features(
        names, chromosomes, excluded_chromosomes=excluded, excluded_gene_sets=excluded_gene_sets
    )
    status[on_sex] = "sex_chromosome"
    status[in_set] = "excluded_gene_set"
    removed_by_chrom = {c: int((chromosomes == c).sum()) for c in excluded}
    idx2 = np.flatnonzero(~on_sex & ~in_set)

    # (3) abundance: CPM over the post-step-2 library sizes
    c2 = counts[:, idx2]
    lib_filter = c2.sum(axis=1)
    cpm = c2 / lib_filter[:, None] * 1e6
    passes = (cpm >= min_cpm).sum(axis=0) >= k
    status[idx2[~passes]] = "abundance"
    idx3 = idx2[passes]

    # (4) TMM on the filtered counts, (5) log2-CPM
    c3 = counts[:, idx3]
    lib = c3.sum(axis=1)
    factors = tmm_factors(c3, lib_size=lib)
    lcpm = log_cpm(c3, lib * factors, prior_count=prior_count)
    lcpm_step5 = lcpm

    # (5a) trend normalization (cyclic loess against the per-feature mean)
    normalization_fit = normalized = None
    if normalization == "cyclic_loess":
        normalization_fit = cyclic_loess(lcpm)
        lcpm = normalized = normalization_fit.pop("values")

    # (5b) residualization on an intercept and the design; means kept
    if residualize is not None:
        lcpm, r2 = residualize_values(lcpm, residualize)
        r_squared = np.full(n_feat, np.nan)
        r_squared[idx3] = r2

    # (6) top-N by the selection rule, ties by name; output keeps input order
    var3 = lcpm.var(axis=0, ddof=0)
    mean3 = lcpm.mean(axis=0)
    variance[idx3] = var3
    mean[idx3] = mean3
    if selection_rule == "variance":
        score3, dropped_status = var3, "variance"
    elif selection_rule == "mean":
        score3, dropped_status = mean3, "selection"
    else:
        trend3, ratio3, trend_bins = trend_ratio(mean3, var3)
        trend, ratio = np.full(n_feat, np.nan), np.full(n_feat, np.nan)
        trend[idx3], ratio[idx3] = trend3, ratio3
        score3, dropped_status = ratio3, "selection"
    score[idx3] = score3
    order = np.lexsort((names[idx3], -score3))
    if idx3.size < top_n:
        findings.append(f"only {idx3.size} features pass the abundance filter (< N = {top_n}); all are kept")
    sel = np.sort(order[:top_n])
    dropped = np.setdiff1d(np.arange(idx3.size), sel)
    status[idx3[dropped]] = dropped_status
    values = lcpm[:, sel]
    if not np.all(np.isfinite(values)):
        raise ValueError("log2-CPM contains non-finite values")

    steps = {"input": int(n_feat), "after_sex_chromosome": int((~on_sex).sum())}
    if excluded_gene_sets:
        steps["after_gene_set"] = int(idx2.size)
    steps["after_abundance"] = int(idx3.size)
    steps["after_variance"] = int(sel.size)
    mean_lcpm = values.mean(axis=0) if values.size else np.array([np.nan])
    return {
        "values": values,
        "kept": names[idx3[sel]],
        "status": status,
        "variance": variance,
        "mean": mean,
        "selection_score": score,
        "trend": trend,
        "ratio": ratio,
        "trend_bins": trend_bins,
        "r_squared": r_squared,
        "removed_by_gene_set": matched,
        "selection_rule": selection_rule,
        "normalization": normalization,
        "normalization_fit": normalization_fit,
        "steps": steps,
        "removed_by_chromosome": removed_by_chrom,
        "lib_size_filter": lib_filter,
        "lib_size": lib,
        "tmm_factors": factors,
        "mean_log_cpm_range": [float(np.min(mean_lcpm)), float(np.max(mean_lcpm))],
        "findings": findings,
        "filtered_counts": c3,
        "filtered_log_cpm": lcpm_step5,
        "filtered_normalized": normalized,
        "filtered_names": names[idx3],
    }


_EXPORT_INPUTS = {
    "rna_counts.h5ad": "pseudobulk",
    "atac_counts.h5ad": "pseudobulk",
    "donor_table.csv": "pseudobulk",
    "manifest.json": "pseudobulk",
    "donor_covariates.csv": "covariates",
}


def check_export_inputs(
    pseudobulk_dir: str | Path, cohorts: Iterable[str], *, cell_type: str = DEFAULT_CELL_TYPE
) -> dict:
    """Validate the pseudobulk directory before exporting; return its manifest.

    Raises
    ------
    FileNotFoundError
        A required file is missing (naming the command that produces it).
    ValueError
        The pseudobulk run is partial, its manifest records another cell
        type than ``cell_type`` (:func:`check_manifest_cell_type`), or a
        cohort is unknown or not analysable.
    """
    pb = Path(pseudobulk_dir)
    for name, command in _EXPORT_INPUTS.items():
        if not (pb / name).exists():
            raise FileNotFoundError(
                f"{pb / name} is missing; run `sead-sexdiff {command}` (and sync) to produce it"
            )
    manifest = json.loads((pb / "manifest.json").read_text())
    check_manifest_cell_type(manifest, pb / "manifest.json", cell_type)
    if manifest.get("partial", True):
        raise ValueError(f"{pb / 'manifest.json'}: the pseudobulk run is partial; export needs a full run")
    for cohort in cohorts:
        info = manifest.get("cohorts", {}).get(cohort)
        if info is None:
            raise ValueError(f"{pb / 'manifest.json'}: no cohort {cohort!r}")
        if not info.get("analysable", False):
            raise ValueError(
                f"Cohort {cohort} is marked not analysable in {pb / 'manifest.json'} "
                f"(empty cells: {info.get('empty_cells')})"
            )
    return manifest


def _features_frame(names, chromosomes, res: dict, gene_ids=None) -> pd.DataFrame:
    df = pd.DataFrame({"feature": names})
    if gene_ids is not None:
        df["gene_id"] = gene_ids
    df["chromosome"] = chromosomes
    df["status"] = res["status"]
    df["variance"] = res["variance"]
    if res.get("selection_rule", "variance") != "variance" or res.get("r_squared") is not None:
        df["mean"] = res["mean"]
    if res.get("trend") is not None:
        df["trend"] = res["trend"]
        df["ratio"] = res["ratio"]
    if res.get("r_squared") is not None:
        df["r_squared"] = res["r_squared"]
    return df


def _load_export_inputs(pseudobulk_dir: str | Path, data_root: str | Path | None, *, log_prefix: str = "export") -> dict:
    """Load everything the per-cohort export reads: annotation, tables, matrices.

    Parameters
    ----------
    pseudobulk_dir:
        Output directory of ``pseudobulk`` and ``covariates`` (already checked
        with :func:`check_export_inputs`).
    data_root:
        Data root holding (or receiving) the gene annotation.
    log_prefix:
        Prefix of the progress lines.

    Returns
    -------
    dict
        ``table`` (``donor_table.csv``), ``covariates`` (``donor_covariates.csv``),
        ``blocks`` (omic → ``{"X", "donors", "names", "gene_ids", "chrom"}``),
        ``annotation`` (provenance plus ``n_rna_genes`` / ``n_rna_unmapped``)
        and ``pb_info`` (the export manifest's ``input`` block: pseudobulk
        directory, manifest hashes and timestamp).
    """
    import scipy.sparse as sp

    pb = Path(pseudobulk_dir)
    pb_manifest = json.loads((pb / "manifest.json").read_text())
    mapping, annotation = resolve_gene_annotation(data_root)
    print(f"[{log_prefix}] Annotation: {annotation['release']} ({annotation['n_genes']} genes, sha256 {annotation['sha256'][:12]}…)")

    table = pd.read_csv(pb / "donor_table.csv", dtype={"donor_id": str, "stage": str, "sex": str})
    covariates = pd.read_csv(pb / "donor_covariates.csv", dtype={"donor_id": str, "education_level": str})
    blocks = {}
    for o in ("rna", "atac"):
        a = ad.read_h5ad(pb / f"{o}_counts.h5ad")
        X = sp.csr_matrix(a.X)
        if o == "rna":
            gene_ids = a.var["gene_ids"].astype(str).to_numpy()
            chrom = rna_chromosomes(gene_ids, mapping)
        else:
            gene_ids = None
            chrom = atac_chromosomes(a.var_names)
        blocks[o] = {
            "X": X, "donors": a.obs_names.astype(str).to_numpy(),
            "names": a.var_names.astype(str).to_numpy(), "gene_ids": gene_ids, "chrom": chrom,
        }
        del a
    annotation["n_rna_genes"] = int(blocks["rna"]["names"].size)
    annotation["n_rna_unmapped"] = int((blocks["rna"]["chrom"] == UNMAPPED_CHROMOSOME).sum())
    cov_manifest = pb / "covariates_manifest.json"
    pb_info = {
        "pseudobulk_dir": str(pb.resolve()),
        "pseudobulk_manifest_sha256": _sha256_file(pb / "manifest.json"),
        "pseudobulk_manifest_created_at": pb_manifest.get("created_at"),
        "covariates_manifest_sha256": _sha256_file(cov_manifest) if cov_manifest.exists() else None,
    }
    return {"table": table, "covariates": covariates, "blocks": blocks, "annotation": annotation, "pb_info": pb_info}


def _cohort_rows(table: pd.DataFrame, covariates: pd.DataFrame, cohort: str, pseudobulk_dir: Path) -> pd.DataFrame:
    """The cohort's rows of the donor table, sorted by ``donor_id`` and checked against the covariates."""
    sub = table[table[cohort].astype(bool)].sort_values("donor_id").reset_index(drop=True)
    missing_cov = sorted(set(sub["donor_id"]) - set(covariates.donor_id))
    if missing_cov:
        raise ValueError(f"{Path(pseudobulk_dir) / 'donor_covariates.csv'}: no row for cohort {cohort} donor(s) {missing_cov}")
    return sub


def _residualization_design(
    donors: np.ndarray,
    counts: dict[str, np.ndarray],
    blocks: dict,
    technical: pd.DataFrame,
    *,
    excluded_chromosomes: Iterable[str],
    excluded_gene_sets: Iterable[str],
    pseudobulk_dir: str | Path | None = None,
) -> pd.DataFrame:
    """The residualization design of one cohort: ``RESIDUALIZATION_COVARIATES`` per donor, in export order.

    Parameters
    ----------
    donors:
        The cohort's donors in export order.
    counts:
        Omic → cohort donors × features raw counts (the rows ``_export_cohort``
        passes to ``preprocess_block``).
    blocks:
        Omic → ``{"names", "chrom", ...}`` (``_load_export_inputs``).
    technical:
        ``donor_technical.csv`` indexed by ``donor_id`` (or with that column).
    excluded_chromosomes, excluded_gene_sets:
        The export's step-2 settings; ``log10_lib_<omic>`` is the log10 of the
        donor's counts over the features kept at step 2 (the step-3 library
        size, ``lib_size_filter``).
    pseudobulk_dir:
        Named in error messages.

    Returns
    -------
    pd.DataFrame
        Donors × ``RESIDUALIZATION_COVARIATES`` (exactly these columns, in
        this order), indexed by ``donor_id``. No sex, stage or analysis
        covariate.

    Raises
    ------
    ValueError
        A technical column is absent, or missing for a cohort donor (the
        message names the donor and the column); a library size is not
        positive.
    """
    tech = technical.set_index("donor_id") if "donor_id" in technical.columns else technical
    label = f"{Path(pseudobulk_dir) / TECHNICAL_OUTPUTS['table']}" if pseudobulk_dir is not None else TECHNICAL_OUTPUTS["table"]
    design = pd.DataFrame(index=pd.Index(donors, name="donor_id"))
    for col in RESIDUALIZATION_COVARIATES:
        if col.startswith("log10_lib_"):
            o = col[len("log10_lib_"):]
            on_chrom, in_set, _ = _excluded_features(
                blocks[o]["names"], blocks[o]["chrom"],
                excluded_chromosomes=excluded_chromosomes, excluded_gene_sets=excluded_gene_sets,
            )
            lib = np.asarray(counts[o], dtype=float)[:, ~on_chrom & ~in_set].sum(axis=1)
            if np.any(lib <= 0):
                raise ValueError(f"{o}: a cohort donor has no counts on the features kept at step 2; log10 library size undefined")
            design[col] = np.log10(lib)
        else:
            if col not in tech.columns:
                raise ValueError(f"{label} lacks column {col!r} (RESIDUALIZATION_COVARIATES)")
            absent = sorted(set(donors) - set(tech.index))
            if absent:
                raise ValueError(f"{label} has no row for cohort donor(s) {absent}")
            values = pd.to_numeric(tech.loc[donors, col], errors="coerce").to_numpy(dtype=float)
            bad = [d for d, v in zip(donors, values) if not np.isfinite(v)]
            if bad:
                raise ValueError(f"{label}: column {col!r} is missing for cohort donor(s) {bad}; the residualized variant needs every covariate for every cohort donor")
            design[col] = values
    return design


def _quantiles(values) -> dict[str, float]:
    """``DIAGNOSTIC_QUANTILES`` of the values, keyed ``q05`` … ``q95``."""
    values = np.asarray(values, dtype=float)
    return {f"q{int(round(q * 100)):02d}": float(np.quantile(values, q)) for q in DIAGNOSTIC_QUANTILES}


def _normalization_trend_table(donors: np.ndarray, fits: dict[str, dict]) -> pd.DataFrame:
    """``normalization_trend.csv``: per normalized omic and donor, the tilt statistics before and after and the largest correction."""
    frames = []
    for o, fit in fits.items():
        frame = pd.DataFrame({"donor_id": donors, "omic": o})
        for stat in NORMALIZATION_TILT_STATISTICS:
            frame[f"{stat}_before"] = fit["tilt_before"][stat]
            frame[f"{stat}_after"] = fit["tilt_after"][stat]
        frame["max_abs_correction"] = fit["max_abs_correction"]
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def _normalization_fit_block(fits: dict[str, dict], *, trend_sha256: str | None) -> dict:
    """The export manifest's ``normalization_fit`` block for the normalized omics."""
    from importlib.metadata import PackageNotFoundError, version

    try:
        statsmodels_version = version("statsmodels")
    except PackageNotFoundError:
        statsmodels_version = None
    block: dict = {}
    for o, fit in fits.items():
        block[o] = {
            "rule": "cyclic_loess",
            "reference": "feature_mean",
            "n_features": fit["n_features"],
            "span": fit["span"],
            "adaptive_span": True,
            "iterations": fit["iterations"],
            "robustness_iterations": fit["robustness_iterations"],
            "delta": fit["delta"],
            "tilt_quantiles": {
                stat: {"before": _quantiles(fit["tilt_before"][stat]), "after": _quantiles(fit["tilt_after"][stat])}
                for stat in NORMALIZATION_TILT_STATISTICS
            },
            "max_abs_correction": {
                "max": float(np.max(fit["max_abs_correction"])),
                "quantiles": _quantiles(fit["max_abs_correction"]),
            },
        }
    block.update({
        "method": "limma normalizeCyclicLoess(method = \"fast\"): at each iteration the reference is the per-feature "
                  "mean over donors (A); each donor's deviation from it (M) is fitted by lowess on A and subtracted",
        "lowess": f"statsmodels.nonparametric.smoothers_lowess.lowess (statsmodels {statsmodels_version}); frac = span, "
                  f"it = {LOESS_ROBUSTNESS_ITERATIONS}, delta = {LOESS_DELTA_FRACTION} x range(A), as R's lowess in limma's loessFit",
        "span_rule": {**LOESS_SPAN_RULE, "formula": "min(min_span + (1 - min_span) * (small_n / n) ** power, 1), "
                      "n = number of abundance-passing features (limma chooseLowessSpan)"},
        "delta_fraction": LOESS_DELTA_FRACTION,
        "position": "after log2-CPM (step 5), before the residualization (when declared) and the selection step (step 6)",
        "fitting_set": "all abundance-passing features of the block (the features TMM is estimated on)",
        "label_free": "reads only the block's log2-CPM; no sex, stage, analysis covariate or technical variable",
        "tilt": {
            "n_bins": NORMALIZATION_TILT_BINS,
            "definition": "per donor, A = per-feature mean over donors and M = value - A of the values before or after "
                          "normalization; slope = OLS slope of M on A; low_bin / high_bin = median M over the lowest / "
                          f"highest of {NORMALIZATION_TILT_BINS} equal-count bins of A",
            "quantiles_over": "donors",
        },
        "trend_file": NORMALIZATION_TREND_FILE,
        "trend_sha256": trend_sha256,
    })
    return block


def _per_omic_setting(value, *, default: str, name: str, label: str, allowed: Iterable[str]) -> dict[str, str]:
    """Resolve an export setting given per omic, as one value for both omics, or ``None`` (the default for both)."""
    if value is None:
        out = {o: default for o in ("rna", "atac")}
    elif isinstance(value, str):
        out = {o: value for o in ("rna", "atac")}
    else:
        out = {o: default for o in ("rna", "atac")}
        unknown = sorted(set(value) - set(out))
        if unknown:
            raise ValueError(f"{name} names unknown omic(s) {unknown}; expected rna and/or atac")
        out.update(value)
    allowed = tuple(allowed)
    for o, rule in out.items():
        if rule not in allowed:
            raise ValueError(f"unknown {label} {rule!r} for {o}; expected one of {allowed}")
    return out


def _export_cohort(
    cohort: str,
    sub: pd.DataFrame,
    blocks: dict,
    covariates: pd.DataFrame,
    annotation: dict,
    pb_info: dict,
    *,
    out_dir: Path,
    final_name: str,
    excluded_chromosomes: Iterable[str],
    excluded_gene_sets: Iterable[str] = (),
    selection_rule: dict[str, str] | str | None = None,
    normalization: dict[str, str] | str | None = None,
    residualize: pd.DataFrame | None = None,
    variant: dict | None = None,
    debug_dir: Path | None = None,
    log_prefix: str = "export",
    cell_type: str = DEFAULT_CELL_TYPE,
    input_findings: Iterable[str] = (),
) -> dict:
    """Preprocess one cohort and write its complete MOTCO input directory.

    Shared by ``run_export`` (the primary export, ``EXPORT_EXCLUDED_CHROMOSOMES``)
    and ``run_preprocessing_sensitivity`` (a registered variant): the two differ
    only in the keyword settings, so a variant export is the primary export
    with one setting changed and nothing else.

    Parameters
    ----------
    cohort:
        Cohort name (manifest and progress lines).
    sub:
        The cohort's rows of ``donor_table.csv``, sorted by ``donor_id``
        (``_cohort_rows``).
    blocks:
        Omic → ``{"X", "donors", "names", "gene_ids", "chrom"}``
        (``_load_export_inputs``).
    covariates:
        ``donor_covariates.csv``.
    annotation:
        Provenance of the gene annotation, recorded in the manifest.
    pb_info:
        The manifest's ``input`` block (pseudobulk directory and manifest hashes).
    out_dir, final_name:
        The directory is written to ``<out_dir>/<final_name>.tmp/`` and renamed
        to ``<out_dir>/<final_name>`` when complete.
    excluded_chromosomes:
        Chromosomes whose features are removed before the abundance filter
        (recorded in ``thresholds.excluded_chromosomes``).
    excluded_gene_sets:
        Names of ``TECHNICAL_GENE_SETS`` removed from the RNA block at step 2
        after the chromosome removal (``preprocess_block``); every set must
        match at least one RNA gene. Empty for the primary export (recorded as
        ``excluded_gene_sets: []``; ``removed_by_gene_set`` is written only when
        sets were excluded).
    selection_rule:
        The step-6 rule per omic (``{"rna": ..., "atac": ...}``), one rule for
        both, or ``None`` for the primary's ``variance`` (recorded as
        ``selection_rule: {rna: variance, atac: variance}``; ``trend_bins`` is
        added for an omic ranked by ``trend_ratio``).
    normalization:
        The trend normalization per omic (``{"rna": ..., "atac": ...}`` of
        ``NORMALIZATION_RULES``), one rule for both, or ``None`` for the
        primary's ``none`` (recorded as ``normalization: {rna: none, atac:
        none}``). When an omic is normalized the manifest gains a
        ``normalization_fit`` block (span, iterations, reference, tilt
        quantiles per normalized omic), ``preprocessing_order`` gains the
        step after log2-CPM, and ``normalization_trend.csv`` (per donor and
        normalized omic, the tilt statistics before and after and the largest
        absolute correction) is written beside the blocks.
    residualize:
        ``None`` (the primary's; recorded as ``residualization: null``) or
        ``donor_technical.csv`` (indexed by ``donor_id`` or with that column):
        the design ``RESIDUALIZATION_COVARIATES`` is built per cohort
        (``_residualization_design``), written as ``residualization.csv`` and
        passed to ``preprocess_block`` for both omics.
    variant:
        For a preprocessing-sensitivity variant, ``{"name", "description",
        "settings", "reference_export_dir"}``, written as the manifest's
        ``variant`` block; ``None`` for the primary export, whose manifest then
        has no such key.
    debug_dir:
        If given, the filtered counts, TMM factors and log2-CPM per omic go
        there, plus the normalized log2-CPM and the loess settings of a
        normalized omic (``_write_debug``).
    log_prefix:
        Prefix of the progress lines, ``[export]`` by default.
    cell_type:
        The cell type of the pseudobulk (checked by the caller), recorded as
        the manifest's ``cell_type``.
    input_findings:
        Findings from the caller's input checks (such as a pseudobulk
        manifest without a ``cell_type``), listed first.

    Returns
    -------
    dict
        ``{"dir": Path, "manifest": dict, "findings": [...], "k": int, "cells": dict}``.
    """
    import shutil

    excluded = tuple(excluded_chromosomes)
    gene_sets = tuple(excluded_gene_sets)
    for set_name in gene_sets:
        if set_name not in TECHNICAL_GENE_SETS:
            raise ValueError(f"unknown gene set {set_name!r}; TECHNICAL_GENE_SETS holds {list(TECHNICAL_GENE_SETS)}")
    rules = _per_omic_setting(selection_rule, default="variance", name="selection_rule", label="selection rule", allowed=SELECTION_RULES)
    norms = _per_omic_setting(normalization, default="none", name="normalization", label="normalization", allowed=NORMALIZATION_RULES)
    normalized_omics = [o for o in ("rna", "atac") if norms[o] != "none"]
    donors = sub["donor_id"].to_numpy()
    k, k_cell, cells = smallest_cell(sub)
    print(f"[{log_prefix}] {cohort}: {len(donors)} donors, k = {k} (smallest cell {k_cell})")
    findings = list(input_findings)
    n_unmapped = annotation.get("n_rna_unmapped", 0)
    if n_unmapped:
        findings.append(
            f"rna: {n_unmapped} of {blocks['rna']['names'].size} genes are absent from "
            f"{annotation['release']}; they are kept (chromosome '{UNMAPPED_CHROMOSOME}')"
        )
    counts_by_omic = {}
    for o, blk in blocks.items():
        pos = pd.Index(blk["donors"]).get_indexer(donors)
        if (pos < 0).any():
            raise ValueError(f"{o}_counts.h5ad lacks cohort {cohort} donor(s) {list(donors[pos < 0])}")
        counts_by_omic[o] = blk["X"][pos].toarray()
    if gene_sets:
        _, _, matched_rna = _excluded_features(
            blocks["rna"]["names"], blocks["rna"]["chrom"], excluded_chromosomes=excluded, excluded_gene_sets=gene_sets
        )
        for set_name, info in matched_rna.items():
            if info["n_matched"] == 0:
                raise ValueError(
                    f"gene set {set_name!r} matches no RNA gene name (patterns "
                    f"{list(TECHNICAL_GENE_SETS[set_name])}); the feature index is not the expected one"
                )
    design = None
    if residualize is not None:
        design = _residualization_design(
            donors, counts_by_omic, blocks, residualize,
            excluded_chromosomes=excluded, excluded_gene_sets=gene_sets, pseudobulk_dir=pb_info.get("pseudobulk_dir"),
        )
        print(f"[{log_prefix}] {cohort}: residualization design {design.shape[0]} donors × "
              f"{design.shape[1]} covariates ({', '.join(design.columns)})")
    res = {}
    for o, blk in blocks.items():
        res[o] = preprocess_block(
            counts_by_omic[o], blk["names"], blk["chrom"], k=k, excluded_chromosomes=excluded,
            excluded_gene_sets=gene_sets, selection_rule=rules[o], normalization=norms[o],
            residualize=None if design is None else design.to_numpy(dtype=float),
        )
        findings += [f"{o}: {msg}" for msg in res[o]["findings"]]
        s = res[o]["steps"]
        print(f"[{log_prefix}] {cohort} {o}: " + " → ".join(str(v) for v in s.values())
              + ("" if rules[o] == "variance" else f" (selection rule {rules[o]})"))
        fit = res[o]["normalization_fit"]
        if fit is not None:
            print(f"[{log_prefix}] {cohort} {o}: {norms[o]} over {fit['n_features']} features (span {fit['span']:.4f}, "
                  f"{fit['iterations']} iterations); |slope of M on A| median "
                  f"{np.median(np.abs(fit['tilt_before']['slope'])):.4g} → {np.median(np.abs(fit['tilt_after']['slope'])):.4g}, "
                  f"largest absolute correction {fit['max_abs_correction'].max():.4g}")
        if debug_dir is not None:
            _write_debug(Path(debug_dir), cohort, o, donors, res[o])
    del counts_by_omic

    tmp = out_dir / f"{final_name}.tmp"
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)
    for o in ("rna", "atac"):
        df = pd.DataFrame(res[o]["values"], index=pd.Index(donors, name="donor_id"), columns=res[o]["kept"])
        df.to_csv(tmp / f"{o}.csv", float_format=EXPORT_FLOAT_FORMAT)
        _features_frame(blocks[o]["names"], blocks[o]["chrom"], res[o], blocks[o]["gene_ids"]).to_csv(
            tmp / f"features_{o}.csv", index=False, float_format=EXPORT_FLOAT_FORMAT
        )
    design_sha256 = None
    if design is not None:
        design.to_csv(tmp / RESIDUALIZATION_FILE, float_format="%.17g")
        design_sha256 = _sha256_file(tmp / RESIDUALIZATION_FILE)
    trend_sha256 = None
    if normalized_omics:
        _normalization_trend_table(donors, {o: res[o]["normalization_fit"] for o in normalized_omics}).to_csv(
            tmp / NORMALIZATION_TREND_FILE, index=False, float_format="%.17g"
        )
        trend_sha256 = _sha256_file(tmp / NORMALIZATION_TREND_FILE)
    meta = pd.DataFrame({
        "donor_id": donors,
        "sex": sub["sex"].to_numpy(),
        "stage": sub["stage"].to_numpy(),
        "stage_label": sub["stage"].map(STAGE_LABELS).to_numpy(),
    })
    meta.to_csv(tmp / "metadata.csv", index=False)
    cov = covariates.set_index("donor_id").loc[donors].reset_index()
    cov["n_rna"] = sub["n_rna"].to_numpy()
    cov["n_atac"] = sub["n_atac"].to_numpy()
    cov.to_csv(tmp / "covariates.csv", index=False)

    cov_cols = [col for col in cov.columns if col not in ("donor_id", "age_capped")]
    missingness = {
        col: sorted(cov.loc[cov[col].isna(), "donor_id"]) for col in cov_cols
    }
    rule_text = {
        "variance": "population variance (ddof 0)",
        "trend_ratio": (
            f"variance over a mean–variance trend (median variance of {TREND_N_BINS} equal-count bins of mean "
            "log2-CPM, linearly interpolated at the feature's mean, flat beyond the end bins)"
        ),
        "mean": "mean log2-CPM",
    }
    order_lines = ["restrict to cohort donors", "remove features on excluded chromosomes"]
    if gene_sets:
        order_lines.append(f"remove RNA features of the excluded gene sets ({', '.join(gene_sets)}; TECHNICAL_GENE_SETS)")
        order_lines.append("keep features with CPM >= min_cpm in >= k donors (raw library sizes after chromosome and gene-set removal)")
    else:
        order_lines.append("keep features with CPM >= min_cpm in >= k donors (raw library sizes after chromosome removal)")
    order_lines += [
        "TMM factors on the filtered counts (library size = filtered column sums)",
        "log2-CPM, prior count scaled by effective library size (edgeR cpm(log = TRUE))",
    ]
    if normalized_omics:
        order_lines.append(
            "trend-normalize the log2-CPM of the abundance-passing features ("
            + "; ".join(f"{o}: {norms[o]}" for o in ("rna", "atac"))
            + f"; limma's fast cyclic loess: {LOESS_ITERATIONS} iterations, each fitting a lowess curve of every donor's "
            "deviation from the per-feature mean over donors against that mean and subtracting it, adaptive span per block)"
        )
    if design is not None:
        order_lines.append(
            "replace each feature by its mean plus the residual of its OLS fit on an intercept and the "
            f"residualization covariates ({', '.join(design.columns)}); units stay log2-CPM"
        )
    if all(rule == "variance" for rule in rules.values()):
        order_lines.append("keep the N features with the highest population variance (ddof 0), ties by name; input order kept")
    else:
        order_lines.append(
            "keep the N features with the highest selection score per omic ("
            + "; ".join(f"{o}: {rules[o]}, {rule_text[rules[o]]}" for o in ("rna", "atac"))
            + "), ties by name; input order kept"
        )
    selection_block: dict = dict(rules)
    trend_bins = {o: res[o]["trend_bins"] for o in ("rna", "atac") if res[o]["trend_bins"] is not None}
    if trend_bins:
        selection_block["trend_bins"] = {
            o: {"n_bins": TREND_N_BINS, "bins": bins, "note": "equal-count bins of mean log2-CPM (np.array_split of the "
                "abundance-passing features sorted by mean); median_mean and median_variance are the interpolation nodes"}
            for o, bins in trend_bins.items()
        }
    residualization_block = None
    if design is not None:
        residualization_block = {
            "covariates": list(design.columns),
            "candidate_fractions": list(RESIDUALIZATION_CANDIDATE_FRACTIONS),
            "cell_state_exclusions": dict(RESIDUALIZATION_CELL_STATE_EXCLUSIONS),
            "design_file": RESIDUALIZATION_FILE,
            "design_sha256": design_sha256,
            "derivation": {
                "log10_lib_<omic>": "log10 of the donor's counts over the features kept at step 2 of this export "
                                    "(the step-3 library size, tmm.<omic>.<donor>.lib_size_filter)",
                "<set>_pseudobulk": f"donor_technical.csv: the UMI-weighted pseudobulk fraction of the gene set over the donor's {cell_type} nuclei",
            },
            "position": (
                "after log2-CPM (step 5) and the trend normalization, before the selection step (step 6)"
                if normalized_omics else "after log2-CPM (step 5), before the selection step (step 6)"
            ),
            "model": "per feature, OLS on an intercept and the covariates; exported value = feature mean + residual; "
                     "variance in features_<omic>.csv and the selection score are those of the residualized values",
            "degrees_of_freedom": int(design.shape[1] + 1),
            "label_free": "no sex, stage or analysis covariate enters the design",
            "r_squared_quantiles": {
                o: {f"q{int(round(q * 100)):02d}": float(np.nanquantile(res[o]["r_squared"], q)) for q in DIAGNOSTIC_QUANTILES}
                for o in ("rna", "atac")
            },
        }
    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "cell_type": cell_type,
        "cohort": cohort,
        "n_donors": int(len(donors)),
        "blocks": {
            o: {"file": f"{o}.csv", "units": "log2-CPM", "motco_transform": "identity",
                "n_features": int(res[o]["values"].shape[1]), "features_file": f"features_{o}.csv"}
            for o in ("rna", "atac")
        },
        "metadata": {
            "file": "metadata.csv",
            "sample_col": "donor_id",
            "group_col": "sex",
            "group_labels": sorted(sub["sex"].unique()),
            "stage_col": "stage",
            "stage_order": STAGES,
            "stage_labels": STAGE_LABELS,
            "cells": cells,
        },
        "preprocessing_order": order_lines,
        "thresholds": {
            "excluded_chromosomes": list(excluded),
            "min_cpm": EXPORT_MIN_CPM,
            "k": k,
            "k_source_cell": k_cell,
            "prior_count": EXPORT_PRIOR_COUNT,
            "N": EXPORT_TOP_N,
            "tmm": {"logratio_trim": TMM_LOGRATIO_TRIM, "sum_trim": TMM_SUM_TRIM,
                    "weighting": True, "reference": "upper-quartile factor closest to the mean"},
            "float_format": EXPORT_FLOAT_FORMAT,
        },
        **({"variant": dict(variant)} if variant is not None else {}),
        "feature_counts": {o: res[o]["steps"] for o in ("rna", "atac")},
        "removed_by_chromosome": {o: res[o]["removed_by_chromosome"] for o in ("rna", "atac")},
        "excluded_gene_sets": list(gene_sets),
        **({"removed_by_gene_set": {
            o: {
                s: {key: val for key, val in info.items() if key != "matched"}
                for s, info in res[o]["removed_by_gene_set"].items()
            }
            for o in ("rna", "atac")
        }} if gene_sets else {}),
        "selection_rule": selection_block,
        "normalization": dict(norms),
        **({"normalization_fit": _normalization_fit_block(
            {o: res[o]["normalization_fit"] for o in normalized_omics}, trend_sha256=trend_sha256,
        )} if normalized_omics else {}),
        "residualization": residualization_block,
        "kept_mean_log_cpm_range": {o: res[o]["mean_log_cpm_range"] for o in ("rna", "atac")},
        "tmm": {
            o: {
                d: {"lib_size_filter": int(res[o]["lib_size_filter"][i]),
                    "lib_size": int(res[o]["lib_size"][i]),
                    "factor": float(res[o]["tmm_factors"][i])}
                for i, d in enumerate(donors)
            }
            for o in ("rna", "atac")
        },
        "annotation": annotation,
        "input": dict(pb_info),
        "covariates": {
            "file": "covariates.csv",
            "applied": False,
            "note": "covariates exported, not applied (adjustment belongs to MOTCO)",
            "missing": missingness,
            "n_age_capped": int(cov["age_capped"].astype(bool).sum()),
        },
        "software": _software_versions(),
        "git": _git_revision(),
        "findings": findings,
    }
    (tmp / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    final = out_dir / final_name
    if final.exists():
        shutil.rmtree(final)
    os.replace(tmp, final)
    print(f"[{log_prefix}] Wrote {final}")
    for msg in findings:
        print(f"[{log_prefix}] FINDING: {cohort}: {msg}")
    return {"dir": final, "manifest": manifest, "findings": findings, "k": k, "cells": cells}


def run_export(
    *,
    cell_type: str = DEFAULT_CELL_TYPE,
    pseudobulk_dir: str | Path | None = None,
    out_dir: str | Path | None = None,
    cohort: str = "all",
    data_root: str | Path | None = None,
    force: bool = False,
    debug_dir: str | Path | None = None,
) -> dict:
    """Export MOTCO-ready RNA and ATAC blocks per cohort.

    Parameters
    ----------
    cell_type:
        Registered cell type (default ``astrocyte``): the default
        directories, and the cell type the pseudobulk manifest must record.
    pseudobulk_dir:
        Output directory of ``pseudobulk`` and ``covariates`` (default the
        cell type's, ``results/pseudobulk`` for ``astrocyte``).
    out_dir:
        Export root; each cohort goes to ``<out_dir>/<cohort>/`` (default
        the cell type's ``motco_input``).
    cohort:
        ``"cohort_30"``, ``"cohort_50"`` or ``"all"``.
    data_root:
        Data root holding (or receiving) the gene annotation.
    force:
        Replace existing cohort directories.
    debug_dir:
        If given, also write the filtered counts, TMM factors and filtered
        log2-CPM per cohort and omic there (for ``scripts/check_tmm_edger.R``).

    Returns
    -------
    dict
        Cohort → ``{"dir": Path, "steps": {omic: counts}, "findings": [...],
        "manifest": dict, "k": int, "cells": dict}``.

    Notes
    -----
    Per cohort and omic: restrict to the cohort → drop chrX/Y/M → CPM ≥ 1 in
    ≥ k donors → TMM → log2-CPM (prior 2) → top N by variance. Covariates are
    exported alongside, not applied. Each cohort is written to
    ``<cohort>.tmp/`` and renamed when complete (``_export_cohort``, shared
    with ``run_preprocessing_sensitivity``, which calls it with one setting
    changed).
    """
    cell_type_entry(cell_type)
    pb = _resolve_dir(pseudobulk_dir, cell_type, "pseudobulk")
    out_dir = _resolve_dir(out_dir, cell_type, "motco_input")
    cohorts = list(EXPORT_COHORTS) if cohort == "all" else [cohort]
    for c in cohorts:
        if c not in EXPORT_COHORTS:
            raise ValueError(f"Unknown cohort {c!r}; expected one of {EXPORT_COHORTS} or 'all'")
    pb_manifest = check_export_inputs(pb, cohorts, cell_type=cell_type)
    input_findings = check_manifest_cell_type(pb_manifest, pb / "manifest.json", cell_type)
    for c in cohorts:
        if (out_dir / c).exists() and not force:
            raise FileExistsError(f"{out_dir / c} exists; rerun with --force to replace it")
    print(f"[export] Cell type: {cell_type}")
    print(f"[export] Pseudobulk directory: {pb}")
    print(f"[export] Output directory: {out_dir}")

    inputs = _load_export_inputs(pb, data_root)
    results = {}
    for c in cohorts:
        sub = _cohort_rows(inputs["table"], inputs["covariates"], c, pb)
        res = _export_cohort(
            c, sub, inputs["blocks"], inputs["covariates"], inputs["annotation"], inputs["pb_info"],
            out_dir=out_dir, final_name=c, excluded_chromosomes=EXPORT_EXCLUDED_CHROMOSOMES,
            debug_dir=None if debug_dir is None else Path(debug_dir),
            cell_type=cell_type, input_findings=input_findings,
        )
        results[c] = {"steps": res["manifest"]["feature_counts"], **res}

    print("[export] Features after each step:")
    rows = [
        {"cohort": c, "omic": o, **steps}
        for c, r in results.items() for o, steps in r["steps"].items()
    ]
    print(pd.DataFrame(rows).to_string(index=False))
    return results


def _write_debug(debug_dir: Path, cohort: str, omic: str, donors: np.ndarray, res: dict) -> None:
    """Write filtered counts, TMM factors and log2-CPM for the edgeR cross-check.

    For a trend-normalized omic also ``<cohort>_<omic>_normalized.csv`` (the
    normalized log2-CPM, before any residualization) and
    ``<cohort>_<omic>_loess.csv`` (span, iterations, number of features) for
    ``scripts/check_cyclic_loess_limma.R``. ``_logcpm.csv`` is always the
    step-5 log2-CPM.
    """
    debug_dir.mkdir(parents=True, exist_ok=True)
    stem = debug_dir / f"{cohort}_{omic}"
    idx = pd.Index(donors, name="donor_id")
    pd.DataFrame(res["filtered_counts"].astype(np.int64), index=idx, columns=res["filtered_names"]).to_csv(
        f"{stem}_counts.csv"
    )
    pd.DataFrame({"lib_size": res["lib_size"], "factor": res["tmm_factors"]}, index=idx).to_csv(
        f"{stem}_tmm.csv", float_format="%.17g"
    )
    pd.DataFrame(res["filtered_log_cpm"], index=idx, columns=res["filtered_names"]).to_csv(
        f"{stem}_logcpm.csv", float_format="%.17g"
    )
    fit = res.get("normalization_fit")
    if fit is not None:
        pd.DataFrame(res["filtered_normalized"], index=idx, columns=res["filtered_names"]).to_csv(
            f"{stem}_normalized.csv", float_format="%.17g"
        )
        pd.DataFrame([{
            "rule": res["normalization"], "span": fit["span"], "iterations": fit["iterations"],
            "robustness_iterations": fit["robustness_iterations"], "n_features": fit["n_features"],
        }]).to_csv(f"{stem}_loess.csv", index=False, float_format="%.17g")


# ---------------------------------------------------------------------------
# MOTCO trajectory analysis (covariate-adjusted, per cohort)
# ---------------------------------------------------------------------------

ANALYSIS_BLOCKS = ("rna", "atac")
# Entered additively, in this order
ANALYSIS_COVARIATES = (
    "age_at_death", "pmi", "rin", "years_education", "log10_n_rna", "log10_n_atac",
)
# Derived column → exported column it is the base-10 logarithm of
ANALYSIS_LOG10_COLUMNS: dict[str, str] = {"log10_n_rna": "n_rna", "log10_n_atac": "n_atac"}
# Exported covariate columns that do not enter the model, with the reason
ANALYSIS_EXCLUDED_COVARIATES: dict[str, str] = {
    "education_level": "not entered by decision: categorical, it would spend one design column per level",
    "age_capped": "constant in this release (no capped ages); MOTCO rejects a constant covariate",
}
ANALYSIS_MISSING_COVARIATES = "error"
ANALYSIS_PERMUTATIONS = 999
ANALYSIS_SEED = 20261019
ANALYSIS_ALPHA = 0.05
# As in MOTCO's covariate Type I check; exact size alpha for (1 + b) / (1 + B) p-values
ANALYSIS_REJECTION_RULE = "p <= alpha"
ANALYSIS_MOTCO_MIN_VERSION = "0.8.0"
ANALYSIS_STATISTICS = ("delta", "angle", "shape")

_ANALYSIS_INPUTS = ("rna.csv", "atac.csv", "metadata.csv", "covariates.csv", "manifest.json")
_ANALYSIS_DESIGN_COLUMNS = ("donor_id", "sex", "stage", "stage_label")


def check_motco_version(
    *,
    required: str = ANALYSIS_MOTCO_MIN_VERSION,
    purpose: str = "the analysis needs",
    feature: str = "covariate-adjusted analyze_trajectories",
) -> str:
    """Return the installed MOTCO version; raise if it is older than the minimum.

    Parameters
    ----------
    required:
        The minimum version; ``ANALYSIS_MOTCO_MIN_VERSION`` by default, and
        ``LATENT_MOTCO_MIN_VERSION`` for ``latent-sensitivity``.
    purpose, feature:
        Named in the error message (what needs the version, and why).

    Raises
    ------
    RuntimeError
        The installed MOTCO is older than ``required``; the message names both
        versions.
    """
    import motco

    installed = str(motco.__version__)
    try:
        from packaging.version import Version

        too_old = Version(installed) < Version(required)
    except ImportError:
        def key(v: str) -> tuple[int, ...]:
            return tuple(int(p) for p in re.findall(r"\d+", v)[:3])

        too_old = key(installed) < key(required)
    if too_old:
        raise RuntimeError(
            f"MOTCO {installed} is installed but {purpose} MOTCO >= {required} "
            f"({feature}); run `uv sync --locked`"
        )
    return installed


def load_analysis_inputs(cohort_dir: str | Path, *, extra_covariates: pd.DataFrame | None = None) -> dict:
    """Read one cohort's export, cross-check it and build the model's metadata table.

    Parameters
    ----------
    cohort_dir:
        ``<export root>/<cohort>/`` as written by ``run_export``.
    extra_covariates:
        Covariates from outside the export (``latent-sensitivity``: the
        technical table's fractions): a table with a ``donor_id`` column and
        one column per covariate, merged by ``donor_id`` after the log10
        derivation and appended after ``ANALYSIS_COVARIATES`` in its column
        order. Every export donor needs a row (extra rows are ignored).
        ``None`` (the default) gives the metadata table of the ``analyze`` step.

    Returns
    -------
    dict
        ``blocks`` (omic → donors × features DataFrame, read as ``motco analyze``
        reads ``--block`` files), ``metadata`` (design columns plus
        ``ANALYSIS_COVARIATES`` and the extra covariates, one row per donor
        sorted by ``donor_id``), ``export_manifest`` and ``inputs`` (file →
        path and sha256 of the five export files).

    Raises
    ------
    FileNotFoundError
        The cohort directory or one of its five files is missing.
    ValueError
        The tables list different donors or disagree with the export manifest;
        or ``extra_covariates`` lacks ``donor_id``, repeats a donor, has no row
        for an export donor, or names a column the metadata table already has.

    Notes
    -----
    Missing covariate values are passed through: MOTCO's
    ``missing_covariates="error"`` rejects them, naming the column and donors.
    """
    d = Path(cohort_dir)
    if not d.is_dir():
        raise FileNotFoundError(f"{d} does not exist; run `sead-sexdiff export` to produce it")
    for name in _ANALYSIS_INPUTS:
        if not (d / name).exists():
            raise FileNotFoundError(f"{d / name} is missing; run `sead-sexdiff export` to produce it")
    inputs = {name: {"path": str(d / name), "sha256": _sha256_file(d / name)} for name in _ANALYSIS_INPUTS}
    manifest = json.loads((d / "manifest.json").read_text())
    meta = pd.read_csv(d / "metadata.csv", dtype={"donor_id": str, "sex": str, "stage": str})
    cov = pd.read_csv(d / "covariates.csv", dtype={"donor_id": str, "education_level": str})
    blocks = {}
    for o in ANALYSIS_BLOCKS:
        values = pd.read_csv(d / f"{o}.csv", index_col=0)
        values.index = values.index.astype(str)
        blocks[o] = values

    tables = {
        "metadata.csv": meta["donor_id"].tolist(),
        "covariates.csv": cov["donor_id"].tolist(),
        **{f"{o}.csv": b.index.tolist() for o, b in blocks.items()},
    }
    for name, ids in tables.items():
        repeated = sorted({i for i in ids if ids.count(i) > 1})
        if repeated:
            raise ValueError(f"{d / name}: repeated donor_id(s) {repeated}")
    reference = set(tables["metadata.csv"])
    for name, ids in tables.items():
        if set(ids) != reference:
            raise ValueError(
                f"{d / 'metadata.csv'} and {d / name} list different donors "
                f"(only in metadata.csv: {sorted(reference - set(ids))}; "
                f"only in {name}: {sorted(set(ids) - reference)})"
            )
    donors = sorted(reference)
    for o, b in blocks.items():
        # analysis_metadata.csv is sorted by donor_id and `motco analyze` needs the
        # block rows in metadata order, so the blocks must already be sorted
        if b.index.tolist() != donors:
            raise ValueError(f"{d / f'{o}.csv'}: rows are not sorted by donor_id (as `export` writes them)")
    if manifest.get("n_donors") != len(donors):
        raise ValueError(
            f"{d / 'manifest.json'} records n_donors = {manifest.get('n_donors')} "
            f"but {d / 'metadata.csv'} lists {len(donors)} donors"
        )
    for o, b in blocks.items():
        n = manifest.get("blocks", {}).get(o, {}).get("n_features")
        if n != b.shape[1]:
            raise ValueError(
                f"{d / 'manifest.json'} records n_features = {n} for {o} but {d / f'{o}.csv'} has {b.shape[1]}"
            )

    needed = [c for c in ANALYSIS_COVARIATES if c not in ANALYSIS_LOG10_COLUMNS]
    needed += list(ANALYSIS_LOG10_COLUMNS.values())
    absent = [c for c in needed if c not in cov.columns]
    if absent:
        raise ValueError(f"{d / 'covariates.csv'} lacks column(s) {absent}")
    absent = [c for c in _ANALYSIS_DESIGN_COLUMNS if c not in meta.columns]
    if absent:
        raise ValueError(f"{d / 'metadata.csv'} lacks column(s) {absent}")
    merged = meta.merge(cov, on="donor_id", how="inner", validate="one_to_one")
    for derived, source in ANALYSIS_LOG10_COLUMNS.items():
        counts = merged[source].astype(float)
        bad = merged.loc[counts <= 0, "donor_id"]
        if len(bad):
            raise ValueError(f"{d / 'covariates.csv'}: {source} must be positive for log10; donor(s) {sorted(bad)}")
        merged[derived] = np.log10(counts)
    merged = merged[list(_ANALYSIS_DESIGN_COLUMNS) + list(ANALYSIS_COVARIATES)]
    if extra_covariates is not None:
        extra = extra_covariates.copy()
        if "donor_id" not in extra.columns:
            raise ValueError("extra_covariates lacks a donor_id column")
        extra["donor_id"] = extra["donor_id"].astype(str)
        clash = [c for c in extra.columns if c != "donor_id" and c in merged.columns]
        if clash:
            raise ValueError(f"extra_covariates column(s) {clash} are already in the metadata table of {d}")
        repeated = sorted(extra.loc[extra["donor_id"].duplicated(), "donor_id"].unique())
        if repeated:
            raise ValueError(f"extra_covariates: repeated donor_id(s) {repeated}")
        absent = sorted(set(donors) - set(extra["donor_id"]))
        if absent:
            raise ValueError(f"extra_covariates has no row for export donor(s) {absent} of {d}")
        merged = merged.merge(extra, on="donor_id", how="left", validate="one_to_one")
    merged = merged.sort_values("donor_id").reset_index(drop=True)
    return {"blocks": blocks, "metadata": merged, "export_manifest": manifest, "inputs": inputs}


def _fmt(value, spec: str = ".6g") -> str:
    return "n/a" if value is None else format(value, spec)


def render_analysis_summary(record: dict, *, cohort: str, alpha: float = ANALYSIS_ALPHA) -> str:
    """Render ``summary.md`` from a MOTCO analysis record (``result.json``).

    Parameters
    ----------
    record:
        The parsed ``result.json``.
    cohort:
        Cohort name for the title.
    alpha:
        Significance level; a statistic is rejected when ``p <= alpha``.

    Returns
    -------
    str
        Markdown: statistics (observed, RRPP p, null mean/sd/q95, CV-selected
        rank in the ``angle`` row), cell sizes, latent block and covariate block.
        For a forced-rank record (``latent.component_selection == "forced"``)
        the ``angle`` row names the forced rank and says the CV-selected
        analysis is primary, and the latent block says cross-validation was
        skipped; the output for CV records is unchanged. A record with a
        ``latent_adjustment`` block (MOTCO ≥ 0.9.0) gets a "Latent adjustment"
        section (conditioning, residual df, smallest CV training df, median
        feature R² and the fraction above 0.5); records without one are
        rendered as before.
    """
    design = record["design"]
    meas = record["measurement"]
    latent = record["latent"]
    forced = latent["component_selection"] == "forced"
    stages = design["stage_levels"]
    lines = [
        f"# MOTCO trajectory analysis: {cohort}",
        "",
        f"MOTCO {record['software']['motco']}; n = {design['n_samples']}; groups "
        f"{' vs '.join(design['group_levels'])}; stages {' → '.join(stages)}; orientation "
        f"`{meas['orientation_estimator']}`; {meas['permutations']} RRPP permutations, seed {meas['seed']}; "
        f"α = {alpha} (rejected when p ≤ α).",
        "",
        "## Statistics",
        "",
        "| statistic | observed | p (RRPP) | rejected | null mean | null sd | null q95 | note |",
        "|---|---:|---:|:---:|---:|---:|---:|---|",
    ]
    for s in ANALYSIS_STATISTICS:
        p = record["p_values"].get(s)
        null = record["null_summary"].get(s, {})
        rejected = "n/a" if p is None else ("yes" if p <= alpha else "no")
        note = ""
        if s == "angle" and forced:
            note = (
                f"degrees; forced rank {latent['selected_lv']} (dimensionality sensitivity; "
                "the CV-selected analysis is primary)"
            )
        elif s == "angle":
            note = f"degrees; CV-selected rank {latent['selected_lv']}"
            if latent["selected_lv"] == 1:
                note += " (at rank 1 the angle is 0° or 180° by construction)"
        lines.append(
            f"| {s} | {_fmt(record['statistics'].get(s))} | {_fmt(p, '.4g')} | {rejected} | "
            f"{_fmt(null.get('mean'))} | {_fmt(null.get('sd'))} | {_fmt(null.get('q95'))} | {note} |"
        )
    lines += [
        "",
        "## Cell sizes",
        "",
        "| group | " + " | ".join(f"stage {s}" for s in stages) + " | total |",
        "|---|" + "---:|" * (len(stages) + 1),
    ]
    for g in design["group_levels"]:
        cells = design["cell_sizes"][g]
        lines.append(f"| {g} | " + " | ".join(str(cells[s]) for s in stages) + f" | {sum(cells.values())} |")
    cv = latent.get("integration_params", {})
    if forced:
        cv_line = "- cross-validation: skipped (rank forced)"
    else:
        cv_line = (
            f"- CV: cv2_splits {cv.get('cv2_splits')}, cv1_splits {cv.get('cv1_splits')}, n_repeats "
            f"{cv.get('n_repeats')} ({latent.get('cv_repeats_completed')} completed), max_components "
            f"{cv.get('max_components')}, random_state {cv.get('random_state')}"
        )
    lines += [
        "",
        "## Latent space",
        "",
        f"- PLS-DA rank: {latent['selected_lv']} (`{latent['component_selection']}`)",
        f"- cv_mean_auroc: {_fmt(latent.get('cv_mean_auroc'), '.4f')}",
        cv_line,
        "- blocks: " + ", ".join(
            f"{b['name']} ({b['n_features']} features, {b['transform']})" for b in record["blocks"]
        ),
    ]
    cov = record.get("covariates")
    lines += ["", "## Covariates", ""]
    if cov is None:
        lines.append("None entered.")
    else:
        dropped = cov["dropped_sample_ids"]
        lines += [
            "- entered (in order): " + ", ".join(f"`{c}` ({cov['kinds'][c]})" for c in cov["columns"]),
            f"- design columns: {len(cov['design_columns'])} (" + ", ".join(cov["design_columns"]) + ")",
            f"- missing policy: `{cov['missing_policy']}`; dropped donors: {cov['n_dropped']}"
            + (f" ({', '.join(dropped)})" if dropped else ""),
            f"- LS-means evaluated at the reference row (`{cov['lsmean_reference']}`):",
            "",
            "| design column | reference value |",
            "|---|---:|",
        ]
        lines += [f"| {c} | {_fmt(v)} |" for c, v in cov["reference_row"].items()]
    adjustment = record.get("latent_adjustment")
    if adjustment is not None:
        r2 = adjustment.get("feature_r2", {})
        lines += [
            "",
            "## Latent adjustment",
            "",
            f"- basis fitted on the {adjustment['feature_scale']} features residualized on the covariates "
            + ", ".join(f"`{c}`" for c in adjustment["covariates"])
            + f", conditional on `{adjustment['conditioning']}`; the trajectory model, LS-means and RRPP are those "
            "of the covariate block above",
            f"- residualization design columns beyond the intercept: {len(adjustment['design_columns'])} ("
            + ", ".join(adjustment["design_columns"]) + f"); n = {adjustment['n_samples']}",
            f"- residual df: {adjustment['residual_df']}; smallest CV training df: "
            + ("n/a (rank forced)" if adjustment.get("cv_min_training_df") is None else str(adjustment["cv_min_training_df"])),
            f"- per-feature partial R² of the covariates given the conditioning: median {_fmt(r2.get('q50'), '.4f')}, "
            f"fraction above 0.5 {_fmt(r2.get('fraction_above_0.5'), '.4f')}",
        ]
    return "\n".join(lines) + "\n"


def _cli_equivalent(
    cohort_dir: Path, metadata_path: Path, params, n_jobs: int, *, covariates: Iterable[str] = ANALYSIS_COVARIATES,
) -> dict:
    """The ``motco analyze`` command that reproduces one analysis directory (CV or forced rank, optionally latent-adjusted)."""
    import shlex

    argv = ["motco", "analyze"]
    for o in ANALYSIS_BLOCKS:
        argv += ["--block", f"{o}={cohort_dir / f'{o}.csv'}"]
    argv += [
        "--metadata", str(metadata_path),
        "--sample-col", "donor_id", "--group-col", "sex", "--stage-col", "stage",
        "--stage-order", ",".join(STAGES),
    ]
    for c in covariates:
        argv += ["--covariate", c]
    argv += ["--missing-covariates", ANALYSIS_MISSING_COVARIATES]
    latent = getattr(params, "covariate_latent_adjustment", None)
    if latent is not None:
        argv += ["--latent-adjustment", latent]
    argv += [
        "--permutations", str(params.permutations), "--seed", str(params.seed), "--n-jobs", str(n_jobs),
        "--orientation-estimator", params.orientation_estimator,
    ]
    if params.forced_components is not None:
        # Cross-validation does not run, so none of its knobs apply
        argv += ["--forced-components", str(params.forced_components)]
    else:
        argv += [
            "--cv1-splits", str(params.cv1_splits), "--cv2-splits", str(params.cv2_splits),
            "--n-repeats", str(params.n_repeats), "--max-components", str(params.max_components),
            "--cv-random-state", str(params.random_state),
        ]
    argv += ["--out-dir", "OUT_DIR"]
    return {
        "command": shlex.join(argv),
        "argv": argv,
        "note": "replace OUT_DIR with a new directory; paths are relative to the directory run_analysis ran in",
    }


def _analyze_cohort(
    cohort: str,
    cohort_dir: Path,
    loaded: dict,
    params,
    *,
    out_dir: Path,
    final_name: str,
    n_jobs: int,
    findings: list[str] | None = None,
    published_dir: Path | None = None,
    log_prefix: str = "analyze",
    covariates: Iterable[str] = ANALYSIS_COVARIATES,
    covariate_sources: dict[str, str] | None = None,
    extra_inputs: dict | None = None,
    cell_type: str | None = None,
) -> dict:
    """Run MOTCO on one loaded cohort and write a complete analysis directory.

    Shared by ``run_analysis`` (CV-selected rank), ``run_rank_sensitivity``
    (forced rank), ``run_preprocessing_sensitivity`` (a variant export) and
    ``run_latent_sensitivity`` (another covariate set, optionally with the
    latent adjustment): they differ only in ``params``, the covariate set and
    where the directory goes, so every directory is a primary directory with
    the named setting(s) changed.

    Parameters
    ----------
    cohort:
        Cohort name (title and manifest).
    cohort_dir:
        The export directory the inputs were read from.
    loaded:
        ``load_analysis_inputs(cohort_dir)`` (with ``extra_covariates`` when
        ``covariates`` goes beyond ``ANALYSIS_COVARIATES``).
    params:
        ``motco.analysis.TrajectoryAnalysisParams``; ``forced_components`` set
        for a forced-rank run, ``None`` for CV rank selection;
        ``covariate_latent_adjustment`` set for a latent-adjusted run (recorded
        as ``measurement.latent_adjustment`` and ``--latent-adjustment`` in
        ``cli_equivalent``; both absent when it is ``None``, as MOTCO leaves
        the setting out of its record).
    out_dir:
        The directory is written to ``<out_dir>/<final_name>.tmp/`` and renamed
        to ``<out_dir>/<final_name>`` when complete.
    final_name:
        Name of the finished directory under ``out_dir``.
    n_jobs:
        Recorded in the manifest and in ``cli_equivalent``.
    findings:
        Findings recorded before the helper's own (capped ages, ``n_jobs``,
        rank 1).
    published_dir:
        The path the directory will have once the caller is done, when that
        differs from ``<out_dir>/<final_name>`` (the sweep writes its ranks
        under ``<cohort>.tmp/`` first); used for ``cli_equivalent``.
    log_prefix:
        Prefix of the progress lines, ``[analyze]`` by default.
    covariates:
        The covariates entered, additively and in this order; columns of
        ``loaded["metadata"]``. ``ANALYSIS_COVARIATES`` by default.
    covariate_sources:
        Covariate → where it comes from, for covariates outside the export;
        recorded as ``covariates.sources`` (absent when ``None``).
    extra_inputs:
        Inputs beside the five export files (``latent-sensitivity``: the
        technical table), merged into the manifest's ``inputs`` block next to
        ``files``, which stays the five export files so ``axis-diagnostics``
        can check the directory against the export.
    cell_type:
        Recorded as the manifest's ``cell_type`` (checked against the export
        manifest by the caller). ``None`` records the export manifest's,
        ``astrocyte`` when it has none.

    Returns
    -------
    dict
        ``{"dir": Path, "record": dict, "manifest": dict, "findings": [...]}``.
    """
    import io
    import shutil
    from importlib.metadata import PackageNotFoundError, version

    from motco.analysis import OmicsBlock, TrajectoryDataset, analyze_trajectories, write_analysis_outputs

    motco_version = check_motco_version()
    covariates = tuple(covariates)
    blocks = loaded["blocks"]
    export_manifest = loaded["export_manifest"]
    if cell_type is None:
        cell_type = export_manifest.get("cell_type", DEFAULT_CELL_TYPE)
    final = out_dir / final_name
    published = Path(published_dir) if published_dir is not None else final
    absent = [c for c in covariates if c not in loaded["metadata"].columns]
    if absent:
        raise ValueError(f"covariate(s) {absent} are not in the metadata table loaded from {cohort_dir}")
    # Parse the table back as `motco analyze --metadata` does, so the
    # CLI-equivalent run sees exactly what this one sees
    metadata_csv = loaded["metadata"].to_csv(index=False)
    metadata = pd.read_csv(io.StringIO(metadata_csv))
    metadata["donor_id"] = metadata["donor_id"].astype(str)
    forced = params.forced_components
    latent = getattr(params, "covariate_latent_adjustment", None)
    print(
        f"[{log_prefix}] {cohort}: {len(metadata)} donors, "
        + " + ".join(f"{b.shape[1]} {o}" for o, b in blocks.items())
        + f" features, covariates {', '.join(covariates)}"
    )

    dataset = TrajectoryDataset(
        blocks=[OmicsBlock(o, blocks[o], "identity") for o in ANALYSIS_BLOCKS],
        metadata=metadata,
        sample_col="donor_id",
        group_col="sex",
        stage_col="stage",
        stage_order=STAGES,
        covariate_cols=covariates,
        missing_covariates=ANALYSIS_MISSING_COVARIATES,
    )
    selection = "CV rank selection" if forced is None else f"forced rank {forced}"
    if latent is not None:
        selection += f" with the {latent}-conditioned covariate latent adjustment"
    print(f"[{log_prefix}] {cohort}: {selection} + {params.permutations} RRPP permutations (n_jobs {n_jobs})...")
    result = analyze_trajectories(dataset, params)
    record = result.record
    print(
        f"[{log_prefix}] {cohort}: rank {record['latent']['selected_lv']} ({record['latent']['component_selection']}), "
        f"{record['measurement']['runtime_seconds']:.1f} s"
    )

    findings = list(findings or [])
    n_capped = export_manifest.get("covariates", {}).get("n_age_capped")
    if n_capped:
        findings.append(
            f"{n_capped} donor(s) have capped ages but age_capped is not entered "
            "(the covariate set assumed it constant); revisit the covariate set"
        )
    if n_jobs != 1:
        findings.append(
            f"n_jobs = {n_jobs}: MOTCO's parallel RRPP seeds one stream per worker, so the permutation "
            "draws, p-values and null summaries differ from an n_jobs = 1 run (statistics, rank, latent "
            "scores and loadings do not); the case-study results are reported at n_jobs = 1"
        )
    if record["latent"]["selected_lv"] == 1:
        findings.append(
            "CV selected rank 1: in a one-dimensional latent space the net displacements are parallel "
            "or antiparallel, so angle is 0° or 180° by construction and carries no orientation information"
        )

    tmp = out_dir / f"{final_name}.tmp"
    if tmp.exists():
        shutil.rmtree(tmp)
    written = write_analysis_outputs(result, tmp)
    (tmp / "analysis_metadata.csv").write_text(metadata_csv)
    (tmp / "summary.md").write_text(render_analysis_summary(record, cohort=cohort))
    software = _software_versions()
    software["motco"] = motco_version
    try:
        software["scikit-learn"] = version("scikit-learn")
    except PackageNotFoundError:
        software["scikit-learn"] = None
    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "cell_type": cell_type,
        "cohort": cohort,
        "n_donors": int(len(metadata)),
        "inputs": {
            "dir": str(cohort_dir),
            "files": loaded["inputs"],
            "export_manifest_sha256": loaded["inputs"]["manifest.json"]["sha256"],
            "export_manifest_created_at": export_manifest.get("created_at"),
            "export_git": export_manifest.get("git"),
            **(extra_inputs or {}),
        },
        "model": {
            "entry_point": "motco.analysis.analyze_trajectories",
            "blocks": [
                {"name": o, "file": f"{o}.csv", "transform": "identity", "n_features": int(blocks[o].shape[1])}
                for o in ANALYSIS_BLOCKS
            ],
            "metadata_file": "analysis_metadata.csv",
            "sample_col": "donor_id",
            "group_col": "sex",
            "group_levels": record["design"]["group_levels"],
            "stage_col": "stage",
            "stage_order": list(STAGES),
            "stage_labels": STAGE_LABELS,
        },
        "preprocessing": {
            "source": "export manifest",
            "thresholds": export_manifest.get("thresholds"),
            "excluded_gene_sets": export_manifest.get("excluded_gene_sets"),
            "selection_rule": {
                o: rule for o, rule in (export_manifest.get("selection_rule") or {}).items() if o in ANALYSIS_BLOCKS
            } or None,
            "normalization": export_manifest.get("normalization"),
            "normalization_fit": (
                None if not export_manifest.get("normalization_fit")
                else {
                    o: {key: fit.get(key) for key in ("rule", "reference", "span", "iterations")}
                    for o, fit in export_manifest["normalization_fit"].items() if o in ANALYSIS_BLOCKS
                }
            ),
            "residualization": (
                None if not export_manifest.get("residualization")
                else {key: export_manifest["residualization"].get(key) for key in ("covariates", "design_sha256", "degrees_of_freedom")}
            ),
            "feature_counts": export_manifest.get("feature_counts"),
            "n_age_capped": n_capped,
            "export_findings": export_manifest.get("findings", []),
        },
        "covariates": {
            "entered": list(covariates),
            "derived": {name: f"log10({src})" for name, src in ANALYSIS_LOG10_COLUMNS.items()},
            **({"sources": dict(covariate_sources)} if covariate_sources else {}),
            "excluded": dict(ANALYSIS_EXCLUDED_COVARIATES),
            "missing_policy": ANALYSIS_MISSING_COVARIATES,
            "imputed": False,
        },
        "measurement": {
            "orientation_estimator": params.orientation_estimator,
            "component_selection": record["latent"]["component_selection"],
            "forced_components": params.forced_components,
            **(
                {"cv": None, "cv_note": "cross-validation skipped: the latent rank was forced"}
                if forced is not None else
                {"cv": {
                    "cv1_splits": params.cv1_splits, "cv2_splits": params.cv2_splits,
                    "n_repeats": params.n_repeats, "max_components": params.max_components,
                    "random_state": params.random_state,
                }}
            ),
            **({"latent_adjustment": latent} if latent is not None else {}),
            "permutations": params.permutations,
            "seed": params.seed,
            "alpha": ANALYSIS_ALPHA,
            "rejection_rule": ANALYSIS_REJECTION_RULE,
            "n_jobs": n_jobs,
            "n_jobs_note": "RRPP draws depend on n_jobs (one seed stream per worker); rerun with the same n_jobs",
            "attribution": False,
        },
        "cli_equivalent": _cli_equivalent(
            cohort_dir, published / "analysis_metadata.csv", params, n_jobs, covariates=covariates
        ),
        "outputs": sorted([p.name for p in written.values()] + [
            "analysis_metadata.csv", "summary.md", "analysis_manifest.json",
        ]),
        "software": software,
        "git": _git_revision(),
        "findings": findings,
    }
    (tmp / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if final.exists():
        shutil.rmtree(final)
    os.replace(tmp, final)
    print(f"[{log_prefix}] Wrote {final}")
    for msg in findings:
        print(f"[{log_prefix}] FINDING: {cohort}: {msg}")
    return {"dir": final, "record": record, "manifest": manifest, "findings": findings}


def run_analysis(
    cohort: str = "all",
    *,
    cell_type: str = DEFAULT_CELL_TYPE,
    input_dir: str | Path | None = None,
    out_dir: str | Path | None = None,
    n_jobs: int = 1,
    force: bool = False,
) -> dict:
    """Run the covariate-adjusted MOTCO trajectory analysis per cohort.

    Parameters
    ----------
    cohort:
        ``"cohort_30"``, ``"cohort_50"`` or ``"all"``.
    cell_type:
        Registered cell type (default ``astrocyte``): the default
        directories, and the cell type every export manifest must record.
    input_dir:
        Export root; each cohort is read from ``<input_dir>/<cohort>/``
        (default the cell type's ``motco_input``).
    out_dir:
        Analysis root; each cohort goes to ``<out_dir>/<cohort>/`` (default
        the cell type's ``motco_analysis``).
    n_jobs:
        Parallel jobs for MOTCO's CV and RRPP (``-1`` = all CPUs). MOTCO 0.8.0
        seeds one RRPP stream per worker, so values other than 1 change the
        permutation draws (p-values, null summaries) but not the statistics,
        rank or latent space; the reported results use 1.
    force:
        Replace existing cohort directories.

    Returns
    -------
    dict
        Cohort → ``{"dir": Path, "record": dict, "manifest": dict, "findings": [...]}``.

    Notes
    -----
    Blocks ``rna`` and ``atac`` in identity units, group ``sex``, stage
    ``stage`` in the order 0, 1, 2, covariates ``ANALYSIS_COVARIATES``
    (missing policy ``error``), MOTCO's default net-displacement orientation
    and double-CV rank, ``ANALYSIS_PERMUTATIONS`` RRPP permutations with
    ``ANALYSIS_SEED``. The model's metadata table is serialized to
    ``analysis_metadata.csv`` and parsed back exactly as ``motco analyze``
    parses ``--metadata``, so the manifest's ``cli_equivalent`` reproduces the
    run. Each cohort is written to ``<cohort>.tmp/`` and renamed when complete
    (``_analyze_cohort``, shared with ``run_rank_sensitivity``).
    """
    motco_version = check_motco_version()
    from motco.analysis import TrajectoryAnalysisParams

    cell_type_entry(cell_type)
    input_dir = _resolve_dir(input_dir, cell_type, "motco_input")
    out_dir = _resolve_dir(out_dir, cell_type, "motco_analysis")
    cohorts = list(EXPORT_COHORTS) if cohort == "all" else [cohort]
    for c in cohorts:
        if c not in EXPORT_COHORTS:
            raise ValueError(f"Unknown cohort {c!r}; expected one of {EXPORT_COHORTS} or 'all'")
    for c in cohorts:
        if (out_dir / c).exists() and not force:
            raise FileExistsError(f"{out_dir / c} exists; rerun with --force to replace it")
    for c in cohorts:
        if not (input_dir / c).is_dir():
            raise FileNotFoundError(f"{input_dir / c} does not exist; run `sead-sexdiff export` to produce it")
    print(f"[analyze] MOTCO {motco_version}")
    print(f"[analyze] Cell type: {cell_type}")
    print(f"[analyze] Input directory: {input_dir}")
    print(f"[analyze] Output directory: {out_dir}")

    # Every cohort's export is read and its cell type checked before anything is written
    loaded_by_cohort = {}
    for c in cohorts:
        loaded = load_analysis_inputs(input_dir / c)
        findings = check_manifest_cell_type(loaded["export_manifest"], input_dir / c / "manifest.json", cell_type)
        loaded_by_cohort[c] = (loaded, findings)
    results = {}
    for c in cohorts:
        cohort_dir = input_dir / c
        loaded, findings = loaded_by_cohort.pop(c)
        params = TrajectoryAnalysisParams(permutations=ANALYSIS_PERMUTATIONS, seed=ANALYSIS_SEED, n_jobs=n_jobs)
        results[c] = _analyze_cohort(
            c, cohort_dir, loaded, params, out_dir=out_dir, final_name=c, n_jobs=n_jobs, findings=findings,
            cell_type=cell_type,
        )
    return results


# ---------------------------------------------------------------------------
# MOTCO rank sensitivity (forced-rank sweep beside the CV-selected analysis)
# ---------------------------------------------------------------------------

# Forced latent ranks, a constant and not a flag: 2 = n_stages − 1 (the lowest
# rank at which `angle` carries orientation information, and the rank MOTCO's
# stage-supervised rule selects at three stages in its operating studies); 3
# the next rank; 4, 6, 9 and 12 the steps of MOTCO's own latent-rank ladder
# (4 = 2(K − 1) is the forced rank with a Type I check at this design). Rank 1
# cannot be forced (MOTCO accepts forced_components >= 2 only): the CV-selected
# analysis is the reference row instead. The sweep never chooses a rank.
ANALYSIS_FORCED_RANKS = (2, 3, 4, 6, 9, 12)
RANK_SENSITIVITY_DIR_TEMPLATE = "rank_{:02d}"
RANK_SENSITIVITY_TABLE = "rank_sensitivity.csv"
RANK_SENSITIVITY_MARKDOWN = "rank_sensitivity.md"
RANK_SENSITIVITY_MANIFEST = "rank_sensitivity_manifest.json"
RANK_1_NOTE = (
    "rank 1 cannot be forced (MOTCO accepts forced_components >= 2 only); "
    "the CV-selected analysis is the reference row"
)
_REFERENCE_FILES = ("result.json", "analysis_manifest.json")
# rank_sensitivity.csv holds only values that reproduce across reruns; runtimes
# go to rank_sensitivity.md and the manifest
_RANK_SENSITIVITY_CSV_EXCLUDED = ("runtime_seconds",)


def _csv_float(value) -> str:
    """Shortest round-trip representation, so a CSV cell equals the JSON value exactly."""
    return repr(float(value))


def _check_reference(cohort_dir: str | Path, analysis_dir: str | Path, inputs: dict, motco_version: str) -> dict:
    """Read the CV-selected analysis of one cohort and check it matches the sweep.

    Parameters
    ----------
    cohort_dir:
        The export directory the sweep reads.
    analysis_dir:
        ``<analysis root>/<cohort>/`` as written by ``run_analysis``.
    inputs:
        ``load_analysis_inputs(cohort_dir)["inputs"]``: file → path and sha256
        of the five export files now on disk.
    motco_version:
        The installed MOTCO version.

    Returns
    -------
    dict
        ``dir``, ``record`` (parsed ``result.json``), ``manifest`` (parsed
        ``analysis_manifest.json``) and ``files`` (the two files' paths and
        sha256).

    Raises
    ------
    FileNotFoundError
        The analysis directory or one of its two files is missing (the message
        names the ``analyze`` command).
    ValueError
        An export file's sha256, the MOTCO version, the seed, the permutation
        count or the component selection recorded by the analysis differs from
        the sweep's; the message names the file or setting.
    """
    d = Path(analysis_dir)
    hint = "run `sead-sexdiff analyze` to produce it"
    if not d.is_dir():
        raise FileNotFoundError(f"{d} does not exist; {hint}")
    for name in _REFERENCE_FILES:
        if not (d / name).exists():
            raise FileNotFoundError(f"{d / name} is missing; {hint}")
    record = json.loads((d / "result.json").read_text())
    manifest = json.loads((d / "analysis_manifest.json").read_text())
    rerun = "rerun `sead-sexdiff analyze` so the reference matches the inputs"
    recorded = manifest.get("inputs", {}).get("files", {})
    for name in _ANALYSIS_INPUTS:
        before = recorded.get(name, {}).get("sha256")
        now = inputs[name]["sha256"]
        if before != now:
            raise ValueError(
                f"{Path(cohort_dir) / name} has changed since the analysis in {d} "
                f"(sha256 {before} recorded, {now} now); {rerun}"
            )
    ref_version = manifest.get("software", {}).get("motco")
    if ref_version != motco_version:
        raise ValueError(
            f"{d / 'analysis_manifest.json'} records MOTCO version {ref_version} but {motco_version} is "
            f"installed; {rerun}"
        )
    meas = manifest.get("measurement", {})
    for key, expected in (("seed", ANALYSIS_SEED), ("permutations", ANALYSIS_PERMUTATIONS)):
        if meas.get(key) != expected:
            raise ValueError(
                f"{d / 'analysis_manifest.json'} records {key} = {meas.get(key)} but the sweep uses {expected}; {rerun}"
            )
    selection = record.get("latent", {}).get("component_selection")
    if selection != "cv":
        raise ValueError(
            f"{d / 'result.json'} records component_selection = {selection!r}; the reference must be the "
            f"CV-selected analysis; {rerun}"
        )
    files = {name: {"path": str(d / name), "sha256": _sha256_file(d / name)} for name in _REFERENCE_FILES}
    return {"dir": d, "record": record, "manifest": manifest, "files": files}


def _sensitivity_row(record: dict, *, cohort: str, component_selection: str, directory: str) -> dict:
    """Flatten one analysis record into a row of the cross-rank table.

    ``directory`` is the reference analysis directory for the ``cv`` row and
    the ``rank_<kk>`` name (relative to the cohort directory) for a forced
    row, so the table does not depend on where the sweep was written.
    ``runtime_seconds`` is rendered in the Markdown and recorded in the
    manifest but left out of the CSV.
    """
    row = {"cohort": cohort, "component_selection": component_selection, "rank": record["latent"]["selected_lv"]}
    for s in ANALYSIS_STATISTICS:
        p = record["p_values"].get(s)
        null = record["null_summary"].get(s, {})
        row[s] = record["statistics"].get(s)
        row[f"{s}_p"] = p
        row[f"{s}_rejected"] = None if p is None else bool(p <= ANALYSIS_ALPHA)
        row[f"{s}_null_mean"] = null.get("mean")
        row[f"{s}_null_sd"] = null.get("sd")
        row[f"{s}_null_q95"] = null.get("q95")
    row["cv_mean_auroc"] = record["latent"].get("cv_mean_auroc")
    row["runtime_seconds"] = record["measurement"]["runtime_seconds"]
    row["dir"] = directory
    return row


def render_rank_sensitivity(
    rows: list[dict],
    *,
    cohort: str,
    seed: int = ANALYSIS_SEED,
    permutations: int = ANALYSIS_PERMUTATIONS,
    alpha: float = ANALYSIS_ALPHA,
) -> str:
    """Render ``rank_sensitivity.md`` from the rows of ``rank_sensitivity.csv``.

    Parameters
    ----------
    rows:
        ``_sensitivity_row`` dicts, the reference (``cv``) row first and the
        forced ranks in ascending order.
    cohort:
        Cohort name for the title.
    seed, permutations, alpha:
        Named in the header; a statistic is rejected when ``p <= alpha``.

    Returns
    -------
    str
        Markdown: a header naming the ladder, seed, permutations, α and the
        reference rank, one table per statistic (rank × observed / p /
        rejected / null mean / sd / q95) and a table of the runs (CV AUROC,
        runtime, directory). Every number comes from the rows; the runtimes
        are the one thing here that is not in ``rank_sensitivity.csv``.
    """
    ref = rows[0]
    ranks = [r["rank"] for r in rows[1:]]
    lines = [
        f"# MOTCO rank sensitivity: {cohort}",
        "",
        f"Forced-rank sweep of the covariate-adjusted trajectory analysis. The reference row (`cv`) is the "
        f"CV-selected analysis at rank {ref['rank']}, which stays the primary result; the forced ranks "
        f"{', '.join(str(k) for k in ranks)} are read as sensitivity to it, and no rank is chosen from the "
        f"results. {permutations} RRPP permutations with seed {seed} at every rank (the same permutation "
        f"indices throughout); α = {alpha} (rejected when p ≤ α). `angle` and `shape` non-rejections are "
        "uninformative at every rank; a `delta` rejection is a difference, not specifically of size.",
    ]
    for s in ANALYSIS_STATISTICS:
        lines += [
            "",
            f"## {s}" + (" (degrees)" if s == "angle" else ""),
            "",
            "| selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |",
            "|---|---:|---:|---:|:---:|---:|---:|---:|",
        ]
        for r in rows:
            rejected = r[f"{s}_rejected"]
            lines.append(
                f"| {r['component_selection']} | {r['rank']} | {_fmt(r[s])} | {_fmt(r[f'{s}_p'], '.4g')} | "
                f"{'n/a' if rejected is None else ('yes' if rejected else 'no')} | "
                f"{_fmt(r[f'{s}_null_mean'])} | {_fmt(r[f'{s}_null_sd'])} | {_fmt(r[f'{s}_null_q95'])} |"
            )
    lines += [
        "",
        "## Runs",
        "",
        "| selection | rank | cv_mean_auroc | runtime (s) | directory |",
        "|---|---:|---:|---:|---|",
    ]
    for r in rows:
        lines.append(
            f"| {r['component_selection']} | {r['rank']} | {_fmt(r['cv_mean_auroc'], '.4f')} | "
            f"{_fmt(r['runtime_seconds'], '.1f')} | `{r['dir']}` |"
        )
    return "\n".join(lines) + "\n"


def run_rank_sensitivity(
    cohort: str = "all",
    *,
    input_dir: str | Path = "results/motco_input",
    analysis_dir: str | Path = "results/motco_analysis",
    out_dir: str | Path = "results/motco_rank_sensitivity",
    n_jobs: int = 1,
    force: bool = False,
) -> dict:
    """Rerun the covariate-adjusted analysis at every forced rank of the ladder.

    Parameters
    ----------
    cohort:
        ``"cohort_30"``, ``"cohort_50"`` or ``"all"``.
    input_dir:
        Export root; each cohort is read from ``<input_dir>/<cohort>/``.
    analysis_dir:
        Analysis root of the ``analyze`` step; ``<analysis_dir>/<cohort>/`` is
        the reference, read and checked, never written.
    out_dir:
        Sweep root; each cohort goes to ``<out_dir>/<cohort>/``.
    n_jobs:
        Parallel jobs for MOTCO's RRPP (``-1`` = all CPUs); as in
        ``run_analysis``, values other than 1 change the permutation draws and
        are recorded as a finding. The reported results use 1.
    force:
        Replace existing cohort directories.

    Returns
    -------
    dict
        Cohort → ``{"dir", "reference", "ranks" (rank → _analyze_cohort
        result), "table" (DataFrame), "manifest", "findings"}``.

    Raises
    ------
    FileExistsError
        A cohort directory exists and ``force`` is false.
    FileNotFoundError
        The export or the analysis of a cohort is missing.
    ValueError
        The reference analysis does not match the inputs or settings, or a rank
        of the ladder is infeasible for a cohort.

    Notes
    -----
    Every check (existing output, export and analysis present, reference
    hashes and settings, feasibility of every rank) runs for every requested
    cohort before anything is written. Each rank is one ``_analyze_cohort``
    call with ``TrajectoryAnalysisParams(permutations=ANALYSIS_PERMUTATIONS,
    seed=ANALYSIS_SEED, n_jobs=n_jobs, forced_components=k)``, so the only
    difference from the ``analyze`` step is the rank; cross-validation does
    not run. A cohort is written under ``<cohort>.tmp/`` (``rank_<kk>/``
    directories, ``rank_sensitivity.csv``, ``rank_sensitivity.md``,
    ``rank_sensitivity_manifest.json``) and renamed when complete.
    """
    import shutil
    from importlib.metadata import PackageNotFoundError, version

    motco_version = check_motco_version()
    from motco.analysis import TrajectoryAnalysisParams

    input_dir = Path(input_dir)
    analysis_dir = Path(analysis_dir)
    out_dir = Path(out_dir)
    ladder = tuple(ANALYSIS_FORCED_RANKS)
    cohorts = list(EXPORT_COHORTS) if cohort == "all" else [cohort]
    for c in cohorts:
        if c not in EXPORT_COHORTS:
            raise ValueError(f"Unknown cohort {c!r}; expected one of {EXPORT_COHORTS} or 'all'")
    for c in cohorts:
        if (out_dir / c).exists() and not force:
            raise FileExistsError(f"{out_dir / c} exists; rerun with --force to replace it")
    for c in cohorts:
        if not (input_dir / c).is_dir():
            raise FileNotFoundError(f"{input_dir / c} does not exist; run `sead-sexdiff export` to produce it")
        if not (analysis_dir / c).is_dir():
            raise FileNotFoundError(
                f"{analysis_dir / c} does not exist; run `sead-sexdiff analyze` to produce the reference analysis"
            )
    print(f"[rank-sensitivity] MOTCO {motco_version}")
    print(f"[rank-sensitivity] Input directory: {input_dir}")
    print(f"[rank-sensitivity] Reference analysis directory: {analysis_dir}")
    print(f"[rank-sensitivity] Output directory: {out_dir}")
    print(f"[rank-sensitivity] Forced ranks: {', '.join(str(k) for k in ladder)} ({RANK_1_NOTE})")

    # Load and check every cohort before writing anything
    prepared = {}
    for c in cohorts:
        cohort_dir = input_dir / c
        loaded = load_analysis_inputs(cohort_dir)
        reference = _check_reference(cohort_dir, analysis_dir / c, loaded["inputs"], motco_version)
        n_donors = len(loaded["metadata"])
        n_features = sum(b.shape[1] for b in loaded["blocks"].values())
        max_rank = min(n_donors, n_features)
        for k in ladder:
            if not (2 <= k <= max_rank):
                raise ValueError(
                    f"forced rank {k} is infeasible for {c}: MOTCO accepts 2 <= forced_components <= "
                    f"min(n_donors = {n_donors}, n_features = {n_features}) = {max_rank}; the ladder "
                    "ANALYSIS_FORCED_RANKS is a constant, nothing was written"
                )
        ref_rank = reference["record"]["latent"]["selected_lv"]
        print(
            f"[rank-sensitivity] {c}: reference {reference['dir']} (CV-selected rank {ref_rank}, "
            f"hashes and settings match)"
        )
        prepared[c] = (loaded, reference)

    results = {}
    for c in cohorts:
        loaded, reference = prepared[c]
        cohort_dir = input_dir / c
        final = out_dir / c
        ref_record = reference["record"]
        ref_rank = ref_record["latent"]["selected_lv"]
        tmp = out_dir / f"{c}.tmp"
        if tmp.exists():
            shutil.rmtree(tmp)
        tmp.mkdir(parents=True)

        runs = {}
        for k in ladder:
            name = RANK_SENSITIVITY_DIR_TEMPLATE.format(k)
            params = TrajectoryAnalysisParams(
                permutations=ANALYSIS_PERMUTATIONS, seed=ANALYSIS_SEED, n_jobs=n_jobs, forced_components=k
            )
            runs[k] = _analyze_cohort(
                c, cohort_dir, loaded, params,
                out_dir=tmp, final_name=name, published_dir=final / name, n_jobs=n_jobs,
                findings=[
                    f"forced rank {k} (dimensionality sensitivity): the CV-selected analysis in "
                    f"{reference['dir']} is primary"
                ],
                log_prefix="rank-sensitivity",
            )
            rec = runs[k]["record"]
            print(
                f"[rank-sensitivity] {c}: rank {k}: "
                + ", ".join(f"{s} {_fmt(rec['statistics'].get(s))} (p {_fmt(rec['p_values'].get(s), '.4g')})"
                            for s in ANALYSIS_STATISTICS)
            )

        rows = [_sensitivity_row(ref_record, cohort=c, component_selection="cv", directory=str(reference["dir"]))]
        rows += [
            _sensitivity_row(runs[k]["record"], cohort=c, component_selection="forced",
                             directory=RANK_SENSITIVITY_DIR_TEMPLATE.format(k))
            for k in ladder
        ]
        table = pd.DataFrame(rows).drop(columns=list(_RANK_SENSITIVITY_CSV_EXCLUDED))
        table.to_csv(tmp / RANK_SENSITIVITY_TABLE, index=False, float_format=_csv_float)
        (tmp / RANK_SENSITIVITY_MARKDOWN).write_text(
            render_rank_sensitivity(rows, cohort=c, seed=ANALYSIS_SEED, permutations=ANALYSIS_PERMUTATIONS)
        )

        findings = []
        at_cv_rank = None
        if ref_rank in ladder:
            forced_rec = runs[ref_rank]["record"]
            stats_equal = all(
                np.isclose(forced_rec["statistics"][s], ref_record["statistics"][s], rtol=1e-9, atol=0.0)
                for s in ANALYSIS_STATISTICS
            )
            p_equal = forced_rec["p_values"] == ref_record["p_values"]
            at_cv_rank = {"rank": ref_rank, "statistics_equal": bool(stats_equal), "p_values_equal": bool(p_equal)}
            if not (stats_equal and p_equal):
                findings.append(
                    f"the forced run at the CV-selected rank {ref_rank} differs from the reference analysis "
                    f"(statistics equal: {stats_equal}; p-values equal: {p_equal}); the same data, rank and "
                    "seed should give the same record"
                )
        elif ref_rank == 1:
            findings.append(f"CV-selected rank 1 is not in the ladder {ladder}: {RANK_1_NOTE}")
        else:
            findings.append(
                f"CV-selected rank {ref_rank} is not in the ladder {ladder}; the reference row stands in for it"
            )
        if n_jobs != 1:
            findings.append(
                f"n_jobs = {n_jobs}: MOTCO's parallel RRPP seeds one stream per worker, so the permutation "
                "draws, p-values and null summaries differ from an n_jobs = 1 run at every rank; the "
                "case-study results are reported at n_jobs = 1"
            )
        for k in ladder:
            for s in ANALYSIS_STATISTICS:
                p_ref, p_k = ref_record["p_values"].get(s), runs[k]["record"]["p_values"].get(s)
                if p_ref is not None and p_k is not None and (p_ref <= ANALYSIS_ALPHA) != (p_k <= ANALYSIS_ALPHA):
                    findings.append(
                        f"{s}: rejection status at forced rank {k} (p = {_fmt(p_k, '.4g')}) differs from the "
                        f"CV-selected analysis (p = {_fmt(p_ref, '.4g')}); rank dependence of the sensitivity "
                        "analysis, the CV-selected result stays primary"
                    )

        software = _software_versions()
        software["motco"] = motco_version
        try:
            software["scikit-learn"] = version("scikit-learn")
        except PackageNotFoundError:
            software["scikit-learn"] = None
        manifest = {
            "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "cohort": c,
            "n_donors": int(len(loaded["metadata"])),
            "forced_ranks": list(ladder),
            "rank_1_note": RANK_1_NOTE,
            "seed": ANALYSIS_SEED,
            "permutations": ANALYSIS_PERMUTATIONS,
            "same_seed_note": (
                "the same seed and permutation count at every rank, so every rank draws the same permutation "
                "indices and between-rank differences are measurement differences"
            ),
            "alpha": ANALYSIS_ALPHA,
            "rejection_rule": ANALYSIS_REJECTION_RULE,
            "n_jobs": n_jobs,
            "n_jobs_note": "RRPP draws depend on n_jobs (one seed stream per worker); rerun with the same n_jobs",
            "inputs": {"dir": str(cohort_dir), "files": loaded["inputs"]},
            "reference": {
                "dir": str(reference["dir"]),
                "files": reference["files"],
                "component_selection": ref_record["latent"]["component_selection"],
                "selected_lv": ref_rank,
                "cv_mean_auroc": ref_record["latent"].get("cv_mean_auroc"),
                "motco": reference["manifest"]["software"]["motco"],
                "seed": reference["manifest"]["measurement"]["seed"],
                "permutations": reference["manifest"]["measurement"]["permutations"],
                "in_ladder": ref_rank in ladder,
                "forced_run_at_cv_rank": at_cv_rank,
                "note": "read and checked (export hashes, MOTCO version, seed, permutations), never rerun or modified",
            },
            "ranks": {
                str(k): {
                    "dir": str(final / RANK_SENSITIVITY_DIR_TEMPLATE.format(k)),
                    "selected_lv": runs[k]["record"]["latent"]["selected_lv"],
                    "component_selection": runs[k]["record"]["latent"]["component_selection"],
                    "runtime_seconds": runs[k]["record"]["measurement"]["runtime_seconds"],
                }
                for k in ladder
            },
            "outputs": sorted(
                [RANK_SENSITIVITY_DIR_TEMPLATE.format(k) + "/" for k in ladder]
                + [RANK_SENSITIVITY_TABLE, RANK_SENSITIVITY_MARKDOWN, RANK_SENSITIVITY_MANIFEST]
            ),
            "software": software,
            "git": _git_revision(),
            "findings": findings,
        }
        (tmp / RANK_SENSITIVITY_MANIFEST).write_text(json.dumps(manifest, indent=2) + "\n")
        if final.exists():
            shutil.rmtree(final)
        os.replace(tmp, final)
        print(f"[rank-sensitivity] Wrote {final}")
        for msg in findings:
            print(f"[rank-sensitivity] FINDING: {c}: {msg}")
        results[c] = {
            "dir": final, "reference": reference, "ranks": runs, "table": table,
            "manifest": manifest, "findings": findings,
        }
    return results


# ---------------------------------------------------------------------------
# MOTCO preprocessing sensitivity (export variants beside the primary analysis)
# ---------------------------------------------------------------------------

# Registered export variants, a constant and not a flag. Each variant is the
# primary export with the listed `_export_cohort` keyword(s) changed, followed
# by the primary analysis unchanged, and is read as sensitivity to the primary
# analysis only. `keep_chrx` keeps chrX (chrY and chrM are still removed): the
# whole export pipeline is recomputed with chrX present, so library sizes, the
# abundance filter, TMM factors, log2-CPM and the variance ranking all change
# and chrX features compete for the N slots like any other feature. Variants
# are never added, dropped or altered from the results; the primary analysis
# stays primary.
PREPROCESSING_VARIANTS: dict[str, dict] = {
    "keep_chrx": {
        "description": "chrX kept; chrY and chrM removed; every other export setting equal to the primary",
        "export": {"excluded_chromosomes": ("chrY", "chrM")},
    },
    "exclude_families": {
        "description": (
            "RNA genes of the ribosomal, mitochondrial and hemoglobin sets (TECHNICAL_GENE_SETS) removed before the "
            "abundance filter; the ATAC block and every other export setting equal to the primary"
        ),
        "export": {"excluded_gene_sets": ("ribosomal", "mitochondrial", "hemoglobin")},
    },
    "rna_trend_ratio": {
        "description": (
            "RNA features ranked at the selection step by the ratio of their log2-CPM variance to a binned-median "
            "mean–variance trend instead of by raw variance; the ATAC block and every other export setting equal to "
            "the primary"
        ),
        "export": {"selection_rule": {"rna": "trend_ratio", "atac": "variance"}},
    },
    "atac_mean_accessibility": {
        "description": (
            "ATAC peaks ranked at the selection step by mean log2-CPM (total accessibility) instead of by variance; "
            "the RNA block and every other export setting equal to the primary"
        ),
        "export": {"selection_rule": {"rna": "variance", "atac": "mean"}},
    },
    "residualized": {
        "description": (
            "after log2-CPM and before the selection step, each abundance-passing feature replaced by its mean plus "
            "the residual of its regression on an intercept and RESIDUALIZATION_COVARIATES (the export's log10 "
            "library sizes and the technical table's gene-set fractions); the selection ranks the residual variance; "
            "every other export setting equal to the primary"
        ),
        "export": {"residualize": True},
    },
    "revised_candidate": {
        "description": (
            "the four settings together in the primary's order: gene sets removed before the abundance filter, "
            "residualization after log2-CPM, trend-ratio ranking for RNA and mean-accessibility ranking for ATAC at "
            "the selection step; the adoption rule ADOPTION_RULE is evaluated for this variant"
        ),
        "export": {
            "excluded_gene_sets": ("ribosomal", "mitochondrial", "hemoglobin"),
            "residualize": True,
            "selection_rule": {"rna": "trend_ratio", "atac": "mean"},
        },
        "adoption_check": True,
    },
    "trend_normalized": {
        "description": (
            "after log2-CPM and before the selection step, each donor's log2-CPM over the abundance-passing features "
            "corrected by limma's fast cyclic loess (a lowess curve of its deviation from the per-feature mean over "
            "donors against that mean, three iterations, adaptive span) in both omics; every other export setting "
            "equal to the primary; the adoption rule ADOPTION_RULE is evaluated for this variant"
        ),
        "export": {"normalization": "cyclic_loess"},
        "adoption_check": True,
    },
    "revised_candidate_trend": {
        "description": (
            "the revised candidate's four settings plus the cyclic-loess trend normalization of both omics, in the "
            "primary's order: gene sets removed before the abundance filter, cyclic loess after log2-CPM, "
            "residualization after the normalization, trend-ratio ranking for RNA and mean-accessibility ranking for "
            "ATAC at the selection step; the adoption rule ADOPTION_RULE is evaluated for this variant"
        ),
        "export": {
            "excluded_gene_sets": ("ribosomal", "mitochondrial", "hemoglobin"),
            "normalization": "cyclic_loess",
            "residualize": True,
            "selection_rule": {"rna": "trend_ratio", "atac": "mean"},
        },
        "adoption_check": True,
    },
}
PREPROCESSING_PRIMARY_LABEL = "primary"
# The primary export's value of every registry setting (named in the comparison)
PREPROCESSING_PRIMARY_SETTINGS: dict = {
    "excluded_chromosomes": list(EXPORT_EXCLUDED_CHROMOSOMES),
    "excluded_gene_sets": [],
    "selection_rule": {"rna": "variance", "atac": "variance"},
    "normalization": {"rna": "none", "atac": "none"},
    "residualize": None,
}
# The adoption rule of the revised candidate (design D7): a constant over
# generated numbers, evaluated by the step from the candidate's latent scores
# and written to comparison/<cohort>/adoption_check.csv. The check variables
# are technical variables that are not regressors of the residualization: the
# primary export's TMM factors (the normalization by-product that proxies the
# composition effect), the nucleus QC means of the technical table, RIN and
# PMI. No p-value, rejection status or CV AUROC is read. Never chosen from
# results; flag and ask before changing.
ADOPTION_CHECK_VARIABLES = (
    "tmm_rna",
    "tmm_atac",
    "fraction_mitochondrial_umis_mean",
    "doublet_score_mean",
    "gex_reads_mapped_confidently_to_intronic_regions_mean",
    "gex_fraction_of_transcriptomic_reads_in_cells_mean",
    "rin",
    "pmi",
)
ADOPTION_RULE = (
    "holds when no latent component of the candidate's CV-selected analysis reaches |Spearman| >= "
    "DIAGNOSTIC_CORRELATION_THRESHOLD against any variable of ADOPTION_CHECK_VARIABLES, in both cohorts; "
    "p-values, rejection statuses and CV AUROC are not read. If it holds, a later change may propose the candidate "
    "as a revised primary with the current primary kept as the reference; if not, the variants stay sensitivity analyses"
)
PREPROCESSING_ADOPTION_TABLE = "adoption_check.csv"
# The precedence rule between the adoption candidates of the trend-normalization
# change, written before any run: applied by the step from the ``holds`` values
# of both candidates' adoption tables in both cohorts and written into each
# candidate's comparison manifest (``adoption_precedence``); ``null`` with the
# reason "incomplete" until both candidates have run in both cohorts.
# ``trend_normalized`` wins ties because it changes the least. Never reads a
# p-value, rank or CV AUROC; flag and ask before changing.
ADOPTION_PRECEDENCE_CANDIDATES = ("trend_normalized", "revised_candidate_trend")
ADOPTION_PRECEDENCE = (
    "if trend_normalized meets ADOPTION_RULE in both cohorts it is the candidate a later revised-primary change "
    "proposes, whatever revised_candidate_trend gives (no gene removed, no degrees of freedom spent, no technical "
    "covariate read); if only revised_candidate_trend meets it in both cohorts, that one is; if neither does, the "
    "variants stay sensitivity analyses. Only the holds values of the adoption tables are read"
)
# Layout under <out_dir>/<variant>/: mirrors results/motco_input and
# results/motco_analysis so the rank sweep can read a variant with no code change
PREPROCESSING_INPUT_DIR = "motco_input"
PREPROCESSING_ANALYSIS_DIR = "motco_analysis"
PREPROCESSING_COMPARISON_DIR = "comparison"
PREPROCESSING_SENSITIVITY_TABLE = "sensitivity.csv"
PREPROCESSING_SENSITIVITY_MARKDOWN = "sensitivity.md"
PREPROCESSING_SENSITIVITY_MANIFEST = "sensitivity_manifest.json"
# X-dosage marker genes whose status is reported beside the feature accounting
PREPROCESSING_NAMED_FEATURES = ("XIST",)


def _check_export_reference(
    primary_export_dir: str | Path,
    pseudobulk_dir: str | Path,
    annotation: dict,
    *,
    k: int,
    cells: dict,
    n_donors: int,
) -> dict:
    """Read the primary export of one cohort and check it was made from the inputs now on disk.

    Parameters
    ----------
    primary_export_dir:
        ``<export root>/<cohort>/`` as written by ``run_export``.
    pseudobulk_dir:
        The pseudobulk directory the sensitivity step reads.
    annotation:
        Provenance of the annotation now resolved (``resolve_gene_annotation``).
    k, cells, n_donors:
        The values ``donor_table.csv`` gives now for the cohort
        (``smallest_cell``), which the primary export must have used.

    Returns
    -------
    dict
        ``dir``, ``manifest`` (parsed ``manifest.json``) and ``files`` (the
        five analysis-input files' paths and sha256).

    Raises
    ------
    FileNotFoundError
        The export directory or one of its files is missing (the message
        names the ``export`` command).
    ValueError
        The pseudobulk manifest, the covariates manifest or the annotation now
        on disk has a sha256 different from the one the export recorded; a
        threshold differs from the export constants; or ``k``, the cells or
        the donor count differ from the donor table now; the message names the
        file or setting and the ``export`` command. Nothing is written.
    """
    d = Path(primary_export_dir)
    pb = Path(pseudobulk_dir)
    hint = "run `sead-sexdiff export` to produce the primary export"
    if not d.is_dir():
        raise FileNotFoundError(f"{d} does not exist; {hint}")
    for name in _ANALYSIS_INPUTS:
        if not (d / name).exists():
            raise FileNotFoundError(f"{d / name} is missing; {hint}")
    manifest_path = d / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    rerun = (
        "the primary export does not match its inputs or the export constants; rerun `sead-sexdiff export` "
        "(then `sead-sexdiff analyze`) before the preprocessing sensitivity analysis"
    )
    recorded_input = manifest.get("input", {})
    for name, key in (
        ("manifest.json", "pseudobulk_manifest_sha256"),
        ("covariates_manifest.json", "covariates_manifest_sha256"),
    ):
        before = recorded_input.get(key)
        path = pb / name
        now = _sha256_file(path) if path.exists() else None
        if before != now:
            raise ValueError(
                f"{path} has sha256 {now} but {manifest_path} recorded {before} ({key}); {rerun}"
            )
    before = manifest.get("annotation", {}).get("sha256")
    if before != annotation.get("sha256"):
        raise ValueError(
            f"annotation file {annotation.get('file')} has sha256 {annotation.get('sha256')} but "
            f"{manifest_path} recorded {before}; {rerun}"
        )
    thr = manifest.get("thresholds", {})
    recorded = {
        "excluded_chromosomes": thr.get("excluded_chromosomes"),
        "min_cpm": thr.get("min_cpm"),
        "prior_count": thr.get("prior_count"),
        "N": thr.get("N"),
        "tmm.logratio_trim": thr.get("tmm", {}).get("logratio_trim"),
        "tmm.sum_trim": thr.get("tmm", {}).get("sum_trim"),
        "float_format": thr.get("float_format"),
        "k": thr.get("k"),
        "cells": manifest.get("metadata", {}).get("cells"),
        "n_donors": manifest.get("n_donors"),
    }
    expected = {
        "excluded_chromosomes": list(EXPORT_EXCLUDED_CHROMOSOMES),
        "min_cpm": EXPORT_MIN_CPM,
        "prior_count": EXPORT_PRIOR_COUNT,
        "N": EXPORT_TOP_N,
        "tmm.logratio_trim": TMM_LOGRATIO_TRIM,
        "tmm.sum_trim": TMM_SUM_TRIM,
        "float_format": EXPORT_FLOAT_FORMAT,
        "k": k,
        "cells": cells,
        "n_donors": n_donors,
    }
    for key, want in expected.items():
        if recorded[key] != want:
            source = "the donor table now gives" if key in ("k", "cells", "n_donors") else "the export constants give"
            raise ValueError(f"{manifest_path} records {key} = {recorded[key]!r} but {source} {want!r}; {rerun}")
    files = {name: {"path": str(d / name), "sha256": _sha256_file(d / name)} for name in _ANALYSIS_INPUTS}
    return {"dir": d, "manifest": manifest, "files": files}


def _variant_needs_technical(spec: dict) -> bool:
    """Whether a registry entry reads the technical table (gene sets excluded, residualized, or the adoption check)."""
    export = spec.get("export", {})
    return bool(export.get("excluded_gene_sets")) or bool(export.get("residualize")) or bool(spec.get("adoption_check"))


def _technical_columns_needed(spec: dict) -> list[str]:
    """The technical-table columns a registry entry needs complete for every cohort donor."""
    columns: list[str] = []
    if spec.get("export", {}).get("residualize"):
        columns += [c for c in RESIDUALIZATION_COVARIATES if not c.startswith("log10_lib_")]
    if spec.get("adoption_check"):
        columns += [c for c in ADOPTION_CHECK_VARIABLES if not c.startswith("tmm_") and c not in ("rin", "pmi")]
    return list(dict.fromkeys(columns))


def _check_technical_reference(pseudobulk_dir: str | Path, *, donors: Iterable[str], columns: Iterable[str] = ()) -> dict:
    """Read the technical table and check it was made from the pseudobulk now on disk and is complete for the donors.

    Parameters
    ----------
    pseudobulk_dir:
        Directory holding ``manifest.json`` (the pseudobulk manifest the export
        reads), ``donor_technical.csv`` and ``technical_manifest.json``.
    donors:
        The cohort's donors; every one needs a row.
    columns:
        Technical columns that must be non-missing for every donor (the
        residualization fractions, the adoption-check QC means).

    Returns
    -------
    dict
        ``table`` (``donor_technical.csv``), ``manifest``
        (``technical_manifest.json``), ``files`` (both paths and sha256) and
        ``pseudobulk_manifest_sha256`` (the pseudobulk manifest's sha256 now).

    Raises
    ------
    FileNotFoundError
        The table or its manifest is missing (the message names the file and
        the ``technical-covariates`` command).
    ValueError
        The technical manifest records a pseudobulk manifest sha256 other than
        the file's now (both named), its gene sets differ from
        ``TECHNICAL_GENE_SETS``, a donor has no row, or a column is missing
        for a donor (donor and column named). Nothing is written.
    """
    pb = Path(pseudobulk_dir)
    table_path = pb / TECHNICAL_OUTPUTS["table"]
    manifest_path = pb / TECHNICAL_OUTPUTS["manifest"]
    hint = (
        "run `sead-sexdiff technical-covariates --omic both` on the cluster and sync results/pseudobulk/ "
        "(scripts/sync_pseudobulk.sh) to produce it"
    )
    for p in (table_path, manifest_path):
        if not p.exists():
            raise FileNotFoundError(f"{p} is missing; {hint}")
    manifest = json.loads(manifest_path.read_text())
    recorded = manifest.get("pseudobulk", {}).get("manifest_sha256")
    now = _sha256_file(pb / "manifest.json")
    if recorded != now:
        raise ValueError(
            f"{manifest_path} records the pseudobulk manifest sha256 {recorded} but {pb / 'manifest.json'} has sha256 "
            f"{now}: the technical table was made from a different pseudobulk than the export reads; rerun "
            "`sead-sexdiff technical-covariates --omic both` on the cluster and sync so both match; nothing was written"
        )
    recorded_sets = {
        s: list(v.get("patterns", [])) if isinstance(v, dict) else list(v)
        for s, v in (manifest.get("gene_sets") or {}).items()
    }
    expected_sets = {s: list(p) for s, p in TECHNICAL_GENE_SETS.items()}
    if recorded_sets != expected_sets:
        raise ValueError(
            f"{manifest_path} records gene sets {recorded_sets} but TECHNICAL_GENE_SETS is {expected_sets}; the "
            "gene-set exclusion and the technical fractions must use the same sets; rerun `sead-sexdiff "
            "technical-covariates --omic both` and sync; nothing was written"
        )
    table = pd.read_csv(table_path, dtype={"donor_id": str})
    if "donor_id" not in table.columns:
        raise ValueError(f"{table_path} lacks a donor_id column")
    tech = table.set_index("donor_id")
    donors = list(donors)
    absent = sorted(set(donors) - set(tech.index))
    if absent:
        raise ValueError(f"{table_path} has no row for cohort donor(s) {absent}; {hint}")
    for col in columns:
        if col not in tech.columns:
            raise ValueError(f"{table_path} lacks column {col!r}; {hint}")
        values = pd.to_numeric(tech.loc[donors, col], errors="coerce").to_numpy(dtype=float)
        bad = [d for d, v in zip(donors, values) if not np.isfinite(v)]
        if bad:
            raise ValueError(
                f"{table_path}: column {col!r} is missing for cohort donor(s) {bad}; a variant or arm that reads the "
                "technical table needs every value for every cohort donor (nothing is imputed); nothing was written"
            )
    files = {
        p.name: {"path": str(p), "sha256": _sha256_file(p)} for p in (table_path, manifest_path)
    }
    return {"table": table, "manifest": manifest, "files": files, "pseudobulk_manifest_sha256": now}


def _variant_settings_record(spec: dict) -> dict:
    """The registry entry's export settings as the manifests record them.

    Tuples become lists, ``residualize`` the covariates and ``normalization``
    the per-omic rules (as the export manifest records them).
    """
    out = {}
    for key, val in spec["export"].items():
        if key == "residualize":
            out[key] = list(RESIDUALIZATION_COVARIATES) if val else None
        elif key == "normalization":
            out[key] = _per_omic_setting(
                val, default="none", name="normalization", label="normalization", allowed=NORMALIZATION_RULES
            )
        elif isinstance(val, (list, tuple)):
            out[key] = list(val)
        else:
            out[key] = val
    return out


def _variant_row(record: dict, *, variant: str, cohort: str, directory: str) -> dict:
    """``_sensitivity_row`` with a leading ``variant`` column (``primary`` or the variant name)."""
    row = _sensitivity_row(
        record, cohort=cohort, component_selection=record["latent"]["component_selection"], directory=directory
    )
    return {"variant": variant, **row}


def _kept_mean_quantiles(export_dir: Path, omic: str) -> dict[str, float]:
    """Quantiles (``DIAGNOSTIC_QUANTILES``) of the kept features' mean log2-CPM, from the block CSV."""
    block = pd.read_csv(export_dir / f"{omic}.csv", index_col=0)
    means = block.mean(axis=0).to_numpy(dtype=float)
    return {f"q{int(round(q * 100)):02d}": float(np.quantile(means, q)) for q in DIAGNOSTIC_QUANTILES}


def _feature_accounting(primary_export_dir: str | Path, variant_export_dir: str | Path) -> dict:
    """Compare the feature filtering of the primary and a variant export.

    Parameters
    ----------
    primary_export_dir, variant_export_dir:
        Export directories of one cohort (``features_<omic>.csv``,
        ``manifest.json`` and, for the kept-mean quantiles, ``<omic>.csv``).

    Returns
    -------
    dict
        Omic → ``feature_counts`` and ``removed_by_chromosome`` of both exports,
        ``selection_rule`` of both, ``kept`` (both kept counts; shared,
        primary-only and variant-only) and, when the variant's settings call
        for them: ``chrx`` and ``named`` (when the excluded chromosomes differ:
        chrX input count, the variant's status counts, kept chrX count and
        names, and the status of ``PREPROCESSING_NAMED_FEATURES`` in both
        exports), ``gene_sets`` (sets excluded: matched, removed and
        already-on-an-excluded-chromosome counts per set), ``selection``
        (rules differ: rule and kept-mean quantiles of both exports),
        ``residualization`` (covariates, degrees of freedom and R² quantiles)
        and ``normalization`` (normalization rules differ: both rules, the
        variant's span, iterations, reference, tilt quantiles over donors
        before and after, the largest absolute correction and the trend file).

    Raises
    ------
    ValueError
        The two exports list different input features; a chrX feature of a
        variant that keeps chrX carries ``sex_chromosome``; the gene-set or
        residualization records disagree with ``features_<omic>.csv``; or the
        normalization record is missing or disagrees with
        ``normalization_trend.csv``.
    """
    p, v = Path(primary_export_dir), Path(variant_export_dir)
    pm = json.loads((p / "manifest.json").read_text())
    vm = json.loads((v / "manifest.json").read_text())
    p_excluded = list(pm.get("thresholds", {}).get("excluded_chromosomes", []))
    v_excluded = list(vm.get("thresholds", {}).get("excluded_chromosomes", []))
    chromosomes_differ = p_excluded != v_excluded
    dtypes = {"feature": str, "chromosome": str, "status": str}
    out = {}
    for o in ANALYSIS_BLOCKS:
        pf = pd.read_csv(p / f"features_{o}.csv", dtype=dtypes)
        vf = pd.read_csv(v / f"features_{o}.csv", dtype=dtypes)
        if pf["feature"].tolist() != vf["feature"].tolist():
            raise ValueError(f"{p / f'features_{o}.csv'} and {v / f'features_{o}.csv'} list different input features")
        p_kept = set(pf.loc[pf["status"] == "kept", "feature"])
        v_kept = set(vf.loc[vf["status"] == "kept", "feature"])
        p_rule = (pm.get("selection_rule") or {}).get(o, "variance")
        v_rule = (vm.get("selection_rule") or {}).get(o, "variance")
        entry = {
            "feature_counts": {"primary": pm["feature_counts"][o], "variant": vm["feature_counts"][o]},
            "removed_by_chromosome": {
                "primary": pm.get("removed_by_chromosome", {}).get(o),
                "variant": vm.get("removed_by_chromosome", {}).get(o),
            },
            "selection_rule": {"primary": p_rule, "variant": v_rule},
            "kept": {
                "primary": len(p_kept),
                "variant": len(v_kept),
                "shared": len(p_kept & v_kept),
                "primary_only": len(p_kept - v_kept),
                "variant_only": len(v_kept - p_kept),
            },
        }
        if chromosomes_differ:
            v_chrx = vf[vf["chromosome"] == "chrX"]
            if "chrX" not in v_excluded:
                bad = v_chrx.loc[v_chrx["status"] == "sex_chromosome", "feature"].tolist()
                if bad:
                    raise ValueError(
                        f"{v / f'features_{o}.csv'}: {len(bad)} chrX feature(s) carry status sex_chromosome "
                        f"(e.g. {bad[0]!r}) although the variant keeps chrX"
                    )
            kept_chrx = sorted(v_chrx.loc[v_chrx["status"] == "kept", "feature"])
            named = {}
            for name in PREPROCESSING_NAMED_FEATURES:
                rows_p, rows_v = pf[pf["feature"] == name], vf[vf["feature"] == name]
                if len(rows_p) == 1 and len(rows_v) == 1:
                    named[name] = {
                        "chromosome": rows_v["chromosome"].iloc[0],
                        "primary": rows_p["status"].iloc[0],
                        "variant": rows_v["status"].iloc[0],
                    }
            entry["chrx"] = {
                "input": int(len(v_chrx)),
                "status_counts_variant": {s: int(n) for s, n in sorted(v_chrx["status"].value_counts().items())},
                "kept": len(kept_chrx),
                "kept_features": kept_chrx,
            }
            entry["named"] = named
        if vm.get("excluded_gene_sets"):
            removed = (vm.get("removed_by_gene_set") or {}).get(o, {})
            n_status = int((vf["status"] == "excluded_gene_set").sum())
            n_removed = sum(int(info.get("n_removed", 0)) for info in removed.values())
            if n_status != n_removed:
                raise ValueError(
                    f"{v / f'features_{o}.csv'} has {n_status} features with status excluded_gene_set but the manifest "
                    f"records {n_removed} removed by gene set"
                )
            entry["gene_sets"] = {
                "excluded": list(vm["excluded_gene_sets"]),
                "removed": {
                    s: {
                        "n_matched": int(info.get("n_matched", 0)),
                        "n_removed": int(info.get("n_removed", 0)),
                        "n_on_excluded_chromosome": len(info.get("on_excluded_chromosome", [])),
                    }
                    for s, info in removed.items()
                },
            }
        if p_rule != v_rule:
            entry["selection"] = {
                "primary": {"rule": p_rule, "kept_mean_quantiles": _kept_mean_quantiles(p, o)},
                "variant": {"rule": v_rule, "kept_mean_quantiles": _kept_mean_quantiles(v, o)},
            }
        if vm.get("residualization"):
            resid = vm["residualization"]
            if "r_squared" not in vf.columns or vf.loc[vf["status"] == "kept", "r_squared"].isna().any():
                raise ValueError(f"{v / f'features_{o}.csv'} lacks an R² for a kept feature although the export residualized")
            entry["residualization"] = {
                "covariates": list(resid.get("covariates", [])),
                "degrees_of_freedom": resid.get("degrees_of_freedom"),
                "design_file": resid.get("design_file"),
                "design_sha256": resid.get("design_sha256"),
                "r_squared_quantiles": (resid.get("r_squared_quantiles") or {}).get(o),
            }
        p_norm = (pm.get("normalization") or {}).get(o, "none")
        v_norm = (vm.get("normalization") or {}).get(o, "none")
        if p_norm != v_norm:
            fit_block = vm.get("normalization_fit") or {}
            fit = fit_block.get(o)
            if fit is None:
                raise ValueError(f"{v / 'manifest.json'} records normalization {v_norm!r} for {o} but no normalization_fit for it")
            trend_path = v / fit_block.get("trend_file", NORMALIZATION_TREND_FILE)
            if not trend_path.exists():
                raise ValueError(f"{trend_path} is missing although the export trend-normalized {o}")
            trend = pd.read_csv(trend_path, dtype={"donor_id": str, "omic": str})
            trend = trend[trend["omic"] == o]
            if len(trend) != vm.get("n_donors"):
                raise ValueError(f"{trend_path} has {len(trend)} rows for {o}, the export has {vm.get('n_donors')} donors")
            for stat in NORMALIZATION_TILT_STATISTICS:
                for when in ("before", "after"):
                    recorded = fit["tilt_quantiles"][stat][when]
                    if any(not np.isclose(q, recorded[key], rtol=0, atol=1e-12) for key, q in _quantiles(trend[f"{stat}_{when}"]).items()):
                        raise ValueError(f"{trend_path}: the {stat}_{when} quantiles of {o} disagree with {v / 'manifest.json'}")
            entry["normalization"] = {
                "primary": p_norm,
                "variant": v_norm,
                "reference": fit.get("reference"),
                "n_features": fit.get("n_features"),
                "span": fit.get("span"),
                "iterations": fit.get("iterations"),
                "robustness_iterations": fit.get("robustness_iterations"),
                "tilt_quantiles": fit.get("tilt_quantiles"),
                "max_abs_correction": fit.get("max_abs_correction"),
                "trend_file": fit_block.get("trend_file"),
                "trend_sha256": fit_block.get("trend_sha256"),
            }
        out[o] = entry
    return out


def _adoption_variables(primary_export_dir: str | Path, technical: pd.DataFrame, donors: list[str]) -> tuple[pd.DataFrame, dict[str, str]]:
    """``ADOPTION_CHECK_VARIABLES`` per donor: the primary export's TMM factors, its RIN and PMI, the technical table's QC means.

    Raises
    ------
    ValueError
        A variable is absent, or missing for a donor (donor and column named);
        no cell is skipped.
    """
    d = Path(primary_export_dir)
    manifest = json.loads((d / "manifest.json").read_text())
    cov = pd.read_csv(d / "covariates.csv", dtype={"donor_id": str, "education_level": str}).set_index("donor_id")
    tech = technical.set_index("donor_id") if "donor_id" in technical.columns else technical
    table = pd.DataFrame(index=pd.Index(donors, name="donor_id"))
    sources: dict[str, str] = {}
    for name in ADOPTION_CHECK_VARIABLES:
        if name.startswith("tmm_"):
            o = name[len("tmm_"):]
            per_donor = manifest.get("tmm", {}).get(o, {})
            absent = [x for x in donors if x not in per_donor]
            if absent:
                raise ValueError(f"{d / 'manifest.json'} tmm[{o!r}] lacks donor(s) {absent}")
            values = np.array([float(per_donor[x]["factor"]) for x in donors])
            sources[name] = "manifest.json (primary export, TMM factor)"
        elif name in ("rin", "pmi"):
            if name not in cov.columns:
                raise ValueError(f"{d / 'covariates.csv'} lacks column {name!r}")
            values = pd.to_numeric(cov.reindex(donors)[name], errors="coerce").to_numpy(dtype=float)
            sources[name] = "covariates.csv (primary export)"
        else:
            if name not in tech.columns:
                raise ValueError(f"{TECHNICAL_OUTPUTS['table']} lacks column {name!r} (ADOPTION_CHECK_VARIABLES)")
            values = pd.to_numeric(tech.reindex(donors)[name], errors="coerce").to_numpy(dtype=float)
            sources[name] = TECHNICAL_OUTPUTS["table"]
        bad = [x for x, val in zip(donors, values) if not np.isfinite(val)]
        if bad:
            raise ValueError(
                f"adoption check variable {name!r} is missing for cohort donor(s) {bad}; every check variable is "
                "needed for every cohort donor (no cell is skipped)"
            )
        table[name] = values
    return table, sources


def adoption_check(
    scores: pd.DataFrame, variables: pd.DataFrame, *, sources: dict[str, str] | None = None, threshold: float | None = None
) -> pd.DataFrame:
    """Evaluate ``ADOPTION_RULE`` for one cohort from the candidate's latent scores and the check variables.

    Parameters
    ----------
    scores:
        ``latent_scores.csv`` indexed by ``donor_id`` (columns ``pls_<k>``).
    variables:
        ``_adoption_variables`` output, same donors.
    sources:
        Variable → source, recorded per row.
    threshold:
        ``DIAGNOSTIC_CORRELATION_THRESHOLD`` by default.

    Returns
    -------
    pd.DataFrame
        One row per component × variable (``component``, ``variable``,
        ``source``, ``n``, ``spearman``, ``flag`` = |rho| ≥ threshold, ``holds``
        empty) and a final row (``component = all``, ``variable = all``) whose
        ``flag`` is whether any cell is flagged and whose ``holds`` is true
        only when none is. Reads no p-value, rejection status or CV AUROC.

    Raises
    ------
    ValueError
        The donors differ, or a correlation is undefined (constant values).
    """
    if threshold is None:
        threshold = DIAGNOSTIC_CORRELATION_THRESHOLD
    sources = sources or {}
    if sorted(scores.index) != sorted(variables.index):
        raise ValueError("adoption check: the latent scores and the check variables list different donors")
    variables = variables.loc[scores.index]
    comps = sorted([col for col in scores.columns if col.startswith("pls_")], key=_component_order)
    rows = []
    for comp in comps:
        s = scores[comp].to_numpy(dtype=float)
        for var in variables.columns:
            rho, n = spearman(s, variables[var].to_numpy(dtype=float))
            if not np.isfinite(rho):
                raise ValueError(f"adoption check: the Spearman correlation of {comp} with {var} is undefined (constant values)")
            rows.append({
                "component": comp, "variable": var, "source": sources.get(var, ""), "n": n,
                "spearman": rho, "flag": bool(abs(rho) >= threshold), "holds": "",
            })
    any_flag = any(r["flag"] for r in rows)
    rows.append({
        "component": "all", "variable": "all", "source": "", "n": int(len(scores)),
        "spearman": float("nan"), "flag": any_flag, "holds": not any_flag,
    })
    return pd.DataFrame(rows, columns=["component", "variable", "source", "n", "spearman", "flag", "holds"])


def _read_adoption_table(path: Path) -> tuple[bool, list[str]]:
    """The ``holds`` value of an ``adoption_check.csv`` and its flagged cells (``component × variable (rho)``)."""
    table = pd.read_csv(path, dtype={"component": str, "variable": str, "flag": str, "holds": str}, keep_default_na=False)
    last = table.iloc[-1]
    if last["component"] != "all" or last["holds"].strip().lower() not in ("true", "false"):
        raise ValueError(f"{path}: the last row is not the overall holds row")
    body = table.iloc[:-1]
    flagged = [
        f"{r['component']} × {r['variable']} ({float(r['spearman']):.2f})"
        for _, r in body.iterrows() if r["flag"].strip().lower() == "true"
    ]
    return last["holds"].strip().lower() == "true", flagged


def adoption_precedence(holds: dict[str, dict[str, bool | None]], flagged: dict | None = None) -> dict:
    """Apply ``ADOPTION_PRECEDENCE`` to the candidates' ``holds`` values.

    Parameters
    ----------
    holds:
        Candidate (``ADOPTION_PRECEDENCE_CANDIDATES``) → cohort
        (``EXPORT_COHORTS``) → ``holds`` of its ``adoption_check.csv``, or
        ``None`` (absent) when that candidate × cohort has not run.
    flagged:
        Candidate → cohort → flagged cells, recorded beside the decision
        (never read by it).

    Returns
    -------
    dict
        ``rule``, ``candidates``, ``holds``, ``candidate`` (the named
        candidate or ``None``), ``reason`` (``"incomplete"`` while a
        candidate × cohort is missing, with ``missing`` listing them) and
        ``flagged``. Only the ``holds`` values decide.
    """
    first, second = ADOPTION_PRECEDENCE_CANDIDATES
    table = {v: {c: (holds.get(v) or {}).get(c) for c in EXPORT_COHORTS} for v in ADOPTION_PRECEDENCE_CANDIDATES}
    out = {
        "rule": ADOPTION_PRECEDENCE,
        "candidates": list(ADOPTION_PRECEDENCE_CANDIDATES),
        "holds": table,
        "flagged": flagged or {},
    }
    missing = [f"{v} {c}" for v in ADOPTION_PRECEDENCE_CANDIDATES for c in EXPORT_COHORTS if table[v][c] is None]
    if missing:
        return {**out, "candidate": None, "reason": "incomplete", "missing": missing}
    meets = {v: all(table[v].values()) for v in ADOPTION_PRECEDENCE_CANDIDATES}
    if meets[first]:
        reason = (
            f"{first} meets ADOPTION_RULE in both cohorts; it is the candidate whatever {second} gives "
            f"({second} {'also meets' if meets[second] else 'does not meet'} it), because it changes the least: no gene "
            "removed, no degrees of freedom spent, no technical covariate read"
        )
        return {**out, "candidate": first, "reason": reason}
    if meets[second]:
        reason = f"only {second} meets ADOPTION_RULE in both cohorts ({first} does not), so {second} is the candidate"
        return {**out, "candidate": second, "reason": reason}
    reason = (
        f"neither {first} nor {second} meets ADOPTION_RULE in both cohorts; no revised-primary change follows and the "
        "variants stay sensitivity analyses"
    )
    return {**out, "candidate": None, "reason": reason}


def _update_adoption_precedence(out_dir: Path) -> dict:
    """Evaluate ``ADOPTION_PRECEDENCE`` from the candidates' adoption tables under ``out_dir`` and record it.

    The result is written as ``adoption_precedence`` into every existing
    comparison manifest of the candidates (``<out_dir>/<candidate>/comparison/
    <cohort>/sensitivity_manifest.json``) and returned.
    """
    holds: dict[str, dict[str, bool | None]] = {}
    flagged: dict[str, dict[str, list[str]]] = {}
    for v in ADOPTION_PRECEDENCE_CANDIDATES:
        holds[v], flagged[v] = {}, {}
        for c in EXPORT_COHORTS:
            path = out_dir / v / PREPROCESSING_COMPARISON_DIR / c / PREPROCESSING_ADOPTION_TABLE
            if path.exists():
                holds[v][c], flagged[v][c] = _read_adoption_table(path)
            else:
                holds[v][c] = None
    result = adoption_precedence(holds, flagged)
    for v in ADOPTION_PRECEDENCE_CANDIDATES:
        for c in EXPORT_COHORTS:
            path = out_dir / v / PREPROCESSING_COMPARISON_DIR / c / PREPROCESSING_SENSITIVITY_MANIFEST
            if path.exists():
                manifest = json.loads(path.read_text())
                manifest["adoption_precedence"] = result
                path.write_text(json.dumps(manifest, indent=2) + "\n")
    return result


def render_preprocessing_sensitivity(
    rows: list[dict],
    features: dict,
    *,
    variant: str,
    cohort: str,
    settings: dict,
    description: str = "",
    seed: int = ANALYSIS_SEED,
    permutations: int = ANALYSIS_PERMUTATIONS,
    alpha: float = ANALYSIS_ALPHA,
    findings: list[str] | None = None,
    adoption: pd.DataFrame | None = None,
) -> str:
    """Render ``sensitivity.md`` from the rows of ``sensitivity.csv`` and the feature accounting.

    Parameters
    ----------
    rows:
        ``_variant_row`` dicts, the primary row first and the variant row second.
    features:
        ``_feature_accounting`` output.
    variant, settings, description:
        The variant's name, the export settings it changes (as the manifests
        record them) and its description (from ``PREPROCESSING_VARIANTS``).
    cohort:
        Cohort name for the title.
    seed, permutations, alpha:
        Named in the header; a statistic is rejected when ``p <= alpha``.
    findings:
        The comparison's findings, listed at the end.
    adoption:
        ``adoption_check`` output for the registry entry with the adoption
        check; ``None`` otherwise.

    Returns
    -------
    str
        Markdown: a header naming the variant, its settings beside the
        primary's, seed, permutations, α and both ranks and the reading rules;
        one table per statistic (two rows); the feature table per omic and
        export and, when the variant calls for them, the chrX, gene-set,
        selection-rule, residualization, normalization and adoption-check
        blocks; a run table
        (rank selection, CV AUROC, runtime, directory); the findings. Every
        number comes from the rows, the accounting and the adoption table; the
        runtimes are the one thing here that is not in ``sensitivity.csv``.
    """
    primary, var = rows[0], rows[1]
    changed = ", ".join(f"`{key}` = {val}" for key, val in settings.items())
    primary_values = ", ".join(f"`{key}` = {PREPROCESSING_PRIMARY_SETTINGS.get(key)}" for key in settings)
    has = lambda key: any(key in acc for acc in features.values())  # noqa: E731
    blocks_measured = (
        "neither chrX-kept, family-excluded, trend- or mean-selected, residualized nor trend-normalized"
        if has("normalization") else "neither chrX-kept, family-excluded, trend- or mean-selected nor residualized"
    )
    lines = [
        f"# MOTCO preprocessing sensitivity: {cohort}, variant `{variant}`",
        "",
        f"{description}. Export setting(s) changed: {changed} (primary: {primary_values}); every other export "
        f"setting and every analysis setting equal the primary's. The `{primary['variant']}` row is the CV-selected "
        f"primary analysis at rank {primary['rank']}, which stays the primary result; the `{variant}` row "
        f"(CV-selected rank {var['rank']}) is read as sensitivity to it, never as its replacement, and no variant is "
        f"chosen from the results. {permutations} RRPP permutations with seed {seed} in both rows (the same "
        f"permutation indices); α = {alpha} (rejected when p ≤ α). `angle` and `shape` non-rejections are "
        "uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection "
        "status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result "
        "of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. "
        f"MOTCO's operating characteristics were measured on none of the variant blocks ({blocks_measured}), so the "
        "variant carries no Type I or power claim of its own.",
    ]
    for s in ANALYSIS_STATISTICS:
        lines += [
            "",
            f"## {s}" + (" (degrees)" if s == "angle" else ""),
            "",
            "| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |",
            "|---|---|---:|---:|---:|:---:|---:|---:|---:|",
        ]
        for r in rows:
            rejected = r[f"{s}_rejected"]
            lines.append(
                f"| {r['variant']} | {r['component_selection']} | {r['rank']} | {_fmt(r[s])} | "
                f"{_fmt(r[f'{s}_p'], '.4g')} | {'n/a' if rejected is None else ('yes' if rejected else 'no')} | "
                f"{_fmt(r[f'{s}_null_mean'])} | {_fmt(r[f'{s}_null_sd'])} | {_fmt(r[f'{s}_null_q95'])} |"
            )

    show_chrx, show_sets = has("chrx"), has("gene_sets")
    header = ["omic", "export", "input", "after chromosome removal"]
    if show_sets:
        header.append("after gene-set removal")
    header += ["after abundance", "kept (after selection)"]
    if show_chrx:
        header.append("kept on chrX")
    header += ["kept shared", "kept only here"]
    align = ["---", "---"] + ["---:"] * (len(header) - 2)
    note = "Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports"
    if show_chrx:
        note += "; `kept on chrX` is per export"
    lines += ["", "## Features", "", note + ".", "", "| " + " | ".join(header) + " |", "|" + "|".join(align) + "|"]
    for o, acc in features.items():
        for label, key, only in (
            (primary["variant"], "primary", acc["kept"]["primary_only"]),
            (variant, "variant", acc["kept"]["variant_only"]),
        ):
            fc = acc["feature_counts"][key]
            cells = [o, label, str(fc["input"]), str(fc["after_sex_chromosome"])]
            if show_sets:
                cells.append(str(fc.get("after_gene_set", fc["after_sex_chromosome"])))
            cells += [str(fc["after_abundance"]), str(fc["after_variance"])]
            if show_chrx:
                cells.append(str(acc["chrx"]["kept"] if key == "variant" else 0))
            cells += [str(acc["kept"]["shared"]), str(only)]
            lines.append("| " + " | ".join(cells) + " |")
    named = [(o, name, info) for o, acc in features.items() for name, info in acc.get("named", {}).items()]
    if named:
        lines.append("")
        for o, name, info in named:
            lines.append(
                f"- `{name}` ({info['chromosome']}, {o}): status `{info['primary']}` in the "
                f"{primary['variant']} export, `{info['variant']}` in the {variant} export."
            )

    if show_sets:
        lines += [
            "",
            "## Gene sets removed",
            "",
            "RNA genes matched by the excluded sets of `TECHNICAL_GENE_SETS` are removed at step 2 after the chromosome "
            "removal (status `excluded_gene_set` in `features_rna.csv`); a matched gene on an excluded chromosome keeps "
            "`sex_chromosome` and is counted under `already on an excluded chromosome`. ATAC names match no set.",
            "",
            "| omic | set | matched | removed here | already on an excluded chromosome |",
            "|---|---|---:|---:|---:|",
        ]
        for o, acc in features.items():
            for s, info in acc.get("gene_sets", {}).get("removed", {}).items():
                lines.append(f"| {o} | {s} | {info['n_matched']} | {info['n_removed']} | {info['n_on_excluded_chromosome']} |")
    if has("selection"):
        lines += [
            "",
            "## Selection rule",
            "",
            f"`variance`: population variance (ddof 0) of log2-CPM, the primary's rule. `trend_ratio`: variance over a "
            f"mean–variance trend, the median variance of {TREND_N_BINS} equal-count bins of mean log2-CPM placed at "
            "each bin's median mean and linearly interpolated (flat beyond the end bins; the bins are in the variant "
            "export's manifest). `mean`: mean log2-CPM. Ties by name under every rule; the quantiles are of the kept "
            "features' mean log2-CPM.",
            "",
            "| omic | export | rule | kept mean q05 | q25 | q50 | q75 | q95 |",
            "|---|---|---|---:|---:|---:|---:|---:|",
        ]
        for o, acc in features.items():
            sel = acc.get("selection")
            if not sel:
                continue
            for label, key in ((primary["variant"], "primary"), (variant, "variant")):
                q = sel[key]["kept_mean_quantiles"]
                lines.append(
                    f"| {o} | {label} | {sel[key]['rule']} | " + " | ".join(_fmt(q[k], ".3f") for k in ("q05", "q25", "q50", "q75", "q95")) + " |"
                )
    if has("normalization"):
        norm_omics = [o for o, acc in features.items() if "normalization" in acc]
        first = features[norm_omics[0]]["normalization"]
        lines += [
            "",
            "## Normalization",
            "",
            f"After log2-CPM and before {'the residualization and ' if has('residualization') else ''}the selection step, "
            "each donor's log2-CPM over the abundance-passing features was corrected by limma's fast cyclic loess "
            f"(`normalizeCyclicLoess(method = \"fast\")`): at each of {first['iterations']} iterations the reference is "
            "the per-feature mean over donors (A), and each donor's deviation from it (M) is fitted by a lowess curve on "
            f"A ({first['robustness_iterations']} robustifying iterations, delta = {LOESS_DELTA_FRACTION:g} of the range "
            "of A) and subtracted. The span is limma's adaptive default from the number of features "
            f"(`LOESS_SPAN_RULE`: {LOESS_SPAN_RULE['min_span']} + (1 − {LOESS_SPAN_RULE['min_span']}) · "
            f"({LOESS_SPAN_RULE['small_n']} / n)^(1/3)). Only the block's values are read (no label, covariate or "
            "technical variable). The tilt statistics are per donor, with A and M taken against the per-feature mean "
            "of the values before or after the normalization: the OLS slope of M on A and the median M of the lowest "
            f"and highest of {NORMALIZATION_TILT_BINS} equal-count bins of A; the quantiles are over donors and the "
            f"per-donor values are in `{first['trend_file']}` beside the variant blocks (sha256 "
            f"`{str(first['trend_sha256'])[:12]}…`). Limitations: the normalization assumes that most features do not "
            "change at any abundance level, so a stage or sex effect shaped like a tilt is removed with the artefact; "
            "the curve is least determined at the rare end, where the selection draws most features; the "
            "lncRNA-versus-coding offset at matched abundance is not corrected; and no Type I check exists for the "
            "normalized blocks.",
            "",
            "| omic | primary | variant | features | span | iterations | largest absolute correction |",
            "|---|---|---|---:|---:|---:|---:|",
        ]
        for o in norm_omics:
            n = features[o]["normalization"]
            lines.append(
                f"| {o} | {n['primary']} | {n['variant']} | {n['n_features']} | {_fmt(n['span'], '.6f')} | "
                f"{n['iterations']} | {_fmt(n['max_abs_correction']['max'], '.4f')} |"
            )
        lines += [
            "",
            "| omic | statistic | values | q05 | q25 | q50 | q75 | q95 |",
            "|---|---|---|---:|---:|---:|---:|---:|",
        ]
        labels = {"slope": "slope of M on A", "low_bin": "median M, lowest bin", "high_bin": "median M, highest bin"}
        for o in norm_omics:
            n = features[o]["normalization"]
            for stat in NORMALIZATION_TILT_STATISTICS:
                for when in ("before", "after"):
                    q = n["tilt_quantiles"][stat][when]
                    lines.append(
                        f"| {o} | {labels[stat]} | {when} | " + " | ".join(_fmt(q[k], '.4f') for k in ("q05", "q25", "q50", "q75", "q95")) + " |"
                    )
            q = n["max_abs_correction"]["quantiles"]
            lines.append(
                f"| {o} | largest absolute correction | per donor | " + " | ".join(_fmt(q[k], '.4f') for k in ("q05", "q25", "q50", "q75", "q95")) + " |"
            )
    if has("residualization"):
        resid = next(acc["residualization"] for acc in features.values() if "residualization" in acc)
        covs = resid["covariates"]
        exclusions = ", ".join(f"`{name}` ({reason})" for name, reason in RESIDUALIZATION_CELL_STATE_EXCLUSIONS.items()) or "none"
        lines += [
            "",
            "## Residualization",
            "",
            f"After log2-CPM{', the trend normalization' if has('normalization') else ''} and before the selection step, "
            "each abundance-passing feature was replaced by its mean plus "
            f"the residual of its OLS fit on an intercept and {len(covs)} covariates: {', '.join(f'`{c}`' for c in covs)} "
            f"(`RESIDUALIZATION_COVARIATES`: the export's own log10 library sizes after step 2 and the technical table's "
            f"gene-set pseudobulk fractions; excluded as cell state under the technical-axis report's rule: {exclusions}). "
            f"The design (`{resid['design_file']}` beside the variant blocks, sha256 `{str(resid['design_sha256'])[:12]}…`) "
            "holds no sex, stage or analysis covariate. Limitations: the fit spends "
            f"{resid['degrees_of_freedom']} degrees of freedom per feature (intercept plus covariates) before the "
            "permutation null, which does not see them; it removes whatever part of the stage and sex signal the "
            "covariates carry; the analysis covariates were left unchanged, so the library-size proxies are adjusted "
            "twice (log10 library sizes here, log10 nuclei counts in the RRPP model); and no Type I check exists for "
            "the residualized blocks.",
            "",
            "| omic | R² q05 | q25 | q50 | q75 | q95 |",
            "|---|---:|---:|---:|---:|---:|",
        ]
        for o, acc in features.items():
            q = acc.get("residualization", {}).get("r_squared_quantiles") or {}
            if q:
                lines.append(f"| {o} | " + " | ".join(_fmt(q[k], ".3f") for k in ("q05", "q25", "q50", "q75", "q95")) + " |")
    if adoption is not None:
        last = adoption.iloc[-1]
        holds = bool(last["holds"])
        body = adoption[adoption["component"] != "all"]
        flagged = body[body["flag"]]
        lines += [
            "",
            "## Adoption check",
            "",
            f"`ADOPTION_RULE`: {ADOPTION_RULE}. Threshold {DIAGNOSTIC_CORRELATION_THRESHOLD} (`DIAGNOSTIC_CORRELATION_THRESHOLD`). "
            f"In this cohort the rule {'holds' if holds else 'does not hold'}"
            + ("" if holds else ": flagged " + ", ".join(f"{r['component']} × `{r['variable']}` ({_fmt(r['spearman'], '.2f')})" for _, r in flagged.iterrows()))
            + ". The rule needs both cohorts; whether the candidate becomes a revised primary is a later change with the "
            "current primary kept as the reference, and the candidate's statistics above are sensitivity, not the result."
            + (
                f" Between {' and '.join(f'`{c}`' for c in ADOPTION_PRECEDENCE_CANDIDATES)}, `ADOPTION_PRECEDENCE` decides "
                f"from both cohorts' `holds` values: {ADOPTION_PRECEDENCE}; its outcome is recorded as `adoption_precedence` "
                "in the comparison manifests."
                if variant in ADOPTION_PRECEDENCE_CANDIDATES else ""
            ),
            "",
            "| component | variable | source | n | Spearman | flagged |",
            "|---|---|---|---:|---:|:---:|",
        ]
        for _, r in body.iterrows():
            lines.append(f"| {r['component']} | `{r['variable']}` | {r['source']} | {r['n']} | {_fmt(r['spearman'], '.3f')} | {'yes' if r['flag'] else 'no'} |")
        lines.append(f"| all | all | | {last['n']} | | {'yes' if last['flag'] else 'no'} → holds = {'true' if holds else 'false'} |")

    lines += [
        "",
        "## Runs",
        "",
        "| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |",
        "|---|---|---:|---:|---:|---|",
    ]
    for r in rows:
        lines.append(
            f"| {r['variant']} | {r['component_selection']} | {r['rank']} | {_fmt(r['cv_mean_auroc'], '.4f')} | "
            f"{_fmt(r['runtime_seconds'], '.1f')} | `{r['dir']}` |"
        )
    lines += ["", "## Findings", ""]
    lines += [f"- {msg}" for msg in (findings or [])] or ["None."]
    return "\n".join(lines) + "\n"


def run_preprocessing_sensitivity(
    variant: str = "all",
    cohort: str = "all",
    *,
    pseudobulk_dir: str | Path = "results/pseudobulk",
    input_dir: str | Path = "results/motco_input",
    analysis_dir: str | Path = "results/motco_analysis",
    out_dir: str | Path = "results/motco_sensitivity",
    data_root: str | Path | None = None,
    n_jobs: int = 1,
    force: bool = False,
) -> dict:
    """Rerun the export with one registered setting changed, then the primary analysis unchanged.

    Parameters
    ----------
    variant:
        A key of ``PREPROCESSING_VARIANTS`` or ``"all"``.
    cohort:
        ``"cohort_30"``, ``"cohort_50"`` or ``"all"``.
    pseudobulk_dir:
        Output directory of ``pseudobulk`` and ``covariates`` (the export's inputs).
    input_dir:
        Export root of the ``export`` step; ``<input_dir>/<cohort>/`` is the
        primary export, read and checked, never written.
    analysis_dir:
        Analysis root of the ``analyze`` step; ``<analysis_dir>/<cohort>/`` is
        the primary analysis, read and checked, never written.
    out_dir:
        Sensitivity root; each variant goes to ``<out_dir>/<variant>/`` with
        ``motco_input/<cohort>/`` (the variant export, plus ``residualization.csv``
        and ``normalization_trend.csv`` when those steps apply),
        ``motco_analysis/<cohort>/`` (the variant analysis) and
        ``comparison/<cohort>/`` (``sensitivity.csv``, ``sensitivity.md``,
        ``sensitivity_manifest.json`` and, for an adoption candidate,
        ``adoption_check.csv``).
    data_root:
        Data root holding (or receiving) the gene annotation.
    n_jobs:
        Parallel jobs for MOTCO's CV and RRPP (``-1`` = all CPUs); as in
        ``run_analysis``, values other than 1 change the permutation draws and
        are recorded as a finding. The reported results use 1.
    force:
        Replace a cohort's three directories under the variant. Without it,
        ``variant = "all"`` skips (with a message) every variant × cohort whose
        comparison directory already exists, so ``all`` is resumable; an
        explicitly named variant with existing output is an error, as is a
        partial tree (export or analysis without the comparison).

    Returns
    -------
    dict
        Variant → cohort → ``{"export", "analysis", "comparison", "rows", "table",
        "manifest", "findings"}`` (skipped variant × cohort pairs are absent).

    Raises
    ------
    ValueError
        An unregistered variant; the primary export or analysis does not match
        the inputs, constants or settings; or the variant's donor tables differ
        from the primary's.
    FileExistsError
        One of a cohort's three directories exists and ``force`` is false.
    FileNotFoundError
        A pseudobulk input, the primary export or the primary analysis is missing.

    Notes
    -----
    Every check (registered variant, existing output, inputs present, primary
    export vs the pseudobulk manifest, covariates manifest, annotation,
    thresholds, ``k``, cells and donor count, primary analysis vs the export
    hashes, MOTCO version, seed, permutations and CV selection; for a variant
    that excludes gene sets, residualizes or runs the adoption check, the
    technical table vs the pseudobulk manifest, ``TECHNICAL_GENE_SETS`` and
    completeness for the cohort donors, ``_check_technical_reference``) runs
    for every requested variant and cohort before anything is written. Per
    cohort the variant export is one ``_export_cohort`` call with the
    variant's keyword(s) (the primary pipeline with the registered setting(s)
    changed: library sizes, the abundance filter, TMM, log2-CPM, the optional
    residualization and the selection recomputed), its ``metadata.csv`` and
    ``covariates.csv`` must be byte-identical to the primary's, and the
    variant analysis is one ``_analyze_cohort`` call with the primary's
    parameters (CV rank selection, ``ANALYSIS_PERMUTATIONS``,
    ``ANALYSIS_SEED``, so RRPP draws the same permutation indices as the
    primary). The comparison holds the primary row first and the variant row
    second; a rejection status that differs is preprocessing dependence of the
    sensitivity analysis and the primary analysis stays primary. For the
    registry entries with ``adoption_check`` the comparison also holds
    ``adoption_check.csv`` (``ADOPTION_RULE`` over ``ADOPTION_CHECK_VARIABLES``);
    after each cohort of an ``ADOPTION_PRECEDENCE_CANDIDATES`` entry,
    ``ADOPTION_PRECEDENCE`` is evaluated from the candidates' adoption tables
    under ``out_dir`` and written as ``adoption_precedence`` into every
    existing candidate comparison manifest (``_update_adoption_precedence``).
    """
    import shutil
    from importlib.metadata import PackageNotFoundError, version

    log = "preprocessing-sensitivity"
    pb = Path(pseudobulk_dir)
    input_dir = Path(input_dir)
    analysis_dir = Path(analysis_dir)
    out_dir = Path(out_dir)
    variants = list(PREPROCESSING_VARIANTS) if variant == "all" else [variant]
    for v in variants:
        if v not in PREPROCESSING_VARIANTS:
            raise ValueError(
                f"Unknown preprocessing variant {v!r}; the registry PREPROCESSING_VARIANTS holds "
                f"{list(PREPROCESSING_VARIANTS)} (or 'all'). Variants are predeclared constants, never chosen "
                "from results; nothing was written"
            )
    cohorts = list(EXPORT_COHORTS) if cohort == "all" else [cohort]
    for c in cohorts:
        if c not in EXPORT_COHORTS:
            raise ValueError(f"Unknown cohort {c!r}; expected one of {EXPORT_COHORTS} or 'all'")
    subdirs = (PREPROCESSING_INPUT_DIR, PREPROCESSING_ANALYSIS_DIR, PREPROCESSING_COMPARISON_DIR)
    skipped: list[tuple[str, str]] = []
    for v in variants:
        for c in cohorts:
            existing = [sub for sub in subdirs if (out_dir / v / sub / c).exists()]
            if not existing or force:
                continue
            if variant == "all" and PREPROCESSING_COMPARISON_DIR in existing:
                skipped.append((v, c))
                continue
            raise FileExistsError(
                f"{out_dir / v / existing[0] / c} exists; rerun with --force to replace the variant's export, "
                f"analysis and comparison of {c}"
            )
    todo = {v: [c for c in cohorts if (v, c) not in skipped] for v in variants}
    for v, c in skipped:
        print(f"[{log}] {v} {c}: comparison directory {out_dir / v / PREPROCESSING_COMPARISON_DIR / c} exists; "
              "skipped (rerun with --force to replace it)")
    variants = [v for v in variants if todo[v]]
    if not variants:
        print(f"[{log}] Nothing to run: every requested variant and cohort is complete under {out_dir}")
        return {}
    cohorts = [c for c in cohorts if any(c in todo[v] for v in variants)]
    check_export_inputs(pb, cohorts)
    for c in cohorts:
        if not (input_dir / c).is_dir():
            raise FileNotFoundError(f"{input_dir / c} does not exist; run `sead-sexdiff export` to produce the primary export")
        if not (analysis_dir / c).is_dir():
            raise FileNotFoundError(
                f"{analysis_dir / c} does not exist; run `sead-sexdiff analyze` to produce the primary analysis"
            )
    motco_version = check_motco_version()
    from motco.analysis import TrajectoryAnalysisParams

    print(f"[{log}] MOTCO {motco_version}")
    print(f"[{log}] Pseudobulk directory: {pb}")
    print(f"[{log}] Primary export directory: {input_dir}")
    print(f"[{log}] Primary analysis directory: {analysis_dir}")
    print(f"[{log}] Output directory: {out_dir}")
    print(f"[{log}] Variants: {', '.join(variants)} (registry PREPROCESSING_VARIANTS; the primary analysis stays primary)")

    inputs = _load_export_inputs(pb, data_root, log_prefix=log)
    annotation = inputs["annotation"]
    # Check every cohort's primary export and analysis (and, when a requested
    # variant reads it, the technical table) before writing anything
    prepared = {}
    for c in cohorts:
        sub = _cohort_rows(inputs["table"], inputs["covariates"], c, pb)
        k, k_cell, cells = smallest_cell(sub)
        ref_export = _check_export_reference(input_dir / c, pb, annotation, k=k, cells=cells, n_donors=len(sub))
        loaded = load_analysis_inputs(input_dir / c)
        reference = _check_reference(input_dir / c, analysis_dir / c, loaded["inputs"], motco_version)
        ref_rank = reference["record"]["latent"]["selected_lv"]
        print(
            f"[{log}] {c}: primary export {ref_export['dir']} ({len(sub)} donors, k = {k} from {k_cell}; pseudobulk, "
            f"covariates, annotation and thresholds match) and primary analysis {reference['dir']} "
            f"(CV-selected rank {ref_rank}; hashes and settings match)"
        )
        needing = [v for v in variants if c in todo[v] and _variant_needs_technical(PREPROCESSING_VARIANTS[v])]
        technical = None
        if needing:
            columns = list(dict.fromkeys(
                col for v in needing for col in _technical_columns_needed(PREPROCESSING_VARIANTS[v])
            ))
            technical = _check_technical_reference(pb, donors=sub["donor_id"].tolist(), columns=columns)
            print(
                f"[{log}] {c}: technical table {pb / TECHNICAL_OUTPUTS['table']} matches the pseudobulk manifest and "
                f"TECHNICAL_GENE_SETS and is complete for {len(sub)} donors × {len(columns)} column(s) "
                f"(needed by {', '.join(needing)})"
            )
        prepared[c] = {
            "sub": sub, "k": k, "cells": cells, "export": ref_export, "analysis": reference, "technical": technical,
        }

    software = _software_versions()
    software["motco"] = motco_version
    try:
        software["scikit-learn"] = version("scikit-learn")
    except PackageNotFoundError:
        software["scikit-learn"] = None

    results: dict[str, dict] = {v: {} for v in variants}
    for v in variants:
        spec = PREPROCESSING_VARIANTS[v]
        v_root = out_dir / v
        recorded_settings = _variant_settings_record(spec)
        for c in todo[v]:
            prep = prepared[c]
            primary_export, primary_analysis = prep["export"], prep["analysis"]
            settings = {"excluded_chromosomes": EXPORT_EXCLUDED_CHROMOSOMES, **spec["export"]}
            if settings.get("residualize"):
                settings["residualize"] = prep["technical"]["table"]
            else:
                settings.pop("residualize", None)
            for sub in subdirs:  # with force, all three are replaced before the first is rewritten
                if (v_root / sub / c).exists():
                    shutil.rmtree(v_root / sub / c)
            variant_block = {
                "name": v,
                "description": spec["description"],
                "settings": recorded_settings,
                "reference_export_dir": str(primary_export["dir"]),
            }
            exp = _export_cohort(
                c, prep["sub"], inputs["blocks"], inputs["covariates"], annotation, inputs["pb_info"],
                out_dir=v_root / PREPROCESSING_INPUT_DIR, final_name=c, variant=variant_block, log_prefix=log,
                **settings,
            )
            for name in ("metadata.csv", "covariates.csv"):
                if _sha256_file(exp["dir"] / name) != primary_export["files"][name]["sha256"]:
                    raise ValueError(
                        f"{exp['dir'] / name} is not byte-identical to {primary_export['dir'] / name}: the variant "
                        "export must use the primary export's donors and covariates, so the pseudobulk inputs are "
                        "not the ones the primary export was made from"
                    )
            loaded_v = load_analysis_inputs(exp["dir"])
            params = TrajectoryAnalysisParams(permutations=ANALYSIS_PERMUTATIONS, seed=ANALYSIS_SEED, n_jobs=n_jobs)
            ana = _analyze_cohort(
                c, exp["dir"], loaded_v, params,
                out_dir=v_root / PREPROCESSING_ANALYSIS_DIR, final_name=c, n_jobs=n_jobs,
                findings=[f"preprocessing variant {v}: the primary analysis in {primary_analysis['dir']} is primary"],
                log_prefix=log,
            )

            # Comparison: the primary row first, the variant row second
            ref_record, var_record = primary_analysis["record"], ana["record"]
            rows = [
                _variant_row(ref_record, variant=PREPROCESSING_PRIMARY_LABEL, cohort=c, directory=str(primary_analysis["dir"])),
                _variant_row(var_record, variant=v, cohort=c, directory=f"{PREPROCESSING_ANALYSIS_DIR}/{c}"),
            ]
            for r in rows:
                print(
                    f"[{log}] {c}: {r['variant']} (rank {r['rank']}): "
                    + ", ".join(f"{s} {_fmt(r[s])} (p {_fmt(r[f'{s}_p'], '.4g')})" for s in ANALYSIS_STATISTICS)
                )
            features = _feature_accounting(primary_export["dir"], exp["dir"])
            for o, acc in features.items():
                extra = ""
                if "chrx" in acc:
                    extra += f"; {acc['chrx']['kept']} of {acc['chrx']['input']} chrX features kept by {v}"
                if "gene_sets" in acc:
                    extra += "; removed by gene set " + ", ".join(f"{s} {i['n_removed']}" for s, i in acc["gene_sets"]["removed"].items())
                if "selection" in acc:
                    extra += (f"; rule {acc['selection']['variant']['rule']} (median kept mean "
                              f"{acc['selection']['variant']['kept_mean_quantiles']['q50']:.3f} vs "
                              f"{acc['selection']['primary']['kept_mean_quantiles']['q50']:.3f})")
                if "residualization" in acc and acc["residualization"].get("r_squared_quantiles"):
                    extra += f"; median R² {acc['residualization']['r_squared_quantiles']['q50']:.3f}"
                if "normalization" in acc:
                    tq = acc["normalization"]["tilt_quantiles"]["slope"]
                    extra += (f"; {acc['normalization']['variant']} span {acc['normalization']['span']:.4f}, median slope of M "
                              f"on A {tq['before']['q50']:.4g} → {tq['after']['q50']:.4g}")
                print(
                    f"[{log}] {c} {o}: kept {acc['kept']['shared']} shared, {acc['kept']['primary_only']} primary-only, "
                    f"{acc['kept']['variant_only']} {v}-only{extra}"
                )
            adoption = None
            if spec.get("adoption_check"):
                scores = pd.read_csv(ana["dir"] / "latent_scores.csv", dtype={"donor_id": str}).set_index("donor_id")
                variables, sources = _adoption_variables(primary_export["dir"], prep["technical"]["table"], scores.index.tolist())
                adoption = adoption_check(scores, variables, sources=sources)
                adoption_holds = bool(adoption.iloc[-1]["holds"])
                flagged = adoption[(adoption["component"] != "all") & adoption["flag"]]
                print(f"[{log}] {c}: adoption check {'holds' if adoption_holds else 'does not hold'} "
                      f"({len(flagged)} of {len(adoption) - 1} cells at |rho| >= {DIAGNOSTIC_CORRELATION_THRESHOLD})")
            findings: list[str] = []
            if adoption is not None:
                if adoption_holds:
                    findings.append(
                        f"adoption check holds in {c}: no candidate component reaches |Spearman| >= "
                        f"{DIAGNOSTIC_CORRELATION_THRESHOLD} against any of {len(ADOPTION_CHECK_VARIABLES)} check variables "
                        "(the rule needs both cohorts; a revised primary is a later change)"
                    )
                else:
                    findings.append(
                        f"adoption check does not hold in {c}: " + ", ".join(
                            f"{r['component']} x {r['variable']} ({_fmt(r['spearman'], '.2f')})" for _, r in flagged.iterrows()
                        ) + "; the variants stay sensitivity analyses"
                    )
            ref_rank, var_rank = ref_record["latent"]["selected_lv"], var_record["latent"]["selected_lv"]
            if ref_rank != var_rank:
                findings.append(
                    f"CV-selected rank differs: {ref_rank} in the primary analysis, {var_rank} in the {v} analysis; "
                    "reported beside angle, not read as better or worse"
                )
            for s in ANALYSIS_STATISTICS:
                p_ref, p_v = ref_record["p_values"].get(s), var_record["p_values"].get(s)
                if p_ref is not None and p_v is not None and (p_ref <= ANALYSIS_ALPHA) != (p_v <= ANALYSIS_ALPHA):
                    findings.append(
                        f"{s}: rejection status in the {v} analysis (p = {_fmt(p_v, '.4g')}) differs from the primary "
                        f"analysis (p = {_fmt(p_ref, '.4g')}); preprocessing dependence of the sensitivity analysis, "
                        "the primary result stays primary"
                    )
            if n_jobs != 1:
                findings.append(
                    f"n_jobs = {n_jobs}: MOTCO's parallel RRPP seeds one stream per worker, so the permutation "
                    "draws, p-values and null summaries differ from an n_jobs = 1 run (statistics, rank, latent "
                    "scores and loadings do not); the case-study results are reported at n_jobs = 1"
                )
            findings += [f"{v} export: {msg}" for msg in exp["findings"]]

            comp_root = v_root / PREPROCESSING_COMPARISON_DIR
            tmp = comp_root / f"{c}.tmp"
            if tmp.exists():
                shutil.rmtree(tmp)
            tmp.mkdir(parents=True)
            table = pd.DataFrame(rows).drop(columns=list(_RANK_SENSITIVITY_CSV_EXCLUDED))
            table.to_csv(tmp / PREPROCESSING_SENSITIVITY_TABLE, index=False, float_format=_csv_float)
            if adoption is not None:
                adoption.to_csv(tmp / PREPROCESSING_ADOPTION_TABLE, index=False, float_format=_csv_float)
            (tmp / PREPROCESSING_SENSITIVITY_MARKDOWN).write_text(
                render_preprocessing_sensitivity(
                    rows, features, variant=v, cohort=c, settings=recorded_settings,
                    description=spec["description"], seed=ANALYSIS_SEED, permutations=ANALYSIS_PERMUTATIONS,
                    alpha=ANALYSIS_ALPHA, findings=findings, adoption=adoption,
                )
            )
            variant_files = {
                name: {"path": str(exp["dir"] / name), "sha256": _sha256_file(exp["dir"] / name)}
                for name in _ANALYSIS_INPUTS
            }
            if (exp["dir"] / RESIDUALIZATION_FILE).exists():
                variant_files[RESIDUALIZATION_FILE] = {
                    "path": str(exp["dir"] / RESIDUALIZATION_FILE), "sha256": _sha256_file(exp["dir"] / RESIDUALIZATION_FILE),
                    "note": "the residualization design; not an analysis input (load_analysis_inputs reads the five files above)",
                }
            if (exp["dir"] / NORMALIZATION_TREND_FILE).exists():
                variant_files[NORMALIZATION_TREND_FILE] = {
                    "path": str(exp["dir"] / NORMALIZATION_TREND_FILE),
                    "sha256": _sha256_file(exp["dir"] / NORMALIZATION_TREND_FILE),
                    "note": "the per-donor tilt accounting of the trend normalization; not an analysis input",
                }
            variant_analysis_files = {
                name: {"path": str(ana["dir"] / name), "sha256": _sha256_file(ana["dir"] / name)}
                for name in _REFERENCE_FILES
            }
            manifest = {
                "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "variant": {"name": v, "description": spec["description"], "settings": recorded_settings},
                "cohort": c,
                "n_donors": int(len(prep["sub"])),
                "seed": ANALYSIS_SEED,
                "permutations": ANALYSIS_PERMUTATIONS,
                "same_seed_note": (
                    "the same seed and permutation count as the primary analysis, so both draw the same "
                    "permutation indices and the difference between the rows is a data difference"
                ),
                "alpha": ANALYSIS_ALPHA,
                "rejection_rule": ANALYSIS_REJECTION_RULE,
                "n_jobs": n_jobs,
                "n_jobs_note": "RRPP draws depend on n_jobs (one seed stream per worker); rerun with the same n_jobs",
                "inputs": {
                    **inputs["pb_info"],
                    "annotation": {key: annotation.get(key) for key in ("release", "file", "sha256")},
                    "technical": None if prep["technical"] is None or not _variant_needs_technical(spec) else {
                        **prep["technical"]["files"],
                        "pseudobulk_manifest_sha256": prep["technical"]["pseudobulk_manifest_sha256"],
                        "note": "checked against manifest.json and TECHNICAL_GENE_SETS, complete for the cohort donors",
                    },
                },
                "primary": {
                    "export": {
                        "dir": str(primary_export["dir"]),
                        "files": primary_export["files"],
                        "manifest_sha256": primary_export["files"]["manifest.json"]["sha256"],
                        "k": prep["k"],
                        "cells": prep["cells"],
                        "thresholds": primary_export["manifest"].get("thresholds"),
                    },
                    "analysis": {
                        "dir": str(primary_analysis["dir"]),
                        "files": primary_analysis["files"],
                        "component_selection": ref_record["latent"]["component_selection"],
                        "selected_lv": ref_rank,
                        "cv_mean_auroc": ref_record["latent"].get("cv_mean_auroc"),
                        "runtime_seconds": ref_record["measurement"]["runtime_seconds"],
                        "motco": primary_analysis["manifest"]["software"]["motco"],
                        "seed": primary_analysis["manifest"]["measurement"]["seed"],
                        "permutations": primary_analysis["manifest"]["measurement"]["permutations"],
                    },
                    "note": (
                        "read and checked (pseudobulk, covariates and annotation hashes, thresholds, k, cells, "
                        "export file hashes, MOTCO version, seed, permutations, CV selection), never rerun or modified"
                    ),
                },
                "variant_run": {
                    "export": {
                        "dir": str(exp["dir"]),
                        "files": variant_files,
                        "manifest_sha256": variant_files["manifest.json"]["sha256"],
                        "thresholds": exp["manifest"]["thresholds"],
                        "excluded_gene_sets": exp["manifest"].get("excluded_gene_sets"),
                        "selection_rule": {
                            o: rule for o, rule in (exp["manifest"].get("selection_rule") or {}).items() if o in ANALYSIS_BLOCKS
                        },
                        "residualization": None if not exp["manifest"].get("residualization") else {
                            key: exp["manifest"]["residualization"].get(key)
                            for key in ("covariates", "cell_state_exclusions", "design_file", "design_sha256", "degrees_of_freedom", "r_squared_quantiles")
                        },
                        "normalization": exp["manifest"].get("normalization"),
                        **({"normalization_fit": {
                            o: {key: fit.get(key) for key in ("rule", "reference", "n_features", "span", "iterations", "robustness_iterations")}
                            for o, fit in exp["manifest"]["normalization_fit"].items() if o in ANALYSIS_BLOCKS
                        }} if exp["manifest"].get("normalization_fit") else {}),
                        "feature_counts": exp["manifest"]["feature_counts"],
                        "findings": exp["findings"],
                    },
                    "analysis": {
                        "dir": str(ana["dir"]),
                        "files": variant_analysis_files,
                        "component_selection": var_record["latent"]["component_selection"],
                        "selected_lv": var_rank,
                        "cv_mean_auroc": var_record["latent"].get("cv_mean_auroc"),
                        "runtime_seconds": var_record["measurement"]["runtime_seconds"],
                        "findings": ana["findings"],
                    },
                },
                "features": features,
                "adoption_check": None if adoption is None else {
                    "rule": ADOPTION_RULE,
                    "threshold": DIAGNOSTIC_CORRELATION_THRESHOLD,
                    "variables": [{"name": name, "source": sources[name]} for name in ADOPTION_CHECK_VARIABLES],
                    "components": sorted({r for r in adoption["component"] if r != "all"}, key=_component_order),
                    "holds": adoption_holds,
                    "n_flagged": int(len(flagged)),
                    "flagged": [
                        {"component": r["component"], "variable": r["variable"], "spearman": float(r["spearman"])}
                        for _, r in flagged.iterrows()
                    ],
                    "file": PREPROCESSING_ADOPTION_TABLE,
                    "sha256": _sha256_file(tmp / PREPROCESSING_ADOPTION_TABLE),
                    "reads": "latent_scores.csv of the candidate analysis, the primary export's manifest (tmm) and "
                             "covariates.csv (rin, pmi), donor_technical.csv; no p-value, rejection status or CV AUROC",
                },
                "outputs": sorted([
                    PREPROCESSING_SENSITIVITY_TABLE, PREPROCESSING_SENSITIVITY_MARKDOWN,
                    PREPROCESSING_SENSITIVITY_MANIFEST,
                ] + ([PREPROCESSING_ADOPTION_TABLE] if adoption is not None else [])),
                "software": software,
                "git": _git_revision(),
                "findings": findings,
            }
            (tmp / PREPROCESSING_SENSITIVITY_MANIFEST).write_text(json.dumps(manifest, indent=2) + "\n")
            final = comp_root / c
            if final.exists():
                shutil.rmtree(final)
            os.replace(tmp, final)
            print(f"[{log}] Wrote {final}")
            for msg in findings:
                print(f"[{log}] FINDING: {v} {c}: {msg}")
            results[v][c] = {
                "export": exp, "analysis": ana, "comparison": final, "rows": rows, "table": table,
                "manifest": manifest, "findings": findings, "adoption": adoption,
            }
            if v in ADOPTION_PRECEDENCE_CANDIDATES:
                precedence = _update_adoption_precedence(out_dir)
                manifest["adoption_precedence"] = precedence
                print(
                    f"[{log}] ADOPTION_PRECEDENCE: candidate {precedence['candidate']} ({precedence['reason']}"
                    + (f"; missing {', '.join(precedence['missing'])}" if precedence.get("missing") else "")
                    + f"), written to the comparison manifests of {', '.join(ADOPTION_PRECEDENCE_CANDIDATES)} under {out_dir}"
                )
    return results


# ---------------------------------------------------------------------------
# MOTCO axis diagnostics (what the latent and block axes correlate with)
# ---------------------------------------------------------------------------

# A reporting flag, not an adoption rule: correlations with |rho| at or above
# it are listed first in diagnostics.md. No p-value is attached to any
# correlation. Adoption rules for preprocessing variants belong to the change
# that runs them and are written there before its runs.
DIAGNOSTIC_CORRELATION_THRESHOLD = 0.3
DIAGNOSTIC_N_BLOCK_PCS = 3
DIAGNOSTIC_QUANTILES = (0.05, 0.25, 0.5, 0.75, 0.95)
DIAGNOSTIC_MEAN_CUTS = (0.0, 1.0)
# Variable sources, a constant and not a flag; the resolved variable list with
# its sources is recorded in the manifest. Never chosen or extended from results.
DIAGNOSTIC_VARIABLES: dict[str, dict] = {
    "technical": {
        "file": "donor_technical.csv",
        "columns": "every numeric column except donor_id",
    },
    "export": {
        "file": "manifest.json (tmm block of the export)",
        "columns": ["tmm_rna", "log10_lib_rna", "tmm_atac", "log10_lib_atac"],
        "derivation": "tmm_<omic> = the TMM factor; log10_lib_<omic> = log10 of the TMM library size (filtered column sums)",
    },
    "covariates": {
        "file": "covariates.csv",
        "columns": list(ANALYSIS_COVARIATES),
        "derivation": {name: f"log10({src})" for name, src in ANALYSIS_LOG10_COLUMNS.items()},
    },
    "design": {
        "file": "metadata.csv",
        "columns": ["stage", "male"],
        "derivation": "stage as the ordinal 0, 1, 2; male = 1 for Male and 0 for Female",
    },
}
DIAGNOSTIC_OUTPUTS: dict[str, str] = {
    "axis": "axis_correlations.csv",
    "block_pcs": "block_pcs.csv",
    "geometry": "geometry_shares.csv",
    "features": "feature_selection.csv",
    "markdown": "diagnostics.md",
    "manifest": "diagnostics_manifest.json",
}
DIAGNOSTIC_PRIMARY_LABEL = "primary"
_DIAGNOSTIC_ANALYSIS_FILES = ("result.json", "analysis_manifest.json", "latent_scores.csv", "stage_means.csv")
_DIAGNOSTIC_GEOMETRY = {
    "net_displacement": "a sex's net displacement (last stage mean minus first stage mean)",
    "sex_difference_by_stage": "the per-stage sex-difference vectors (first group minus second), squares summed over stages",
    "sex_difference_centred": "the per-stage sex-difference vectors minus their mean over stages, squares summed over stages",
}


def spearman(x, y) -> tuple[float, int]:
    """Spearman correlation over the pairwise-complete values (NaN with fewer than three or a constant)."""
    import warnings

    from scipy.stats import spearmanr

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    n = int(ok.sum())
    if n < 3 or np.ptp(x[ok]) == 0 or np.ptp(y[ok]) == 0:
        return float("nan"), n
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rho = spearmanr(x[ok], y[ok]).statistic
    return float(rho), n


def partial_spearman(x, y, z) -> tuple[float, int]:
    """Spearman correlation of the OLS residuals of ``x`` and ``y`` on an intercept and ``z``.

    Parameters
    ----------
    x, y:
        The two variables.
    z:
        Conditioning variables, samples × covariates.

    Returns
    -------
    (float, int)
        The partial correlation and the number of rows complete in all three;
        NaN when fewer than ``z.shape[1] + 3`` rows are complete or when a
        variable lies in the span of the conditioning set (its residual is
        numerically zero, as for one of the covariates themselves).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    if z.ndim == 1:
        z = z[:, None]
    ok = np.isfinite(x) & np.isfinite(y) & np.all(np.isfinite(z), axis=1)
    n = int(ok.sum())
    if n < z.shape[1] + 3:
        return float("nan"), n
    design = np.column_stack([np.ones(n), z[ok]])
    res = []
    for v in (x[ok], y[ok]):
        beta = np.linalg.lstsq(design, v, rcond=None)[0]
        r = v - design @ beta
        if np.std(r) <= 1e-8 * max(np.std(v), np.finfo(float).tiny):
            return float("nan"), n
        res.append(r)
    return spearman(res[0], res[1])[0], n


def block_principal_components(values: pd.DataFrame, *, n_components: int = DIAGNOSTIC_N_BLOCK_PCS) -> dict:
    """Principal components of a z-scored block (population SD).

    Parameters
    ----------
    values:
        Donors × features (one exported block).
    n_components:
        Components returned.

    Returns
    -------
    dict
        ``scores`` (donors × ``pc1..``, the SVD scores ``U S``), ``variance_share``
        per component (``s² / Σ s²`` over all components) and ``n_features``.
        Each component's sign is fixed so that its largest-magnitude loading is
        positive.
    """
    x = np.asarray(values, dtype=float)
    mu = x.mean(axis=0)
    sd = x.std(axis=0, ddof=0)
    sd = np.where(sd > 0, sd, 1.0)
    zed = (x - mu) / sd
    u, s, vt = np.linalg.svd(zed, full_matrices=False)
    k = min(n_components, s.size)
    scores = u[:, :k] * s[:k]
    for j in range(k):
        lead = np.argmax(np.abs(vt[j]))
        if vt[j, lead] < 0:
            scores[:, j] = -scores[:, j]
    share = (s ** 2) / float((s ** 2).sum()) if s.size else np.array([])
    cols = [f"pc{j + 1}" for j in range(k)]
    return {
        "scores": pd.DataFrame(scores, index=values.index, columns=cols),
        "variance_share": {c: float(share[j]) for j, c in enumerate(cols)},
        "n_features": int(x.shape[1]),
    }


def geometry_shares(stage_means: pd.DataFrame, *, groups: list[str], stages: list[str]) -> list[dict]:
    """Share of each trajectory-geometry vector's squared length along each latent component.

    Parameters
    ----------
    stage_means:
        MOTCO's ``stage_means.csv`` (``group``, ``stage``, ``pls_k`` columns).
    groups, stages:
        The group and stage levels in order (``result.json`` ``design``).

    Returns
    -------
    list[dict]
        One row per quantity and component: ``quantity``
        (``net_displacement_<group>``, ``sex_difference_by_stage``,
        ``sex_difference_centred``), ``component``, ``value`` (the signed
        component of a single vector; NaN for the summed quantities),
        ``along`` (the squared length along the component), ``total`` (the
        squared length over all components) and ``share`` = ``along / total``.
    """
    comps = [c for c in stage_means.columns if c.startswith("pls_")]
    comps.sort(key=lambda c: int(c.split("_")[1]))
    sm = stage_means.copy()
    sm["group"] = sm["group"].astype(str)
    sm["stage"] = sm["stage"].astype(str)
    means = {
        (g, st): sm.loc[(sm["group"] == g) & (sm["stage"] == st), comps].to_numpy(dtype=float).ravel()
        for g in groups for st in stages
    }
    for key, v in means.items():
        if v.size != len(comps):
            raise ValueError(f"stage_means.csv lacks the row for group {key[0]!r}, stage {key[1]!r}")
    quantities: dict[str, list[np.ndarray]] = {}
    for g in groups:
        quantities[f"net_displacement_{g}"] = [means[(g, stages[-1])] - means[(g, stages[0])]]
    diffs = [means[(groups[0], st)] - means[(groups[1], st)] for st in stages]
    quantities["sex_difference_by_stage"] = diffs
    centre = np.mean(diffs, axis=0)
    quantities["sex_difference_centred"] = [d - centre for d in diffs]
    rows = []
    for name, vectors in quantities.items():
        total = float(sum((v ** 2).sum() for v in vectors))
        for j, c in enumerate(comps):
            along = float(sum(v[j] ** 2 for v in vectors))
            rows.append({
                "quantity": name,
                "component": c,
                "value": float(vectors[0][j]) if len(vectors) == 1 else float("nan"),
                "along": along,
                "total": total,
                "share": along / total if total > 0 else float("nan"),
            })
    return rows


def feature_selection_accounting(
    export_dir: str | Path, omic: str, first_component: pd.Series,
    *, threshold: float = DIAGNOSTIC_CORRELATION_THRESHOLD,
) -> dict:
    """Account for what the export's feature selection kept in one omic.

    Parameters
    ----------
    export_dir:
        The cohort's export directory (``features_<omic>.csv`` and ``<omic>.csv``).
    omic:
        ``"rna"`` or ``"atac"``.
    first_component:
        The primary analysis's ``pls_0`` scores indexed by ``donor_id``.
    threshold:
        A kept feature is "correlated" when ``|Spearman(feature, pls_0)|`` is
        at least this value.

    Returns
    -------
    dict
        ``omic``, ``n_input``, ``n_abundance_passing`` (features with a
        recorded variance), ``n_kept``; the quantiles of the kept features'
        mean log2-CPM (``mean_q05`` … ``mean_q95``); ``share_mean_below_0``,
        ``share_mean_below_1``; ``variance_share_kept`` (variance of the kept
        features over the variance of all abundance-passing features);
        ``variance_kept_min``, ``variance_dropped_median`` (status
        ``variance``); ``share_correlated_pls_0`` and ``n_correlated_pls_0``.
    """
    d = Path(export_dir)
    feats = pd.read_csv(d / f"features_{omic}.csv", dtype={"feature": str, "chromosome": str, "status": str})
    block = pd.read_csv(d / f"{omic}.csv", index_col=0)
    block.index = block.index.astype(str)
    kept_names = feats.loc[feats["status"] == "kept", "feature"].tolist()
    if sorted(kept_names) != sorted(block.columns.astype(str)):
        raise ValueError(f"{d / f'features_{omic}.csv'} kept features differ from the columns of {d / f'{omic}.csv'}")
    passing = feats["variance"].notna()
    var_kept = feats.loc[feats["status"] == "kept", "variance"].to_numpy(dtype=float)
    var_dropped = feats.loc[feats["status"] == "variance", "variance"].to_numpy(dtype=float)
    means = block.mean(axis=0).to_numpy(dtype=float)
    scores = first_component.reindex(block.index).to_numpy(dtype=float)
    if not np.all(np.isfinite(scores)):
        raise ValueError(f"{d / f'{omic}.csv'} lists donors without a pls_0 score")
    from scipy.stats import rankdata

    rank_scores = rankdata(scores)
    rank_block = np.apply_along_axis(rankdata, 0, block.to_numpy(dtype=float))
    rs = (rank_scores - rank_scores.mean()) / rank_scores.std()
    rb = (rank_block - rank_block.mean(axis=0)) / rank_block.std(axis=0)
    rho = (rb * rs[:, None]).mean(axis=0)
    correlated = int((np.abs(rho) >= threshold).sum())
    out = {
        "omic": omic,
        "n_input": int(len(feats)),
        "n_abundance_passing": int(passing.sum()),
        "n_kept": int(len(kept_names)),
    }
    for q in DIAGNOSTIC_QUANTILES:
        out[f"mean_q{int(round(q * 100)):02d}"] = float(np.quantile(means, q))
    for cut in DIAGNOSTIC_MEAN_CUTS:
        out[f"share_mean_below_{int(cut)}"] = float((means < cut).mean())
    out["variance_share_kept"] = float(var_kept.sum() / feats.loc[passing, "variance"].sum())
    out["variance_kept_min"] = float(var_kept.min())
    out["variance_dropped_median"] = float(np.median(var_dropped)) if var_dropped.size else float("nan")
    out["share_correlated_pls_0"] = correlated / len(kept_names)
    out["n_correlated_pls_0"] = correlated
    return out


def _check_analysis_inputs(export_dir: str | Path, analysis_dir: str | Path) -> dict:
    """Read one analysis directory and check its recorded input hashes against the export files.

    Returns
    -------
    dict
        ``dir``, ``record`` (``result.json``), ``manifest``
        (``analysis_manifest.json``) and ``files`` (path and sha256 of the
        analysis files the diagnostic reads).

    Raises
    ------
    FileNotFoundError
        The directory or one of its files is missing.
    ValueError
        A recorded input hash differs from the sha256 of the export file now
        on disk (the message names the file and both hashes).
    """
    d = Path(analysis_dir)
    e = Path(export_dir)
    if not d.is_dir():
        raise FileNotFoundError(f"{d} does not exist; run `sead-sexdiff analyze` (or `rank-sensitivity`) to produce it")
    for name in _DIAGNOSTIC_ANALYSIS_FILES:
        if not (d / name).exists():
            raise FileNotFoundError(f"{d / name} is missing; run `sead-sexdiff analyze` (or `rank-sensitivity`) to produce it")
    manifest = json.loads((d / "analysis_manifest.json").read_text())
    record = json.loads((d / "result.json").read_text())
    recorded = manifest.get("inputs", {}).get("files", {})
    if not recorded:
        raise ValueError(f"{d / 'analysis_manifest.json'} records no input hashes")
    for name, info in recorded.items():
        path = e / name
        now = _sha256_file(path) if path.exists() else None
        if info.get("sha256") != now:
            raise ValueError(
                f"{path} has changed since the analysis in {d} (sha256 {info.get('sha256')} recorded, {now} now); "
                "the diagnostic reads an export and an analysis made from it; rerun `sead-sexdiff analyze` "
                "(and `rank-sensitivity`) or restore the export"
            )
    files = {name: {"path": str(d / name), "sha256": _sha256_file(d / name)} for name in _DIAGNOSTIC_ANALYSIS_FILES}
    return {"dir": d, "record": record, "manifest": manifest, "files": files}


def _diagnostic_variables(export_dir: Path, technical: pd.DataFrame, donors: list[str], meta: pd.DataFrame, pseudobulk_dir: Path) -> tuple[pd.DataFrame, dict[str, dict]]:
    """Assemble the diagnostic variables (``DIAGNOSTIC_VARIABLES``) for the cohort donors.

    Returns the donors × variables table (sorted by ``donor_id``) and the
    variable → ``{"source", "derivation"}`` record for the manifest.
    """
    export_manifest = json.loads((export_dir / "manifest.json").read_text())
    cov = pd.read_csv(export_dir / "covariates.csv", dtype={"donor_id": str, "education_level": str})
    cov = cov.set_index("donor_id").loc[donors]
    table = pd.DataFrame(index=pd.Index(donors, name="donor_id"))
    sources: dict[str, dict] = {}

    tech = technical.set_index("donor_id")
    missing = sorted(set(donors) - set(tech.index))
    if missing:
        raise ValueError(
            f"{pseudobulk_dir / 'donor_technical.csv'} has no row for cohort donor(s) {missing}; run "
            "`sead-sexdiff technical-covariates` on the cluster and sync"
        )
    tech = tech.loc[donors]
    for col in tech.columns:
        if pd.api.types.is_bool_dtype(tech[col]) or not pd.api.types.is_numeric_dtype(tech[col]):
            continue
        table[col] = tech[col].to_numpy(dtype=float)
        sources[col] = {"source": "donor_technical.csv", "derivation": None}

    tmm = export_manifest.get("tmm", {})
    for o in ANALYSIS_BLOCKS:
        per_donor = tmm.get(o, {})
        absent = [d for d in donors if d not in per_donor]
        if absent:
            raise ValueError(f"{export_dir / 'manifest.json'} tmm[{o!r}] lacks donor(s) {absent}")
        table[f"tmm_{o}"] = [float(per_donor[d]["factor"]) for d in donors]
        table[f"log10_lib_{o}"] = [float(np.log10(per_donor[d]["lib_size"])) for d in donors]
        sources[f"tmm_{o}"] = {"source": "manifest.json (export, tmm)", "derivation": "TMM factor"}
        sources[f"log10_lib_{o}"] = {"source": "manifest.json (export, tmm)", "derivation": "log10(lib_size), the TMM library size"}

    for name in ANALYSIS_COVARIATES:
        if name in ANALYSIS_LOG10_COLUMNS:
            src = ANALYSIS_LOG10_COLUMNS[name]
            values = cov[src].to_numpy(dtype=float)
            if np.any(values <= 0):
                raise ValueError(f"{export_dir / 'covariates.csv'}: {src} must be positive for log10")
            table[name] = np.log10(values)
            sources[name] = {"source": "covariates.csv", "derivation": f"log10({src})"}
        else:
            if name not in cov.columns:
                raise ValueError(f"{export_dir / 'covariates.csv'} lacks column {name!r}")
            table[name] = cov[name].to_numpy(dtype=float)
            sources[name] = {"source": "covariates.csv", "derivation": None}

    m = meta.set_index("donor_id").loc[donors]
    table["stage"] = m["stage"].astype(int).to_numpy(dtype=float)
    table["male"] = (m["sex"].astype(str) == "Male").to_numpy(dtype=float)
    sources["stage"] = {"source": "metadata.csv", "derivation": "ordinal stage 0, 1, 2"}
    sources["male"] = {"source": "metadata.csv", "derivation": "1 for Male, 0 for Female"}
    return table, sources


def _correlation_rows(scores: pd.DataFrame, variables: pd.DataFrame, z: np.ndarray, *, threshold: float) -> list[dict]:
    """Raw and partial Spearman rows of every score column against every variable (variables sorted by name)."""
    rows = []
    for comp in scores.columns:
        s = scores[comp].to_numpy(dtype=float)
        sd = float(np.std(s, ddof=1)) if s.size > 1 else float("nan")
        for var in sorted(variables.columns):
            v = variables[var].to_numpy(dtype=float)
            rho, n = spearman(s, v)
            prho, n_partial = partial_spearman(s, v, z)
            rows.append({
                "component": comp,
                "component_sd": sd,
                "variable": var,
                "n": n,
                "spearman": rho,
                "n_partial": n_partial,
                "partial_spearman": prho,
                "flag_raw": bool(np.isfinite(rho) and abs(rho) >= threshold),
                "flag_partial": bool(np.isfinite(prho) and abs(prho) >= threshold),
            })
    return rows


def _component_order(name: str) -> int:
    return int(name.rsplit("_", 1)[1]) if "_" in name else int(name.lstrip("pc"))


def render_axis_diagnostics(
    tables: dict[str, pd.DataFrame],
    *,
    cohort: str,
    analyses: list[dict],
    sources: dict[str, dict],
    findings: list[str],
    threshold: float = DIAGNOSTIC_CORRELATION_THRESHOLD,
) -> str:
    """Render ``diagnostics.md`` from the four diagnostic tables.

    Parameters
    ----------
    tables:
        ``axis``, ``block_pcs``, ``geometry``, ``features`` as written to the CSVs.
    cohort:
        Cohort name for the title.
    analyses:
        One dict per analysis read, in order: ``label``, ``rank``,
        ``selection``, ``dir``.
    sources:
        Variable → ``{"source", "derivation"}``.
    findings:
        The run's findings, listed at the end.
    threshold:
        The flag threshold named in the header.

    Returns
    -------
    str
        Markdown: a header stating what the tables are and are not (no
        p-values, no choice made here); the flagged correlations first (latent
        components, then block PCs, by |rho|); the primary analysis's full
        correlation table; a per-component summary of every forced-rank
        analysis; the block PC tables; the geometry shares; the
        feature-selection accounting; the findings. Every number comes from
        the tables.
    """
    axis, pcs, geo, feats = tables["axis"], tables["block_pcs"], tables["geometry"], tables["features"]
    primary = next(a for a in analyses if a["label"] == DIAGNOSTIC_PRIMARY_LABEL)
    forced = [a for a in analyses if a["label"] != DIAGNOSTIC_PRIMARY_LABEL]
    n_var = len(sources)
    lines = [
        f"# MOTCO axis diagnostics: {cohort}",
        "",
        f"Descriptive measurement-space diagnostics of the covariate-adjusted trajectory analysis: Spearman "
        f"correlations of every latent component of the primary analysis (CV-selected rank {primary['rank']})"
        + (f" and of the forced-rank analyses ({', '.join(str(a['rank']) for a in forced)})" if forced else "")
        + f" with {n_var} predeclared variables (technical table, export TMM factors and library sizes, the six "
        f"analysis covariates, stage and sex), the partial correlations given the six analysis covariates (Spearman "
        f"correlation of the OLS residuals of both variables on an intercept and the covariates), the same for the "
        f"first {DIAGNOSTIC_N_BLOCK_PCS} principal components of each z-scored block, the share of the trajectory "
        f"geometry along each component, and the feature-selection accounting. Rows with |rho| ≥ {threshold} are "
        "flagged and listed first; the flag is a reporting threshold, not an adoption rule. No p-value is attached to "
        "any correlation, the tables support no inferential claim, and nothing here chooses a preprocessing variant: "
        "an axis correlated with a technical variable is a finding about the measurement space of the primary "
        "analysis, not a reinterpretation of its statistics, and the primary analysis stays primary with its "
        "non-rejection reading rules.",
        "",
        "## Flagged correlations",
        "",
    ]
    flagged = axis[axis["flag_raw"] | axis["flag_partial"]].copy()
    if len(flagged):
        flagged["abs"] = flagged["spearman"].abs()
        flagged = flagged.sort_values(["analysis", "component", "abs"], ascending=[True, True, False], key=lambda s: s.map(_analysis_sort_key) if s.name == "analysis" else s)
        lines += [
            "| analysis | rank | component | sd | variable | n | rho | rho \\| covariates | flag |",
            "|---|---:|---|---:|---|---:|---:|---:|---|",
        ]
        for r in flagged.itertuples(index=False):
            flag = "+".join(f for f, on in (("raw", r.flag_raw), ("partial", r.flag_partial)) if on)
            lines.append(
                f"| {r.analysis} | {r.rank} | {r.component} | {_fmt(r.component_sd, '.4g')} | `{r.variable}` | {r.n} | "
                f"{_fmt(r.spearman, '.3f')} | {_fmt_nan(r.partial_spearman, '.3f')} | {flag} |"
            )
    else:
        lines.append(f"No latent-component correlation reaches |rho| ≥ {threshold}.")
    pflag = pcs[pcs["flag_raw"] | pcs["flag_partial"]].copy()
    lines += ["", "Block principal components:", ""]
    if len(pflag):
        pflag["abs"] = pflag["spearman"].abs()
        pflag = pflag.sort_values(["block", "component", "abs"], ascending=[True, True, False])
        lines += [
            "| block | component | variance share | variable | n | rho | rho \\| covariates | flag |",
            "|---|---|---:|---|---:|---:|---:|---|",
        ]
        for r in pflag.itertuples(index=False):
            flag = "+".join(f for f, on in (("raw", r.flag_raw), ("partial", r.flag_partial)) if on)
            lines.append(
                f"| {r.block} | {r.component} | {_fmt(r.variance_share, '.3f')} | `{r.variable}` | {r.n} | "
                f"{_fmt(r.spearman, '.3f')} | {_fmt_nan(r.partial_spearman, '.3f')} | {flag} |"
            )
    else:
        lines.append(f"No block principal component correlation reaches |rho| ≥ {threshold}.")

    # Primary analysis: variables × components
    prim = axis[axis["analysis"] == DIAGNOSTIC_PRIMARY_LABEL]
    comps = sorted(prim["component"].unique(), key=_component_order)
    lines += [
        "",
        f"## Primary analysis (`{primary['dir']}`, {primary['selection']} rank {primary['rank']})",
        "",
        "Component score SD: " + ", ".join(
            f"{c} {_fmt(prim.loc[prim['component'] == c, 'component_sd'].iloc[0], '.4g')}" for c in comps
        ) + ". Columns per component: Spearman rho, then rho given the six analysis covariates (n/a for a covariate "
        "itself). Variables sorted by name; `n` is the number of donors with the variable.",
        "",
        "| variable | source | n | " + " | ".join(f"{c} rho | {c} rho \\| cov" for c in comps) + " |",
        "|---|---|---:|" + "---:|---:|" * len(comps),
    ]
    for var in sorted(prim["variable"].unique()):
        sub = prim[prim["variable"] == var].set_index("component")
        cells = []
        for c in comps:
            r = sub.loc[c]
            cells.append(f"{_fmt(r['spearman'], '.3f')}{'*' if r['flag_raw'] else ''}")
            cells.append(f"{_fmt_nan(r['partial_spearman'], '.3f')}{'*' if r['flag_partial'] else ''}")
        n = int(sub["n"].iloc[0])
        lines.append(f"| `{var}` | {sources[var]['source']} | {n} | " + " | ".join(cells) + " |")
    lines.append("")
    lines.append(f"`*` marks |rho| ≥ {threshold}.")

    # Forced-rank analyses: per-component summary
    if forced:
        lines += [
            "",
            "## Forced-rank analyses",
            "",
            "Per component: score SD, the number of variables flagged on the raw correlation, the variable with "
            "the largest |rho| and that rho, the share of each sex's net displacement and of the per-stage "
            "sex-difference vectors along the component. The full rows are in `axis_correlations.csv` and "
            "`geometry_shares.csv`.",
            "",
        ]
        groups = sorted({q.replace("net_displacement_", "") for q in geo["quantity"].unique() if q.startswith("net_displacement_")})
        lines += [
            "| analysis | rank | component | sd | flagged | top variable | rho | " + " | ".join(f"share net {g}" for g in groups) + " | share sex diff | share sex diff centred |",
            "|---|---:|---|---:|---:|---|---:|" + "---:|" * (len(groups) + 2),
        ]
        for a in forced:
            sub = axis[axis["analysis"] == a["label"]]
            g = geo[geo["analysis"] == a["label"]]
            for c in sorted(sub["component"].unique(), key=_component_order):
                rows = sub[sub["component"] == c]
                top = rows.loc[rows["spearman"].abs().idxmax()]
                shares = [
                    _fmt(g.loc[(g["component"] == c) & (g["quantity"] == f"net_displacement_{grp}"), "share"].iloc[0], ".3f")
                    for grp in groups
                ]
                s1 = _fmt(g.loc[(g["component"] == c) & (g["quantity"] == "sex_difference_by_stage"), "share"].iloc[0], ".3f")
                s2 = _fmt(g.loc[(g["component"] == c) & (g["quantity"] == "sex_difference_centred"), "share"].iloc[0], ".3f")
                lines.append(
                    f"| {a['label']} | {a['rank']} | {c} | {_fmt(top['component_sd'], '.4g')} | {int(rows['flag_raw'].sum())} | "
                    f"`{top['variable']}` | {_fmt(top['spearman'], '.3f')} | " + " | ".join(shares) + f" | {s1} | {s2} |"
                )

    # Block principal components
    lines += ["", "## Block principal components", ""]
    for block in sorted(pcs["block"].unique()):
        sub = pcs[pcs["block"] == block]
        bcomps = sorted(sub["component"].unique(), key=_component_order)
        shares = ", ".join(f"{c} {_fmt(sub.loc[sub['component'] == c, 'variance_share'].iloc[0], '.3f')}" for c in bcomps)
        lines += [
            f"### {block} (variance share: {shares})",
            "",
            "| variable | n | " + " | ".join(f"{c} rho | {c} rho \\| cov" for c in bcomps) + " |",
            "|---|---:|" + "---:|---:|" * len(bcomps),
        ]
        for var in sorted(sub["variable"].unique()):
            vs = sub[sub["variable"] == var].set_index("component")
            cells = []
            for c in bcomps:
                r = vs.loc[c]
                cells.append(f"{_fmt(r['spearman'], '.3f')}{'*' if r['flag_raw'] else ''}")
                cells.append(f"{_fmt_nan(r['partial_spearman'], '.3f')}{'*' if r['flag_partial'] else ''}")
            lines.append(f"| `{var}` | {int(vs['n'].iloc[0])} | " + " | ".join(cells) + " |")
        lines.append("")

    # Geometry shares
    lines += [
        "## Geometry shares",
        "",
        "Share of the squared length of each vector of the trajectory geometry (from `stage_means.csv`) that lies "
        "along each component; `value` is the signed component of a single vector (the net displacements), empty "
        "for the summed quantities. At rank 1 every share is 1 by construction.",
        "",
        "| analysis | rank | quantity | component | value | along | total | share |",
        "|---|---:|---|---|---:|---:|---:|---:|",
    ]
    for r in geo.itertuples(index=False):
        lines.append(
            f"| {r.analysis} | {r.rank} | {r.quantity} | {r.component} | {_fmt_nan(r.value, '.4g')} | "
            f"{_fmt(r.along, '.4g')} | {_fmt(r.total, '.4g')} | {_fmt(r.share, '.3f')} |"
        )

    # Feature selection
    qcols = [f"mean_q{int(round(q * 100)):02d}" for q in DIAGNOSTIC_QUANTILES]
    lines += [
        "",
        "## Feature selection",
        "",
        "Per omic, from the export's `features_<omic>.csv` and block: input features, features passing the "
        "abundance filter, kept features; quantiles of the kept features' mean log2-CPM; the share of kept features "
        "with mean below 0 and below 1; the share of the abundance-passing variance held by the kept features; the "
        "variance at the cut (kept minimum, dropped median); the share of kept features whose |Spearman| with the "
        f"primary analysis's `pls_0` is ≥ {threshold}.",
        "",
        "| omic | input | abundance-passing | kept | " + " | ".join(q.replace("mean_q", "mean q") for q in qcols)
        + " | share mean < 0 | share mean < 1 | variance share kept | kept min variance | dropped median variance | share \\|rho(pls_0)\\| ≥ thr |",
        "|---|---:|---:|---:|" + "---:|" * (len(qcols) + 6),
    ]
    for r in feats.to_dict("records"):
        lines.append(
            f"| {r['omic']} | {r['n_input']} | {r['n_abundance_passing']} | {r['n_kept']} | "
            + " | ".join(_fmt(r[q], ".3f") for q in qcols)
            + f" | {_fmt(r['share_mean_below_0'], '.3f')} | {_fmt(r['share_mean_below_1'], '.3f')} | "
            f"{_fmt(r['variance_share_kept'], '.3f')} | {_fmt(r['variance_kept_min'], '.4g')} | "
            f"{_fmt_nan(r['variance_dropped_median'], '.4g')} | {_fmt(r['share_correlated_pls_0'], '.3f')} ({r['n_correlated_pls_0']}) |"
        )
    lines += ["", "## Analyses read", "", "| analysis | selection | rank | directory |", "|---|---|---:|---|"]
    for a in analyses:
        lines.append(f"| {a['label']} | {a['selection']} | {a['rank']} | `{a['dir']}` |")
    lines += ["", "## Findings", ""]
    lines += [f"- {msg}" for msg in findings] or ["None."]
    return "\n".join(lines) + "\n"


def _fmt_nan(value, spec: str = ".6g") -> str:
    """``_fmt`` that also renders NaN as ``n/a``."""
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return "n/a"
    return format(value, spec)


def _analysis_sort_key(label: str) -> tuple[int, int]:
    if label == DIAGNOSTIC_PRIMARY_LABEL:
        return (0, 0)
    m = re.search(r"(\d+)$", label)
    return (1, int(m.group(1)) if m else 0)


def run_axis_diagnostics(
    cohort: str = "all",
    *,
    cell_type: str = DEFAULT_CELL_TYPE,
    input_dir: str | Path | None = None,
    analysis_dir: str | Path | None = None,
    rank_dir: str | Path | None = None,
    pseudobulk_dir: str | Path | None = None,
    out_dir: str | Path | None = None,
    force: bool = False,
) -> dict:
    """Measure what the MOTCO latent axes and block axes correlate with, per cohort.

    Every directory defaults to the cell type's (``cell_type_dirs``), which
    for ``astrocyte`` are the ``results/`` directories.

    Parameters
    ----------
    cohort:
        ``"cohort_30"``, ``"cohort_50"`` or ``"all"``.
    cell_type:
        Registered cell type (default ``astrocyte``); every export manifest
        read must record it.
    input_dir:
        Export root; ``<input_dir>/<cohort>/`` is read (blocks, features,
        metadata, covariates, manifest).
    analysis_dir:
        Analysis root of the ``analyze`` step; ``<analysis_dir>/<cohort>/`` is
        the primary analysis, read and checked, never written.
    rank_dir:
        Sweep root of ``rank-sensitivity``; every ``<rank_dir>/<cohort>/rank_*/``
        found is read and checked; an absent cohort directory is a finding.
    pseudobulk_dir:
        Directory holding ``donor_technical.csv`` (and ``technical_manifest.json``).
    out_dir:
        Diagnostics root; each cohort goes to ``<out_dir>/<cohort>/``.
    force:
        Replace existing cohort directories.

    Returns
    -------
    dict
        Cohort → ``{"dir", "tables", "analyses", "findings", "manifest"}``.

    Raises
    ------
    FileExistsError
        A cohort directory exists and ``force`` is false.
    FileNotFoundError
        A required input is missing (the message names the file and the
        command that produces it).
    ValueError
        An analysis manifest's input hashes differ from the export files
        (``_check_analysis_inputs``), or the tables disagree.

    Notes
    -----
    Reads only kept outputs and writes nothing into the directories it reads.
    Its outputs are read by no export or analysis command: preprocessing
    variants are predeclared in a change proposal, never derived from a
    diagnostic value. Every check runs for every requested cohort before
    anything is written; each cohort is written to ``<cohort>.tmp/`` and
    renamed. The CSVs are byte-identical across reruns on the same inputs.
    """
    import shutil

    log = "axis-diagnostics"
    cell_type_entry(cell_type)
    input_dir = _resolve_dir(input_dir, cell_type, "motco_input")
    analysis_dir = _resolve_dir(analysis_dir, cell_type, "motco_analysis")
    rank_dir = _resolve_dir(rank_dir, cell_type, "motco_rank_sensitivity")
    pb = _resolve_dir(pseudobulk_dir, cell_type, "pseudobulk")
    out_dir = _resolve_dir(out_dir, cell_type, "motco_diagnostics")
    cohorts = list(EXPORT_COHORTS) if cohort == "all" else [cohort]
    for c in cohorts:
        if c not in EXPORT_COHORTS:
            raise ValueError(f"Unknown cohort {c!r}; expected one of {EXPORT_COHORTS} or 'all'")
    for c in cohorts:
        if (out_dir / c).exists() and not force:
            raise FileExistsError(f"{out_dir / c} exists; rerun with --force to replace it")
    technical_path = pb / TECHNICAL_OUTPUTS["table"]
    if not technical_path.exists():
        raise FileNotFoundError(
            f"{technical_path} is missing; run `sead-sexdiff technical-covariates` on the cluster and sync "
            "results/pseudobulk/ to produce it"
        )
    for c in cohorts:
        if not (input_dir / c).is_dir():
            raise FileNotFoundError(f"{input_dir / c} does not exist; run `sead-sexdiff export` to produce it")
        for name in _ANALYSIS_INPUTS + tuple(f"features_{o}.csv" for o in ANALYSIS_BLOCKS):
            if not (input_dir / c / name).exists():
                raise FileNotFoundError(f"{input_dir / c / name} is missing; run `sead-sexdiff export` to produce it")
        if not (analysis_dir / c).is_dir():
            raise FileNotFoundError(f"{analysis_dir / c} does not exist; run `sead-sexdiff analyze` to produce it")
    print(f"[{log}] Cell type: {cell_type}")
    print(f"[{log}] Export directory: {input_dir}")
    print(f"[{log}] Analysis directory: {analysis_dir}")
    print(f"[{log}] Rank-sensitivity directory: {rank_dir}")
    print(f"[{log}] Technical table: {technical_path}")
    print(f"[{log}] Output directory: {out_dir}")

    technical = pd.read_csv(technical_path, dtype={"donor_id": str})
    technical_manifest_path = pb / TECHNICAL_OUTPUTS["manifest"]
    technical_manifest = json.loads(technical_manifest_path.read_text()) if technical_manifest_path.exists() else None

    # Check every cohort before writing anything
    prepared = {}
    for c in cohorts:
        export_dir = input_dir / c
        rank_findings = check_manifest_cell_type(
            json.loads((export_dir / "manifest.json").read_text()), export_dir / "manifest.json", cell_type
        )
        analyses = [(DIAGNOSTIC_PRIMARY_LABEL, _check_analysis_inputs(export_dir, analysis_dir / c))]
        if (rank_dir / c).is_dir():
            rank_dirs = sorted(p for p in (rank_dir / c).glob("rank_*") if p.is_dir())
            for p in rank_dirs:
                analyses.append((p.name, _check_analysis_inputs(export_dir, p)))
            if not rank_dirs:
                rank_findings.append(f"no forced-rank analyses were read: {rank_dir / c} holds no rank_* directory")
        else:
            rank_findings.append(f"no forced-rank analyses were read: {rank_dir / c} does not exist")
        for label, a in analyses:
            selection = a["record"]["latent"]["component_selection"]
            if label == DIAGNOSTIC_PRIMARY_LABEL and selection != "cv":
                raise ValueError(f"{a['dir'] / 'result.json'} records component_selection = {selection!r}; the primary analysis must be CV-selected")
        print(f"[{log}] {c}: {len(analyses)} analysis directories read and checked against the export hashes "
              f"({', '.join(label for label, _ in analyses)})")
        prepared[c] = (analyses, rank_findings)

    results = {}
    for c in cohorts:
        analyses, findings = prepared[c]
        findings = list(findings)
        export_dir = input_dir / c
        meta = pd.read_csv(export_dir / "metadata.csv", dtype={"donor_id": str, "sex": str, "stage": str})
        donors = sorted(meta["donor_id"].tolist())
        variables, sources = _diagnostic_variables(export_dir, technical, donors, meta, pb)
        z = variables[list(ANALYSIS_COVARIATES)].to_numpy(dtype=float)
        n_missing = {col: int(variables[col].isna().sum()) for col in variables.columns if variables[col].isna().any()}
        if n_missing:
            findings.append(f"variables with missing donors (pairwise-complete correlations): {n_missing}")
        primary_record = analyses[0][1]["record"]
        groups = [str(g) for g in primary_record["design"]["group_levels"]]
        stages = [str(s) for s in primary_record["design"]["stage_levels"]]

        axis_rows, geo_rows, analysis_info, scores_by_label = [], [], [], {}
        for label, a in analyses:
            scores = pd.read_csv(a["dir"] / "latent_scores.csv", dtype={"donor_id": str}).set_index("donor_id")
            if sorted(scores.index) != donors:
                raise ValueError(f"{a['dir'] / 'latent_scores.csv'} lists different donors from {export_dir / 'metadata.csv'}")
            scores = scores.loc[donors]
            comps = sorted([col for col in scores.columns if col.startswith("pls_")], key=_component_order)
            scores = scores[comps]
            scores_by_label[label] = scores
            rank = int(a["record"]["latent"]["selected_lv"])
            if rank != len(comps):
                raise ValueError(f"{a['dir']}: result.json records rank {rank} but latent_scores.csv has {len(comps)} components")
            info = {
                "label": label, "rank": rank, "selection": a["record"]["latent"]["component_selection"],
                "dir": str(a["dir"]), "files": a["files"],
                "component_sd": {col: float(np.std(scores[col].to_numpy(dtype=float), ddof=1)) for col in comps},
            }
            analysis_info.append(info)
            for row in _correlation_rows(scores, variables, z, threshold=DIAGNOSTIC_CORRELATION_THRESHOLD):
                axis_rows.append({"analysis": label, "rank": rank, **row, "source": sources[row["variable"]]["source"]})
            stage_means = pd.read_csv(a["dir"] / "stage_means.csv", dtype={"group": str, "stage": str})
            for row in geometry_shares(stage_means, groups=groups, stages=stages):
                geo_rows.append({"analysis": label, "rank": rank, **row})
            print(f"[{log}] {c}: {label} (rank {rank}): "
                  + ", ".join(f"{col} sd {info['component_sd'][col]:.4g}" for col in comps[:4])
                  + (" …" if len(comps) > 4 else ""))

        # Nestedness of pls_0 across forced ranks
        p0 = scores_by_label[DIAGNOSTIC_PRIMARY_LABEL]["pls_0"].to_numpy(dtype=float)
        nested = {}
        for label, scores in scores_by_label.items():
            if label == DIAGNOSTIC_PRIMARY_LABEL:
                continue
            other = scores["pls_0"].to_numpy(dtype=float)
            nested[label] = {
                "allclose": bool(np.allclose(other, p0)),
                "max_abs_diff": float(np.max(np.abs(other - p0))),
                "spearman_with_primary": spearman(other, p0)[0],
            }
        if nested:
            same = [k for k, v in nested.items() if v["allclose"]]
            differ = [k for k, v in nested.items() if not v["allclose"]]
            if differ:
                findings.append(
                    f"pls_0 is not identical across ranks: differs from the primary's in {differ} (max |diff| "
                    + ", ".join(f"{k} {nested[k]['max_abs_diff']:.3g}" for k in differ) + f"); identical in {same}"
                )
            else:
                findings.append(
                    f"pls_0 is identical (np.allclose) across the primary and every forced rank {same}: PLS "
                    "components nest, so every forced rank keeps the primary's first axis"
                )

        # Block principal components
        pc_rows, pc_info = [], {}
        for o in ANALYSIS_BLOCKS:
            block = pd.read_csv(export_dir / f"{o}.csv", index_col=0)
            block.index = block.index.astype(str)
            if sorted(block.index) != donors:
                raise ValueError(f"{export_dir / f'{o}.csv'} lists different donors from {export_dir / 'metadata.csv'}")
            block = block.loc[donors]
            pcs = block_principal_components(block)
            pc_info[o] = {"variance_share": pcs["variance_share"], "n_features": pcs["n_features"]}
            for row in _correlation_rows(pcs["scores"], variables, z, threshold=DIAGNOSTIC_CORRELATION_THRESHOLD):
                pc_rows.append({
                    "block": o, **{k: v for k, v in row.items() if k != "component_sd"},
                    "variance_share": pcs["variance_share"][row["component"]],
                    "source": sources[row["variable"]]["source"],
                })
            print(f"[{log}] {c} {o}: PC variance shares " + ", ".join(f"{k} {v:.3f}" for k, v in pcs["variance_share"].items()))

        # Feature selection accounting
        feat_rows = [
            feature_selection_accounting(export_dir, o, scores_by_label[DIAGNOSTIC_PRIMARY_LABEL]["pls_0"])
            for o in ANALYSIS_BLOCKS
        ]
        for r in feat_rows:
            print(f"[{log}] {c} {r['omic']}: kept {r['n_kept']} of {r['n_abundance_passing']} abundance-passing; "
                  f"median kept mean {r['mean_q50']:.3f}; share mean < 1 {r['share_mean_below_1']:.3f}; variance share "
                  f"{r['variance_share_kept']:.3f}; |rho(pls_0)| ≥ {DIAGNOSTIC_CORRELATION_THRESHOLD}: {r['share_correlated_pls_0']:.3f}")

        # Technical table provenance against the export's pseudobulk
        export_manifest = json.loads((export_dir / "manifest.json").read_text())
        pb_hash_export = export_manifest.get("input", {}).get("pseudobulk_manifest_sha256")
        pb_hash_technical = (technical_manifest or {}).get("pseudobulk", {}).get("manifest_sha256")
        if technical_manifest is None:
            findings.append(f"{technical_manifest_path} not found: the technical table's pseudobulk provenance was not checked")
        elif pb_hash_export != pb_hash_technical:
            findings.append(
                f"the technical table was made from a pseudobulk manifest with sha256 {pb_hash_technical} but the export "
                f"records {pb_hash_export}; the technical variables and the blocks may come from different pseudobulks"
            )
        flagged_primary = [
            r for r in axis_rows
            if r["analysis"] == DIAGNOSTIC_PRIMARY_LABEL and r["component"] == "pls_0" and (r["flag_raw"] or r["flag_partial"])
        ]
        findings.append(
            f"primary pls_0: {len(flagged_primary)} of {len(variables.columns)} variables flagged at |rho| ≥ "
            f"{DIAGNOSTIC_CORRELATION_THRESHOLD} (raw or partial): "
            + ", ".join(f"{r['variable']} ({_fmt(r['spearman'], '.2f')}; partial {_fmt_nan(r['partial_spearman'], '.2f')})"
                        for r in sorted(flagged_primary, key=lambda r: -abs(r["spearman"])))
        )

        axis_cols = ["analysis", "rank", "component", "component_sd", "variable", "source", "n", "spearman",
                     "n_partial", "partial_spearman", "flag_raw", "flag_partial"]
        pc_cols = ["block", "component", "variance_share", "variable", "source", "n", "spearman", "n_partial",
                   "partial_spearman", "flag_raw", "flag_partial"]
        geo_cols = ["analysis", "rank", "quantity", "component", "value", "along", "total", "share"]
        tables = {
            "axis": pd.DataFrame(axis_rows, columns=axis_cols),
            "block_pcs": pd.DataFrame(pc_rows, columns=pc_cols),
            "geometry": pd.DataFrame(geo_rows, columns=geo_cols),
            "features": pd.DataFrame(feat_rows),
        }
        tmp = out_dir / f"{c}.tmp"
        if tmp.exists():
            shutil.rmtree(tmp)
        tmp.mkdir(parents=True)
        for key in ("axis", "block_pcs", "geometry", "features"):
            tables[key].to_csv(tmp / DIAGNOSTIC_OUTPUTS[key], index=False, float_format=_csv_float)
        (tmp / DIAGNOSTIC_OUTPUTS["markdown"]).write_text(
            render_axis_diagnostics(
                tables, cohort=c, analyses=analysis_info, sources=sources, findings=findings,
                threshold=DIAGNOSTIC_CORRELATION_THRESHOLD,
            )
        )
        export_files = {
            name: {"path": str(export_dir / name), "sha256": _sha256_file(export_dir / name)}
            for name in _ANALYSIS_INPUTS + tuple(f"features_{o}.csv" for o in ANALYSIS_BLOCKS)
        }
        manifest = {
            "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "cell_type": cell_type,
            "cohort": c,
            "n_donors": int(len(donors)),
            "threshold": DIAGNOSTIC_CORRELATION_THRESHOLD,
            "threshold_note": "a reporting flag, not an adoption rule; no p-value is attached to any correlation",
            "inputs": {
                "export": {"dir": str(export_dir), "files": export_files, "pseudobulk_manifest_sha256": pb_hash_export},
                "analyses": [
                    {"label": a["label"], "dir": a["dir"], "selection": a["selection"], "rank": a["rank"], "files": a["files"]}
                    for a in analysis_info
                ],
                "technical": {
                    "table": {"path": str(technical_path), "sha256": _sha256_file(technical_path)},
                    "manifest": {
                        "path": str(technical_manifest_path),
                        "sha256": _sha256_file(technical_manifest_path) if technical_manifest is not None else None,
                        "pseudobulk_manifest_sha256": pb_hash_technical,
                        "matches_export": None if technical_manifest is None else bool(pb_hash_export == pb_hash_technical),
                    },
                },
                "rank_dir": str(rank_dir / c),
            },
            "variables": [{"name": v, **sources[v], "n_missing": n_missing.get(v, 0)} for v in sorted(variables.columns)],
            "variable_sources": DIAGNOSTIC_VARIABLES,
            "conditioning_set": list(ANALYSIS_COVARIATES),
            "definitions": {
                "spearman": "scipy.stats.spearmanr over the donors with both values present (n recorded)",
                "partial_spearman": (
                    "Spearman correlation of the OLS residuals of the component and the variable on an intercept and the "
                    "six analysis covariates, over the donors complete in all; n/a for a variable in the conditioning set"
                ),
                "component_sd": "sample SD (ddof 1) of the component's latent scores over the cohort donors",
                "block_pcs": (
                    f"SVD of the z-scored block (population SD); first {DIAGNOSTIC_N_BLOCK_PCS} components' scores U S, "
                    "sign fixed so the largest-magnitude loading is positive; variance_share = s^2 / sum s^2"
                ),
                "geometry": _DIAGNOSTIC_GEOMETRY,
                "geometry_share": "squared length along the component over the squared length over all components",
                "feature_selection": (
                    "kept features' mean log2-CPM quantiles and shares below the cuts; variance_share_kept = variance of the "
                    "kept features over the variance of all abundance-passing features; variance at the cut = kept minimum "
                    "and dropped (status variance) median; share_correlated_pls_0 = share of kept features with "
                    f"|Spearman(feature, primary pls_0)| >= {DIAGNOSTIC_CORRELATION_THRESHOLD}"
                ),
                "nestedness": "np.allclose of each forced rank's pls_0 with the primary's",
            },
            "block_pcs": pc_info,
            "nestedness": nested,
            "n_flagged": {
                "axis_raw": int(tables["axis"]["flag_raw"].sum()),
                "axis_partial": int(tables["axis"]["flag_partial"].sum()),
                "block_pcs_raw": int(tables["block_pcs"]["flag_raw"].sum()),
                "block_pcs_partial": int(tables["block_pcs"]["flag_partial"].sum()),
            },
            "downstream": "read by no export or analysis command; preprocessing variants are predeclared in a change proposal, never derived from a diagnostic value",
            "outputs": {
                DIAGNOSTIC_OUTPUTS[k]: _sha256_file(tmp / DIAGNOSTIC_OUTPUTS[k])
                for k in ("axis", "block_pcs", "geometry", "features", "markdown")
            },
            "software": _software_versions(),
            "git": _git_revision(),
            "findings": findings,
        }
        (tmp / DIAGNOSTIC_OUTPUTS["manifest"]).write_text(json.dumps(manifest, indent=2) + "\n")
        final = out_dir / c
        if final.exists():
            shutil.rmtree(final)
        os.replace(tmp, final)
        print(f"[{log}] Wrote {final}")
        for msg in findings:
            print(f"[{log}] FINDING: {c}: {msg}")
        results[c] = {"dir": final, "tables": tables, "analyses": analysis_info, "findings": findings, "manifest": manifest}
    return results


# ---------------------------------------------------------------------------
# MOTCO latent-adjustment sensitivity (covariate set and latent basis beside the primary analysis)
# ---------------------------------------------------------------------------

# MOTCO's covariate_latent_adjustment (covariate-aware latent construction)
LATENT_MOTCO_MIN_VERSION = "0.9.0"
# Added after ANALYSIS_COVARIATES, in this order, on the raw fraction scale of
# donor_technical.csv (the scale of the `residualized` variant): the two
# fractions the technical-axis report classified as residualization candidates
# (between-donor share dominant, within-donor supertype share small), i.e.
# RESIDUALIZATION_COVARIATES without the export's log10 library sizes (the
# nuclei counts already enter; MOTCO's Type I check used eight covariates).
# A constant, never a flag, never chosen from results; flag and ask before changing.
LATENT_SENSITIVITY_ADDED_COVARIATES = ("ribosomal_pseudobulk", "mitochondrial_pseudobulk")
LATENT_SENSITIVITY_COVARIATES = ANALYSIS_COVARIATES + LATENT_SENSITIVITY_ADDED_COVARIATES
if LATENT_SENSITIVITY_ADDED_COVARIATES != tuple(c for c in RESIDUALIZATION_COVARIATES if not c.startswith("log10_lib_")):
    raise AssertionError(
        f"LATENT_SENSITIVITY_ADDED_COVARIATES {LATENT_SENSITIVITY_ADDED_COVARIATES} must equal RESIDUALIZATION_COVARIATES "
        f"without log10_lib_* ({RESIDUALIZATION_COVARIATES}); a change to one set is a change to both"
    )
# MOTCO's admissible conditioning (Type I adopt_stage); group_stage is not run
LATENT_SENSITIVITY_CONDITIONING = "stage"
# The arms, a registry and not a flag, each run once per cohort on the primary
# export with every other analysis setting equal to the primary's. The primary
# analysis (six covariates, covariate-blind basis) is the reference row and
# stays primary; primary → extended_covariates is the effect of the covariate
# set, extended_covariates → latent_adjusted the effect of the basis alone.
LATENT_SENSITIVITY_ARMS: dict[str, dict] = {
    "extended_covariates": {
        "description": (
            "the six analysis covariates plus the ribosomal and mitochondrial pseudobulk fractions of the technical "
            "table, entered additively in the trajectory model; covariate-blind latent basis (no latent adjustment)"
        ),
        "covariates": LATENT_SENSITIVITY_COVARIATES,
        "latent_adjustment": None,
    },
    "latent_adjusted": {
        "description": (
            "the same eight covariates, with the latent basis fitted on features residualized on them conditional on "
            "stage (MOTCO's covariate_latent_adjustment = 'stage'); trajectory model, LS-means and RRPP as in "
            "extended_covariates"
        ),
        "covariates": LATENT_SENSITIVITY_COVARIATES,
        "latent_adjustment": LATENT_SENSITIVITY_CONDITIONING,
    },
}
LATENT_SENSITIVITY_REFERENCE_LABEL = "primary"
# Layout under <out_dir>: <arm>/motco_analysis/<cohort>/ mirrors results/motco_analysis
# so axis-diagnostics and rank-sensitivity can read an arm; comparison/<cohort>/
LATENT_SENSITIVITY_ANALYSIS_DIR = "motco_analysis"
LATENT_SENSITIVITY_COMPARISON_DIR = "comparison"
LATENT_SENSITIVITY_TABLE = "latent_sensitivity.csv"
LATENT_SENSITIVITY_MARKDOWN = "latent_sensitivity.md"
LATENT_SENSITIVITY_MANIFEST = "latent_sensitivity_manifest.json"
LATENT_SENSITIVITY_MEDIATOR_TABLE = "mediator_check.csv"
# MOTCO's latent_adjustment.feature_r2 key → latent_sensitivity.csv column suffix
LATENT_SENSITIVITY_FEATURE_R2: dict[str, str] = {
    "q05": "q05", "q25": "q25", "q50": "q50", "q75": "q75", "q95": "q95", "mean": "mean",
    "fraction_above_0.5": "frac_gt_0_5",
}
_MEDIATOR_SUBSETS: dict[str, str | None] = {"all": None, "F": "Female", "M": "Male"}


def mediator_check(
    technical: pd.DataFrame,
    metadata: pd.DataFrame,
    columns: Iterable[str] = LATENT_SENSITIVITY_ADDED_COVARIATES,
    *,
    cohort: str,
    threshold: float = DIAGNOSTIC_CORRELATION_THRESHOLD,
) -> pd.DataFrame:
    """Spearman correlation of each added covariate with stage, over all donors and within each sex.

    Parameters
    ----------
    technical:
        A table with ``donor_id`` and the ``columns`` (``donor_technical.csv``).
    metadata:
        The cohort's design table (``donor_id``, ``sex`` as ``Female`` /
        ``Male``, ``stage`` as ``0`` / ``1`` / ``2``); its donors are the rows used.
    columns:
        The covariates checked, in output order.
    cohort:
        Written in the ``cohort`` column.
    threshold:
        A cell with ``|ρ| >= threshold`` is flagged as a mediator candidate
        (``DIAGNOSTIC_CORRELATION_THRESHOLD``).

    Returns
    -------
    pandas.DataFrame
        One row per covariate × subset (``all``, ``F``, ``M``): ``cohort``,
        ``covariate``, ``subset``, ``n`` (pairwise-complete donors),
        ``spearman`` (NaN with fewer than three donors or a constant) and
        ``flagged``.

    Notes
    -----
    Predeclared and written before any arm runs. No p-value is computed, and
    the correlation is not partial on the analysis covariates: the question is
    whether the fraction moves with stage at all. A flag removes no covariate,
    changes no arm and stops no run; it binds the report to say that the
    ``latent_adjusted`` basis may have removed stage signal carried by that
    fraction (no simulation can tell a confounder from a mediator).
    """
    meta = metadata[["donor_id", "sex", "stage"]].copy()
    meta["donor_id"] = meta["donor_id"].astype(str)
    unknown = sorted(set(meta["sex"].astype(str)) - {s for s in _MEDIATOR_SUBSETS.values() if s})
    if unknown:
        raise ValueError(f"mediator_check: sex label(s) {unknown} are neither Female nor Male")
    tech = technical.copy()
    tech["donor_id"] = tech["donor_id"].astype(str)
    tech = tech.set_index("donor_id")
    absent = sorted(set(meta["donor_id"]) - set(tech.index))
    if absent:
        raise ValueError(f"mediator_check: the technical table has no row for donor(s) {absent}")
    stage = meta["stage"].astype(int).to_numpy(dtype=float)
    sex = meta["sex"].astype(str).to_numpy()
    rows = []
    for col in columns:
        values = pd.to_numeric(tech.loc[meta["donor_id"], col], errors="coerce").to_numpy(dtype=float)
        for subset, label in _MEDIATOR_SUBSETS.items():
            keep = np.ones(len(meta), dtype=bool) if label is None else sex == label
            rho, n = spearman(values[keep], stage[keep])
            rows.append({
                "cohort": cohort,
                "covariate": col,
                "subset": subset,
                "n": n,
                "spearman": rho,
                "flagged": bool(np.isfinite(rho) and abs(rho) >= threshold),
            })
    return pd.DataFrame(rows, columns=["cohort", "covariate", "subset", "n", "spearman", "flagged"])


def _latent_row(record: dict, *, arm: str, cohort: str, directory: str) -> dict:
    """``_sensitivity_row`` with the arm, covariate set, conditioning and latent-adjustment block.

    The latent-adjustment columns are ``None`` (empty in the CSV) for a record
    without MOTCO's ``latent_adjustment`` block; ``conditioning`` is ``none``.
    """
    base = _sensitivity_row(
        record, cohort=cohort, component_selection=record["latent"]["component_selection"], directory=directory
    )
    adjustment = record.get("latent_adjustment")
    row = {
        "arm": arm,
        "cohort": base.pop("cohort"),
        "component_selection": base.pop("component_selection"),
        "rank": base.pop("rank"),
        "covariates": ";".join((record.get("covariates") or {}).get("columns", [])),
        "conditioning": "none" if adjustment is None else adjustment["conditioning"],
    }
    directory = base.pop("dir")
    row.update(base)
    row["residual_df"] = None if adjustment is None else adjustment["residual_df"]
    row["cv_min_training_df"] = None if adjustment is None else adjustment["cv_min_training_df"]
    r2 = {} if adjustment is None else adjustment.get("feature_r2", {})
    for key, suffix in LATENT_SENSITIVITY_FEATURE_R2.items():
        row[f"feature_r2_{suffix}"] = r2.get(key)
    row["dir"] = directory
    return row


def render_latent_sensitivity(
    rows: list[dict],
    mediator: pd.DataFrame,
    *,
    cohort: str,
    seed: int = ANALYSIS_SEED,
    permutations: int = ANALYSIS_PERMUTATIONS,
    alpha: float = ANALYSIS_ALPHA,
    threshold: float = DIAGNOSTIC_CORRELATION_THRESHOLD,
    findings: list[str] | None = None,
) -> str:
    """Render ``latent_sensitivity.md`` from the rows of ``latent_sensitivity.csv`` and the mediator check.

    Parameters
    ----------
    rows:
        ``_latent_row`` dicts: the primary (reference) row first, then the arms
        in ``LATENT_SENSITIVITY_ARMS`` order.
    mediator:
        ``mediator_check`` output for the cohort.
    cohort:
        Cohort name for the title.
    seed, permutations, alpha, threshold:
        Named in the text; a statistic is rejected when ``p <= alpha``.
    findings:
        The comparison's findings, listed at the end.

    Returns
    -------
    str
        Markdown: a header with the arms and the reading rules, one table per
        statistic (observed, p, rejected, null mean / sd / q95 per row), a run
        table (selection, rank, CV AUROC, covariates, conditioning, runtime,
        directory), the latent-adjustment block of the adjusted arm(s), the
        mediator check with flagged rows first, and the findings. Every number
        comes from the rows and the mediator table; the runtimes are the one
        thing here that is not in ``latent_sensitivity.csv``.
    """
    ref = rows[0]
    added = ", ".join(f"`{c}`" for c in LATENT_SENSITIVITY_ADDED_COVARIATES)
    lines = [
        f"# MOTCO latent-adjustment sensitivity: {cohort}",
        "",
        f"Covariate-set and latent-basis sensitivity of the covariate-adjusted trajectory analysis. The "
        f"`{ref['arm']}` row is the primary analysis (CV-selected rank {ref['rank']}, the six analysis covariates, "
        "covariate-blind basis), read from its `result.json` and never rerun; it stays the primary result, and "
        "MOTCO's power comparison (`candidate_fails`) excludes the latent adjustment as its replacement. "
        f"`extended_covariates` adds {added} (raw fractions from `donor_technical.csv`) with the covariate-blind "
        f"basis; `latent_adjusted` enters the same eight covariates and fits the basis on features residualized on "
        f"them conditional on `{LATENT_SENSITIVITY_CONDITIONING}`, with the trajectory model, LS-means and RRPP of "
        "`extended_covariates`, so those two rows differ through the basis and the CV-selected rank only. Both arms "
        "are read as sensitivity, never as the result: a rejection status that differs between `primary` and "
        "`extended_covariates` is covariate-set dependence, one that differs between `extended_covariates` and "
        "`latent_adjusted` is basis dependence; a CV-selected rank that differs is reported beside `angle`, not read "
        "as better or worse; `angle` and `shape` non-rejections are uninformative and a `delta` rejection is a "
        f"difference, not specifically of size, in every row. {permutations} RRPP permutations with seed {seed} in "
        f"every row (the same permutation indices); α = {alpha} (rejected when p ≤ α). MOTCO measured the latent "
        "adjustment's operating characteristics at its Phase 6 design only (Type I `adopt_stage`; power "
        "`candidate_fails`, with `shape` power lost under a nuisance axis and small-magnitude power lost without "
        "one), on neither these blocks nor this eight-covariate set, so neither arm carries a Type I or power claim "
        "of its own.",
    ]
    for s in ANALYSIS_STATISTICS:
        lines += [
            "",
            f"## {s}" + (" (degrees)" if s == "angle" else ""),
            "",
            "| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |",
            "|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|",
        ]
        for r in rows:
            rejected = r[f"{s}_rejected"]
            n_cov = len(r["covariates"].split(";")) if r["covariates"] else 0
            lines.append(
                f"| {r['arm']} | {n_cov} | {r['conditioning']} | {r['rank']} | {_fmt(r[s])} | "
                f"{_fmt(r[f'{s}_p'], '.4g')} | {'n/a' if rejected is None else ('yes' if rejected else 'no')} | "
                f"{_fmt(r[f'{s}_null_mean'])} | {_fmt(r[f'{s}_null_sd'])} | {_fmt(r[f'{s}_null_q95'])} |"
            )
    lines += [
        "",
        "## Runs",
        "",
        "| arm | selection | rank | cv_mean_auroc | covariates (in order) | conditioning | runtime (s) | directory |",
        "|---|---|---:|---:|---|---|---:|---|",
    ]
    for r in rows:
        covs = ", ".join(f"`{c}`" for c in r["covariates"].split(";")) if r["covariates"] else "none"
        lines.append(
            f"| {r['arm']} | {r['component_selection']} | {r['rank']} | {_fmt(r['cv_mean_auroc'], '.4f')} | {covs} | "
            f"{r['conditioning']} | {_fmt(r['runtime_seconds'], '.1f')} | `{r['dir']}` |"
        )
    adjusted = [r for r in rows if r["conditioning"] != "none"]
    if adjusted:
        suffixes = list(LATENT_SENSITIVITY_FEATURE_R2.values())
        labels = {"frac_gt_0_5": "fraction > 0.5"}
        lines += [
            "",
            "## Latent adjustment",
            "",
            "MOTCO's `latent_adjustment` block: the residualization's degrees of freedom on the full sample "
            "(`n − 1 − k − q`) and the smallest over the CV training sets, and the quantiles over features of the "
            "partial R² of the covariates given the intercept and the conditioning indicators (pooled-standardized "
            "features).",
            "",
            "| arm | conditioning | residual df | smallest CV training df | "
            + " | ".join(f"R² {labels.get(x, x)}" for x in suffixes) + " |",
            "|---|---|---:|---:|" + "---:|" * len(suffixes),
        ]
        for r in adjusted:
            lines.append(
                f"| {r['arm']} | {r['conditioning']} | {_fmt(r['residual_df'], 'd')} | "
                f"{'n/a' if r['cv_min_training_df'] is None else _fmt(r['cv_min_training_df'], 'd')} | "
                + " | ".join(_fmt(r[f"feature_r2_{x}"], ".4f") for x in suffixes) + " |"
            )
    flagged = mediator[mediator["flagged"]]
    ordered = pd.concat([flagged, mediator[~mediator["flagged"]]])
    lines += [
        "",
        "## Mediator check",
        "",
        "Spearman correlation of each added covariate with stage (ordinal 0, 1, 2) over all cohort donors and within "
        f"each sex, written before the arms ran; |ρ| ≥ {threshold} (`DIAGNOSTIC_CORRELATION_THRESHOLD`) flags a "
        "mediator candidate (flagged rows first). No p-value is computed, and a flag removes no covariate and changes "
        "no arm: no simulation can tell a confounder from a mediator, so a flagged fraction may carry stage signal "
        "that the `latent_adjusted` basis removed (and that the covariate adjustment of both arms absorbed).",
        "",
        "| covariate | subset | n | Spearman | flagged |",
        "|---|---|---:|---:|:---:|",
    ]
    for _, r in ordered.iterrows():
        lines.append(
            f"| `{r['covariate']}` | {r['subset']} | {r['n']} | {_fmt_nan(r['spearman'], '.3f')} | "
            f"{'yes' if r['flagged'] else 'no'} |"
        )
    lines += ["", "## Findings", ""]
    lines += [f"- {msg}" for msg in (findings or [])] or ["None."]
    return "\n".join(lines) + "\n"


def _same_dir(a: Path, b: Path) -> bool:
    """Whether ``a`` is ``b`` or lies inside it (resolved paths)."""
    a, b = a.resolve(), b.resolve()
    return a == b or b in a.parents


def run_latent_sensitivity(
    cohort: str = "all",
    *,
    cell_type: str = DEFAULT_CELL_TYPE,
    input_dir: str | Path | None = None,
    analysis_dir: str | Path | None = None,
    pseudobulk_dir: str | Path | None = None,
    out_dir: str | Path | None = None,
    n_jobs: int = 1,
    force: bool = False,
) -> dict:
    """Rerun the primary analysis with the extended covariate set, without and with the latent adjustment.

    Parameters
    ----------
    cohort:
        ``"cohort_30"``, ``"cohort_50"`` or ``"all"``.
    cell_type:
        A registered cell type (``CELL_TYPES``, default ``astrocyte``). The
        directory defaults below come from its results root
        (``cell_type_dirs``); every manifest read (export, analysis,
        technical, pseudobulk) must record it, a manifest without the key
        reading as ``astrocyte`` with a finding (``check_manifest_cell_type``),
        and the arm and comparison manifests record it. The registry of arms,
        the covariate set, the conditioning, the mediator check and every
        analysis setting are the same for every cell type.
    input_dir:
        Export root of the ``export`` step (default ``<root>/motco_input``);
        ``<input_dir>/<cohort>/`` is the primary export, read, never written.
        Both arms run on it.
    analysis_dir:
        Analysis root of the ``analyze`` step (default
        ``<root>/motco_analysis``); ``<analysis_dir>/<cohort>/`` is the
        primary analysis (the reference row), read and checked, never rerun
        or written.
    pseudobulk_dir:
        Directory holding ``manifest.json``, ``donor_technical.csv`` and
        ``technical_manifest.json`` (default ``<root>/pseudobulk``; read,
        never written).
    out_dir:
        Output root (default ``<root>/motco_latent_sensitivity``):
        ``<out_dir>/<arm>/motco_analysis/<cohort>/`` per arm of
        ``LATENT_SENSITIVITY_ARMS`` (the files ``analyze`` writes) and
        ``<out_dir>/comparison/<cohort>/`` (``latent_sensitivity.csv``,
        ``latent_sensitivity.md``, ``mediator_check.csv``,
        ``latent_sensitivity_manifest.json``). It may not lie inside an input
        directory.
    n_jobs:
        Parallel jobs for MOTCO's CV and RRPP (``-1`` = all CPUs); as in
        ``run_analysis``, values other than 1 change the permutation draws and
        are recorded as a finding. The reported results use 1.
    force:
        Replace a cohort's existing arm and comparison directories.

    Returns
    -------
    dict
        Cohort → ``{"dir", "arms" (arm → _analyze_cohort result), "reference",
        "technical", "rows", "table", "mediator", "manifest", "findings"}``.

    Raises
    ------
    RuntimeError
        The installed MOTCO is older than ``LATENT_MOTCO_MIN_VERSION`` (raised
        before any input is read).
    FileExistsError
        A cohort's arm or comparison directory exists and ``force`` is false.
    FileNotFoundError
        The primary export, the primary analysis, the pseudobulk manifest or
        the technical table (or its manifest) is missing; the message names
        the file and the command that produces it.
    ValueError
        The primary analysis does not match the export, the installed MOTCO,
        the seed, the permutations or the CV selection (``_check_reference``);
        the technical table does not match the pseudobulk manifest or
        ``TECHNICAL_GENE_SETS``, or lacks an added covariate for a cohort donor
        (``_check_technical_reference``); a manifest read records another
        cell type (``check_manifest_cell_type``); or ``out_dir`` lies inside
        an input directory.

    Notes
    -----
    Every check runs for every requested cohort before anything is written.
    The cell-type findings of the manifests read (a manifest written before
    the cell-type registry) go to the comparison manifest's ``findings``,
    not to ``latent_sensitivity.md``, whose findings are the comparison's.
    Per cohort, ``mediator_check.csv`` is written to
    ``comparison/<cohort>.tmp/`` before the first arm runs; each arm is one
    ``_analyze_cohort`` call on the primary export with
    ``covariates = LATENT_SENSITIVITY_COVARIATES``, the fractions merged from
    the technical table (``load_analysis_inputs(extra_covariates=...)``) and
    ``TrajectoryAnalysisParams(permutations=ANALYSIS_PERMUTATIONS,
    seed=ANALYSIS_SEED, n_jobs=n_jobs, covariate_latent_adjustment=...)``, so
    everything else equals ``analyze`` (blocks, stage order, orientation, CV
    rank selection, seed, permutations). The arm manifests record the five
    export files as ``inputs.files`` (so ``axis-diagnostics`` can read an arm
    against the primary export) and the technical table under
    ``inputs.technical``. The arms are written as ``<cohort>.tmp/`` and the
    three directories renamed once the comparison is complete; on any error
    the cohort's temporaries are removed.
    """
    import shutil
    from importlib.metadata import PackageNotFoundError, version

    log = "latent-sensitivity"
    motco_version = check_motco_version(
        required=LATENT_MOTCO_MIN_VERSION, purpose="latent-sensitivity needs", feature="covariate_latent_adjustment"
    )
    from motco.analysis import TrajectoryAnalysisParams

    cell_type_entry(cell_type)
    flag = "" if cell_type == DEFAULT_CELL_TYPE else f" --cell-type {cell_type}"
    input_dir = _resolve_dir(input_dir, cell_type, "motco_input")
    analysis_dir = _resolve_dir(analysis_dir, cell_type, "motco_analysis")
    pb = _resolve_dir(pseudobulk_dir, cell_type, "pseudobulk")
    out_dir = _resolve_dir(out_dir, cell_type, "motco_latent_sensitivity")
    cohorts = list(EXPORT_COHORTS) if cohort == "all" else [cohort]
    for c in cohorts:
        if c not in EXPORT_COHORTS:
            raise ValueError(f"Unknown cohort {c!r}; expected one of {EXPORT_COHORTS} or 'all'")
    for name, d in (("--input-dir", input_dir), ("--analysis-dir", analysis_dir), ("--pseudobulk-dir", pb)):
        if _same_dir(out_dir, d):
            raise ValueError(
                f"the output directory {out_dir} lies inside the input directory {d} ({name}); the step never "
                "writes under its inputs; nothing was written"
            )
    arms = list(LATENT_SENSITIVITY_ARMS)
    for c in cohorts:
        targets = [out_dir / LATENT_SENSITIVITY_COMPARISON_DIR / c] + [
            out_dir / arm / LATENT_SENSITIVITY_ANALYSIS_DIR / c for arm in arms
        ]
        for t in targets:
            if t.exists() and not force:
                raise FileExistsError(
                    f"{t} exists; rerun with --force to replace the arms and the comparison of {c}; nothing was written"
                )
    for c in cohorts:
        if not (input_dir / c).is_dir():
            raise FileNotFoundError(
                f"{input_dir / c} does not exist; run `sead-sexdiff export{flag}` to produce the primary export"
            )
        if not (analysis_dir / c).is_dir():
            raise FileNotFoundError(
                f"{analysis_dir / c} does not exist; run `sead-sexdiff analyze{flag}` to produce the primary analysis"
            )
    if not (pb / "manifest.json").exists():
        raise FileNotFoundError(
            f"{pb / 'manifest.json'} is missing; run `sead-sexdiff pseudobulk{flag}` on the cluster and sync {pb}/ "
            f"(scripts/sync_pseudobulk.sh{flag}) to produce it"
        )
    print(f"[{log}] MOTCO {motco_version} (latent adjustment needs >= {LATENT_MOTCO_MIN_VERSION})")
    print(f"[{log}] Cell type: {cell_type}")
    print(f"[{log}] Primary export directory: {input_dir}")
    print(f"[{log}] Primary analysis directory: {analysis_dir}")
    print(f"[{log}] Pseudobulk directory: {pb}")
    print(f"[{log}] Output directory: {out_dir}")
    print(
        f"[{log}] Arms: {', '.join(arms)} (registry LATENT_SENSITIVITY_ARMS; covariates "
        f"{', '.join(LATENT_SENSITIVITY_COVARIATES)}; the primary analysis stays primary)"
    )

    params = {
        arm: TrajectoryAnalysisParams(
            permutations=ANALYSIS_PERMUTATIONS, seed=ANALYSIS_SEED, n_jobs=n_jobs,
            covariate_latent_adjustment=spec["latent_adjustment"],
        )
        for arm, spec in LATENT_SENSITIVITY_ARMS.items()
    }
    # Check every cohort's primary export, primary analysis and technical table before writing anything
    prepared = {}
    for c in cohorts:
        cohort_dir = input_dir / c
        loaded = load_analysis_inputs(cohort_dir)
        cell_type_findings = check_manifest_cell_type(loaded["export_manifest"], cohort_dir / "manifest.json", cell_type)
        reference = _check_reference(cohort_dir, analysis_dir / c, loaded["inputs"], motco_version)
        cell_type_findings += check_manifest_cell_type(
            reference["manifest"], reference["dir"] / "analysis_manifest.json", cell_type
        )
        meas = reference["manifest"].get("measurement", {})
        for arm, p in params.items():
            used = {
                "orientation_estimator": p.orientation_estimator,
                "cv": {
                    "cv1_splits": p.cv1_splits, "cv2_splits": p.cv2_splits, "n_repeats": p.n_repeats,
                    "max_components": p.max_components, "random_state": p.random_state,
                },
            }
            for key, want in used.items():
                if meas.get(key) != want:
                    raise ValueError(
                        f"{reference['dir'] / 'analysis_manifest.json'} records {key} = {meas.get(key)!r} but the "
                        f"{arm} arm would use {want!r}; the arms must use the primary analysis's settings; rerun "
                        "`sead-sexdiff analyze`; nothing was written"
                    )
        donors = loaded["metadata"]["donor_id"].tolist()
        technical = _check_technical_reference(pb, donors=donors, columns=LATENT_SENSITIVITY_ADDED_COVARIATES)
        cell_type_findings += check_manifest_cell_type(
            technical["manifest"], pb / TECHNICAL_OUTPUTS["manifest"], cell_type
        )
        cell_type_findings += check_manifest_cell_type(
            json.loads((pb / "manifest.json").read_text()), pb / "manifest.json", cell_type
        )
        extra = technical["table"][["donor_id", *LATENT_SENSITIVITY_ADDED_COVARIATES]].copy()
        extra["donor_id"] = extra["donor_id"].astype(str)
        extra = extra[extra["donor_id"].isin(donors)]
        for col in LATENT_SENSITIVITY_ADDED_COVARIATES:
            extra[col] = pd.to_numeric(extra[col], errors="raise").astype(float)
        loaded_arms = load_analysis_inputs(cohort_dir, extra_covariates=extra)
        if loaded_arms["inputs"] != loaded["inputs"]:
            raise ValueError(f"the export files under {cohort_dir} changed while they were read; nothing was written")
        ref_rank = reference["record"]["latent"]["selected_lv"]
        print(
            f"[{log}] {c}: primary analysis {reference['dir']} (CV-selected rank {ref_rank}; hashes and settings "
            f"match); technical table {pb / TECHNICAL_OUTPUTS['table']} matches the pseudobulk manifest and "
            f"TECHNICAL_GENE_SETS and is complete for {len(donors)} donors × "
            f"{len(LATENT_SENSITIVITY_ADDED_COVARIATES)} added covariates"
        )
        prepared[c] = {
            "loaded": loaded_arms, "reference": reference, "technical": technical,
            "cell_type_findings": cell_type_findings,
        }

    software = _software_versions()
    software["motco"] = motco_version
    try:
        software["scikit-learn"] = version("scikit-learn")
    except PackageNotFoundError:
        software["scikit-learn"] = None
    sources = {
        col: f"{pb / TECHNICAL_OUTPUTS['table']} (raw fraction scale; checked by _check_technical_reference)"
        for col in LATENT_SENSITIVITY_ADDED_COVARIATES
    }

    results = {}
    for c in cohorts:
        prep = prepared[c]
        loaded, reference, technical = prep["loaded"], prep["reference"], prep["technical"]
        cohort_dir = input_dir / c
        comp_root = out_dir / LATENT_SENSITIVITY_COMPARISON_DIR
        comp_tmp = comp_root / f"{c}.tmp"
        arm_roots = {arm: out_dir / arm / LATENT_SENSITIVITY_ANALYSIS_DIR for arm in arms}
        temporaries = [comp_tmp] + [root / f"{c}.tmp" for root in arm_roots.values()]
        technical_inputs = {
            "technical": {
                **technical["files"],
                "pseudobulk_manifest": {
                    "path": str(pb / "manifest.json"), "sha256": technical["pseudobulk_manifest_sha256"],
                },
                "columns": list(LATENT_SENSITIVITY_ADDED_COVARIATES),
                "note": (
                    "the added covariates' source, checked against manifest.json and TECHNICAL_GENE_SETS and complete "
                    "for the cohort donors; not an export file (inputs.files holds the five export files)"
                ),
            },
        }
        try:
            for t in temporaries:
                if t.exists():
                    shutil.rmtree(t)
            comp_tmp.mkdir(parents=True)
            mediator = mediator_check(technical["table"], loaded["metadata"], LATENT_SENSITIVITY_ADDED_COVARIATES, cohort=c)
            mediator.to_csv(comp_tmp / LATENT_SENSITIVITY_MEDIATOR_TABLE, index=False, float_format=_csv_float)
            flagged = mediator[mediator["flagged"]]
            print(
                f"[{log}] {c}: mediator check written before the arms ({len(flagged)} of {len(mediator)} cells at "
                f"|rho| >= {DIAGNOSTIC_CORRELATION_THRESHOLD}"
                + ("" if flagged.empty else ": " + ", ".join(
                    f"{r['covariate']} {r['subset']} {_fmt(r['spearman'], '.3f')}" for _, r in flagged.iterrows()
                )) + ")"
            )

            runs = {}
            for arm, spec in LATENT_SENSITIVITY_ARMS.items():
                runs[arm] = _analyze_cohort(
                    c, cohort_dir, loaded, params[arm],
                    out_dir=arm_roots[arm], final_name=f"{c}.tmp", published_dir=arm_roots[arm] / c, n_jobs=n_jobs,
                    findings=[f"latent-sensitivity arm {arm}: the primary analysis in {reference['dir']} is primary"],
                    log_prefix=log, covariates=spec["covariates"], covariate_sources=sources,
                    extra_inputs=technical_inputs, cell_type=cell_type,
                )
            for arm, spec in LATENT_SENSITIVITY_ARMS.items():
                rec = runs[arm]["record"]
                entered = (rec.get("covariates") or {}).get("columns")
                if entered != list(spec["covariates"]):
                    raise ValueError(f"the {arm} record lists covariates {entered}, not {list(spec['covariates'])}")
                block = rec.get("latent_adjustment")
                conditioning = None if block is None else block.get("conditioning")
                if conditioning != spec["latent_adjustment"]:
                    raise ValueError(
                        f"the {arm} record has latent-adjustment conditioning {conditioning!r}, not "
                        f"{spec['latent_adjustment']!r}"
                    )
            covariate_blocks = {
                arm: {k: runs[arm]["record"]["covariates"][k] for k in ("columns", "kinds", "design_columns", "reference_row")}
                for arm in arms
            }
            if any(b != covariate_blocks[arms[0]] for b in covariate_blocks.values()):
                raise ValueError(
                    f"the arms' covariate blocks differ ({covariate_blocks}); the latent adjustment must change only "
                    "the basis"
                )

            rows = [_latent_row(
                reference["record"], arm=LATENT_SENSITIVITY_REFERENCE_LABEL, cohort=c, directory=str(reference["dir"])
            )]
            rows += [
                _latent_row(runs[arm]["record"], arm=arm, cohort=c, directory=f"{arm}/{LATENT_SENSITIVITY_ANALYSIS_DIR}/{c}")
                for arm in arms
            ]
            for r in rows:
                print(
                    f"[{log}] {c}: {r['arm']} (rank {r['rank']}, conditioning {r['conditioning']}): "
                    + ", ".join(f"{s} {_fmt(r[s])} (p {_fmt(r[f'{s}_p'], '.4g')})" for s in ANALYSIS_STATISTICS)
                )

            findings: list[str] = []
            by_arm = {r["arm"]: r for r in rows}
            ref_label = LATENT_SENSITIVITY_REFERENCE_LABEL
            pairs = (
                (ref_label, "extended_covariates", "covariate-set dependence"),
                ("extended_covariates", "latent_adjusted", "basis dependence"),
            )
            for a, b, kind in pairs:
                if by_arm[a]["rank"] != by_arm[b]["rank"]:
                    findings.append(
                        f"CV-selected rank differs: {by_arm[a]['rank']} in {a}, {by_arm[b]['rank']} in {b}; reported "
                        "beside angle, not read as better or worse"
                    )
                for s in ANALYSIS_STATISTICS:
                    ra, rb = by_arm[a][f"{s}_rejected"], by_arm[b][f"{s}_rejected"]
                    if ra is not None and rb is not None and ra != rb:
                        findings.append(
                            f"{s}: rejection status differs between {a} (p = {_fmt(by_arm[a][f'{s}_p'], '.4g')}) and "
                            f"{b} (p = {_fmt(by_arm[b][f'{s}_p'], '.4g')}); {kind} of the sensitivity analysis, not a "
                            "result of the arm; the primary result stays primary"
                        )
            for _, r in flagged.iterrows():
                findings.append(
                    f"mediator candidate: {r['covariate']} has Spearman {_fmt(r['spearman'], '.3f')} with stage "
                    f"({r['subset']}, n = {r['n']}); the latent_adjusted basis may have removed stage signal carried "
                    "by it (no simulation can tell a confounder from a mediator; the flag changes no arm)"
                )
            if n_jobs != 1:
                findings.append(
                    f"n_jobs = {n_jobs}: MOTCO's parallel RRPP seeds one stream per worker, so the permutation "
                    "draws, p-values and null summaries differ from an n_jobs = 1 run (statistics, rank, latent "
                    "scores and loadings do not); the case-study results are reported at n_jobs = 1"
                )

            table = pd.DataFrame(rows).drop(columns=list(_RANK_SENSITIVITY_CSV_EXCLUDED))
            for col in ("residual_df", "cv_min_training_df"):
                table[col] = table[col].astype("Int64")
            table.to_csv(comp_tmp / LATENT_SENSITIVITY_TABLE, index=False, float_format=_csv_float)
            (comp_tmp / LATENT_SENSITIVITY_MARKDOWN).write_text(
                render_latent_sensitivity(
                    rows, mediator, cohort=c, seed=ANALYSIS_SEED, permutations=ANALYSIS_PERMUTATIONS,
                    alpha=ANALYSIS_ALPHA, findings=findings,
                )
            )
            reference_now = {name: _sha256_file(Path(info["path"])) for name, info in reference["files"].items()}
            if reference_now != {name: info["sha256"] for name, info in reference["files"].items()}:
                raise ValueError(f"the primary analysis files under {reference['dir']} changed during the run")
            ref_record = reference["record"]
            manifest_findings = [*prep["cell_type_findings"], *findings]
            manifest = {
                "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "cell_type": cell_type,
                "cohort": c,
                "n_donors": int(len(loaded["metadata"])),
                "covariates": {
                    "entered": list(LATENT_SENSITIVITY_COVARIATES),
                    "primary": list(ANALYSIS_COVARIATES),
                    "added": list(LATENT_SENSITIVITY_ADDED_COVARIATES),
                    "sources": sources,
                    "scale": "raw fraction (as the residualized variant)",
                    "missing_policy": ANALYSIS_MISSING_COVARIATES,
                    "imputed": False,
                },
                "conditioning": LATENT_SENSITIVITY_CONDITIONING,
                "arms": {
                    arm: {
                        "description": spec["description"],
                        "covariates": list(spec["covariates"]),
                        "latent_adjustment": spec["latent_adjustment"],
                        "dir": str(arm_roots[arm] / c),
                        "files": {
                            name: {"path": str(arm_roots[arm] / c / name), "sha256": _sha256_file(runs[arm]["dir"] / name)}
                            for name in _REFERENCE_FILES
                        },
                        "component_selection": runs[arm]["record"]["latent"]["component_selection"],
                        "selected_lv": runs[arm]["record"]["latent"]["selected_lv"],
                        "cv_mean_auroc": runs[arm]["record"]["latent"].get("cv_mean_auroc"),
                        "runtime_seconds": runs[arm]["record"]["measurement"]["runtime_seconds"],
                        "latent_adjustment_block": runs[arm]["record"].get("latent_adjustment"),
                        "findings": runs[arm]["findings"],
                    }
                    for arm, spec in LATENT_SENSITIVITY_ARMS.items()
                },
                "seed": ANALYSIS_SEED,
                "permutations": ANALYSIS_PERMUTATIONS,
                "same_seed_note": (
                    "the same seed and permutation count as the primary analysis in every arm, so every row draws the "
                    "same permutation indices"
                ),
                "alpha": ANALYSIS_ALPHA,
                "rejection_rule": ANALYSIS_REJECTION_RULE,
                "n_jobs": n_jobs,
                "n_jobs_note": "RRPP draws depend on n_jobs (one seed stream per worker); rerun with the same n_jobs",
                "inputs": {
                    "export": {"dir": str(cohort_dir), "files": loaded["inputs"]},
                    **technical_inputs,
                },
                "reference": {
                    "label": LATENT_SENSITIVITY_REFERENCE_LABEL,
                    "dir": str(reference["dir"]),
                    "files": reference["files"],
                    "covariates": (ref_record.get("covariates") or {}).get("columns"),
                    "component_selection": ref_record["latent"]["component_selection"],
                    "selected_lv": ref_record["latent"]["selected_lv"],
                    "cv_mean_auroc": ref_record["latent"].get("cv_mean_auroc"),
                    "runtime_seconds": ref_record["measurement"]["runtime_seconds"],
                    "motco": reference["manifest"]["software"]["motco"],
                    "seed": reference["manifest"]["measurement"]["seed"],
                    "permutations": reference["manifest"]["measurement"]["permutations"],
                    "note": (
                        "read and checked (export hashes, MOTCO version, seed, permutations, CV selection, orientation "
                        "and CV settings), never rerun or modified; its files hash identically after the run"
                    ),
                },
                "mediator_check": {
                    "file": LATENT_SENSITIVITY_MEDIATOR_TABLE,
                    "sha256": _sha256_file(comp_tmp / LATENT_SENSITIVITY_MEDIATOR_TABLE),
                    "threshold": DIAGNOSTIC_CORRELATION_THRESHOLD,
                    "statistic": "Spearman with stage (ordinal 0, 1, 2), all donors and per sex; no p-value",
                    "written_before_arms": True,
                    "note": "a flag removes no covariate, changes no arm and stops no run",
                },
                "mediator_flags": [
                    {"covariate": r["covariate"], "subset": r["subset"], "n": int(r["n"]), "spearman": float(r["spearman"])}
                    for _, r in flagged.iterrows()
                ],
                "outputs": sorted([
                    LATENT_SENSITIVITY_TABLE, LATENT_SENSITIVITY_MARKDOWN, LATENT_SENSITIVITY_MEDIATOR_TABLE,
                    LATENT_SENSITIVITY_MANIFEST,
                ]),
                "software": software,
                "git": _git_revision(),
                "findings": manifest_findings,
            }
            (comp_tmp / LATENT_SENSITIVITY_MANIFEST).write_text(json.dumps(manifest, indent=2) + "\n")
            # Arms first, the comparison last: a complete comparison directory marks a complete cohort
            finals = [(arm_roots[arm] / f"{c}.tmp", arm_roots[arm] / c) for arm in arms]
            finals.append((comp_tmp, comp_root / c))
            for tmp, final in finals:
                if final.exists():
                    shutil.rmtree(final)
                os.replace(tmp, final)
        except BaseException:
            for t in temporaries:
                if t.exists():
                    shutil.rmtree(t, ignore_errors=True)
            raise
        for arm in arms:
            runs[arm]["dir"] = arm_roots[arm] / c
        final = comp_root / c
        print(f"[{log}] Wrote {final} and " + ", ".join(str(arm_roots[arm] / c) for arm in arms))
        for msg in manifest_findings:
            print(f"[{log}] FINDING: {c}: {msg}")
        results[c] = {
            "dir": final, "arms": runs, "reference": reference, "technical": technical, "rows": rows,
            "table": table, "mediator": mediator, "manifest": manifest, "findings": manifest_findings,
        }
    return results


# ---------------------------------------------------------------------------
# Pseudobulk stage DE (per-feature DESeq2 of stage and sex × stage beside the primary analysis)
# ---------------------------------------------------------------------------

# pydeseq2's Benjamini–Hochberg padj with its defaults (independent filtering,
# Cook's outlier handling); a feature is significant when padj <= α. A
# constant, never a flag, never chosen from results; flag and ask before changing.
STAGE_DE_ALPHA = 0.05
# Treatment coding: the reference level is the first category of each factor
STAGE_DE_REFERENCE_LEVELS: dict[str, str] = {"sex": "Female", "stage": "0"}
_STAGE_DE_LEVELS: dict[str, tuple[str, ...]] = {"sex": ("Female", "Male"), "stage": tuple(STAGES)}
if any(_STAGE_DE_LEVELS[v][0] != ref for v, ref in STAGE_DE_REFERENCE_LEVELS.items()):
    raise AssertionError("STAGE_DE_REFERENCE_LEVELS must name the first level of each factor")
# features_<omic>.csv statuses of the tested features: everything that passed
# the export's chromosome and abundance filters, in the export's order;
# STAGE_DE_EXPORTED_STATUS marks MOTCO's block (the `exported` column)
STAGE_DE_FEATURE_STATUSES = ("variance", "kept")
STAGE_DE_EXPORTED_STATUS = "kept"
# The statuses a primary export writes; any other stops the step
_STAGE_DE_EXPORT_STATUSES = ("sex_chromosome", "abundance", "variance", "kept")
# The summary's reporting sets: `exported` (MOTCO's features, the primary
# reading) and `abundance` (every tested feature). padj is DESeq2's over all
# tested features in both, never recomputed on the exported subset.
STAGE_DE_FEATURE_SETS = ("exported", "abundance")
# Continuous covariates per omic on their recorded scale (covariates.csv; the
# omic's own nuclei count on the log10 scale, derived as the analysis derives
# it): ANALYSIS_COVARIATES without the other omic's nuclei count. Nothing is
# imputed. Flag and ask before changing.
STAGE_DE_COVARIATES: dict[str, tuple[str, ...]] = {
    "rna": ("age_at_death", "pmi", "rin", "years_education", "log10_n_rna"),
    "atac": ("age_at_death", "pmi", "rin", "years_education", "log10_n_atac"),
}
if any(
    cols != tuple(c for c in ANALYSIS_COVARIATES if c not in ANALYSIS_LOG10_COLUMNS or c == f"log10_n_{o}")
    for o, cols in STAGE_DE_COVARIATES.items()
):
    raise AssertionError(
        f"STAGE_DE_COVARIATES {STAGE_DE_COVARIATES} must be ANALYSIS_COVARIATES {ANALYSIS_COVARIATES} without the "
        "other omic's log10 nuclei count"
    )
# The predeclared models. `formula` holds the design terms (the omic's
# covariates are appended by stage_de_formula); `tested` lists (label, design
# column) per Wald-tested coefficient, the column resolved by name in the
# fitted design matrix. A factor model's two tests are also reported as their
# union, which carries no joint error-rate claim (pydeseq2 0.5.3 has no
# likelihood-ratio test). factor_interaction is the primary interaction model
# (MOTCO's ordered-categorical stages); trend_interaction is a second
# predeclared model, not a fallback. A registry, never a flag, never chosen
# from results; flag and ask before changing.
STAGE_DE_MODELS: dict[str, dict] = {
    "factor_additive": {
        "formula": "~ sex + stage",
        "stage_coding": "factor",
        "interaction": False,
        "tested": (("stage_1_vs_0", "stage[T.1]"), ("stage_2_vs_0", "stage[T.2]")),
    },
    "factor_interaction": {
        "formula": "~ sex + stage + sex:stage",
        "stage_coding": "factor",
        "interaction": True,
        "tested": (("male_x_stage_1", "sex[T.Male]:stage[T.1]"), ("male_x_stage_2", "sex[T.Male]:stage[T.2]")),
    },
    "trend_additive": {
        "formula": "~ sex + stage_num",
        "stage_coding": "numeric",
        "interaction": False,
        "tested": (("stage_num", "stage_num"),),
    },
    "trend_interaction": {
        "formula": "~ sex + stage_num + sex:stage_num",
        "stage_coding": "numeric",
        "interaction": True,
        "tested": (("male_x_stage_num", "sex[T.Male]:stage_num"),),
    },
}
STAGE_DE_UNION = "union"
# Written before either count existed (proposal); evaluated per cohort × omic
# from the generated files. The cell reads a count and the primary analysis's
# delta rejection status (ANALYSIS_REJECTION_RULE, as every comparison table
# derives it), never a p-value magnitude; it changes no MOTCO output. A
# constant; flag and ask before changing.
STAGE_DE_CONCORDANCE_RULE: dict[str, str] = {
    "model": "factor_interaction",
    "coefficient": STAGE_DE_UNION,
    "feature_set": "exported",
    "interaction_detected": (
        "at least one exported feature (status kept) with padj <= STAGE_DE_ALPHA on either interaction coefficient "
        "of factor_interaction"
    ),
    "delta_rejected": "the primary analysis rejects delta (p <= ANALYSIS_ALPHA, ANALYSIS_REJECTION_RULE; result.json)",
    "cell": "concordant when both are true or both are false; discordant otherwise, written as a finding",
}
if not (
    STAGE_DE_MODELS[STAGE_DE_CONCORDANCE_RULE["model"]]["interaction"]
    and STAGE_DE_MODELS[STAGE_DE_CONCORDANCE_RULE["model"]]["stage_coding"] == "factor"
):
    raise AssertionError("STAGE_DE_CONCORDANCE_RULE must name the factor interaction model")
# Reproduction of the export from the counts now on disk (design D8; the export
# records no count-file hash): statuses and library sizes exact, TMM factors
# and the exported values to rounding (the export writes %.10g)
STAGE_DE_TMM_RTOL = 1e-10
STAGE_DE_VALUES_RTOL = 1e-9
# Layout under <out_dir>/<cohort>/: <omic>/<model>__<coefficient>.csv per
# tested coefficient and <omic>/<model>__fit.json per fit (the resume record),
# then the summary, the Markdown and the manifest
STAGE_DE_TABLE_TEMPLATE = "{model}__{coefficient}.csv"
STAGE_DE_FIT_TEMPLATE = "{model}__fit.json"
STAGE_DE_SUMMARY = "stage_de_summary.csv"
STAGE_DE_MARKDOWN = "stage_de.md"
STAGE_DE_MANIFEST = "stage_de_manifest.json"
STAGE_DE_TABLE_COLUMNS = ("feature", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "exported")
STAGE_DE_SUMMARY_COLUMNS = (
    "cohort", "omic", "model", "coefficient", "kind", "design_column", "feature_set", "n_tested", "n_padj_missing",
    "n_significant", "min_padj", "interaction_detected", "delta_rejected", "concordance",
)
_STAGE_DE_OMICS = ("rna", "atac")
_STAGE_DE_PSEUDOBULK_INPUTS = ("rna_counts.h5ad", "atac_counts.h5ad", "donor_table.csv", "manifest.json")
_STAGE_DE_EXPORT_INPUTS = (
    "metadata.csv", "covariates.csv", "features_rna.csv", "features_atac.csv", "manifest.json", "rna.csv", "atac.csv",
)
_STAGE_DE_IGNORED_WARNINGS = (r"Setting an item of incompatible dtype",)
# The parts of a fit record the summary and the manifest read; a record
# lacking one (an earlier version's, or a damaged one) is stale
_STAGE_DE_RECORD_KEYS = (
    ("tables",), ("export_reproduction", "variance_max_rel_diff"), ("fit", "warnings"), ("fit", "n_cpus"),
    ("fit", "convergence"), ("fit", "dispersion_trend", "type"), ("fit", "settings", "fit_type_used"),
    ("fit", "coefficients"), ("fit", "cooks"), ("fit", "tested"),
)
def _stage_de_warnings_scope(n_cpus: int) -> str:
    """What a fit record's ``warnings`` cover: with ``n_cpus > 1`` pydeseq2's worker processes raise their own."""
    if n_cpus == 1:
        return "all warnings of the fit (n_cpus = 1)"
    return (
        f"main process only (n_cpus = {n_cpus}): warnings raised in pydeseq2's worker processes, such as the IRLS "
        "fallback's RuntimeWarnings, are not captured; the tables do not depend on n_cpus"
    )


_STAGE_DE_PADJ_NOTE = (
    "pydeseq2's Benjamini-Hochberg over the tested features that pass its independent filter (a baseMean cutoff "
    "chosen per coefficient at alpha); filtered features have no padj and count as not significant; never "
    "recomputed on the exported subset"
)


def stage_de_formula(model: str, omic: str) -> str:
    """The full design formula of one registered model for one omic: the design terms, then the omic's covariates."""
    if model not in STAGE_DE_MODELS:
        raise ValueError(f"unknown stage-DE model {model!r}; expected one of {tuple(STAGE_DE_MODELS)}")
    if omic not in STAGE_DE_COVARIATES:
        raise ValueError(f"unknown omic {omic!r}; expected one of {tuple(STAGE_DE_COVARIATES)}")
    return " + ".join([STAGE_DE_MODELS[model]["formula"], *STAGE_DE_COVARIATES[omic]])


def stage_de_design(metadata: pd.DataFrame, omic: str, *, source: str | Path = "covariates.csv") -> pd.DataFrame:
    """The donors × design-variable table of the stage-DE models for one omic.

    Parameters
    ----------
    metadata:
        ``load_analysis_inputs(...)["metadata"]``: ``donor_id``, ``sex``
        (``Female`` / ``Male``), ``stage`` (``0`` / ``1`` / ``2``) and the
        analysis covariates (the log10 nuclei counts derived).
    omic:
        ``rna`` or ``atac``; selects ``STAGE_DE_COVARIATES[omic]``.
    source:
        Named in the error for a missing covariate value.

    Returns
    -------
    pandas.DataFrame
        Indexed by ``donor_id`` in the input order: ``sex`` and ``stage`` as
        categoricals whose first category is the reference level of
        ``STAGE_DE_REFERENCE_LEVELS``, ``stage_num`` (0, 1, 2 as float) and the
        omic's covariates as floats on their recorded scale.

    Raises
    ------
    ValueError
        A sex or stage label outside the expected levels, a level without
        donors, a covariate column absent, or a missing covariate value (donors
        and column named; nothing is imputed).
    """
    if omic not in STAGE_DE_COVARIATES:
        raise ValueError(f"unknown omic {omic!r}; expected one of {tuple(STAGE_DE_COVARIATES)}")
    donors = metadata["donor_id"].astype(str).to_numpy()
    design = pd.DataFrame(index=pd.Index(donors, name="donor_id"))
    for var, levels in _STAGE_DE_LEVELS.items():
        values = metadata[var].astype(str).to_numpy()
        unknown = sorted(set(values) - set(levels))
        if unknown:
            raise ValueError(f"stage-de: {var} label(s) {unknown} are not among the levels {list(levels)}")
        empty = [lv for lv in levels if lv not in set(values)]
        if empty:
            raise ValueError(f"stage-de: no donor has {var} = {empty}; every model needs every level")
        design[var] = pd.Categorical(values, categories=list(levels))
    design["stage_num"] = metadata["stage"].astype(str).astype(int).to_numpy(dtype=float)
    for col in STAGE_DE_COVARIATES[omic]:
        if col not in metadata.columns:
            raise ValueError(f"stage-de: the metadata table lacks covariate {col!r} ({source})")
        values = pd.to_numeric(metadata[col], errors="coerce").to_numpy(dtype=float)
        bad = donors[~np.isfinite(values)].tolist()
        if bad:
            raise ValueError(
                f"stage-de: covariate {col!r} is missing for donor(s) {bad} ({source}); nothing is imputed; nothing "
                "was written"
            )
        design[col] = values
    return design


def _resolve_design_column(columns: list[str], expected: str, *, model: str, formula: str) -> str:
    """The design-matrix column of a tested coefficient, found by name, never by position."""
    if expected in columns:
        return expected
    raise ValueError(
        f"stage-de: the {model} design matrix ({formula}) has no column {expected!r} (columns: {columns}); the tested "
        "coefficient cannot be resolved by name"
    )


def stage_de_design_check(design: pd.DataFrame, model: str, omic: str) -> dict:
    """Build one model's design matrix as pydeseq2 does, check its rank and resolve the tested columns.

    Run before anything is written, so a design that is not full rank or a
    coefficient whose column cannot be found by name stops the step before
    any fit (``stage_de_fit`` repeats both checks on the fitted matrix).

    Returns
    -------
    dict
        ``formula``, ``design_columns``, ``rank`` and ``resolved``
        (coefficient label → design column).

    Raises
    ------
    ValueError
        The design matrix is not full rank, or a tested column is absent (the
        name and the available columns are given).
    """
    from formulaic_contrasts import FormulaicContrasts

    formula = stage_de_formula(model, omic)
    matrix = FormulaicContrasts(design, formula).design_matrix
    columns = [str(c) for c in matrix.columns]
    rank = int(np.linalg.matrix_rank(np.asarray(matrix, dtype=float)))
    if rank < len(columns):
        raise ValueError(
            f"stage-de: the {model} design ({formula}) has rank {rank} < {len(columns)} columns for these "
            f"{len(design)} donors; the model cannot be fitted; nothing was written"
        )
    resolved = {
        label: _resolve_design_column(columns, expected, model=model, formula=formula)
        for label, expected in STAGE_DE_MODELS[model]["tested"]
    }
    return {"formula": formula, "design_columns": columns, "rank": rank, "resolved": resolved}


def stage_de_fit(
    counts: pd.DataFrame,
    design: pd.DataFrame,
    model: str,
    *,
    omic: str,
    exported: Iterable[bool] | None = None,
    n_cpus: int = 1,
) -> dict:
    """Fit one registered stage-DE model with pydeseq2 and Wald-test its coefficients.

    Parameters
    ----------
    counts:
        Donors × features raw counts (integer-valued), rows in ``design``
        order, columns the tested feature set in the export's order. Size
        factors and dispersions are fitted on exactly these features.
    design:
        ``stage_de_design`` output for the same donors.
    model:
        A key of ``STAGE_DE_MODELS``.
    omic:
        Selects the covariates of the formula (``stage_de_formula``).
    exported:
        Per feature (``counts`` column order), whether it is one of MOTCO's
        exported features; written as the ``exported`` column. ``None``
        leaves the column out.
    n_cpus:
        pydeseq2's worker count (runtime only; recorded).

    Returns
    -------
    dict
        ``tables`` (coefficient label → DataFrame with ``feature``,
        ``baseMean``, ``log2FoldChange``, ``lfcSE``, ``stat``, ``pvalue``,
        ``padj`` and ``exported``, one row per feature in ``counts`` order),
        ``fit`` (formula, design-matrix columns, each tested coefficient's
        resolved column, size factors, dispersion trend requested and used,
        Cook's counts, non-converged dispersion and LFC fits, per-coefficient
        missing-padj and significant counts, the warnings pydeseq2 raised,
        settings, runtime) and ``dispersions`` (per feature: gene-wise, trend,
        MAP and final dispersion and the three convergence flags).

    Raises
    ------
    ValueError
        ``counts`` and ``design`` list different donors, a count is negative
        or not an integer, the design matrix is not full rank, or a tested
        coefficient's design column cannot be found by name (the expected name
        and the available columns are given).

    Notes
    -----
    Each tested coefficient is a unit contrast vector on the design-matrix
    column found by name, with ``DeseqStats(alpha=STAGE_DE_ALPHA)`` and
    pydeseq2's defaults (Cook's refit, Cook's filtering, independent
    filtering); ``log2FoldChange`` is the coefficient on the log2 scale.
    With continuous covariates no donor shares its design row with another,
    so pydeseq2 (as R's DESeq2) replaces no outlier and filters no Cook's
    outlier; the counts are recorded.
    """
    import warnings

    if model not in STAGE_DE_MODELS:
        raise ValueError(f"unknown stage-DE model {model!r}; expected one of {tuple(STAGE_DE_MODELS)}")
    spec = STAGE_DE_MODELS[model]
    formula = stage_de_formula(model, omic)
    if counts.index.astype(str).tolist() != design.index.astype(str).tolist():
        raise ValueError("stage-de: the count rows and the design rows list different donors (or orders)")
    values = counts.to_numpy()
    if not np.issubdtype(values.dtype, np.integer):
        values = np.asarray(values, dtype=float)
        if not np.all(np.isfinite(values)) or not np.array_equal(values, np.round(values)):
            raise ValueError("stage-de: the counts are not integer-valued")
    if (values < 0).any():
        raise ValueError("stage-de: the counts contain negative values")
    feature_names = counts.columns.astype(str)
    counts_df = pd.DataFrame(values.astype(np.int64), index=design.index.astype(str), columns=feature_names)
    meta = design.copy()
    meta.index = meta.index.astype(str)
    mask = None if exported is None else np.asarray(list(exported), dtype=bool)
    if mask is not None and mask.size != counts_df.shape[1]:
        raise ValueError(f"stage-de: exported has {mask.size} entries for {counts_df.shape[1]} features")

    t0 = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        for pattern in _STAGE_DE_IGNORED_WARNINGS:
            warnings.filterwarnings("ignore", message=pattern)
        dds = DeseqDataSet(
            counts=counts_df, metadata=meta, design=formula, refit_cooks=True, n_cpus=n_cpus, quiet=True,
        )
        if any("not full rank" in str(w.message) for w in caught):
            raise ValueError(
                f"stage-de: the {model} design ({formula}) is not full rank for these {len(meta)} donors; the model "
                "cannot be fitted"
            )
        dds.deseq2()
        columns = [str(c) for c in dds.obsm["design_matrix"].columns]
        resolved = {
            label: _resolve_design_column(columns, expected, model=model, formula=formula)
            for label, expected in spec["tested"]
        }
        tables = {}
        for label, column in resolved.items():
            contrast = np.zeros(len(columns))
            contrast[columns.index(column)] = 1.0
            ds = DeseqStats(
                dds, contrast=contrast, alpha=STAGE_DE_ALPHA, cooks_filter=True, independent_filter=True,
                quiet=True, n_cpus=n_cpus,
            )
            ds.summary()
            res = ds.results_df.reindex(feature_names)
            table = pd.DataFrame({"feature": np.asarray(feature_names, dtype=object)})
            for col in STAGE_DE_TABLE_COLUMNS[1:-1]:
                table[col] = res[col].to_numpy(dtype=float)
            if mask is not None:
                table["exported"] = mask
            tables[label] = table
    runtime = time.perf_counter() - t0
    messages = sorted({f"{w.category.__name__}: {w.message}" for w in caught})

    uns = dds.uns
    trend_type = uns.get("disp_function_type")
    trend = {"type": trend_type}
    if trend_type == "parametric":
        trend["coefficients"] = {str(k): float(v) for k, v in uns["trend_coeffs"].items()}
    elif trend_type == "mean":
        trend["mean_dispersion"] = float(uns["mean_disp"])

    def _count(frame, key) -> int:
        return int(np.asarray(frame[key], dtype=bool).sum()) if key in frame else 0

    def _not_converged(key) -> int:
        if key not in dds.var:
            return 0
        flags = pd.to_numeric(dds.var[key], errors="coerce").to_numpy(dtype=float)
        return int((np.isfinite(flags) & (flags == 0)).sum())

    var = dds.var.reindex(feature_names)
    dispersions = pd.DataFrame({"feature": np.asarray(feature_names, dtype=object)})
    for out_col, key in (
        ("genewise", "genewise_dispersions"), ("trend", "fitted_dispersions"), ("map", "MAP_dispersions"),
        ("dispersion", "dispersions"),
    ):
        dispersions[out_col] = pd.to_numeric(var[key], errors="coerce").to_numpy(dtype=float) if key in var else np.nan
    for out_col, key in (
        ("genewise_converged", "_genewise_converged"), ("map_converged", "_MAP_converged"),
        ("lfc_converged", "_LFC_converged"),
    ):
        dispersions[out_col] = pd.to_numeric(var[key], errors="coerce").to_numpy(dtype=float) if key in var else np.nan

    fit = {
        "model": model,
        "omic": omic,
        "formula": formula,
        "stage_coding": spec["stage_coding"],
        "interaction": spec["interaction"],
        "design_columns": columns,
        "tested": [
            {"coefficient": label, "expected_column": expected, "design_column": resolved[label]}
            for label, expected in spec["tested"]
        ],
        "n_donors": int(counts_df.shape[0]),
        "n_features": int(counts_df.shape[1]),
        "n_exported": None if mask is None else int(mask.sum()),
        "size_factors": {str(d): float(v) for d, v in dds.obs["size_factors"].items()},
        "dispersion_trend": trend,
        "cooks": {
            "replaceable_donors": _count(dds.obs, "replaceable"),
            "replaced_features": _count(dds.var, "replaced"),
            "refitted_features": _count(dds.var, "refitted"),
            "outlier_features": int(np.asarray(dds.cooks_outlier(), dtype=bool).sum()),
            "note": (
                "Cook's refit and filtering need donors that share a design row (>= 7 for the refit, >= 3 for the "
                "filter); with continuous covariates none does, so both are inactive by construction (as in R's "
                "DESeq2) and the per-feature results carry no outlier protection"
            ),
        },
        "convergence": {
            "genewise_not_converged": _not_converged("_genewise_converged"),
            "map_not_converged": _not_converged("_MAP_converged"),
            "lfc_not_converged": _not_converged("_LFC_converged"),
            "note": (
                "pydeseq2's L-BFGS-B dispersion fits fall back to a grid search and its IRLS to an optimizer when "
                "they do not converge; these features are the ones least reproducible across implementations"
            ),
        },
        "coefficients": {
            label: {
                "n_pvalue_missing": int(t["pvalue"].isna().sum()),
                "n_padj_missing": int(t["padj"].isna().sum()),
                "n_significant": int((t["padj"] <= STAGE_DE_ALPHA).sum()),
            }
            for label, t in tables.items()
        },
        "settings": {
            "alpha": STAGE_DE_ALPHA,
            "padj": _STAGE_DE_PADJ_NOTE,
            "independent_filter": True,
            "cooks_filter": True,
            "refit_cooks": True,
            "min_replicates": int(dds.min_replicates),
            "fit_type_requested": str(dds.fit_type),
            "fit_type_used": trend_type,
            "size_factors_fit_type": str(dds.size_factors_fit_type),
            "lfc_shrinkage": False,
        },
        "n_cpus": n_cpus,
        "warnings": messages,
        "warnings_scope": _stage_de_warnings_scope(n_cpus),
        "runtime_seconds": runtime,
    }
    del dds
    gc.collect()
    return {"tables": tables, "fit": fit, "dispersions": dispersions}


def _stage_de_export_settings(manifest: dict, path: Path) -> None:
    """Stop unless the export was made with the primary settings (the reproduction reruns those)."""
    problems = []
    if manifest.get("variant"):
        problems.append(f"variant {manifest['variant'].get('name', manifest['variant'])!r}")
    if manifest.get("excluded_gene_sets"):
        problems.append(f"excluded_gene_sets {manifest['excluded_gene_sets']}")
    for o in _STAGE_DE_OMICS:
        rule = (manifest.get("selection_rule") or {}).get(o, "variance")
        if rule != "variance":
            problems.append(f"selection_rule {o} = {rule!r}")
        norm = (manifest.get("normalization") or {}).get(o, "none")
        if norm != "none":
            problems.append(f"normalization {o} = {norm!r}")
    if manifest.get("residualization"):
        problems.append("residualization")
    if problems:
        raise ValueError(
            f"{path} records {', '.join(problems)}; stage-de binds to an export made with the primary settings "
            "(chromosome and abundance filters, variance selection, no normalization or residualization); "
            "nothing was written"
        )


def _stage_de_inputs(cohort: str, pseudobulk_dir: Path, input_dir: Path, analysis_dir: Path) -> dict:
    """Check the input files exist (naming the producing command); return group → name → path."""
    groups = {
        "pseudobulk": (
            pseudobulk_dir, _STAGE_DE_PSEUDOBULK_INPUTS,
            "run `sead-sexdiff pseudobulk` on the cluster and sync results/pseudobulk/ (scripts/sync_pseudobulk.sh) "
            "to produce it",
        ),
        "export": (input_dir / cohort, _STAGE_DE_EXPORT_INPUTS, "run `sead-sexdiff export` to produce it"),
        "analysis": (analysis_dir / cohort, _REFERENCE_FILES, "run `sead-sexdiff analyze` to produce it"),
    }
    paths: dict[str, dict[str, Path]] = {}
    for group, (d, names, hint) in groups.items():
        paths[group] = {}
        for name in names:
            p = d / name
            if not p.exists():
                raise FileNotFoundError(f"{p} is missing; {hint}; nothing was written")
            paths[group][name] = p
    return paths


def _stage_de_reproduce(
    counts: np.ndarray,
    features: pd.DataFrame,
    exported_values: pd.DataFrame,
    export_manifest: dict,
    *,
    omic: str,
    donors: list[str],
    counts_path: Path,
    export_dir: Path,
) -> dict:
    """Rerun the export's preprocessing on the cohort counts and compare with what the export recorded.

    Raises a ``ValueError`` naming ``counts_path`` and the ``export`` command
    when a status, a library size, a TMM factor, a tested feature's log2-CPM
    variance or an exported value differs.
    """
    rerun = (
        "the counts on disk do not reproduce the primary export; they changed since the export; rerun "
        "`sead-sexdiff export` (then `sead-sexdiff analyze`); nothing was written"
    )
    thr = export_manifest["thresholds"]
    # preprocess_block takes k, min_cpm, prior_count, N and the chromosomes from
    # the record; the TMM trims and the float format are the module constants
    for key, recorded, constant in (
        ("thresholds.tmm.logratio_trim", (thr.get("tmm") or {}).get("logratio_trim"), TMM_LOGRATIO_TRIM),
        ("thresholds.tmm.sum_trim", (thr.get("tmm") or {}).get("sum_trim"), TMM_SUM_TRIM),
        ("thresholds.float_format", thr.get("float_format"), EXPORT_FLOAT_FORMAT),
    ):
        if recorded != constant:
            raise ValueError(
                f"{export_dir / 'manifest.json'} records {key} = {recorded!r} but the export constant is "
                f"{constant!r}; the reproduction reruns the export with the constants; rerun `sead-sexdiff export` "
                "(then `sead-sexdiff analyze`); nothing was written"
            )
    names = features["feature"].to_numpy(dtype=str)
    res = preprocess_block(
        counts, names, features["chromosome"].to_numpy(dtype=object), k=int(thr["k"]), min_cpm=float(thr["min_cpm"]),
        prior_count=float(thr["prior_count"]), top_n=int(thr["N"]), excluded_chromosomes=tuple(thr["excluded_chromosomes"]),
    )
    status_now = np.asarray(res["status"], dtype=str)
    status_then = features["status"].to_numpy(dtype=str)
    diff = np.flatnonzero(status_now != status_then)
    if diff.size:
        examples = ", ".join(f"{names[i]} ({status_then[i]} recorded, {status_now[i]} now)" for i in diff[:3])
        raise ValueError(
            f"{counts_path}: {diff.size} feature status(es) differ from {export_dir / f'features_{omic}.csv'} for the "
            f"cohort donors (e.g. {examples}); {rerun}"
        )
    tmm = export_manifest.get("tmm", {}).get(omic, {})
    absent = [d for d in donors if d not in tmm]
    if absent:
        raise ValueError(f"{export_dir / 'manifest.json'} records no TMM entry for {omic} donor(s) {absent}; {rerun}")
    for key, now in (("lib_size_filter", res["lib_size_filter"]), ("lib_size", res["lib_size"])):
        then = np.array([tmm[d][key] for d in donors], dtype=float)
        if not np.array_equal(np.asarray(now, dtype=float), then):
            bad = [d for d, a, b in zip(donors, now, then) if a != b]
            raise ValueError(
                f"{counts_path}: the {key} of donor(s) {bad[:5]} differs from {export_dir / 'manifest.json'}; {rerun}"
            )
    factors_then = np.array([tmm[d]["factor"] for d in donors], dtype=float)
    factor_rel = float(np.max(np.abs(res["tmm_factors"] - factors_then) / np.abs(factors_then)))
    if not factor_rel <= STAGE_DE_TMM_RTOL:
        raise ValueError(
            f"{counts_path}: the TMM factors differ from {export_dir / 'manifest.json'} by up to {factor_rel:.3g} "
            f"(relative; tolerance {STAGE_DE_TMM_RTOL}); {rerun}"
        )
    # Every tested feature's log2-CPM variance (features_<omic>.csv), so a count
    # change confined to tested-but-not-exported features is caught as well
    tested = np.isin(status_then, STAGE_DE_FEATURE_STATUSES)
    variance_then = pd.to_numeric(features["variance"], errors="coerce").to_numpy(dtype=float)[tested]
    variance_now = np.asarray(res["variance"], dtype=float)[tested]
    if not np.all(np.isfinite(variance_then)):
        raise ValueError(f"{export_dir / f'features_{omic}.csv'} lacks the variance of a tested feature; {rerun}")
    variance_diff = np.abs(variance_now - variance_then)
    variance_bad = ~(variance_diff <= STAGE_DE_VALUES_RTOL * np.abs(variance_then) + 1e-12)
    if variance_bad.any():
        examples = ", ".join(names[tested][variance_bad][:3])
        raise ValueError(
            f"{counts_path}: the log2-CPM variance of {int(variance_bad.sum())} tested feature(s) differs from "
            f"{export_dir / f'features_{omic}.csv'} (e.g. {examples}; up to {float(np.max(variance_diff)):.3g}); {rerun}"
        )
    variance_rel = float(np.max(variance_diff / np.maximum(np.abs(variance_then), 1e-300))) if variance_diff.size else 0.0
    kept = list(res["kept"])
    if kept != exported_values.columns.astype(str).tolist():
        raise ValueError(f"{counts_path}: the reproduced kept features differ from {export_dir / f'{omic}.csv'}; {rerun}")
    then_values = exported_values.loc[donors].to_numpy(dtype=float)
    abs_diff = np.abs(res["values"] - then_values)
    value_rel = float(np.max(abs_diff / np.maximum(np.abs(then_values), 1e-300)))
    value_abs = float(np.max(abs_diff))
    if not np.all(abs_diff <= STAGE_DE_VALUES_RTOL * np.abs(then_values) + 1e-12):
        raise ValueError(
            f"{counts_path}: the reproduced log2-CPM values differ from {export_dir / f'{omic}.csv'} by up to "
            f"{value_abs:.3g} (absolute; tolerance {STAGE_DE_VALUES_RTOL} relative); {rerun}"
        )
    return {
        "statuses_equal": True,
        "lib_sizes_equal": True,
        "tmm_factor_max_rel_diff": factor_rel,
        "tmm_rtol": STAGE_DE_TMM_RTOL,
        "variance_max_rel_diff": variance_rel,
        "values_max_abs_diff": value_abs,
        "values_max_rel_diff": value_rel,
        "values_rtol": STAGE_DE_VALUES_RTOL,
        "thresholds": {
            "excluded_chromosomes": list(thr["excluded_chromosomes"]), "min_cpm": thr["min_cpm"], "k": thr["k"],
            "prior_count": thr["prior_count"], "N": thr["N"],
        },
        "steps": res["steps"],
        "note": (
            "preprocess_block rerun on the cohort donors' counts with the chromosomes of features_<omic>.csv and the "
            "recorded thresholds; statuses and library sizes equal exactly, the TMM factors, every tested feature's "
            "log2-CPM variance and the exported values to rounding. Not a hash: a count change that leaves all of "
            "these unchanged would pass"
        ),
    }


def _stage_de_prepare(
    cohort: str,
    *,
    pseudobulk_dir: Path,
    input_dir: Path,
    analysis_dir: Path,
    omics: list[str],
    motco_version: str,
    log: str,
    cell_type: str = DEFAULT_CELL_TYPE,
) -> dict:
    """Every check of one cohort before anything is written; the loaded counts, designs and references.

    The export manifest must record ``cell_type`` (:func:`check_manifest_cell_type`).

    Returns
    -------
    dict
        ``paths``, ``hashes`` (group → name → sha256 of every input file),
        ``cell_type_findings``,
        ``loaded`` (``load_analysis_inputs``), ``reference`` (``_check_reference``),
        ``delta_rejected``, ``donors``, ``cells``, ``k``, ``features`` (omic →
        ``features_<omic>.csv``), ``feature_sets`` (omic → abundance / exported
        counts), ``designs`` and, for the requested omics, ``counts`` (donors ×
        tested features), ``exported`` (mask) and ``reproduction``.
    """
    c = cohort
    pb, ex_dir, an_dir = pseudobulk_dir, input_dir / c, analysis_dir / c
    paths = _stage_de_inputs(c, pb, input_dir, analysis_dir)
    rerun_export = "rerun `sead-sexdiff export` (then `sead-sexdiff analyze`); nothing was written"
    # Hash every file before it is read (load_analysis_inputs hashes the five
    # analysis inputs of the export itself before reading them)
    hashes = {
        group: {
            name: _sha256_file(p) for name, p in files.items() if not (group == "export" and name in _ANALYSIS_INPUTS)
        }
        for group, files in paths.items()
    }
    pb_manifest = json.loads((pb / "manifest.json").read_text())
    if pb_manifest.get("partial", True):
        raise ValueError(f"{pb / 'manifest.json'}: the pseudobulk run is partial; stage-de needs a full run")
    info = pb_manifest.get("cohorts", {}).get(c)
    if info is None or not info.get("analysable", False):
        raise ValueError(f"{pb / 'manifest.json'}: cohort {c} is missing or not analysable")

    loaded = load_analysis_inputs(ex_dir)
    cell_type_findings = check_manifest_cell_type(loaded["export_manifest"], ex_dir / "manifest.json", cell_type)
    hashes["export"].update({name: info["sha256"] for name, info in loaded["inputs"].items()})
    meta = loaded["metadata"]
    # The designs first, so a missing covariate is reported by donor and column
    designs = {o: stage_de_design(meta, o, source=ex_dir / "covariates.csv") for o in _STAGE_DE_OMICS}
    design_checks = {o: {m: stage_de_design_check(designs[o], m, o) for m in STAGE_DE_MODELS} for o in _STAGE_DE_OMICS}
    reference = _check_reference(ex_dir, an_dir, loaded["inputs"], motco_version)
    for name, info in reference["files"].items():
        if info["sha256"] != hashes["analysis"][name]:
            raise ValueError(f"{an_dir / name} changed while it was read; nothing was written")
    export_manifest = loaded["export_manifest"]
    recorded = export_manifest.get("input", {}).get("pseudobulk_manifest_sha256")
    now = hashes["pseudobulk"]["manifest.json"]
    if recorded != now:
        raise ValueError(
            f"{pb / 'manifest.json'} has sha256 {now} but {ex_dir / 'manifest.json'} recorded {recorded} "
            f"(input.pseudobulk_manifest_sha256): the export was made from a different pseudobulk; {rerun_export}"
        )
    _stage_de_export_settings(export_manifest, ex_dir / "manifest.json")

    donors = meta["donor_id"].astype(str).tolist()
    table = pd.read_csv(pb / "donor_table.csv", dtype={"donor_id": str, "stage": str, "sex": str})
    if c not in table.columns:
        raise ValueError(f"{pb / 'donor_table.csv'} has no {c} column")
    sub = table[table[c].astype(bool)].set_index("donor_id")
    if set(sub.index) != set(donors):
        raise ValueError(
            f"{pb / 'donor_table.csv'} and {ex_dir / 'metadata.csv'} list different {c} donors (only in the donor "
            f"table: {sorted(set(sub.index) - set(donors))}; only in the export: {sorted(set(donors) - set(sub.index))}); "
            f"{rerun_export}"
        )
    for col in ("sex", "stage"):
        bad = [d for d, v in zip(donors, meta[col].astype(str)) if str(sub.loc[d, col]) != v]
        if bad:
            raise ValueError(
                f"{pb / 'donor_table.csv'} and {ex_dir / 'metadata.csv'} disagree on {col} for donor(s) {bad}; "
                f"{rerun_export}"
            )
    k, _, cells = smallest_cell(sub.reset_index())
    thr = export_manifest.get("thresholds", {})
    if thr.get("k") != k or export_manifest.get("metadata", {}).get("cells") != cells:
        raise ValueError(
            f"{ex_dir / 'manifest.json'} records k = {thr.get('k')} and cells "
            f"{export_manifest.get('metadata', {}).get('cells')} but {pb / 'donor_table.csv'} gives k = {k} and cells "
            f"{cells}; {rerun_export}"
        )

    features, feature_sets = {}, {}
    for o in _STAGE_DE_OMICS:
        path = ex_dir / f"features_{o}.csv"
        f = pd.read_csv(path, dtype={"feature": str, "chromosome": str, "status": str}, keep_default_na=False)
        unknown = sorted(set(f["status"]) - set(_STAGE_DE_EXPORT_STATUSES))
        if unknown:
            raise ValueError(f"{path} has status(es) {unknown}; stage-de reads a primary export ({_STAGE_DE_EXPORT_STATUSES})")
        if not f["feature"].is_unique:
            raise ValueError(f"{path}: feature names repeat")
        tested = f["status"].isin(STAGE_DE_FEATURE_STATUSES).to_numpy()
        kept = (f["status"] == STAGE_DE_EXPORTED_STATUS).to_numpy()
        steps = export_manifest.get("feature_counts", {}).get(o, {})
        expected = {
            "input": len(f), "after_sex_chromosome": int((f["status"] != "sex_chromosome").sum()),
            "after_abundance": int(tested.sum()), "after_variance": int(kept.sum()),
        }
        for key, val in expected.items():
            if steps.get(key) != val:
                raise ValueError(
                    f"{path} gives {key} = {val} but {ex_dir / 'manifest.json'} records {steps.get(key)}; {rerun_export}"
                )
        if f.loc[kept, "feature"].tolist() != loaded["blocks"][o].columns.astype(str).tolist():
            raise ValueError(f"{path}: the kept features differ from the columns of {ex_dir / f'{o}.csv'}; {rerun_export}")
        features[o] = f
        feature_sets[o] = {"abundance": int(tested.sum()), "exported": int(kept.sum())}

    counts, exported, reproduction = {}, {}, {}
    for o in omics:
        path = pb / f"{o}_counts.h5ad"
        f = features[o]
        a = ad.read_h5ad(path)
        names = a.var_names.astype(str).to_numpy()
        if not np.array_equal(names, f["feature"].to_numpy(dtype=str)):
            raise ValueError(
                f"{path} lists features other than, or in another order than, {ex_dir / f'features_{o}.csv'}; "
                f"{rerun_export}"
            )
        obs = a.obs.copy()
        obs.index = obs.index.astype(str)
        absent = sorted(set(donors) - set(obs.index))
        if absent:
            raise ValueError(f"{path} has no row for {c} donor(s) {absent}; {rerun_export}")
        for col in ("sex", "stage"):
            if col not in obs.columns:
                raise ValueError(f"{path}: obs lacks {col!r}")
            bad = [d for d, v in zip(donors, meta[col].astype(str)) if str(obs.loc[d, col]) != v]
            if bad:
                raise ValueError(f"{path} and {ex_dir / 'metadata.csv'} disagree on {col} for donor(s) {bad}; {rerun_export}")
        rows = obs.index.get_indexer(donors)
        X = a.X[rows]
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        del a
        if not np.array_equal(X, np.round(X)) or (X < 0).any():
            raise ValueError(f"{path}: the counts are not non-negative integers")
        reproduction[o] = _stage_de_reproduce(
            X, f, loaded["blocks"][o], export_manifest, omic=o, donors=donors, counts_path=path, export_dir=ex_dir,
        )
        tested = f["status"].isin(STAGE_DE_FEATURE_STATUSES).to_numpy()
        counts[o] = pd.DataFrame(
            np.asarray(X[:, tested], dtype=np.int64), index=pd.Index(donors, name="donor_id"),
            columns=f.loc[tested, "feature"].to_numpy(dtype=str),
        )
        exported[o] = (f.loc[tested, "status"] == STAGE_DE_EXPORTED_STATUS).to_numpy()
        del X
        gc.collect()
        print(
            f"[{log}] {c} {o}: {path.name} reproduces the export (statuses, library sizes; TMM factors to "
            f"{reproduction[o]['tmm_factor_max_rel_diff']:.2g}, tested-feature variances to "
            f"{reproduction[o]['variance_max_rel_diff']:.2g}, values to {reproduction[o]['values_max_abs_diff']:.2g}); "
            f"{feature_sets[o]['abundance']} features tested, {feature_sets[o]['exported']} exported"
        )

    delta_rejected = _sensitivity_row(
        reference["record"], cohort=c, component_selection="cv", directory=str(an_dir)
    )["delta_rejected"]
    if delta_rejected is None:
        raise ValueError(f"{an_dir / 'result.json'} has no delta p-value; rerun `sead-sexdiff analyze`")
    return {
        "paths": paths, "hashes": hashes, "cell_type_findings": cell_type_findings, "loaded": loaded,
        "reference": reference, "delta_rejected": delta_rejected,
        "donors": donors, "cells": cells, "k": k, "features": features, "feature_sets": feature_sets,
        "designs": designs, "design_checks": design_checks, "counts": counts, "exported": exported,
        "reproduction": reproduction,
    }


def _stage_de_fingerprint(prep: dict, cohort: str, omic: str, model: str) -> dict:
    """What a fit depends on: inputs (sha256), formula, tested coefficients, covariates, levels, α, pydeseq2."""
    from importlib.metadata import PackageNotFoundError, version

    try:
        pydeseq2_version = version("pydeseq2")
    except PackageNotFoundError:
        pydeseq2_version = None
    h = prep["hashes"]
    return {
        "cohort": cohort,
        "omic": omic,
        "model": model,
        "formula": stage_de_formula(model, omic),
        "tested": [list(t) for t in STAGE_DE_MODELS[model]["tested"]],
        "covariates": list(STAGE_DE_COVARIATES[omic]),
        "reference_levels": dict(STAGE_DE_REFERENCE_LEVELS),
        "alpha": STAGE_DE_ALPHA,
        "feature_statuses": list(STAGE_DE_FEATURE_STATUSES),
        "donors": list(prep["donors"]),
        "inputs": {
            f"pseudobulk/{omic}_counts.h5ad": h["pseudobulk"][f"{omic}_counts.h5ad"],
            "pseudobulk/manifest.json": h["pseudobulk"]["manifest.json"],
            "pseudobulk/donor_table.csv": h["pseudobulk"]["donor_table.csv"],
            **{f"export/{n}": h["export"][n] for n in ("metadata.csv", "covariates.csv", f"features_{omic}.csv",
                                                      "manifest.json", f"{omic}.csv")},
            **{f"analysis/{n}": h["analysis"][n] for n in _REFERENCE_FILES},
        },
        "pydeseq2": pydeseq2_version,
    }


def _stage_de_model_state(omic_dir: Path, model: str, fingerprint: dict) -> str:
    """``complete`` (fit record and tables match), ``missing`` (to fit) or ``stale``.

    ``stale`` covers a record made from other inputs or settings, an
    unreadable record, and a record lacking a part the summary and manifest
    read (``_STAGE_DE_RECORD_KEYS``), such as one written by an earlier
    version of the step.
    """
    record_path = omic_dir / STAGE_DE_FIT_TEMPLATE.format(model=model)
    if not record_path.exists():
        return "missing"
    try:
        record = json.loads(record_path.read_text())
    except (OSError, json.JSONDecodeError, UnicodeDecodeError):
        return "stale"
    if not isinstance(record, dict) or record.get("fingerprint") != fingerprint:
        return "stale"
    for path in _STAGE_DE_RECORD_KEYS:
        node = record
        for key in path:
            if not isinstance(node, dict) or key not in node:
                return "stale"
            node = node[key]
    tables = record["tables"]
    if not isinstance(tables, dict) or not all(
        isinstance(info, dict) and isinstance(info.get("file"), str) and isinstance(info.get("sha256"), str)
        for info in tables.values()
    ):
        return "stale"
    for label, info in tables.items():
        p = omic_dir / info["file"]
        if not p.is_file() or _sha256_file(p) != info["sha256"]:
            return "missing"
    if set(tables) != {label for label, _ in STAGE_DE_MODELS[model]["tested"]}:
        return "missing"
    return "complete"


def _stage_de_summary_rows(cohort: str, omic: str, tables: dict[str, dict[str, pd.DataFrame]], delta_rejected: bool) -> list[dict]:
    """Summary rows of one omic: per model × coefficient (and the union of a two-coefficient model) × feature set."""
    rule = STAGE_DE_CONCORDANCE_RULE
    rows = []
    for model, spec in STAGE_DE_MODELS.items():
        model_tables = tables[model]
        labels = [label for label, _ in spec["tested"]]
        features = model_tables[labels[0]]["feature"].tolist()
        for label in labels[1:]:
            if model_tables[label]["feature"].tolist() != features:
                raise ValueError(f"stage-de: the {model} tables of {cohort} {omic} list different features")
        exported = model_tables[labels[0]]["exported"].to_numpy(dtype=bool)
        entries = [(label, "wald", col, model_tables[label]["padj"].to_numpy(dtype=float)[:, None])
                   for label, col in spec["tested"]]
        if len(labels) > 1:
            entries.append((
                STAGE_DE_UNION, "union", ";".join(col for _, col in spec["tested"]),
                np.column_stack([model_tables[label]["padj"].to_numpy(dtype=float) for label in labels]),
            ))
        for label, kind, column, padj in entries:
            for feature_set in STAGE_DE_FEATURE_SETS:
                keep = exported if feature_set == "exported" else np.ones(len(features), dtype=bool)
                p = padj[keep]
                finite = np.isfinite(p)
                significant = (np.where(finite, p, np.inf) <= STAGE_DE_ALPHA).any(axis=1)
                row = {
                    "cohort": cohort, "omic": omic, "model": model, "coefficient": label, "kind": kind,
                    "design_column": column, "feature_set": feature_set, "n_tested": int(keep.sum()),
                    "n_padj_missing": int((~finite).all(axis=1).sum()), "n_significant": int(significant.sum()),
                    "min_padj": float(np.min(p[finite])) if finite.any() else float("nan"),
                    "interaction_detected": None, "delta_rejected": None, "concordance": None,
                }
                if (model, label, feature_set) == (rule["model"], rule["coefficient"], rule["feature_set"]):
                    detected = row["n_significant"] > 0
                    row["interaction_detected"] = detected
                    row["delta_rejected"] = bool(delta_rejected)
                    row["concordance"] = "concordant" if detected == bool(delta_rejected) else "discordant"
                rows.append(row)
    return rows


def render_stage_de(summary: pd.DataFrame, manifest: dict) -> str:
    """Render ``stage_de.md`` from ``stage_de_summary.csv`` and ``stage_de_manifest.json``.

    Parameters
    ----------
    summary:
        The cohort's summary table (``STAGE_DE_SUMMARY_COLUMNS``).
    manifest:
        The cohort's manifest.

    Returns
    -------
    str
        Markdown: the reading rules, the concordance table, one count table
        per omic (per model × coefficient and union: significant / tested on
        the exported and the abundance-passing set, missing padj, smallest
        padj), the models with their formulas and resolved design columns, the
        covariates, the feature sets, the fits (runtime, CPUs, Cook's counts,
        dispersion trend, warnings), the inputs with their sha256 and the
        reproduction of the export, and the findings. Every number comes from
        the two arguments.
    """
    c = manifest["cohort"]
    alpha = manifest["alpha"]
    ref = manifest["reference"]
    lines = [
        f"# Pseudobulk stage DE: {c}",
        "",
        f"Per-feature DESeq2 (pydeseq2 {manifest['software'].get('pydeseq2')}) on the raw pseudobulk counts of the "
        f"{manifest['n_donors']} {c} donors, covariate-adjusted models of ADNC stage and of sex × stage, on the "
        "features that passed the primary export's chromosome and abundance filters (size factors and dispersions "
        "fitted on that set; DESeq2's own normalization, not the export's TMM). A feature is significant when "
        f"padj ≤ {alpha}: Benjamini–Hochberg over the tested features that pass pydeseq2's independent filter (a "
        "baseMean cutoff chosen per coefficient; a filtered feature has no padj and counts as not significant); no "
        "permutation null. The FDR is approximate: DESeq2's Wald test plugs in estimated dispersions, and in "
        "simulations under DESeq2's own model the realized FDR at padj ≤ 0.05 was about 0.06–0.09, higher when few "
        "features are non-null. With continuous covariates no two donors share a design row, so Cook's refit and Cook's "
        "filtering are inactive by construction (as in R's DESeq2): the per-feature results carry no outlier "
        "protection. The `exported` counts (the features MOTCO measured) are the primary reading; padj is never "
        "recomputed on that subset. A factor model's `union` counts features significant on either Wald-tested "
        "coefficient and carries no joint error-rate claim (pydeseq2 has no likelihood-ratio test). "
        "Reading rules: the stage counts describe the stage signal in the pseudobulks; the interaction counts and the "
        "concordance cells cross-check the primary analysis and are not a second test of its hypothesis; a "
        "discordance in either direction is reported as such, since per-feature FDR over thousands of features and "
        "one latent-space permutation test ask related, not identical, questions; a zero interaction count means none "
        "detected at this n, not no interaction; the primary analysis stays primary with its non-rejection rules "
        f"(`angle` and `shape` non-rejections are uninformative). The step read sex and stage, so it is a check beside "
        "the analysis, never a filter: no export, analysis or sensitivity run reads its outputs.",
        "",
        "## Concordance",
        "",
        f"Rule (`STAGE_DE_CONCORDANCE_RULE`, predeclared): `interaction_detected` = "
        f"{manifest['concordance_rule']['interaction_detected']}; `delta_rejected` = "
        f"{manifest['concordance_rule']['delta_rejected']}; {manifest['concordance_rule']['cell']}. The primary "
        f"analysis (`{ref['dir']}`, CV-selected rank {ref['selected_lv']}, MOTCO {ref['motco']}) is read, never "
        "rerun or modified; a DE result changes no MOTCO reading.",
        "",
        "| omic | interaction detected | exported features with interaction padj ≤ α | delta rejected | cell |",
        "|---|:---:|---:|:---:|---|",
    ]
    for o, cell in manifest["concordance"].items():
        lines.append(
            f"| {o} | {'yes' if cell['interaction_detected'] else 'no'} | {cell['n_interaction_exported']} | "
            f"{'yes' if cell['delta_rejected'] else 'no'} | {cell['cell']} |"
        )
    for o in manifest["omics"]:
        sub = summary[summary["omic"] == o]
        lines += [
            "",
            f"## Counts: {o}",
            "",
            f"{manifest['feature_sets'][o]['exported']} exported and {manifest['feature_sets'][o]['abundance']} "
            "abundance-passing features tested; significant / tested at padj ≤ α, missing padj (removed by the "
            "independent filter; for `union`, missing on both coefficients) and the smallest padj per set.",
            "",
            "| model | coefficient | design column | exported: significant / tested | abundance: significant / tested | "
            "padj missing (exported / abundance) | min padj (exported) | min padj (abundance) |",
            "|---|---|---|---:|---:|---:|---:|---:|",
        ]
        keys = sub[["model", "coefficient"]].drop_duplicates().itertuples(index=False)
        for model, coefficient in keys:
            pair = sub[(sub["model"] == model) & (sub["coefficient"] == coefficient)].set_index("feature_set")
            e, a = pair.loc["exported"], pair.loc["abundance"]
            label = f"{coefficient} (no joint error rate)" if coefficient == STAGE_DE_UNION else coefficient
            lines.append(
                f"| {model} | {label} | `{e['design_column']}` | {int(e['n_significant'])} / {int(e['n_tested'])} | "
                f"{int(a['n_significant'])} / {int(a['n_tested'])} | {int(e['n_padj_missing'])} / "
                f"{int(a['n_padj_missing'])} | {_fmt_nan(float(e['min_padj']), '.3g')} | "
                f"{_fmt_nan(float(a['min_padj']), '.3g')} |"
            )
    lines += [
        "",
        "## Models",
        "",
        f"Reference levels: {', '.join(f'{k} = `{v}`' for k, v in manifest['reference_levels'].items())} "
        f"(treatment coding). Primary interaction model: `{manifest['primary_interaction_model']}`.",
        "",
        "| model | stage coding | formula (RNA) | formula (ATAC) | tested coefficient → design column |",
        "|---|---|---|---|---|",
    ]
    for model, info in manifest["models"].items():
        tested = "; ".join(f"`{t['coefficient']}` → `{t['design_column']}`" for t in info["tested"])
        lines.append(
            f"| {model} | {info['stage_coding']} | `{info['formulas'].get('rna', 'n/a')}` | "
            f"`{info['formulas'].get('atac', 'n/a')}` | {tested} |"
        )
    cov = manifest["covariates"]
    lines += [
        "",
        "## Covariates",
        "",
        f"Continuous, {cov['scale']}; source `{cov['source']}`; missing values: {cov['missing_policy']} (imputed: "
        f"{'yes' if cov['imputed'] else 'no'}).",
        "",
    ]
    for o in _STAGE_DE_OMICS:
        lines.append(f"- {o}: {', '.join(f'`{x}`' for x in cov['per_omic'][o])}")
    lines += [
        "",
        "## Feature sets",
        "",
        "| omic | abundance-passing (tested) | exported (`kept`) | statuses tested |",
        "|---|---:|---:|---|",
    ]
    for o, fs in manifest["feature_sets"].items():
        lines.append(f"| {o} | {fs['abundance']} | {fs['exported']} | {', '.join(fs['statuses'])} |")
    lines += [
        "",
        "## Fits",
        "",
        "Dispersion trend used (pydeseq2 switches from the parametric fit to a mean trend when the parametric fit "
        "does not converge), Cook's counts (replaceable donors, refitted and outlier features; zero by construction "
        "here) and the features whose gene-wise or MAP dispersion or LFC fit did not converge. The warnings column "
        "counts the warnings captured in the main process; a fit made with more than one CPU does not capture those "
        "of pydeseq2's worker processes.",
        "",
        "| omic | model | features × donors | CPUs | runtime (s) | dispersion trend | Cook's replaceable / refitted / "
        "outliers | not converged (gene-wise / MAP / LFC) | warnings |",
        "|---|---|---:|---:|---:|---|---:|---:|---:|",
    ]
    for o, fits in manifest["fits"].items():
        for model, fit in fits.items():
            conv = fit.get("convergence", {})
            lines.append(
                f"| {o} | {model} | {fit['n_features']} × {fit['n_donors']} | {fit['n_cpus']} | "
                f"{_fmt(fit['runtime_seconds'], '.1f')} | {fit['dispersion_trend']['type']} | "
                f"{fit['cooks']['replaceable_donors']} / {fit['cooks']['refitted_features']} / "
                f"{fit['cooks']['outlier_features']} | {conv.get('genewise_not_converged', 'n/a')} / "
                f"{conv.get('map_not_converged', 'n/a')} / {conv.get('lfc_not_converged', 'n/a')} | "
                f"{len(fit['warnings'])} |"
            )
    lines += ["", "## Inputs", ""]
    for group, files in manifest["inputs"]["files"].items():
        for name, info in files.items():
            lines.append(f"- {group}: `{info['path']}` sha256 `{info['sha256']}`")
    lines += ["", "Reproduction of the export from the counts:", ""]
    for o, rep in manifest["export_reproduction"].items():
        lines.append(
            f"- {o}: statuses and library sizes equal; TMM factors to {rep['tmm_factor_max_rel_diff']:.3g} (relative), "
            f"tested-feature log2-CPM variances to {rep.get('variance_max_rel_diff', float('nan')):.3g} (relative), "
            f"exported values to {rep['values_max_abs_diff']:.3g} (absolute)"
        )
    lines += ["", "## Findings", ""]
    lines += [f"- {msg}" for msg in manifest["findings"]] or ["None."]
    return "\n".join(lines) + "\n"


def _write_stage_de_debug(debug_dir: Path, cohort: str, omic: str, counts: pd.DataFrame, design: pd.DataFrame) -> None:
    """Write the tested counts and the design table for ``scripts/check_stage_de_deseq2.R``."""
    debug_dir.mkdir(parents=True, exist_ok=True)
    counts.to_csv(debug_dir / f"{cohort}_{omic}_counts.csv")
    out = design.copy()
    for col in ("sex", "stage"):
        out[col] = out[col].astype(str)
    out.to_csv(debug_dir / f"{cohort}_{omic}_design.csv", float_format=_csv_float)


def run_stage_de(
    cohort: str = "all",
    *,
    omic: str = "both",
    cell_type: str = DEFAULT_CELL_TYPE,
    pseudobulk_dir: str | Path | None = None,
    input_dir: str | Path | None = None,
    analysis_dir: str | Path | None = None,
    out_dir: str | Path | None = None,
    n_cpus: int = 1,
    force: bool = False,
    debug_dir: str | Path | None = None,
) -> dict:
    """Per-feature DESeq2 of ADNC stage and sex × stage on the cohort pseudobulks, beside the primary analysis.

    Parameters
    ----------
    cohort:
        ``"cohort_30"``, ``"cohort_50"`` or ``"all"``.
    omic:
        ``"rna"``, ``"atac"`` or ``"both"``: the omics fitted in this run.
        The summary covers every omic whose fits are complete under the
        cohort directory (an omic fitted earlier from the same inputs is kept).
    cell_type:
        Registered cell type (default ``astrocyte``): the default directories
        (``cell_type_dirs``), and the cell type every export manifest read
        must record. It is recorded in the manifest, not in the fit records'
        fingerprints (their input hashes already differ between cell types).
    pseudobulk_dir:
        Directory holding ``rna_counts.h5ad``, ``atac_counts.h5ad``,
        ``donor_table.csv`` and ``manifest.json`` (read, never written).
    input_dir:
        Export root; ``<input_dir>/<cohort>/`` is the primary export (read,
        never written).
    analysis_dir:
        Analysis root; ``<analysis_dir>/<cohort>/`` is the primary analysis
        (read and checked, never rerun or written).
    out_dir:
        Output root: ``<out_dir>/<cohort>/`` with ``<omic>/<model>__<coefficient>.csv``
        (one DESeq2 table per tested coefficient), ``<omic>/<model>__fit.json``
        (the fit record that makes the step resumable), ``stage_de_summary.csv``,
        ``stage_de.md`` and ``stage_de_manifest.json``. It may not lie inside
        an input directory.
    n_cpus:
        pydeseq2's worker count; runtime only (recorded).
    force:
        Refit the requested omics of a cohort even when their fits are
        complete (their directories and the cohort's summary files are
        replaced); without it a complete fit is skipped and a fit made from
        other inputs or settings is an error.
    debug_dir:
        Python API only: write ``<cohort>_<omic>_counts.csv`` and
        ``<cohort>_<omic>_design.csv`` (the tested counts and the design
        table) and, per model, ``<cohort>_<omic>_<model>_trend.csv``,
        ``_dispersions.csv`` and ``_size_factors.csv`` there for
        ``scripts/check_stage_de_deseq2.R``. The per-model files come from the
        fits of this call, so with ``debug_dir`` every requested model must
        be fitted: a complete model without ``force`` is an error.

    Returns
    -------
    dict
        Cohort → ``{"dir", "summary", "manifest", "concordance", "findings", "fitted", "skipped"}``.

    Raises
    ------
    FileNotFoundError
        An input file is missing; the message names it and the command that
        produces it.
    ValueError
        A reference check fails (pseudobulk manifest hash, reproduction of the
        export from the counts, donors / sex / stage, k and cells, feature
        statuses, the analysis's export hashes, MOTCO version, seed,
        permutations or CV selection), a covariate value is missing, a fit
        record was made from other inputs (without ``force``), a design is not
        full rank or a coefficient cannot be resolved, or ``out_dir`` lies
        inside an input directory.

    Notes
    -----
    Every check runs for every requested cohort before anything is written.
    Fits run one model at a time; each writes its tables and then its fit
    record (``.tmp`` files renamed), so an interrupted run resumes at the first
    model without a matching record. After a cohort's fits the input files
    are hashed again and must be unchanged. The step reads sex and stage, so
    its outputs are an input of no export, analysis or sensitivity step.
    """
    import shutil
    from importlib.metadata import PackageNotFoundError, version

    log = "stage-de"
    cohorts = list(EXPORT_COHORTS) if cohort == "all" else [cohort]
    for c in cohorts:
        if c not in EXPORT_COHORTS:
            raise ValueError(f"Unknown cohort {c!r}; expected one of {EXPORT_COHORTS} or 'all'")
    omics = list(_STAGE_DE_OMICS) if omic == "both" else [omic]
    for o in omics:
        if o not in _STAGE_DE_OMICS:
            raise ValueError(f"Unknown omic {o!r}; expected one of {_STAGE_DE_OMICS} or 'both'")
    if not isinstance(n_cpus, int) or n_cpus < 1:
        raise ValueError(f"n_cpus must be a positive integer, not {n_cpus!r}")
    cell_type_entry(cell_type)
    pb = _resolve_dir(pseudobulk_dir, cell_type, "pseudobulk")
    input_dir = _resolve_dir(input_dir, cell_type, "motco_input")
    analysis_dir = _resolve_dir(analysis_dir, cell_type, "motco_analysis")
    out_dir = _resolve_dir(out_dir, cell_type, "stage_de")
    debug = None if debug_dir is None else Path(debug_dir)
    for target_name, target in (("--out-dir", out_dir), ("debug_dir", debug)):
        if target is None:
            continue
        for name, d in (("--pseudobulk-dir", pb), ("--input-dir", input_dir), ("--analysis-dir", analysis_dir)):
            if _same_dir(target, d):
                raise ValueError(
                    f"{target_name} {target} lies inside the input directory {d} ({name}); the step never writes "
                    "under its inputs; nothing was written"
                )
    motco_version = check_motco_version()
    print(f"[{log}] Cell type: {cell_type}")
    print(f"[{log}] Pseudobulk directory: {pb}")
    print(f"[{log}] Primary export directory: {input_dir}")
    print(f"[{log}] Primary analysis directory: {analysis_dir} (read and checked; it stays primary)")
    print(f"[{log}] Output directory: {out_dir}")
    print(
        f"[{log}] Models: {', '.join(STAGE_DE_MODELS)} (registry STAGE_DE_MODELS); alpha {STAGE_DE_ALPHA}; "
        f"omics fitted: {', '.join(omics)}; n_cpus {n_cpus}"
    )

    # Every check of every cohort before anything is written
    prepared = {}
    for c in cohorts:
        prep = _stage_de_prepare(
            c, pseudobulk_dir=pb, input_dir=input_dir, analysis_dir=analysis_dir, omics=omics,
            motco_version=motco_version, log=log, cell_type=cell_type,
        )
        prep["fingerprints"] = {
            o: {m: _stage_de_fingerprint(prep, c, o, m) for m in STAGE_DE_MODELS} for o in _STAGE_DE_OMICS
        }
        states = {}
        for o in _STAGE_DE_OMICS:
            for m in STAGE_DE_MODELS:
                states[(o, m)] = _stage_de_model_state(out_dir / c / o, m, prep["fingerprints"][o][m])
                if states[(o, m)] == "stale" and not (force and o in omics):
                    hint = "rerun with --force" + ("" if o in omics else f" --omic {o} (or both)")
                    raise ValueError(
                        f"{out_dir / c / o / STAGE_DE_FIT_TEMPLATE.format(model=m)} is unreadable or was made from "
                        f"other inputs or settings than the files now on disk; {hint}; nothing was written"
                    )
        prep["states"] = states
        if debug is not None and not force:
            done = [f"{o}/{m}" for o in omics for m in STAGE_DE_MODELS if states[(o, m)] == "complete"]
            if done:
                raise ValueError(
                    f"debug_dir needs every requested model fitted in this call, but {', '.join(done)} of {c} are "
                    f"complete under {out_dir / c}; pass force=True or another out_dir; nothing was written"
                )
        print(
            f"[{log}] {c}: primary analysis {prep['reference']['dir']} matches the export (hashes, MOTCO "
            f"{motco_version}, seed, permutations, CV selection); delta rejected: "
            f"{'yes' if prep['delta_rejected'] else 'no'}"
        )
        prepared[c] = prep

    software = _software_versions()
    for pkg in ("pydeseq2", "formulaic", "formulaic-contrasts", "statsmodels", "scikit-learn"):
        try:
            software[pkg] = version(pkg)
        except PackageNotFoundError:
            software[pkg] = None
    software["motco"] = motco_version

    results = {}
    for c in cohorts:
        prep = prepared[c]
        cohort_dir = out_dir / c
        if force:
            for o in omics:
                if (cohort_dir / o).exists():
                    shutil.rmtree(cohort_dir / o)
                for m in STAGE_DE_MODELS:
                    prep["states"][(o, m)] = "missing"
            for name in (STAGE_DE_SUMMARY, STAGE_DE_MARKDOWN, STAGE_DE_MANIFEST):
                if (cohort_dir / name).exists():
                    (cohort_dir / name).unlink()
        fitted, skipped = [], []
        for o in omics:
            omic_dir = cohort_dir / o
            omic_dir.mkdir(parents=True, exist_ok=True)
            if debug is not None:
                _write_stage_de_debug(debug, c, o, prep["counts"][o], prep["designs"][o])
            for m, spec in STAGE_DE_MODELS.items():
                if prep["states"][(o, m)] == "complete":
                    print(f"[{log}] {c} {o} {m}: complete fit record found; skipped (rerun with --force to refit)")
                    skipped.append((o, m))
                    continue
                out = stage_de_fit(
                    prep["counts"][o], prep["designs"][o], m, omic=o, exported=prep["exported"][o], n_cpus=n_cpus,
                )
                checked = prep["design_checks"][o][m]
                if out["fit"]["design_columns"] != checked["design_columns"]:
                    raise ValueError(
                        f"stage-de: the fitted {m} design columns {out['fit']['design_columns']} differ from the "
                        f"pre-checked ones {checked['design_columns']}"
                    )
                if debug is not None:
                    out["dispersions"].to_csv(debug / f"{c}_{o}_{m}_dispersions.csv", index=False, float_format=_csv_float)
                    pd.Series(out["fit"]["size_factors"], name="size_factor").rename_axis("donor_id").to_csv(
                        debug / f"{c}_{o}_{m}_size_factors.csv", float_format=_csv_float
                    )
                    trend = out["fit"]["dispersion_trend"]
                    pd.DataFrame([{
                        "type": trend["type"],
                        "a0": (trend.get("coefficients") or {}).get("a0"),
                        "a1": (trend.get("coefficients") or {}).get("a1"),
                        "mean_dispersion": trend.get("mean_dispersion"),
                    }]).to_csv(debug / f"{c}_{o}_{m}_trend.csv", index=False, float_format=_csv_float)
                table_files = {}
                for label, table in out["tables"].items():
                    name = STAGE_DE_TABLE_TEMPLATE.format(model=m, coefficient=label)
                    tmp = omic_dir / f"{name}.tmp"
                    table.to_csv(tmp, index=False, float_format=_csv_float)
                    os.replace(tmp, omic_dir / name)
                    table_files[label] = {"file": name, "sha256": _sha256_file(omic_dir / name)}
                record = {
                    "fingerprint": prep["fingerprints"][o][m],
                    "fit": out["fit"],
                    "export_reproduction": prep["reproduction"][o],
                    "tables": table_files,
                    "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                }
                record_path = omic_dir / STAGE_DE_FIT_TEMPLATE.format(model=m)
                tmp = record_path.with_name(record_path.name + ".tmp")
                tmp.write_text(json.dumps(record, indent=2) + "\n")
                os.replace(tmp, record_path)
                prep["states"][(o, m)] = "complete"
                fitted.append((o, m))
                print(
                    f"[{log}] {c} {o} {m}: {out['fit']['n_features']} features × {out['fit']['n_donors']} donors in "
                    f"{out['fit']['runtime_seconds']:.1f} s; padj <= {STAGE_DE_ALPHA}: "
                    + ", ".join(
                        f"{label} {info['n_significant']}" for label, info in out["fit"]["coefficients"].items()
                    )
                )
            prep["counts"].pop(o, None)
            gc.collect()

        # Omics not fitted in this run: their state again now (another run may have completed them)
        for o in _STAGE_DE_OMICS:
            if o in omics:
                continue
            for m in STAGE_DE_MODELS:
                state = _stage_de_model_state(cohort_dir / o, m, prep["fingerprints"][o][m])
                if state == "stale":
                    raise ValueError(
                        f"{cohort_dir / o / STAGE_DE_FIT_TEMPLATE.format(model=m)} became unreadable or was replaced "
                        f"by a fit from other inputs during this run; rerun with --force --omic {o}"
                    )
                prep["states"][(o, m)] = state
        included = [o for o in _STAGE_DE_OMICS if all(prep["states"][(o, m)] == "complete" for m in STAGE_DE_MODELS)]
        records, tables = {}, {}
        for o in included:
            records[o], tables[o] = {}, {}
            for m, spec in STAGE_DE_MODELS.items():
                records[o][m] = json.loads((cohort_dir / o / STAGE_DE_FIT_TEMPLATE.format(model=m)).read_text())
                tables[o][m] = {
                    label: pd.read_csv(
                        cohort_dir / o / STAGE_DE_TABLE_TEMPLATE.format(model=m, coefficient=label),
                        dtype={"feature": str}, keep_default_na=False, na_values=[""], float_precision="round_trip",
                    )
                    for label, _ in spec["tested"]
                }
        rows = []
        for o in included:
            rows += _stage_de_summary_rows(c, o, tables[o], prep["delta_rejected"])
        summary = pd.DataFrame(rows, columns=list(STAGE_DE_SUMMARY_COLUMNS))
        rule = STAGE_DE_CONCORDANCE_RULE
        concordance, findings = {}, list(prep["cell_type_findings"])
        for o in included:
            cell = summary[
                (summary["omic"] == o) & (summary["model"] == rule["model"])
                & (summary["coefficient"] == rule["coefficient"]) & (summary["feature_set"] == rule["feature_set"])
            ].iloc[0]
            concordance[o] = {
                "interaction_detected": bool(cell["interaction_detected"]),
                "n_interaction_exported": int(cell["n_significant"]),
                "delta_rejected": bool(cell["delta_rejected"]),
                "cell": str(cell["concordance"]),
            }
            if concordance[o]["cell"] == "discordant":
                if concordance[o]["interaction_detected"]:
                    what = (
                        f"{concordance[o]['n_interaction_exported']} exported {o} feature(s) have padj <= "
                        f"{STAGE_DE_ALPHA} on an interaction coefficient of {rule['model']} and the primary analysis "
                        "did not reject delta"
                    )
                else:
                    what = (
                        f"no exported {o} feature has padj <= {STAGE_DE_ALPHA} on an interaction coefficient of "
                        f"{rule['model']} and the primary analysis rejected delta"
                    )
                findings.append(
                    f"{c} {o}: discordant: {what}; reported beside both (per-feature FDR and one latent-space "
                    "permutation test ask related, not identical, questions); the primary analysis stays primary and "
                    "is unchanged"
                )
        missing = [o for o in _STAGE_DE_OMICS if o not in included]
        for o in missing:
            findings.append(f"{c} {o}: not fitted yet (run `sead-sexdiff stage-de --omic {o}`); the summary lacks it")
        for o in included:
            for m in STAGE_DE_MODELS:
                for msg in records[o][m]["fit"]["warnings"]:
                    findings.append(f"{c} {o} {m}: pydeseq2 warning: {msg}")
        fit_cpus = {
            f"{o}/{m}": records[o][m]["fit"]["n_cpus"] for o in included for m in STAGE_DE_MODELS
            if records[o][m]["fit"]["n_cpus"] != 1
        }
        if fit_cpus:
            findings.append(
                f"fits made with n_cpus != 1: {', '.join(f'{k} ({v})' for k, v in fit_cpus.items())}; pydeseq2's "
                "per-feature fits are expected not to depend on n_cpus (the report states whether this was checked), "
                "but warnings raised in pydeseq2's worker processes were not captured for these fits, so their "
                "warning records and findings are incomplete"
            )

        hashes_now = {
            group: {name: _sha256_file(p) for name, p in files.items()} for group, files in prep["paths"].items()
        }
        if hashes_now != prep["hashes"]:
            changed = [f"{g}/{n}" for g in hashes_now for n in hashes_now[g] if hashes_now[g][n] != prep["hashes"][g][n]]
            raise ValueError(f"input file(s) {changed} changed during the stage-de run of {c}")

        reference = prep["reference"]
        ref_record = reference["record"]
        manifest = {
            "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "cell_type": cell_type,
            "cohort": c,
            "n_donors": len(prep["donors"]),
            "cells": prep["cells"],
            "k": prep["k"],
            "omics": included,
            "complete": not missing,
            "inputs": {
                "files": {
                    group: {name: {"path": str(p), "sha256": prep["hashes"][group][name]} for name, p in files.items()}
                    for group, files in prep["paths"].items()
                },
                "pseudobulk_manifest_sha256": prep["hashes"]["pseudobulk"]["manifest.json"],
                "export_manifest_sha256": prep["hashes"]["export"]["manifest.json"],
                "export_input": prep["loaded"]["export_manifest"].get("input"),
                "note": (
                    "every file hashed before it was read and again after the run; the export records the pseudobulk "
                    "manifest's sha256, not the count files', so the counts are bound to the export by reproduction "
                    "(export_reproduction: statuses, library sizes, TMM factors, tested-feature log2-CPM variances, "
                    "exported values), not by a hash"
                ),
            },
            "export_reproduction": {
                o: prep["reproduction"].get(o, records[o][next(iter(STAGE_DE_MODELS))]["export_reproduction"])
                for o in included
            },
            "reference": {
                "dir": str(reference["dir"]),
                "files": reference["files"],
                "result_sha256": reference["files"]["result.json"]["sha256"],
                "delta_rejected": prep["delta_rejected"],
                "alpha": ANALYSIS_ALPHA,
                "rejection_rule": ANALYSIS_REJECTION_RULE,
                "component_selection": ref_record["latent"]["component_selection"],
                "selected_lv": ref_record["latent"]["selected_lv"],
                "motco": reference["manifest"]["software"]["motco"],
                "seed": reference["manifest"]["measurement"]["seed"],
                "permutations": reference["manifest"]["measurement"]["permutations"],
                "note": (
                    "read and checked (export hashes, MOTCO version, seed, permutations, CV selection), never rerun "
                    "or modified; only its delta rejection status enters the concordance rule"
                ),
            },
            "models": {
                m: {
                    "formula": spec["formula"],
                    "formulas": {o: stage_de_formula(m, o) for o in _STAGE_DE_OMICS},
                    "stage_coding": spec["stage_coding"],
                    "interaction": spec["interaction"],
                    "tested": [
                        {
                            "coefficient": label,
                            "design_column": (
                                records[included[0]][m]["fit"]["tested"][i]["design_column"] if included else column
                            ),
                        }
                        for i, (label, column) in enumerate(spec["tested"])
                    ],
                }
                for m, spec in STAGE_DE_MODELS.items()
            },
            "primary_interaction_model": "factor_interaction",
            "union": (
                "a factor model's two Wald tests are also counted as their union (significant on either coefficient; "
                "its padj missing = missing on both); the union carries no joint error-rate claim (pydeseq2 "
                f"{software.get('pydeseq2')} has no likelihood-ratio test)"
            ),
            "covariates": {
                "per_omic": {o: list(cols) for o, cols in STAGE_DE_COVARIATES.items()},
                "scale": (
                    "on their recorded scale, not centered or scaled (the MLE of the tested coefficients is invariant "
                    "to affine covariate transforms; features whose dispersion optimizer fell back to a grid search "
                    "may differ slightly, counted under fits.*.convergence); nuclei counts on the log10 scale as the "
                    "analysis derives them"
                ),
                "source": str(input_dir / c / "covariates.csv"),
                "missing_policy": ANALYSIS_MISSING_COVARIATES,
                "imputed": False,
            },
            "reference_levels": dict(STAGE_DE_REFERENCE_LEVELS),
            "alpha": STAGE_DE_ALPHA,
            "significance": f"padj <= alpha; padj is {_STAGE_DE_PADJ_NOTE}; no permutation null",
            "cooks": (
                "Cook's refit and filtering are on (pydeseq2's defaults) but inactive by construction: with continuous "
                "covariates no two donors share a design row, so no donor is replaceable and no Cook's outlier is "
                "filtered (as in R's DESeq2); the per-feature results carry no outlier protection"
            ),
            "design_checks": {
                o: {m: {k: v for k, v in chk.items() if k != "formula"} for m, chk in prep["design_checks"][o].items()}
                for o in _STAGE_DE_OMICS
            },
            "normalization": (
                "DESeq2 median-of-ratios size factors and dispersions fitted on the abundance-passing set; the "
                "export's TMM factors are not reused"
            ),
            "feature_sets": {
                o: {
                    **prep["feature_sets"][o],
                    "statuses": list(STAGE_DE_FEATURE_STATUSES),
                    "exported_status": STAGE_DE_EXPORTED_STATUS,
                    "source": str(input_dir / c / f"features_{o}.csv"),
                }
                for o in _STAGE_DE_OMICS
            },
            "fits": {
                o: {
                    m: {
                        **records[o][m]["fit"],
                        "warnings_scope": records[o][m]["fit"].get(
                            "warnings_scope", _stage_de_warnings_scope(records[o][m]["fit"]["n_cpus"])
                        ),
                    }
                    for m in STAGE_DE_MODELS
                }
                for o in included
            },
            "tables": {o: {m: records[o][m]["tables"] for m in STAGE_DE_MODELS} for o in included},
            "concordance_rule": dict(STAGE_DE_CONCORDANCE_RULE),
            "concordance": concordance,
            "run": {"omics_fitted": [f"{o}/{m}" for o, m in fitted], "skipped": [f"{o}/{m}" for o, m in skipped],
                    "n_cpus": n_cpus},
            "outputs": sorted([STAGE_DE_SUMMARY, STAGE_DE_MARKDOWN, STAGE_DE_MANIFEST]),
            "software": software,
            "git": _git_revision(),
            "findings": findings,
        }
        cohort_dir.mkdir(parents=True, exist_ok=True)
        tmp = cohort_dir / f"{STAGE_DE_SUMMARY}.tmp"
        summary.to_csv(tmp, index=False, float_format=_csv_float)
        os.replace(tmp, cohort_dir / STAGE_DE_SUMMARY)
        written = pd.read_csv(cohort_dir / STAGE_DE_SUMMARY, dtype={"design_column": str}, float_precision="round_trip")
        tmp = cohort_dir / f"{STAGE_DE_MARKDOWN}.tmp"
        tmp.write_text(render_stage_de(written, manifest))
        os.replace(tmp, cohort_dir / STAGE_DE_MARKDOWN)
        tmp = cohort_dir / f"{STAGE_DE_MANIFEST}.tmp"
        tmp.write_text(json.dumps(manifest, indent=2) + "\n")
        os.replace(tmp, cohort_dir / STAGE_DE_MANIFEST)
        print(f"[{log}] Wrote {cohort_dir}")
        for msg in findings:
            print(f"[{log}] FINDING: {msg}")
        results[c] = {
            "dir": cohort_dir, "summary": summary, "manifest": manifest, "concordance": concordance,
            "findings": findings, "fitted": fitted, "skipped": skipped,
        }
        del prepared[c]
        gc.collect()
    return results


# ---------------------------------------------------------------------------
# Pseudobulk stage DE, limma-voom (a second predeclared method beside stage-de)
# ---------------------------------------------------------------------------

# The method, fixed before any limma run (design D3–D5): edgeR TMM on the
# tested counts with the export's trims (library size = their column sums) →
# voomWithQualityWeights (voom's precision weights times one quality weight per
# donor from arrayWeights) → lmFit → eBayes, each with limma 3.68.4's defaults.
# The R script passes every argument it can as written here and stops when an
# installed default differs (arrayWeights' prior.n, which
# voomWithQualityWeights does not pass on, is only checked). A constant, never
# a flag, never changed from results; flag and ask before changing.
STAGE_DE_LIMMA_METHOD: dict = {
    "normalization": {
        "function": "edgeR::calcNormFactors",
        "method": "TMM",
        "logratioTrim": TMM_LOGRATIO_TRIM,
        "sumTrim": TMM_SUM_TRIM,
        "doWeighting": True,
        "Acutoff": -1e10,
        "refColumn": "automatic (upper-quartile factor closest to the mean)",
        "lib_size": "column sums of the tested counts",
    },
    "weights": {
        "function": "limma::voomWithQualityWeights",
        "normalize.method": "none",
        "span": 0.5,
        "adaptive.span": True,
        "var.design": None,
        "var.group": None,
        "method": "genebygene",
        "maxiter": 50,
        "tol": 1e-5,
        "prior.n": 10,
    },
    "fit": {"function": "limma::lmFit", "method": "ls"},
    "moderation": {
        "function": "limma::eBayes",
        "proportion": 0.01,
        "stdev.coef.lim": [0.1, 4.0],
        "trend": False,
        "robust": False,
        "legacy": None,
    },
    "adjust_method": "BH",
    "bh_family": "all tested features (limma has no independent filter; the export's abundance filter is the only one)",
}
if (
    STAGE_DE_LIMMA_METHOD["normalization"]["logratioTrim"] != TMM_LOGRATIO_TRIM
    or STAGE_DE_LIMMA_METHOD["normalization"]["sumTrim"] != TMM_SUM_TRIM
):
    raise AssertionError("STAGE_DE_LIMMA_METHOD's TMM trims must be TMM_LOGRATIO_TRIM / TMM_SUM_TRIM (the export's)")
# One joint moderated F per factor model over its two tested coefficients
# (label, (coefficient labels)); a trend model has one coefficient, whose t is
# its test. Tied to STAGE_DE_MODELS; flag and ask before changing.
STAGE_DE_LIMMA_JOINT_TESTS: dict[str, tuple[str, tuple[str, str]]] = {
    "factor_additive": ("stage_joint", ("stage_1_vs_0", "stage_2_vs_0")),
    "factor_interaction": ("interaction_joint", ("male_x_stage_1", "male_x_stage_2")),
}
if set(STAGE_DE_LIMMA_JOINT_TESTS) != {m for m, spec in STAGE_DE_MODELS.items() if len(spec["tested"]) > 1} or any(
    labels != tuple(label for label, _ in STAGE_DE_MODELS[m]["tested"])
    or joint in labels or joint == STAGE_DE_UNION
    for m, (joint, labels) in STAGE_DE_LIMMA_JOINT_TESTS.items()
) or len({joint for joint, _ in STAGE_DE_LIMMA_JOINT_TESTS.values()}) != len(STAGE_DE_LIMMA_JOINT_TESTS):
    raise AssertionError(
        "STAGE_DE_LIMMA_JOINT_TESTS must name, per factor model of STAGE_DE_MODELS, exactly its two tested labels "
        "in order, under a label of its own"
    )
# Written before any limma run (proposal); evaluated per cohort × omic beside
# the DESeq2 cell (STAGE_DE_CONCORDANCE_RULE), which stays the reported cell of
# stage-de. Reads a count and the primary analysis's delta status, never a
# p-value magnitude of MOTCO's; changes no MOTCO output. Flag and ask before changing.
STAGE_DE_LIMMA_CONCORDANCE_RULE: dict[str, str] = {
    "model": "factor_interaction",
    "test": "interaction_joint",
    "feature_set": "exported",
    "interaction_detected": (
        "at least one exported feature (status kept) with adj.P.Val <= STAGE_DE_ALPHA on the joint moderated F of "
        "factor_interaction's two interaction coefficients (interaction_joint)"
    ),
    "delta_rejected": STAGE_DE_CONCORDANCE_RULE["delta_rejected"],
    "cell": STAGE_DE_CONCORDANCE_RULE["cell"],
    "beside": (
        "reported beside the DESeq2 cell (STAGE_DE_CONCORDANCE_RULE, the union of two Wald tests), which stays the "
        "reported cell of stage-de; methods_agree compares the two interaction_detected values and a disagreement "
        "is a method-dependence finding; neither cell replaces the other or changes a MOTCO reading"
    ),
}
if not (
    STAGE_DE_LIMMA_CONCORDANCE_RULE["model"] == STAGE_DE_CONCORDANCE_RULE["model"]
    and STAGE_DE_LIMMA_CONCORDANCE_RULE["test"] == STAGE_DE_LIMMA_JOINT_TESTS[STAGE_DE_CONCORDANCE_RULE["model"]][0]
    and STAGE_DE_LIMMA_CONCORDANCE_RULE["feature_set"] == STAGE_DE_CONCORDANCE_RULE["feature_set"]
):
    raise AssertionError(
        "STAGE_DE_LIMMA_CONCORDANCE_RULE must name the factor interaction model's joint test on the exported set"
    )
# The step's TMM factors (edgeR, in R) against the export manifest's, absolute;
# the tolerance of scripts/check_tmm_edger.R. Beyond it the cohort stops.
STAGE_DE_LIMMA_TMM_ATOL = 1e-8
# Layout under <out_dir>/<cohort>/ (stage-de's): <omic>/<model>__<label>.csv per
# tested coefficient and joint test, <omic>/<model>__fit.json per fit (the
# resume record), then the summary, the method comparison, the Markdown and
# the manifest. The R script ships beside this module.
STAGE_DE_LIMMA_TABLE_TEMPLATE = "{model}__{label}.csv"
STAGE_DE_LIMMA_FIT_TEMPLATE = "{model}__fit.json"
STAGE_DE_LIMMA_SUMMARY = "stage_de_limma_summary.csv"
STAGE_DE_LIMMA_COMPARISON = "method_comparison.csv"
STAGE_DE_LIMMA_MARKDOWN = "stage_de_limma.md"
STAGE_DE_LIMMA_MANIFEST = "stage_de_limma_manifest.json"
STAGE_DE_LIMMA_SCRIPT = "stage_de_limma.R"
STAGE_DE_LIMMA_TABLE_COLUMNS = ("feature", "AveExpr", "logFC", "SE", "t", "P.Value", "adj.P.Val", "B", "exported")
# Summary rows: stage-de's columns; kind is moderated_t, union or joint_f
STAGE_DE_LIMMA_SUMMARY_COLUMNS = STAGE_DE_SUMMARY_COLUMNS
STAGE_DE_LIMMA_COMPARISON_COLUMNS = (
    "cohort", "omic", "model", "test", "kind", "deseq2_test", "limma_test", "feature_set", "n_tested",
    "n_significant_deseq2", "n_significant_limma", "n_both", "n_only_deseq2", "n_only_limma", "detected_deseq2",
    "detected_limma", "detection_agrees", "spearman_log2fc", "spearman_stat", "n_compared", "delta_rejected",
    "concordance_deseq2", "concordance_limma", "methods_agree",
)
_STAGE_DE_LIMMA_THREAD_ENV = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")


def stage_de_limma_script() -> Path:
    """The packaged R script (``stage_de_limma.R``, shipped beside this module)."""
    return Path(__file__).with_name(STAGE_DE_LIMMA_SCRIPT)


def _run_stage_de_limma_r(args: list[str], *, what: str) -> subprocess.CompletedProcess:
    """Run the packaged R script with one BLAS thread; raise with the end of R's stderr on failure."""
    import shutil

    exe = shutil.which("Rscript")
    if exe is None:
        raise RuntimeError(
            "Rscript is not on PATH; stage-de-limma runs R (limma, edgeR) from the project's Nix shell: enter it with "
            "direnv or `nix develop`, or run `nix develop -c uv run sead-sexdiff stage-de-limma ...`; nothing was written"
        )
    env = dict(os.environ)
    for name in _STAGE_DE_LIMMA_THREAD_ENV:
        env[name] = "1"
    proc = subprocess.run(
        [exe, "--vanilla", str(stage_de_limma_script()), *args], capture_output=True, text=True, env=env, check=False,
    )
    if proc.returncode != 0:
        tail = "\n".join(proc.stderr.strip().splitlines()[-15:])
        raise RuntimeError(
            f"stage-de-limma: R failed ({what}; exit status {proc.returncode}); nothing of it was written. End of R's "
            f"error output:\n{tail}"
        )
    return proc


def stage_de_limma_versions() -> dict:
    """R, limma, edgeR and jsonlite versions as the packaged script reports them, plus the script's sha256."""
    proc = _run_stage_de_limma_r(["--versions"], what="--versions")
    out = json.loads(proc.stdout)
    out["script_sha256"] = _sha256_file(stage_de_limma_script())
    return out


def _stage_de_limma_settings(model: str, omic: str, check: dict) -> dict:
    """``settings.json`` of one model: what the R script fits, tests and checks."""
    resolved = check["resolved"]
    joint = []
    if model in STAGE_DE_LIMMA_JOINT_TESTS:
        label, coefficients = STAGE_DE_LIMMA_JOINT_TESTS[model]
        joint.append({"label": label, "coefficients": list(coefficients), "columns": [resolved[c] for c in coefficients]})
    return {
        "model": model,
        "omic": omic,
        "formula": stage_de_formula(model, omic),
        "levels": {v: list(levels) for v, levels in _STAGE_DE_LEVELS.items()},
        "tested": [{"label": label, "column": resolved[label]} for label, _ in STAGE_DE_MODELS[model]["tested"]],
        "joint": joint,
        "alpha": STAGE_DE_ALPHA,
        "tmm_atol": STAGE_DE_LIMMA_TMM_ATOL,
        "method": STAGE_DE_LIMMA_METHOD,
    }


def write_stage_de_limma_exchange(
    directory: str | Path,
    counts: pd.DataFrame,
    design: pd.DataFrame,
    *,
    omic: str,
    export_tmm: pd.Series,
    models: Iterable[str] | None = None,
) -> dict[str, Path]:
    """Write one omic's exchange files for ``stage_de_limma.R``.

    Parameters
    ----------
    directory:
        The omic directory: ``counts.csv`` (features × donors, integers) and
        ``export_tmm.csv`` (``donor_id``, ``factor``) are written once there,
        and one subdirectory per model holds links to both plus
        ``design_matrix.csv``, ``design_table.csv`` (``%.17g``) and
        ``settings.json``, so each model directory can be rerun by hand.
    counts:
        Donors × tested features raw counts, rows in ``design`` order.
    design:
        ``stage_de_design`` output for the same donors.
    omic:
        Selects the covariates of the formulas.
    export_tmm:
        The export manifest's TMM factor per donor (indexed by donor id).
    models:
        Keys of ``STAGE_DE_MODELS`` (default: all).

    Returns
    -------
    dict
        Model → its directory.

    Raises
    ------
    ValueError
        Donors differ between the counts, the design and the TMM factors, a
        count is not a non-negative integer, a name holds a character the CSV
        exchange cannot carry, or a design fails ``stage_de_design_check``.
    """
    import shutil

    from formulaic_contrasts import FormulaicContrasts

    d = Path(directory)
    models = list(STAGE_DE_MODELS) if models is None else list(models)
    donors = design.index.astype(str).tolist()
    if counts.index.astype(str).tolist() != donors:
        raise ValueError("stage-de-limma: the count rows and the design rows list different donors (or orders)")
    tmm = export_tmm.copy()
    tmm.index = tmm.index.astype(str)
    if sorted(tmm.index) != sorted(donors):
        raise ValueError("stage-de-limma: the export TMM factors list other donors than the counts")
    names = [*counts.columns.astype(str), *donors]
    bad = [n for n in names if any(ch in n for ch in ',"\n\r')]
    if bad:
        raise ValueError(f"stage-de-limma: name(s) {bad[:3]} hold a comma, quote or line break; the exchange cannot carry them")
    values = counts.to_numpy()
    if not np.issubdtype(values.dtype, np.integer):
        values = np.asarray(values, dtype=float)
        if not np.all(np.isfinite(values)) or not np.array_equal(values, np.round(values)):
            raise ValueError("stage-de-limma: the counts are not integer-valued")
    if (values < 0).any():
        raise ValueError("stage-de-limma: the counts contain negative values")
    d.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(values.astype(np.int64).T, index=pd.Index(counts.columns.astype(str), name="feature"),
                         columns=donors)
    frame.to_csv(d / "counts.csv")
    pd.DataFrame({"donor_id": donors, "factor": tmm.loc[donors].to_numpy(dtype=float)}).to_csv(
        d / "export_tmm.csv", index=False, float_format="%.17g"
    )
    table = design.copy()
    table.index = table.index.astype(str)
    table.index.name = "donor_id"
    for col in _STAGE_DE_LEVELS:
        table[col] = table[col].astype(str)
    dirs = {}
    for model in models:
        check = stage_de_design_check(design, model, omic)
        matrix = FormulaicContrasts(design, stage_de_formula(model, omic)).design_matrix
        md = d / model
        md.mkdir(exist_ok=True)
        for name in ("counts.csv", "export_tmm.csv"):
            target = md / name
            if target.exists():
                target.unlink()
            try:
                os.link(d / name, target)
            except OSError:
                shutil.copyfile(d / name, target)
        pd.DataFrame(np.asarray(matrix, dtype=float), index=pd.Index(donors, name="donor_id"),
                     columns=[str(c) for c in matrix.columns]).to_csv(md / "design_matrix.csv", float_format="%.17g")
        table.to_csv(md / "design_table.csv", float_format="%.17g")
        (md / "settings.json").write_text(json.dumps(_stage_de_limma_settings(model, omic, check), indent=2) + "\n")
        dirs[model] = md
    return dirs


def stage_de_limma_check(directory: str | Path, model_dirs: dict[str, Path]) -> dict:
    """Run the R script's ``--check``: donor ids, TMM against the export, every model's design rebuilt in R.

    Returns ``check.json`` (TMM result, design checks per model, versions);
    raises with R's message when a check fails.
    """
    d = Path(directory)
    _run_stage_de_limma_r(["--check", str(d), *[str(p) for p in model_dirs.values()]], what=f"checks in {d}")
    return json.loads((d / "check.json").read_text())


def _read_r_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype={"feature": str, "donor_id": str}, keep_default_na=False, na_values=["NA"],
                       float_precision="round_trip")


def stage_de_limma_fit(
    counts: pd.DataFrame,
    design: pd.DataFrame,
    model: str,
    *,
    omic: str,
    export_tmm: pd.Series,
    exported: Iterable[bool] | None = None,
    model_dir: str | Path | None = None,
) -> dict:
    """Fit one registered stage-DE model with limma-voom (quality weights) in R and test its coefficients.

    Parameters
    ----------
    counts:
        Donors × tested features raw counts (integers), rows in ``design``
        order, columns in the export's order.
    design:
        ``stage_de_design`` output for the same donors.
    model:
        A key of ``STAGE_DE_MODELS``.
    omic:
        Selects the covariates of the formula.
    export_tmm:
        The export manifest's TMM factor per donor; the R side's edgeR TMM
        must equal it to ``STAGE_DE_LIMMA_TMM_ATOL``.
    exported:
        Per feature, whether it is one of MOTCO's exported features (the
        ``exported`` column); ``None`` leaves the column out.
    model_dir:
        A model directory already written by ``write_stage_de_limma_exchange``
        (kept, with R's outputs, for a rerun by hand). ``None`` writes the
        exchange files into a temporary directory removed afterwards.

    Returns
    -------
    dict
        ``tables`` (label → DataFrame: per tested coefficient ``feature``,
        ``AveExpr``, ``logFC``, ``SE``, ``t``, ``P.Value``, ``adj.P.Val``,
        ``B``; per joint test ``feature``, ``AveExpr``, the two
        ``logFC_<coefficient>``, ``F``, ``P.Value``, ``adj.P.Val``; each with
        ``exported``, one row per feature in ``counts`` order) and ``fit``
        (formula, design columns and resolved tests, the TMM check and
        factors, the R-side design check, voom's span, the sample weights per
        donor with their minimum, median, maximum and per-cell means, the
        empirical-Bayes hyperparameters, per-test counts, the exact-F
        diagnostic, versions, the script's sha256, runtime).

    Raises
    ------
    RuntimeError
        ``Rscript`` is not on ``PATH`` (the message names the Nix shell) or
        the R process fails (the message ends with R's stderr).
    ValueError
        The inputs are inconsistent (see ``write_stage_de_limma_exchange``)
        or R's tables do not list the tested features in order or lack an
        adjusted p-value.
    """
    import shutil
    import tempfile

    if model not in STAGE_DE_MODELS:
        raise ValueError(f"unknown stage-DE model {model!r}; expected one of {tuple(STAGE_DE_MODELS)}")
    spec = STAGE_DE_MODELS[model]
    features = counts.columns.astype(str).tolist()
    mask = None if exported is None else np.asarray(list(exported), dtype=bool)
    if mask is not None and mask.size != len(features):
        raise ValueError(f"stage-de-limma: exported has {mask.size} entries for {len(features)} features")
    check = stage_de_design_check(design, model, omic)
    t0 = time.perf_counter()
    tmp = None
    try:
        if model_dir is None:
            tmp = Path(tempfile.mkdtemp(prefix="stage_de_limma_"))
            md = write_stage_de_limma_exchange(tmp, counts, design, omic=omic, export_tmm=export_tmm, models=[model])[model]
        else:
            md = Path(model_dir)
        for p in md.glob("*__*.csv"):
            p.unlink()
        for name in ("fit.json", "sample_weights.csv", "tmm.csv"):
            if (md / name).exists():
                (md / name).unlink()
        _run_stage_de_limma_r([str(md)], what=f"{omic} {model}")
        r_fit = json.loads((md / "fit.json").read_text())
        tables = {}
        labels = [(label, f"coef__{label}.csv", False) for label, _ in spec["tested"]]
        if model in STAGE_DE_LIMMA_JOINT_TESTS:
            labels.append((STAGE_DE_LIMMA_JOINT_TESTS[model][0], f"joint__{STAGE_DE_LIMMA_JOINT_TESTS[model][0]}.csv", True))
        for label, name, joint in labels:
            table = _read_r_table(md / name)
            if table["feature"].tolist() != features:
                raise ValueError(f"stage-de-limma: R's {name} for {omic} {model} does not list the tested features in order")
            if table["adj.P.Val"].isna().any():
                raise ValueError(f"stage-de-limma: R's {name} for {omic} {model} lacks an adjusted p-value")
            if mask is not None:
                table["exported"] = mask
            tables[label] = table
        weights = _read_r_table(md / "sample_weights.csv").set_index("donor_id")["weight"]
        tmm = _read_r_table(md / "tmm.csv").set_index("donor_id")
    finally:
        if tmp is not None:
            shutil.rmtree(tmp, ignore_errors=True)
    runtime = time.perf_counter() - t0

    sex = design["sex"].astype(str)
    stage = design["stage"].astype(str)
    cells = {}
    for s in _STAGE_DE_LEVELS["sex"]:
        for t in _STAGE_DE_LEVELS["stage"]:
            members = design.index[(sex == s).to_numpy() & (stage == t).to_numpy()].astype(str)
            cells[f"{s}/stage {t}"] = {"n": int(len(members)), "mean_weight": float(weights.loc[members].mean())}
    versions = dict(r_fit["versions"])
    versions["script_sha256"] = _sha256_file(stage_de_limma_script())
    joint_info = {}
    if model in STAGE_DE_LIMMA_JOINT_TESTS:
        label, coefficients = STAGE_DE_LIMMA_JOINT_TESTS[model]
        joint_info[label] = {
            "coefficients": list(coefficients),
            "design_columns": [check["resolved"][c] for c in coefficients],
            **r_fit["joint"][label],
        }
    fit = {
        "model": model,
        "omic": omic,
        "formula": check["formula"],
        "stage_coding": spec["stage_coding"],
        "interaction": spec["interaction"],
        "design_columns": check["design_columns"],
        "tested": [
            {"coefficient": label, "expected_column": expected, "design_column": check["resolved"][label]}
            for label, expected in spec["tested"]
        ],
        "joint": joint_info,
        "n_donors": int(r_fit["n_donors"]),
        "n_features": int(r_fit["n_features"]),
        "n_exported": None if mask is None else int(mask.sum()),
        "tmm": {**r_fit["tmm"], "factors": {str(k): float(v) for k, v in tmm["factor"].items()}},
        "design_check": r_fit["design_check"],
        "installed_defaults_equal_method": r_fit["installed_defaults_equal_method"],
        "voom": r_fit["voom"],
        "sample_weights": {
            "per_donor": {str(k): float(v) for k, v in weights.items()},
            "min": float(weights.min()),
            "median": float(weights.median()),
            "max": float(weights.max()),
            "per_cell_mean": cells,
        },
        "ebayes": r_fit["ebayes"],
        "coefficients": {
            label: {
                "n_padj_missing": int(t["adj.P.Val"].isna().sum()),
                "n_significant": int((t["adj.P.Val"] <= STAGE_DE_ALPHA).sum()),
            }
            for label, t in tables.items()
        },
        "settings": STAGE_DE_LIMMA_METHOD,
        "versions": versions,
        "runtime_seconds": runtime,
        "r_runtime_seconds": r_fit["runtime_seconds"],
    }
    return {"tables": tables, "fit": fit}


# The parts of a fit record the summary and the manifest read; a record lacking
# one (an earlier version's, or a damaged one) is stale
_STAGE_DE_LIMMA_RECORD_KEYS = (
    ("tables",), ("fit", "tmm", "max_abs_diff"), ("fit", "sample_weights", "per_cell_mean"), ("fit", "ebayes"),
    ("fit", "coefficients"), ("fit", "versions", "script_sha256"), ("fit", "tested"), ("fit", "joint"),
    ("fit", "design_check"), ("fit", "voom"),
)


def _stage_de_limma_labels(model: str) -> list[tuple[str, str]]:
    """(label, kind) of every table of one model: its moderated t's, then its joint F."""
    labels = [(label, "moderated_t") for label, _ in STAGE_DE_MODELS[model]["tested"]]
    if model in STAGE_DE_LIMMA_JOINT_TESTS:
        labels.append((STAGE_DE_LIMMA_JOINT_TESTS[model][0], "joint_f"))
    return labels


def _stage_de_limma_fingerprint(prep: dict, cohort: str, omic: str, model: str, versions: dict) -> dict:
    """What a limma fit depends on: stage-de's fingerprint without pydeseq2, plus the method, joint tests, script and R versions."""
    fingerprint = _stage_de_fingerprint(prep, cohort, omic, model)
    fingerprint.pop("pydeseq2")
    fingerprint.update({
        "joint": [STAGE_DE_LIMMA_JOINT_TESTS[model][0], list(STAGE_DE_LIMMA_JOINT_TESTS[model][1])]
        if model in STAGE_DE_LIMMA_JOINT_TESTS else None,
        "method": STAGE_DE_LIMMA_METHOD,
        "tmm_atol": STAGE_DE_LIMMA_TMM_ATOL,
        "script_sha256": versions["script_sha256"],
        "r_versions": {k: versions[k] for k in ("R", "limma", "edgeR", "jsonlite")},
    })
    return fingerprint


def _stage_de_limma_model_state(omic_dir: Path, model: str, fingerprint: dict) -> str:
    """``complete`` (fit record and tables match), ``missing`` (to fit) or ``stale`` (stage-de's resume semantics)."""
    record_path = omic_dir / STAGE_DE_LIMMA_FIT_TEMPLATE.format(model=model)
    if not record_path.exists():
        return "missing"
    try:
        record = json.loads(record_path.read_text())
    except (OSError, json.JSONDecodeError, UnicodeDecodeError):
        return "stale"
    if not isinstance(record, dict) or record.get("fingerprint") != fingerprint:
        return "stale"
    for path in _STAGE_DE_LIMMA_RECORD_KEYS:
        node = record
        for key in path:
            if not isinstance(node, dict) or key not in node:
                return "stale"
            node = node[key]
    tables = record["tables"]
    if not isinstance(tables, dict) or not all(
        isinstance(info, dict) and isinstance(info.get("file"), str) and isinstance(info.get("sha256"), str)
        for info in tables.values()
    ):
        return "stale"
    for info in tables.values():
        p = omic_dir / info["file"]
        if not p.is_file() or _sha256_file(p) != info["sha256"]:
            return "missing"
    if set(tables) != {label for label, _ in _stage_de_limma_labels(model)}:
        return "missing"
    return "complete"


def _read_stage_de_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype={"feature": str}, keep_default_na=False, na_values=[""], float_precision="round_trip")


def _stage_de_limma_deseq2_reference(
    deseq2_dir: Path, cohort: str, prep: dict, *, cell_type: str = DEFAULT_CELL_TYPE
) -> dict:
    """Read the DESeq2 stage-DE outputs of one cohort and check they belong to the inputs now on disk (design D10).

    The DESeq2 manifest must record ``cell_type`` (checked first, so a
    DESeq2 directory of another cell type is named as such). Raises a
    ``FileNotFoundError`` or ``ValueError`` naming the file and the
    ``stage-de`` command; nothing is written under ``deseq2_dir``.
    """
    d = deseq2_dir / cohort
    rerun = "rerun `sead-sexdiff stage-de --force` (DESeq2) on the inputs now on disk; nothing was written"
    manifest_path = d / STAGE_DE_MANIFEST
    for p in (manifest_path, d / STAGE_DE_SUMMARY):
        if not p.exists():
            raise FileNotFoundError(f"{p} is missing; run `sead-sexdiff stage-de` (DESeq2) to produce it; nothing was written")
    files = {str(p.relative_to(d)): _sha256_file(p) for p in sorted(d.rglob("*")) if p.is_file()}
    manifest = json.loads(manifest_path.read_text())
    cell_type_findings = check_manifest_cell_type(manifest, manifest_path, cell_type)
    if not manifest.get("complete") or sorted(manifest.get("omics", [])) != sorted(_STAGE_DE_OMICS):
        raise ValueError(
            f"{manifest_path} is not complete for both omics (complete = {manifest.get('complete')}, omics "
            f"{manifest.get('omics')}); run `sead-sexdiff stage-de --omic both`; nothing was written"
        )
    recorded = manifest.get("inputs", {}).get("files", {})
    for group, names in prep["hashes"].items():
        for name, sha in names.items():
            then = recorded.get(group, {}).get(name, {}).get("sha256")
            if then != sha:
                raise ValueError(
                    f"{manifest_path} records sha256 {then} for {group}/{name} but the file now on disk has {sha}: the "
                    f"DESeq2 run was made from other inputs; {rerun}"
                )
    extra = sorted(f"{g}/{n}" for g, names in recorded.items() for n in names if n not in prep["hashes"].get(g, {}))
    if extra:
        raise ValueError(f"{manifest_path} records input(s) {extra} that this step does not read; {rerun}")
    tables = {}
    for o in _STAGE_DE_OMICS:
        tables[o] = {}
        for m, spec in STAGE_DE_MODELS.items():
            entries = manifest.get("tables", {}).get(o, {}).get(m, {})
            tables[o][m] = {}
            for label, _ in spec["tested"]:
                info = entries.get(label)
                if not isinstance(info, dict):
                    raise ValueError(f"{manifest_path} records no table for {o} {m} {label}; {rerun}")
                p = d / o / info["file"]
                if not p.is_file() or files.get(str(p.relative_to(d))) != info["sha256"]:
                    raise ValueError(
                        f"{p} does not hash as {manifest_path} records (sha256 {info['sha256']}); the DESeq2 table was "
                        f"edited or replaced; {rerun}"
                    )
                tables[o][m][label] = _read_stage_de_table(p)
    summary = pd.read_csv(d / STAGE_DE_SUMMARY, dtype={"design_column": str}, float_precision="round_trip")
    rows = []
    for o in _STAGE_DE_OMICS:
        rows += _stage_de_summary_rows(cohort, o, tables[o], prep["delta_rejected"])
    recomputed = pd.DataFrame(rows, columns=list(STAGE_DE_SUMMARY_COLUMNS))
    keys = ["omic", "model", "coefficient", "feature_set"]
    merged = recomputed.merge(summary, on=keys, how="outer", suffixes=("_now", "_file"), indicator=True)
    bad = merged[merged["_merge"] != "both"]
    if len(bad):
        raise ValueError(f"{d / STAGE_DE_SUMMARY} lists other rows than its tables give ({len(bad)} row(s) differ); {rerun}")
    for col in ("n_tested", "n_significant", "n_padj_missing"):
        diff = merged[merged[f"{col}_now"].astype(int) != merged[f"{col}_file"].astype(int)]
        if len(diff):
            r = diff.iloc[0]
            raise ValueError(
                f"{d / STAGE_DE_SUMMARY} gives {col} = {int(r[f'{col}_file'])} for {r['omic']} {r['model']} "
                f"{r['coefficient']} ({r['feature_set']}) but its tables give {int(r[f'{col}_now'])}; {rerun}"
            )
    cells = {}
    rule = STAGE_DE_CONCORDANCE_RULE
    for o in _STAGE_DE_OMICS:
        sel = (merged["omic"] == o) & (merged["model"] == rule["model"]) & (merged["coefficient"] == rule["coefficient"]) & (
            merged["feature_set"] == rule["feature_set"])
        r = merged[sel].iloc[0]
        now = {"interaction_detected": bool(r["interaction_detected_now"]), "delta_rejected": bool(r["delta_rejected_now"]),
               "cell": str(r["concordance_now"])}
        file = {"interaction_detected": str(r["interaction_detected_file"]) == "True",
                "delta_rejected": str(r["delta_rejected_file"]) == "True", "cell": str(r["concordance_file"])}
        recorded_cell = manifest.get("concordance", {}).get(o, {})
        if now != file or recorded_cell.get("cell") != now["cell"] or recorded_cell.get(
                "n_interaction_exported") != int(r["n_significant_now"]):
            raise ValueError(
                f"{d / STAGE_DE_SUMMARY} or {manifest_path} gives the {o} concordance cell {file} but the tables give "
                f"{now}; {rerun}"
            )
        cells[o] = {**now, "n_interaction_exported": int(r["n_significant_now"])}
    return {
        "dir": d, "manifest": manifest, "files": files, "tables": tables, "summary": recomputed, "concordance": cells,
        "manifest_sha256": files[STAGE_DE_MANIFEST], "summary_sha256": files[STAGE_DE_SUMMARY],
        "software": manifest.get("software", {}), "cell_type_findings": cell_type_findings,
    }


def _stage_de_limma_summary_rows(
    cohort: str, omic: str, tables: dict[str, dict[str, pd.DataFrame]], delta_rejected: bool,
) -> list[dict]:
    """Summary rows of one omic: per model × moderated t, joint F and union × feature set."""
    rule = STAGE_DE_LIMMA_CONCORDANCE_RULE
    rows = []
    for model, spec in STAGE_DE_MODELS.items():
        model_tables = tables[model]
        labels = _stage_de_limma_labels(model)
        features = model_tables[labels[0][0]]["feature"].tolist()
        for label, _ in labels[1:]:
            if model_tables[label]["feature"].tolist() != features:
                raise ValueError(f"stage-de-limma: the {model} tables of {cohort} {omic} list different features")
        exported = model_tables[labels[0][0]]["exported"].to_numpy(dtype=bool)
        columns = dict(spec["tested"])
        entries = []
        for label, kind in labels:
            column = columns[label] if kind == "moderated_t" else ";".join(
                columns[c] for c in STAGE_DE_LIMMA_JOINT_TESTS[model][1])
            entries.append((label, kind, column, model_tables[label]["adj.P.Val"].to_numpy(dtype=float)[:, None]))
        if len(spec["tested"]) > 1:
            entries.append((
                STAGE_DE_UNION, "union", ";".join(col for _, col in spec["tested"]),
                np.column_stack([model_tables[label]["adj.P.Val"].to_numpy(dtype=float) for label, _ in spec["tested"]]),
            ))
        for label, kind, column, padj in entries:
            for feature_set in STAGE_DE_FEATURE_SETS:
                keep = exported if feature_set == "exported" else np.ones(len(features), dtype=bool)
                p = padj[keep]
                finite = np.isfinite(p)
                significant = (np.where(finite, p, np.inf) <= STAGE_DE_ALPHA).any(axis=1)
                row = {
                    "cohort": cohort, "omic": omic, "model": model, "coefficient": label, "kind": kind,
                    "design_column": column, "feature_set": feature_set, "n_tested": int(keep.sum()),
                    "n_padj_missing": int((~finite).all(axis=1).sum()), "n_significant": int(significant.sum()),
                    "min_padj": float(np.min(p[finite])) if finite.any() else float("nan"),
                    "interaction_detected": None, "delta_rejected": None, "concordance": None,
                }
                if (model, label, feature_set) == (rule["model"], rule["test"], rule["feature_set"]):
                    detected = row["n_significant"] > 0
                    row["interaction_detected"] = detected
                    row["delta_rejected"] = bool(delta_rejected)
                    row["concordance"] = "concordant" if detected == bool(delta_rejected) else "discordant"
                rows.append(row)
    return rows


def _stage_de_limma_comparison_rows(
    cohort: str, omic: str, limma: dict[str, dict[str, pd.DataFrame]], deseq2: dict[str, dict[str, pd.DataFrame]],
    concordance: dict,
) -> list[dict]:
    """Method-comparison rows of one omic (design D9): coefficient, union and joint-F-against-union rows per feature set."""
    from scipy.stats import spearmanr

    def significant(frame: pd.DataFrame, col: str) -> np.ndarray:
        p = frame[col].to_numpy(dtype=float)
        return np.where(np.isfinite(p), p, np.inf) <= STAGE_DE_ALPHA

    rows = []
    for model, spec in STAGE_DE_MODELS.items():
        labels = [label for label, _ in spec["tested"]]
        features = deseq2[model][labels[0]]["feature"].tolist()
        for label, _ in _stage_de_limma_labels(model):
            if limma[model][label]["feature"].tolist() != features:
                raise ValueError(
                    f"stage-de-limma: the {cohort} {omic} {model} {label} table lists other features than the DESeq2 "
                    "tables; both must come from the same export"
                )
        exported = deseq2[model][labels[0]]["exported"].to_numpy(dtype=bool)
        if not np.array_equal(exported, limma[model][labels[0]]["exported"].to_numpy(dtype=bool)):
            raise ValueError(f"stage-de-limma: the {cohort} {omic} {model} tables disagree on `exported`")
        entries = [
            (label, "coefficient", label, label, significant(deseq2[model][label], "padj"),
             significant(limma[model][label], "adj.P.Val"), label)
            for label in labels
        ]
        if len(labels) > 1:
            d_union = np.any([significant(deseq2[model][label], "padj") for label in labels], axis=0)
            l_union = np.any([significant(limma[model][label], "adj.P.Val") for label in labels], axis=0)
            entries.append((STAGE_DE_UNION, "union", STAGE_DE_UNION, STAGE_DE_UNION, d_union, l_union, None))
            joint = STAGE_DE_LIMMA_JOINT_TESTS[model][0]
            entries.append((joint, "joint_f_vs_union", STAGE_DE_UNION, joint, d_union,
                            significant(limma[model][joint], "adj.P.Val"), None))
        for test, kind, d_test, l_test, d_sig, l_sig, coef in entries:
            for feature_set in STAGE_DE_FEATURE_SETS:
                keep = exported if feature_set == "exported" else np.ones(len(features), dtype=bool)
                ds, ls = d_sig[keep], l_sig[keep]
                row = {
                    "cohort": cohort, "omic": omic, "model": model, "test": test, "kind": kind, "deseq2_test": d_test,
                    "limma_test": l_test, "feature_set": feature_set, "n_tested": int(keep.sum()),
                    "n_significant_deseq2": int(ds.sum()), "n_significant_limma": int(ls.sum()),
                    "n_both": int((ds & ls).sum()), "n_only_deseq2": int((ds & ~ls).sum()),
                    "n_only_limma": int((~ds & ls).sum()), "detected_deseq2": bool(ds.any()),
                    "detected_limma": bool(ls.any()), "detection_agrees": bool(ds.any()) == bool(ls.any()),
                    "spearman_log2fc": None, "spearman_stat": None, "n_compared": None,
                    "delta_rejected": None, "concordance_deseq2": None, "concordance_limma": None, "methods_agree": None,
                }
                if coef is not None:
                    dt, lt = deseq2[model][coef], limma[model][coef]
                    pairs = {
                        "spearman_log2fc": (dt["log2FoldChange"], lt["logFC"]),
                        "spearman_stat": (dt["stat"], lt["t"]),
                    }
                    finite = keep.copy()
                    for a, b in pairs.values():
                        finite &= np.isfinite(a.to_numpy(dtype=float)) & np.isfinite(b.to_numpy(dtype=float))
                    row["n_compared"] = int(finite.sum())
                    for name, (a, b) in pairs.items():
                        row[name] = float(spearmanr(a.to_numpy(dtype=float)[finite], b.to_numpy(dtype=float)[finite])[0]) \
                            if finite.sum() > 2 else float("nan")
                rule = STAGE_DE_LIMMA_CONCORDANCE_RULE
                if (model, test, feature_set) == (rule["model"], rule["test"], rule["feature_set"]):
                    row.update({
                        "delta_rejected": concordance["delta_rejected"],
                        "concordance_deseq2": concordance["deseq2"]["cell"],
                        "concordance_limma": concordance["limma"]["cell"],
                        "methods_agree": concordance["methods_agree"],
                    })
                rows.append(row)
    return rows


def render_stage_de_limma(summary: pd.DataFrame, comparison: pd.DataFrame, manifest: dict) -> str:
    """Render ``stage_de_limma.md`` from the summary, the method comparison and the manifest.

    Parameters
    ----------
    summary:
        The cohort's ``stage_de_limma_summary.csv``.
    comparison:
        The cohort's ``method_comparison.csv``.
    manifest:
        The cohort's ``stage_de_limma_manifest.json``.

    Returns
    -------
    str
        Markdown: the reading rules, the side-by-side concordance table, one
        count table per omic with the DESeq2 counts beside limma's, the
        comparison rows, the sample weights per model and cell, the TMM
        check, the exact-F diagnostic, the inputs and the findings. Every
        number comes from the three arguments.
    """
    c = manifest["cohort"]
    alpha = manifest["alpha"]
    ref = manifest["reference"]
    sw = manifest["software"]
    lines = [
        f"# Pseudobulk stage DE, limma-voom: {c}",
        "",
        f"A second predeclared method beside DESeq2 (`stage-de`), on the same {manifest['n_donors']} {c} donors, raw "
        "pseudobulk counts, registered models, design matrices and abundance-passing features. Method: edgeR TMM on "
        "the tested counts with the export's trims (checked against the export manifest's factors) → "
        "voomWithQualityWeights (voom's precision weights times one quality weight per donor) → lmFit → eBayes, with "
        f"limma's defaults (R {sw['R']}, limma {sw['limma']}, edgeR {sw['edgeR']}; script sha256 "
        f"`{sw['stage_de_limma_script_sha256'][:12]}…`). Moderated t per tested coefficient; a joint moderated F per "
        "factor model over its two coefficients (`stage_joint`, `interaction_joint`); a feature is significant when "
        f"its Benjamini–Hochberg adjusted p ≤ {alpha} over all tested features (limma has no independent filter; "
        "DESeq2's BH family is smaller where its filter removes features). voom's log-CPM "
        "(`log2((count + 0.5) / (lib.size × factor + 1) × 1e6)`) is not the export's prior-2 log2-CPM. The joint F "
        "combines the moderated t's through the unweighted design's coefficient correlation, so with voom's weights "
        "its claim is approximate (the exact-F diagnostic below measures it); a factor model's `union` counts features "
        "significant on either moderated t and carries no joint error-rate claim.",
        "",
        "Reading rules: this change was proposed after the DESeq2 results were known, and the DESeq2 results stay the "
        "reported ones; a call that differs between the methods is method dependence, not a result of either; the "
        "same donors and counts enter both methods, so agreement is not replication; a zero count means none detected "
        "at this n, not no effect; the concordance cells cross-check the primary analysis and are not a second test "
        "of its hypothesis, whose non-rejection rules are untouched (`angle` and `shape` non-rejections are "
        "uninformative). The step read sex and stage, so it is a check beside the analysis, never a filter: no "
        "export, analysis or sensitivity run reads its outputs.",
        "",
        "## Concordance",
        "",
        f"Rule (`STAGE_DE_LIMMA_CONCORDANCE_RULE`, predeclared): `interaction_detected` = "
        f"{manifest['concordance_rule']['interaction_detected']}; `delta_rejected` = "
        f"{manifest['concordance_rule']['delta_rejected']}; {manifest['concordance_rule']['cell']}. "
        f"{manifest['concordance_rule']['beside'][:1].upper()}{manifest['concordance_rule']['beside'][1:]}. The primary "
        f"analysis (`{ref['dir']}`, CV-selected rank "
        f"{ref['selected_lv']}, MOTCO {ref['motco']}) is read, never rerun or modified.",
        "",
        "| omic | DESeq2: exported features (union of two Wald tests) | DESeq2 cell | limma: exported features "
        "(interaction_joint) | limma cell | delta rejected | methods agree |",
        "|---|---:|---|---:|---|:---:|:---:|",
    ]
    for o, cell in manifest["concordance"].items():
        lines.append(
            f"| {o} | {cell['deseq2']['n_interaction_exported']} | {cell['deseq2']['cell']} | "
            f"{cell['limma']['n_interaction_exported']} | {cell['limma']['cell']} | "
            f"{'yes' if cell['delta_rejected'] else 'no'} | {'yes' if cell['methods_agree'] else 'no'} |"
        )
    kind_label = {"moderated_t": "moderated t", "joint_f": "joint F", "union": "union (no joint error rate)"}
    for o in manifest["omics"]:
        sub = summary[summary["omic"] == o]
        comp = comparison[comparison["omic"] == o]
        lines += [
            "",
            f"## Counts: {o}",
            "",
            f"{manifest['feature_sets'][o]['exported']} exported and {manifest['feature_sets'][o]['abundance']} "
            f"abundance-passing features tested; significant / tested at adjusted p ≤ α. The DESeq2 columns are the "
            "same coefficient's Wald test, the union for a union row, and DESeq2's union of the model's two Wald tests "
            "for a joint-F row (a different test, `joint_f_vs_union`).",
            "",
            "| model | test | kind | design column | limma exported | limma abundance | DESeq2 exported | DESeq2 abundance "
            "| min adj. p limma (exported / abundance) |",
            "|---|---|---|---|---:|---:|---:|---:|---:|",
        ]
        for (model, test), pair in sub.groupby(["model", "coefficient"], sort=False):
            pair = pair.set_index("feature_set")
            e, a = pair.loc["exported"], pair.loc["abundance"]
            cmp = comp[(comp["model"] == model) & (comp["test"] == test)].set_index("feature_set")
            de, da = cmp.loc["exported"], cmp.loc["abundance"]
            lines.append(
                f"| {model} | {test} | {kind_label[e['kind']]} | `{e['design_column']}` | "
                f"{int(e['n_significant'])} / {int(e['n_tested'])} | {int(a['n_significant'])} / {int(a['n_tested'])} | "
                f"{int(de['n_significant_deseq2'])} / {int(de['n_tested'])} | "
                f"{int(da['n_significant_deseq2'])} / {int(da['n_tested'])} | "
                f"{_fmt_nan(float(e['min_padj']), '.3g')} / {_fmt_nan(float(a['min_padj']), '.3g')} |"
            )
        lines += [
            "",
            f"## Method comparison: {o}",
            "",
            "Features called by both methods and by one only, whether each method detects anything, and, per "
            "coefficient, the Spearman correlations of DESeq2's log2 fold change with limma's logFC and of the Wald "
            "statistic with the moderated t (over the features with both values finite). Descriptive: no p-value, "
            "no gate; a row whose detection differs is a method-dependence finding.",
            "",
            "| model | test | kind | set | DESeq2 | limma | both | only DESeq2 | only limma | detection agrees | "
            "ρ log2FC | ρ statistic | compared |",
            "|---|---|---|---|---:|---:|---:|---:|---:|:---:|---:|---:|---:|",
        ]
        for _, r in comp.iterrows():
            rho = [_fmt_nan(None if pd.isna(r[k]) else float(r[k]), ".3f") for k in ("spearman_log2fc", "spearman_stat")]
            n_comp = "" if pd.isna(r["n_compared"]) else str(int(r["n_compared"]))
            lines.append(
                f"| {r['model']} | {r['test']} | {r['kind']} | {r['feature_set']} | {int(r['n_significant_deseq2'])} | "
                f"{int(r['n_significant_limma'])} | {int(r['n_both'])} | {int(r['n_only_deseq2'])} | "
                f"{int(r['n_only_limma'])} | {'yes' if bool(r['detection_agrees']) else 'no'} | {rho[0]} | {rho[1]} | "
                f"{n_comp} |"
            )
    lines += [
        "",
        "## Sample weights",
        "",
        "One quality weight per donor and model (arrayWeights inside voomWithQualityWeights, shrunk toward 1 with "
        "`prior.n = 10`): minimum, median, maximum and the mean per sex × stage cell. Shown without a gate; a cell "
        "with a low mean weight counts less in that model's fit.",
        "",
    ]
    cell_names = None
    for o, fits in manifest["fits"].items():
        for model, fit in fits.items():
            w = fit["sample_weights"]
            if cell_names is None:
                cell_names = list(w["per_cell_mean"])
                lines += [
                    "| omic | model | min | median | max | " + " | ".join(cell_names) + " |",
                    "|---|---|---:|---:|---:|" + "---:|" * len(cell_names),
                ]
            lines.append(
                f"| {o} | {model} | {w['min']:.3f} | {w['median']:.3f} | {w['max']:.3f} | "
                + " | ".join(f"{w['per_cell_mean'][k]['mean_weight']:.3f}" for k in cell_names) + " |"
            )
    lines += [
        "",
        "## TMM check",
        "",
        f"edgeR's TMM factors on the tested counts (library size = their column sums, the export's trims) against "
        f"the export manifest's, absolute tolerance {manifest['tmm_atol']}:",
        "",
    ]
    for o, chk in manifest["tmm_check"].items():
        lines.append(
            f"- {o}: largest absolute difference {chk['max_abs_diff']:.3g} over {chk['n_donors']} donors and "
            f"{chk['n_features']} features ({'passed' if chk['passed'] else 'FAILED'})"
        )
    lines += [
        "",
        "## Joint F diagnostic",
        "",
        "limma's F is recomputed from the moderated t's and the unweighted design's coefficient correlation (it must "
        "equal topTable's), and compared with the exact F from each feature's own weighted coefficient covariance "
        "(same variance and degrees of freedom; BH over all tested features). limma's F stays the reported test.",
        "",
        "| omic | model | test | df1 / df2 | limma F recomputed (max rel. diff) | exact F: max / median rel. diff | "
        "significant: limma / exact | changed calls |",
        "|---|---|---|---|---:|---:|---:|---|",
    ]
    for o, fits in manifest["fits"].items():
        for model, fit in fits.items():
            for label, j in fit["joint"].items():
                ex = j["exact"]
                changed = ", ".join(f"`{x['feature']}`" for x in ex["changed_calls"]) or "none"
                lines.append(
                    f"| {o} | {model} | {label} | {j['df1']} / {_fmt(j['df2'], '.4g')} | "
                    f"{j['limma_f_recomputed_max_rel_diff']:.2g} | {ex['max_rel_diff']:.3g} / {ex['median_rel_diff']:.3g} | "
                    f"{j['n_significant']} / {ex['n_significant_exact']} | {ex['n_changed_calls']}: {changed} |"
                )
    lines += [
        "",
        "## Fits",
        "",
        "| omic | model | features × donors | voom span | df prior | s² prior | runtime (s) |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for o, fits in manifest["fits"].items():
        for model, fit in fits.items():
            lines.append(
                f"| {o} | {model} | {fit['n_features']} × {fit['n_donors']} | {fit['voom']['span']:.4f} | "
                f"{_fmt(fit['ebayes']['df_prior'], '.4g')} | {_fmt(fit['ebayes']['s2_prior'], '.4g')} | "
                f"{_fmt(fit['runtime_seconds'], '.1f')} |"
            )
    lines += ["", "## Inputs", ""]
    for group, files in manifest["inputs"]["files"].items():
        for name, info in files.items():
            lines.append(f"- {group}: `{info['path']}` sha256 `{info['sha256']}`")
    dq = manifest["inputs"]["deseq2"]
    lines.append(
        f"- DESeq2 (`stage-de`, read and checked, never written): `{dq['dir']}`, manifest sha256 "
        f"`{dq['manifest_sha256']}`, summary sha256 `{dq['summary_sha256']}`, {len(dq['files'])} files hashed before "
        "and after the run"
    )
    lines += ["", "## Findings", ""]
    lines += [f"- {msg}" for msg in manifest["findings"]] or ["None."]
    return "\n".join(lines) + "\n"


def run_stage_de_limma(
    cohort: str = "all",
    *,
    omic: str = "both",
    cell_type: str = DEFAULT_CELL_TYPE,
    pseudobulk_dir: str | Path | None = None,
    input_dir: str | Path | None = None,
    analysis_dir: str | Path | None = None,
    deseq2_dir: str | Path | None = None,
    out_dir: str | Path | None = None,
    force: bool = False,
    debug_dir: str | Path | None = None,
) -> dict:
    """limma-voom fits of the registered stage-DE models beside the DESeq2 step, with the method comparison.

    Parameters
    ----------
    cohort:
        ``"cohort_30"``, ``"cohort_50"`` or ``"all"``.
    omic:
        ``"rna"``, ``"atac"`` or ``"both"``: the omics fitted in this run.
        The summary covers every omic whose fits are complete under the
        cohort directory.
    pseudobulk_dir, input_dir, analysis_dir:
        As for ``run_stage_de`` (read and checked through the same code
        path, ``_stage_de_prepare``; never written).
    deseq2_dir:
        The DESeq2 step's output root; ``<deseq2_dir>/<cohort>/`` must be a
        complete ``stage-de`` run made from the inputs now on disk (read and
        checked, never written).
    out_dir:
        Output root: ``<out_dir>/<cohort>/`` with ``<omic>/<model>__<label>.csv``
        (per tested coefficient and joint test), ``<omic>/<model>__fit.json``
        (the resume record), ``stage_de_limma_summary.csv``,
        ``method_comparison.csv``, ``stage_de_limma.md`` and
        ``stage_de_limma_manifest.json``. It may not lie inside an input
        directory or the DESeq2 directory.
    force:
        Refit the requested omics of a cohort even when their fits are
        complete; without it a complete fit is skipped and a fit made from
        other inputs, settings, script or R versions is an error.
    debug_dir:
        Python API only: keep the exchange files and R's outputs under
        ``<debug_dir>/<cohort>_<omic>/<model>/`` so the R script can be rerun
        by hand (``Rscript --vanilla stage_de_limma.R <model_dir>``). Every
        requested model must then be fitted in the call.

    Returns
    -------
    dict
        Cohort → ``{"dir", "summary", "comparison", "manifest", "concordance", "findings", "fitted", "skipped"}``.

    Raises
    ------
    FileNotFoundError
        An input file or a DESeq2 output is missing (the message names it and
        the command that produces it).
    RuntimeError
        ``Rscript`` is missing (the message names the Nix shell) or R fails
        (the message ends with R's error output).
    ValueError
        A reference check of the DESeq2 step fails, the DESeq2 outputs do not
        belong to the inputs now on disk, the TMM factors or a design differ
        in R, a fit record is stale (without ``force``), or ``out_dir`` lies
        inside an input or the DESeq2 directory.

    Notes
    -----
    Every check of every cohort runs before anything is written: the DESeq2
    step's preparation and reproduction, the DESeq2 outputs against the
    inputs, and R's donor-id, TMM and design checks. Fits run one model at a
    time, tables first and then the record (``.tmp`` files renamed). After a
    cohort the inputs and the DESeq2 outputs are hashed again and must be
    unchanged. The step reads sex and stage, so its outputs are an input of no
    export, analysis or sensitivity step.
    """
    import shutil
    import tempfile
    from importlib.metadata import PackageNotFoundError, version

    log = "stage-de-limma"
    cohorts = list(EXPORT_COHORTS) if cohort == "all" else [cohort]
    for c in cohorts:
        if c not in EXPORT_COHORTS:
            raise ValueError(f"Unknown cohort {c!r}; expected one of {EXPORT_COHORTS} or 'all'")
    omics = list(_STAGE_DE_OMICS) if omic == "both" else [omic]
    for o in omics:
        if o not in _STAGE_DE_OMICS:
            raise ValueError(f"Unknown omic {o!r}; expected one of {_STAGE_DE_OMICS} or 'both'")
    cell_type_entry(cell_type)
    pb = _resolve_dir(pseudobulk_dir, cell_type, "pseudobulk")
    input_dir = _resolve_dir(input_dir, cell_type, "motco_input")
    analysis_dir = _resolve_dir(analysis_dir, cell_type, "motco_analysis")
    deseq2_dir = _resolve_dir(deseq2_dir, cell_type, "stage_de")
    out_dir = _resolve_dir(out_dir, cell_type, "stage_de_limma")
    debug = None if debug_dir is None else Path(debug_dir)
    for target_name, target in (("--out-dir", out_dir), ("debug_dir", debug)):
        if target is None:
            continue
        for name, d in (("--pseudobulk-dir", pb), ("--input-dir", input_dir), ("--analysis-dir", analysis_dir),
                        ("--deseq2-dir", deseq2_dir)):
            if _same_dir(target, d):
                raise ValueError(
                    f"{target_name} {target} lies inside the input directory {d} ({name}); the step never writes "
                    "under its inputs; nothing was written"
                )
    versions = stage_de_limma_versions()
    motco_version = check_motco_version()
    print(f"[{log}] Cell type: {cell_type}")
    print(f"[{log}] Pseudobulk directory: {pb}")
    print(f"[{log}] Primary export directory: {input_dir}")
    print(f"[{log}] Primary analysis directory: {analysis_dir} (read and checked; it stays primary)")
    print(f"[{log}] DESeq2 stage-DE directory: {deseq2_dir} (read and checked; its results stay the reported ones)")
    print(f"[{log}] Output directory: {out_dir}")
    print(
        f"[{log}] R {versions['R']}, limma {versions['limma']}, edgeR {versions['edgeR']}; script "
        f"{stage_de_limma_script().name} sha256 {versions['script_sha256'][:12]}…; models: {', '.join(STAGE_DE_MODELS)}; "
        f"joint tests: {', '.join(j for j, _ in STAGE_DE_LIMMA_JOINT_TESTS.values())}; alpha {STAGE_DE_ALPHA}; omics "
        f"fitted: {', '.join(omics)}"
    )

    work_root = debug if debug is not None else Path(tempfile.mkdtemp(prefix="stage_de_limma_"))
    try:
        # Every check of every cohort before anything is written
        prepared = {}
        for c in cohorts:
            prep = _stage_de_prepare(
                c, pseudobulk_dir=pb, input_dir=input_dir, analysis_dir=analysis_dir, omics=omics,
                motco_version=motco_version, log=log, cell_type=cell_type,
            )
            prep["deseq2"] = _stage_de_limma_deseq2_reference(deseq2_dir, c, prep, cell_type=cell_type)
            print(
                f"[{log}] {c}: DESeq2 outputs {prep['deseq2']['dir']} match the inputs (hashes, tables, summary, "
                "concordance cells)"
            )
            prep["fingerprints"] = {
                o: {m: _stage_de_limma_fingerprint(prep, c, o, m, versions) for m in STAGE_DE_MODELS}
                for o in _STAGE_DE_OMICS
            }
            states = {}
            for o in _STAGE_DE_OMICS:
                for m in STAGE_DE_MODELS:
                    states[(o, m)] = _stage_de_limma_model_state(out_dir / c / o, m, prep["fingerprints"][o][m])
                    if states[(o, m)] == "stale" and not (force and o in omics):
                        hint = "rerun with --force" + ("" if o in omics else f" --omic {o} (or both)")
                        raise ValueError(
                            f"{out_dir / c / o / STAGE_DE_LIMMA_FIT_TEMPLATE.format(model=m)} is unreadable or was "
                            "made from other inputs, settings, script or R versions than those now present; "
                            f"{hint}; nothing was written"
                        )
            prep["states"] = states
            if debug is not None and not force:
                done = [f"{o}/{m}" for o in omics for m in STAGE_DE_MODELS if states[(o, m)] == "complete"]
                if done:
                    raise ValueError(
                        f"debug_dir needs every requested model fitted in this call, but {', '.join(done)} of {c} are "
                        f"complete under {out_dir / c}; pass force=True or another out_dir; nothing was written"
                    )
            # R's checks (donor ids, TMM against the export, every design rebuilt) on the exchange files
            prep["exchange"], prep["r_checks"] = {}, {}
            for o in omics:
                tmm = pd.Series({d: v["factor"] for d, v in prep["loaded"]["export_manifest"]["tmm"][o].items()})
                omic_dir = work_root / f"{c}_{o}"
                if omic_dir.exists():
                    shutil.rmtree(omic_dir)
                prep["exchange"][o] = write_stage_de_limma_exchange(
                    omic_dir, prep["counts"][o], prep["designs"][o], omic=o, export_tmm=tmm.loc[prep["donors"]],
                )
                prep["r_checks"][o] = stage_de_limma_check(omic_dir, prep["exchange"][o])
                print(
                    f"[{log}] {c} {o}: edgeR TMM equals the export's to {prep['r_checks'][o]['tmm']['max_abs_diff']:.2g} "
                    f"(tolerance {STAGE_DE_LIMMA_TMM_ATOL}); every design rebuilt by model.matrix equals the formulaic one"
                )
            print(
                f"[{log}] {c}: primary analysis {prep['reference']['dir']} matches the export (hashes, MOTCO "
                f"{motco_version}, seed, permutations, CV selection); delta rejected: "
                f"{'yes' if prep['delta_rejected'] else 'no'}"
            )
            prepared[c] = prep

        software = _software_versions()
        for pkg in ("pydeseq2", "formulaic", "formulaic-contrasts", "scipy"):
            try:
                software[pkg] = version(pkg)
            except PackageNotFoundError:
                software[pkg] = None
        software.update({"motco": motco_version, "R": versions["R"], "limma": versions["limma"],
                         "edgeR": versions["edgeR"], "jsonlite": versions["jsonlite"],
                         "stage_de_limma_script_sha256": versions["script_sha256"]})

        results = {}
        for c in cohorts:
            prep = prepared[c]
            cohort_dir = out_dir / c
            if force:
                for o in omics:
                    if (cohort_dir / o).exists():
                        shutil.rmtree(cohort_dir / o)
                    for m in STAGE_DE_MODELS:
                        prep["states"][(o, m)] = "missing"
                for name in (STAGE_DE_LIMMA_SUMMARY, STAGE_DE_LIMMA_COMPARISON, STAGE_DE_LIMMA_MARKDOWN,
                             STAGE_DE_LIMMA_MANIFEST):
                    if (cohort_dir / name).exists():
                        (cohort_dir / name).unlink()
            fitted, skipped = [], []
            for o in omics:
                omic_dir = cohort_dir / o
                omic_dir.mkdir(parents=True, exist_ok=True)
                tmm = pd.Series({d: v["factor"] for d, v in prep["loaded"]["export_manifest"]["tmm"][o].items()})
                for m in STAGE_DE_MODELS:
                    if prep["states"][(o, m)] == "complete":
                        print(f"[{log}] {c} {o} {m}: complete fit record found; skipped (rerun with --force to refit)")
                        skipped.append((o, m))
                        continue
                    out = stage_de_limma_fit(
                        prep["counts"][o], prep["designs"][o], m, omic=o, export_tmm=tmm.loc[prep["donors"]],
                        exported=prep["exported"][o], model_dir=prep["exchange"][o][m],
                    )
                    table_files = {}
                    for label, table in out["tables"].items():
                        name = STAGE_DE_LIMMA_TABLE_TEMPLATE.format(model=m, label=label)
                        tmp = omic_dir / f"{name}.tmp"
                        table.to_csv(tmp, index=False, float_format=_csv_float)
                        os.replace(tmp, omic_dir / name)
                        table_files[label] = {"file": name, "sha256": _sha256_file(omic_dir / name)}
                    record = {
                        "fingerprint": prep["fingerprints"][o][m],
                        "fit": out["fit"],
                        "export_reproduction": prep["reproduction"][o],
                        "tables": table_files,
                        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                    }
                    record_path = omic_dir / STAGE_DE_LIMMA_FIT_TEMPLATE.format(model=m)
                    tmp = record_path.with_name(record_path.name + ".tmp")
                    tmp.write_text(json.dumps(record, indent=2) + "\n")
                    os.replace(tmp, record_path)
                    prep["states"][(o, m)] = "complete"
                    fitted.append((o, m))
                    print(
                        f"[{log}] {c} {o} {m}: {out['fit']['n_features']} features × {out['fit']['n_donors']} donors in "
                        f"{out['fit']['runtime_seconds']:.1f} s; adj. p <= {STAGE_DE_ALPHA}: "
                        + ", ".join(f"{label} {info['n_significant']}" for label, info in out["fit"]["coefficients"].items())
                    )
                prep["counts"].pop(o, None)
                gc.collect()

            # Omics not fitted in this run: their state again now
            for o in _STAGE_DE_OMICS:
                if o in omics:
                    continue
                for m in STAGE_DE_MODELS:
                    state = _stage_de_limma_model_state(cohort_dir / o, m, prep["fingerprints"][o][m])
                    if state == "stale":
                        raise ValueError(
                            f"{cohort_dir / o / STAGE_DE_LIMMA_FIT_TEMPLATE.format(model=m)} became unreadable or was "
                            f"replaced by a fit from other inputs during this run; rerun with --force --omic {o}"
                        )
                    prep["states"][(o, m)] = state
            included = [o for o in _STAGE_DE_OMICS if all(prep["states"][(o, m)] == "complete" for m in STAGE_DE_MODELS)]
            records, tables = {}, {}
            for o in included:
                records[o], tables[o] = {}, {}
                for m in STAGE_DE_MODELS:
                    records[o][m] = json.loads((cohort_dir / o / STAGE_DE_LIMMA_FIT_TEMPLATE.format(model=m)).read_text())
                    tables[o][m] = {
                        label: _read_stage_de_table(cohort_dir / o / STAGE_DE_LIMMA_TABLE_TEMPLATE.format(model=m, label=label))
                        for label, _ in _stage_de_limma_labels(m)
                    }
            rows = []
            for o in included:
                rows += _stage_de_limma_summary_rows(c, o, tables[o], prep["delta_rejected"])
            summary = pd.DataFrame(rows, columns=list(STAGE_DE_LIMMA_SUMMARY_COLUMNS))
            rule = STAGE_DE_LIMMA_CONCORDANCE_RULE
            concordance, findings = {}, list(prep["cell_type_findings"]) + list(prep["deseq2"]["cell_type_findings"])
            for o in included:
                cell = summary[
                    (summary["omic"] == o) & (summary["model"] == rule["model"])
                    & (summary["coefficient"] == rule["test"]) & (summary["feature_set"] == rule["feature_set"])
                ].iloc[0]
                limma_cell = {
                    "interaction_detected": bool(cell["interaction_detected"]),
                    "n_interaction_exported": int(cell["n_significant"]),
                    "cell": str(cell["concordance"]),
                }
                deseq2_cell = prep["deseq2"]["concordance"][o]
                concordance[o] = {
                    "delta_rejected": bool(cell["delta_rejected"]),
                    "deseq2": {k: deseq2_cell[k] for k in ("interaction_detected", "n_interaction_exported", "cell")},
                    "limma": limma_cell,
                    "methods_agree": deseq2_cell["interaction_detected"] == limma_cell["interaction_detected"],
                }
                if limma_cell["cell"] == "discordant":
                    if limma_cell["interaction_detected"]:
                        what = (
                            f"{limma_cell['n_interaction_exported']} exported {o} feature(s) have adj. p <= "
                            f"{STAGE_DE_ALPHA} on {rule['test']} and the primary analysis did not reject delta"
                        )
                    else:
                        what = f"no exported {o} feature has adj. p <= {STAGE_DE_ALPHA} on {rule['test']} and the primary analysis rejected delta"
                    findings.append(
                        f"{c} {o}: limma cell discordant: {what}; reported beside both (per-feature FDR and one "
                        "latent-space permutation test ask related, not identical, questions); the primary analysis "
                        "stays primary and is unchanged"
                    )
                if not concordance[o]["methods_agree"]:
                    findings.append(
                        f"{c} {o}: method dependence of the concordance cells: limma's {rule['test']} finds "
                        f"{limma_cell['n_interaction_exported']} exported feature(s), DESeq2's union of two Wald tests "
                        f"{deseq2_cell['n_interaction_exported']}; neither cell replaces the other, the DESeq2 cell stays "
                        "the reported one, and the DESeq2 outputs and the primary analysis are unchanged"
                    )
            comp_rows = []
            for o in included:
                comp_rows += _stage_de_limma_comparison_rows(
                    c, o, tables[o], prep["deseq2"]["tables"][o], concordance[o],
                )
            comparison = pd.DataFrame(comp_rows, columns=list(STAGE_DE_LIMMA_COMPARISON_COLUMNS))
            for _, r in comparison[~comparison["detection_agrees"].astype(bool)].iterrows():
                findings.append(
                    f"{c} {r['omic']} {r['model']} {r['test']} ({r['kind']}, {r['feature_set']}): method dependence: "
                    f"DESeq2 {r['deseq2_test']} detects {int(r['n_significant_deseq2'])}, limma {r['limma_test']} detects "
                    f"{int(r['n_significant_limma'])}; a result of neither method"
                )
            for o in included:
                for m in STAGE_DE_MODELS:
                    for label, j in records[o][m]["fit"]["joint"].items():
                        ex = j["exact"]
                        if ex["n_changed_calls"]:
                            features = ", ".join(x["feature"] for x in ex["changed_calls"])
                            findings.append(
                                f"{c} {o} {m} {label}: the exact F changes {ex['n_changed_calls']} call(s) at adj. p <= "
                                f"{STAGE_DE_ALPHA} ({features}); limma's F stays the reported test"
                            )
            missing = [o for o in _STAGE_DE_OMICS if o not in included]
            for o in missing:
                findings.append(f"{c} {o}: not fitted yet (run `sead-sexdiff stage-de-limma --omic {o}`); the summary lacks it")

            hashes_now = {
                group: {name: _sha256_file(p) for name, p in files.items()} for group, files in prep["paths"].items()
            }
            if hashes_now != prep["hashes"]:
                changed = [f"{g}/{n}" for g in hashes_now for n in hashes_now[g] if hashes_now[g][n] != prep["hashes"][g][n]]
                raise ValueError(f"input file(s) {changed} changed during the stage-de-limma run of {c}")
            dq = prep["deseq2"]
            deseq2_now = {str(p.relative_to(dq["dir"])): _sha256_file(p) for p in sorted(dq["dir"].rglob("*")) if p.is_file()}
            if deseq2_now != dq["files"]:
                raise ValueError(f"the DESeq2 outputs under {dq['dir']} changed during the stage-de-limma run of {c}")

            reference = prep["reference"]
            ref_record = reference["record"]
            any_fit = {o: records[o][next(iter(STAGE_DE_MODELS))]["fit"] for o in included}
            manifest = {
                "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "cell_type": cell_type,
                "cohort": c,
                "n_donors": len(prep["donors"]),
                "cells": prep["cells"],
                "k": prep["k"],
                "omics": included,
                "complete": not missing,
                "inputs": {
                    "files": {
                        group: {name: {"path": str(p), "sha256": prep["hashes"][group][name]} for name, p in files.items()}
                        for group, files in prep["paths"].items()
                    },
                    "pseudobulk_manifest_sha256": prep["hashes"]["pseudobulk"]["manifest.json"],
                    "export_manifest_sha256": prep["hashes"]["export"]["manifest.json"],
                    "export_input": prep["loaded"]["export_manifest"].get("input"),
                    "deseq2": {
                        "dir": str(dq["dir"]),
                        "manifest_sha256": dq["manifest_sha256"],
                        "summary_sha256": dq["summary_sha256"],
                        "pydeseq2": dq["software"].get("pydeseq2"),
                        "files": dq["files"],
                        "note": (
                            "complete for both omics, made from the inputs now on disk (every inputs.files sha256), "
                            "every table hashing as recorded and the summary and concordance cells equal to a "
                            "recomputation from the tables; hashed before and after the run, never written"
                        ),
                    },
                    "note": (
                        "every file hashed before it was read and again after the run; the counts are bound to the "
                        "export by stage-de's reproduction (export_reproduction) and by edgeR's TMM against the export "
                        "manifest's factors (tmm_check)"
                    ),
                },
                "export_reproduction": {
                    o: prep["reproduction"].get(o, records[o][next(iter(STAGE_DE_MODELS))]["export_reproduction"])
                    for o in included
                },
                "reference": {
                    "dir": str(reference["dir"]),
                    "files": reference["files"],
                    "result_sha256": reference["files"]["result.json"]["sha256"],
                    "delta_rejected": prep["delta_rejected"],
                    "alpha": ANALYSIS_ALPHA,
                    "rejection_rule": ANALYSIS_REJECTION_RULE,
                    "component_selection": ref_record["latent"]["component_selection"],
                    "selected_lv": ref_record["latent"]["selected_lv"],
                    "motco": reference["manifest"]["software"]["motco"],
                    "seed": reference["manifest"]["measurement"]["seed"],
                    "permutations": reference["manifest"]["measurement"]["permutations"],
                    "note": (
                        "read and checked (export hashes, MOTCO version, seed, permutations, CV selection), never "
                        "rerun or modified; only its delta rejection status enters the concordance rules"
                    ),
                },
                "models": {
                    m: {
                        "formula": spec["formula"],
                        "formulas": {o: stage_de_formula(m, o) for o in _STAGE_DE_OMICS},
                        "stage_coding": spec["stage_coding"],
                        "interaction": spec["interaction"],
                        "tested": [{"coefficient": label, "design_column": column} for label, column in spec["tested"]],
                        "joint": (
                            {"test": STAGE_DE_LIMMA_JOINT_TESTS[m][0], "coefficients": list(STAGE_DE_LIMMA_JOINT_TESTS[m][1])}
                            if m in STAGE_DE_LIMMA_JOINT_TESTS else None
                        ),
                    }
                    for m, spec in STAGE_DE_MODELS.items()
                },
                "method": STAGE_DE_LIMMA_METHOD,
                "log_cpm": (
                    "voom's log2((count + 0.5) / (lib.size × TMM factor + 1) × 1e6), not the export's log2-CPM with "
                    "prior count 2 scaled by the effective library size"
                ),
                "joint_tests": {m: {"test": j, "coefficients": list(labels)} for m, (j, labels) in STAGE_DE_LIMMA_JOINT_TESTS.items()},
                "joint_f": (
                    "limma's moderated F (topTable over both coefficients): the moderated t's combined through the "
                    "correlation of the unweighted design's coefficients (cov.coefficients), df1 = 2, df2 = df.prior + "
                    "df.residual; approximate with voom's observation weights, measured per feature by the exact-F "
                    "diagnostic (fits.*.joint.*.exact), which never replaces it"
                ),
                "union": (
                    "a factor model's two moderated t tests are also counted as their union (significant on either); "
                    "the union carries no joint error-rate claim and is reported for the like-for-like comparison with "
                    "the DESeq2 union only"
                ),
                "covariates": {
                    "per_omic": {o: list(cols) for o, cols in STAGE_DE_COVARIATES.items()},
                    "scale": "on their recorded scale, not centered or scaled; nuclei counts on the log10 scale",
                    "source": str(input_dir / c / "covariates.csv"),
                    "missing_policy": ANALYSIS_MISSING_COVARIATES,
                    "imputed": False,
                },
                "reference_levels": dict(STAGE_DE_REFERENCE_LEVELS),
                "alpha": STAGE_DE_ALPHA,
                "significance": (
                    "adj.P.Val <= alpha; Benjamini-Hochberg over all tested features (no independent filter); DESeq2's "
                    "BH family is the tested features passing pydeseq2's per-coefficient independent filter"
                ),
                "tmm_atol": STAGE_DE_LIMMA_TMM_ATOL,
                "tmm_check": {
                    o: {
                        "max_abs_diff": any_fit[o]["tmm"]["max_abs_diff"], "atol": STAGE_DE_LIMMA_TMM_ATOL,
                        "passed": bool(any_fit[o]["tmm"]["passed"]), "n_donors": any_fit[o]["n_donors"],
                        "n_features": any_fit[o]["n_features"],
                    }
                    for o in included
                },
                "design_checks": {
                    o: {
                        m: {
                            "python": {k: v for k, v in prep["design_checks"][o][m].items() if k != "formula"},
                            "r": records[o][m]["fit"]["design_check"],
                        }
                        for m in STAGE_DE_MODELS
                    }
                    for o in included
                },
                "feature_sets": {
                    o: {
                        **prep["feature_sets"][o],
                        "statuses": list(STAGE_DE_FEATURE_STATUSES),
                        "exported_status": STAGE_DE_EXPORTED_STATUS,
                        "source": str(input_dir / c / f"features_{o}.csv"),
                    }
                    for o in _STAGE_DE_OMICS
                },
                "fits": {o: {m: records[o][m]["fit"] for m in STAGE_DE_MODELS} for o in included},
                "tables": {o: {m: records[o][m]["tables"] for m in STAGE_DE_MODELS} for o in included},
                "concordance_rule": dict(STAGE_DE_LIMMA_CONCORDANCE_RULE),
                "deseq2_concordance_rule": dict(STAGE_DE_CONCORDANCE_RULE),
                "concordance": concordance,
                "run": {"omics_fitted": [f"{o}/{m}" for o, m in fitted], "skipped": [f"{o}/{m}" for o, m in skipped]},
                "outputs": sorted([STAGE_DE_LIMMA_SUMMARY, STAGE_DE_LIMMA_COMPARISON, STAGE_DE_LIMMA_MARKDOWN,
                                   STAGE_DE_LIMMA_MANIFEST]),
                "software": software,
                "git": _git_revision(),
                "findings": findings,
            }
            cohort_dir.mkdir(parents=True, exist_ok=True)
            for name, frame in ((STAGE_DE_LIMMA_SUMMARY, summary), (STAGE_DE_LIMMA_COMPARISON, comparison)):
                tmp = cohort_dir / f"{name}.tmp"
                frame.to_csv(tmp, index=False, float_format=_csv_float)
                os.replace(tmp, cohort_dir / name)
            written = pd.read_csv(cohort_dir / STAGE_DE_LIMMA_SUMMARY, dtype={"design_column": str}, float_precision="round_trip")
            written_cmp = pd.read_csv(cohort_dir / STAGE_DE_LIMMA_COMPARISON, float_precision="round_trip")
            tmp = cohort_dir / f"{STAGE_DE_LIMMA_MARKDOWN}.tmp"
            tmp.write_text(render_stage_de_limma(written, written_cmp, manifest))
            os.replace(tmp, cohort_dir / STAGE_DE_LIMMA_MARKDOWN)
            tmp = cohort_dir / f"{STAGE_DE_LIMMA_MANIFEST}.tmp"
            tmp.write_text(json.dumps(manifest, indent=2) + "\n")
            os.replace(tmp, cohort_dir / STAGE_DE_LIMMA_MANIFEST)
            print(f"[{log}] Wrote {cohort_dir}")
            for msg in findings:
                print(f"[{log}] FINDING: {msg}")
            results[c] = {
                "dir": cohort_dir, "summary": summary, "comparison": comparison, "manifest": manifest,
                "concordance": concordance, "findings": findings, "fitted": fitted, "skipped": skipped,
            }
            del prepared[c]
            gc.collect()
    finally:
        if debug is None:
            shutil.rmtree(work_root, ignore_errors=True)
    return results


# ---------------------------------------------------------------------------
# Cell-type positive control (microglia against astrocytes)
# ---------------------------------------------------------------------------

# Roles of the two cell types: the control has a documented ADNC-stage
# response; the reference is the case-study cell type, whose outputs are read
# and never rerun.
POSITIVE_CONTROL_CELL_TYPES: dict[str, str] = {"control": "microglia", "reference": "astrocyte"}
POSITIVE_CONTROL_PRIMARY_COHORT = "cohort_30"
# ``source`` values of axis_correlations.csv that are technical variables for
# E3: the technical table's columns and the export's TMM factors and library
# sizes (DIAGNOSTIC_VARIABLES' ``technical`` and ``export`` groups). The six
# analysis covariates (covariates.csv) and ``male`` (metadata.csv) are not.
POSITIVE_CONTROL_TECHNICAL_SOURCES = ("donor_technical.csv", "manifest.json (export, tmm)")
# Written before any microglia run (proposal, design D7) and never changed
# from results. No expectation reads a p-value, a rejection status or a
# trajectory statistic; the values are compared, not tested.
POSITIVE_CONTROL_EXPECTATIONS: dict[str, dict] = {
    "E1": {
        "name": "stage signal",
        "statement": (
            "under each stage-DE method, the number of exported RNA features with a stage effect at FDR "
            f"{STAGE_DE_ALPHA} is larger for the control than for the reference, compared within the method"
        ),
        "methods": {
            "deseq2": {
                "command": "stage-de", "subdir": CELL_TYPE_SUBDIRS["stage_de"], "summary": STAGE_DE_SUMMARY,
                "manifest": STAGE_DE_MANIFEST, "omic": "rna", "model": "factor_additive",
                "coefficient": STAGE_DE_UNION, "feature_set": "exported", "value": "n_significant",
            },
            "limma": {
                "command": "stage-de-limma", "subdir": CELL_TYPE_SUBDIRS["stage_de_limma"],
                "summary": STAGE_DE_LIMMA_SUMMARY, "manifest": STAGE_DE_LIMMA_MANIFEST, "omic": "rna",
                "model": "factor_additive", "coefficient": STAGE_DE_LIMMA_JOINT_TESTS["factor_additive"][0],
                "feature_set": "exported", "value": "n_significant",
            },
        },
        "comparison": "control > reference, within each method",
        "rule": (
            "holds only when it holds under both methods; when it holds under one only it is method_dependent, "
            "which counts as not holding, and the method under which it held is named"
        ),
    },
    "E2": {
        "name": "latent separation",
        "statement": "the control's primary analysis has a larger CV mean AUROC than the reference's",
        "subdir": CELL_TYPE_SUBDIRS["motco_analysis"],
        "file": "result.json",
        "value": "latent.cv_mean_auroc",
        "comparison": "control > reference",
    },
    "E3": {
        "name": "axis identity",
        "statement": (
            "along the first latent component of the control's primary analysis, |Spearman| with stage exceeds "
            "|Spearman| with every technical variable of the diagnostic's variable set"
        ),
        "subdir": CELL_TYPE_SUBDIRS["motco_diagnostics"],
        "file": DIAGNOSTIC_OUTPUTS["axis"],
        "analysis": DIAGNOSTIC_PRIMARY_LABEL,
        "component": "pls_0",
        "value": "spearman (raw, not partial)",
        "stage_variable": "stage",
        "technical_sources": list(POSITIVE_CONTROL_TECHNICAL_SOURCES),
        "comparison": "|rho(stage)| > max |rho(technical)|, control only",
    },
}
POSITIVE_CONTROL_RULE = (
    "the control holds when E1, E2 and E3 all hold in the primary cohort (POSITIVE_CONTROL_PRIMARY_COHORT); if it "
    "holds, the pipeline resolves stage structure on real data and the reference's null is not a pipeline failure; "
    "if not, the latent construction is in question for any cell type"
)
# The analysis and the diagnostic E2 and E3 read, per cell type: the primary
# analysis (the default; today's inputs, output path and bytes) or the
# latent_adjusted arm of latent-sensitivity and the diagnostic made on it
# (outputs under <out_dir>/latent_adjusted/<cohort>/). E1 is read unchanged
# on every basis (basis-independent). A registry, never a free directory,
# never extended from results; the expectations, their comparisons, E3's
# technical sources and the all-three rule do not change with the basis.
# The restated control's reading was written before any microglia arm ran
# (proposal of mtg-microglia-latent-sensitivity) and is never changed from
# results. Flag and ask before changing.
POSITIVE_CONTROL_DEFAULT_BASIS = "primary"
POSITIVE_CONTROL_BASES: dict[str, dict] = {
    "primary": {
        "description": "the primary analysis of each cell type and the diagnostic made from it",
        "analysis_subdir": CELL_TYPE_SUBDIRS["motco_analysis"],
        "diagnostics_subdir": CELL_TYPE_SUBDIRS["motco_diagnostics"],
        "arm": None,
        "reading": None,
    },
    "latent_adjusted": {
        "description": "the latent_adjusted arm of latent-sensitivity of each cell type and the diagnostic made from it",
        "analysis_subdir": f"{CELL_TYPE_SUBDIRS['motco_latent_sensitivity']}/latent_adjusted/{LATENT_SENSITIVITY_ANALYSIS_DIR}",
        "diagnostics_subdir": f"{CELL_TYPE_SUBDIRS['motco_latent_sensitivity']}/diagnostics/latent_adjusted",
        "arm": "latent_adjusted",
        "reading": (
            "a diagnosis of the primary-basis control, which stays the control's outcome: if the restated control "
            "holds in the primary cohort, the measured quality gradient explains the primary control's failure (once "
            "the library-quality direction is conditioned out of the basis, the pipeline resolves stage in the cell "
            "type with the stronger stage signal), without promoting the arm (MOTCO's power comparison, "
            "candidate_fails, excludes the latent adjustment as a replacement of the primary), changing the "
            "primary-basis outcome or changing the reference cell type's reading; if it does not hold, the latent "
            "construction is in question beyond the quality gradient (the composition of the control population, the "
            "rank-1 selection and the control's ATAC nuclei per donor are the stated candidates) and the diagnosis "
            "'technical gradient' is not supported. E3 on this basis carries the mediator caveat in both directions: "
            "the adjustment removes the quality direction and, with it, any stage signal the added fractions carry, "
            "so the stage |rho| may rise or fall; a failure for that reason is still a failure and the expectation "
            "is not re-tuned. The sensitivity cohort is descriptive, as on the primary basis"
        ),
    },
}
# Only an arm with a latent adjustment of its own can be a basis: with the
# covariate-blind basis (extended_covariates) the latent space, CV rank and
# AUROC equal the primary's by MOTCO's contract, so E2 and E3 would repeat
# the primary's values; that equality is stated and checked in the report.
for _name, _entry in POSITIVE_CONTROL_BASES.items():
    _arm = _entry["arm"]
    if _name == POSITIVE_CONTROL_DEFAULT_BASIS:
        if _arm is not None:
            raise AssertionError("the default positive-control basis is the primary analysis, not an arm")
        continue
    if _arm not in LATENT_SENSITIVITY_ARMS or not LATENT_SENSITIVITY_ARMS[_arm]["latent_adjustment"]:
        raise AssertionError(
            f"positive-control basis {_name!r} must be an arm of LATENT_SENSITIVITY_ARMS with a latent adjustment; "
            f"got {_arm!r}"
        )
    if _entry["analysis_subdir"] != (
        f"{CELL_TYPE_SUBDIRS['motco_latent_sensitivity']}/{_arm}/{LATENT_SENSITIVITY_ANALYSIS_DIR}"
    ) or _entry["diagnostics_subdir"] != f"{CELL_TYPE_SUBDIRS['motco_latent_sensitivity']}/diagnostics/{_arm}":
        raise AssertionError(f"positive-control basis {_name!r} must point at the {_arm} arm's layout")
del _name, _entry, _arm
POSITIVE_CONTROL_TABLE = "positive_control.csv"
POSITIVE_CONTROL_MARKDOWN = "positive_control.md"
POSITIVE_CONTROL_MANIFEST = "positive_control_manifest.json"
POSITIVE_CONTROL_COLUMNS = (
    "cohort", "primary_cohort", "expectation", "method", "quantity", "control", "control_value", "reference",
    "reference_value", "comparison", "holds", "outcome", "detail",
)


def _positive_control_read(
    root: Path, cohort: str, role: str, *, diagnostics: bool, basis: str = POSITIVE_CONTROL_DEFAULT_BASIS
) -> dict:
    """Read and check one root's files of one cohort for the positive control (design D7).

    Checks the cell type of every manifest read (``POSITIVE_CONTROL_CELL_TYPES[role]``)
    and the chain: both stage-DE steps were made from the primary analysis
    (``result.json`` sha256) and the limma step from the DESeq2 summary read
    for E1; the diagnostic read for E3 from the analysis read for E2. On a
    non-primary ``basis`` (``POSITIVE_CONTROL_BASES``) E2 and E3 read the
    arm analysis and the diagnostic made on it, and the chain also checks
    that the arm was made from the primary export now on disk (its recorded
    input hashes, the check ``axis-diagnostics`` runs), that the arm's
    comparison manifest names the primary analysis now on disk as its
    reference (so E1, E2 and E3 hang off one primary analysis) and that it
    lists this arm directory, and that the arm's conditioning is the
    registered one.

    Returns
    -------
    dict
        ``files`` (label → path and sha256), ``findings``, ``basis``,
        ``arm`` (``None`` on the primary basis; the arm's directory, hashes
        and adjustment block otherwise), ``e1`` (method → the summary row's
        ``n_significant`` and ``n_tested``), ``e2`` (``cv_mean_auroc``) and,
        with ``diagnostics``, ``e3`` (the ``pls_0`` rows of
        ``axis_correlations.csv`` labelled ``primary`` by the diagnostic).
    """
    cell_type = POSITIVE_CONTROL_CELL_TYPES[role]
    flag = "" if cell_type == DEFAULT_CELL_TYPE else f" --cell-type {cell_type}"
    entry = POSITIVE_CONTROL_BASES[basis]
    files, findings = {}, []

    def need(path: Path, command: str) -> Path:
        if not path.exists():
            raise FileNotFoundError(
                f"{path} is missing ({role} root, cell type {cell_type}); run `sead-sexdiff {command}{flag}` to "
                "produce it; nothing was written"
            )
        files[str(path.relative_to(root))] = {"path": str(path), "sha256": _sha256_file(path)}
        return path

    def manifest(path: Path, command: str) -> dict:
        m = json.loads(need(path, command).read_text())
        findings.extend(check_manifest_cell_type(m, path, cell_type))
        return m

    # The primary analysis: E1's stage-DE steps are chained to it on every basis
    an_dir = root / CELL_TYPE_SUBDIRS["motco_analysis"] / cohort
    primary_path = need(an_dir / "result.json", "analyze")
    an_manifest = manifest(an_dir / "analysis_manifest.json", "analyze")
    primary_record = json.loads(primary_path.read_text())
    if primary_record.get("latent", {}).get("component_selection") != "cv":
        raise ValueError(f"{primary_path} is not a CV-selected analysis; the positive control reads the primary analysis")
    primary_sha = files[str(primary_path.relative_to(root))]["sha256"]
    export_sha = an_manifest.get("inputs", {}).get("export_manifest_sha256")
    arm_info = None
    if entry["arm"] is None:
        result_path, record, result_sha = primary_path, primary_record, primary_sha
        diag_command = "axis-diagnostics"
    else:
        arm = entry["arm"]
        arm_command = "latent-sensitivity"
        arm_dir = root / entry["analysis_subdir"] / cohort
        result_path = need(arm_dir / "result.json", arm_command)
        arm_manifest = manifest(arm_dir / "analysis_manifest.json", arm_command)
        record = json.loads(result_path.read_text())
        if record.get("latent", {}).get("component_selection") != "cv":
            raise ValueError(f"{result_path} is not a CV-selected analysis; the {arm} arm keeps the primary's CV rank selection")
        conditioning = (record.get("latent_adjustment") or {}).get("conditioning")
        if conditioning != LATENT_SENSITIVITY_ARMS[arm]["latent_adjustment"]:
            raise ValueError(
                f"{result_path} records latent-adjustment conditioning {conditioning!r}, not "
                f"{LATENT_SENSITIVITY_ARMS[arm]['latent_adjustment']!r} (the {arm} arm of LATENT_SENSITIVITY_ARMS); "
                "nothing was written"
            )
        result_sha = files[str(result_path.relative_to(root))]["sha256"]
        # The arm was made from the primary export now on disk (the check axis-diagnostics runs) ...
        export_dir = root / CELL_TYPE_SUBDIRS["motco_input"] / cohort
        recorded = arm_manifest.get("inputs", {}).get("files", {})
        if not recorded:
            raise ValueError(f"{arm_dir / 'analysis_manifest.json'} records no input hashes; nothing was written")
        for name, info in recorded.items():
            path = export_dir / name
            now = _sha256_file(path) if path.exists() else None
            if info.get("sha256") != now:
                raise ValueError(
                    f"{path} differs from the file the {arm} arm in {arm_dir} was made from (sha256 "
                    f"{info.get('sha256')} recorded, {now} now); rerun `sead-sexdiff {arm_command}{flag}` or restore "
                    "the export; nothing was written"
                )
            files[str(path.relative_to(root))] = {"path": str(path), "sha256": now}
        if arm_manifest.get("inputs", {}).get("export_manifest_sha256") != export_sha:
            raise ValueError(
                f"{arm_dir / 'analysis_manifest.json'} and {an_dir / 'analysis_manifest.json'} record different export "
                "manifests; nothing was written"
            )
        # ... and compared with the primary analysis now on disk, the one E1's stage-DE steps are checked against
        comp_dir = root / CELL_TYPE_SUBDIRS["motco_latent_sensitivity"] / LATENT_SENSITIVITY_COMPARISON_DIR / cohort
        comp_path = comp_dir / LATENT_SENSITIVITY_MANIFEST
        comp_manifest = manifest(comp_path, arm_command)
        comp_ref = comp_manifest.get("reference", {}).get("files", {}).get("result.json", {}).get("sha256")
        if comp_ref != primary_sha:
            raise ValueError(
                f"{comp_path} compared the arms with a primary analysis whose result.json sha256 is {comp_ref}, not "
                f"{primary_path} ({primary_sha}), the analysis E1's stage-DE steps are checked against; rerun "
                f"`sead-sexdiff {arm_command}{flag}`; nothing was written"
            )
        comp_arm = comp_manifest.get("arms", {}).get(arm, {}).get("files", {}).get("result.json", {}).get("sha256")
        if comp_arm != result_sha:
            raise ValueError(
                f"{comp_path} records result.json sha256 {comp_arm} for the {arm} arm, not {result_path} "
                f"({result_sha}); rerun `sead-sexdiff {arm_command}{flag}`; nothing was written"
            )
        diag_command = (
            f"axis-diagnostics --analysis-dir {root / entry['analysis_subdir']} "
            f"--out-dir {root / entry['diagnostics_subdir']}"
        )
        arm_info = {
            "arm": arm,
            "dir": str(arm_dir),
            "result_sha256": result_sha,
            "selected_lv": record["latent"]["selected_lv"],
            "covariates": (record.get("covariates") or {}).get("columns"),
            "latent_adjustment": record.get("latent_adjustment"),
            "comparison_manifest": {"path": str(comp_path), "sha256": files[str(comp_path.relative_to(root))]["sha256"]},
            "primary": {"dir": str(an_dir), "result_sha256": primary_sha},
        }
    auroc = record["latent"].get("cv_mean_auroc")
    if auroc is None:
        raise ValueError(f"{result_path} records no latent.cv_mean_auroc")

    e1, summaries = {}, {}
    for method, spec in POSITIVE_CONTROL_EXPECTATIONS["E1"]["methods"].items():
        d = root / spec["subdir"] / cohort
        summary_path = need(d / spec["summary"], spec["command"])
        m = manifest(d / spec["manifest"], spec["command"])
        if not m.get("complete") or spec["omic"] not in m.get("omics", []):
            raise ValueError(f"{d / spec['manifest']} is not complete for {spec['omic']}; run `sead-sexdiff {spec['command']}{flag}`")
        ref_sha = m.get("reference", {}).get("result_sha256")
        if ref_sha != primary_sha:
            raise ValueError(
                f"{d / spec['manifest']} was made against an analysis with result.json sha256 {ref_sha}, not "
                f"{primary_path} ({primary_sha}); rerun `sead-sexdiff {spec['command']}{flag}`"
            )
        if m.get("inputs", {}).get("export_manifest_sha256") != export_sha:
            raise ValueError(f"{d / spec['manifest']} and {an_dir / 'analysis_manifest.json'} record different export manifests")
        summary = pd.read_csv(summary_path, dtype={"design_column": str})
        summaries[method] = files[str(summary_path.relative_to(root))]["sha256"]
        row = summary[
            (summary["omic"] == spec["omic"]) & (summary["model"] == spec["model"])
            & (summary["coefficient"] == spec["coefficient"]) & (summary["feature_set"] == spec["feature_set"])
        ]
        if len(row) != 1:
            raise ValueError(
                f"{summary_path} has {len(row)} row(s) for {spec['omic']} {spec['model']} {spec['coefficient']} "
                f"({spec['feature_set']}); expected one"
            )
        e1[method] = {"n_significant": int(row[spec["value"]].iloc[0]), "n_tested": int(row["n_tested"].iloc[0])}
        if method == "limma":
            recorded = m.get("inputs", {}).get("deseq2", {}).get("summary_sha256")
            if recorded != summaries["deseq2"]:
                raise ValueError(
                    f"{d / spec['manifest']} checked a DESeq2 summary with sha256 {recorded}, not the one read for E1 "
                    f"({summaries['deseq2']}); rerun `sead-sexdiff stage-de-limma{flag}`"
                )
    out = {"files": files, "findings": findings, "basis": basis, "arm": arm_info, "e1": e1, "e2": float(auroc)}

    if diagnostics:
        spec = POSITIVE_CONTROL_EXPECTATIONS["E3"]
        d = root / entry["diagnostics_subdir"] / cohort
        axis_path = need(d / spec["file"], diag_command)
        m = manifest(d / DIAGNOSTIC_OUTPUTS["manifest"], diag_command)
        primary = [a for a in m.get("inputs", {}).get("analyses", []) if a.get("label") == spec["analysis"]]
        if len(primary) != 1 or primary[0].get("files", {}).get("result.json", {}).get("sha256") != result_sha:
            raise ValueError(
                f"{d / DIAGNOSTIC_OUTPUTS['manifest']} was not made from {result_path}; rerun "
                f"`sead-sexdiff {diag_command}{flag}`"
            )
        if m.get("inputs", {}).get("technical", {}).get("manifest", {}).get("matches_export") is not True:
            raise ValueError(
                f"{d / DIAGNOSTIC_OUTPUTS['manifest']}: the technical table was not made from the export's pseudobulk "
                "(inputs.technical.manifest.matches_export is not true); E3 needs the cell type's own technical table"
            )
        axis = pd.read_csv(axis_path, dtype={"variable": str, "source": str, "analysis": str, "component": str})
        out["e3"] = axis[(axis["analysis"] == spec["analysis"]) & (axis["component"] == spec["component"])].reset_index(drop=True)
    return out


def _positive_control_rows(cohort: str, control: dict, reference: dict) -> tuple[list[dict], dict]:
    """The rows of ``positive_control.csv`` for one cohort and the outcomes (E1 per method, E1, E2, E3, overall)."""
    types = POSITIVE_CONTROL_CELL_TYPES
    primary = cohort == POSITIVE_CONTROL_PRIMARY_COHORT
    rows, outcomes, findings = [], {}, []
    base = {"cohort": cohort, "primary_cohort": primary}
    arm = POSITIVE_CONTROL_BASES[control["basis"]]["arm"]
    if reference["basis"] != control["basis"]:
        raise ValueError(f"the control was read on basis {control['basis']!r} and the reference on {reference['basis']!r}")
    analysis_name = "the primary analysis" if arm is None else f"the {arm} arm analysis (basis {control['basis']})"

    held = {}
    for method, spec in POSITIVE_CONTROL_EXPECTATIONS["E1"]["methods"].items():
        a, b = control["e1"][method]["n_significant"], reference["e1"][method]["n_significant"]
        held[method] = a > b
        rows.append({
            **base, "expectation": "E1", "method": method,
            "quantity": f"exported {spec['omic'].upper()} features at FDR {STAGE_DE_ALPHA}: {spec['model']} {spec['coefficient']}",
            "control": types["control"], "control_value": a, "reference": types["reference"], "reference_value": b,
            "comparison": ">", "holds": held[method], "outcome": "holds" if held[method] else "fails",
            "detail": f"tested {control['e1'][method]['n_tested']} (control), {reference['e1'][method]['n_tested']} (reference)",
        })
    if all(held.values()):
        e1 = "holds"
        detail = "holds under both methods"
    elif any(held.values()):
        e1 = "method_dependent"
        detail = "held under " + ", ".join(m for m, h in held.items() if h) + " only; counts as not holding"
    else:
        e1 = "fails"
        detail = "holds under neither method"
    outcomes["E1"] = e1
    outcomes["E1_methods"] = held
    rows.append({
        **base, "expectation": "E1", "method": "both", "quantity": "E1 under both stage-DE methods",
        "control": "", "control_value": None, "reference": "", "reference_value": None,
        "comparison": "both hold", "holds": e1 == "holds", "outcome": e1, "detail": detail,
    })

    a, b = control["e2"], reference["e2"]
    outcomes["E2"] = "holds" if a > b else "fails"
    rows.append({
        **base, "expectation": "E2", "method": "", "quantity": f"CV mean AUROC of {analysis_name} (latent.cv_mean_auroc)",
        "control": types["control"], "control_value": a, "reference": types["reference"], "reference_value": b,
        "comparison": ">", "holds": a > b, "outcome": outcomes["E2"], "detail": "",
    })

    spec = POSITIVE_CONTROL_EXPECTATIONS["E3"]
    t = control["e3"]
    stage = t[t["variable"] == spec["stage_variable"]]
    if len(stage) != 1 or not np.isfinite(stage["spearman"].iloc[0]):
        raise ValueError(f"axis_correlations.csv has no single finite {spec['stage_variable']} row for {cohort} {spec['component']}")
    technical = t[t["source"].isin(spec["technical_sources"])]
    if technical.empty:
        raise ValueError(f"axis_correlations.csv has no technical variable row for {cohort} {spec['component']}")
    finite = technical[np.isfinite(technical["spearman"].to_numpy(dtype=float))]
    if len(finite) < len(technical):
        findings.append(
            f"{cohort}: E3 technical variable(s) with an undefined Spearman correlation are not compared: "
            f"{sorted(technical.loc[~technical.index.isin(finite.index), 'variable'])}"
        )
    abs_rho = finite["spearman"].abs()
    top = finite.loc[abs_rho.idxmax()]
    a, b = float(abs(stage["spearman"].iloc[0])), float(abs(top["spearman"]))
    outcomes["E3"] = "holds" if a > b else "fails"
    outcomes["E3_technical_max"] = {"variable": str(top["variable"]), "source": str(top["source"]), "abs_spearman": b}
    rows.append({
        **base, "expectation": "E3", "method": "",
        "quantity": f"|Spearman| of {spec['component']} ("
                    + (spec["analysis"] if arm is None else f"{arm} arm, labelled {spec['analysis']} by its diagnostic")
                    + ") with stage vs the largest with a technical variable",
        "control": spec["stage_variable"], "control_value": a, "reference": str(top["variable"]), "reference_value": b,
        "comparison": ">", "holds": a > b, "outcome": outcomes["E3"],
        "detail": f"{len(finite)} technical variables compared ({', '.join(spec['technical_sources'])})",
    })

    failing = [e for e in ("E1", "E2", "E3") if outcomes[e] != "holds"]
    overall = not failing
    outcomes["all_hold"] = overall
    outcomes["failing"] = failing
    rows.append({
        **base, "expectation": "overall", "method": "", "quantity": "E1, E2 and E3 all hold",
        "control": "", "control_value": None, "reference": "", "reference_value": None, "comparison": "all hold",
        "holds": overall, "outcome": "holds" if overall else "does_not_hold",
        "detail": ("the primary cohort: this row is the control's outcome" if primary
                   else f"not the primary cohort ({POSITIVE_CONTROL_PRIMARY_COHORT}); descriptive")
                  + ("" if overall else "; not holding: " + ", ".join(failing)),
    })
    return rows, {**outcomes, "findings": findings}


def render_positive_control(table: pd.DataFrame, manifest: dict) -> str:
    """Render ``positive_control.md`` from the written table and the manifest."""
    c = manifest["cohort"]
    types = manifest["cell_types"]
    basis = manifest.get("basis", POSITIVE_CONTROL_DEFAULT_BASIS)
    entry = POSITIVE_CONTROL_BASES[basis]
    overall = table[table["expectation"] == "overall"].iloc[0]
    lines = [
        f"# Positive control: {types['control']} against {types['reference']} ({c})"
        + ("" if entry["arm"] is None else f", basis `{basis}`"),
        "",
        f"Predeclared expectations (POSITIVE_CONTROL_EXPECTATIONS), evaluated from the generated files of both cell "
        f"types. {POSITIVE_CONTROL_RULE[0].upper() + POSITIVE_CONTROL_RULE[1:]}. No expectation reads a p-value, a "
        "rejection status or a trajectory statistic; the values are compared, not tested.",
        "",
        f"This is a positive control for the stage signal, not for the sex × stage hypothesis: the "
        f"{types['control']} trajectory result is read under the same binding rules as the {types['reference']} one, "
        f"and nothing here changes any {types['reference']} output or reading.",
        "",
    ]
    if entry["arm"] is not None:
        lines += [
            f"**Basis `{basis}`** (POSITIVE_CONTROL_BASES): E2 and E3 read {entry['description']}; E1 is "
            "basis-independent and equals the primary-basis rows. The expectations, their comparisons, E3's technical "
            f"variables and the all-three rule are the primary basis's. Read as {entry['reading']}.",
            "",
        ]
    lines += [
        ("**Primary cohort.** " if manifest["primary_cohort"] else f"**Not the primary cohort** "
         f"({POSITIVE_CONTROL_PRIMARY_COHORT} is): descriptive. ")
        + f"Overall: **{overall['outcome'].replace('_', ' ')}**"
        + ("" if bool(overall["holds"]) else f" ({', '.join(manifest['outcomes']['failing'])} not holding)") + ".",
        "",
        "| Expectation | Method | Quantity | Control | Value | Compared with | Value | Holds | Outcome | Detail |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]

    def cell(v, spec=".6g"):
        if v is None or (isinstance(v, float) and not np.isfinite(v)):
            return ""
        return format(v, spec) if isinstance(v, float) else str(v).replace("|", "\\|")

    for r in table.to_dict("records"):
        lines.append(
            f"| {r['expectation']} | {cell(r['method'])} | {cell(r['quantity'])} | {cell(r['control'])} | "
            f"{cell(r['control_value'])} | {cell(r['reference'])} | {cell(r['reference_value'])} | "
            f"{'yes' if r['holds'] else 'no'} | {r['outcome']} | {cell(r['detail'])} |"
        )
    lines += ["", "## Expectations", ""]
    for key, spec in POSITIVE_CONTROL_EXPECTATIONS.items():
        lines.append(f"- **{key} ({spec['name']})**: {spec['statement']}." + (f" {spec['rule'][0].upper() + spec['rule'][1:]}." if "rule" in spec else ""))
    lines += ["", "## Inputs", ""]
    for role in ("control", "reference"):
        info = manifest["inputs"][role]
        lines.append(f"- {role} ({types[role]}): `{info['root']}`")
        for rel, f in info["files"].items():
            lines.append(f"  - `{rel}` sha256 `{f['sha256'][:12]}…`")
    if manifest["findings"]:
        lines += ["", "## Findings", ""] + [f"- {msg}" for msg in manifest["findings"]]
    return "\n".join(lines) + "\n"


def run_positive_control(
    cohort: str = "all",
    *,
    basis: str = POSITIVE_CONTROL_DEFAULT_BASIS,
    control_root: str | Path | None = None,
    reference_root: str | Path | None = None,
    out_dir: str | Path | None = None,
    force: bool = False,
) -> dict:
    """Evaluate the predeclared positive-control expectations E1–E3 per cohort (design D7).

    Parameters
    ----------
    cohort:
        ``"cohort_30"``, ``"cohort_50"`` or ``"all"``.
    basis:
        A name of ``POSITIVE_CONTROL_BASES``: ``primary`` (the default; E2
        and E3 from each cell type's primary analysis and the diagnostic made
        from it, outputs under ``<out_dir>/<cohort>/``) or ``latent_adjusted``
        (E2 and E3 from each cell type's ``latent_adjusted`` arm of
        ``latent-sensitivity`` and the diagnostic made on that arm, outputs
        under ``<out_dir>/<basis>/<cohort>/``). E1 is read unchanged on every
        basis. The expectations, their comparisons, E3's technical sources
        and the all-three rule do not change with the basis.
    control_root:
        Results root of the control cell type (``POSITIVE_CONTROL_CELL_TYPES["control"]``,
        default its registry root ``results/microglia``): ``stage_de/``,
        ``stage_de_limma/``, ``motco_analysis/`` and ``motco_diagnostics/``
        are read; every manifest must record the control cell type.
    reference_root:
        Results root of the reference cell type (default ``results``):
        ``stage_de/``, ``stage_de_limma/`` and ``motco_analysis/`` are read;
        every manifest must record the reference cell type (a manifest with
        no ``cell_type`` reads as astrocyte, with a finding).
    out_dir:
        Output root; each cohort goes to ``<out_dir>/<cohort>/`` on the
        primary basis and to ``<out_dir>/<basis>/<cohort>/`` otherwise
        (default ``<control_root>/positive_control``).
    force:
        Replace existing cohort directories.

    Returns
    -------
    dict
        Cohort → ``{"dir", "table", "manifest", "outcomes"}``.

    Raises
    ------
    FileNotFoundError
        A summary, analysis or diagnostic file is missing; the message names
        it and the command that produces it.
    ValueError
        ``basis`` is not registered; a manifest records another cell type
        than its root's role; or the chain does not hold (the stage-DE steps
        or the diagnostic were made from another analysis, the limma step
        checked another DESeq2 summary, the technical table belongs to
        another pseudobulk; on a non-primary basis also an arm made from
        another export, compared with another primary analysis or with
        another conditioning than the registered one).

    Notes
    -----
    E1 compares the control's and the reference's exported RNA stage counts
    within each stage-DE method (DESeq2 ``factor_additive`` union, limma
    ``stage_joint``) and holds only under both; E2 compares the primary
    analyses' CV mean AUROC; E3 compares, along the control's first latent
    component, |Spearman| with stage against the largest |Spearman| with a
    technical variable (``POSITIVE_CONTROL_TECHNICAL_SOURCES``). The overall
    outcome holds when all three hold; it is the control's outcome in
    ``POSITIVE_CONTROL_PRIMARY_COHORT``. Every check of every cohort runs
    before anything is written; each cohort is written to ``<cohort>.tmp/``
    and renamed.
    """
    import shutil

    log = "positive-control"
    types = POSITIVE_CONTROL_CELL_TYPES
    if basis not in POSITIVE_CONTROL_BASES:
        raise ValueError(f"Unknown basis {basis!r}; registered bases (POSITIVE_CONTROL_BASES): {list(POSITIVE_CONTROL_BASES)}")
    entry = POSITIVE_CONTROL_BASES[basis]
    control_root = Path(control_root) if control_root is not None else Path(cell_type_entry(types["control"])["results_root"])
    reference_root = (
        Path(reference_root) if reference_root is not None else Path(cell_type_entry(types["reference"])["results_root"])
    )
    out_dir = Path(out_dir) if out_dir is not None else control_root / CELL_TYPE_SUBDIRS["positive_control"]
    cohort_root = out_dir if entry["arm"] is None else out_dir / basis
    cohorts = list(EXPORT_COHORTS) if cohort == "all" else [cohort]
    for c in cohorts:
        if c not in EXPORT_COHORTS:
            raise ValueError(f"Unknown cohort {c!r}; expected one of {EXPORT_COHORTS} or 'all'")
        if (cohort_root / c).exists() and not force:
            raise FileExistsError(f"{cohort_root / c} exists; rerun with --force to replace it")
    print(f"[{log}] Control ({types['control']}): {control_root}")
    print(f"[{log}] Reference ({types['reference']}): {reference_root} (read, never rerun)")
    print(f"[{log}] Basis: {basis} ({entry['description']})")
    print(f"[{log}] Output directory: {cohort_root}")

    prepared = {}
    for c in cohorts:
        control = _positive_control_read(control_root, c, "control", diagnostics=True, basis=basis)
        reference = _positive_control_read(reference_root, c, "reference", diagnostics=False, basis=basis)
        rows, outcomes = _positive_control_rows(c, control, reference)
        prepared[c] = (control, reference, rows, outcomes)

    results = {}
    for c in cohorts:
        control, reference, rows, outcomes = prepared[c]
        # object dtype: counts stay integers and floats keep their shortest round-trip form in the CSV
        table = pd.DataFrame(rows, columns=list(POSITIVE_CONTROL_COLUMNS), dtype=object)
        findings = (
            [f"control: {m}" for m in control["findings"]] + [f"reference: {m}" for m in reference["findings"]]
            + outcomes.pop("findings")
        )
        manifest = {
            "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "cell_type": types["control"],
            "cell_types": dict(types),
            "cohort": c,
            "primary_cohort": c == POSITIVE_CONTROL_PRIMARY_COHORT,
            "basis": basis,
            "basis_entry": dict(entry),
            "e1_basis_independent": True,
            "expectations": POSITIVE_CONTROL_EXPECTATIONS,
            "rule": POSITIVE_CONTROL_RULE,
            "inputs": {
                "control": {"root": str(control_root), "files": control["files"], "arm": control["arm"]},
                "reference": {"root": str(reference_root), "files": reference["files"], "arm": reference["arm"]},
            },
            "values": {
                "E1": {m: {"control": control["e1"][m], "reference": reference["e1"][m]} for m in control["e1"]},
                "E2": {"control": control["e2"], "reference": reference["e2"]},
                "E3": {
                    "stage_abs_spearman": float(table.loc[table["expectation"] == "E3", "control_value"].iloc[0]),
                    "technical_max": outcomes["E3_technical_max"],
                },
            },
            "outcomes": outcomes,
            "reads": "counts (n_significant), latent.cv_mean_auroc and raw Spearman correlations only; no p-value, "
                     "rejection status or trajectory statistic",
            "outputs": [POSITIVE_CONTROL_TABLE, POSITIVE_CONTROL_MARKDOWN, POSITIVE_CONTROL_MANIFEST],
            "software": _software_versions(),
            "git": _git_revision(),
            "findings": findings,
        }
        tmp = cohort_root / f"{c}.tmp"
        if tmp.exists():
            shutil.rmtree(tmp)
        tmp.mkdir(parents=True)
        table.to_csv(tmp / POSITIVE_CONTROL_TABLE, index=False, float_format=_csv_float)
        (tmp / POSITIVE_CONTROL_MARKDOWN).write_text(render_positive_control(table, manifest))
        (tmp / POSITIVE_CONTROL_MANIFEST).write_text(json.dumps(manifest, indent=2) + "\n")
        final = cohort_root / c
        if final.exists():
            shutil.rmtree(final)
        os.replace(tmp, final)
        print(f"[{log}] Wrote {final}")
        for r in rows:
            print(
                f"[{log}] {c} {r['expectation']}{' ' + r['method'] if r['method'] else ''}: {r['outcome']}"
                + (f" ({r['control']} {r['control_value']} vs {r['reference']} {r['reference_value']})"
                   if r["control_value"] is not None else f" ({r['detail']})")
            )
        for msg in findings:
            print(f"[{log}] FINDING: {c}: {msg}")
        results[c] = {"dir": final, "table": table, "manifest": manifest, "outcomes": outcomes}
    return results


__all__ = [
    "DEFAULT_BUCKET",
    "download_data",
    "subset_adata",
    "subset_and_concat_folder",
    "write_subset_from_folder",
    "prep_anndata",
    "pseudobulk_and_filter",
    "differential_expression",
    "resolve_data_root",
    "MTG_PREFIXES",
    "TAXONOMY",
    "DEFAULT_CELL_TYPE",
    "CELL_TYPES",
    "PREREGISTERED_COHORTS",
    "cell_type_dirs",
    "check_manifest_cell_type",
    "taxonomy_for",
    "select_nuclei",
    "select_astrocyte_nuclei",
    "pseudobulk_donor",
    "build_donor_table",
    "summarize_cohorts",
    "run_pseudobulk",
    "run_covariates",
    "tmm_factors",
    "log_cpm",
    "preprocess_block",
    "run_export",
    "load_analysis_inputs",
    "render_analysis_summary",
    "run_analysis",
    "render_rank_sensitivity",
    "run_rank_sensitivity",
    "PREPROCESSING_VARIANTS",
    "RESIDUALIZATION_COVARIATES",
    "ADOPTION_CHECK_VARIABLES",
    "ADOPTION_RULE",
    "trend_ratio",
    "residualize_values",
    "adoption_check",
    "render_preprocessing_sensitivity",
    "run_preprocessing_sensitivity",
    "TECHNICAL_GENE_SETS",
    "match_gene_sets",
    "technical_donor",
    "run_technical_covariates",
    "DIAGNOSTIC_VARIABLES",
    "render_axis_diagnostics",
    "run_axis_diagnostics",
    "LATENT_SENSITIVITY_ARMS",
    "LATENT_SENSITIVITY_COVARIATES",
    "mediator_check",
    "render_latent_sensitivity",
    "run_latent_sensitivity",
    "STAGE_DE_MODELS",
    "STAGE_DE_COVARIATES",
    "STAGE_DE_CONCORDANCE_RULE",
    "stage_de_formula",
    "stage_de_design",
    "stage_de_fit",
    "render_stage_de",
    "run_stage_de",
    "STAGE_DE_LIMMA_METHOD",
    "STAGE_DE_LIMMA_JOINT_TESTS",
    "STAGE_DE_LIMMA_CONCORDANCE_RULE",
    "stage_de_limma_script",
    "stage_de_limma_versions",
    "write_stage_de_limma_exchange",
    "stage_de_limma_check",
    "stage_de_limma_fit",
    "render_stage_de_limma",
    "run_stage_de_limma",
    "POSITIVE_CONTROL_EXPECTATIONS",
    "POSITIVE_CONTROL_PRIMARY_COHORT",
    "POSITIVE_CONTROL_BASES",
    "render_positive_control",
    "run_positive_control",
]
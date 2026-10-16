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
_SHARED_SUBCLASSES = [
    "L2/3 IT", "L4 IT", "L5 IT", "L6 IT", "L6 IT Car3", "L5 ET", "L5/6 NP", "L6 CT",
    "L6b", "Lamp5", "Lamp5 Lhx6", "Pax6", "Sncg", "Vip", "Sst", "Sst Chodl", "Pvalb",
    "Chandelier", "OPC", "Oligodendrocyte", "Endothelial",
]

# Taxonomy mapping per release, keyed by the date in the file name
# (``..._final-nuclei.<date>.h5ad``). Labels confirmed by the spike in
# ``openspec/changes/mtg-astrocyte-donor-pseudobulks/spike.md``. Every
# ``Subclass`` value of a donor object must be listed (astrocyte or other) and
# every astrocyte ``Supertype`` must be mapped, otherwise the run fails.
TAXONOMY: dict[str, dict] = {
    "2026-06-22": {
        "omic": "rna",
        "astrocyte_subclass": ["Astrocyte"],
        "other_subclasses": _SHARED_SUBCLASSES + ["Immune", "VLMC & Perivascular"],
        "supertypes": {s: s for s in _ASTRO_SUPERTYPES},
    },
    "2024-12-06": {
        "omic": "atac",
        "astrocyte_subclass": ["Astrocyte"],
        "other_subclasses": _SHARED_SUBCLASSES + ["Microglia-PVM", "VLMC"],
        "supertypes": {s: s for s in _ASTRO_SUPERTYPES},
    },
}

# Raw count matrix per omic (h5ad element path)
MATRIX_SOURCE: dict[str, str] = {"rna": "layers/UMIs", "atac": "X"}
RNA_MIN_GENES = 200
ADNC_STAGES: dict[str, str] = {"Not AD": "0", "Low": "0", "Intermediate": "1", "High": "2"}
COHORT_THRESHOLDS: dict[str, int] = {"cohort_30": 30, "cohort_50": 50}
EXPECTED_DONORS: dict[str, int] = {"rna": 84, "atac": 86}
PREREGISTERED_COHORTS: dict[str, dict] = {
    "cohort_30": {"n": 80, "Female": [11, 10, 28], "Male": [9, 10, 12]},
    "cohort_50": {"n": 74},
}
STAGES = ["0", "1", "2"]

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


def taxonomy_for(path: str | Path) -> dict:
    """Return the taxonomy mapping of the release a donor object belongs to."""
    release = release_of(path)
    if release not in TAXONOMY:
        raise ValueError(f"{Path(path).name}: no taxonomy mapping for release {release}")
    return TAXONOMY[release]


def validate_taxonomy(
    path: str | Path, subclasses: Iterable[str], astro_supertypes: Iterable[str]
) -> None:
    """Fail on any ``Subclass`` or astrocyte ``Supertype`` label absent from the mapping.

    Parameters
    ----------
    path:
        Donor object, used for the release and the error message.
    subclasses:
        All ``Subclass`` values of the object.
    astro_supertypes:
        ``Supertype`` values of the astrocyte nuclei.
    """
    tax = taxonomy_for(path)
    known = set(tax["astrocyte_subclass"]) | set(tax["other_subclasses"])
    for label in sorted(set(map(str, subclasses)) - known):
        raise ValueError(
            f"{Path(path).name}: Subclass label {label!r} is not in the taxonomy "
            f"mapping for release {release_of(path)}"
        )
    for label in sorted(set(map(str, astro_supertypes)) - set(tax["supertypes"])):
        raise ValueError(
            f"{Path(path).name}: astrocyte Supertype label {label!r} is not in the "
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


def pseudobulk_donor(
    path: str | Path,
    omic: str,
    *,
    min_genes: int = RNA_MIN_GENES,
    reference_features: np.ndarray | None = None,
) -> dict:
    """Sum raw counts over the astrocyte nuclei of one donor object.

    Parameters
    ----------
    path:
        Donor ``.h5ad`` object (one donor per file).
    omic:
        ``"rna"`` (sums ``layers/UMIs``) or ``"atac"`` (sums ``X``).
    min_genes:
        RNA only: astrocyte nuclei with fewer detected genes are dropped
        before summing. ATAC applies no nucleus QC.
    reference_features:
        Ordered feature names every donor of the omic must match. ``None``
        skips the check (first donor).

    Returns
    -------
    dict
        ``donor_id``, ``omic``, ``release``, ``sex``, ``adnc``, ``excluded``
        (``True`` for ``Reference`` donors), ``n_astro`` (before QC),
        ``n_nuclei`` (summed), ``supertype_counts`` (canonical → count, all
        mapped supertypes), ``row`` (1 × features int64 CSR, ``None`` when
        excluded) and ``var`` (feature annotations).

    Notes
    -----
    Only ``obs`` columns, ``var`` and the astrocyte rows of the count matrix
    are read, so the normalized RNA ``X`` is never loaded.
    """
    import h5py
    import scipy.sparse as sp
    from anndata.io import read_elem

    path = Path(path)
    if omic not in MATRIX_SOURCE:
        raise ValueError(f"Unknown omic {omic!r}")
    tax = taxonomy_for(path)
    if tax["omic"] != omic:
        raise ValueError(f"{path.name}: release {release_of(path)} is {tax['omic']}, not {omic}")
    source = MATRIX_SOURCE[omic]

    with h5py.File(path, "r") as f:
        obs = _read_obs_columns(
            f, ["Donor ID", "Sex", ADNC_COLUMN, "Subclass", "Supertype"], path
        )
        donor = _single_value(obs["Donor ID"], path, "Donor ID")
        sex = _single_value(obs["Sex"], path, "Sex")
        adnc = _single_value(obs[ADNC_COLUMN], path, ADNC_COLUMN)

        subclass = obs["Subclass"].astype(str).to_numpy()
        rows = np.flatnonzero(np.isin(subclass, tax["astrocyte_subclass"]))
        supertype_raw = obs["Supertype"].astype(str).to_numpy()[rows]
        validate_taxonomy(path, subclass, supertype_raw)

        var = read_elem(f["var"])
        if reference_features is not None:
            _check_features(var.index.to_numpy(), reference_features, path)

        record = {
            "donor_id": donor,
            "omic": omic,
            "release": release_of(path),
            "sex": sex,
            "adnc": adnc,
            "excluded": adnc == "Reference",
            "n_astro": int(rows.size),
            "n_nuclei": 0,
            "supertype_counts": {c: 0 for c in tax["supertypes"].values()},
            "row": None,
            "var": var,
        }
        if record["excluded"]:
            return record

        if source not in f:
            raise KeyError(f"{path.name}: count matrix {source!r} not found")
        m = _read_rows(f[source], rows)

    _check_counts(m, path, source)
    m = m.astype(np.int64)
    if omic == "rna":
        n_genes = np.asarray((m > 0).sum(axis=1)).ravel()
        keep = n_genes >= min_genes
        m = m[keep]
        supertype_raw = supertype_raw[keep]

    record["n_nuclei"] = int(m.shape[0])
    record["row"] = sp.csr_matrix(np.asarray(m.sum(axis=0), dtype=np.int64))
    for raw, n in zip(*np.unique(supertype_raw, return_counts=True)):
        record["supertype_counts"][tax["supertypes"][str(raw)]] += int(n)
    del m
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


def compare_preregistered(info: dict) -> list[str]:
    """List deviations of the observed cohorts from the preregistered counts."""
    findings = []
    for cohort, expected in PREREGISTERED_COHORTS.items():
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
    data_root: str | Path | None = None,
    out_dir: str | Path = "results/pseudobulk",
    limit: int | None = None,
    donors: list[str] | None = None,
    force: bool = False,
) -> dict:
    """Build donor-level astrocyte pseudobulks, the donor table and the manifest.

    Parameters
    ----------
    omic:
        ``"rna"``, ``"atac"`` or ``"both"``.
    data_root:
        Data root (resolved with :func:`resolve_data_root`). Donor objects
        are read from ``<root>/<MTG_PREFIXES[omic]>``.
    out_dir:
        Output directory (default ``results/pseudobulk``).
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
    Outputs: ``rna_counts.h5ad``, ``atac_counts.h5ad`` (donors × features,
    int64 CSR), ``donor_table.csv``, ``supertype_composition.csv``,
    ``cohort_summary.csv`` and ``manifest.json``.
    """
    import scipy.sparse as sp

    root = resolve_data_root(data_root)
    out_dir = Path(out_dir)
    omics = ["rna", "atac"] if omic == "both" else [omic]
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
                print(f"{tag} cached ({rec['n_nuclei']} nuclei)")
            else:
                raw = pseudobulk_donor(
                    path, o,
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
                    print(f"{tag} {rec['n_astro']} astrocytes, {rec['n_nuclei']} summed")
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
    findings = compare_preregistered(cohort_info) if not partial else []

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
        "taxonomy": TAXONOMY,
        "qc": {
            "rna": {"min_genes": RNA_MIN_GENES, "applied_to": "astrocyte nuclei before summing"},
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
        "preregistered_cohorts": PREREGISTERED_COHORTS,
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
    data_root: str | Path | None = None,
    out_dir: str | Path = "results/pseudobulk",
) -> dict:
    """Write the donor covariate table and its manifest.

    Parameters
    ----------
    omic:
        ``"both"`` (RNA values, cross-checked against ATAC) or ``"rna"``
        (no cross-check).
    data_root:
        Data root (resolved with :func:`resolve_data_root`).
    out_dir:
        Output directory (default ``results/pseudobulk``).

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
    root = resolve_data_root(data_root)
    out_dir = Path(out_dir)
    omics = ["rna", "atac"] if omic == "both" else ["rna"]
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

    Returns
    -------
    dict
        ``values`` (donors × kept features, log2-CPM, input feature order),
        ``kept`` (names), ``status`` and ``variance`` per input feature,
        ``steps`` (feature counts after each step), ``removed_by_chromosome``,
        ``lib_size_filter``, ``lib_size``, ``tmm_factors``,
        ``mean_log_cpm_range`` of the kept features, ``findings``, and the
        filtered ``filtered_counts`` / ``filtered_log_cpm`` / ``filtered_names``.

    Notes
    -----
    Labels are never used: ``k`` is the only cohort-derived input.
    """
    counts = np.asarray(counts, dtype=float)
    names = np.asarray(names, dtype=str)
    chromosomes = np.asarray(chromosomes, dtype=object)
    n_feat = counts.shape[1]
    excluded = tuple(excluded_chromosomes)
    status = np.full(n_feat, "kept", dtype=object)
    variance = np.full(n_feat, np.nan)
    findings: list[str] = []

    # (2) sex chromosomes
    on_sex = np.isin(chromosomes, excluded)
    status[on_sex] = "sex_chromosome"
    removed_by_chrom = {c: int((chromosomes == c).sum()) for c in excluded}
    idx2 = np.flatnonzero(~on_sex)

    # (3) abundance: CPM over the post-chromosome library sizes
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

    # (6) top-N by population variance, ties by name; output keeps input order
    var3 = lcpm.var(axis=0, ddof=0)
    variance[idx3] = var3
    order = np.lexsort((names[idx3], -var3))
    if idx3.size < top_n:
        findings.append(f"only {idx3.size} features pass the abundance filter (< N = {top_n}); all are kept")
    sel = np.sort(order[:top_n])
    dropped = np.setdiff1d(np.arange(idx3.size), sel)
    status[idx3[dropped]] = "variance"
    values = lcpm[:, sel]
    if not np.all(np.isfinite(values)):
        raise ValueError("log2-CPM contains non-finite values")

    mean_lcpm = values.mean(axis=0) if values.size else np.array([np.nan])
    return {
        "values": values,
        "kept": names[idx3[sel]],
        "status": status,
        "variance": variance,
        "steps": {
            "input": int(n_feat),
            "after_sex_chromosome": int(idx2.size),
            "after_abundance": int(idx3.size),
            "after_variance": int(sel.size),
        },
        "removed_by_chromosome": removed_by_chrom,
        "lib_size_filter": lib_filter,
        "lib_size": lib,
        "tmm_factors": factors,
        "mean_log_cpm_range": [float(np.min(mean_lcpm)), float(np.max(mean_lcpm))],
        "findings": findings,
        "filtered_counts": c3,
        "filtered_log_cpm": lcpm,
        "filtered_names": names[idx3],
    }


_EXPORT_INPUTS = {
    "rna_counts.h5ad": "pseudobulk",
    "atac_counts.h5ad": "pseudobulk",
    "donor_table.csv": "pseudobulk",
    "manifest.json": "pseudobulk",
    "donor_covariates.csv": "covariates",
}


def check_export_inputs(pseudobulk_dir: str | Path, cohorts: Iterable[str]) -> dict:
    """Validate the pseudobulk directory before exporting; return its manifest.

    Raises
    ------
    FileNotFoundError
        A required file is missing (naming the command that produces it).
    ValueError
        The pseudobulk run is partial, or a cohort is unknown or not analysable.
    """
    pb = Path(pseudobulk_dir)
    for name, command in _EXPORT_INPUTS.items():
        if not (pb / name).exists():
            raise FileNotFoundError(
                f"{pb / name} is missing; run `sead-sexdiff {command}` (and sync) to produce it"
            )
    manifest = json.loads((pb / "manifest.json").read_text())
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
    return df


def run_export(
    *,
    pseudobulk_dir: str | Path = "results/pseudobulk",
    out_dir: str | Path = "results/motco_input",
    cohort: str = "all",
    data_root: str | Path | None = None,
    force: bool = False,
    debug_dir: str | Path | None = None,
) -> dict:
    """Export MOTCO-ready RNA and ATAC blocks per cohort.

    Parameters
    ----------
    pseudobulk_dir:
        Output directory of ``pseudobulk`` and ``covariates``.
    out_dir:
        Export root; each cohort goes to ``<out_dir>/<cohort>/``.
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
        Cohort → ``{"dir": Path, "steps": {omic: counts}, "findings": [...]}``.

    Notes
    -----
    Per cohort and omic: restrict to the cohort → drop chrX/Y/M → CPM ≥ 1 in
    ≥ k donors → TMM → log2-CPM (prior 2) → top N by variance. Covariates are
    exported alongside, not applied. Each cohort is written to
    ``<cohort>.tmp/`` and renamed when complete.
    """
    import shutil
    import scipy.sparse as sp

    pb = Path(pseudobulk_dir)
    out_dir = Path(out_dir)
    cohorts = list(EXPORT_COHORTS) if cohort == "all" else [cohort]
    for c in cohorts:
        if c not in EXPORT_COHORTS:
            raise ValueError(f"Unknown cohort {c!r}; expected one of {EXPORT_COHORTS} or 'all'")
    pb_manifest = check_export_inputs(pb, cohorts)
    for c in cohorts:
        if (out_dir / c).exists() and not force:
            raise FileExistsError(f"{out_dir / c} exists; rerun with --force to replace it")
    print(f"[export] Pseudobulk directory: {pb}")
    print(f"[export] Output directory: {out_dir}")

    mapping, annotation = resolve_gene_annotation(data_root)
    print(f"[export] Annotation: {annotation['release']} ({annotation['n_genes']} genes, sha256 {annotation['sha256'][:12]}…)")

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
    n_unmapped = int((blocks["rna"]["chrom"] == UNMAPPED_CHROMOSOME).sum())
    annotation["n_rna_genes"] = int(blocks["rna"]["names"].size)
    annotation["n_rna_unmapped"] = n_unmapped
    global_findings = []
    if n_unmapped:
        global_findings.append(
            f"rna: {n_unmapped} of {blocks['rna']['names'].size} genes are absent from "
            f"{annotation['release']}; they are kept (chromosome '{UNMAPPED_CHROMOSOME}')"
        )
    pb_manifest_sha = _sha256_file(pb / "manifest.json")
    cov_manifest = pb / "covariates_manifest.json"

    results = {}
    for c in cohorts:
        sub = table[table[c].astype(bool)].sort_values("donor_id").reset_index(drop=True)
        donors = sub["donor_id"].to_numpy()
        missing_cov = sorted(set(donors) - set(covariates.donor_id))
        if missing_cov:
            raise ValueError(f"{pb / 'donor_covariates.csv'}: no row for cohort {c} donor(s) {missing_cov}")
        k, k_cell, cells = smallest_cell(sub)
        print(f"[export] {c}: {len(donors)} donors, k = {k} (smallest cell {k_cell})")
        findings = list(global_findings)
        res = {}
        for o, blk in blocks.items():
            pos = pd.Index(blk["donors"]).get_indexer(donors)
            if (pos < 0).any():
                raise ValueError(f"{o}_counts.h5ad lacks cohort {c} donor(s) {list(donors[pos < 0])}")
            counts = blk["X"][pos].toarray()
            res[o] = preprocess_block(counts, blk["names"], blk["chrom"], k=k)
            findings += [f"{o}: {msg}" for msg in res[o]["findings"]]
            s = res[o]["steps"]
            print(f"[export] {c} {o}: {s['input']} → {s['after_sex_chromosome']} → "
                  f"{s['after_abundance']} → {s['after_variance']}")
            if debug_dir is not None:
                _write_debug(Path(debug_dir), c, o, donors, res[o])

        tmp = out_dir / f"{c}.tmp"
        if tmp.exists():
            shutil.rmtree(tmp)
        tmp.mkdir(parents=True)
        for o in ("rna", "atac"):
            df = pd.DataFrame(res[o]["values"], index=pd.Index(donors, name="donor_id"), columns=res[o]["kept"])
            df.to_csv(tmp / f"{o}.csv", float_format=EXPORT_FLOAT_FORMAT)
            _features_frame(blocks[o]["names"], blocks[o]["chrom"], res[o], blocks[o]["gene_ids"]).to_csv(
                tmp / f"features_{o}.csv", index=False, float_format=EXPORT_FLOAT_FORMAT
            )
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
        manifest = {
            "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "cohort": c,
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
            "preprocessing_order": [
                "restrict to cohort donors",
                "remove features on excluded chromosomes",
                "keep features with CPM >= min_cpm in >= k donors (raw library sizes after chromosome removal)",
                "TMM factors on the filtered counts (library size = filtered column sums)",
                "log2-CPM, prior count scaled by effective library size (edgeR cpm(log = TRUE))",
                "keep the N features with the highest population variance (ddof 0), ties by name; input order kept",
            ],
            "thresholds": {
                "excluded_chromosomes": list(EXPORT_EXCLUDED_CHROMOSOMES),
                "min_cpm": EXPORT_MIN_CPM,
                "k": k,
                "k_source_cell": k_cell,
                "prior_count": EXPORT_PRIOR_COUNT,
                "N": EXPORT_TOP_N,
                "tmm": {"logratio_trim": TMM_LOGRATIO_TRIM, "sum_trim": TMM_SUM_TRIM,
                        "weighting": True, "reference": "upper-quartile factor closest to the mean"},
                "float_format": EXPORT_FLOAT_FORMAT,
            },
            "feature_counts": {o: res[o]["steps"] for o in ("rna", "atac")},
            "removed_by_chromosome": {o: res[o]["removed_by_chromosome"] for o in ("rna", "atac")},
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
            "input": {
                "pseudobulk_dir": str(pb.resolve()),
                "pseudobulk_manifest_sha256": pb_manifest_sha,
                "pseudobulk_manifest_created_at": pb_manifest.get("created_at"),
                "covariates_manifest_sha256": _sha256_file(cov_manifest) if cov_manifest.exists() else None,
            },
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
        final = out_dir / c
        if final.exists():
            shutil.rmtree(final)
        os.replace(tmp, final)
        print(f"[export] Wrote {final}")
        for msg in findings:
            print(f"[export] FINDING: {c}: {msg}")
        results[c] = {"dir": final, "steps": manifest["feature_counts"], "findings": findings}

    print("[export] Features after each step:")
    rows = [
        {"cohort": c, "omic": o, **steps}
        for c, r in results.items() for o, steps in r["steps"].items()
    ]
    print(pd.DataFrame(rows).to_string(index=False))
    return results


def _write_debug(debug_dir: Path, cohort: str, omic: str, donors: np.ndarray, res: dict) -> None:
    """Write filtered counts, TMM factors and log2-CPM for the edgeR cross-check."""
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
    "pseudobulk_donor",
    "build_donor_table",
    "summarize_cohorts",
    "run_pseudobulk",
    "run_covariates",
    "tmm_factors",
    "log_cpm",
    "preprocess_block",
    "run_export",
]
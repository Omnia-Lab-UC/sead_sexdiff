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
]
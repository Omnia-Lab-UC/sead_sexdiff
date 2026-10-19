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


def select_astrocyte_nuclei(
    f,
    path: str | Path,
    omic: str,
    *,
    min_genes: int = RNA_MIN_GENES,
    reference_features: np.ndarray | None = None,
    read_counts: bool = True,
) -> dict:
    """Select the astrocyte nuclei of one open donor object and read their raw counts.

    The one nucleus-selection code path of the MTG steps: ``pseudobulk_donor``
    sums what it returns and ``technical_donor`` summarizes the same nuclei,
    so "the same nuclei as the pseudobulk" is a property of the code, not a
    claim.

    Parameters
    ----------
    f:
        The donor object, opened with ``h5py.File(path, "r")``.
    path:
        Its path (release, taxonomy mapping and error messages).
    omic:
        ``"rna"`` (counts in ``layers/UMIs``) or ``"atac"`` (counts in ``X``).
    min_genes:
        RNA only: astrocyte nuclei with fewer detected genes are dropped.
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
        ``donor_id``, ``omic``, ``release``, ``sex``, ``adnc``, ``excluded``
        (``True`` for ``Reference`` donors), ``n_astro`` (astrocyte nuclei
        before QC), ``var`` (feature annotations), ``rows`` (``obs`` positions
        of the selected nuclei, after QC), ``supertypes`` (their raw
        ``Supertype`` labels) and ``counts`` (their raw counts, nuclei ×
        features int64 CSR; ``None`` when not read). For an excluded donor
        ``rows`` and ``supertypes`` are empty and ``counts`` is ``None``.

    Notes
    -----
    Only ``obs`` columns, ``var`` and the selected rows of the count matrix
    are read, so the normalized RNA ``X`` is never loaded.
    """
    from anndata.io import read_elem

    path = Path(path)
    if omic not in MATRIX_SOURCE:
        raise ValueError(f"Unknown omic {omic!r}")
    tax = taxonomy_for(path)
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
    rows = np.flatnonzero(np.isin(subclass, tax["astrocyte_subclass"]))
    supertype_raw = obs["Supertype"].astype(str).to_numpy()[rows]
    validate_taxonomy(path, subclass, supertype_raw)

    var = read_elem(f["var"])
    if reference_features is not None:
        _check_features(var.index.to_numpy(), reference_features, path)

    selection = {
        "donor_id": donor,
        "omic": omic,
        "release": release_of(path),
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
    The nuclei are selected by :func:`select_astrocyte_nuclei`, the code path
    shared with the technical-covariates step; this function only sums them.
    Only ``obs`` columns, ``var`` and the astrocyte rows of the count matrix
    are read, so the normalized RNA ``X`` is never loaded.
    """
    import h5py
    import scipy.sparse as sp

    path = Path(path)
    if omic not in MATRIX_SOURCE:
        raise ValueError(f"Unknown omic {omic!r}")
    tax = taxonomy_for(path)
    if tax["omic"] != omic:
        raise ValueError(f"{path.name}: release {release_of(path)} is {tax['omic']}, not {omic}")

    with h5py.File(path, "r") as f:
        sel = select_astrocyte_nuclei(
            f, path, omic, min_genes=min_genes, reference_features=reference_features
        )

    record = {
        "donor_id": sel["donor_id"],
        "omic": omic,
        "release": sel["release"],
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
    gene_sets: dict[str, list[str]] | None = None,
    min_genes: int = RNA_MIN_GENES,
    reference_features: np.ndarray | None = None,
) -> dict:
    """Technical quality summaries of the astrocyte nuclei of one donor object.

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
    min_genes, reference_features:
        Passed to :func:`select_astrocyte_nuclei`, so the nuclei are exactly
        the pseudobulk's.

    Returns
    -------
    dict
        JSON-serializable. Header (``donor_id``, ``omic``, ``release``,
        ``sex``, ``adnc``, ``excluded``, ``n_astro``, ``n_nuclei``,
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
    spanning several libraries become astrocyte-weighted means. ATAC reads
    ``obs`` only; nothing reads the normalized RNA ``X``.
    """
    import h5py

    path = Path(path)
    if omic not in TECHNICAL_OBS_COLUMNS:
        raise ValueError(f"Unknown omic {omic!r}")
    columns = list(TECHNICAL_OBS_COLUMNS[omic])
    with h5py.File(path, "r") as f:
        sel = select_astrocyte_nuclei(
            f, path, omic, min_genes=min_genes, reference_features=reference_features,
            read_counts=omic == "rna",
        )
        qc = _read_obs_columns(f, columns, path)
    tax = taxonomy_for(path)
    record = {
        "donor_id": sel["donor_id"],
        "omic": omic,
        "release": sel["release"],
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
    donor: str, rec: dict, pb_fractions: dict, table: pd.DataFrame, pb_dir: Path
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
        n_cache = int(_load_cache(cache)["n_nuclei"])
        out["pseudobulk_cache_n_nuclei"] = n_cache
        if n_cache != rec["n_nuclei"]:
            raise ValueError(
                f"{donor}: {rec['n_nuclei']} astrocyte nuclei selected but the pseudobulk cache {cache} "
                f"records n_nuclei = {n_cache}"
            )
    hit = table.loc[table["donor_id"] == donor, "n_rna"]
    out["donor_table_n_rna"] = int(hit.iloc[0]) if len(hit) else None
    if len(hit) and int(hit.iloc[0]) != rec["n_nuclei"]:
        raise ValueError(
            f"{donor}: {rec['n_nuclei']} astrocyte nuclei selected but {pb_dir / 'donor_table.csv'} "
            f"records n_rna = {int(hit.iloc[0])}"
        )
    return out


def run_technical_covariates(
    omic: str = "both",
    *,
    data_root: str | Path | None = None,
    pseudobulk_dir: str | Path = "results/pseudobulk",
    out_dir: str | Path = "results/pseudobulk",
    force: bool = False,
) -> dict:
    """Write the donor technical table of the astrocyte nuclei and its manifest.

    Parameters
    ----------
    omic:
        ``"both"`` (RNA objects plus the ATAC library metrics) or ``"rna"``.
    data_root:
        Data root (resolved with :func:`resolve_data_root`) holding the MTG
        donor objects.
    pseudobulk_dir:
        Output directory of ``pseudobulk`` (``rna_counts.h5ad``,
        ``donor_table.csv``, ``manifest.json`` are read; the per-donor cache
        is read when present).
    out_dir:
        Output directory (default ``results/pseudobulk``).
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
        made with different gene sets, QC columns or ``min_genes``, or a
        donor object's ``Donor ID`` differs from its file name.

    Notes
    -----
    Per RNA donor object the astrocyte nuclei are selected by
    :func:`select_astrocyte_nuclei` (the pseudobulk's code path), summarized
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
    root = resolve_data_root(data_root)
    pb = Path(pseudobulk_dir)
    out_dir = Path(out_dir)
    omics = ["rna", "atac"] if omic == "both" else ["rna"]
    for name in ("rna_counts.h5ad", "donor_table.csv", "manifest.json"):
        if not (pb / name).exists():
            raise FileNotFoundError(
                f"{pb / name} is missing; run `sead-sexdiff pseudobulk` (and sync) to produce it"
            )
    cov_files = {name: pb / name for name in ("donor_covariates.csv", "covariates_manifest.json")}
    cov_before = {name: _sha256_file(p) if p.exists() else None for name, p in cov_files.items()}
    print(f"[{log}] Data root: {root.resolve()}")
    print(f"[{log}] Pseudobulk directory: {pb}")
    print(f"[{log}] Output directory: {out_dir}")

    pb_manifest = json.loads((pb / "manifest.json").read_text())
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
                n_cached[o] += 1
                print(f"{tag} cached ({rec['n_nuclei']} nuclei)")
            else:
                rec = technical_donor(
                    path, o,
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
                    print(f"{tag} {rec['n_astro']} astrocytes, {rec['n_nuclei']} nuclei, {rec['n_libraries']} libraries; {frac}")
                else:
                    print(f"{tag} {rec['n_nuclei']} nuclei")
            records[o][donor] = rec
            if o == "rna":
                cross[donor] = _cross_check_technical(donor, rec, pb_fractions, table, pb)
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

    findings: list[str] = []
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
            "function": "select_astrocyte_nuclei (shared with pseudobulk_donor)",
            "rna": {"matrix": MATRIX_SOURCE["rna"], "min_genes": RNA_MIN_GENES},
            "atac": {"matrix": None, "note": "obs columns only; the ATAC count matrix is not read"},
            "reference_excluded": True,
            "taxonomy": TAXONOMY,
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
            "<qc>_mean, <qc>_sd": "nucleus-weighted mean and sample SD of the obs column over the selected nuclei (library-level columns become astrocyte-weighted means)",
            "n_libraries_rna": f"distinct values of obs['{TECHNICAL_LIBRARY_COLUMN}'] among the selected nuclei",
            "atac <qc>_mean": "nucleus-weighted mean over the selected ATAC astrocyte nuclei; 'constant' per donor and column in 'atac_constancy'",
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
    variant: dict | None = None,
    debug_dir: Path | None = None,
    log_prefix: str = "export",
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
    variant:
        For a preprocessing-sensitivity variant, ``{"name", "description",
        "settings", "reference_export_dir"}``, written as the manifest's
        ``variant`` block; ``None`` for the primary export, whose manifest then
        has no such key.
    debug_dir:
        If given, the filtered counts, TMM factors and log2-CPM per omic go
        there (``_write_debug``).
    log_prefix:
        Prefix of the progress lines, ``[export]`` by default.

    Returns
    -------
    dict
        ``{"dir": Path, "manifest": dict, "findings": [...], "k": int, "cells": dict}``.
    """
    import shutil

    excluded = tuple(excluded_chromosomes)
    donors = sub["donor_id"].to_numpy()
    k, k_cell, cells = smallest_cell(sub)
    print(f"[{log_prefix}] {cohort}: {len(donors)} donors, k = {k} (smallest cell {k_cell})")
    findings = []
    n_unmapped = annotation.get("n_rna_unmapped", 0)
    if n_unmapped:
        findings.append(
            f"rna: {n_unmapped} of {blocks['rna']['names'].size} genes are absent from "
            f"{annotation['release']}; they are kept (chromosome '{UNMAPPED_CHROMOSOME}')"
        )
    res = {}
    for o, blk in blocks.items():
        pos = pd.Index(blk["donors"]).get_indexer(donors)
        if (pos < 0).any():
            raise ValueError(f"{o}_counts.h5ad lacks cohort {cohort} donor(s) {list(donors[pos < 0])}")
        counts = blk["X"][pos].toarray()
        res[o] = preprocess_block(counts, blk["names"], blk["chrom"], k=k, excluded_chromosomes=excluded)
        findings += [f"{o}: {msg}" for msg in res[o]["findings"]]
        s = res[o]["steps"]
        print(f"[{log_prefix}] {cohort} {o}: {s['input']} → {s['after_sex_chromosome']} → "
              f"{s['after_abundance']} → {s['after_variance']}")
        if debug_dir is not None:
            _write_debug(Path(debug_dir), cohort, o, donors, res[o])

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
        "preprocessing_order": [
            "restrict to cohort donors",
            "remove features on excluded chromosomes",
            "keep features with CPM >= min_cpm in >= k donors (raw library sizes after chromosome removal)",
            "TMM factors on the filtered counts (library size = filtered column sums)",
            "log2-CPM, prior count scaled by effective library size (edgeR cpm(log = TRUE))",
            "keep the N features with the highest population variance (ddof 0), ties by name; input order kept",
        ],
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
    pb = Path(pseudobulk_dir)
    out_dir = Path(out_dir)
    cohorts = list(EXPORT_COHORTS) if cohort == "all" else [cohort]
    for c in cohorts:
        if c not in EXPORT_COHORTS:
            raise ValueError(f"Unknown cohort {c!r}; expected one of {EXPORT_COHORTS} or 'all'")
    check_export_inputs(pb, cohorts)
    for c in cohorts:
        if (out_dir / c).exists() and not force:
            raise FileExistsError(f"{out_dir / c} exists; rerun with --force to replace it")
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


def check_motco_version() -> str:
    """Return the installed MOTCO version; raise if it is older than the minimum.

    Raises
    ------
    RuntimeError
        The installed MOTCO is older than ``ANALYSIS_MOTCO_MIN_VERSION``.
    """
    import motco

    installed = str(motco.__version__)
    required = ANALYSIS_MOTCO_MIN_VERSION
    try:
        from packaging.version import Version

        too_old = Version(installed) < Version(required)
    except ImportError:
        def key(v: str) -> tuple[int, ...]:
            return tuple(int(p) for p in re.findall(r"\d+", v)[:3])

        too_old = key(installed) < key(required)
    if too_old:
        raise RuntimeError(
            f"MOTCO {installed} is installed but the analysis needs MOTCO >= {required} "
            "(covariate-adjusted analyze_trajectories); run `uv sync --locked`"
        )
    return installed


def load_analysis_inputs(cohort_dir: str | Path) -> dict:
    """Read one cohort's export, cross-check it and build the model's metadata table.

    Parameters
    ----------
    cohort_dir:
        ``<export root>/<cohort>/`` as written by ``run_export``.

    Returns
    -------
    dict
        ``blocks`` (omic → donors × features DataFrame, read as ``motco analyze``
        reads ``--block`` files), ``metadata`` (design columns plus
        ``ANALYSIS_COVARIATES``, one row per donor sorted by ``donor_id``),
        ``export_manifest`` and ``inputs`` (file → path and sha256).

    Raises
    ------
    FileNotFoundError
        The cohort directory or one of its five files is missing.
    ValueError
        The tables list different donors or disagree with the export manifest.

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
        skipped; the output for CV records is unchanged.
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
    return "\n".join(lines) + "\n"


def _cli_equivalent(cohort_dir: Path, metadata_path: Path, params, n_jobs: int) -> dict:
    """The ``motco analyze`` command that reproduces one analysis directory (CV or forced rank)."""
    import shlex

    argv = ["motco", "analyze"]
    for o in ANALYSIS_BLOCKS:
        argv += ["--block", f"{o}={cohort_dir / f'{o}.csv'}"]
    argv += [
        "--metadata", str(metadata_path),
        "--sample-col", "donor_id", "--group-col", "sex", "--stage-col", "stage",
        "--stage-order", ",".join(STAGES),
    ]
    for c in ANALYSIS_COVARIATES:
        argv += ["--covariate", c]
    argv += [
        "--missing-covariates", ANALYSIS_MISSING_COVARIATES,
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
) -> dict:
    """Run MOTCO on one loaded cohort and write a complete analysis directory.

    Shared by ``run_analysis`` (CV-selected rank) and ``run_rank_sensitivity``
    (forced rank): the two differ only in ``params`` and in where the directory
    goes, so a forced-rank directory is a CV directory with one setting changed.

    Parameters
    ----------
    cohort:
        Cohort name (title and manifest).
    cohort_dir:
        The export directory the inputs were read from.
    loaded:
        ``load_analysis_inputs(cohort_dir)``.
    params:
        ``motco.analysis.TrajectoryAnalysisParams``; ``forced_components`` set
        for a forced-rank run, ``None`` for CV rank selection.
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
    blocks = loaded["blocks"]
    export_manifest = loaded["export_manifest"]
    final = out_dir / final_name
    published = Path(published_dir) if published_dir is not None else final
    # Parse the table back as `motco analyze --metadata` does, so the
    # CLI-equivalent run sees exactly what this one sees
    metadata_csv = loaded["metadata"].to_csv(index=False)
    metadata = pd.read_csv(io.StringIO(metadata_csv))
    metadata["donor_id"] = metadata["donor_id"].astype(str)
    forced = params.forced_components
    print(
        f"[{log_prefix}] {cohort}: {len(metadata)} donors, "
        + " + ".join(f"{b.shape[1]} {o}" for o, b in blocks.items())
        + f" features, covariates {', '.join(ANALYSIS_COVARIATES)}"
    )

    dataset = TrajectoryDataset(
        blocks=[OmicsBlock(o, blocks[o], "identity") for o in ANALYSIS_BLOCKS],
        metadata=metadata,
        sample_col="donor_id",
        group_col="sex",
        stage_col="stage",
        stage_order=STAGES,
        covariate_cols=ANALYSIS_COVARIATES,
        missing_covariates=ANALYSIS_MISSING_COVARIATES,
    )
    selection = "CV rank selection" if forced is None else f"forced rank {forced}"
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
        "cohort": cohort,
        "n_donors": int(len(metadata)),
        "inputs": {
            "dir": str(cohort_dir),
            "files": loaded["inputs"],
            "export_manifest_sha256": loaded["inputs"]["manifest.json"]["sha256"],
            "export_manifest_created_at": export_manifest.get("created_at"),
            "export_git": export_manifest.get("git"),
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
            "feature_counts": export_manifest.get("feature_counts"),
            "n_age_capped": n_capped,
            "export_findings": export_manifest.get("findings", []),
        },
        "covariates": {
            "entered": list(ANALYSIS_COVARIATES),
            "derived": {name: f"log10({src})" for name, src in ANALYSIS_LOG10_COLUMNS.items()},
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
            "permutations": params.permutations,
            "seed": params.seed,
            "alpha": ANALYSIS_ALPHA,
            "rejection_rule": ANALYSIS_REJECTION_RULE,
            "n_jobs": n_jobs,
            "n_jobs_note": "RRPP draws depend on n_jobs (one seed stream per worker); rerun with the same n_jobs",
            "attribution": False,
        },
        "cli_equivalent": _cli_equivalent(cohort_dir, published / "analysis_metadata.csv", params, n_jobs),
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
    input_dir: str | Path = "results/motco_input",
    out_dir: str | Path = "results/motco_analysis",
    n_jobs: int = 1,
    force: bool = False,
) -> dict:
    """Run the covariate-adjusted MOTCO trajectory analysis per cohort.

    Parameters
    ----------
    cohort:
        ``"cohort_30"``, ``"cohort_50"`` or ``"all"``.
    input_dir:
        Export root; each cohort is read from ``<input_dir>/<cohort>/``.
    out_dir:
        Analysis root; each cohort goes to ``<out_dir>/<cohort>/``.
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

    input_dir = Path(input_dir)
    out_dir = Path(out_dir)
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
    print(f"[analyze] Input directory: {input_dir}")
    print(f"[analyze] Output directory: {out_dir}")

    results = {}
    for c in cohorts:
        cohort_dir = input_dir / c
        loaded = load_analysis_inputs(cohort_dir)
        params = TrajectoryAnalysisParams(permutations=ANALYSIS_PERMUTATIONS, seed=ANALYSIS_SEED, n_jobs=n_jobs)
        results[c] = _analyze_cohort(c, cohort_dir, loaded, params, out_dir=out_dir, final_name=c, n_jobs=n_jobs)
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
}
PREPROCESSING_PRIMARY_LABEL = "primary"
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


def _variant_row(record: dict, *, variant: str, cohort: str, directory: str) -> dict:
    """``_sensitivity_row`` with a leading ``variant`` column (``primary`` or the variant name)."""
    row = _sensitivity_row(
        record, cohort=cohort, component_selection=record["latent"]["component_selection"], directory=directory
    )
    return {"variant": variant, **row}


def _feature_accounting(primary_export_dir: str | Path, variant_export_dir: str | Path) -> dict:
    """Compare the feature filtering of the primary and a variant export.

    Parameters
    ----------
    primary_export_dir, variant_export_dir:
        Export directories of one cohort (``features_<omic>.csv`` and
        ``manifest.json`` are read).

    Returns
    -------
    dict
        Omic → ``feature_counts`` and ``removed_by_chromosome`` of both exports,
        ``chrx`` (input count, the variant's status counts, kept count and
        kept feature names), ``kept`` (both kept counts; shared, primary-only
        and variant-only) and ``named`` (status of
        ``PREPROCESSING_NAMED_FEATURES`` in both exports, when present).

    Raises
    ------
    ValueError
        The two exports list different input features, or a chrX feature of
        the variant carries the ``sex_chromosome`` status.
    """
    p, v = Path(primary_export_dir), Path(variant_export_dir)
    pm = json.loads((p / "manifest.json").read_text())
    vm = json.loads((v / "manifest.json").read_text())
    dtypes = {"feature": str, "chromosome": str, "status": str}
    out = {}
    for o in ANALYSIS_BLOCKS:
        pf = pd.read_csv(p / f"features_{o}.csv", dtype=dtypes)
        vf = pd.read_csv(v / f"features_{o}.csv", dtype=dtypes)
        if pf["feature"].tolist() != vf["feature"].tolist():
            raise ValueError(f"{p / f'features_{o}.csv'} and {v / f'features_{o}.csv'} list different input features")
        p_kept = set(pf.loc[pf["status"] == "kept", "feature"])
        v_kept = set(vf.loc[vf["status"] == "kept", "feature"])
        v_chrx = vf[vf["chromosome"] == "chrX"]
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
        out[o] = {
            "feature_counts": {"primary": pm["feature_counts"][o], "variant": vm["feature_counts"][o]},
            "removed_by_chromosome": {
                "primary": pm.get("removed_by_chromosome", {}).get(o),
                "variant": vm.get("removed_by_chromosome", {}).get(o),
            },
            "chrx": {
                "input": int(len(v_chrx)),
                "status_counts_variant": {s: int(n) for s, n in sorted(v_chrx["status"].value_counts().items())},
                "kept": len(kept_chrx),
                "kept_features": kept_chrx,
            },
            "kept": {
                "primary": len(p_kept),
                "variant": len(v_kept),
                "shared": len(p_kept & v_kept),
                "primary_only": len(p_kept - v_kept),
                "variant_only": len(v_kept - p_kept),
            },
            "named": named,
        }
    return out


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
) -> str:
    """Render ``sensitivity.md`` from the rows of ``sensitivity.csv`` and the feature accounting.

    Parameters
    ----------
    rows:
        ``_variant_row`` dicts, the primary row first and the variant row second.
    features:
        ``_feature_accounting`` output.
    variant, settings, description:
        The variant's name, the export settings it changes and its description
        (from ``PREPROCESSING_VARIANTS``).
    cohort:
        Cohort name for the title.
    seed, permutations, alpha:
        Named in the header; a statistic is rejected when ``p <= alpha``.
    findings:
        The comparison's findings, listed at the end.

    Returns
    -------
    str
        Markdown: a header naming the variant, its settings, seed,
        permutations, α and both ranks and the reading rules; one table per
        statistic (two rows); the feature table per omic and export; a run
        table (rank selection, CV AUROC, runtime, directory); the findings.
        Every number comes from the rows and the accounting; the runtimes are
        the one thing here that is not in ``sensitivity.csv``.
    """
    primary, var = rows[0], rows[1]
    changed = ", ".join(f"`{key}` = {list(val) if isinstance(val, (list, tuple)) else val}" for key, val in settings.items())
    lines = [
        f"# MOTCO preprocessing sensitivity: {cohort}, variant `{variant}`",
        "",
        f"{description}. Export setting(s) changed: {changed} (primary: `excluded_chromosomes` = "
        f"{list(EXPORT_EXCLUDED_CHROMOSOMES)}); every other export setting and every analysis setting equal "
        f"the primary's. The `{primary['variant']}` row is the CV-selected primary analysis at rank "
        f"{primary['rank']}, which stays the primary result; the `{variant}` row (CV-selected rank "
        f"{var['rank']}) is read as sensitivity to it, never as its replacement, and no variant is chosen from "
        f"the results. {permutations} RRPP permutations with seed {seed} in both rows (the same permutation "
        f"indices); α = {alpha} (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in "
        "both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that "
        "differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of "
        "the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse.",
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
    lines += [
        "",
        "## Features",
        "",
        "Counts after each export step per omic; `kept on chrX` and `kept only here` are per export, "
        "`kept shared` is the number of kept features common to both exports.",
        "",
        "| omic | export | input | after chromosome removal | after abundance | kept (after variance) | "
        "kept on chrX | kept shared | kept only here |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for o, acc in features.items():
        for label, key, chrx, only in (
            (primary["variant"], "primary", 0, acc["kept"]["primary_only"]),
            (variant, "variant", acc["chrx"]["kept"], acc["kept"]["variant_only"]),
        ):
            fc = acc["feature_counts"][key]
            lines.append(
                f"| {o} | {label} | {fc['input']} | {fc['after_sex_chromosome']} | {fc['after_abundance']} | "
                f"{fc['after_variance']} | {chrx} | {acc['kept']['shared']} | {only} |"
            )
    named = [(o, name, info) for o, acc in features.items() for name, info in acc["named"].items()]
    if named:
        lines.append("")
        for o, name, info in named:
            lines.append(
                f"- `{name}` ({info['chromosome']}, {o}): status `{info['primary']}` in the "
                f"{primary['variant']} export, `{info['variant']}` in the {variant} export."
            )
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
        ``motco_input/<cohort>/`` (the variant export), ``motco_analysis/<cohort>/``
        (the variant analysis) and ``comparison/<cohort>/`` (``sensitivity.csv``,
        ``sensitivity.md``, ``sensitivity_manifest.json``).
    data_root:
        Data root holding (or receiving) the gene annotation.
    n_jobs:
        Parallel jobs for MOTCO's CV and RRPP (``-1`` = all CPUs); as in
        ``run_analysis``, values other than 1 change the permutation draws and
        are recorded as a finding. The reported results use 1.
    force:
        Replace a cohort's three directories under the variant.

    Returns
    -------
    dict
        Variant → cohort → ``{"export", "analysis", "comparison", "rows", "table",
        "manifest", "findings"}``.

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
    hashes, MOTCO version, seed, permutations and CV selection) runs for every
    requested variant and cohort before anything is written. Per cohort the
    variant export is one ``_export_cohort`` call with the variant's keyword(s)
    (the primary pipeline with one setting changed: library sizes, the
    abundance filter, TMM, log2-CPM and the variance ranking recomputed), its
    ``metadata.csv`` and ``covariates.csv`` must be byte-identical to the
    primary's, and the variant analysis is one ``_analyze_cohort`` call with
    the primary's parameters (CV rank selection, ``ANALYSIS_PERMUTATIONS``,
    ``ANALYSIS_SEED``, so RRPP draws the same permutation indices as the
    primary). The comparison holds the primary row first and the variant row
    second; a rejection status that differs is preprocessing dependence of the
    sensitivity analysis and the primary analysis stays primary.
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
    for v in variants:
        for c in cohorts:
            for sub in subdirs:
                if (out_dir / v / sub / c).exists() and not force:
                    raise FileExistsError(
                        f"{out_dir / v / sub / c} exists; rerun with --force to replace the variant's export, "
                        f"analysis and comparison of {c}"
                    )
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
    # Check every cohort's primary export and analysis before writing anything
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
        prepared[c] = {"sub": sub, "k": k, "cells": cells, "export": ref_export, "analysis": reference}

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
        settings = {"excluded_chromosomes": EXPORT_EXCLUDED_CHROMOSOMES, **spec["export"]}
        recorded_settings = {
            key: list(val) if isinstance(val, (list, tuple)) else val for key, val in spec["export"].items()
        }
        for c in cohorts:
            prep = prepared[c]
            primary_export, primary_analysis = prep["export"], prep["analysis"]
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
                print(
                    f"[{log}] {c} {o}: kept {acc['kept']['shared']} shared, {acc['kept']['primary_only']} primary-only, "
                    f"{acc['kept']['variant_only']} {v}-only; {acc['chrx']['kept']} of {acc['chrx']['input']} chrX "
                    f"features kept by {v}"
                )
            findings: list[str] = []
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
            (tmp / PREPROCESSING_SENSITIVITY_MARKDOWN).write_text(
                render_preprocessing_sensitivity(
                    rows, features, variant=v, cohort=c, settings=recorded_settings,
                    description=spec["description"], seed=ANALYSIS_SEED, permutations=ANALYSIS_PERMUTATIONS,
                    alpha=ANALYSIS_ALPHA, findings=findings,
                )
            )
            variant_files = {
                name: {"path": str(exp["dir"] / name), "sha256": _sha256_file(exp["dir"] / name)}
                for name in _ANALYSIS_INPUTS
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
                "outputs": sorted([
                    PREPROCESSING_SENSITIVITY_TABLE, PREPROCESSING_SENSITIVITY_MARKDOWN,
                    PREPROCESSING_SENSITIVITY_MANIFEST,
                ]),
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
                "manifest": manifest, "findings": findings,
            }
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
    input_dir: str | Path = "results/motco_input",
    analysis_dir: str | Path = "results/motco_analysis",
    rank_dir: str | Path = "results/motco_rank_sensitivity",
    pseudobulk_dir: str | Path = "results/pseudobulk",
    out_dir: str | Path = "results/motco_diagnostics",
    force: bool = False,
) -> dict:
    """Measure what the MOTCO latent axes and block axes correlate with, per cohort.

    Parameters
    ----------
    cohort:
        ``"cohort_30"``, ``"cohort_50"`` or ``"all"``.
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
    input_dir, analysis_dir, rank_dir = Path(input_dir), Path(analysis_dir), Path(rank_dir)
    pb, out_dir = Path(pseudobulk_dir), Path(out_dir)
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
        analyses = [(DIAGNOSTIC_PRIMARY_LABEL, _check_analysis_inputs(export_dir, analysis_dir / c))]
        rank_findings = []
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
    "render_preprocessing_sensitivity",
    "run_preprocessing_sensitivity",
    "TECHNICAL_GENE_SETS",
    "match_gene_sets",
    "technical_donor",
    "run_technical_covariates",
    "DIAGNOSTIC_VARIABLES",
    "render_axis_diagnostics",
    "run_axis_diagnostics",
]
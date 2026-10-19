"""Command-line interface for sead-sexdiff.

Provides a subcommand to download data from the public SEA-AD S3 bucket
using the `download_data` utility.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List

from .data import (
    download_data,
    DEFAULT_BUCKET,
    DATA_ROOT_ENV,
    MTG_PREFIXES,
    resolve_data_root,
    EXPORT_COHORTS,
    ANALYSIS_STATISTICS,
    ANALYSIS_FORCED_RANKS,
    PREPROCESSING_VARIANTS,
    DIAGNOSTIC_CORRELATION_THRESHOLD,
    run_analysis,
    run_rank_sensitivity,
    run_preprocessing_sensitivity,
    run_covariates,
    run_technical_covariates,
    run_axis_diagnostics,
    run_export,
    run_pseudobulk,
    write_subset_from_folder,
    pseudobulk_and_filter,
    differential_expression,
)


def _parse_suffixes(values: list[str] | None) -> list[str] | None:
    """Normalize suffix arguments.

    Supports passing multiple `--suffix` flags or a single comma-separated
    value. Returns `None` if no suffixes are provided so that `download_data`
    uses its own default behavior.
    """
    if not values:
        return None
    out: list[str] = []
    for v in values:
        if not v:
            continue
        # allow comma-separated
        parts = [p.strip() for p in v.split(",") if p.strip()]
        out.extend(parts)
    return out or None


def _parse_floats_csv(value: str | None) -> list[float] | None:
    """Parse a comma- or space-separated string of numbers into a list of floats.

    Returns None if the input is falsy.
    """
    if not value:
        return None
    parts = []
    # Allow commas and/or whitespace
    for token in value.replace(",", " ").split():
        if not token:
            continue
        parts.append(float(token))
    return parts or None


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="sead-sexdiff",
        description="Utilities for SEA-AD sex-difference analyses",
    )

    sub = parser.add_subparsers(dest="command", required=True)

    dl = sub.add_parser(
        "download",
        help="Download data from the SEA-AD public S3 bucket",
        description=(
            "Download public SEA-AD data files via anonymous S3 access. "
            "By default downloads CSV and HTML files under the given prefix. "
            "You may also specify explicit files to download, which overrides suffix filtering."
        ),
    )
    dl.add_argument(
        "--prefix",
        default="PFC/RNAseq/",
        help="S3 prefix (folder path) to search within the bucket (default: %(default)s)",
    )
    dl.add_argument(
        "--bucket",
        default=DEFAULT_BUCKET,
        help="S3 bucket name (default: %(default)s)",
    )
    dl.add_argument(
        "--suffix",
        action="append",
        dest="suffixes",
        metavar="SUFFIX",
        help=(
            "File suffix or pattern to include (e.g., .csv, .html, *.csv). "
            "Repeat the flag or use comma-separated values to provide multiple."
        ),
    )
    dl.add_argument(
        "--file",
        action="append",
        dest="files",
        metavar="FILE",
        help=(
            "Specific file under the prefix to download. Repeat the flag or use "
            "comma-separated values to provide multiple. When provided, this "
            "overrides --suffix filtering. You can pass keys relative to the "
            "prefix (e.g., public_datasets/SEA-AD-cell-annotation.2024-08-27.csv) "
            "or full keys that already start with the prefix."
        ),
    )
    dl.add_argument(
        "--omic",
        choices=["rna", "atac", "both"],
        default=None,
        help=(
            "Download the MTG per-donor objects for this omic "
            f"(rna: {MTG_PREFIXES['rna']}, atac: {MTG_PREFIXES['atac']}). "
            "Overrides --prefix and downloads every file under the prefix unless "
            "--suffix or --file is given. Without it, --prefix is used as before."
        ),
    )
    dl.add_argument(
        "--dest-root",
        default=None,
        type=Path,
        help=f"Local data root directory (default: ${DATA_ROOT_ENV}, else data)",
    )
    dl.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable progress bar with percentage and ETA during downloads",
    )

    # subset command: build a concatenated subset and write to disk
    sub_subset = sub.add_parser(
        "subset",
        help="Subset donor .h5ad files by subclass and write a single output",
        description=(
            "Load each donor .h5ad in a folder, subset by obs['Subclass'] == label, "
            "concatenate them, and write one compressed H5AD file."
        ),
    )
    sub_subset.add_argument(
        "--folder",
        default="data/PFC/RNAseq/donor_objects",
        type=Path,
        help=(
            "Input folder that contains donor .h5ad files (default: %(default)s)"
        ),
    )
    sub_subset.add_argument(
        "--subclass-label",
        default=None,
        help=(
            "Value to match in adata.obs['Subclass']. If omitted, a subset file "
            "is produced for EACH Subclass present across the input files."
        ),
    )
    sub_subset.add_argument(
        "--out-dir",
        default="results/cleaned_files",
        type=Path,
        help="Directory to write the output file (default: %(default)s)",
    )
    sub_subset.add_argument(
        "--out-name",
        default=None,
        help=(
            "Output filename (used only when --subclass-label is provided). "
            "If omitted, the filename is inferred from the subclass label "
            "(e.g., 'Microglia-PVM_subset.h5ad'). Ignored when producing all subclasses."
        ),
    )
    sub_subset.add_argument(
        "--compression",
        default="lzf",
        help=(
            "Compression for AnnData.write (e.g., lzf, gzip). Use 'none' to disable "
            "compression. (default: %(default)s)"
        ),
    )

    # differential expression command
    sub_de = sub.add_parser(
        "de",
        help="Run differential expression with pydeseq2 on pseudobulk data",
        description=(
            "Load a pseudobulk .h5ad file and, by default, run DE for ALL supertypes. "
            "Fits a DESeq2 model with design '~Sex + ADNC + Sex:ADNC', computes a "
            "contrast, writes one results table and volcano plot per supertype."
        ),
    )
    sub_de.add_argument(
        "--pseudobulk-file",
        default=Path("results/cleaned_files/microglia_pseudobulk.h5ad"),
        type=Path,
        help="Path to pseudobulk .h5ad file (default: %(default)s)",
    )
    sub_de.add_argument(
        "--supertype",
        default=None,
        help=(
            "Value in obs.Supertype to subset. If omitted, runs DE for ALL "
            "supertypes present in the data."
        ),
    )
    sub_de.add_argument(
        "--results-dir",
        default=Path("results/DE"),
        type=Path,
        help=(
            "Directory to write outputs. A CSV and a volcano PDF are produced per "
            "supertype (default: %(default)s)"
        ),
    )
    sub_de.add_argument(
        "--contrast",
        default=None,
        help=(
            "Contrast vector as comma-separated numbers (e.g., '0,0,0,0,0,0,0,1'). "
            "If omitted, the default contrast (0,0,0,0,0,0,0,1) is used."
        ),
    )

    # pseudobulk command: MTG astrocyte donor-level RNA + ATAC pseudobulks
    sub_pb = sub.add_parser(
        "pseudobulk",
        help="Build donor-level astrocyte pseudobulks from the MTG donor objects",
        description=(
            "Read the MTG RNA and/or ATAC donor objects one at a time, sum raw "
            "counts over astrocyte nuclei per donor, and write count matrices, "
            "the donor table with cohort flags, supertype composition, a cohort "
            "summary and a run manifest."
        ),
    )
    sub_pb.add_argument(
        "--omic",
        choices=["rna", "atac", "both"],
        default="both",
        help="Omic(s) to process (default: %(default)s)",
    )
    sub_pb.add_argument(
        "--data-root",
        default=None,
        type=Path,
        help=f"Local data root holding the MTG donor objects (default: ${DATA_ROOT_ENV}, else data)",
    )
    sub_pb.add_argument(
        "--out-dir",
        default=Path("results/pseudobulk"),
        type=Path,
        help="Output directory (default: %(default)s)",
    )
    sub_pb.add_argument(
        "--limit",
        type=int,
        default=None,
        help=(
            "Process only the first N donors (sorted, present in every selected "
            "omic); marks the run partial"
        ),
    )
    sub_pb.add_argument(
        "--donors",
        default=None,
        help="Comma-separated donor IDs to process; marks the run partial",
    )
    sub_pb.add_argument(
        "--force",
        action="store_true",
        help="Recompute donors already in the per-donor cache",
    )

    # covariates command: donor-level covariates from the MTG donor objects
    sub_cov = sub.add_parser(
        "covariates",
        help="Write the donor covariate table from the MTG donor object metadata",
        description=(
            "Read Age at Death, PMI, RIN and education from the obs of every MTG "
            "donor object (no matrices), code them numerically, cross-check RNA "
            "against ATAC, and write donor_covariates.csv plus "
            "covariates_manifest.json. Needs the raw donor objects (cluster)."
        ),
    )
    sub_cov.add_argument(
        "--omic",
        choices=["rna", "both"],
        default="both",
        help="Read RNA only, or RNA plus the ATAC cross-check (default: %(default)s)",
    )
    sub_cov.add_argument(
        "--data-root",
        default=None,
        type=Path,
        help=f"Local data root holding the MTG donor objects (default: ${DATA_ROOT_ENV}, else data)",
    )
    sub_cov.add_argument(
        "--out-dir",
        default=Path("results/pseudobulk"),
        type=Path,
        help="Output directory (default: %(default)s)",
    )

    # technical-covariates command: donor technical table of the astrocyte nuclei
    sub_tc = sub.add_parser(
        "technical-covariates",
        help="Write the donor technical table of the astrocyte nuclei (gene-set UMI fractions, nucleus and library QC)",
        description=(
            "For every MTG RNA donor object, select the astrocyte nuclei exactly as the "
            "pseudobulk step does (shared code path), read their rows of the raw UMIs layer "
            "and their obs QC columns (never the normalized X), and write per donor the "
            "nucleus means and SDs and the UMI-weighted pseudobulk fractions of the predeclared "
            "gene sets (ribosomal, mitochondrial, hemoglobin), the within-donor supertype share "
            "of each fraction, the nucleus and library QC means and the number of libraries; "
            "from the ATAC objects the library-level ATAC_* metrics (obs only). The pseudobulk "
            "fractions are cross-checked against rna_counts.h5ad and the nucleus counts against "
            "the pseudobulk cache and donor_table.csv. Writes donor_technical.csv, "
            "nucleus_technical_summary.csv and technical_manifest.json; never writes "
            "donor_covariates.csv or covariates_manifest.json. Resumable per donor "
            "(per_donor_technical/). Needs the raw donor objects (cluster)."
        ),
    )
    sub_tc.add_argument(
        "--omic",
        choices=["rna", "both"],
        default="both",
        help="RNA objects only, or RNA plus the ATAC library metrics (default: %(default)s)",
    )
    sub_tc.add_argument(
        "--data-root",
        default=None,
        type=Path,
        help=f"Local data root holding the MTG donor objects (default: ${DATA_ROOT_ENV}, else data)",
    )
    sub_tc.add_argument(
        "--pseudobulk-dir",
        default=Path("results/pseudobulk"),
        type=Path,
        help="Directory with the pseudobulk outputs, read for the cross-check (default: %(default)s)",
    )
    sub_tc.add_argument(
        "--out-dir",
        default=Path("results/pseudobulk"),
        type=Path,
        help="Output directory (default: %(default)s)",
    )
    sub_tc.add_argument(
        "--force",
        action="store_true",
        help="Recompute donors already in the per-donor cache",
    )

    # export command: MOTCO-ready blocks per cohort
    sub_ex = sub.add_parser(
        "export",
        help="Preprocess the pseudobulks into MOTCO input blocks per cohort",
        description=(
            "Per cohort and omic: restrict to the cohort, drop chrX/Y/M, keep "
            "features with CPM >= 1 in >= k donors (k = smallest sex x stage "
            "cell), TMM log2-CPM (prior count 2), keep the 5,000 most variable "
            "features. Writes rna.csv, atac.csv, metadata.csv, covariates.csv, "
            "features_<omic>.csv and manifest.json. Reads only the pseudobulk "
            "outputs plus the GENCODE annotation (downloaded once)."
        ),
    )
    sub_ex.add_argument(
        "--pseudobulk-dir",
        default=Path("results/pseudobulk"),
        type=Path,
        help="Directory with the pseudobulk and covariates outputs (default: %(default)s)",
    )
    sub_ex.add_argument(
        "--out-dir",
        default=Path("results/motco_input"),
        type=Path,
        help="Export root; one subdirectory per cohort (default: %(default)s)",
    )
    sub_ex.add_argument(
        "--cohort",
        choices=list(EXPORT_COHORTS) + ["all"],
        default="all",
        help="Cohort(s) to export (default: %(default)s)",
    )
    sub_ex.add_argument(
        "--data-root",
        default=None,
        type=Path,
        help=(
            "Data root holding <root>/annotation/ for the GENCODE GTF "
            f"(default: ${DATA_ROOT_ENV}, else data)"
        ),
    )
    sub_ex.add_argument(
        "--force",
        action="store_true",
        help="Replace existing cohort directories",
    )

    # analyze command: covariate-adjusted MOTCO trajectory analysis per cohort
    sub_an = sub.add_parser(
        "analyze",
        help="Run the covariate-adjusted MOTCO trajectory analysis per cohort",
        description=(
            "Per cohort, from the export outputs only: compare the female and male "
            "astrocyte trajectories across ADNC stages 0, 1, 2 with MOTCO's "
            "analyze_trajectories (RNA + ATAC blocks, net-displacement orientation, "
            "CV-selected rank, 999 RRPP permutations, fixed seed), adjusting for "
            "age_at_death, pmi, rin, years_education, log10_n_rna and log10_n_atac. "
            "Writes MOTCO's artifact plus analysis_metadata.csv, "
            "analysis_manifest.json and summary.md."
        ),
    )
    sub_an.add_argument(
        "--cohort",
        choices=list(EXPORT_COHORTS) + ["all"],
        default="all",
        help="Cohort(s) to analyse (default: %(default)s)",
    )
    sub_an.add_argument(
        "--input-dir",
        default=Path("results/motco_input"),
        type=Path,
        help="Export root with one subdirectory per cohort (default: %(default)s)",
    )
    sub_an.add_argument(
        "--out-dir",
        default=Path("results/motco_analysis"),
        type=Path,
        help="Analysis root; one subdirectory per cohort (default: %(default)s)",
    )
    sub_an.add_argument(
        "--n-jobs",
        type=int,
        default=1,
        help=(
            "Parallel jobs for MOTCO's CV and RRPP, -1 = all CPUs. Values other than 1 "
            "change the RRPP draws (p-values); results are reported at 1 (default: %(default)s)"
        ),
    )
    sub_an.add_argument(
        "--force",
        action="store_true",
        help="Replace existing cohort directories",
    )

    # rank-sensitivity command: forced-rank sweep beside the CV-selected analysis
    ladder = ", ".join(str(k) for k in ANALYSIS_FORCED_RANKS)
    sub_rs = sub.add_parser(
        "rank-sensitivity",
        help="Rerun the MOTCO trajectory analysis at a fixed ladder of forced latent ranks",
        description=(
            "Per cohort, from the export outputs and the analyze step's outputs: rerun the "
            f"covariate-adjusted trajectory analysis with the latent rank forced to {ladder} "
            "(a constant, not a flag), everything else equal to the analyze step (same blocks, "
            "covariates, orientation, permutations and seed at every rank). The CV-selected "
            "analysis is read, checked against the inputs and settings, and placed as the "
            "reference row of the cross-rank table; it is never rerun or modified and stays the "
            "primary result. Rank 1 cannot be forced. Writes rank_<kk>/ analysis directories, "
            "rank_sensitivity.csv, rank_sensitivity.md and rank_sensitivity_manifest.json."
        ),
    )
    sub_rs.add_argument(
        "--cohort",
        choices=list(EXPORT_COHORTS) + ["all"],
        default="all",
        help="Cohort(s) to sweep (default: %(default)s)",
    )
    sub_rs.add_argument(
        "--input-dir",
        default=Path("results/motco_input"),
        type=Path,
        help="Export root with one subdirectory per cohort (default: %(default)s)",
    )
    sub_rs.add_argument(
        "--analysis-dir",
        default=Path("results/motco_analysis"),
        type=Path,
        help="Analysis root of the analyze step; read as the reference, never written (default: %(default)s)",
    )
    sub_rs.add_argument(
        "--out-dir",
        default=Path("results/motco_rank_sensitivity"),
        type=Path,
        help="Sweep root; one subdirectory per cohort (default: %(default)s)",
    )
    sub_rs.add_argument(
        "--n-jobs",
        type=int,
        default=1,
        help=(
            "Parallel jobs for MOTCO's RRPP, -1 = all CPUs. Values other than 1 change the RRPP "
            "draws (p-values) at every rank; results are reported at 1 (default: %(default)s)"
        ),
    )
    sub_rs.add_argument(
        "--force",
        action="store_true",
        help="Replace existing cohort directories under --out-dir",
    )

    # preprocessing-sensitivity command: export variants beside the primary analysis
    variants = ", ".join(PREPROCESSING_VARIANTS)
    sub_ps = sub.add_parser(
        "preprocessing-sensitivity",
        help="Rerun the export with one registered setting changed, then the analysis, beside the primary",
        description=(
            "Per variant of the registry PREPROCESSING_VARIANTS (a constant, not a flag; now: "
            f"{variants}) and per cohort: rerun the export from the pseudobulk outputs with the variant's "
            "setting changed and everything else equal to the export step (keep_chrx removes chrY and chrM "
            "only, so chrX features are kept and the abundance filter, TMM, log2-CPM and the variance ranking "
            "are recomputed with them present), then rerun the covariate-adjusted trajectory analysis on the "
            "variant blocks with the analyze step's settings (CV-selected rank, 999 RRPP permutations, same "
            "seed), and compare it with the primary analysis. The primary export and analysis are read and "
            "checked (pseudobulk, covariates and annotation hashes, thresholds, k, cells; export hashes, MOTCO "
            "version, seed, permutations, CV selection), never rerun or modified, and stay the primary result. "
            "Writes <out-dir>/<variant>/motco_input/<cohort>/ (the variant export), "
            "<out-dir>/<variant>/motco_analysis/<cohort>/ (the variant analysis) and "
            "<out-dir>/<variant>/comparison/<cohort>/ (sensitivity.csv, sensitivity.md, sensitivity_manifest.json)."
        ),
    )
    sub_ps.add_argument(
        "--variant",
        choices=list(PREPROCESSING_VARIANTS) + ["all"],
        default="all",
        help="Registered variant(s) to run (default: %(default)s)",
    )
    sub_ps.add_argument(
        "--cohort",
        choices=list(EXPORT_COHORTS) + ["all"],
        default="all",
        help="Cohort(s) to run (default: %(default)s)",
    )
    sub_ps.add_argument(
        "--pseudobulk-dir",
        default=Path("results/pseudobulk"),
        type=Path,
        help="Directory with the pseudobulk and covariates outputs, the export's inputs (default: %(default)s)",
    )
    sub_ps.add_argument(
        "--input-dir",
        default=Path("results/motco_input"),
        type=Path,
        help="Export root of the export step; read as the primary export, never written (default: %(default)s)",
    )
    sub_ps.add_argument(
        "--analysis-dir",
        default=Path("results/motco_analysis"),
        type=Path,
        help="Analysis root of the analyze step; read as the primary analysis, never written (default: %(default)s)",
    )
    sub_ps.add_argument(
        "--out-dir",
        default=Path("results/motco_sensitivity"),
        type=Path,
        help="Sensitivity root; one subdirectory per variant (default: %(default)s)",
    )
    sub_ps.add_argument(
        "--data-root",
        default=None,
        type=Path,
        help=(
            "Data root holding <root>/annotation/ for the GENCODE GTF "
            f"(default: ${DATA_ROOT_ENV}, else data)"
        ),
    )
    sub_ps.add_argument(
        "--n-jobs",
        type=int,
        default=1,
        help=(
            "Parallel jobs for MOTCO's CV and RRPP, -1 = all CPUs. Values other than 1 change the RRPP "
            "draws (p-values); results are reported at 1 (default: %(default)s)"
        ),
    )
    sub_ps.add_argument(
        "--force",
        action="store_true",
        help="Replace a cohort's existing variant export, analysis and comparison directories under --out-dir",
    )
    # axis-diagnostics command: what the latent and block axes correlate with
    sub_ad = sub.add_parser(
        "axis-diagnostics",
        help="Measure what the MOTCO latent axes and the exported blocks' principal axes correlate with",
        description=(
            "Per cohort, from kept outputs only (the export, the analyze step's analysis, the technical table "
            "and, when present, the rank-sensitivity analyses): Spearman correlations of every latent component "
            "with the predeclared variables (the technical table's donor columns, the export's TMM factors and "
            "log10 library sizes, the six analysis covariates, stage and sex) and the partial correlations given "
            "the six analysis covariates; the same for the first three principal components of each z-scored "
            "block; the share of the trajectory geometry (net displacements, sex differences) along each "
            "component; and the feature-selection accounting (kept features' mean log2-CPM quantiles, variance "
            "share, variance at the cut, share correlated with pls_0). Each analysis manifest's input hashes are "
            "checked against the export files first; nothing under the input directories is written. Rows with "
            f"|rho| >= {DIAGNOSTIC_CORRELATION_THRESHOLD} are flagged (a reporting flag, not an adoption rule; no "
            "p-values). Writes axis_correlations.csv, block_pcs.csv, geometry_shares.csv, feature_selection.csv, "
            "diagnostics.md and diagnostics_manifest.json. Its outputs feed no export or analysis command: "
            "preprocessing variants are predeclared in a change proposal, never derived from a diagnostic value."
        ),
    )
    sub_ad.add_argument(
        "--cohort",
        choices=list(EXPORT_COHORTS) + ["all"],
        default="all",
        help="Cohort(s) to diagnose (default: %(default)s)",
    )
    sub_ad.add_argument(
        "--input-dir",
        default=Path("results/motco_input"),
        type=Path,
        help="Export root with one subdirectory per cohort (default: %(default)s)",
    )
    sub_ad.add_argument(
        "--analysis-dir",
        default=Path("results/motco_analysis"),
        type=Path,
        help="Analysis root of the analyze step; read as the primary analysis, never written (default: %(default)s)",
    )
    sub_ad.add_argument(
        "--rank-dir",
        default=Path("results/motco_rank_sensitivity"),
        type=Path,
        help="Sweep root of the rank-sensitivity step; its rank_*/ analyses are read when present (default: %(default)s)",
    )
    sub_ad.add_argument(
        "--pseudobulk-dir",
        default=Path("results/pseudobulk"),
        type=Path,
        help="Directory holding donor_technical.csv from the technical-covariates step (default: %(default)s)",
    )
    sub_ad.add_argument(
        "--out-dir",
        default=Path("results/motco_diagnostics"),
        type=Path,
        help="Diagnostics root; one subdirectory per cohort (default: %(default)s)",
    )
    sub_ad.add_argument(
        "--force",
        action="store_true",
        help="Replace existing cohort directories under --out-dir",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "download":
        suffixes = _parse_suffixes(args.suffixes)
        files = _parse_suffixes(args.files)
        root = resolve_data_root(args.dest_root)
        print(f"[download] Data root: {root.resolve()}")
        if args.omic is None:
            prefixes = [args.prefix]
            if suffixes is None:
                suffixes = [".csv", ".html"]
        else:
            omics = ["rna", "atac"] if args.omic == "both" else [args.omic]
            prefixes = [MTG_PREFIXES[o] for o in omics]
        n_down = n_skip = n_failed = 0
        for prefix in prefixes:
            prefix_files = files
            if files and len(prefixes) > 1:
                # With --omic both, route each file to the omic named in it
                token = "RNAseq" if prefix == MTG_PREFIXES["rna"] else "ATACseq"
                other = "ATACseq" if token == "RNAseq" else "RNAseq"
                prefix_files = [f for f in files if token in f or other not in f]
                if not prefix_files:
                    continue
            print(f"[download] Prefix: {prefix}")
            paths = download_data(
                prefix=prefix,
                bucket=args.bucket,
                suffixes=suffixes,
                files=prefix_files,
                dest_root=root,
                show_progress=not args.no_progress,
            )
            for p in paths:
                print(p)
            if paths.manifest_path is not None:
                print(f"[download] Manifest: {paths.manifest_path}")
            n_down += len(paths)
            n_skip += len(paths.skipped)
            n_failed += len(paths.failed)
        if len(prefixes) > 1:
            print(f"[download] Total: {n_down} downloaded, {n_skip} skipped, {n_failed} failed")
        return 1 if n_failed else 0

    if args.command == "subset":
        compression = None if (str(args.compression).lower() in {"none", "no", "false", "0"}) else args.compression

        print("[subset] Starting subset and concatenate step...")
        print(f"[subset] Folder: {args.folder}")
        print(f"[subset] Subclass label: {args.subclass_label if args.subclass_label is not None else 'ALL (one file per Subclass)'}")
        print(f"[subset] Output directory: {args.out_dir}")
        if args.subclass_label is not None:
            print(f"[subset] Output name: {args.out_name or '(auto from subclass)'}")
        else:
            print("[subset] Output name: (auto: '<Subclass>_subset.h5ad' per subclass)")
        print(f"[subset] Compression: {compression or 'none'}")
        print("[subset] Building concatenated subset (this may take a while)...")
        out_paths = write_subset_from_folder(
            folder=args.folder,
            subclass_label=args.subclass_label,
            out_dir=args.out_dir,
            out_name=args.out_name,
            compression=compression,
        )
        if isinstance(out_paths, list):
            print("[subset] Wrote the following subset files:")
            for p in out_paths:
                print(f"  - {p}")
            print("[subset] Note: Renaming obs column 'Overall AD neuropathological Change' to 'ADNC' is applied before writing when present.")
            print("[subset] Pseudobulk files were generated alongside subset files.")
        else:
            print(f"[subset] Wrote subset file to: {out_paths}")
            print("[subset] Note: Renaming obs column 'Overall AD neuropathological Change' to 'ADNC' is applied before writing when present.")
            print("[subset] Pseudobulk file was generated alongside subset file.")
        return 0

    if args.command == "de":
        print("[de] Running differential expression analysis...")
        contrast_list = _parse_floats_csv(args.contrast)
        if contrast_list is None:
            outputs = differential_expression(
                pseudobulk_file=args.pseudobulk_file,
                supertype=args.supertype,
                results_dir=args.results_dir,
            )
        else:
            outputs = differential_expression(
                pseudobulk_file=args.pseudobulk_file,
                supertype=args.supertype,
                results_dir=args.results_dir,
                contrast=contrast_list,
            )
        if isinstance(outputs, dict):
            # Print a line per supertype
            for st, paths in outputs.items():
                csv_p = paths.get("results_csv")
                pdf_p = paths.get("volcano_pdf")
                if csv_p:
                    print(f"[de] {st} results: {csv_p}")
                if pdf_p:
                    print(f"[de] {st} volcano: {pdf_p}")
        print("[de] Done.")
        return 0

    if args.command == "pseudobulk":
        donors = _parse_suffixes([args.donors]) if args.donors else None
        print("[pseudobulk] Building astrocyte donor pseudobulks...")
        outputs = run_pseudobulk(
            omic=args.omic,
            data_root=args.data_root,
            out_dir=args.out_dir,
            limit=args.limit,
            donors=donors,
            force=args.force,
        )
        for name in ("rna_counts", "atac_counts", "donor_table", "supertype_composition", "cohort_summary", "manifest"):
            if name in outputs:
                print(f"[pseudobulk] Wrote {outputs[name]}")
        print("[pseudobulk] Done.")
        return 0

    if args.command == "covariates":
        print("[covariates] Reading donor covariates...")
        outputs = run_covariates(omic=args.omic, data_root=args.data_root, out_dir=args.out_dir)
        for name in ("covariates", "manifest"):
            print(f"[covariates] Wrote {outputs[name]}")
        print("[covariates] Done.")
        return 0

    if args.command == "technical-covariates":
        print("[technical-covariates] Summarizing the astrocyte nuclei of every donor object...")
        outputs = run_technical_covariates(
            args.omic,
            data_root=args.data_root,
            pseudobulk_dir=args.pseudobulk_dir,
            out_dir=args.out_dir,
            force=args.force,
        )
        for name in ("table", "summary", "manifest"):
            print(f"[technical-covariates] Wrote {outputs[name]}")
        print("[technical-covariates] Done.")
        return 0

    if args.command == "export":
        print("[export] Exporting MOTCO input blocks...")
        run_export(
            pseudobulk_dir=args.pseudobulk_dir,
            out_dir=args.out_dir,
            cohort=args.cohort,
            data_root=args.data_root,
            force=args.force,
        )
        print("[export] Done.")
        return 0

    if args.command == "analyze":
        import pandas as pd

        print("[analyze] Running the MOTCO trajectory analysis...")
        results = run_analysis(
            args.cohort,
            input_dir=args.input_dir,
            out_dir=args.out_dir,
            n_jobs=args.n_jobs,
            force=args.force,
        )
        rows = []
        for cohort, res in results.items():
            record = res["record"]
            for s in ANALYSIS_STATISTICS:
                rows.append({
                    "cohort": cohort,
                    "statistic": s,
                    "observed": record["statistics"].get(s),
                    "p": record["p_values"].get(s),
                    "rank": record["latent"]["selected_lv"] if s == "angle" else "",
                })
        print("[analyze] Statistics (p from RRPP; rank = CV-selected latent rank):")
        print(pd.DataFrame(rows).to_string(index=False))
        print("[analyze] Done.")
        return 0

    if args.command == "rank-sensitivity":
        import pandas as pd

        print("[rank-sensitivity] Running the forced-rank sweep...")
        results = run_rank_sensitivity(
            args.cohort,
            input_dir=args.input_dir,
            analysis_dir=args.analysis_dir,
            out_dir=args.out_dir,
            n_jobs=args.n_jobs,
            force=args.force,
        )
        columns = ["cohort", "component_selection", "rank"]
        for s in ANALYSIS_STATISTICS:
            columns += [s, f"{s}_p", f"{s}_rejected"]
        table = pd.concat([res["table"][columns] for res in results.values()], ignore_index=True)
        print("[rank-sensitivity] Cross-rank table (cv = the CV-selected analysis, primary; forced = sensitivity):")
        print(table.to_string(index=False))
        print("[rank-sensitivity] Done.")
        return 0

    if args.command == "preprocessing-sensitivity":
        import pandas as pd

        print("[preprocessing-sensitivity] Running the preprocessing sensitivity analysis...")
        results = run_preprocessing_sensitivity(
            args.variant,
            args.cohort,
            pseudobulk_dir=args.pseudobulk_dir,
            input_dir=args.input_dir,
            analysis_dir=args.analysis_dir,
            out_dir=args.out_dir,
            data_root=args.data_root,
            n_jobs=args.n_jobs,
            force=args.force,
        )
        columns = ["variant", "cohort", "component_selection", "rank"]
        for s in ANALYSIS_STATISTICS:
            columns += [s, f"{s}_p", f"{s}_rejected"]
        table = pd.concat(
            [res["table"][columns] for per_cohort in results.values() for res in per_cohort.values()],
            ignore_index=True,
        )
        print("[preprocessing-sensitivity] Comparison table (primary = the analyze step's result, primary; "
              "the variant rows are sensitivity):")
        print(table.to_string(index=False))
        print("[preprocessing-sensitivity] Done.")
        return 0

    if args.command == "axis-diagnostics":
        print("[axis-diagnostics] Measuring the latent and block axes...")
        results = run_axis_diagnostics(
            args.cohort,
            input_dir=args.input_dir,
            analysis_dir=args.analysis_dir,
            rank_dir=args.rank_dir,
            pseudobulk_dir=args.pseudobulk_dir,
            out_dir=args.out_dir,
            force=args.force,
        )
        for cohort, res in results.items():
            axis = res["tables"]["axis"]
            flagged = axis[(axis["analysis"] == "primary") & (axis["component"] == "pls_0") & (axis["flag_raw"] | axis["flag_partial"])]
            print(f"[axis-diagnostics] {cohort}: {len(flagged)} variables flagged against the primary pls_0 "
                  f"(|rho| >= {DIAGNOSTIC_CORRELATION_THRESHOLD}, raw or partial; no p-values):")
            if len(flagged):
                print(flagged[["variable", "n", "spearman", "partial_spearman"]].to_string(index=False))
        print("[axis-diagnostics] Done.")
        return 0

    parser.error("Unknown command")
    return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())

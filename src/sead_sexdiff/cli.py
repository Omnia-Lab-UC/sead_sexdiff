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
    CELL_TYPES,
    CELL_TYPE_SUBDIRS,
    DEFAULT_CELL_TYPE,
    cell_type_dirs,
    resolve_data_root,
    EXPORT_COHORTS,
    ANALYSIS_STATISTICS,
    ANALYSIS_FORCED_RANKS,
    PREPROCESSING_VARIANTS,
    DIAGNOSTIC_CORRELATION_THRESHOLD,
    LATENT_SENSITIVITY_ARMS,
    LATENT_SENSITIVITY_ADDED_COVARIATES,
    LATENT_SENSITIVITY_CONDITIONING,
    LATENT_MOTCO_MIN_VERSION,
    STAGE_DE_ALPHA,
    STAGE_DE_MODELS,
    STAGE_DE_LIMMA_JOINT_TESTS,
    POSITIVE_CONTROL_BASES,
    POSITIVE_CONTROL_CELL_TYPES,
    POSITIVE_CONTROL_DEFAULT_BASIS,
    POSITIVE_CONTROL_PRIMARY_COHORT,
    run_analysis,
    run_rank_sensitivity,
    run_preprocessing_sensitivity,
    run_covariates,
    run_technical_covariates,
    run_axis_diagnostics,
    run_latent_sensitivity,
    run_stage_de,
    run_stage_de_limma,
    run_positive_control,
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


def _add_cell_type_argument(parser: argparse.ArgumentParser) -> None:
    """``--cell-type``: a registered name of ``CELL_TYPES`` (the mapping itself is not an option)."""
    parser.add_argument(
        "--cell-type",
        choices=list(CELL_TYPES),
        default=DEFAULT_CELL_TYPE,
        help=(
            "Registered cell type (CELL_TYPES): the population and the default directories under its results "
            "root; explicit directory flags take precedence (default: %(default)s)"
        ),
    )


def _dir_default(key: str) -> str:
    """Help text of a directory flag that defaults to the cell type's directory."""
    return (
        f"(default: <results root of --cell-type>/{CELL_TYPE_SUBDIRS[key]}, "
        f"{cell_type_dirs(DEFAULT_CELL_TYPE)[key]} for {DEFAULT_CELL_TYPE})"
    )


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
        help="Build donor-level pseudobulks of a registered cell type (default astrocyte) from the MTG donor objects",
        description=(
            "Read the MTG RNA and/or ATAC donor objects one at a time, sum raw "
            "counts over the nuclei of the cell type (--cell-type, default astrocyte) "
            "per donor, and write count matrices, "
            "the donor table with cohort flags, supertype composition, a cohort "
            "summary and a run manifest."
        ),
    )
    _add_cell_type_argument(sub_pb)
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
        default=None,
        type=Path,
        help="Output directory " + _dir_default("pseudobulk"),
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
    _add_cell_type_argument(sub_cov)
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
        default=None,
        type=Path,
        help="Output directory " + _dir_default("pseudobulk"),
    )

    # technical-covariates command: donor technical table of the cell type's nuclei
    sub_tc = sub.add_parser(
        "technical-covariates",
        help="Write the donor technical table of the cell type's nuclei (gene-set UMI fractions, nucleus and library QC)",
        description=(
            "For every MTG RNA donor object, select the nuclei of the cell type (--cell-type, "
            "default astrocyte) exactly as the pseudobulk step does (shared code path; the "
            "pseudobulk manifest must record the same cell type), read their rows of the raw UMIs layer "
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
    _add_cell_type_argument(sub_tc)
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
        default=None,
        type=Path,
        help="Directory with the pseudobulk outputs, read for the cross-check " + _dir_default("pseudobulk"),
    )
    sub_tc.add_argument(
        "--out-dir",
        default=None,
        type=Path,
        help="Output directory " + _dir_default("pseudobulk"),
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
    _add_cell_type_argument(sub_ex)
    sub_ex.add_argument(
        "--pseudobulk-dir",
        default=None,
        type=Path,
        help="Directory with the pseudobulk and covariates outputs " + _dir_default("pseudobulk"),
    )
    sub_ex.add_argument(
        "--out-dir",
        default=None,
        type=Path,
        help="Export root; one subdirectory per cohort " + _dir_default("motco_input"),
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
            "trajectories of the cell type (--cell-type, default astrocyte; the export "
            "manifest must record it) across ADNC stages 0, 1, 2 with MOTCO's "
            "analyze_trajectories (RNA + ATAC blocks, net-displacement orientation, "
            "CV-selected rank, 999 RRPP permutations, fixed seed), adjusting for "
            "age_at_death, pmi, rin, years_education, log10_n_rna and log10_n_atac. "
            "Writes MOTCO's artifact plus analysis_metadata.csv, "
            "analysis_manifest.json and summary.md."
        ),
    )
    _add_cell_type_argument(sub_an)
    sub_an.add_argument(
        "--cohort",
        choices=list(EXPORT_COHORTS) + ["all"],
        default="all",
        help="Cohort(s) to analyse (default: %(default)s)",
    )
    sub_an.add_argument(
        "--input-dir",
        default=None,
        type=Path,
        help="Export root with one subdirectory per cohort " + _dir_default("motco_input"),
    )
    sub_an.add_argument(
        "--out-dir",
        default=None,
        type=Path,
        help="Analysis root; one subdirectory per cohort " + _dir_default("motco_analysis"),
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
            "registered setting(s) changed and everything else equal to the export step (keep_chrx removes chrY "
            "and chrM only; exclude_families removes the ribosomal, mitochondrial and hemoglobin RNA genes before "
            "the abundance filter; rna_trend_ratio ranks RNA features by variance over a binned mean-variance "
            "trend; atac_mean_accessibility ranks ATAC peaks by mean log2-CPM; residualized regresses each "
            "feature on the export's log10 library sizes and the technical table's gene-set fractions after "
            "log2-CPM; revised_candidate composes the four; trend_normalized corrects each donor's log2-CPM by "
            "limma's fast cyclic loess against the per-feature mean after log2-CPM; revised_candidate_trend adds it "
            "to revised_candidate's four), then rerun the covariate-adjusted trajectory "
            "analysis on the variant blocks with the analyze step's settings (CV-selected rank, 999 RRPP "
            "permutations, same seed), and compare it with the primary analysis. The primary export and analysis "
            "are read and checked (pseudobulk, covariates and annotation hashes, thresholds, k, cells; export "
            "hashes, MOTCO version, seed, permutations, CV selection; for the gene-set, residualized, "
            "trend-normalized and candidate variants also the technical table against the pseudobulk manifest), "
            "never rerun or modified, and stay the primary result. For revised_candidate, trend_normalized and "
            "revised_candidate_trend the predeclared adoption rule ADOPTION_RULE is evaluated from the generated "
            "files into adoption_check.csv (no p-value is read), and ADOPTION_PRECEDENCE, applied to the holds "
            "values of trend_normalized and revised_candidate_trend in both cohorts, is recorded in their comparison "
            "manifests. Writes <out-dir>/<variant>/motco_input/<cohort>/ (the variant export, plus "
            "residualization.csv when it residualizes and normalization_trend.csv when it trend-normalizes), "
            "<out-dir>/<variant>/motco_analysis/<cohort>/ (the variant analysis) and "
            "<out-dir>/<variant>/comparison/<cohort>/ (sensitivity.csv, sensitivity.md, sensitivity_manifest.json "
            "and, for the candidates, adoption_check.csv). With --variant all, a variant x cohort whose comparison "
            "directory exists is skipped with a message unless --force is given, so all is resumable."
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
        help=(
            "Replace a cohort's existing variant export, analysis and comparison directories under --out-dir "
            "(without it, --variant all skips complete variant x cohort pairs; a named variant with existing output is an error)"
        ),
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
    _add_cell_type_argument(sub_ad)
    sub_ad.add_argument(
        "--cohort",
        choices=list(EXPORT_COHORTS) + ["all"],
        default="all",
        help="Cohort(s) to diagnose (default: %(default)s)",
    )
    sub_ad.add_argument(
        "--input-dir",
        default=None,
        type=Path,
        help="Export root with one subdirectory per cohort " + _dir_default("motco_input"),
    )
    sub_ad.add_argument(
        "--analysis-dir",
        default=None,
        type=Path,
        help="Analysis root of the analyze step; read as the primary analysis, never written " + _dir_default("motco_analysis"),
    )
    sub_ad.add_argument(
        "--rank-dir",
        default=None,
        type=Path,
        help="Sweep root of the rank-sensitivity step; its rank_*/ analyses are read when present " + _dir_default("motco_rank_sensitivity"),
    )
    sub_ad.add_argument(
        "--pseudobulk-dir",
        default=None,
        type=Path,
        help="Directory holding donor_technical.csv from the technical-covariates step " + _dir_default("pseudobulk"),
    )
    sub_ad.add_argument(
        "--out-dir",
        default=None,
        type=Path,
        help="Diagnostics root; one subdirectory per cohort " + _dir_default("motco_diagnostics"),
    )
    sub_ad.add_argument(
        "--force",
        action="store_true",
        help="Replace existing cohort directories under --out-dir",
    )

    # latent-sensitivity command: covariate set and latent basis beside the primary analysis
    arms = ", ".join(LATENT_SENSITIVITY_ARMS)
    added = ", ".join(LATENT_SENSITIVITY_ADDED_COVARIATES)
    sub_ls = sub.add_parser(
        "latent-sensitivity",
        help="Rerun the analysis with the extended covariate set, without and with MOTCO's latent adjustment",
        description=(
            f"Per cohort, from the export, the analyze step's analysis and the technical table: run the arms of the "
            f"registry LATENT_SENSITIVITY_ARMS (a constant, not a flag; now: {arms}) on the unchanged primary export. "
            f"Both enter the six analysis covariates plus {added} (raw fractions from donor_technical.csv; the set is "
            "a constant); extended_covariates keeps the covariate-blind latent basis, latent_adjusted fits the basis "
            f"on features residualized on the eight covariates conditional on {LATENT_SENSITIVITY_CONDITIONING} "
            f"(MOTCO >= {LATENT_MOTCO_MIN_VERSION}). Everything else equals the analyze step (blocks, stage order, "
            "orientation, CV-selected rank, 999 RRPP permutations, same seed). The primary analysis is read and "
            "checked (export hashes, MOTCO version, seed, permutations, CV selection), never rerun or modified, and "
            "stays the primary result; the technical table is checked against the pseudobulk manifest, "
            "TECHNICAL_GENE_SETS and completeness. Before the arms run, mediator_check.csv gives each added "
            "covariate's Spearman correlation with stage (all donors and per sex; flagged at |rho| >= "
            f"{DIAGNOSTIC_CORRELATION_THRESHOLD}, no p-value, never a gate). Writes <out-dir>/<arm>/motco_analysis/<cohort>/ "
            "(the files analyze writes) and <out-dir>/comparison/<cohort>/ (latent_sensitivity.csv, "
            "latent_sensitivity.md, mediator_check.csv, latent_sensitivity_manifest.json). The arms run per registered "
            "cell type (--cell-type; the registry, the covariate set and every setting are the same for every cell "
            "type); every manifest read must record the cell type."
        ),
    )
    _add_cell_type_argument(sub_ls)
    sub_ls.add_argument(
        "--cohort",
        choices=list(EXPORT_COHORTS) + ["all"],
        default="all",
        help="Cohort(s) to run (default: %(default)s)",
    )
    sub_ls.add_argument(
        "--input-dir",
        default=None,
        type=Path,
        help=f"Export root of the export step; read as the primary export, never written {_dir_default('motco_input')}",
    )
    sub_ls.add_argument(
        "--analysis-dir",
        default=None,
        type=Path,
        help=(
            "Analysis root of the analyze step; read as the primary analysis, never written "
            f"{_dir_default('motco_analysis')}"
        ),
    )
    sub_ls.add_argument(
        "--pseudobulk-dir",
        default=None,
        type=Path,
        help=(
            "Directory holding manifest.json, donor_technical.csv and technical_manifest.json; read, never written "
            f"{_dir_default('pseudobulk')}"
        ),
    )
    sub_ls.add_argument(
        "--out-dir",
        default=None,
        type=Path,
        help=f"Output root; one subdirectory per arm plus comparison/ {_dir_default('motco_latent_sensitivity')}",
    )
    sub_ls.add_argument(
        "--n-jobs",
        type=int,
        default=1,
        help=(
            "Parallel jobs for MOTCO's CV and RRPP, -1 = all CPUs. Values other than 1 change the RRPP "
            "draws (p-values); results are reported at 1 (default: %(default)s)"
        ),
    )
    sub_ls.add_argument(
        "--force",
        action="store_true",
        help="Replace a cohort's existing arm and comparison directories under --out-dir",
    )

    # stage-de command: per-feature DESeq2 of stage and sex x stage beside the primary analysis
    models = ", ".join(STAGE_DE_MODELS)
    sub_sd = sub.add_parser(
        "stage-de",
        help="Per-feature DESeq2 of ADNC stage and sex x stage on the cohort pseudobulks, beside the primary analysis",
        description=(
            "Per cohort and omic, from the pseudobulk counts, the primary export and the primary analysis: fit the "
            f"covariate-adjusted DESeq2 models of the registry STAGE_DE_MODELS (a constant, not a flag; now: {models}) "
            "with pydeseq2 on the raw counts of the cohort donors, on the features that passed the export's "
            "chromosome and abundance filters (features_<omic>.csv status variance or kept; size factors and "
            "dispersions fitted on that set), with sex (reference Female), stage (reference 0) and the omic's "
            "covariates (age_at_death, pmi, rin, years_education, log10_n_<omic>). One Wald test per tested "
            f"coefficient, Benjamini-Hochberg padj at alpha = {STAGE_DE_ALPHA} with pydeseq2's defaults; counts on the "
            "exported set (MOTCO's features, the primary reading) and on the abundance-passing set. The predeclared "
            "concordance rule compares the exported interaction count of factor_interaction with the primary "
            "analysis's delta rejection status. Before anything is fitted the export is checked against the "
            "pseudobulk (manifest hash, and the export's preprocessing reproduced from the counts) and the analysis "
            "against the export (hashes, MOTCO version, seed, permutations, CV selection). The primary analysis "
            "stays primary and is never rerun or written. Writes <out-dir>/<cohort>/ (<omic>/<model>__<coefficient>.csv, "
            "<omic>/<model>__fit.json, stage_de_summary.csv, stage_de.md, stage_de_manifest.json); resumable per "
            "model."
        ),
    )
    _add_cell_type_argument(sub_sd)
    sub_sd.add_argument(
        "--cohort",
        choices=list(EXPORT_COHORTS) + ["all"],
        default="all",
        help="Cohort(s) to run (default: %(default)s)",
    )
    sub_sd.add_argument(
        "--omic",
        choices=["rna", "atac", "both"],
        default="both",
        help=(
            "Omic(s) to fit; the summary covers every omic whose fits are complete in the cohort directory "
            "(default: %(default)s)"
        ),
    )
    sub_sd.add_argument(
        "--pseudobulk-dir",
        default=None,
        type=Path,
        help=(
            "Directory holding rna_counts.h5ad, atac_counts.h5ad, donor_table.csv and manifest.json; read, never "
            "written " + _dir_default("pseudobulk")
        ),
    )
    sub_sd.add_argument(
        "--input-dir",
        default=None,
        type=Path,
        help="Export root of the export step; read as the primary export, never written " + _dir_default("motco_input"),
    )
    sub_sd.add_argument(
        "--analysis-dir",
        default=None,
        type=Path,
        help="Analysis root of the analyze step; read as the primary analysis, never written " + _dir_default("motco_analysis"),
    )
    sub_sd.add_argument(
        "--out-dir",
        default=None,
        type=Path,
        help="Output root; one subdirectory per cohort " + _dir_default("stage_de"),
    )
    sub_sd.add_argument(
        "--n-cpus",
        type=int,
        default=1,
        help="pydeseq2 worker processes; runtime only, recorded in the manifest (default: %(default)s)",
    )
    sub_sd.add_argument(
        "--force",
        action="store_true",
        help="Refit the requested omics even when their fits are complete, replacing them and the cohort's summary",
    )

    # stage-de-limma command: the same models with limma-voom, beside the DESeq2 step
    joints = ", ".join(f"{j} ({m})" for m, (j, _) in STAGE_DE_LIMMA_JOINT_TESTS.items())
    sub_sl = sub.add_parser(
        "stage-de-limma",
        help="limma-voom fits of the stage-DE models beside the DESeq2 step (second predeclared method)",
        description=(
            "Per cohort and omic, on the same donors, raw counts, registered models (STAGE_DE_MODELS: "
            f"{models}), design matrices and abundance-passing features as `stage-de`: edgeR TMM with the export's "
            "trims (checked against the export manifest's factors to 1e-8) -> voomWithQualityWeights -> lmFit -> "
            "eBayes with limma's defaults (STAGE_DE_LIMMA_METHOD, a constant), run in R through the packaged script "
            "stage_de_limma.R (Rscript from the Nix shell, one BLAS thread). A moderated t per tested coefficient and "
            f"a joint moderated F per factor model ({joints}); Benjamini-Hochberg over all tested features at "
            f"alpha = {STAGE_DE_ALPHA}. The predeclared limma concordance cell (exported interaction_joint count "
            "against the primary analysis's delta rejection status) is reported beside the DESeq2 cell, which stays "
            "the reported one; method_comparison.csv compares the two methods per coefficient, union and joint test. "
            "Every check of `stage-de` runs first, and its outputs (--deseq2-dir) must be complete and made from the "
            "inputs now on disk; nothing is written under the inputs or the DESeq2 directory. Writes "
            "<out-dir>/<cohort>/ (<omic>/<model>__<label>.csv, <omic>/<model>__fit.json, stage_de_limma_summary.csv, "
            "method_comparison.csv, stage_de_limma.md, stage_de_limma_manifest.json); resumable per model."
        ),
    )
    _add_cell_type_argument(sub_sl)
    sub_sl.add_argument(
        "--cohort",
        choices=list(EXPORT_COHORTS) + ["all"],
        default="all",
        help="Cohort(s) to run (default: %(default)s)",
    )
    sub_sl.add_argument(
        "--omic",
        choices=["rna", "atac", "both"],
        default="both",
        help=(
            "Omic(s) to fit; the summary covers every omic whose fits are complete in the cohort directory "
            "(default: %(default)s)"
        ),
    )
    sub_sl.add_argument(
        "--pseudobulk-dir",
        default=None,
        type=Path,
        help=(
            "Directory holding rna_counts.h5ad, atac_counts.h5ad, donor_table.csv and manifest.json; read, never "
            "written " + _dir_default("pseudobulk")
        ),
    )
    sub_sl.add_argument(
        "--input-dir",
        default=None,
        type=Path,
        help="Export root of the export step; read as the primary export, never written " + _dir_default("motco_input"),
    )
    sub_sl.add_argument(
        "--analysis-dir",
        default=None,
        type=Path,
        help="Analysis root of the analyze step; read as the primary analysis, never written " + _dir_default("motco_analysis"),
    )
    sub_sl.add_argument(
        "--deseq2-dir",
        default=None,
        type=Path,
        help="Output root of `stage-de` (DESeq2); read and checked, never written " + _dir_default("stage_de"),
    )
    sub_sl.add_argument(
        "--out-dir",
        default=None,
        type=Path,
        help="Output root; one subdirectory per cohort " + _dir_default("stage_de_limma"),
    )
    sub_sl.add_argument(
        "--force",
        action="store_true",
        help="Refit the requested omics even when their fits are complete, replacing them and the cohort's summary",
    )

    # positive-control command: predeclared expectations of the control cell type against the reference
    control, reference = POSITIVE_CONTROL_CELL_TYPES["control"], POSITIVE_CONTROL_CELL_TYPES["reference"]
    sub_pc = sub.add_parser(
        "positive-control",
        help=f"Evaluate the predeclared positive-control expectations E1-E3 ({control} against {reference})",
        description=(
            f"Per cohort, from the generated files of both cell types only: E1 (stage signal) compares the number "
            f"of exported RNA features with a stage effect at FDR {STAGE_DE_ALPHA} between {control} and {reference} "
            "within each stage-DE method (DESeq2 factor_additive union, limma factor_additive stage_joint) and holds "
            "only under both (under one only: method_dependent, not holding); E2 (latent separation) compares the "
            "primary analyses' CV mean AUROC; E3 (axis identity) compares, along the first latent component of the "
            f"{control} primary analysis, |Spearman| with stage against the largest |Spearman| with a technical "
            "variable (technical table, TMM factors, library sizes) from axis-diagnostics. The control holds when "
            f"E1, E2 and E3 all hold in {POSITIVE_CONTROL_PRIMARY_COHORT}. No p-value, rejection status or trajectory "
            "statistic is read; nothing is tested. Every manifest read must record its root's cell type, and the "
            "stage-DE steps and the diagnostic must have been made from the analysis read for E2. Writes "
            "<out-dir>/<cohort>/ (positive_control.csv, positive_control.md, positive_control_manifest.json). With "
            "--basis latent_adjusted, E2 and E3 read each cell type's latent_adjusted arm of latent-sensitivity and "
            "the diagnostic made on it (E1 unchanged; the arm chain checked), writing <out-dir>/latent_adjusted/<cohort>/; "
            "the restated control is read as a diagnosis of the primary-basis control, never as adopting the arm."
        ),
    )
    sub_pc.add_argument(
        "--cohort",
        choices=list(EXPORT_COHORTS) + ["all"],
        default="all",
        help="Cohort(s) to evaluate (default: %(default)s)",
    )
    sub_pc.add_argument(
        "--basis",
        choices=list(POSITIVE_CONTROL_BASES),
        default=POSITIVE_CONTROL_DEFAULT_BASIS,
        help=(
            "The analysis and diagnostic E2 and E3 read (registry POSITIVE_CONTROL_BASES): primary, or latent_adjusted "
            "(the latent_adjusted arm of latent-sensitivity of each cell type and the diagnostic made on it; outputs "
            "under <out-dir>/latent_adjusted/<cohort>/) (default: %(default)s)"
        ),
    )
    sub_pc.add_argument(
        "--control-root",
        default=None,
        type=Path,
        help=f"Results root of the control cell type, {control} (default: {cell_type_dirs(control)['pseudobulk'].parent})",
    )
    sub_pc.add_argument(
        "--reference-root",
        default=None,
        type=Path,
        help=(
            f"Results root of the reference cell type, {reference}; read, never rerun "
            f"(default: {cell_type_dirs(reference)['pseudobulk'].parent})"
        ),
    )
    sub_pc.add_argument(
        "--out-dir",
        default=None,
        type=Path,
        help="Output root; one subdirectory per cohort (default: <control root>/positive_control)",
    )
    sub_pc.add_argument(
        "--force",
        action="store_true",
        help="Replace existing cohort directories",
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
        print(f"[pseudobulk] Building {args.cell_type} donor pseudobulks...")
        outputs = run_pseudobulk(
            omic=args.omic,
            cell_type=args.cell_type,
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
        outputs = run_covariates(
            omic=args.omic, cell_type=args.cell_type, data_root=args.data_root, out_dir=args.out_dir
        )
        for name in ("covariates", "manifest"):
            print(f"[covariates] Wrote {outputs[name]}")
        print("[covariates] Done.")
        return 0

    if args.command == "technical-covariates":
        print(f"[technical-covariates] Summarizing the {args.cell_type} nuclei of every donor object...")
        outputs = run_technical_covariates(
            args.omic,
            cell_type=args.cell_type,
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
            cell_type=args.cell_type,
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
            cell_type=args.cell_type,
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
        tables = [res["table"][columns] for per_cohort in results.values() for res in per_cohort.values()]
        if tables:
            table = pd.concat(tables, ignore_index=True)
            print("[preprocessing-sensitivity] Comparison table (primary = the analyze step's result, primary; "
                  "the variant rows are sensitivity):")
            print(table.to_string(index=False))
            for per_cohort in results.values():
                for res in per_cohort.values():
                    if res.get("adoption") is not None:
                        holds = bool(res["adoption"].iloc[-1]["holds"])
                        print(f"[preprocessing-sensitivity] {res['manifest']['variant']['name']} {res['manifest']['cohort']}: "
                              f"adoption check {'holds' if holds else 'does not hold'} (adoption_check.csv)")
        print("[preprocessing-sensitivity] Done.")
        return 0

    if args.command == "axis-diagnostics":
        print("[axis-diagnostics] Measuring the latent and block axes...")
        results = run_axis_diagnostics(
            args.cohort,
            cell_type=args.cell_type,
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

    if args.command == "latent-sensitivity":
        import pandas as pd

        print("[latent-sensitivity] Running the covariate-set and latent-basis sensitivity analysis...")
        results = run_latent_sensitivity(
            args.cohort,
            cell_type=args.cell_type,
            input_dir=args.input_dir,
            analysis_dir=args.analysis_dir,
            pseudobulk_dir=args.pseudobulk_dir,
            out_dir=args.out_dir,
            n_jobs=args.n_jobs,
            force=args.force,
        )
        columns = ["arm", "cohort", "conditioning", "rank"]
        for s in ANALYSIS_STATISTICS:
            columns += [s, f"{s}_p", f"{s}_rejected"]
        table = pd.concat([res["table"][columns] for res in results.values()], ignore_index=True)
        print("[latent-sensitivity] Comparison table (primary = the analyze step's result, primary; the arms are "
              "sensitivity):")
        print(table.to_string(index=False))
        for cohort, res in results.items():
            flags = res["manifest"]["mediator_flags"]
            print(f"[latent-sensitivity] {cohort}: mediator check ({res['dir'] / 'mediator_check.csv'}): "
                  + ("no flag" if not flags else ", ".join(
                      f"{f['covariate']} {f['subset']} rho {f['spearman']:.3f}" for f in flags)))
        print("[latent-sensitivity] Done.")
        return 0

    if args.command == "stage-de":
        import pandas as pd

        print("[stage-de] Running the per-feature stage and sex x stage DESeq2 models...")
        results = run_stage_de(
            args.cohort,
            omic=args.omic,
            cell_type=args.cell_type,
            pseudobulk_dir=args.pseudobulk_dir,
            input_dir=args.input_dir,
            analysis_dir=args.analysis_dir,
            out_dir=args.out_dir,
            n_cpus=args.n_cpus,
            force=args.force,
        )
        columns = ["cohort", "omic", "model", "coefficient", "feature_set", "n_tested", "n_significant", "min_padj"]
        table = pd.concat([res["summary"][columns] for res in results.values()], ignore_index=True)
        print(f"[stage-de] Counts at padj <= {STAGE_DE_ALPHA} (union = either coefficient, no joint error-rate claim):")
        print(table.to_string(index=False))
        for cohort, res in results.items():
            for omic, cell in res["concordance"].items():
                print(
                    f"[stage-de] {cohort} {omic}: concordance {cell['cell']} (interaction detected: "
                    f"{'yes' if cell['interaction_detected'] else 'no'}, {cell['n_interaction_exported']} exported "
                    f"features; primary delta rejected: {'yes' if cell['delta_rejected'] else 'no'})"
                )
        print("[stage-de] Done.")
        return 0

    if args.command == "stage-de-limma":
        import pandas as pd

        print("[stage-de-limma] Running the limma-voom stage and sex x stage models beside the DESeq2 step...")
        results = run_stage_de_limma(
            args.cohort,
            omic=args.omic,
            cell_type=args.cell_type,
            pseudobulk_dir=args.pseudobulk_dir,
            input_dir=args.input_dir,
            analysis_dir=args.analysis_dir,
            deseq2_dir=args.deseq2_dir,
            out_dir=args.out_dir,
            force=args.force,
        )
        columns = ["cohort", "omic", "model", "coefficient", "kind", "feature_set", "n_tested", "n_significant", "min_padj"]
        table = pd.concat([res["summary"][columns] for res in results.values()], ignore_index=True)
        print(
            f"[stage-de-limma] Counts at adj. p <= {STAGE_DE_ALPHA} (joint_f = moderated F over both coefficients; "
            "union = either coefficient, no joint error-rate claim):"
        )
        print(table.to_string(index=False))
        for cohort, res in results.items():
            for omic, cell in res["concordance"].items():
                print(
                    f"[stage-de-limma] {cohort} {omic}: limma cell {cell['limma']['cell']} "
                    f"({cell['limma']['n_interaction_exported']} exported features on the joint F), DESeq2 cell "
                    f"{cell['deseq2']['cell']} ({cell['deseq2']['n_interaction_exported']} on the union); primary delta "
                    f"rejected: {'yes' if cell['delta_rejected'] else 'no'}; methods agree: "
                    f"{'yes' if cell['methods_agree'] else 'no'}"
                )
        print("[stage-de-limma] Done.")
        return 0

    if args.command == "positive-control":
        print("[positive-control] Evaluating the predeclared positive-control expectations...")
        results = run_positive_control(
            args.cohort,
            basis=args.basis,
            control_root=args.control_root,
            reference_root=args.reference_root,
            out_dir=args.out_dir,
            force=args.force,
        )
        for cohort, res in results.items():
            o = res["outcomes"]
            print(
                f"[positive-control] {cohort}{' (primary cohort)' if res['manifest']['primary_cohort'] else ''}"
                f"{'' if args.basis == POSITIVE_CONTROL_DEFAULT_BASIS else f' (basis {args.basis})'}: "
                f"E1 {o['E1']}, E2 {o['E2']}, E3 {o['E3']}; all three hold: {'yes' if o['all_hold'] else 'no'}"
            )
        print("[positive-control] Done.")
        return 0

    parser.error("Unknown command")
    return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())

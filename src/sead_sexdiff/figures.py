"""Case-study figures: ``sead-sexdiff figures``.

Draws the registered figure set (``FIGURES``) of the SEA-AD MTG MOTCO case
study from the generated outputs of both cell types only (the default
directories of ``CELL_TYPES``, plus the astrocyte preprocessing variants),
after a consistency check of every manifest those outputs carry, with
captions generated from the files and a manifest (``figures_manifest.json``)
that ties every image to the files it was made from.

``results/figures/`` is a working copy: a run does nothing when every figure
is fresh and redraws the whole set otherwise; ``--check`` lists the stale
figures; ``--report <stem>`` freezes the set under ``reports/figures/<stem>/``
for a committed report. The figures are a display step: they test nothing and
are an input of no other step.

Plotting lives here rather than in :mod:`sead_sexdiff.data` (the exception to
"all logic in data.py", decided 2026-11-02); this module imports from
``data.py`` and never the other way round.
"""

from __future__ import annotations

import hashlib
import inspect
import io
import json
import os
import re
import shutil
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path, PurePosixPath

import matplotlib
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.collections import PathCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

from .data import (
    ANALYSIS_ALPHA,
    ANALYSIS_BLOCKS,
    ANALYSIS_FORCED_RANKS,
    ANALYSIS_INTERACTION_SCAN,
    ANALYSIS_PERMUTATIONS,
    ANALYSIS_RANK_FLOOR,
    ANALYSIS_RESIDUAL_CORRECTION,
    ANALYSIS_SEED,
    ANALYSIS_SPAN_SPLIT_TEST,
    ANALYSIS_STATISTICS,
    CELL_TYPES,
    COHORT_THRESHOLDS,
    DEFAULT_CELL_TYPE,
    DIAGNOSTIC_CORRELATION_THRESHOLD,
    DIAGNOSTIC_OUTPUTS,
    DIAGNOSTIC_PRIMARY_LABEL,
    DONOR_EXCLUSION_ANALYSIS_DIR,
    DONOR_EXCLUSION_COMPARISON_DIR,
    DONOR_EXCLUSION_MANIFEST,
    DONOR_EXCLUSION_VARIANTS,
    EXPORT_COHORTS,
    INTERACTION_SCAN_POWER_STATEMENT,
    INTERACTION_SCAN_STATISTIC,
    INTERACTION_SCAN_TABLE,
    LATENT_SENSITIVITY_ANALYSIS_DIR,
    LATENT_SENSITIVITY_ARMS,
    LATENT_SENSITIVITY_COMPARISON_DIR,
    LATENT_SENSITIVITY_MANIFEST,
    LATENT_SENSITIVITY_SPLIT_UNAVAILABLE,
    POSITIVE_CONTROL_BASES,
    POSITIVE_CONTROL_CELL_TYPES,
    POSITIVE_CONTROL_DEFAULT_BASIS,
    POSITIVE_CONTROL_MANIFEST,
    POSITIVE_CONTROL_PRIMARY_COHORT,
    POSITIVE_CONTROL_TABLE,
    POSITIVE_CONTROL_TECHNICAL_SOURCES,
    PREPROCESSING_ANALYSIS_DIR,
    PREPROCESSING_COMPARISON_DIR,
    PREPROCESSING_INPUT_DIR,
    PREPROCESSING_SENSITIVITY_MANIFEST,
    PREPROCESSING_VARIANTS,
    RANK_SENSITIVITY_DIR_TEMPLATE,
    RANK_SENSITIVITY_MANIFEST,
    SPAN_SPLIT_CONCORDANCE_PART,
    SPAN_SPLIT_CONFIRMATORY_PART,
    SPAN_SPLIT_PARTS,
    SPAN_SPLIT_ROLES,
    SPAN_SPLIT_STATISTIC,
    STAGE_DE_ALPHA,
    STAGE_DE_LIMMA_COMPARISON,
    STAGE_DE_LIMMA_JOINT_TESTS,
    STAGE_DE_LIMMA_MANIFEST,
    STAGE_DE_LIMMA_SUMMARY,
    STAGE_DE_LIMMA_TABLE_TEMPLATE,
    STAGE_DE_MANIFEST,
    STAGE_DE_MODELS,
    STAGE_DE_SUMMARY,
    STAGE_DE_TABLE_TEMPLATE,
    STAGE_DE_UNION,
    STAGE_LABELS,
    TECHNICAL_GENE_SETS,
    TECHNICAL_OUTPUTS,
    _git_revision,
    _interaction_scan_status,
    _sha256_file,
    _span_split_status,
    cell_type_dirs,
    check_manifest_cell_type,
    run_preprocessing_sensitivity,
    spearman,
)

# ---------------------------------------------------------------------------
# Display constants (predeclared in the proposal, never flags, never chosen
# from results; recorded in the manifest). Flag and ask before changing.
# ---------------------------------------------------------------------------

# The cohort of the single-cohort figures (03, 07, 08, 10): the primary cohort;
# every figure that draws both cohorts draws cohort_50 beside it
FIGURE_COHORT = "cohort_30"
# Figure 10's volcano and method comparison: the contrast between the extreme
# stages on RNA; the caption gives the counts of the coefficients not drawn
FIGURE_STAGE_DE: dict[str, str] = {"omic": "rna", "model": "factor_additive", "coefficient": "stage_2_vs_0"}
# Figure 08's technical variable: E3's rule (POSITIVE_CONTROL_EXPECTATIONS),
# applied to every cell type
FIGURE_FIRST_AXIS_RULE = (
    "the variable with the largest |Spearman| with pls_0 of the primary analysis among the axis_correlations.csv rows "
    "whose source is in POSITIVE_CONTROL_TECHNICAL_SOURCES (the technical table's columns, the export's TMM factors "
    "and library sizes), as E3 picks it"
)
# rank-sensitivity and preprocessing-sensitivity take no --cell-type and run on
# astrocytes only; the preprocessing variants live outside cell_type_dirs
FIGURE_SENSITIVITY_CELL_TYPE = DEFAULT_CELL_TYPE
FIGURE_PREPROCESSING_ROOT = Path("results/motco_sensitivity")
# Scatter layers with more offsets than this are embedded as images in the
# PDF and SVG files (at FIGURE_DPI); axes, text, lines and smaller layers stay vector
FIGURE_RASTER_MIN_POINTS = 500
# Figure 12 (a): the smallest adjusted p-values drawn per block and primary
FIGURE_SCAN_TOP_N = 50
FIGURE_DPI = 200
FIGURE_FORMATS = ("png", "pdf", "svg")
FIGURE_SVG_HASHSALT = "sead-sexdiff"
FIGURE_FONT = "DejaVu Sans"
# Jitter of the strip plots (numpy default_rng seeds), one per figure
FIGURE_JITTER_SEEDS: dict[str, int] = {"components_by_stage": 0, "first_axis": 1}
# Reproduction tolerances (spec: Recomputed values reproduce their source)
FIGURE_STAGE_MEANS_ATOL = 1e-8
FIGURE_SPEARMAN_ATOL = 1e-9
# Palette: sex and cell-type pairs validated all-pairs (dataviz reference
# instance); the aqua cell-type colour is below 3:1 on white, so every
# cell-type mark carries a legend entry or a direct label
FIGURE_PALETTE: dict = {
    "ink": "#0b0b0b",
    "ink2": "#52514e",
    "muted": "#898781",
    "grid": "#e1e0d9",
    "axis": "#c3c2b7",
    "neutral": "#f0efec",
    "sex": {"Female": "#2a78d6", "Male": "#eb6834"},
    "cell_type": {"astrocyte": "#1baf7a", "microglia": "#6250d6"},
    "good": "#0ca30c",
    "critical": "#d03b3b",
    "alpha_band": "#f4b6b5",
    "diverging": ["#104281", "#3987e5", "#b7d3f6", "#f0efec", "#f4b6b5", "#e34948", "#8f1f1f"],
}

# Paths and names
FIGURE_RESULTS_ROOT = "results"
FIGURE_SUPERSEDED_DIR = "superseded"
FIGURE_OUT_DIR = Path("results/figures")
FIGURE_REPORTS_DIR = Path("reports")
FIGURE_REPORT_SUBDIR = "figures"
FIGURE_REPORT_STEM_RE = re.compile(r"^[a-z0-9][a-z0-9.-]*$")
FIGURE_FROZEN_FORMAT = "svg"
FIGURE_MANIFEST = "figures_manifest.json"
FIGURE_CAPTIONS = "figures.md"
_RESULT = "result.json"
_ANALYSIS_MANIFEST = "analysis_manifest.json"
_PSEUDOBULK_MANIFEST = "manifest.json"
_EXPORT_MANIFEST = "manifest.json"

# Import-time ties
if FIGURE_COHORT != POSITIVE_CONTROL_PRIMARY_COHORT:
    raise AssertionError(
        f"FIGURE_COHORT ({FIGURE_COHORT!r}) must be the primary cohort, POSITIVE_CONTROL_PRIMARY_COHORT "
        f"({POSITIVE_CONTROL_PRIMARY_COHORT!r})"
    )
if FIGURE_STAGE_DE["model"] not in STAGE_DE_MODELS or FIGURE_STAGE_DE["coefficient"] not in [
    label for label, _ in STAGE_DE_MODELS[FIGURE_STAGE_DE["model"]]["tested"]
]:
    raise AssertionError(
        f"FIGURE_STAGE_DE names model {FIGURE_STAGE_DE['model']!r} and coefficient {FIGURE_STAGE_DE['coefficient']!r}, "
        "which STAGE_DE_MODELS does not test"
    )
if FIGURE_STAGE_DE["model"] not in STAGE_DE_LIMMA_JOINT_TESTS or FIGURE_STAGE_DE["omic"] not in ANALYSIS_BLOCKS:
    raise AssertionError("FIGURE_STAGE_DE must name a factor model with a limma joint test and an analysis block")
if inspect.signature(run_preprocessing_sensitivity).parameters["out_dir"].default != str(FIGURE_PREPROCESSING_ROOT):
    raise AssertionError("FIGURE_PREPROCESSING_ROOT must be run_preprocessing_sensitivity's default out_dir")

# The reading rules the captions restate (AGENTS.md, the binding MOTCO rules of
# ../MOTCO/docs/roadmap.md, and each step's registered reading). A rule
# changed in AGENTS.md by a later change changes its text here; flag and ask
# before changing.
FIGURE_READING_RULES: dict[str, str] = {
    "display_only": (
        "The figure displays outputs already reported; it tests nothing, is an input of no step and changes no "
        "reading of the reports that made the outputs."
    ),
    "non_rejection": (
        "angle and shape non-rejections are uninformative; a delta rejection is a difference, not specifically of "
        "size; a delta non-rejection is weak evidence (delta is conservative under the df correction)."
    ),
    "motco_binding": (
        "MOTCO results are reported under the binding rules: the CV-selected rank beside angle; angle and shape "
        "non-rejections uninformative; a delta rejection a difference, not specifically of size; the covariate "
        "block and the feature filtering beside the results. The rank floor defines angle and shape in at least two "
        "dimensions without making them powerful (two dimensions, not a second axis with signal), so each "
        "analysis's configuration eigengap and null summary belong beside angle."
    ),
    "df_size": (
        "Under the df residual correction the p-values are size-controlled at the configuration MOTCO measured (six "
        "covariates, both cohort_30 tables), with delta conservative there, so a delta non-rejection stays weak "
        "evidence; the size is inferred from the mechanism (n and the reduced rank), not measured, for cohort_50, the "
        "eight-covariate latent arms, the exclusion variants, the forced ranks and the preprocessing variants."
    ),
    "sensitivity_dependence": (
        "Every analysis other than the primary is sensitivity only: a rejection status that differs from the "
        "primary's is rank, preprocessing, covariate-set, basis or cohort dependence, not a result of the variant; a "
        "CV-selected rank that differs is reported beside angle, not as better or worse; the non-rejection rules hold "
        "in every variant; the primary analysis stays primary."
    ),
    "forced_rank": (
        "The forced-rank sweep is sensitivity only: the CV-selected analysis stays primary, and a rejection status "
        "that differs between ranks is rank dependence, not a result at that rank. MOTCO's forced-rank operating "
        "characteristics are carried as limitations (Type I checked at forced 4 only, unadjusted; shape power loss "
        "and angle leakage at forced 2(K - 1); a wider angle null above n_stages - 1), and the adjusted analysis has "
        "no Type I check at any forced rank."
    ),
    "latent_adjusted_not_read": (
        "The latent_adjusted arm's angle and shape p-values are reported but not read (MOTCO's diagnosis "
        "residualized_basis: the stage latent adjustment is liberal on both at this configuration, an excess of the "
        "basis the residual correction does not address, so neither a rejection nor a non-rejection is "
        "interpretable); its delta is read under the rules above."
    ),
    "measurement_space": (
        "An axis correlated with a technical variable is a measurement-space finding about the analysis, not a "
        "reinterpretation of its statistics; the primary analysis stays primary with its non-rejection rules."
    ),
    "diagnostic": (
        "The diagnostics read sex and stage, so they are diagnostics and never a filter: their outputs feed no "
        "export, analysis or sensitivity run, and no preprocessing variant is derived from a diagnostic value."
    ),
    "stage_de_counts": (
        "The stage counts describe the stage signal in the pseudobulks; the MOTCO non-rejection rules are untouched. "
        "A zero count is none detected at this n; the union counts carry no joint error-rate claim and limma's joint "
        "F claim is approximate; a call that differs between DESeq2 and limma is method dependence, not a result of "
        "either, and agreement is not replication (same donors and counts). The FDR of padj <= alpha is approximate "
        "and Cook's handling is inactive with continuous covariates; the covariates, the feature set and each "
        "method's normalization (DESeq2's median of ratios; limma-voom's log-CPM on the export's TMM) are stated "
        "beside the counts."
    ),
    "positive_control": (
        "The positive control is about the stage signal, not the sex x stage hypothesis: it holds when E1, E2 and E3 "
        "all hold in cohort_30, and the primary-basis cohort_30 row is its outcome; cohort_50 is descriptive. If it "
        "holds, the reference's null is not a pipeline failure; if not, the latent construction is in question for "
        "any cell type. E1 does "
        "not depend on the basis; the latent_adjusted basis is a diagnosis of the primary-basis control, never a "
        "replacement (the arm is not promoted). No expectation reads a p-value."
    ),
    "span_split": (
        "MOTCO's split interaction test splits the group x stage interaction energy of the analysed features into the "
        "fitted stage-PLS span, its orthogonal complement and all of it. Only the complement is confirmatory: its "
        "pseudo-F p-value is read at alpha = 0.05, unadjusted (no adjustment across delta, angle, shape, the "
        "complement and the interaction scan, so each primary carries five tests read at alpha), and each cell type's "
        "cohort_30 primary gives "
        "the result; a complement status that differs from the primary's elsewhere is rank, preprocessing, "
        "covariate-set or cohort dependence. The span part is a concordance cell against any of delta, angle and "
        "shape rejected in the same analysis, a discordance a finding and never a result; total, the in-span share "
        "and the interaction-SS p-values are reported and not read. A complement rejection is a sex x stage "
        "difference with interaction energy outside the fitted span, not a different molecular program: it may hold "
        "leakage from an in-span difference and what the fitted span missed of the stage program, and it carries no "
        "direction, feature list or attribution. A complement non-rejection is uninformative: no power was shown at "
        "SEA-AD's stage signal, nor for sparse differences at any signal. The in-span share grows with the fitted rank "
        "and is read against no null here. Size was measured on real noise only for the astrocyte cohort_30 export at "
        "this configuration (complement 0.012); for cohort_50, microglia, the eight-covariate arm, the exclusion "
        "variants, the forced ranks and the preprocessing variants it is inferred. The test is not available with a "
        "latent adjustment."
    ),
    "interaction_scan": (
        "MOTCO's per-feature interaction scan tests each analysed feature's group x stage interaction F with the "
        "covariate-adjusted models, combined by the maximum F. Its global p-value is the one confirmatory reading, "
        "read at alpha = 0.05, unadjusted (no adjustment across delta, angle, shape, the complement and the scan, so "
        "each primary carries five tests read at alpha); a hit is a feature with step-down max-T adjusted p <= 0.05 "
        "(familywise over the scanned features), and a rejection comes with at least one hit. Each cell type's "
        "cohort_30 primary gives the result; a scan status that differs elsewhere is preprocessing, covariate-set or "
        "cohort dependence, and the forced ranks and the latent-adjusted arms equal their reference by construction "
        "(the scan reads neither the latent space nor the bootstrap). A hit names a feature whose covariate-adjusted "
        "sex difference changes across stages: a detection, not a driver (correlated neighbours can be hits too), with "
        "no geometry, trajectory or attribution, changing no reading of delta, angle, shape or the complement; the "
        "per-feature F values, the per-stage differences and the null-maximum shares are reported and not read. A "
        "non-rejection of the astrocyte cohort_30 primary carries MOTCO's power on that export (a 20-feature ramp in "
        "one sex detected in 1.000 of replicates at 3 within-stage SDs, 0.692 at 2), which excludes neither smaller "
        "effects, nor effects spread thinly over many features, nor other shapes; every other non-rejection is "
        "uninformative. Size was measured on real noise only for the astrocyte cohort_30 export at six covariates "
        "(global 0.042, familywise 0.042); elsewhere it is inferred. The null maxima lean on ATAC there (0.618 of them "
        "with half the features), and the F is the plain, unmoderated per-feature F. The scan concordance cells beside "
        "the stage DE are findings when discordant, never results."
    ),
}


# ---------------------------------------------------------------------------
# The recording reader (design D3)
# ---------------------------------------------------------------------------

def _stat_key(path: Path) -> tuple[int, int]:
    st = path.stat()
    return st.st_size, st.st_mtime_ns


class FigureReader:
    """The one reader of every file the figures use (design D3).

    Every read records ``(path, sha256)`` under the figure being drawn, so a
    figure's input list is the reader's log and cannot omit a file the figure
    used. Each file is hashed once per run (keyed by size and mtime, which
    also catches a file changed during the run). Absolute paths, paths
    outside ``results/`` and paths under ``results/superseded/`` are refused.
    """

    def __init__(self) -> None:
        self._sha: dict[str, tuple[tuple[int, int], str]] = {}
        self._bytes: dict[str, bytes] = {}
        self.n_hashed = 0
        self.figure: str | None = None
        self.inputs: dict[str, dict[str, str]] = {}
        self.manifests: dict[str, list[str]] = {}
        self.analyses: dict[str, dict[str, dict]] = {}
        self.findings: dict[str, list[str]] = {}

    @staticmethod
    def key(path: str | Path) -> str:
        """The repository-relative POSIX path of a readable output; raise for any other path."""
        p = PurePosixPath(Path(path).as_posix())
        if p.is_absolute():
            raise ValueError(f"{path}: absolute path; the figures read repository-relative paths under results/ only")
        if not p.parts or p.parts[0] != FIGURE_RESULTS_ROOT or ".." in p.parts:
            raise ValueError(f"{path}: outside results/; the figures read generated outputs under results/ only")
        if len(p.parts) > 1 and p.parts[1] == FIGURE_SUPERSEDED_DIR:
            raise ValueError(f"{path}: under results/superseded/; superseded outputs are never read")
        return str(p)

    def begin(self, figure: str) -> None:
        self.figure = figure
        self.inputs[figure] = {}
        self.manifests[figure] = []
        self.analyses[figure] = {}
        self.findings[figure] = []

    def end(self) -> None:
        self.figure = None

    def _hash(self, key: str, stat: tuple[int, int], digest: str | None = None, data: bytes | None = None) -> str:
        cached = self._sha.get(key)
        if cached is not None:
            if cached[0] != stat:
                raise ValueError(f"{key} changed during the run (size or mtime differs from when it was hashed)")
            return cached[1]
        if digest is None:
            digest = hashlib.sha256(data).hexdigest() if data is not None else _sha256_file(Path(key))
        self._sha[key] = (stat, digest)
        self.n_hashed += 1
        return digest

    def sha256(self, path: str | Path) -> str:
        """The file's sha256 (hashed once per run); not logged under a figure."""
        key = self.key(path)
        p = Path(key)
        if not p.is_file():
            raise FileNotFoundError(f"{key}: missing")
        return self._hash(key, _stat_key(p))

    def read_bytes(self, path: str | Path) -> bytes:
        """The file's bytes, logged with their sha256 under the figure being drawn."""
        key = self.key(path)
        if key not in self._bytes:
            p = Path(key)
            if not p.is_file():
                raise FileNotFoundError(f"{key}: missing; the figures read generated outputs only")
            before = _stat_key(p)
            data = p.read_bytes()
            if _stat_key(p) != before or len(data) != before[0]:
                raise ValueError(f"{key} changed while it was read")
            self._hash(key, before, data=data)
            self._bytes[key] = data
        if self.figure is not None:
            self.inputs[self.figure][key] = self._sha[key][1]
        return self._bytes[key]

    def csv(self, path: str | Path, **kwargs) -> pd.DataFrame:
        return pd.read_csv(io.BytesIO(self.read_bytes(path)), **kwargs)

    def json(self, path: str | Path) -> dict:
        return json.loads(self.read_bytes(path))

    def manifest(self, path: str | Path, cell_type: str) -> dict:
        """A step's manifest, its cell type checked (``check_manifest_cell_type``) and recorded as a manifest read."""
        key = self.key(path)
        manifest = self.json(key)
        findings = check_manifest_cell_type(manifest, key, cell_type)
        if self.figure is not None:
            if key not in self.manifests[self.figure]:
                self.manifests[self.figure].append(key)
            for msg in findings:
                if msg not in self.findings[self.figure]:
                    self.findings[self.figure].append(msg)
        return manifest

    def result(self, path: str | Path) -> dict:
        """A MOTCO ``result.json``, recorded with its MOTCO version and measurement under the figure."""
        key = self.key(path)
        record = self.json(key)
        if self.figure is not None:
            self.analyses[self.figure][key] = _measurement(record)
        return record


def _measurement(record: dict) -> dict:
    params = record.get("params", {})
    rc = record.get("measurement", {}).get("residual_correction")
    return {
        "motco": record.get("software", {}).get("motco"),
        "residual_correction": params.get("residual_correction"),
        "kappa": rc.get("kappa") if isinstance(rc, dict) else None,
        "min_components": params.get("min_components"),
        "forced_components": params.get("forced_components"),
        "component_selection": record.get("latent", {}).get("component_selection"),
        "selected_lv": record.get("latent", {}).get("selected_lv"),
        "permutations": params.get("permutations"),
        "seed": params.get("seed"),
        "span_split_test": bool(params.get("span_split_test", False)),
        "span_split_present": isinstance(record.get("span_split"), dict),
        "interaction_scan": bool(params.get("interaction_scan", False)),
        "interaction_scan_present": isinstance(record.get("interaction_scan"), dict),
    }


def _manifest_entries(obj) -> list[tuple[str, str]]:
    """Every nested ``{"path", "sha256"}`` entry of a manifest, in document order."""
    out: list[tuple[str, str]] = []
    if isinstance(obj, dict):
        if isinstance(obj.get("path"), str) and isinstance(obj.get("sha256"), str):
            out.append((obj["path"], obj["sha256"]))
        for value in obj.values():
            out.extend(_manifest_entries(value))
    elif isinstance(obj, list):
        for value in obj:
            out.extend(_manifest_entries(value))
    return out


def _check_made_from(manifest: dict, manifest_path: str | Path, result_path: str | Path) -> None:
    """Raise unless the manifest records ``result_path`` (the chain: the output was made from the analysis drawn)."""
    key = FigureReader.key(result_path)
    recorded = {PurePosixPath(Path(p).as_posix()).as_posix() for p, _ in _manifest_entries(manifest)}
    if key not in recorded:
        raise ValueError(f"{manifest_path} does not record {key}: it was not made from the analysis drawn beside it")


# ---------------------------------------------------------------------------
# The registered analyses and manifests, and the consistency pre-pass (D4)
# ---------------------------------------------------------------------------

def registered_analyses() -> list[dict]:
    """Every MOTCO analysis of the registries, one dict per (analysis, cohort), in drawing order.

    Keys: ``cell_type``, ``family`` (primary, forced rank, preprocessing,
    latent, exclusion), ``label``, ``cohort``, ``dir``, ``forced`` (the forced
    rank or None) and ``command`` (the step that writes it).
    """
    rows = []
    for cell_type in CELL_TYPES:
        d = cell_type_dirs(cell_type)
        for cohort in EXPORT_COHORTS:
            def add(family: str, label: str, directory: Path, command: str, forced: int | None = None) -> None:
                rows.append({"cell_type": cell_type, "family": family, "label": label, "cohort": cohort,
                             "dir": directory, "forced": forced, "command": command})

            add("primary", "primary", d["motco_analysis"] / cohort, "analyze")
            if cell_type == FIGURE_SENSITIVITY_CELL_TYPE:
                for rank in ANALYSIS_FORCED_RANKS:
                    add("forced rank", f"forced rank {rank}",
                        d["motco_rank_sensitivity"] / cohort / RANK_SENSITIVITY_DIR_TEMPLATE.format(rank),
                        "rank-sensitivity", forced=rank)
                for variant in PREPROCESSING_VARIANTS:
                    add("preprocessing", variant, FIGURE_PREPROCESSING_ROOT / variant / PREPROCESSING_ANALYSIS_DIR / cohort,
                        "preprocessing-sensitivity")
            for arm in LATENT_SENSITIVITY_ARMS:
                add("latent", arm, d["motco_latent_sensitivity"] / arm / LATENT_SENSITIVITY_ANALYSIS_DIR / cohort,
                    "latent-sensitivity")
            for variant in DONOR_EXCLUSION_VARIANTS:
                add("exclusion", variant,
                    d["motco_exclusion_sensitivity"] / variant / DONOR_EXCLUSION_ANALYSIS_DIR / cohort,
                    "exclusion-sensitivity")
    return rows


def _positive_control_dir(basis: str, cohort: str) -> Path:
    root = cell_type_dirs(POSITIVE_CONTROL_CELL_TYPES["control"])["positive_control"]
    return root / cohort if basis == POSITIVE_CONTROL_DEFAULT_BASIS else root / basis / cohort


def registered_manifests() -> list[tuple[Path, str, str]]:
    """The manifests the registries say exist: ``(path, cell type of its root, step that writes it)``."""
    out: list[tuple[Path, str, str]] = []
    for cell_type in CELL_TYPES:
        d = cell_type_dirs(cell_type)
        out.append((d["pseudobulk"] / _PSEUDOBULK_MANIFEST, cell_type, "pseudobulk"))
        out.append((d["pseudobulk"] / TECHNICAL_OUTPUTS["manifest"], cell_type, "technical-covariates"))
        for cohort in EXPORT_COHORTS:
            out.append((d["motco_input"] / cohort / _EXPORT_MANIFEST, cell_type, "export"))
            out.append((d["motco_diagnostics"] / cohort / DIAGNOSTIC_OUTPUTS["manifest"], cell_type, "axis-diagnostics"))
            out.append((d["stage_de"] / cohort / STAGE_DE_MANIFEST, cell_type, "stage-de"))
            out.append((d["stage_de_limma"] / cohort / STAGE_DE_LIMMA_MANIFEST, cell_type, "stage-de-limma"))
            out.append((d["motco_latent_sensitivity"] / LATENT_SENSITIVITY_COMPARISON_DIR / cohort
                        / LATENT_SENSITIVITY_MANIFEST, cell_type, "latent-sensitivity"))
            for variant in DONOR_EXCLUSION_VARIANTS:
                out.append((d["motco_exclusion_sensitivity"] / variant / DONOR_EXCLUSION_COMPARISON_DIR / cohort
                            / DONOR_EXCLUSION_MANIFEST, cell_type, "exclusion-sensitivity"))
            if cell_type == FIGURE_SENSITIVITY_CELL_TYPE:
                out.append((d["motco_rank_sensitivity"] / cohort / RANK_SENSITIVITY_MANIFEST, cell_type,
                            "rank-sensitivity"))
                for variant in PREPROCESSING_VARIANTS:
                    root = FIGURE_PREPROCESSING_ROOT / variant
                    out.append((root / PREPROCESSING_INPUT_DIR / cohort / _EXPORT_MANIFEST, cell_type,
                                "preprocessing-sensitivity"))
                    out.append((root / PREPROCESSING_COMPARISON_DIR / cohort / PREPROCESSING_SENSITIVITY_MANIFEST,
                                cell_type, "preprocessing-sensitivity"))
            if cell_type == POSITIVE_CONTROL_CELL_TYPES["control"]:
                for basis in POSITIVE_CONTROL_BASES:
                    out.append((_positive_control_dir(basis, cohort) / POSITIVE_CONTROL_MANIFEST, cell_type,
                                "positive-control"))
    for a in registered_analyses():
        out.append((a["dir"] / _ANALYSIS_MANIFEST, a["cell_type"], a["command"]))
    return out


def _registered_roots() -> list[tuple[Path, str]]:
    """The directories whose manifests are all checked: every cell type's default directories, the variants' root."""
    roots = [(directory, cell_type) for cell_type in CELL_TYPES for directory in cell_type_dirs(cell_type).values()]
    roots.append((FIGURE_PREPROCESSING_ROOT, FIGURE_SENSITIVITY_CELL_TYPE))
    return roots


def check_consistency(reader: FigureReader) -> dict:
    """The pre-pass before anything is drawn (design D4); raise listing every failure, nothing written.

    Checks, for the manifests of the registries (each must exist) and every
    other ``*manifest.json`` under the registered directories: the cell type
    of the root it was read from, and that every recorded ``{"path",
    "sha256"}`` entry exists under ``results/`` outside ``superseded/`` and
    hashes as recorded. Checks every registered analysis's ``result.json``
    for one MOTCO version (the astrocyte primary's) and the registered
    measurement (``ANALYSIS_RESIDUAL_CORRECTION``, ``ANALYSIS_RANK_FLOOR``,
    the forced rank under rank-sensitivity and none elsewhere,
    ``ANALYSIS_PERMUTATIONS``, ``ANALYSIS_SEED``, and the split test
    ``ANALYSIS_SPAN_SPLIT_TEST`` with its ``span_split`` block, except in the
    arms of ``LATENT_SENSITIVITY_SPLIT_UNAVAILABLE``, which must carry
    neither).

    Returns
    -------
    dict
        Counts (manifests from the registries and found, entries, distinct
        files, analyses), the MOTCO version and measurement, and the
        cell-type findings of pre-registry manifests.
    """
    analyses = registered_analyses()
    required = registered_manifests()
    missing = [f"{path} (written by `sead-sexdiff {command}`)" for path, _, command in required if not path.is_file()]
    missing += [f"{a['dir'] / _RESULT} (written by `sead-sexdiff {a['command']}`)" for a in analyses
                if not (a["dir"] / _RESULT).is_file()]
    if missing:
        raise ValueError(
            "registered output(s) missing; the registries say they exist, so the figures wait for them (nothing was "
            "written): " + "; ".join(missing)
        )
    manifests: dict[str, str] = {}
    for path, cell_type, _ in required:
        manifests.setdefault(FigureReader.key(path), cell_type)
    n_registered = len(manifests)
    for root, cell_type in _registered_roots():
        if root.is_dir():
            for path in sorted(root.rglob("*manifest.json")):
                manifests.setdefault(FigureReader.key(path), cell_type)

    errors: list[str] = []
    findings: list[str] = []
    n_entries = 0
    files: set[str] = set()
    for path, cell_type in manifests.items():
        manifest = reader.json(path)
        try:
            findings.extend(check_manifest_cell_type(manifest, path, cell_type))
        except ValueError as exc:
            errors.append(str(exc))
        for recorded, sha in _manifest_entries(manifest):
            n_entries += 1
            try:
                key = FigureReader.key(recorded)
            except ValueError as exc:
                errors.append(f"{path} records {exc}")
                continue
            files.add(key)
            if not Path(key).is_file():
                errors.append(f"{path} records {key} (sha256 {sha}), which is missing")
                continue
            now = reader.sha256(key)
            if now != sha:
                errors.append(f"{path} records {key} with sha256 {sha}; it now hashes {now}")

    reference = analyses[0]
    reference_path = FigureReader.key(reference["dir"] / _RESULT)
    reference_motco = reader.json(reference_path).get("software", {}).get("motco")
    forced_seen = set()
    for a in analyses:
        path = FigureReader.key(a["dir"] / _RESULT)
        m = _measurement(reader.json(path))
        if m["motco"] != reference_motco:
            errors.append(f"{path} records MOTCO {m['motco']}; {reference_path} records MOTCO {reference_motco}; the "
                          "figures draw analyses of one MOTCO version")
        if m["residual_correction"] != ANALYSIS_RESIDUAL_CORRECTION:
            errors.append(f"{path}: params.residual_correction is {m['residual_correction']!r} (absent means none); "
                          f"expected {ANALYSIS_RESIDUAL_CORRECTION!r} (ANALYSIS_RESIDUAL_CORRECTION)")
        if m["min_components"] != ANALYSIS_RANK_FLOOR:
            errors.append(f"{path}: params.min_components is {m['min_components']!r}; expected {ANALYSIS_RANK_FLOOR} "
                          "(ANALYSIS_RANK_FLOOR)")
        if m["forced_components"] != a["forced"]:
            expected = f"forced rank {a['forced']}" if a["forced"] is not None else "None (CV-selected rank)"
            errors.append(f"{path}: params.forced_components is {m['forced_components']!r}; expected {expected}")
        if m["forced_components"] is not None:
            forced_seen.add(m["forced_components"])
        if m["permutations"] != ANALYSIS_PERMUTATIONS:
            errors.append(f"{path}: params.permutations is {m['permutations']!r}; expected {ANALYSIS_PERMUTATIONS}")
        if m["seed"] != ANALYSIS_SEED:
            errors.append(f"{path}: params.seed is {m['seed']!r}; expected {ANALYSIS_SEED}")
        unavailable = a["family"] == "latent" and a["label"] in LATENT_SENSITIVITY_SPLIT_UNAVAILABLE
        want = ANALYSIS_SPAN_SPLIT_TEST and not unavailable
        if (m["span_split_test"], m["span_split_present"]) != (want, want):
            expected = (f"the split test (params.span_split_test = true and a span_split block; ANALYSIS_SPAN_SPLIT_TEST)"
                        if want else "no split test (not available with a latent adjustment)")
            errors.append(f"{path}: params.span_split_test is {m['span_split_test']}, span_split block "
                          f"{'present' if m['span_split_present'] else 'absent'}; expected {expected}")
        if (m["interaction_scan"], m["interaction_scan_present"]) != (ANALYSIS_INTERACTION_SCAN, ANALYSIS_INTERACTION_SCAN):
            errors.append(f"{path}: params.interaction_scan is {m['interaction_scan']}, interaction_scan block "
                          f"{'present' if m['interaction_scan_present'] else 'absent'}; expected the interaction scan "
                          "(params.interaction_scan = true and an interaction_scan block; ANALYSIS_INTERACTION_SCAN)")
    if errors:
        shown = errors[:40]
        more = f"; ... and {len(errors) - len(shown)} more" if len(errors) > len(shown) else ""
        raise ValueError(f"consistency check failed ({len(errors)} problem(s); nothing was written):\n  "
                         + "\n  ".join(shown) + more)
    return {
        "manifests": len(manifests),
        "manifests_from_registries": n_registered,
        "manifests_found": len(manifests) - n_registered,
        "entries": n_entries,
        "distinct_files": len(files),
        "analyses": len(analyses),
        "motco": reference_motco,
        "measurement": {
            "residual_correction": ANALYSIS_RESIDUAL_CORRECTION,
            "min_components": ANALYSIS_RANK_FLOOR,
            "forced_components": sorted(forced_seen),
            "permutations": ANALYSIS_PERMUTATIONS,
            "seed": ANALYSIS_SEED,
            "span_split_test": ANALYSIS_SPAN_SPLIT_TEST,
            "span_split_unavailable": list(LATENT_SENSITIVITY_SPLIT_UNAVAILABLE),
        },
        "findings": findings,
    }


# ---------------------------------------------------------------------------
# Shared drawing helpers
# ---------------------------------------------------------------------------

_P = FIGURE_PALETTE
_INK, _INK2, _MUTED, _GRID, _AXIS, _NEUTRAL = (_P[k] for k in ("ink", "ink2", "muted", "grid", "axis", "neutral"))
_SEX: dict[str, str] = _P["sex"]
_CELL: dict[str, str] = _P["cell_type"]
_STAGE_MARKERS = {0: "o", 1: "s", 2: "^"}
_STAT_TITLE = {"delta": "delta (size)", "angle": "angle (degrees)", "shape": "shape"}
_DIVERGING = LinearSegmentedColormap.from_list("figure_diverging", _P["diverging"])
# axis_correlations.csv sources → heatmap groups, in drawing order
_SOURCE_GROUPS = {
    "metadata.csv": "design",
    "covariates.csv": "analysis covariates",
    "manifest.json (export, tmm)": "export (TMM, library size)",
    "donor_technical.csv": "technical table",
}

FIGURE_RCPARAMS: dict = {
    "font.family": FIGURE_FONT,
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.titleweight": "semibold",
    "axes.titlelocation": "left",
    "axes.labelsize": 8,
    "axes.labelcolor": _INK2,
    "axes.edgecolor": _AXIS,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": _GRID,
    "grid.linewidth": 0.6,
    "axes.axisbelow": True,
    "xtick.color": _MUTED,
    "ytick.color": _MUTED,
    "xtick.labelcolor": _INK2,
    "ytick.labelcolor": _INK2,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "legend.frameon": False,
    "text.color": _INK,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.bbox": "tight",
    "lines.linewidth": 2,
    "lines.solid_capstyle": "round",
    "svg.fonttype": "path",
    "svg.hashsalt": FIGURE_SVG_HASHSALT,
}


def _figure(nrows: int, ncols: int, figsize: tuple[float, float], **kwargs) -> tuple[Figure, np.ndarray]:
    fig = Figure(figsize=figsize)
    axes = fig.subplots(nrows, ncols, squeeze=False, **kwargs)
    return fig, axes


def _suptitle(fig: Figure, text: str) -> None:
    fig.suptitle(text, x=0.01, ha="left", fontsize=10, fontweight="semibold")


def _p(value: float) -> str:
    return f"{value:.3f}"


def _group(values: dict[str, str]) -> str:
    """One value when all labels share it, else ``value (labels); ...`` in first-seen order."""
    groups: dict[str, list[str]] = {}
    for label, value in values.items():
        groups.setdefault(value, [])
        if label not in groups[value]:
            groups[value].append(label)
    if len(groups) == 1:
        return next(iter(groups))
    return "; ".join(f"{value} ({', '.join(labels)})" for value, labels in groups.items())


def _filtering(manifest: dict) -> str:
    """The export's feature filtering, read from an analysis manifest's preprocessing block."""
    pre = manifest["preprocessing"]
    th = pre["thresholds"]
    rule = pre.get("selection_rule") or {}
    selection = _group({b: rule.get(b, "variance") for b in ANALYSIS_BLOCKS})
    return (f"{', '.join(th['excluded_chromosomes'])} removed, CPM ≥ {th['min_cpm']:g} in ≥ k donors, TMM, log2-CPM "
            f"(prior count {th['prior_count']:g}), top {th['N']:,} per block by {selection}")


def _context(items: list[tuple[str, dict, dict]], *, extra: str = "") -> str:
    """The binding reporting context of the analyses drawn, read from their records and manifests.

    ``items`` holds ``(label, result.json record, analysis manifest)`` per analysis.
    """
    covariates = _group({lab: ", ".join(rec["covariates"]["columns"]) for lab, rec, _ in items})
    blocks = _group({lab: ", ".join(f"{b['name']} {b['n_features']:,}" for b in rec["blocks"]) for lab, rec, _ in items})
    filtering = _group({lab: _filtering(man) for lab, _, man in items})
    k = _group({lab: str(man["preprocessing"]["thresholds"]["k"]) for lab, _, man in items})
    rrpp = _group({lab: f"{rec['params']['permutations']} permutations (seed {rec['params']['seed']})"
                   for lab, rec, _ in items})
    correction = _group({lab: f"{rec['params'].get('residual_correction')}" for lab, rec, _ in items})
    kappa = _group({lab: f"{rec['measurement']['residual_correction']['kappa']:.3f}" for lab, rec, _ in items})
    return (f"Covariates (additive): {covariates}. Blocks: {blocks} features. Feature filtering (export): {filtering}; "
            f"k = {k}. RRPP null: {rrpp}, residual correction {correction}, κ = {kappa}. "
            f"Reading: {FIGURE_READING_RULES['non_rejection']}{extra}")


def _footnote(fig: Figure, text: str, *, y: float = -0.01) -> None:
    fig.text(0.01, y, text, fontsize=6.6, color=_INK2, ha="left", va="top", wrap=True)


def _stat_line(record: dict) -> str:
    s, p = record["statistics"], record["p_values"]
    return (f"delta {s['delta']:.2f} (p {_p(p['delta'])}) · angle {s['angle']:.1f}° (p {_p(p['angle'])}, "
            f"CV rank {record['latent']['selected_lv']}) · shape {s['shape']:.3f} (p {_p(p['shape'])})")


def _load_analysis(r: FigureReader, directory: Path, cell_type: str) -> dict:
    """A MOTCO analysis directory with covariate-adjusted donor scores, checked against ``stage_means.csv``.

    Adjusted score = raw PLS score − (x − x_ref) · β, the covariate
    contribution removed at MOTCO's LS-mean reference row; its sex × stage
    cell means must equal ``stage_means.csv`` to ``FIGURE_STAGE_MEANS_ATOL``.
    """
    record = r.result(directory / _RESULT)
    manifest = r.manifest(directory / _ANALYSIS_MANIFEST, cell_type)
    scores = r.csv(directory / "latent_scores.csv", index_col=0, dtype={"donor_id": str})
    meta = r.csv(directory / "analysis_metadata.csv", index_col=0, dtype={"donor_id": str}).loc[scores.index]
    means_path = directory / "stage_means.csv"
    means = r.csv(means_path)
    cov = record["covariates"]
    columns = cov["design_columns"]
    beta = pd.DataFrame(cov["coefficients"]).T.loc[columns, scores.columns].to_numpy()
    offset = (meta[columns] - pd.Series(cov["reference_row"])[columns]).to_numpy() @ beta
    adjusted = scores - offset
    cells = adjusted.groupby([meta["sex"], meta["stage"]]).mean()
    for _, row in means.iterrows():
        cell = (row["group"], int(row["stage"]))
        if cell not in cells.index:
            raise ValueError(f"{means_path}: cell {cell} has no donors in latent_scores.csv")
        for comp in scores.columns:
            ours, theirs = float(cells.loc[cell, comp]), float(row[comp])
            if abs(ours - theirs) > FIGURE_STAGE_MEANS_ATOL:
                raise ValueError(f"{means_path}: adjusted cell mean of {comp} for {cell[0]} stage {cell[1]} is {ours!r}; "
                                 f"the file holds {theirs!r} (tolerance {FIGURE_STAGE_MEANS_ATOL})")
    if len(cells) != len(means):
        raise ValueError(f"{means_path} holds {len(means)} cells; the donors fall in {len(cells)}")
    return {"record": record, "manifest": manifest, "scores": scores, "adjusted": adjusted, "meta": meta,
            "means": means, "dir": directory}


def _geometry(r: FigureReader, cell_type: str, cohort: str, analysis_dir: Path) -> pd.DataFrame:
    """The primary rows of ``geometry_shares.csv``, the diagnostic checked to be made from ``analysis_dir``."""
    d = cell_type_dirs(cell_type)["motco_diagnostics"] / cohort
    manifest = r.manifest(d / DIAGNOSTIC_OUTPUTS["manifest"], cell_type)
    _check_made_from(manifest, d / DIAGNOSTIC_OUTPUTS["manifest"], analysis_dir / _RESULT)
    shares = r.csv(d / DIAGNOSTIC_OUTPUTS["geometry"])
    return shares[shares["analysis"] == DIAGNOSTIC_PRIMARY_LABEL]


def _share(shares: pd.DataFrame, quantity: str, components) -> float:
    sub = shares[(shares["quantity"] == quantity) & (shares["component"].isin(list(components)))]
    return float(sub["share"].sum())


_LABEL_OFFSETS = (
    (4, 2, "left", "bottom"), (-4, 2, "right", "bottom"), (4, -2, "left", "top"), (-4, -2, "right", "top"),
    (4, 9, "left", "bottom"), (-4, 9, "right", "bottom"), (4, -9, "left", "top"), (-4, -9, "right", "top"),
    (0, 12, "center", "bottom"), (0, -12, "center", "top"),
)


def _place_labels(fig: Figure, ax, labels: list[tuple[str, float, float]], avoid: np.ndarray, *,
                  fontsize: float = 6.3, marker_pt: float = 3.0) -> None:
    """Label points without collisions: each label takes the first of ``_LABEL_OFFSETS`` whose box stays in the
    axes and overlaps no label placed before it and no point of ``avoid`` (data coordinates) but its own, or the
    least-overlapping one; deterministic, in the order given. Call once the axes positions and limits are final.
    """
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    renderer = FigureCanvasAgg(fig).get_renderer()
    axbox = ax.get_window_extent(renderer)
    pts = ax.transData.transform(np.asarray(avoid, dtype=float).reshape(-1, 2))
    pad = marker_pt * fig.dpi / 72
    placed = []
    for text, x, y in labels:
        own = ax.transData.transform((x, y))
        others = pts[np.hypot(pts[:, 0] - own[0], pts[:, 1] - own[1]) > 1e-6]
        best = None
        for dx, dy, ha, va in _LABEL_OFFSETS:
            ann = ax.annotate(text, (x, y), xytext=(dx, dy), textcoords="offset points", ha=ha, va=va,
                              fontsize=fontsize, color=_INK)
            box = ann.get_window_extent(renderer)
            hits = int(np.sum((others[:, 0] > box.x0 - pad) & (others[:, 0] < box.x1 + pad)
                              & (others[:, 1] > box.y0 - pad) & (others[:, 1] < box.y1 + pad)))
            area = sum(max(0.0, min(box.x1, b.x1) - max(box.x0, b.x0)) * max(0.0, min(box.y1, b.y1) - max(box.y0, b.y0))
                       for b in placed)
            outside = not (axbox.x0 <= box.x0 and box.x1 <= axbox.x1 and axbox.y0 <= box.y0 and box.y1 <= axbox.y1)
            score = (outside, hits + (area > 0), area)
            ann.remove()
            if best is None or score < best[0]:
                best = (score, (dx, dy, ha, va), box)
            if score == (False, 0, 0):
                break
        dx, dy, ha, va = best[1]
        leader = {"arrowstyle": "-", "color": _MUTED, "lw": 0.5, "shrinkA": 0, "shrinkB": 2} if abs(dy) > 2 else None
        ax.annotate(text, (x, y), xytext=(dx, dy), textcoords="offset points", ha=ha, va=va, fontsize=fontsize,
                    color=_INK, arrowprops=leader)
        placed.append(best[2])


# ---------------------------------------------------------------------------
# The figures. Each returns the figure and the values its caption reads.
# ---------------------------------------------------------------------------

def _draw_cohorts(r: FigureReader) -> tuple[Figure, dict]:
    ncols = 1 + len(EXPORT_COHORTS)
    fig, axes = _figure(len(CELL_TYPES), ncols, (10, 6.2), gridspec_kw={"width_ratios": [1.25] + [1] * len(EXPORT_COHORTS)})
    first = EXPORT_COHORTS[0]
    stages = list(STAGE_LABELS)
    values: dict = {}
    counts_all: dict = {}
    for row, cell_type in enumerate(CELL_TYPES):
        pb = cell_type_dirs(cell_type)["pseudobulk"]
        manifest = r.manifest(pb / _PSEUDOBULK_MANIFEST, cell_type)
        table = r.csv(pb / "donor_table.csv", dtype={"donor_id": str, "stage": str})
        values[cell_type] = {"n_donors": len(table), "cohorts": {}}
        ax = axes[row, 0]
        for sex, color in _SEX.items():
            for inside in (True, False):
                sub = table[(table["sex"] == sex) & (table[first] == inside)]
                ax.scatter(sub["n_rna"], sub["n_atac"], s=22, facecolor=color if inside else "white",
                           edgecolor=color if not inside else "white", linewidth=1.0 if not inside else 0.8,
                           label=sex if inside else None, zorder=3)
        ax.set_xscale("log")
        ax.set_yscale("log")
        for threshold in COHORT_THRESHOLDS.values():
            ax.axvline(threshold, color=_MUTED, linewidth=0.8, zorder=2)
            ax.axhline(threshold, color=_MUTED, linewidth=0.8, zorder=2)
            ax.annotate(f"≥ {threshold}", (threshold, table["n_atac"].min()), xytext=(3, 0),
                        textcoords="offset points", color=_INK2, fontsize=7, va="bottom")
        ax.set_xlabel(f"{cell_type} nuclei, RNA (log)")
        ax.set_ylabel(f"{cell_type} nuclei, ATAC (log)")
        ax.set_title(f"{cell_type}: nuclei per donor (n = {len(table)} donors)")
        if row == 0:
            ax.scatter([], [], s=22, facecolor="white", edgecolor=_MUTED, linewidth=1.0, label=f"outside {first}")
            ax.legend(loc="upper left")
        for col, cohort in enumerate(EXPORT_COHORTS, start=1):
            sub = table[table[cohort]]
            counts = sub.groupby(["stage", "sex"]).size().unstack(fill_value=0).reindex(
                index=stages, columns=list(_SEX), fill_value=0)
            cells = {sex: [int(counts.loc[s, sex]) for s in stages] for sex in _SEX}
            recorded = manifest["cohorts"][cohort]
            if len(sub) != recorded["n"] or cells != {sex: list(recorded["cells"][sex]) for sex in _SEX}:
                raise ValueError(
                    f"{pb / 'donor_table.csv'}: {cell_type} {cohort} has n = {len(sub)} and cells {cells}; "
                    f"{pb / _PSEUDOBULK_MANIFEST} records n = {recorded['n']} and cells {recorded['cells']}"
                )
            k = int(counts.to_numpy().min())
            values[cell_type]["cohorts"][cohort] = {"n": len(sub), "cells": cells, "k": k}
            counts_all[(row, col)] = counts
    top = max(int(c.to_numpy().max()) for c in counts_all.values())
    for (row, col), counts in counts_all.items():
        ax = axes[row, col]
        cohort = EXPORT_COHORTS[col - 1]
        x = np.arange(len(counts))
        width = 0.36
        for k, sex in enumerate(_SEX):
            pos = x + (k - 0.5) * (width + 0.04)
            bars = ax.bar(pos, counts[sex], width=width, color=_SEX[sex], label=sex)
            for b, v in zip(bars, counts[sex]):
                ax.text(b.get_x() + b.get_width() / 2, v + 0.4, str(v), ha="center", va="bottom", fontsize=7,
                        color=_INK2)
        ax.set_xticks(x, [f"{s}\n{STAGE_LABELS[s]}" for s in counts.index])
        ax.set_ylim(0, top * 1.1 + 1)
        ax.grid(axis="x", visible=False)
        ax.set_ylabel("donors")
        v = values[list(CELL_TYPES)[row]]["cohorts"][cohort]
        ax.set_title(f"{cohort}: n = {v['n']}, smallest cell k = {v['k']}")
    _suptitle(fig, "Donor cohorts: nuclei thresholds in both omics, sex × ADNC stage cells")
    fig.tight_layout()
    return fig, values


def _caption_cohorts(v: dict) -> str:
    parts = []
    for cell_type, entry in v.items():
        cohorts = "; ".join(
            f"{cohort} n = {c['n']} (" + ", ".join(f"{sex} {' / '.join(map(str, cells))}" for sex, cells in c["cells"].items())
            + f" by stage {' / '.join(STAGE_LABELS)}; smallest cell k = {c['k']})"
            for cohort, c in entry["cohorts"].items()
        )
        parts.append(f"{cell_type}: {entry['n_donors']} donors in the donor table; {cohorts}.")
    return " ".join(parts)


def _draw_trajectories(r: FigureReader) -> tuple[Figure, dict]:
    fig, axes = _figure(len(CELL_TYPES), len(EXPORT_COHORTS), (10, 9))
    comps = ("pls_0", "pls_1")
    values: dict = {}
    items = []
    for row, cell_type in enumerate(CELL_TYPES):
        for col, cohort in enumerate(EXPORT_COHORTS):
            ax = axes[row, col]
            a = _load_analysis(r, cell_type_dirs(cell_type)["motco_analysis"] / cohort, cell_type)
            adj, meta, means, record = a["adjusted"], a["meta"], a["means"], a["record"]
            items.append((f"{cell_type} {cohort}", record, a["manifest"]))
            for sex, color in _SEX.items():
                for stage, marker in _STAGE_MARKERS.items():
                    m = (meta["sex"] == sex) & (meta["stage"] == stage)
                    ax.scatter(adj.loc[m, comps[0]], adj.loc[m, comps[1]], s=14, marker=marker, color=color,
                               alpha=0.28, linewidth=0, zorder=2)
                traj = means[means["group"] == sex].sort_values("stage")
                for (_, a0), (_, a1) in zip(traj.iloc[:-1].iterrows(), traj.iloc[1:].iterrows()):
                    ax.annotate("", xy=(a1[comps[0]], a1[comps[1]]), xytext=(a0[comps[0]], a0[comps[1]]),
                                arrowprops={"arrowstyle": "-|>", "color": color, "lw": 2, "mutation_scale": 11,
                                            "shrinkA": 7, "shrinkB": 9}, zorder=5)
                for _, t in traj.iterrows():
                    ax.scatter(t[comps[0]], t[comps[1]], s=70, marker=_STAGE_MARKERS[int(t["stage"])], color=color,
                               edgecolor="white", linewidth=2, zorder=6)
                    ax.annotate(str(int(t["stage"])), (t[comps[0]], t[comps[1]]), xytext=(6, 5),
                                textcoords="offset points", fontsize=7.5, color=_INK, fontweight="semibold", zorder=7)
            geometry = _geometry(r, cell_type, cohort, a["dir"])
            shares = {sex: _share(geometry, f"net_displacement_{sex}", comps) for sex in _SEX}
            values[(cell_type, cohort)] = {
                "n": len(meta), "rank": record["latent"]["selected_lv"], "statistics": record["statistics"],
                "p_values": record["p_values"], "plane_shares": shares,
            }
            share_text = ", ".join(f"{sex[0]} {share:.0%}" for sex, share in shares.items())
            ax.set_title(f"{cell_type} · {cohort} (n = {len(meta)}, rank {record['latent']['selected_lv']})\n"
                         f"{_stat_line(record)}", fontsize=8)
            ax.text(0.01, 0.01, f"share of net displacement in this plane: {share_text}", transform=ax.transAxes,
                    fontsize=7, color=_INK2, va="bottom")
            ax.set_xlabel("pls_0 (covariate-adjusted)")
            ax.set_ylabel("pls_1 (covariate-adjusted)")
            ax.axhline(0, color=_AXIS, linewidth=0.6, zorder=1)
            ax.axvline(0, color=_AXIS, linewidth=0.6, zorder=1)
    handles = [Line2D([], [], color=c, label=s) for s, c in _SEX.items()] + [
        Line2D([], [], marker=m, color=_INK2, linestyle="", markersize=6, label=f"stage {s} {STAGE_LABELS[str(s)]}")
        for s, m in _STAGE_MARKERS.items()
    ]
    _suptitle(fig, "Sex × stage trajectories in the first two latent components: LS-means (large) and "
                   "covariate-adjusted donors (faint)")
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    fig.legend(handles=handles, loc="lower center", ncol=5, bbox_to_anchor=(0.5, 0.0))
    _footnote(fig, _context(items, extra=" Donors: raw PLS score minus the covariate contribution at the LS-mean "
                                         "reference row; their cell means equal MOTCO's stage_means.csv."))
    return fig, values


def _caption_trajectories(v: dict) -> str:
    parts = []
    for (cell_type, cohort), e in v.items():
        s, p = e["statistics"], e["p_values"]
        parts.append(
            f"{cell_type} {cohort} (n = {e['n']}, CV-selected rank {e['rank']}): delta {s['delta']:.2f} "
            f"(p {_p(p['delta'])}), angle {s['angle']:.1f}° (p {_p(p['angle'])}, rank {e['rank']}), shape "
            f"{s['shape']:.3f} (p {_p(p['shape'])}); share of each sex's net displacement in the pls_0 × pls_1 plane: "
            + ", ".join(f"{sex} {share:.0%}" for sex, share in e["plane_shares"].items()) + "."
        )
    return " ".join(parts)


def _draw_components(r: FigureReader) -> tuple[Figure, dict]:
    cohort = FIGURE_COHORT
    analyses = {ct: _load_analysis(r, cell_type_dirs(ct)["motco_analysis"] / cohort, ct) for ct in CELL_TYPES}
    ncols = max(a["adjusted"].shape[1] for a in analyses.values())
    fig, axes = _figure(len(CELL_TYPES), ncols, (2.15 * ncols, 5.4))
    rng = np.random.default_rng(FIGURE_JITTER_SEEDS["components_by_stage"])
    values: dict = {}
    empty = []
    for row, (cell_type, a) in enumerate(analyses.items()):
        adj, meta, means = a["adjusted"], a["meta"], a["means"]
        shares = _geometry(r, cell_type, cohort, a["dir"])
        values[cell_type] = {"rank": a["record"]["latent"]["selected_lv"], "n": len(meta), "components": {}}
        for col in range(ncols):
            ax = axes[row, col]
            if col >= adj.shape[1]:
                ax.axis("off")
                empty.append(ax)
                continue
            comp = adj.columns[col]
            for k, (sex, color) in enumerate(_SEX.items()):
                dx = (k - 0.5) * 0.24
                m = meta["sex"] == sex
                xs = meta.loc[m, "stage"] + dx + rng.uniform(-0.06, 0.06, int(m.sum()))
                ax.scatter(xs, adj.loc[m, comp], s=9, color=color, alpha=0.3, linewidth=0)
                traj = means[means["group"] == sex].sort_values("stage")
                ax.plot(traj["stage"] + dx, traj[comp], color=color, marker="o", markersize=5,
                        markeredgecolor="white", markeredgewidth=1.5, label=sex)
            nd = {sex: _share(shares, f"net_displacement_{sex}", [comp]) for sex in _SEX}
            sd = _share(shares, "sex_difference_centred", [comp])
            values[cell_type]["components"][comp] = {"net_displacement": nd, "sex_difference_centred": sd}
            ax.set_title(f"{comp}\nnet displacement F {nd['Female']:.0%} · M {nd['Male']:.0%}\n"
                         f"centred sex difference {sd:.0%}", fontsize=7.5)
            ax.set_xticks([0, 1, 2])
            ax.set_xlim(-0.5, 2.5)
            ax.grid(axis="x", visible=False)
            if col == 0:
                ax.set_ylabel(f"{cell_type}\nadjusted score")
            below = [b for b in range(row + 1, len(CELL_TYPES))
                     if col < list(analyses.values())[b]["adjusted"].shape[1]]
            if not below:
                ax.set_xlabel("ADNC stage")
    handles = [Line2D([], [], color=c, marker="o", markersize=5, markeredgecolor="white", label=s)
               for s, c in _SEX.items()]
    _suptitle(fig, f"{cohort}: LS-mean trajectory per latent component, with the share of each sex's net "
                   "displacement and of the centred sex difference along it")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    if empty:
        empty[0].legend(handles=handles, loc="center", title="LS-means", title_fontsize=7.5)
    else:
        fig.legend(handles=handles, loc="lower center", ncol=len(handles), title="LS-means", bbox_to_anchor=(0.5, -0.04))
    return fig, values


def _caption_components(v: dict) -> str:
    parts = []
    for cell_type, e in v.items():
        comps = "; ".join(
            f"{comp} F {c['net_displacement']['Female']:.0%}, M {c['net_displacement']['Male']:.0%}, sex difference "
            f"{c['sex_difference_centred']:.0%}"
            for comp, c in e["components"].items()
        )
        parts.append(f"{cell_type} {FIGURE_COHORT} (n = {e['n']}, CV-selected rank {e['rank']}): shares of the net "
                     f"displacement and of the centred sex difference per component: {comps}.")
    return " ".join(parts)


def _draw_nulls(r: FigureReader) -> tuple[Figure, dict]:
    rows = [(ct, c) for ct in CELL_TYPES for c in EXPORT_COHORTS]
    records, items = {}, []
    for ct, c in rows:
        d = cell_type_dirs(ct)["motco_analysis"] / c
        records[(ct, c)] = r.result(d / _RESULT)
        items.append((f"{ct} {c}", records[(ct, c)], r.manifest(d / _ANALYSIS_MANIFEST, ct)))
    fig, axes = _figure(1, len(ANALYSIS_STATISTICS), (11, 3.4), sharey=True)
    axes = axes[0]
    for ax, stat in zip(axes, ANALYSIS_STATISTICS):
        for y, key in enumerate(rows):
            rec = records[key]
            ns = rec["null_summary"][stat]
            ax.plot([0, ns["q95"]], [y, y], color=_NEUTRAL, linewidth=9, solid_capstyle="butt", zorder=1)
            ax.plot([ns["q95"], ns["q99"]], [y, y], color=_GRID, linewidth=3, solid_capstyle="butt", zorder=1)
            ax.plot([ns["q50"]] * 2, [y - 0.22, y + 0.22], color=_MUTED, linewidth=1.2, zorder=2)
            ax.plot([ns["q95"]] * 2, [y - 0.22, y + 0.22], color=_INK2, linewidth=1.2, zorder=2)
            obs = rec["statistics"][stat]
            ax.scatter(obs, y, s=46, color=_CELL[key[0]], edgecolor="white", linewidth=2, zorder=3)
            note = f"p {_p(rec['p_values'][stat])}"
            if stat == "angle":
                cs = rec["config_spectrum"]
                pe = cs["permutation_pooled_eigengap"]
                note += (f" · rank {rec['latent']['selected_lv']} · eigengap {cs['pooled']['relative_eigengap']:.2f}"
                         f" (null {pe['q05']:.2f}–{pe['q95']:.2f})")
            ax.annotate(note, (max(obs, ns["q99"]), y), xytext=(6, 0), textcoords="offset points", va="center",
                        fontsize=6.8, color=_INK2)
        ax.set_title(_STAT_TITLE[stat])
        ax.set_xlim(left=0)
        ax.set_xlim(right=ax.get_xlim()[1] * (2.15 if stat == "angle" else 1.45))
        ax.grid(axis="y", visible=False)
    axes[0].set_yticks(range(len(rows)), [f"{ct} · {c}" for ct, c in rows])
    axes[0].invert_yaxis()
    legend = [
        Line2D([], [], color=_NEUTRAL, linewidth=7, label="RRPP null: 0 to 95th pct"),
        Line2D([], [], color=_GRID, linewidth=3, label="95th to 99th pct"),
        Line2D([], [], color=_MUTED, linewidth=1.2, label="null median"),
        Line2D([], [], color=_INK2, linewidth=1.2, label="null 95th pct"),
    ] + [Line2D([], [], marker="o", linestyle="", color=_CELL[ct], markersize=6, label=f"observed, {ct}")
         for ct in CELL_TYPES]
    rejected = sum(rec["p_values"][s] <= ANALYSIS_ALPHA for rec in records.values() for s in ANALYSIS_STATISTICS)
    n_tests = len(ANALYSIS_STATISTICS) * len(records)
    _suptitle(fig, f"Primary analyses: observed statistics against their RRPP null; {rejected} of {n_tests} tests "
                   f"rejected at p ≤ {ANALYSIS_ALPHA}")
    fig.tight_layout(rect=(0, 0.07, 1, 0.94))
    fig.legend(handles=legend, loc="lower center", ncol=len(legend), bbox_to_anchor=(0.5, 0.0))
    _footnote(fig, _context(items), y=0.0)
    values = {"rows": {}, "rejected": rejected, "n_tests": n_tests}
    for key, rec in records.items():
        cs = rec["config_spectrum"]
        values["rows"][key] = {
            "rank": rec["latent"]["selected_lv"],
            "eigengap": cs["pooled"]["relative_eigengap"],
            "eigengap_null": (cs["permutation_pooled_eigengap"]["q05"], cs["permutation_pooled_eigengap"]["q95"]),
            "stats": {s: {"observed": rec["statistics"][s], "q50": rec["null_summary"][s]["q50"],
                          "q95": rec["null_summary"][s]["q95"], "p": rec["p_values"][s]} for s in ANALYSIS_STATISTICS},
        }
    return fig, values


def _caption_nulls(v: dict) -> str:
    parts = [f"{v['rejected']} of {v['n_tests']} tests rejected at p <= {ANALYSIS_ALPHA}."]
    for (cell_type, cohort), e in v["rows"].items():
        stats = "; ".join(
            f"{s} {x['observed']:.{1 if s == 'angle' else 2 if s == 'delta' else 3}f} against null median "
            f"{x['q50']:.{1 if s == 'angle' else 2 if s == 'delta' else 3}f} and 95th percentile "
            f"{x['q95']:.{1 if s == 'angle' else 2 if s == 'delta' else 3}f} (p {_p(x['p'])})"
            for s, x in e["stats"].items()
        )
        parts.append(f"{cell_type} {cohort} (CV-selected rank {e['rank']}; configuration eigengap "
                     f"{e['eigengap']:.2f}, permutation null {e['eigengap_null'][0]:.2f}–{e['eigengap_null'][1]:.2f}): "
                     f"{stats}.")
    return " ".join(parts)


def _draw_pvalue_overview(r: FigureReader) -> tuple[Figure, dict]:
    rows = []
    primaries = []
    added: dict[str, list[str]] = {}
    for a in registered_analyses():
        rec = r.result(a["dir"] / _RESULT)
        man = r.manifest(a["dir"] / _ANALYSIS_MANIFEST, a["cell_type"])
        not_read = a["family"] == "latent" and LATENT_SENSITIVITY_ARMS[a["label"]]["latent_adjustment"] is not None
        rows.append({**{k: a[k] for k in ("cell_type", "family", "label", "cohort")},
                     "rank": rec["latent"]["selected_lv"], "selection": rec["latent"]["component_selection"],
                     "not_read": not_read, **{f"{s}_p": rec["p_values"][s] for s in ANALYSIS_STATISTICS},
                     "columns": tuple(rec["covariates"]["columns"])})
        if a["family"] == "primary":
            primaries.append((f"{a['cell_type']} {a['cohort']}", rec, man))
    table = pd.DataFrame(rows)
    primary_columns = set(table.loc[table["family"] == "primary", "columns"].iloc[0])
    for _, row in table.iterrows():
        extra = [c for c in row["columns"] if c not in primary_columns]
        if extra:
            added.setdefault(", ".join(extra), [])
            if row["label"] not in added[", ".join(extra)]:
                added[", ".join(extra)].append(row["label"])
    keys = table[["cell_type", "family", "label"]].drop_duplicates().reset_index(drop=True)
    n = len(keys)
    fig, axes = _figure(1, len(ANALYSIS_STATISTICS), (11, 0.27 * n + 1.6), sharey=True)
    axes = axes[0]
    for ax, stat in zip(axes, ANALYSIS_STATISTICS):
        ax.axvspan(0, ANALYSIS_ALPHA, color=_P["alpha_band"], alpha=0.45, linewidth=0, zorder=0)
        ax.axvline(ANALYSIS_ALPHA, color=_P["critical"], linewidth=0.8, zorder=1)
        for y, key in keys.iterrows():
            for cohort in EXPORT_COHORTS:
                t = table[(table[["cell_type", "family", "label"]] == key.values).all(axis=1) & (table["cohort"] == cohort)].iloc[0]
                not_read = bool(t["not_read"]) and stat in ("angle", "shape")
                color = _MUTED if not_read else _CELL[t["cell_type"]]
                face = "white" if cohort != EXPORT_COHORTS[0] else color
                ax.scatter(t[f"{stat}_p"], y, s=30, marker="o", facecolor=face, edgecolor=color, linewidth=1.4,
                           zorder=3)
        ax.set_xlim(-0.015, 1.015)
        ax.set_xticks([0, ANALYSIS_ALPHA, 0.25, 0.5, 0.75, 1], ["0", "α", "0.25", "0.5", "0.75", "1"])
        ax.set_title(f"{_STAT_TITLE[stat]}: RRPP p-value")
        ax.grid(axis="y", visible=False)
    labels = []
    for _, key in keys.iterrows():
        ranks = table[(table[["cell_type", "family", "label"]] == key.values).all(axis=1)].set_index("cohort")["rank"]
        labels.append(f"{key['cell_type']} · {key['label']}  [{' / '.join(str(ranks[c]) for c in EXPORT_COHORTS)}]")
    axes[0].set_yticks(range(n), labels)
    axes[0].invert_yaxis()
    previous = None
    for y, key in keys.iterrows():
        group = (key["cell_type"], key["family"])
        if previous is not None and group != previous:
            for ax in axes:
                ax.axhline(y - 0.5, color=_AXIS, linewidth=0.6)
        previous = group
    not_read_arms = [arm for arm, e in LATENT_SENSITIVITY_ARMS.items() if e["latent_adjustment"] is not None]
    legend = [
        Line2D([], [], marker="o", linestyle="", color=_INK2, markersize=6, label=EXPORT_COHORTS[0]),
        *[Line2D([], [], marker="o", linestyle="", markerfacecolor="white", markeredgecolor=_INK2, markersize=6,
                 label=c) for c in EXPORT_COHORTS[1:]],
        *[Line2D([], [], marker="o", linestyle="", color=_CELL[ct], markersize=6, label=ct) for ct in CELL_TYPES],
        Line2D([], [], marker="o", linestyle="", color=_MUTED, markersize=6,
               label=f"{', '.join(not_read_arms)} angle / shape: reported, not read (residualized_basis)"),
    ]
    p_columns = [f"{s}_p" for s in ANALYSIS_STATISTICS]
    rejected = int((table[p_columns] <= ANALYSIS_ALPHA).to_numpy().sum())
    _suptitle(fig, f"All {len(table)} MOTCO analyses (primary and sensitivity): {rejected} of "
                   f"{len(ANALYSIS_STATISTICS) * len(table)} tests rejected at p ≤ {ANALYSIS_ALPHA}; "
                   f"[latent rank {' / '.join(EXPORT_COHORTS)}]")
    fig.tight_layout(rect=(0, 0.035, 1, 0.97))
    fig.legend(handles=legend, loc="lower center", ncol=len(legend), bbox_to_anchor=(0.5, 0.0))
    added_text = "; ".join(f"+ {cols} ({', '.join(labels)})" for cols, labels in added.items())
    _footnote(fig, "Every analysis except the primaries is sensitivity only: a status that differs from the primary's "
                   "is rank, preprocessing, covariate-set, basis or cohort dependence, not a result of the variant; "
                   "each changes only its registered setting (forced rank, preprocessing variant, covariates or latent "
                   f"basis: {added_text}, donor exclusion). Primaries: " + _context(primaries), y=0.0)
    rejections = [
        {"cell_type": t["cell_type"], "label": t["label"], "cohort": t["cohort"], "statistic": s, "p": float(t[f"{s}_p"])}
        for _, t in table.iterrows() for s in ANALYSIS_STATISTICS if t[f"{s}_p"] <= ANALYSIS_ALPHA
    ]
    values = {"n_analyses": len(table), "n_p_values": len(ANALYSIS_STATISTICS) * len(table), "rejected": rejected,
              "rejections": rejections, "n_rows": n,
              "ranks": {f"{t['cell_type']} {t['label']} {t['cohort']}": int(t["rank"]) for _, t in table.iterrows()}}
    return fig, values


def _caption_pvalue_overview(v: dict) -> str:
    listed = "; ".join(f"{x['cell_type']} {x['label']} {x['cohort']} {x['statistic']} p {_p(x['p'])}"
                       for x in v["rejections"]) or "none"
    ranks = ", ".join(f"{name} {rank}" for name, rank in v["ranks"].items())
    return (f"{v['n_analyses']} analyses in {v['n_rows']} rows, {v['n_p_values']} p-values; {v['rejected']} at "
            f"p <= {ANALYSIS_ALPHA}: {listed}. Latent rank per analysis (CV-selected, or forced): {ranks}.")


def _draw_rank_sensitivity(r: FigureReader) -> tuple[Figure, dict]:
    cell_type = FIGURE_SENSITIVITY_CELL_TYPE
    d = cell_type_dirs(cell_type)
    fig, axes = _figure(len(EXPORT_COHORTS), len(ANALYSIS_STATISTICS), (11, 5.6))
    ranks = np.array(ANALYSIS_FORCED_RANKS)
    values: dict = {}
    items = []
    for row, cohort in enumerate(EXPORT_COHORTS):
        base = d["motco_rank_sensitivity"] / cohort
        manifest = r.manifest(base / RANK_SENSITIVITY_MANIFEST, cell_type)
        cv_dir = d["motco_analysis"] / cohort
        _check_made_from(manifest, base / RANK_SENSITIVITY_MANIFEST, cv_dir / _RESULT)
        forced = {}
        for rank in ANALYSIS_FORCED_RANKS:
            rd = base / RANK_SENSITIVITY_DIR_TEMPLATE.format(rank)
            forced[rank] = r.result(rd / _RESULT)
            items.append((cohort, forced[rank], r.manifest(rd / _ANALYSIS_MANIFEST, cell_type)))
        cv = r.result(cv_dir / _RESULT)
        items.append((cohort, cv, r.manifest(cv_dir / _ANALYSIS_MANIFEST, cell_type)))
        cv_rank = cv["latent"]["selected_lv"]
        values[cohort] = {"cv_rank": cv_rank, "cv_p": cv["p_values"],
                          "forced_p": {rank: forced[rank]["p_values"] for rank in ANALYSIS_FORCED_RANKS}}
        for col, stat in enumerate(ANALYSIS_STATISTICS):
            ax = axes[row, col]
            q50 = [forced[k]["null_summary"][stat]["q50"] for k in ranks]
            q95 = [forced[k]["null_summary"][stat]["q95"] for k in ranks]
            q99 = [forced[k]["null_summary"][stat]["q99"] for k in ranks]
            ax.fill_between(ranks, 0, q95, color=_NEUTRAL, linewidth=0, zorder=1, label="null 0–95th pct")
            ax.plot(ranks, q50, color=_MUTED, linewidth=1.2, zorder=2, label="null median")
            ax.plot(ranks, q95, color=_INK2, linewidth=1.2, zorder=2, label="null 95th pct")
            ax.plot(ranks, q99, color=_AXIS, linewidth=1.0, zorder=2, label="null 99th pct")
            obs = [forced[k]["statistics"][stat] for k in ranks]
            ax.plot(ranks, obs, color=_CELL[cell_type], marker="o", markersize=5, markeredgecolor="white",
                    markeredgewidth=1.5, zorder=3, label="observed, forced rank")
            ax.scatter(cv_rank, cv["statistics"][stat], s=80, marker="D", color=_INK, edgecolor="white", linewidth=2,
                       zorder=4, label="observed, CV-selected rank")
            top = max(max(q99), max(obs), cv["statistics"][stat]) * 1.18
            ax.set_ylim(0, top)
            for k in ranks:
                ax.text(k, top * 0.985, f"{forced[k]['p_values'][stat]:.2f}", ha="center", va="top", fontsize=6.3,
                        color=_INK2)
            ax.text(ranks[0] - 0.9, top * 0.985, "p", ha="center", va="top", fontsize=6.3, color=_INK2,
                    fontstyle="italic")
            ax.set_xticks(ranks)
            ax.set_xlim(ranks[0] - 1.2, ranks[-1] + 0.6)
            ax.set_title(f"{cohort} · {_STAT_TITLE[stat]} · CV rank {cv_rank}: p {_p(cv['p_values'][stat])}")
            if row == len(EXPORT_COHORTS) - 1:
                ax.set_xlabel("latent rank")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=len(handles), bbox_to_anchor=(0.5, -0.03))
    _suptitle(fig, f"{cell_type}: dimensionality sensitivity (forced-rank p-values along the top); the CV-selected "
                   "analysis stays primary")
    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    _footnote(fig, _context(items, extra=" A status that differs between ranks is rank dependence, not a result at "
                                         "that rank."), y=-0.05)
    return fig, values


def _caption_rank_sensitivity(v: dict) -> str:
    parts = []
    for cohort, e in v.items():
        cv = ", ".join(f"{s} {_p(e['cv_p'][s])}" for s in ANALYSIS_STATISTICS)
        forced = "; ".join(
            f"{s} " + ", ".join(f"{_p(p[s])}" for p in e["forced_p"].values()) for s in ANALYSIS_STATISTICS
        )
        parts.append(f"{FIGURE_SENSITIVITY_CELL_TYPE} {cohort}: CV-selected rank {e['cv_rank']} (p {cv}); forced ranks "
                     f"{', '.join(map(str, e['forced_p']))}: p {forced}.")
    return " ".join(parts)


def _draw_axis_heatmap(r: FigureReader) -> tuple[Figure, dict]:
    cohort = FIGURE_COHORT
    tables = {}
    for ct in CELL_TYPES:
        dd = cell_type_dirs(ct)
        d = dd["motco_diagnostics"] / cohort
        manifest = r.manifest(d / DIAGNOSTIC_OUTPUTS["manifest"], ct)
        _check_made_from(manifest, d / DIAGNOSTIC_OUTPUTS["manifest"], dd["motco_analysis"] / cohort / _RESULT)
        t = r.csv(d / DIAGNOSTIC_OUTPUTS["axis"], dtype={"variable": str, "source": str, "analysis": str})
        tables[ct] = t[t["analysis"] == DIAGNOSTIC_PRIMARY_LABEL]
    order = list(dict.fromkeys(_SOURCE_GROUPS.values()))
    variables: list[str] = []
    groups: dict[str, str] = {}
    for t in tables.values():
        for var, source in t.drop_duplicates("variable")[["variable", "source"]].itertuples(index=False):
            if source not in _SOURCE_GROUPS:
                raise ValueError(f"axis_correlations.csv: variable {var!r} has an unregistered source {source!r}")
            if var not in groups:
                groups[var] = _SOURCE_GROUPS[source]
                variables.append(var)
    variables = sorted(variables, key=lambda v: order.index(groups[v]))
    group_list = [groups[v] for v in variables]
    widths = [tables[ct]["component"].nunique() for ct in CELL_TYPES]
    fig, axes = _figure(1, len(CELL_TYPES), (9.2, 0.165 * len(variables) + 1.6), gridspec_kw={"width_ratios": widths},
                        sharey=True)
    axes = axes[0]
    values: dict = {}
    im = None
    for ax, ct in zip(axes, CELL_TYPES):
        mat = tables[ct].pivot(index="variable", columns="component", values="spearman").reindex(variables)
        mat = mat[sorted(mat.columns, key=lambda c: int(c.rsplit("_", 1)[1]))]
        im = ax.imshow(mat.to_numpy(), cmap=_DIVERGING, vmin=-1, vmax=1, aspect="auto")
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                v = mat.iat[i, j]
                if np.isfinite(v) and abs(v) >= DIAGNOSTIC_CORRELATION_THRESHOLD:
                    ax.text(j, i, f"{v:.2f}".replace("0.", "."), ha="center", va="center", fontsize=5.8,
                            color="white" if abs(v) >= 0.6 else _INK)
        ax.set_xticks(range(mat.shape[1]), mat.columns, rotation=0)
        ax.xaxis.tick_top()
        ax.set_title(f"{ct} (rank {mat.shape[1]})", pad=16)
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_visible(False)
        for i in range(1, len(group_list)):
            if group_list[i] != group_list[i - 1]:
                ax.axhline(i - 0.5, color="white", linewidth=2.5)
        finite = mat.stack().dropna()
        top = {}
        for comp in mat.columns:
            col = mat[comp].dropna()
            var = col.abs().idxmax()
            top[comp] = (var, float(col[var]))
        values[ct] = {"rank": mat.shape[1], "n_cells": int(finite.size),
                      "n_flagged": int((finite.abs() >= DIAGNOSTIC_CORRELATION_THRESHOLD).sum()), "top": top}
    axes[0].set_yticks(range(len(variables)), variables, fontsize=6.3)
    y = 0
    for g in order:
        n = group_list.count(g)
        if n:
            axes[-1].annotate(g, xy=(1.04, y - 0.4), xycoords=("axes fraction", "data"), fontsize=7, color=_INK2,
                              va="top", ha="left", fontweight="semibold")
        y += n
    cbar = fig.colorbar(im, ax=list(axes), orientation="horizontal", fraction=0.025, pad=0.01, aspect=40)
    cbar.set_label("Spearman ρ (latent score vs donor variable)", color=_INK2)
    cbar.outline.set_visible(False)
    _suptitle(fig, f"{cohort}: Spearman ρ of each latent component of the primary analysis with each diagnostic "
                   f"variable;\nvalues shown where |ρ| ≥ {DIAGNOSTIC_CORRELATION_THRESHOLD} (reporting flag)")
    fig.subplots_adjust(top=1 - 1.45 / fig.get_figheight())
    return fig, values


def _caption_axis_heatmap(v: dict) -> str:
    parts = []
    for ct, e in v.items():
        top = "; ".join(f"{comp} {var} {rho:.2f}" for comp, (var, rho) in e["top"].items())
        parts.append(f"{ct} {FIGURE_COHORT} (rank {e['rank']}): {e['n_flagged']} of {e['n_cells']} correlations at "
                     f"|rho| >= {DIAGNOSTIC_CORRELATION_THRESHOLD}; largest |rho| per component: {top}.")
    return " ".join(parts)


def _draw_first_axis(r: FigureReader) -> tuple[Figure, dict]:
    cohort = FIGURE_COHORT
    fig, axes = _figure(len(CELL_TYPES), 2, (9, 7))
    rng = np.random.default_rng(FIGURE_JITTER_SEEDS["first_axis"])
    values: dict = {}
    stages = sorted(int(s) for s in STAGE_LABELS)
    for row, ct in enumerate(CELL_TYPES):
        d = cell_type_dirs(ct)
        ad = d["motco_analysis"] / cohort
        r.result(ad / _RESULT)
        analysis_manifest = r.manifest(ad / _ANALYSIS_MANIFEST, ct)
        scores = r.csv(ad / "latent_scores.csv", index_col=0, dtype={"donor_id": str})
        meta = r.csv(ad / "analysis_metadata.csv", index_col=0, dtype={"donor_id": str}).loc[scores.index]
        dd = d["motco_diagnostics"] / cohort
        diag_manifest = r.manifest(dd / DIAGNOSTIC_OUTPUTS["manifest"], ct)
        _check_made_from(diag_manifest, dd / DIAGNOSTIC_OUTPUTS["manifest"], ad / _RESULT)
        axis_path = dd / DIAGNOSTIC_OUTPUTS["axis"]
        corr = r.csv(axis_path, dtype={"variable": str, "source": str, "analysis": str})
        corr = corr[(corr["analysis"] == DIAGNOSTIC_PRIMARY_LABEL) & (corr["component"] == "pls_0")].reset_index(drop=True)
        stage_rho = float(corr.loc[corr["variable"] == "stage", "spearman"].iloc[0])
        technical = corr[corr["source"].isin(POSITIVE_CONTROL_TECHNICAL_SOURCES)]
        technical = technical[np.isfinite(technical["spearman"].to_numpy(dtype=float))]
        top = technical.loc[technical["spearman"].abs().idxmax()]
        variable, source, top_rho = str(top["variable"]), str(top["source"]), float(top["spearman"])
        if source == "donor_technical.csv":
            tech = r.csv(d["pseudobulk"] / TECHNICAL_OUTPUTS["table"], dtype={"donor_id": str}).set_index("donor_id")
            x = tech.loc[scores.index, variable].astype(float)
        else:
            export = r.manifest(FigureReader.key(analysis_manifest["inputs"]["dir"]) + "/" + _EXPORT_MANIFEST, ct)
            kind, omic = variable.split("_", 1) if variable.startswith("tmm_") else ("log10_lib", variable.rsplit("_", 1)[1])
            tmm = export["tmm"][omic]
            x = pd.Series([float(tmm[dn]["factor"]) if kind == "tmm" else float(np.log10(tmm[dn]["lib_size"]))
                           for dn in scores.index], index=scores.index)
        for name, ours, theirs in (("stage", spearman(scores["pls_0"], meta["stage"])[0], stage_rho),
                                   (variable, spearman(scores["pls_0"], x)[0], top_rho)):
            if not abs(ours - theirs) <= FIGURE_SPEARMAN_ATOL:
                raise ValueError(f"{axis_path}: Spearman of pls_0 with {name} recomputes to {ours!r}; the file holds "
                                 f"{theirs!r} (tolerance {FIGURE_SPEARMAN_ATOL})")
        values[ct] = {"stage_rho": stage_rho, "variable": variable, "source": source, "rho": top_rho,
                      "n": len(scores), "n_technical": len(technical)}
        ax = axes[row, 0]
        for k, (sex, color) in enumerate(_SEX.items()):
            m = meta["sex"] == sex
            ax.scatter(meta.loc[m, "stage"] + (k - 0.5) * 0.22 + rng.uniform(-0.07, 0.07, int(m.sum())),
                       scores.loc[m, "pls_0"], s=18, color=color, edgecolor="white", linewidth=0.6, label=sex)
        med = scores["pls_0"].groupby(meta["stage"]).median()
        ax.hlines(med.values, med.index - 0.3, med.index + 0.3, color=_INK, linewidth=1.5, label="median, all donors",
                  zorder=4)
        ax.set_xticks(stages, [f"{s}\n{STAGE_LABELS[str(s)]}" for s in stages])
        ax.grid(axis="x", visible=False)
        ax.set_ylabel(f"{ct}: pls_0 score (raw)")
        ax.set_title(f"{ct}: pls_0 vs ADNC stage, Spearman ρ = {stage_rho:.2f}")
        if row == 0:
            ax.legend(loc="upper left")
        ax = axes[row, 1]
        for sex, color in _SEX.items():
            m = meta["sex"] == sex
            ax.scatter(x[m], scores.loc[m, "pls_0"], s=18, color=color, edgecolor="white", linewidth=0.6)
        fraction = variable.split("_")[0] in TECHNICAL_GENE_SETS and variable.endswith(("_mean", "_pseudobulk"))
        log = source == "donor_technical.csv" and bool((x > 0).all()) and float(x.max() / x.min()) > 10
        if log:
            ax.set_xscale("log")
        if fraction:
            ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{100 * v:g}%"))
            if log:
                ax.xaxis.set_minor_formatter(mticker.FuncFormatter(
                    lambda v, _: f"{100 * v:g}%" if f"{v:.0e}"[0] in "25" else ""))
        ax.set_xlabel(f"{variable}{' (log scale)' if log else ''}")
        ax.set_title(f"pls_0 vs {variable}, Spearman ρ = {top_rho:.2f}")
    _suptitle(fig, f"{cohort}: the first latent component (raw scores) against ADNC stage and against the\n"
                   "technical variable with the largest |ρ| (E3's rule: technical table, TMM factors, library sizes)")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return fig, values


def _caption_first_axis(v: dict) -> str:
    return " ".join(
        f"{ct} {FIGURE_COHORT} (n = {e['n']}): Spearman of pls_0 with stage {e['stage_rho']:.2f}; the largest |rho| "
        f"among {e['n_technical']} technical variables is {e['variable']} ({e['source']}) at {e['rho']:.2f}."
        for ct, e in v.items()
    )


def _draw_positive_control(r: FigureReader) -> tuple[Figure, dict]:
    control, reference = POSITIVE_CONTROL_CELL_TYPES["control"], POSITIVE_CONTROL_CELL_TYPES["reference"]
    blocks = [(basis, cohort) for basis in POSITIVE_CONTROL_BASES for cohort in EXPORT_COHORTS]
    panels = [("E1", "deseq2", "E1 · DESeq2: exported RNA stage\nfeatures at FDR 0.05 (union)"),
              ("E1", "limma", "E1 · limma: exported RNA stage\nfeatures at FDR 0.05 (joint F)"),
              ("E2", "", "E2 · CV mean AUROC\nof the analysis"),
              ("E3", "", "E3 · control pls_0: |ρ| with stage vs\nthe largest |ρ| with a technical variable")]
    fig, axes = _figure(1, len(panels), (13, 3.9), sharey=True, gridspec_kw={"width_ratios": [1, 1, 1, 1.5]})
    axes = axes[0]
    ylabels = []
    values: dict = {}
    e1_values = []
    e3_names: list = []
    for y, (basis, cohort) in enumerate(blocks):
        d = _positive_control_dir(basis, cohort)
        r.manifest(d / POSITIVE_CONTROL_MANIFEST, control)
        t = r.csv(d / POSITIVE_CONTROL_TABLE, dtype={"method": str})
        t["method"] = t["method"].fillna("")
        outcome = t.loc[t["expectation"] == "overall", "outcome"].iloc[0]
        ylabels.append(f"{basis} · {cohort}\n{'✓' if outcome == 'holds' else '✗'} {outcome.replace('_', ' ')}")
        values[(basis, cohort)] = {"outcome": outcome, "rows": {}}
        for ax, (exp, method, _) in zip(axes, panels):
            row = t[(t["expectation"] == exp) & (t["method"] == method)].iloc[0]
            a, b = float(row["control_value"]), float(row["reference_value"])
            holds = bool(row["holds"])
            values[(basis, cohort)]["rows"][(exp, method)] = {"control": a, "reference": b, "holds": holds,
                                                             "reference_name": str(row["reference"])}
            if exp == "E1":
                e1_values += [a, b]
            cb = _MUTED if exp == "E3" else _CELL[reference]
            fmt = "{:.0f}" if exp == "E1" else ("{:.3f}" if exp == "E2" else "{:.2f}")
            ax.plot([a, b], [y, y], color=_GRID, linewidth=2, zorder=1)
            for value, color, other in ((a, _CELL[control], b), (b, cb, a)):
                ax.scatter([value], [y], s=48, color=color, edgecolor="white", linewidth=2, zorder=3)
                left = value < other or (value == other and color == cb)
                ax.annotate(fmt.format(value), (value, y), xytext=(-1 if left else 1, 6), textcoords="offset points",
                            ha="right" if left else "left", fontsize=6.5, color=_INK2)
            if exp == "E3":
                e3_names.append((ax, str(row["reference"]), b, y))
            ax.annotate("✓ holds" if holds else "✗ fails", (1.0, y), xycoords=("axes fraction", "data"),
                        xytext=(-2, 0), textcoords="offset points", ha="right", va="center", fontsize=7,
                        color=_P["good"] if holds else _P["critical"], fontweight="semibold")
    e1_top = max(10.0, max(e1_values))
    for ax, (exp, _, title) in zip(axes, panels):
        ax.set_title(title, fontsize=8)
        ax.grid(axis="y", visible=False)
        if exp == "E1":
            ax.set_xscale("symlog", linthresh=1)
            ax.set_xlim(-0.4, e1_top * 30)
            ticks = [0, 1] + [10 ** k for k in range(1, int(np.floor(np.log10(e1_top * 30))) + 1)]
            ax.set_xticks(ticks, [f"{t:g}" for t in ticks])
        else:
            lo, hi = ax.get_xlim()
            pad = 0.45 if exp == "E3" else 0.4
            ax.set_xlim(lo - (hi - lo) * 0.15, hi + (hi - lo) * pad)
    for ax, name, b, y in e3_names:
        ax.annotate(name, (b, y), xytext=(0, -8), textcoords="offset points", va="top", fontsize=6, color=_INK2,
                    ha="left" if b <= 0.5 * sum(ax.get_xlim()) else "right")
    axes[0].set_yticks(range(len(blocks)), ylabels, fontsize=7)
    axes[0].set_ylim(len(blocks) - 0.5, -0.6)
    legend = [Line2D([], [], marker="o", linestyle="", color=_CELL[control], markersize=6, label=f"{control} (control)"),
              Line2D([], [], marker="o", linestyle="", color=_CELL[reference], markersize=6,
                     label=f"{reference} (reference)"),
              Line2D([], [], marker="o", linestyle="", color=_MUTED, markersize=6,
                     label="E3: largest technical |ρ| (named below the dot)")]
    _suptitle(fig, f"Positive control ({control} vs {reference}): the control holds when E1–E3 all hold in "
                   f"{POSITIVE_CONTROL_PRIMARY_COHORT}; the {POSITIVE_CONTROL_DEFAULT_BASIS}-basis "
                   f"{POSITIVE_CONTROL_PRIMARY_COHORT} row is its outcome")
    fig.tight_layout(rect=(0, 0.07, 1, 0.94))
    fig.legend(handles=legend, loc="lower center", ncol=len(legend), bbox_to_anchor=(0.5, 0.0))
    other = [b for b in POSITIVE_CONTROL_BASES if b != POSITIVE_CONTROL_DEFAULT_BASIS]
    _footnote(fig, f"E1 does not depend on the basis. {', '.join(other)}: E2 and E3 read from the arm and its "
                   "diagnostic, a diagnosis of the primary-basis control, not a replacement; "
                   f"{', '.join(c for c in EXPORT_COHORTS if c != POSITIVE_CONTROL_PRIMARY_COHORT)} is descriptive. "
                   "No expectation reads a p-value.", y=0.0)
    return fig, values


def _caption_positive_control(v: dict) -> str:
    control, reference = POSITIVE_CONTROL_CELL_TYPES["control"], POSITIVE_CONTROL_CELL_TYPES["reference"]
    parts = []
    for (basis, cohort), e in v.items():
        rows = e["rows"]
        e1 = "; ".join(
            f"{method} {rows[('E1', method)]['control']:.0f} vs {rows[('E1', method)]['reference']:.0f} "
            f"({'holds' if rows[('E1', method)]['holds'] else 'fails'})" for method in ("deseq2", "limma")
        )
        e2, e3 = rows[("E2", "")], rows[("E3", "")]
        parts.append(
            f"{basis} basis, {cohort} ({control} vs {reference}): E1 {e1}; E2 CV mean AUROC {e2['control']:.3f} vs "
            f"{e2['reference']:.3f} ({'holds' if e2['holds'] else 'fails'}); E3 |rho| with stage {e3['control']:.2f} vs "
            f"{e3['reference_name']} {e3['reference']:.2f} ({'holds' if e3['holds'] else 'fails'}); overall "
            f"{e['outcome'].replace('_', ' ')}."
        )
    return " ".join(parts)


def _draw_stage_de(r: FigureReader) -> tuple[Figure, dict]:
    cohort = FIGURE_COHORT
    omic, model, coefficient = FIGURE_STAGE_DE["omic"], FIGURE_STAGE_DE["model"], FIGURE_STAGE_DE["coefficient"]
    joint = STAGE_DE_LIMMA_JOINT_TESTS[model][0]
    tested = [label for label, _ in STAGE_DE_MODELS[model]["tested"]]
    fig, axes = _figure(len(CELL_TYPES), 2, (10, 8.4))
    values: dict = {}
    labels: list = []
    for row, ct in enumerate(CELL_TYPES):
        d = cell_type_dirs(ct)
        de_dir, lm_dir = d["stage_de"] / cohort, d["stage_de_limma"] / cohort
        result = d["motco_analysis"] / cohort / _RESULT
        de_manifest = r.manifest(de_dir / STAGE_DE_MANIFEST, ct)
        _check_made_from(de_manifest, de_dir / STAGE_DE_MANIFEST, result)
        lm_manifest = r.manifest(lm_dir / STAGE_DE_LIMMA_MANIFEST, ct)
        de_table = STAGE_DE_TABLE_TEMPLATE.format(model=model, coefficient=coefficient)
        lm_table = STAGE_DE_LIMMA_TABLE_TEMPLATE.format(model=model, label=coefficient)
        de = r.csv(de_dir / omic / de_table, dtype={"feature": str})
        lm = r.csv(lm_dir / omic / lm_table, dtype={"feature": str})
        de_summary = r.csv(de_dir / STAGE_DE_SUMMARY)
        # The tables drawn are the ones each step recorded, and the limma step was made from the DESeq2
        # outputs drawn beside it (its inputs.deseq2 record, as positive-control checks it)
        dq = lm_manifest["inputs"]["deseq2"]
        for path, recorded, where in (
            (de_dir / omic / de_table, de_manifest["tables"][omic][model][coefficient]["sha256"],
             f"{de_dir / STAGE_DE_MANIFEST} tables"),
            (lm_dir / omic / lm_table, lm_manifest["tables"][omic][model][coefficient]["sha256"],
             f"{lm_dir / STAGE_DE_LIMMA_MANIFEST} tables"),
            (de_dir / omic / de_table, dq["files"][f"{omic}/{de_table}"], f"{lm_dir / STAGE_DE_LIMMA_MANIFEST} inputs.deseq2"),
            (de_dir / STAGE_DE_SUMMARY, dq["summary_sha256"], f"{lm_dir / STAGE_DE_LIMMA_MANIFEST} inputs.deseq2"),
            (de_dir / STAGE_DE_MANIFEST, dq["manifest_sha256"], f"{lm_dir / STAGE_DE_LIMMA_MANIFEST} inputs.deseq2"),
        ):
            if r.sha256(path) != recorded:
                raise ValueError(f"{where} records {path} with sha256 {recorded}; the file read hashes {r.sha256(path)}")
        lm_summary = r.csv(lm_dir / STAGE_DE_LIMMA_SUMMARY)
        comparison_path = lm_dir / STAGE_DE_LIMMA_COMPARISON
        comparison = r.csv(comparison_path)

        def summary(table: pd.DataFrame, coef: str, feature_set: str = "exported") -> pd.Series:
            sub = table[(table["omic"] == omic) & (table["model"] == model) & (table["coefficient"] == coef)
                        & (table["feature_set"] == feature_set)]
            if len(sub) != 1:
                raise ValueError(f"{ct}: {len(sub)} summary rows for {omic} {model} {coef} {feature_set}")
            return sub.iloc[0]

        crow = comparison[(comparison["omic"] == omic) & (comparison["model"] == model)
                          & (comparison["test"] == coefficient) & (comparison["feature_set"] == "exported")].iloc[0]
        sig_all = de["padj"].le(STAGE_DE_ALPHA)
        exported = de["exported"].astype(bool)
        checks = {
            f"{STAGE_DE_SUMMARY} n_significant ({coefficient}, exported)": (int((sig_all & exported).sum()),
                                                                          int(summary(de_summary, coefficient)["n_significant"])),
            f"{STAGE_DE_SUMMARY} n_significant ({coefficient}, abundance)": (int(sig_all.sum()),
                                                                           int(summary(de_summary, coefficient, "abundance")["n_significant"])),
            f"{STAGE_DE_SUMMARY} n_padj_missing ({coefficient}, exported)": (int((de["padj"].isna() & exported).sum()),
                                                                           int(summary(de_summary, coefficient)["n_padj_missing"])),
        }
        m = de[exported].merge(lm[lm["exported"].astype(bool)], on="feature", suffixes=("", "_limma"),
                               validate="one_to_one")
        s_de, s_lm = m["padj"].le(STAGE_DE_ALPHA), m["adj.P.Val"].le(STAGE_DE_ALPHA)
        recomputed = {"n_tested": len(m), "n_significant_deseq2": int(s_de.sum()), "n_significant_limma": int(s_lm.sum()),
                      "n_both": int((s_de & s_lm).sum()), "n_only_deseq2": int((s_de & ~s_lm).sum()),
                      "n_only_limma": int((~s_de & s_lm).sum())}
        for name, ours in recomputed.items():
            checks[f"{STAGE_DE_LIMMA_COMPARISON} {name}"] = (ours, int(crow[name]))
        for name, (ours, theirs) in checks.items():
            if ours != theirs:
                raise ValueError(f"{comparison_path if STAGE_DE_LIMMA_COMPARISON in name else de_dir / STAGE_DE_SUMMARY}: "
                                 f"{ct} {name} recomputes to {ours}; the file holds {theirs}")
        values[ct] = {
            "counts": {coef: int(summary(de_summary, coef)["n_significant"]) for coef in tested + [STAGE_DE_UNION]},
            "limma": {coef: int(summary(lm_summary, coef)["n_significant"]) for coef in tested + [joint]},
            "overlap": {k: recomputed[k] for k in ("n_both", "n_only_deseq2", "n_only_limma")},
            "n_tested": recomputed["n_tested"], "n_abundance": int(summary(de_summary, coefficient, "abundance")["n_tested"]),
            "n_significant_abundance": int(sig_all.sum()),
            "n_padj_missing": int(summary(de_summary, coefficient)["n_padj_missing"]),
            "spearman_log2fc": float(crow["spearman_log2fc"]),
            "covariates": list(de_manifest["covariates"]["per_omic"][omic]),
        }
        ax = axes[row, 0]
        dv = de.dropna(subset=["pvalue"])  # pvalue stays missing where DESeq2 tests nothing; padj also where the independent filter drops
        sig = dv["padj"].le(STAGE_DE_ALPHA)
        ex = dv["exported"].astype(bool)
        yv = -np.log10(dv["pvalue"])
        ax.scatter(dv.loc[~ex, "log2FoldChange"], yv[~ex], s=3, color=_GRID, linewidth=0, label="abundance set only")
        ax.scatter(dv.loc[ex & ~sig, "log2FoldChange"], yv[ex & ~sig], s=4, color=_AXIS, linewidth=0, label="exported")
        ax.scatter(dv.loc[~ex & sig, "log2FoldChange"], yv[~ex & sig], s=16, facecolor="white", edgecolor=_CELL[ct],
                   linewidth=1, label=f"padj ≤ {STAGE_DE_ALPHA}, not exported")
        ax.scatter(dv.loc[ex & sig, "log2FoldChange"], yv[ex & sig], s=16, color=_CELL[ct], edgecolor="white",
                   linewidth=0.6, label=f"padj ≤ {STAGE_DE_ALPHA}, exported")
        top = dv[ex & sig].nsmallest(8, "padj")
        labels.append((ax, [(t["feature"], t["log2FoldChange"], -np.log10(t["pvalue"])) for _, t in top.iterrows()],
                       np.column_stack([dv.loc[sig, "log2FoldChange"], yv[sig]])))
        lim = max(np.nanpercentile(np.abs(dv["log2FoldChange"]), 99.9),
                  dv.loc[sig, "log2FoldChange"].abs().max() if sig.any() else 0.0)
        ax.set_xlim(-lim * 1.1, lim * 1.1)
        ax.set_xlabel(f"log2 fold change, {coefficient} (DESeq2)")
        ax.set_ylabel("−log10 p")
        ax.set_title(f"{ct}: {values[ct]['counts'][coefficient]} exported / {values[ct]['n_significant_abundance']} "
                     f"abundance features at padj ≤ {STAGE_DE_ALPHA}")
        ax.legend(loc="upper left", markerscale=1.5, fontsize=6.5)
        ax = axes[row, 1]
        ax.scatter(m.loc[~s_de & ~s_lm, "log2FoldChange"], m.loc[~s_de & ~s_lm, "logFC"], s=4, color=_AXIS,
                   linewidth=0, label="neither")
        for mask, label in [(s_de & ~s_lm, "DESeq2 only"), (~s_de & s_lm, "limma only"), (s_de & s_lm, "both")]:
            if label == "limma only":
                ax.scatter(m.loc[mask, "log2FoldChange"], m.loc[mask, "logFC"], s=22, marker="s", facecolor="white",
                           edgecolor=_INK, linewidth=1, label=f"{label} ({int(mask.sum())})")
            else:
                ax.scatter(m.loc[mask, "log2FoldChange"], m.loc[mask, "logFC"], s=18,
                           facecolor="white" if label == "DESeq2 only" else _CELL[ct], edgecolor=_CELL[ct],
                           linewidth=1, label=f"{label} ({int(mask.sum())})")
        lo = min(m["log2FoldChange"].quantile(0.001), m["logFC"].quantile(0.001))
        hi = max(m["log2FoldChange"].quantile(0.999), m["logFC"].quantile(0.999))
        ax.plot([lo, hi], [lo, hi], color=_MUTED, linewidth=0.8, zorder=0)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_xlabel("log2 fold change (DESeq2)")
        ax.set_ylabel("logFC (limma-voom)")
        ax.set_title(f"{ct}: DESeq2 vs limma, exported {omic.upper()}; Spearman of log2FC {crow['spearman_log2fc']:.2f}")
        ax.legend(loc="upper left", fontsize=6.5)
    _suptitle(fig, f"{cohort}: pseudobulk stage DE, {model} {coefficient} ({STAGE_DE_MODELS[model]['formula']} "
                   f"+ covariates), {omic.upper()}")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    for ax, texts, avoid in labels:
        _place_labels(fig, ax, texts, avoid)
    return fig, values


def _caption_stage_de(v: dict) -> str:
    model, coefficient = FIGURE_STAGE_DE["model"], FIGURE_STAGE_DE["coefficient"]
    joint = STAGE_DE_LIMMA_JOINT_TESTS[model][0]
    tested = [label for label, _ in STAGE_DE_MODELS[model]["tested"]]
    parts = []
    for ct, e in v.items():
        counts = ", ".join(f"{coef} {e['counts'][coef]}" for coef in tested + [STAGE_DE_UNION])
        limma = ", ".join(f"{coef} {e['limma'][coef]}" for coef in tested + [joint])
        o = e["overlap"]
        parts.append(
            f"{ct} {FIGURE_COHORT}, {FIGURE_STAGE_DE['omic'].upper()} ({e['n_tested']:,} exported of "
            f"{e['n_abundance']:,} abundance features; covariates {', '.join(e['covariates'])}): exported features at "
            f"FDR {STAGE_DE_ALPHA}, DESeq2 {model} {counts}; limma {limma}. {coefficient}: both methods {o['n_both']}, "
            f"DESeq2 only {o['n_only_deseq2']}, limma only {o['n_only_limma']} (DESeq2 padj missing for "
            f"{e['n_padj_missing']:,} exported features, read as not significant); Spearman of log2 fold changes "
            f"{e['spearman_log2fc']:.2f}."
        )
    return " ".join(parts)


_SPLIT_P_TICKS = (0.001, 0.01, ANALYSIS_ALPHA, 0.2, 1.0)


def _split_p_axis(ax) -> None:
    """A log pseudo-F p axis from the smallest attainable p, with the α band."""
    low = 1.0 / (1 + ANALYSIS_PERMUTATIONS)
    ax.axvspan(low * 0.7, ANALYSIS_ALPHA, color=_P["alpha_band"], alpha=0.45, linewidth=0, zorder=0)
    ax.axvline(ANALYSIS_ALPHA, color=_P["critical"], linewidth=0.8, zorder=1)
    ax.set_xscale("log")
    ax.set_xlim(low * 0.7, 1.4)
    ax.set_xticks(_SPLIT_P_TICKS, [f"{x:g}" if x != ANALYSIS_ALPHA else "α" for x in _SPLIT_P_TICKS])
    ax.xaxis.set_minor_locator(mticker.NullLocator())
    ax.grid(axis="y", visible=False)


def _draw_span_split(r: FigureReader) -> tuple[Figure, dict]:
    rows, primaries = [], []
    for a in registered_analyses():
        rec = r.result(a["dir"] / _RESULT)
        man = r.manifest(a["dir"] / _ANALYSIS_MANIFEST, a["cell_type"])
        status = _span_split_status(rec, alpha=ANALYSIS_ALPHA)
        row = {**{k: a[k] for k in ("cell_type", "family", "label", "cohort")},
               "rank": rec["latent"]["selected_lv"], "available": status is not None}
        if status is not None:
            split = rec["span_split"]
            row.update({f"{part}_p": split["parts"][part][SPAN_SPLIT_STATISTIC] for part in SPAN_SPLIT_PARTS})
            row.update({f"{part}_f": split["parts"][part]["f"] for part in SPAN_SPLIT_PARTS})
            row.update({f"{part}_dim": split["parts"][part]["dimension"] for part in SPAN_SPLIT_PARTS})
            row.update(status, share=split["span_share"])
        rows.append(row)
        if a["family"] == "primary":
            primaries.append((f"{a['cell_type']} {a['cohort']}", rec, man))
    table = pd.DataFrame(rows)
    available = table[table["available"]]
    keys = table[["cell_type", "family", "label"]].drop_duplicates().reset_index(drop=True)
    n = len(keys)
    fig = Figure(figsize=(12, 0.27 * n + 2.2))
    gs = fig.add_gridspec(2, 2, width_ratios=(1, 1.2), height_ratios=(1, 1.45))
    ax_a, ax_c, ax_b = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[:, 1])

    # (a) the primaries: every part's pseudo-F p-value, marked by its role
    prim = table[table["family"] == "primary"].reset_index(drop=True)
    _split_p_axis(ax_a)
    offsets = {"total": -0.24, "span": 0.0, "complement": 0.24}
    for y, t in prim.iterrows():
        color = _CELL[t["cell_type"]]
        ax_a.scatter(t["total_p"], y + offsets["total"], s=26, marker="s", color=_MUTED, edgecolor="white",
                     linewidth=1, zorder=3)
        ax_a.scatter(t["span_p"], y + offsets["span"], s=34, marker="D", facecolor="white", edgecolor=color,
                     linewidth=1.5, zorder=3)
        ax_a.scatter(t["complement_p"], y + offsets["complement"], s=40, marker="o", color=color, edgecolor="white",
                     linewidth=1.2, zorder=4)
        ax_a.annotate(f"span {'concordant' if t['span_concordant'] else 'discordant'}", (1.0, y), xytext=(8, 0),
                      textcoords="offset points", va="center", fontsize=6.6, color=_INK2, annotation_clip=False)
    ax_a.set_ylim(len(prim) - 0.5, -0.5)
    ax_a.set_yticks(range(len(prim)), [f"{t['cell_type']} · {t['cohort']}  [rank {t['rank']}]" for _, t in prim.iterrows()])
    ax_a.set_title("(a) primaries: pseudo-F p per part")

    # (b) the complement of every registered analysis that carries the test
    _split_p_axis(ax_b)
    not_available = []
    for y, key in keys.iterrows():
        sub = table[(table[["cell_type", "family", "label"]] == key.values).all(axis=1)].set_index("cohort")
        color = _CELL[key["cell_type"]]
        if not sub["available"].any():
            not_available += [f"{key['cell_type']} {key['label']} {c}" for c in sub.index]
            ax_b.text(0.0012, y, "not available: no split test with a latent adjustment", va="center", ha="left",
                      fontsize=6.6, color=_MUTED, fontstyle="italic")
            continue
        for i, cohort in enumerate(EXPORT_COHORTS):
            t = sub.loc[cohort]
            face = color if cohort == EXPORT_COHORTS[0] else "white"
            ax_b.scatter(t["complement_p"], y + (i - (len(EXPORT_COHORTS) - 1) / 2) * 0.3, s=26, marker="o",
                         facecolor=face, edgecolor=color, linewidth=1.3, zorder=3)
    labels = []
    for _, key in keys.iterrows():
        ranks = table[(table[["cell_type", "family", "label"]] == key.values).all(axis=1)].set_index("cohort")["rank"]
        labels.append(f"{key['cell_type']} · {key['label']}  [{' / '.join(str(ranks[c]) for c in EXPORT_COHORTS)}]")
    ax_b.set_yticks(range(n), labels)
    ax_b.set_ylim(n - 0.5, -0.5)
    previous = None
    for y, key in keys.iterrows():
        group = (key["cell_type"], key["family"])
        if previous is not None and group != previous:
            ax_b.axhline(y - 0.5, color=_AXIS, linewidth=0.6)
        previous = group
    ax_b.set_title(f"(b) complement pseudo-F p per analysis [rank {' / '.join(EXPORT_COHORTS)}]")

    # (c) the in-span share against the fitted rank, no reference line
    shift = {c: (i - (len(EXPORT_COHORTS) - 1) / 2) * 0.22 for i, c in enumerate(EXPORT_COHORTS)}
    for _, t in available.iterrows():
        if pd.isna(t["share"]):
            continue
        color = _CELL[t["cell_type"]]
        face = color if t["cohort"] == EXPORT_COHORTS[0] else "white"
        ax_c.scatter(t["rank"] + shift[t["cohort"]], t["share"], s=26, marker="o", facecolor=face, edgecolor=color,
                     linewidth=1.3, zorder=3)
    ranks = sorted(available["rank"].unique())
    ax_c.set_xticks(ranks)
    ax_c.set_xlim(min(ranks) - 0.8, max(ranks) + 0.8)
    ax_c.set_ylim(bottom=0)
    ax_c.set_xlabel("fitted latent rank (span dimension)")
    ax_c.set_ylabel("in-span share of the interaction SS")
    ax_c.set_title("(c) in-span share by fitted rank (not read)")

    legend = [
        Line2D([], [], marker="o", linestyle="", color=_INK2, markersize=6, label="complement (confirmatory; a)"),
        Line2D([], [], marker="D", linestyle="", markerfacecolor="white", markeredgecolor=_INK2, markersize=5,
               label="span (concordance cell; a)"),
        Line2D([], [], marker="s", linestyle="", color=_MUTED, markersize=6, label="total (not read; a)"),
        Line2D([], [], marker="o", linestyle="", color=_INK2, markersize=6, label=f"{EXPORT_COHORTS[0]} (b, c)"),
        *[Line2D([], [], marker="o", linestyle="", markerfacecolor="white", markeredgecolor=_INK2, markersize=6,
                 label=f"{c} (b, c)") for c in EXPORT_COHORTS[1:]],
        *[Line2D([], [], marker="o", linestyle="", color=_CELL[ct], markersize=6, label=ct) for ct in CELL_TYPES],
    ]
    rejected = available[available["complement_rejected"].astype(bool)]
    _suptitle(fig, f"Split interaction test: complement p ≤ {ANALYSIS_ALPHA} in {len(rejected)} of {len(available)} "
                   f"analyses carrying it; only each cell type's {FIGURE_COHORT} primary complement is confirmatory")
    fig.tight_layout(rect=(0, 0.05, 1, 0.965))
    fig.legend(handles=legend, loc="lower center", ncol=len(legend), bbox_to_anchor=(0.5, 0.012))
    _footnote(fig, "Only the complement is confirmatory (p ≤ α, unadjusted; each primary carries five tests read at α, "
                   "with the interaction scan); "
                   "the span part is a concordance cell against any of delta, angle and shape rejected; total and the "
                   "in-span share are reported and not read; every analysis other than the primaries is sensitivity, a "
                   "differing complement status being rank, preprocessing, covariate-set or cohort dependence; a "
                   "complement non-rejection is uninformative. Primaries: " + _context(primaries), y=0.0)
    values = {
        "n_analyses": len(table),
        "n_drawn": len(available),
        "not_available": not_available,
        "primaries": {
            f"{t['cell_type']} {t['cohort']}": {
                "rank": int(t["rank"]),
                "dims": {part: int(t[f"{part}_dim"]) for part in SPAN_SPLIT_PARTS},
                "f": {part: float(t[f"{part}_f"]) for part in SPAN_SPLIT_PARTS},
                "p": {part: float(t[f"{part}_p"]) for part in SPAN_SPLIT_PARTS},
                "share": None if pd.isna(t["share"]) else float(t["share"]),
                "complement_rejected": bool(t["complement_rejected"]),
                "span_rejected": bool(t["span_rejected"]),
                "production_rejected": bool(t["production_rejected"]),
                "span_concordant": bool(t["span_concordant"]),
            }
            for _, t in prim.iterrows()
        },
        "complement_rejected": [f"{t['cell_type']} {t['label']} {t['cohort']} p {_p(t['complement_p'])}"
                                for _, t in rejected.iterrows()],
        "discordant": [f"{t['cell_type']} {t['label']} {t['cohort']}" for _, t in available.iterrows()
                       if not bool(t["span_concordant"])],
        "share_range": (float(pd.to_numeric(available["share"]).min()), float(pd.to_numeric(available["share"]).max())),
    }
    return fig, values


def _caption_span_split(v: dict) -> str:
    parts = [f"{v['n_drawn']} of {v['n_analyses']} analyses carry the split test (not available: "
             f"{', '.join(v['not_available']) or 'none'})."]
    for name, e in v["primaries"].items():
        dims = ", ".join(f"{part} {e['dims'][part]:,}" for part in SPAN_SPLIT_PARTS)
        ps = ", ".join(f"{part} {e['p'][part]:.3f} (F {e['f'][part]:.4g}; {SPAN_SPLIT_ROLES[part]})"
                       for part in SPAN_SPLIT_PARTS)
        share = "n/a" if e["share"] is None else f"{e['share']:.4f}"
        parts.append(
            f"{name} (rank {e['rank']}; dimensions {dims}): pseudo-F p {ps}; complement "
            f"{'rejected' if e['complement_rejected'] else 'not rejected'}; span cell "
            f"{'concordant' if e['span_concordant'] else 'discordant'} (span rejected "
            f"{'yes' if e['span_rejected'] else 'no'}, any of delta, angle, shape rejected "
            f"{'yes' if e['production_rejected'] else 'no'}); in-span share {share}."
        )
    parts.append(f"Complement p <= {ANALYSIS_ALPHA}: {'; '.join(v['complement_rejected']) or 'none'}. Discordant span "
                 f"cells: {'; '.join(v['discordant']) or 'none'}. In-span share over the analyses drawn "
                 f"{v['share_range'][0]:.4f}–{v['share_range'][1]:.4f}.")
    return " ".join(parts)


def _scan_reference(a: dict) -> tuple[str, str] | None:
    """The (family, label) a registered analysis's scan equals by construction, or ``None``.

    A forced rank equals its cell type and cohort's primary (the scan does not
    read the latent rank); a latent arm with a latent adjustment equals the arm
    without one (the scan reads neither the latent basis nor the bootstrap).
    """
    if a["family"] == "forced rank":
        return ("primary", "primary")
    if a["family"] == "latent" and LATENT_SENSITIVITY_ARMS[a["label"]].get("latent_adjustment") is not None:
        plain = [arm for arm, spec in LATENT_SENSITIVITY_ARMS.items() if spec.get("latent_adjustment") is None]
        return ("latent", plain[0])
    return None


def _scan_table(r: FigureReader, directory: Path, record: dict) -> pd.DataFrame:
    """An analysis's ``interaction_scan.csv``, bound to its record (rows and hits) or the run stops."""
    path = directory / INTERACTION_SCAN_TABLE
    table = r.csv(path, dtype={"feature": str, "block": str})
    scan = record["interaction_scan"]
    hits = table[(~table["degenerate"].astype(bool)) & (table["p_value_adjusted"] <= ANALYSIS_ALPHA)]
    if len(table) != scan["n_features"] or len(hits) != len(scan.get("hits", [])) or set(
            zip(hits["block"], hits["feature"])) != {(h["block"], h["feature"]) for h in scan.get("hits", [])}:
        raise ValueError(f"{path} has {len(table)} rows and {len(hits)} hit(s) but {directory / _RESULT} records "
                         f"{scan['n_features']} features and {len(scan.get('hits', []))} hit(s); nothing was written")
    return table


def _draw_interaction_scan(r: FigureReader) -> tuple[Figure, dict]:
    rows, primaries, records, tables = [], [], {}, {}
    analyses = registered_analyses()
    for a in analyses:
        rec = r.result(a["dir"] / _RESULT)
        man = r.manifest(a["dir"] / _ANALYSIS_MANIFEST, a["cell_type"])
        status = _interaction_scan_status(rec, alpha=ANALYSIS_ALPHA)
        if status is None:
            raise ValueError(f"{a['dir'] / _RESULT} carries no interaction scan; run `sead-sexdiff {a['command']}`")
        records[(a["cell_type"], a["family"], a["label"], a["cohort"])] = rec
        rows.append({**{k: a[k] for k in ("cell_type", "family", "label", "cohort")},
                     "rank": rec["latent"]["selected_lv"], "p": status["p"], "rejected": status["rejected"],
                     **{f"hits_{b}": status["hits"][b] for b in ANALYSIS_BLOCKS},
                     "identical_to": _scan_reference(a)})
        if a["family"] == "primary":
            primaries.append((f"{a['cell_type']} {a['cohort']}", rec, man))
            tables[(a["cell_type"], a["cohort"])] = _scan_table(r, a["dir"], rec)
    # The identities the construction implies, reproduced from the records
    for row in rows:
        ref = row["identical_to"]
        if ref is None:
            continue
        a_rec = records[(row["cell_type"], row["family"], row["label"], row["cohort"])]
        b_rec = records[(row["cell_type"], *ref, row["cohort"])]
        if a_rec["interaction_scan"] != b_rec["interaction_scan"]:
            raise ValueError(f"the {row['cell_type']} {row['label']} {row['cohort']} scan differs from its reference "
                             f"({ref[1]}); the scan reads the features and models only; nothing was written")
    table = pd.DataFrame(rows)
    keys = table[["cell_type", "family", "label"]].drop_duplicates().reset_index(drop=True)
    n = len(keys)
    fig = Figure(figsize=(12, 0.27 * n + 2.2))
    gs = fig.add_gridspec(2, 2, width_ratios=(1, 1.2), height_ratios=(1.45, 1))
    ax_a, ax_c, ax_b = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[:, 1])
    marker = {b: m for b, m in zip(ANALYSIS_BLOCKS, ("o", "^"))}
    offsets = {b: (i - (len(ANALYSIS_BLOCKS) - 1) / 2) * 0.3 for i, b in enumerate(ANALYSIS_BLOCKS)}

    # (a) the primaries: the smallest adjusted p-values per block, hits filled
    prim = table[table["family"] == "primary"].reset_index(drop=True)
    _split_p_axis(ax_a)
    smallest = {}
    for y, t in prim.iterrows():
        color = _CELL[t["cell_type"]]
        st = tables[(t["cell_type"], t["cohort"])]
        st = st[~st["degenerate"].astype(bool)]
        for b in ANALYSIS_BLOCKS:
            top = st[st["block"] == b].nsmallest(FIGURE_SCAN_TOP_N, "p_value_adjusted", keep="first")
            smallest[(t["cell_type"], t["cohort"], b)] = float(top["p_value_adjusted"].min())
            hit = top["p_value_adjusted"] <= ANALYSIS_ALPHA
            jitter = np.linspace(-0.08, 0.08, len(top)) if len(top) > 1 else np.zeros(len(top))
            ax_a.scatter(top.loc[~hit, "p_value_adjusted"], y + offsets[b] + jitter[~hit.to_numpy()], s=14,
                         marker=marker[b], facecolor="white", edgecolor=color, linewidth=0.9, zorder=3)
            ax_a.scatter(top.loc[hit, "p_value_adjusted"], y + offsets[b] + jitter[hit.to_numpy()], s=22,
                         marker=marker[b], color=color, edgecolor="white", linewidth=0.8, zorder=4)
        ax_a.annotate(f"global p {_p(t['p'])}", (1.0, y), xytext=(8, 0), textcoords="offset points", va="center",
                      fontsize=6.6, color=_INK2, annotation_clip=False)
    ax_a.set_ylim(len(prim) - 0.5, -0.5)
    ax_a.set_yticks(range(len(prim)), [f"{t['cell_type']} · {t['cohort']}" for _, t in prim.iterrows()])
    ax_a.set_title(f"(a) primaries: {FIGURE_SCAN_TOP_N} smallest adjusted p per block")

    # (b) the global p-value of every registered analysis
    _split_p_axis(ax_b)
    for y, key in keys.iterrows():
        sub = table[(table[["cell_type", "family", "label"]] == key.values).all(axis=1)].set_index("cohort")
        color = _CELL[key["cell_type"]]
        for i, cohort in enumerate(EXPORT_COHORTS):
            t = sub.loc[cohort]
            yy = y + (i - (len(EXPORT_COHORTS) - 1) / 2) * 0.3
            if t["identical_to"] is not None:
                ax_b.scatter(t["p"], yy, s=40, marker="|", color=_MUTED, linewidth=1.4, zorder=3)
                continue
            face = color if cohort == EXPORT_COHORTS[0] else "white"
            ax_b.scatter(t["p"], yy, s=26, marker="o", facecolor=face, edgecolor=color, linewidth=1.3, zorder=3)
    labels = []
    for _, key in keys.iterrows():
        sub = table[(table[["cell_type", "family", "label"]] == key.values).all(axis=1)]
        ref = sub["identical_to"].iloc[0]
        suffix = "" if ref is None else f"  (= {ref[1]})"
        labels.append(f"{key['cell_type']} · {key['label']}{suffix}")
    ax_b.set_yticks(range(n), labels)
    ax_b.set_ylim(n - 0.5, -0.5)
    previous = None
    for y, key in keys.iterrows():
        group = (key["cell_type"], key["family"])
        if previous is not None and group != previous:
            ax_b.axhline(y - 0.5, color=_AXIS, linewidth=0.6)
        previous = group
    ax_b.set_title("(b) global scan p per analysis (grey ticks: equal to their reference by construction)")

    # (c) the null-maximum share per block beside the block's share of features, no threshold
    shares = {}
    for y, (label, rec, _) in enumerate(primaries):
        color = _CELL[rec_cell(label)]
        scan = rec["interaction_scan"]
        for b in ANALYSIS_BLOCKS:
            share = scan["blocks"][b]["null_max_share"]
            feat = scan["blocks"][b]["n_features"] / scan["n_features"]
            shares[(label, b)] = (share, feat)
            ax_c.scatter(share, y + offsets[b], s=30, marker=marker[b], color=color, edgecolor="white", linewidth=0.8,
                         zorder=3)
            ax_c.scatter(feat, y + offsets[b], s=40, marker="|", color=_INK2, linewidth=1.2, zorder=2)
    ax_c.set_xlim(0, 1)
    ax_c.set_ylim(len(primaries) - 0.5, -0.5)
    ax_c.set_yticks(range(len(primaries)), [label.replace(" ", " · ", 1) for label, _, _ in primaries])
    ax_c.set_xlabel("share of the permutation draws' maximum F (marker) and of the features (tick), per block")
    ax_c.set_title("(c) null-maximum share per block (not read)")

    legend = [
        *[Line2D([], [], marker=marker[b], linestyle="", color=_INK2, markersize=6, label=f"{b} (a, c)")
          for b in ANALYSIS_BLOCKS],
        Line2D([], [], marker="o", linestyle="", markerfacecolor="white", markeredgecolor=_INK2, markersize=6,
               label="not a hit (a)"),
        Line2D([], [], marker="o", linestyle="", color=_INK2, markersize=6, label=f"{EXPORT_COHORTS[0]} (b)"),
        *[Line2D([], [], marker="o", linestyle="", markerfacecolor="white", markeredgecolor=_INK2, markersize=6,
                 label=f"{c} (b)") for c in EXPORT_COHORTS[1:]],
        *[Line2D([], [], marker="o", linestyle="", color=_CELL[ct], markersize=6, label=ct) for ct in CELL_TYPES],
    ]
    distinct = table[table["identical_to"].isna()]
    rejected = distinct[distinct["rejected"].astype(bool)]
    _suptitle(fig, f"Per-feature interaction scan: global p ≤ {ANALYSIS_ALPHA} in {len(rejected)} of {len(distinct)} "
                   f"distinct analyses ({len(table)} drawn); only each cell type's {FIGURE_COHORT} primary is "
                   "confirmatory")
    fig.tight_layout(rect=(0, 0.05, 1, 0.965))
    fig.legend(handles=legend, loc="lower center", ncol=len(legend), bbox_to_anchor=(0.5, 0.012))
    power = INTERACTION_SCAN_POWER_STATEMENT
    _footnote(fig, "The global p-value is confirmatory (p ≤ α, unadjusted; each primary carries five tests read at α); "
                   "hits are features with step-down max-T adjusted p ≤ α (familywise over the scanned features), "
                   "detections, not drivers; F values, per-stage differences and null-maximum shares are reported and "
                   "not read; every analysis other than the primaries is sensitivity, a differing status being "
                   "preprocessing, covariate-set or cohort dependence. A non-rejection of the "
                   f"{power['cell_type']} {power['cohort']} primary carries MOTCO's measured power for a 20-feature "
                   f"ramp in one sex ({power['power']['3.0']:.3f} at 3 SDs, {power['power']['2.0']:.3f} at 2); every "
                   "other non-rejection is uninformative. Primaries: " + _context(primaries), y=0.0)

    # The scan concordance cells beside the stage DE, read for the caption (each made from the primary drawn)
    concordance = {}
    for cell_type in CELL_TYPES:
        d = cell_type_dirs(cell_type)
        for cohort in EXPORT_COHORTS:
            primary_dir = d["motco_analysis"] / cohort
            de_man_path = d["stage_de"] / cohort / STAGE_DE_MANIFEST
            de_man = r.manifest(de_man_path, cell_type)
            _check_made_from(de_man, de_man_path, primary_dir / _RESULT)
            li_man_path = d["stage_de_limma"] / cohort / STAGE_DE_LIMMA_MANIFEST
            li_man = r.manifest(li_man_path, cell_type)
            _check_made_from(li_man, li_man_path, primary_dir / _RESULT)
            for omic, cell in li_man["concordance"].items():
                scan = cell.get("scan")
                if scan is None:
                    raise ValueError(f"{li_man_path} carries no scan concordance cell for {omic}; run "
                                     "`sead-sexdiff stage-de-limma`")
                d_scan = de_man["concordance"][omic].get("scan")
                if d_scan is None or d_scan["cell"] != scan["deseq2_cell"]:
                    raise ValueError(f"{de_man_path} and {li_man_path} disagree on the {omic} DESeq2 scan cell")
                concordance[f"{cell_type} {cohort} {omic}"] = {
                    "hits": len(scan["hits"]), "deseq2": scan["deseq2_cell"], "limma": scan["limma_cell"],
                    "deseq2_n": cell["deseq2"]["n_interaction_exported"], "limma_n": cell["limma"]["n_interaction_exported"],
                }
    values = {
        "n_analyses": len(table),
        "n_distinct": len(distinct),
        "primaries": {
            label: {
                "p": float(rec["interaction_scan"][INTERACTION_SCAN_STATISTIC]),
                "rejected": bool(rec["interaction_scan"][INTERACTION_SCAN_STATISTIC] <= ANALYSIS_ALPHA),
                "max_f": float(rec["interaction_scan"]["max_f"]),
                "max_feature": rec["interaction_scan"]["max_feature"],
                "max_block": rec["interaction_scan"]["max_block"],
                "n_features": int(rec["interaction_scan"]["n_features"]),
                "df": (int(rec["interaction_scan"]["df_interaction"]), int(rec["interaction_scan"]["df_residual"])),
                "hits": {b: int(rec["interaction_scan"]["blocks"][b]["hits"]) for b in ANALYSIS_BLOCKS},
                "hit_names": [f"{h['block']}:{h['feature']}" for h in rec["interaction_scan"].get("hits", [])],
                "shares": {b: float(shares[(label, b)][0]) for b in ANALYSIS_BLOCKS},
            }
            for label, rec, _ in primaries
        },
        "rejected": [f"{t['cell_type']} {t['label']} {t['cohort']} p {_p(t['p'])}" for _, t in rejected.iterrows()],
        "concordance": concordance,
    }
    return fig, values


def rec_cell(label: str) -> str:
    """The cell type of a ``"<cell type> <cohort>"`` label."""
    return label.split(" ", 1)[0]


def _caption_interaction_scan(v: dict) -> str:
    parts = [f"{v['n_analyses']} analyses carry the scan, {v['n_distinct']} of them distinct (the forced ranks and the "
             "latent-adjusted arms equal their reference by construction)."]
    for name, e in v["primaries"].items():
        hits = ", ".join(f"{b} {e['hits'][b]}" for b in ANALYSIS_BLOCKS)
        shares = ", ".join(f"{b} {e['shares'][b]:.4f}" for b in ANALYSIS_BLOCKS)
        names = f" ({'; '.join(e['hit_names'])})" if e["hit_names"] else ""
        parts.append(
            f"{name} ({e['n_features']:,} features, df {e['df'][0]}, {e['df'][1]}): global p {e['p']:.3f}, "
            f"{'rejected' if e['rejected'] else 'not rejected'}; maximum F {e['max_f']:.4g} ({e['max_feature']}, "
            f"{e['max_block']}); hits {hits}{names}; null-maximum shares {shares}."
        )
    cells = "; ".join(
        f"{k}: {c['hits']} hit(s), DESeq2 {c['deseq2_n']} ({c['deseq2']}), limma {c['limma_n']} ({c['limma']})"
        for k, c in v["concordance"].items()
    )
    parts.append(f"Global p <= {ANALYSIS_ALPHA}: {'; '.join(v['rejected']) or 'none'}. Scan concordance cells "
                 f"(scan hits against the exported interaction counts): {cells}.")
    return " ".join(parts)


# ---------------------------------------------------------------------------
# The registry. Entries are added or changed by a change proposal only; a
# change that adds a step adds its figure here. Flag and ask before changing.
# ---------------------------------------------------------------------------

FIGURES: dict[str, dict] = {
    "cohorts": {
        "number": 1,
        "title": "Donor cohorts",
        "draw": _draw_cohorts,
        "caption": _caption_cohorts,
        "construction": (
            "Left: nuclei of the cell type per donor in RNA against ATAC (log scales), from each cell type's donor "
            "table, with the cohort thresholds in both omics; hollow markers lie outside cohort_30. Right: donors per "
            "sex x ADNC stage cell in each cohort; the cells are checked against the pseudobulk manifest."
        ),
        "reading": ("display_only",),
    },
    "trajectories": {
        "number": 2,
        "title": "Sex x stage trajectories",
        "draw": _draw_trajectories,
        "caption": _caption_trajectories,
        "construction": (
            "Per cell type and cohort, the primary analysis in its first two latent components: MOTCO's LS-means per "
            "sex and stage (stage_means.csv, joined by arrows in stage order) and the covariate-adjusted donor scores "
            "(raw PLS score minus the covariate contribution at the LS-mean reference row, whose cell means are checked "
            "against stage_means.csv); the share of each sex's net displacement in the plane is read from "
            "geometry_shares.csv. The panel titles give each statistic with its p-value and the CV-selected rank; the "
            "footnote gives the covariate block, the feature filtering and the RRPP null of the records."
        ),
        "reading": ("display_only", "motco_binding", "df_size"),
    },
    "components_by_stage": {
        "number": 3,
        "title": "Latent components by stage",
        "draw": _draw_components,
        "caption": _caption_components,
        "construction": (
            "For cohort_30, every latent component of each cell type's primary analysis against ADNC stage: "
            "covariate-adjusted donor scores (jittered) and the LS-mean trajectory per sex, with the shares of each "
            "sex's net displacement and of the centred sex difference along the component (geometry_shares.csv of "
            "axis-diagnostics)."
        ),
        "reading": ("display_only", "motco_binding", "diagnostic"),
    },
    "statistics_vs_null": {
        "number": 4,
        "title": "Statistics against their RRPP null",
        "draw": _draw_nulls,
        "caption": _caption_nulls,
        "construction": (
            "Per primary analysis (both cell types and cohorts), each observed statistic against the summary of its "
            "RRPP null from result.json (band to the 95th percentile, thinner band to the 99th, ticks at the median and "
            "the 95th percentile); angle carries the CV-selected rank, the pooled configuration eigengap and the 5th to "
            "95th percentile of its permutation null."
        ),
        "reading": ("display_only", "motco_binding", "df_size"),
    },
    "pvalue_overview": {
        "number": 5,
        "title": "Every MOTCO analysis",
        "draw": _draw_pvalue_overview,
        "caption": _caption_pvalue_overview,
        "construction": (
            "Every MOTCO analysis of the registries outside results/superseded/ (primary, forced ranks, preprocessing "
            "variants, latent-sensitivity arms, donor-exclusion variants; per cohort), one row per analysis with the "
            "latent rank per cohort, and its RRPP p-value per statistic (filled: cohort_30, hollow: cohort_50); the "
            "latent_adjusted angle and shape markers are grey (reported, not read)."
        ),
        "reading": ("display_only", "motco_binding", "df_size", "sensitivity_dependence", "forced_rank",
                    "latent_adjusted_not_read"),
    },
    "rank_sensitivity": {
        "number": 6,
        "title": "Rank sensitivity",
        "draw": _draw_rank_sensitivity,
        "caption": _caption_rank_sensitivity,
        "construction": (
            "For astrocytes (rank-sensitivity runs on astrocytes only), per cohort and statistic, the observed value at "
            "each forced rank of ANALYSIS_FORCED_RANKS against its RRPP null (median, 95th and 99th percentile), the "
            "CV-selected analysis at its rank (diamond), and the forced-rank p-values along the top."
        ),
        "reading": ("display_only", "motco_binding", "forced_rank", "df_size"),
    },
    "axis_correlations": {
        "number": 7,
        "title": "Axis correlations",
        "draw": _draw_axis_heatmap,
        "caption": _caption_axis_heatmap,
        "construction": (
            "For cohort_30, the Spearman correlation of every latent component of each cell type's primary analysis "
            "with every diagnostic variable (axis_correlations.csv of axis-diagnostics: design, analysis covariates, "
            "the export's TMM factors and library sizes, the technical table), grouped by source; values are printed "
            "where |rho| reaches the reporting flag DIAGNOSTIC_CORRELATION_THRESHOLD."
        ),
        "reading": ("display_only", "measurement_space", "diagnostic"),
    },
    "first_axis": {
        "number": 8,
        "title": "The first latent axis",
        "draw": _draw_first_axis,
        "caption": _caption_first_axis,
        "construction": (
            "For cohort_30, the raw pls_0 scores of each cell type's primary analysis against ADNC stage (jittered, with "
            "the median per stage) and against the technical variable picked by E3's rule (the largest |Spearman| among "
            "the technical table's columns and the export's TMM factors and library sizes); both Spearman correlations "
            "are recomputed and checked against axis_correlations.csv."
        ),
        "reading": ("display_only", "measurement_space", "diagnostic"),
    },
    "positive_control": {
        "number": 9,
        "title": "Positive control",
        "draw": _draw_positive_control,
        "caption": _caption_positive_control,
        "construction": (
            "The positive control's expectations E1 (exported RNA stage counts at FDR 0.05 under DESeq2 and limma), E2 "
            "(CV mean AUROC) and E3 (|Spearman| of the control's pls_0 with stage against the largest with a technical "
            "variable) per basis of POSITIVE_CONTROL_BASES and cohort, control against reference, from "
            "positive_control.csv, with each expectation's and the overall outcome as recorded."
        ),
        "reading": ("display_only", "positive_control"),
    },
    "stage_de": {
        "number": 10,
        "title": "Pseudobulk stage DE",
        "draw": _draw_stage_de,
        "caption": _caption_stage_de,
        "construction": (
            "For cohort_30, RNA, factor_additive stage_2_vs_0 (the contrast between the extreme stages): left, the "
            "DESeq2 volcano over the abundance-passing features (exported features darker, padj <= 0.05 in the cell "
            "type's colour, the eight smallest exported padj labelled); right, the DESeq2 log2 fold change against "
            "limma-voom's logFC over the exported features, by significance under each method. The counts are checked "
            "against stage_de_summary.csv and method_comparison.csv, a missing DESeq2 padj counting as not significant."
        ),
        "reading": ("display_only", "stage_de_counts"),
    },
    "span_split": {
        "number": 11,
        "title": "Split interaction test",
        "draw": _draw_span_split,
        "caption": _caption_span_split,
        "construction": (
            "MOTCO's split interaction test from each record's span_split block: (a) per primary analysis (both cell "
            "types and cohorts) the pseudo-F p-value of each part on a log axis with alpha, the complement a filled "
            "circle (confirmatory), the span a hollow diamond with its concordance cell beside the production "
            "statistics, total a grey square (descriptive); (b) the complement p-value of every registered analysis that carries the test, enumerated "
            "as in figure 5 (filled: cohort_30, hollow: cohort_50), the latent_adjusted arms listed as not available; "
            "(c) the in-span share of the same analyses against their fitted rank, with no reference line: the share "
            "is descriptive, grows with the fitted rank, and is read against no null here."
        ),
        "reading": ("display_only", "span_split", "sensitivity_dependence", "forced_rank"),
    },
    "interaction_scan": {
        "number": 12,
        "title": "Per-feature interaction scan",
        "draw": _draw_interaction_scan,
        "caption": _caption_interaction_scan,
        "construction": (
            "MOTCO's per-feature interaction scan from each record's interaction_scan block and each primary's "
            "interaction_scan.csv: (a) per primary analysis (both cell types and cohorts) the "
            f"{FIGURE_SCAN_TOP_N} smallest step-down max-T adjusted p-values per block on a log axis with alpha, "
            "hits filled, the global p-value beside; (b) the global p-value of every registered analysis, enumerated "
            "as in figure 5 (filled: cohort_30, hollow: cohort_50), the forced ranks and the latent-adjusted arms "
            "drawn as grey ticks because their scan equals their reference's by construction; (c) per primary, each "
            "block's share of the permutation draws' maximum F beside its share of the features, with no threshold: "
            "the shares are descriptive."
        ),
        "reading": ("display_only", "interaction_scan", "sensitivity_dependence"),
    },
}
for _name, _entry in FIGURES.items():
    if not set(_entry["reading"]) <= set(FIGURE_READING_RULES):
        raise AssertionError(f"FIGURES[{_name!r}] names reading rules outside FIGURE_READING_RULES")
if [e["number"] for e in FIGURES.values()] != list(range(1, len(FIGURES) + 1)):
    raise AssertionError("FIGURES must be numbered 1..n in registry order")
del _name, _entry


def figure_stem(name: str) -> str:
    """The file stem of a registered figure, ``<nn>_<name>``."""
    return f"{FIGURES[name]['number']:02d}_{name}"


# ---------------------------------------------------------------------------
# Files: determinism, rasterization, fingerprint (D6, D7)
# ---------------------------------------------------------------------------

def _module_sha256() -> str:
    return _sha256_file(Path(__file__))


def _package_version() -> str | None:
    try:
        return version("sead-sexdiff")
    except PackageNotFoundError:
        return None


def _check_font() -> dict[str, str]:
    """The font files the figures resolve, each required to be matplotlib's bundled DejaVu Sans."""
    data = Path(matplotlib.get_data_path()).resolve()
    out = {}
    for weight in ("normal", "semibold"):
        for style in ("normal", "italic"):
            path = Path(font_manager.findfont(font_manager.FontProperties(family=FIGURE_FONT, weight=weight, style=style),
                                              fallback_to_default=False)).resolve()
            if data not in path.parents:
                raise ValueError(f"{FIGURE_FONT} ({weight}, {style}) resolves to {path}, not matplotlib's bundled copy "
                                 f"under {data}; the figures use the bundled font only")
            out[f"{weight}/{style}"] = path.relative_to(data).as_posix()
    return out


def fingerprint(stem: str, inputs: dict[str, str], module_sha256: str, matplotlib_version: str) -> str:
    """The sha256 of the canonical JSON of the figure, its inputs, the module's sha256 and the matplotlib version."""
    payload = {"figure": stem, "inputs": sorted([path, sha] for path, sha in inputs.items()),
               "module_sha256": module_sha256, "matplotlib": matplotlib_version}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _rasterize_dense(fig: Figure) -> int:
    n = 0
    for ax in fig.axes:
        for coll in ax.collections:
            if isinstance(coll, PathCollection) and len(coll.get_offsets()) > FIGURE_RASTER_MIN_POINTS:
                coll.set_rasterized(True)
                n += 1
    return n


def _metadata(fmt: str, stem: str, fp: str) -> dict:
    description = f"sead-sexdiff figure {stem}; fingerprint sha256:{fp}"
    if fmt == "png":
        return {"Title": stem, "Description": description}
    if fmt == "pdf":
        return {"Title": stem, "Subject": description, "CreationDate": None}
    return {"Title": stem, "Description": description, "Date": None}


def _save(fig: Figure, out_dir: Path, stem: str, fp: str) -> dict[str, dict]:
    _rasterize_dense(fig)
    files = {}
    for fmt in FIGURE_FORMATS:
        path = out_dir / f"{stem}.{fmt}"
        fig.savefig(path, format=fmt, dpi=FIGURE_DPI, metadata=_metadata(fmt, stem, fp))
        files[fmt] = {"file": path.name, "sha256": _sha256_file(path)}
    return files


# ---------------------------------------------------------------------------
# Captions, manifest, staleness, run modes (D5, D8, D9)
# ---------------------------------------------------------------------------

def _caption(name: str, values: dict) -> dict:
    entry = FIGURES[name]
    return {
        "construction": entry["construction"],
        "values": entry["caption"](values),
        "reading": [{"rule": key, "text": FIGURE_READING_RULES[key]} for key in entry["reading"]],
    }


def render_captions(figures: dict[str, dict]) -> str:
    """``figures.md``: per figure its heading, the SVG link and the three caption parts as paragraphs."""
    lines = [
        "# Case-study figures",
        "",
        "Drawn by `sead-sexdiff figures` from the generated outputs listed in `figures_manifest.json`. Each caption "
        "has three parts: how the figure is made, the values read from the files when it was drawn, and the "
        "registered reading rules that apply to it.",
        "",
    ]
    for stem, f in figures.items():
        c = f["caption"]
        lines += [
            f"## Figure {f['number']}. {f['title']}",
            "",
            f"![{f['title']}]({stem}.{FIGURE_FROZEN_FORMAT})",
            "",
            f"*Construction.* {c['construction']}",
            "",
            f"*Values.* {c['values']}",
            "",
            "*Reading.* " + " ".join(rule["text"] for rule in c["reading"]),
            "",
        ]
    return "\n".join(lines)


def _jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, Path):
        return obj.as_posix()
    if isinstance(obj, np.generic):
        return obj.item()
    return obj


def _display_constants(fonts: dict[str, str]) -> dict:
    return _jsonable({
        "cohort": FIGURE_COHORT,
        "stage_de": FIGURE_STAGE_DE,
        "first_axis_rule": FIGURE_FIRST_AXIS_RULE,
        "technical_sources": POSITIVE_CONTROL_TECHNICAL_SOURCES,
        "correlation_threshold": DIAGNOSTIC_CORRELATION_THRESHOLD,
        "sensitivity_cell_type": FIGURE_SENSITIVITY_CELL_TYPE,
        "preprocessing_root": FIGURE_PREPROCESSING_ROOT,
        "palette": FIGURE_PALETTE,
        "dpi": FIGURE_DPI,
        "formats": FIGURE_FORMATS,
        "raster_min_points": FIGURE_RASTER_MIN_POINTS,
        "svg_hashsalt": FIGURE_SVG_HASHSALT,
        "font": {"family": FIGURE_FONT, "files": fonts},
        "jitter_seeds": FIGURE_JITTER_SEEDS,
        "tolerances": {"stage_means_atol": FIGURE_STAGE_MEANS_ATOL, "spearman_atol": FIGURE_SPEARMAN_ATOL},
    })


def draw_figures(out_dir: str | Path, *, consistency: dict, reader: FigureReader) -> dict:
    """Draw the whole set into ``<out_dir>.tmp/`` and rename it into place (design D8); return the manifest.

    The manifest is written last; the previous set becomes ``<out_dir>.old/``
    and is removed only after the rename succeeded. On any error the
    temporary directory is removed and the previous set is left unchanged.
    """
    out_dir = Path(out_dir)
    tmp = out_dir.with_name(out_dir.name + ".tmp")
    old = out_dir.with_name(out_dir.name + ".old")
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)
    module_sha, mpl = _module_sha256(), matplotlib.__version__
    try:
        fonts = _check_font()
        figures: dict[str, dict] = {}
        for name, entry in FIGURES.items():
            stem = figure_stem(name)
            reader.begin(stem)
            with matplotlib.rc_context(FIGURE_RCPARAMS):
                fig, values = entry["draw"](reader)
                reader.end()
                inputs = reader.inputs[stem]
                fp = fingerprint(stem, inputs, module_sha, mpl)
                files = _save(fig, tmp, stem, fp)
            figures[stem] = {
                "number": entry["number"],
                "name": name,
                "title": entry["title"],
                "inputs": [{"path": p, "sha256": s} for p, s in sorted(inputs.items())],
                "manifests": list(reader.manifests[stem]),
                "analyses": [{"result": p, **m} for p, m in sorted(reader.analyses[stem].items())],
                "fingerprint": fp,
                "caption": _caption(name, values),
                "files": files,
                "findings": list(reader.findings[stem]),
            }
            print(f"[figures] Drew {stem} ({', '.join(FIGURE_FORMATS)}; {len(inputs)} inputs, "
                  f"{len(reader.manifests[stem])} manifests)")
        (tmp / FIGURE_CAPTIONS).write_text(render_captions(figures))
        manifest = {
            "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "command": "figures",
            "module": "sead_sexdiff/figures.py",
            "module_sha256": module_sha,
            "package_version": _package_version(),
            "git": _git_revision(),
            "software": {"matplotlib": mpl, "numpy": np.__version__, "pandas": pd.__version__},
            "display": _display_constants(fonts),
            "reading_rules": FIGURE_READING_RULES,
            "consistency": consistency,
            "captions": FIGURE_CAPTIONS,
            "figures": figures,
        }
        (tmp / FIGURE_MANIFEST).write_text(json.dumps(_jsonable(manifest), indent=2) + "\n")
    except BaseException:
        reader.end()
        shutil.rmtree(tmp, ignore_errors=True)
        raise
    if old.exists():
        shutil.rmtree(old)
    if out_dir.exists():
        os.replace(out_dir, old)
    os.replace(tmp, out_dir)
    if old.exists():
        shutil.rmtree(old)
    print(f"[figures] Wrote {out_dir} ({len(figures)} figures, {FIGURE_CAPTIONS}, {FIGURE_MANIFEST})")
    return manifest


def figures_status(out_dir: str | Path) -> dict:
    """Compare the working copy's manifest with the files on disk (spec: Staleness check); write nothing.

    A figure is stale when a file it read, or a file recorded with a sha256 by
    a manifest it read, hashes differently or is missing, when the module's
    sha256 or the matplotlib version differs, when it is not in the manifest,
    or when one of its images is missing or differs from its recorded sha256.

    Returns
    -------
    dict
        ``state`` (``fresh``, ``stale`` or ``absent``), ``exit_code`` (0, 1,
        2), ``stale`` (stem → reasons) and ``fresh`` (stems).
    """
    out_dir = Path(out_dir)
    path = out_dir / FIGURE_MANIFEST
    stems = [figure_stem(name) for name in FIGURES]
    if not path.is_file():
        return {"state": "absent", "exit_code": 2, "stale": {}, "fresh": [], "manifest": None}
    manifest = json.loads(path.read_text())
    reader = FigureReader()
    module_sha, mpl = _module_sha256(), matplotlib.__version__
    recorded_figures = manifest.get("figures", {})
    stale: dict[str, list[str]] = {}
    fresh: list[str] = []
    for stem in [s for s in recorded_figures if s not in stems]:
        stale[stem] = ["no longer in the registry FIGURES"]
    for stem in stems:
        f = recorded_figures.get(stem)
        if f is None:
            stale[stem] = ["not in the manifest"]
            continue
        reasons = []
        if manifest.get("module_sha256") != module_sha:
            reasons.append(f"module changed (sha256 {manifest.get('module_sha256')} when drawn, {module_sha} now)")
        if manifest.get("software", {}).get("matplotlib") != mpl:
            reasons.append(f"matplotlib changed ({manifest.get('software', {}).get('matplotlib')} when drawn, {mpl} now)")
        current: dict[str, str] = {}
        unchanged: set[str] = set()
        for entry in f["inputs"]:
            key = FigureReader.key(entry["path"])
            if not Path(key).is_file():
                reasons.append(f"input missing: {key}")
                continue
            current[key] = reader.sha256(key)
            if current[key] != entry["sha256"]:
                reasons.append(f"input changed: {key}")
            else:
                unchanged.add(key)
        for mp in f["manifests"]:
            if mp not in unchanged:
                continue
            for recorded, sha in _manifest_entries(reader.json(mp)):
                try:
                    key = FigureReader.key(recorded)
                except ValueError as exc:
                    reasons.append(f"recorded by {mp}: {exc}")
                    continue
                if not Path(key).is_file():
                    reasons.append(f"recorded by {mp}: {key} missing")
                elif reader.sha256(key) != sha:
                    reasons.append(f"recorded by {mp}: {key} changed")
        for img in f["files"].values():
            p = out_dir / img["file"]
            if not p.is_file():
                reasons.append(f"image {img['file']} missing")
            elif _sha256_file(p) != img["sha256"]:
                reasons.append(f"image {img['file']} differs from its recorded sha256")
        if not reasons and fingerprint(stem, current, module_sha, mpl) != f["fingerprint"]:
            reasons.append("fingerprint differs from the recorded one")
        if reasons:
            stale[stem] = list(dict.fromkeys(reasons))
        else:
            fresh.append(stem)
    state = "stale" if stale else "fresh"
    return {"state": state, "exit_code": 1 if stale else 0, "stale": stale, "fresh": fresh, "manifest": manifest}


def _print_status(out_dir: Path, status: dict) -> None:
    if status["state"] == "absent":
        print(f"[figures] {out_dir}: no {FIGURE_MANIFEST} (no working copy)")
        return
    print(f"[figures] {out_dir}: {len(status['fresh'])} fresh, {len(status['stale'])} stale")
    for stem, reasons in status["stale"].items():
        shown = reasons[:6]
        more = f"; ... and {len(reasons) - len(shown)} more" if len(reasons) > len(shown) else ""
        print(f"[figures] STALE {stem}: " + "; ".join(shown) + more)


def report_destination(stem: str, reports_dir: str | Path = FIGURE_REPORTS_DIR) -> Path:
    """``<reports_dir>/figures/<stem>/``; raise for an invalid stem or an existing destination."""
    if not FIGURE_REPORT_STEM_RE.match(stem):
        raise ValueError(f"report stem {stem!r} must match {FIGURE_REPORT_STEM_RE.pattern}")
    dest = Path(reports_dir) / FIGURE_REPORT_SUBDIR / stem
    if dest.exists():
        raise ValueError(f"{dest} exists; a frozen set is never regenerated or overwritten (nothing was copied)")
    return dest


def freeze_figures(out_dir: str | Path, dest: str | Path) -> Path:
    """Copy the SVGs, ``figures.md`` and the manifest to ``dest`` through ``<dest>.tmp/``, each SVG hash-checked."""
    out_dir, dest = Path(out_dir), Path(dest)
    manifest = json.loads((out_dir / FIGURE_MANIFEST).read_text())
    tmp = dest.with_name(dest.name + ".tmp")
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)
    try:
        for stem, f in manifest["figures"].items():
            img = f["files"][FIGURE_FROZEN_FORMAT]
            shutil.copyfile(out_dir / img["file"], tmp / img["file"])
            sha = _sha256_file(tmp / img["file"])
            if sha != img["sha256"]:
                raise ValueError(f"{tmp / img['file']} hashes {sha}; the manifest records {img['sha256']}")
        for name in (FIGURE_CAPTIONS, FIGURE_MANIFEST):
            shutil.copyfile(out_dir / name, tmp / name)
        if dest.exists():
            raise ValueError(f"{dest} exists; a frozen set is never overwritten")
        os.replace(tmp, dest)
    except BaseException:
        shutil.rmtree(tmp, ignore_errors=True)
        raise
    print(f"[figures] Froze {len(manifest['figures'])} {FIGURE_FROZEN_FORMAT.upper()}s, {FIGURE_CAPTIONS} and "
          f"{FIGURE_MANIFEST} into {dest}")
    return dest


def run_figures(
    *,
    out_dir: str | Path = FIGURE_OUT_DIR,
    check: bool = False,
    report: str | None = None,
    reports_dir: str | Path = FIGURE_REPORTS_DIR,
) -> dict:
    """Draw, check or freeze the registered case-study figures (``sead-sexdiff figures``).

    Parameters
    ----------
    out_dir:
        The working copy (default ``results/figures``). Inputs are read from
        ``results/`` relative to the working directory.
    check:
        Only compare the working copy with the files on disk (exit code 0
        fresh, 1 stale, 2 no manifest); nothing is written.
    report:
        A report stem: run the consistency check, redraw the working copy if
        it is stale or absent, then freeze its SVGs, ``figures.md`` and the
        manifest under ``<reports_dir>/figures/<report>/`` (never overwritten).
    reports_dir:
        The reports root (default ``reports``).

    Returns
    -------
    dict
        ``exit_code``, ``status`` (the staleness check before any drawing),
        ``drawn`` (whether the set was redrawn), ``consistency`` (the
        pre-pass summary when it ran) and ``frozen`` (the frozen directory).

    Notes
    -----
    A plain run does nothing when every figure is fresh; otherwise it runs the
    consistency check and redraws the whole set (design D8). The working copy
    is redrawn freely: the record of a figure is the frozen set a report
    carries.
    """
    out_dir = Path(out_dir)
    if check and report is not None:
        raise ValueError("--check and --report are exclusive")
    if check:
        status = figures_status(out_dir)
        _print_status(out_dir, status)
        return {"exit_code": status["exit_code"], "status": status, "drawn": False, "consistency": None, "frozen": None}
    dest = report_destination(report, reports_dir) if report is not None else None
    status = figures_status(out_dir)
    _print_status(out_dir, status)
    consistency = None
    reader = FigureReader()
    if status["state"] != "fresh" or dest is not None:
        consistency = check_consistency(reader)
        print(f"[figures] Consistency check passed: {consistency['manifests']} manifests "
              f"({consistency['manifests_from_registries']} from the registries, {consistency['manifests_found']} more "
              f"under the registered directories), {consistency['entries']} recorded entries "
              f"({consistency['distinct_files']} files) hashing as recorded; {consistency['analyses']} analyses on "
              f"MOTCO {consistency['motco']}, residual correction {ANALYSIS_RESIDUAL_CORRECTION}, rank floor "
              f"{ANALYSIS_RANK_FLOOR}, {ANALYSIS_PERMUTATIONS} permutations, seed {ANALYSIS_SEED}, split test "
              f"(not available in {', '.join(LATENT_SENSITIVITY_SPLIT_UNAVAILABLE)})")
        if consistency["findings"]:
            print(f"[figures] FINDING: {len(consistency['findings'])} manifests record no cell_type (written before "
                  "the cell-type registry; read as their root's cell type); listed in the manifest")
    drawn = False
    if status["state"] == "fresh":
        print(f"[figures] The working copy is fresh; nothing drawn")
    else:
        draw_figures(out_dir, consistency=consistency, reader=reader)
        drawn = True
    frozen = freeze_figures(out_dir, dest) if dest is not None else None
    return {"exit_code": 0, "status": status, "drawn": drawn, "consistency": consistency, "frozen": frozen}

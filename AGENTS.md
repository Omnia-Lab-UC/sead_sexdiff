# AGENTS.md

Guidance for AI coding agents working in this repository. This file is the single
source of truth; tool-specific files (e.g. `CLAUDE.md`) only reference it.

## Project

Python package + CLI (`sead-sexdiff`) for sex-difference analyses on the public
SEA-AD single-nucleus data. Two pipelines share the `download` step.

Legacy PFC (prefrontal cortex, RNA, DESeq2):

1. `download` – fetch files from the public S3 bucket `sea-ad-single-cell-profiling`
   (anonymous/unsigned boto3 access) into the data root, preserving key
   structure; size-checked and idempotent, writes `download_manifest.json`
   under the prefix.
2. `subset` – read donor `.h5ad` files one at a time, filter by `obs['Subclass']`,
   apply basic QC, concatenate, and write `<Subclass>_subset.h5ad` plus a
   `decoupler` pseudobulk `<Subclass>_pseudobulk.h5ad` to `results/cleaned_files/`.
3. `de` – run pydeseq2 per `obs.Supertype` with design `~Sex + ADNC + Sex:ADNC`,
   writing a CSV and volcano PDF per supertype to `results/DE/`.

MTG astrocyte case study (RNA + ATAC, input to MOTCO):

1. `download --omic {rna,atac,both}` – fetch the MTG per-donor objects
   (`MTG/RNAseq/donors_objects/`, 84 donors, release 2026-06-22;
   `MTG/ATACseq/donor_objects/`, 86 donors, release 2024-12-06).
2. `pseudobulk` – per donor, sum raw counts over astrocyte nuclei into
   `results/pseudobulk/` (`rna_counts.h5ad`, `atac_counts.h5ad`,
   `donor_table.csv`, `supertype_composition.csv`, `cohort_summary.csv`,
   `manifest.json`), with a resumable per-donor cache in `per_donor/`.
   Counts stay raw: transforms and feature filters happen in `export`.
3. `covariates` – read donor-level `obs` columns (age, PMI, RIN, education)
   from the raw donor objects without loading matrices, code them
   numerically (`90+` → 90, flagged), cross-check RNA against ATAC, and write
   `results/pseudobulk/donor_covariates.csv` + `covariates_manifest.json`.
   Runs on the cluster login node (needs the raw objects; seconds per file).
4. `export` – from `results/pseudobulk/` only (any workstation), per cohort:
   drop chrX/Y/M → CPM ≥ 1 in ≥ k donors → TMM → log2-CPM → top 5,000 by
   variance per omic, writing `results/motco_input/<cohort>/` (`rna.csv`,
   `atac.csv`, `metadata.csv`, `covariates.csv`, `features_<omic>.csv`,
   `manifest.json`). Covariates are exported, not applied: adjustment is a
   MOTCO change.
5. `analyze` – from `results/motco_input/` only, per cohort: the
   covariate-adjusted sex × stage trajectory comparison through MOTCO's
   `analyze_trajectories` (≥ 0.8.0), writing `results/motco_analysis/<cohort>/`
   (MOTCO's `result.json`, `latent_scores.csv`, `pls_loadings.csv`,
   `stage_means.csv`, plus `analysis_metadata.csv`, `analysis_manifest.json`,
   `summary.md`). Never assemble it from `motco plsr` + `motco de`.

`README.md` documents CLI flags and the Python API in detail; keep it in sync
when behavior or defaults change.

## Layout

- `src/sead_sexdiff/cli.py` – argparse CLI (`build_parser`, `main`); thin layer
  that calls into `data.py` and prints `[command]`-prefixed progress messages.
- `src/sead_sexdiff/data.py` – all logic: S3 download, subsetting/QC,
  pseudobulk, differential expression, covariates, MOTCO export and analysis.
- `src/sead_sexdiff/DEG.R` – exploratory R script (MAST mixed model); not wired
  into the CLI.
- `scripts/` – `pseudobulk.sbatch` (cluster job), `sync_pseudobulk.sh`
  (rsync of `results/pseudobulk/` to a workstation) and `check_tmm_edger.R`
  (one-off edgeR cross-check of the export's TMM / log2-CPM).
- `data/`, `results/` – gitignored; never commit their contents.
- `reports/` – committed case-study reports (e.g. `motco-astrocyte-<date>.md`).
  Their numbers are pasted from generated outputs (`summary.md`, manifests),
  never retyped, and checked against `result.json` before committing.
- `openspec/` – OpenSpec change proposals/specs (see below).

## Environment & commands

- Python `>=3.12,<3.14`, managed with **uv** (`uv.lock` is committed).
- Setup: the Nix flake (`flake.nix`) provides Python 3.12, uv, R, and the R
  packages used by `DEG.R` and `check_tmm_edger.R`; its shell hook runs `uv sync --locked`. Enter it
  with direnv (`.envrc` → `use flake`) or `nix develop`. Python packages come
  from uv, not Nix; R packages come from Nix (add them to `flake.nix`).
- Run the CLI: `uv run sead-sexdiff <download|pseudobulk|covariates|export|analyze|subset|de> [options]`.
- `motco` is a git dependency pinned to a release tag of
  `https://github.com/tomszar/MOTCO` (`[tool.uv.sources]`; `uv.lock` holds
  the commit). Move it only with `uv add "motco @ git+https://github.com/tomszar/MOTCO@vX.Y.Z"`,
  never as a path dependency on `../MOTCO`. It pins scikit-learn exactly; do
  not loosen that pin to resolve a conflict. `run_analysis` refuses a MOTCO
  older than `ANALYSIS_MOTCO_MIN_VERSION`.
  Outside the Nix shell numpy cannot load its native libraries; run commands
  through `nix develop -c ...` there.
- Data root: explicit flag (`download --dest-root`, `--data-root` elsewhere)
  → `SEAD_DATA_ROOT` → `data/`. Machines differ only by this value. `export`
  uses it only for `<root>/annotation/` (GENCODE GTF, downloaded once).
- Add dependencies with `uv add <pkg>` so `pyproject.toml` and `uv.lock` stay
  consistent.
- There is no test suite or linter configured yet. Verify changes by running
  the relevant CLI command or importing the function; say so explicitly when
  something could not be verified (e.g. real data not available locally).

## Data & domain notes

- Donor `.h5ad` files are large. Process them one at a time and free
  intermediates; do not load the whole folder into memory at once.
- Expected `obs` columns: `Subclass`, `Supertype`, `Donor ID`, `Sex`, and
  `Overall AD neuropathological Change`, which is renamed to `ADNC` before
  writing subsets.
- Raw counts live in the `UMIs` layer; pseudobulk uses it with `mode='sum'`.
  In the MTG RNA objects `X` is normalized; in the MTG ATAC objects raw peak
  counts are in `X` (float32, integer-valued; no layers).
- QC defaults: `min_genes=200` per cell, `min_cells=3` per gene.
- MTG pseudobulk scientific choices (constants in `data.py`; flag and ask
  before changing any of them):
  - `TAXONOMY`: per-release label mapping keyed by the release date in the
    file name. Every `Subclass` value and every astrocyte `Supertype` must be
    listed, otherwise the run fails. RNA calls microglia `Immune`, ATAC
    `Microglia-PVM`; astrocytes are `Astrocyte` in both. Spike notes:
    `openspec/changes/mtg-astrocyte-donor-pseudobulks/spike.md`.
  - `MATRIX_SOURCE`, `RNA_MIN_GENES` (RNA nuclei ≥ 200 genes; ATAC no extra QC).
  - `Reference` ADNC donors excluded; `ADNC_STAGES` (`Not AD`+`Low` → `0`,
    `Intermediate` → `1`, `High` → `2`).
  - `COHORT_THRESHOLDS` (≥ 30 / ≥ 50 astrocyte nuclei in both omics) and
    `PREREGISTERED_COHORTS` (≥ 30: n = 80, F 11/10/28, M 9/10/12; ≥ 50:
    n = 74). Deviations go to `manifest.json` `findings`, never silently.
- Covariate coding (`COVARIATE_COLUMNS`, `COVARIATE_CROSS_CHECK`): RNA values
  are used; `N+` ages are capped and flagged; non-numeric values stay missing
  and are listed; nothing is imputed.
- MOTCO export constants (flag and ask before changing any of them):
  `EXPORT_EXCLUDED_CHROMOSOMES` (chrX/Y/M; sex is the group),
  `EXPORT_MIN_CPM` (1, in ≥ k donors, k = smallest sex × stage cell, computed
  per cohort), TMM (edgeR defaults `TMM_LOGRATIO_TRIM` / `TMM_SUM_TRIM`),
  `EXPORT_PRIOR_COUNT` (2), `EXPORT_TOP_N` (5,000 per block, equal N because
  MOTCO z-scores features without block weights), `EXPORT_COHORTS`,
  `STAGE_LABELS`, `GENCODE_RELEASE` (v32 basic, matching 10x GRCh38-2020-A).
  No step may use sex, stage or covariates. Changes to the numpy TMM must be
  re-checked with `scripts/check_tmm_edger.R`.
- MOTCO analysis constants (flag and ask before changing any of them):
  `ANALYSIS_COVARIATES` (age_at_death, pmi, rin, years_education,
  log10_n_rna, log10_n_atac; additive, in this order),
  `ANALYSIS_LOG10_COLUMNS` (nuclei counts enter on the log10 scale, derived
  at analysis time), `ANALYSIS_EXCLUDED_COVARIATES` (`education_level` by
  decision, `age_capped` constant in this release),
  `ANALYSIS_MISSING_COVARIATES` (`error`, never impute),
  `ANALYSIS_PERMUTATIONS` (999), `ANALYSIS_SEED` (20261019),
  `ANALYSIS_ALPHA` (0.05) with `ANALYSIS_REJECTION_RULE` (p ≤ α),
  `ANALYSIS_MOTCO_MIN_VERSION` (0.8.0). MOTCO's own defaults (net-displacement
  orientation, double-CV rank) are not overridden, and there is no
  unadjusted run. `--n-jobs` changes MOTCO's RRPP draws (one seed stream per
  worker), so report results from `n_jobs = 1`.
- Reporting MOTCO results follows the binding rules in
  `../MOTCO/docs/roadmap.md` (Phase 6, Next changes item 1): CV-selected rank
  beside `angle`; `angle`/`shape` non-rejections are uninformative; a
  `delta` rejection is a difference, not specifically of size; report the
  covariate block and the feature filtering beside the results.
- DE: supertypes whose design matrix is not full rank are skipped (pydeseq2
  `UserWarning`), not treated as errors. Output filenames replace `/` with `-`.
- Statistical choices (design formula, contrasts, QC thresholds, pseudobulk
  parameters) are scientific decisions. Don't change them silently; flag them
  and ask.

## Cluster (`cluster.ing.uc.cl`)

- Raw data is canonical on the cluster at `/home1/tgonza/sead_sexdiff/data`
  (repo clone at `/home1/tgonza/sead_sexdiff`, `uv sync --locked` on the login
  node; uv in `~/.local/bin`; clone over HTTPS, the cluster has no GitHub SSH
  key). Partition `512x1024`, no scratch. Avoid the `full` partition: node n14
  has glibc 2.17 and cannot import the numpy wheels installed from the login
  node (glibc 2.34).
- Compute nodes have no internet access: run `download --omic both` on the
  login node (nohup), then `scripts/pseudobulk.sbatch` runs
  `pseudobulk --omic both` as one sequential job; resubmitting resumes.
  `covariates --omic both` reads metadata only and runs directly on the login
  node.
- Workstations pull only `results/pseudobulk/` with
  `scripts/sync_pseudobulk.sh`; never copy raw donor objects around.

## Code conventions

- `from __future__ import annotations`, modern type hints (`str | Path`,
  `list[str]`), keyword-only args (`*`) for optional parameters.
- NumPy-style docstrings (`Parameters`, `Returns`, `Notes` sections).
- Accept `str | Path` for paths and convert with `Path(...)` internally.
- Use `print` with a `[command]` prefix for CLI progress output (no logging
  framework yet).

## OpenSpec workflow

Non-trivial features go through OpenSpec (`openspec/changes/`, specs in
`openspec/specs/`): propose → apply → archive. Agent skills/commands for this
live in `.claude/` and `.junie/`.

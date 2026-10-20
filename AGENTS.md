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
4. `technical-covariates` – per RNA donor object, select the astrocyte
   nuclei through the pseudobulk's code path (`select_astrocyte_nuclei`),
   read their raw `UMIs` rows and `obs` QC columns (never the normalized
   `X`) and write per donor the nucleus means and SDs and the UMI-weighted
   pseudobulk fractions of the predeclared gene sets (`TECHNICAL_GENE_SETS`:
   ribosomal, mitochondrial, hemoglobin), the within-donor supertype share
   (η²) of each fraction, the nucleus and library QC means
   (`TECHNICAL_OBS_COLUMNS`) and `n_libraries_rna`; from the ATAC objects
   the library-level `ATAC_*` metrics (`obs` only). Fractions are
   cross-checked against `rna_counts.h5ad` and nucleus counts against the
   pseudobulk cache and `donor_table.csv` (any mismatch fails). Writes
   `results/pseudobulk/donor_technical.csv`, `nucleus_technical_summary.csv`
   and `technical_manifest.json`, never `donor_covariates.csv` or
   `covariates_manifest.json` (hash-bound by the reference checks), with a
   resumable per-donor JSON cache in `per_donor_technical/`. Runs on the
   cluster (`scripts/technical.sbatch`, or nohup on the login node).
5. `export` – from `results/pseudobulk/` only (any workstation), per cohort:
   drop chrX/Y/M → CPM ≥ 1 in ≥ k donors → TMM → log2-CPM → top 5,000 by
   variance per omic, writing `results/motco_input/<cohort>/` (`rna.csv`,
   `atac.csv`, `metadata.csv`, `covariates.csv`, `features_<omic>.csv`,
   `manifest.json`). Covariates are exported, not applied: adjustment is a
   MOTCO change.
6. `analyze` – from `results/motco_input/` only, per cohort: the
   covariate-adjusted sex × stage trajectory comparison through MOTCO's
   `analyze_trajectories` (≥ 0.8.0), writing `results/motco_analysis/<cohort>/`
   (MOTCO's `result.json`, `latent_scores.csv`, `pls_loadings.csv`,
   `stage_means.csv`, plus `analysis_metadata.csv`, `analysis_manifest.json`,
   `summary.md`). Never assemble it from `motco plsr` + `motco de`.
7. `rank-sensitivity` – from `results/motco_input/` and
   `results/motco_analysis/` only, per cohort: MOTCO's dimensionality
   sensitivity analysis, the same `analyze_trajectories` call with the latent
   rank forced to each value of `ANALYSIS_FORCED_RANKS`, everything else
   (inputs, covariates, orientation, permutations, seed) equal to `analyze`.
   The CV-selected analysis is read and checked (export hashes, MOTCO
   version, seed, permutations), never rerun or modified, and is the
   reference row. Writes `results/motco_rank_sensitivity/<cohort>/`
   (`rank_<kk>/` analysis directories, `rank_sensitivity.csv`,
   `rank_sensitivity.md`, `rank_sensitivity_manifest.json`).
8. `preprocessing-sensitivity` – from `results/pseudobulk/` (incl.
   `donor_technical.csv` + `technical_manifest.json` for the gene-set,
   residualized, trend-normalized and candidate variants), the annotation,
   `results/motco_input/` and `results/motco_analysis/`, per variant of the
   registry `PREPROCESSING_VARIANTS` (`keep_chrx`, `exclude_families`,
   `rna_trend_ratio`, `atac_mean_accessibility`, `residualized`,
   `revised_candidate`, `trend_normalized`, `revised_candidate_trend`) and
   per cohort: the export rerun with the variant's registered setting(s)
   changed (`_export_cohort`, shared with `export`: `excluded_chromosomes`,
   `excluded_gene_sets`, `selection_rule` per omic, `normalization` per omic,
   `residualize`) and the analysis rerun unchanged on the variant blocks
   (`_analyze_cohort` with the `analyze` parameters), compared with the
   primary analysis. The primary export and analysis are read and checked
   (pseudobulk, covariates and annotation hashes, thresholds, `k`, cells;
   export hashes, MOTCO version, seed, permutations, CV selection; the
   technical table against the pseudobulk manifest, `TECHNICAL_GENE_SETS`
   and completeness for the cohort donors), never rerun or modified. For
   `revised_candidate`, `trend_normalized` and `revised_candidate_trend` the
   predeclared `ADOPTION_RULE` is evaluated from the generated files into
   `adoption_check.csv` (no p-value read), and `ADOPTION_PRECEDENCE` between
   the last two is written into their comparison manifests. Writes
   `results/motco_sensitivity/<variant>/` (`motco_input/<cohort>/` plus
   `residualization.csv` when residualizing and `normalization_trend.csv`
   when trend-normalizing, `motco_analysis/<cohort>/`,
   mirroring the primary layout so `rank-sensitivity` and `axis-diagnostics`
   can read a variant, and `comparison/<cohort>/` with `sensitivity.csv`,
   `sensitivity.md`, `sensitivity_manifest.json` and, for the candidates,
   `adoption_check.csv`). `--variant all` skips complete variant × cohort
   pairs without `--force` (resumable).
9. `axis-diagnostics` – from `results/motco_input/`, `results/motco_analysis/`,
   `results/pseudobulk/donor_technical.csv` and, when present,
   `results/motco_rank_sensitivity/<cohort>/rank_*/` only (any workstation),
   per cohort: Spearman and covariate-partial Spearman correlations of every
   latent component of every analysis read with the predeclared variables
   (`DIAGNOSTIC_VARIABLES`: the technical table's numeric columns, the
   export's TMM factors and log10 library sizes, the six analysis
   covariates, stage, sex), the same for the first three PCs of each
   z-scored block, the geometry shares along each component from
   `stage_means.csv`, and the feature-selection accounting from
   `features_<omic>.csv`. Every analysis manifest's input hashes are checked
   against the export files first; nothing under the inputs is written.
   Writes `results/motco_diagnostics/<cohort>/` (`axis_correlations.csv`,
   `block_pcs.csv`, `geometry_shares.csv`, `feature_selection.csv`,
   `diagnostics.md` with flagged rows first, `diagnostics_manifest.json`).
   Its outputs feed no export or analysis command.

`README.md` documents CLI flags and the Python API in detail; keep it in sync
when behavior or defaults change.

## Layout

- `src/sead_sexdiff/cli.py` – argparse CLI (`build_parser`, `main`); thin layer
  that calls into `data.py` and prints `[command]`-prefixed progress messages.
- `src/sead_sexdiff/data.py` – all logic: S3 download, subsetting/QC,
  pseudobulk, differential expression, covariates, technical covariates,
  MOTCO export, analysis, the rank / preprocessing sensitivity analyses and
  the axis diagnostics.
- `src/sead_sexdiff/DEG.R` – exploratory R script (MAST mixed model); not wired
  into the CLI.
- `scripts/` – `pseudobulk.sbatch` and `technical.sbatch` (cluster jobs),
  `sync_pseudobulk.sh` (rsync of `results/pseudobulk/` to a workstation,
  per-donor caches excluded), `check_tmm_edger.R`
  (one-off edgeR cross-check of the export's TMM / log2-CPM) and
  `check_cyclic_loess_limma.R` (limma cross-check of the variants' cyclic
  loess on `_export_cohort` debug files).
- `data/`, `results/` – gitignored; never commit their contents.
- `reports/` – committed case-study reports (e.g. `motco-astrocyte-<date>.md`).
  Their numbers are pasted from generated outputs (`summary.md`, manifests),
  never retyped, and checked against `result.json` before committing.
- `openspec/` – OpenSpec change proposals/specs (see below).

## Environment & commands

- Python `>=3.12,<3.14`, managed with **uv** (`uv.lock` is committed).
- Setup: the Nix flake (`flake.nix`) provides Python 3.12, uv, R, and the R
  packages used by `DEG.R`, `check_tmm_edger.R` and `check_cyclic_loess_limma.R`; its shell hook runs `uv sync --locked`. Enter it
  with direnv (`.envrc` → `use flake`) or `nix develop`. Python packages come
  from uv, not Nix; R packages come from Nix (add them to `flake.nix`).
- Run the CLI: `uv run sead-sexdiff <download|pseudobulk|covariates|technical-covariates|export|analyze|rank-sensitivity|preprocessing-sensitivity|axis-diagnostics|subset|de> [options]`.
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
- The nucleus selection of the MTG steps is one code path,
  `select_astrocyte_nuclei` (taxonomy mapping, `Reference` excluded,
  `RNA_MIN_GENES` on the `UMIs` layer; ATAC `read_counts=False` reads `obs`
  only): `pseudobulk_donor` sums what it returns and `technical_donor`
  summarizes the same nuclei. Changing it changes both steps; the pseudobulk
  outputs must then be re-checked against the per-donor cache.
- Technical-covariate constants (flag and ask before changing any of them):
  `TECHNICAL_GENE_SETS` (ribosomal `^RP[LS]\d+[A-Z]?$` / `^RPLP[0-2]$` /
  `RPSA`, so pseudogenes and kinases are excluded; the 13 `MT-` genes; `HBB`,
  `HBA1`, `HBA2`; recorded as matched lists, an empty match is an error),
  `TECHNICAL_OBS_COLUMNS` per omic, `TECHNICAL_LIBRARY_COLUMN` (`load_name`),
  `TECHNICAL_SD_DDOF` (1) and `TECHNICAL_MIN_SUPERTYPE_NUCLEI` (2, the η²
  definition), nucleus-weighted means beside UMI-weighted pseudobulk
  fractions, the within-donor supertype share as the biology-versus-
  preparation statistic (read beside the manifest's between-donor share,
  never alone), the `COVARIATE_RTOL` cross-check against `rna_counts.h5ad`,
  and the rule that the technical files stay apart from
  `donor_covariates.csv` / `covariates_manifest.json`.
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
- Rank-sensitivity constants (flag and ask before changing any of them):
  `ANALYSIS_FORCED_RANKS` (2, 3, 4, 6, 9, 12: 2 = `n_stages − 1`, 3 the next
  rank, 4/6/9/12 the steps of MOTCO's latent-rank ladder with 4 = 2(K − 1)
  the forced rank with a Type I check; a constant, never a flag, never
  chosen from results), the same `ANALYSIS_SEED` and `ANALYSIS_PERMUTATIONS`
  at every rank (so every rank draws the same permutation indices), rank 1
  not forceable (MOTCO accepts `forced_components ≥ 2`; the CV-selected
  analysis is the reference row), `n_jobs = 1` for reported results, and the
  rule that the CV-selected analysis stays primary.
- Reporting MOTCO results follows the binding rules in
  `../MOTCO/docs/roadmap.md` (Phase 6, Next changes item 1): CV-selected rank
  beside `angle`; `angle`/`shape` non-rejections are uninformative; a
  `delta` rejection is a difference, not specifically of size; report the
  covariate block and the feature filtering beside the results. The
  forced-rank sweep is read as sensitivity only: the CV-selected analysis
  stays primary, a rejection status that differs between ranks is rank
  dependence (not a result at that rank), the non-rejection rules hold at
  every rank, MOTCO's forced-rank operating characteristics (Type I at
  forced 4 only, unadjusted; `shape` power loss and `angle` leakage at
  forced 2(K − 1); wider `angle` null above `n_stages − 1`) are carried as
  limitations, and the adjusted analysis has no Type I check at any forced
  rank.
- Preprocessing-sensitivity constants (flag and ask before changing any of
  them): `PREPROCESSING_VARIANTS`, a registry keyed by variant name with a
  `description` and the `_export_cohort` keyword(s) it changes (`keep_chrx`:
  `excluded_chromosomes = ("chrY", "chrM")`, so chrX is kept and the whole
  export pipeline is recomputed with it present; `exclude_families`:
  `excluded_gene_sets = ("ribosomal", "mitochondrial", "hemoglobin")` of
  `TECHNICAL_GENE_SETS`, removed at step 2 after the chromosome removal with
  status `excluded_gene_set`, RNA only; `rna_trend_ratio`: `selection_rule =
  {rna: trend_ratio, atac: variance}`; `atac_mean_accessibility`:
  `selection_rule = {rna: variance, atac: mean}`; `residualized`:
  `residualize = True`, the technical table substituted at run time;
  `revised_candidate`: the four together, with `adoption_check: True`;
  `trend_normalized`: `normalization = cyclic_loess` for both omics, with
  `adoption_check: True`; `revised_candidate_trend`: `revised_candidate`'s
  four plus the cyclic loess, with `adoption_check: True`; never a flag,
  never chosen from results, entries only added by a change proposal),
  `SELECTION_RULES` and `TREND_N_BINS` (20 equal-count bins of
  mean log2-CPM via `np.array_split`, median variance at the median mean,
  `np.interp` flat beyond the end bins, ratio = variance / trend, ties by
  name; dropped status `selection` under a non-variance rule, `variance`
  under the primary's), the residualization position (after log2-CPM, before
  the selection; exported value = mean + residual; `np.linalg.lstsq` once
  for all features; R² per feature) and set (`RESIDUALIZATION_COVARIATES` =
  the export's own `log10_lib_<omic>` after step 2 plus
  `RESIDUALIZATION_CANDIDATE_FRACTIONS` minus
  `RESIDUALIZATION_CELL_STATE_EXCLUSIONS`, now `hemoglobin_pseudobulk` with
  the technical-axis report's section 2 as the reason; sex, stage and the
  analysis covariates never enter; the design written as `residualization.csv`,
  `RESIDUALIZATION_FILE`, never an analysis input), the trend normalization
  (`NORMALIZATION_RULES` = `none` (the primary's) / `cyclic_loess`: limma's
  fast `normalizeCyclicLoess` reimplemented with statsmodels' `lowess`;
  position after log2-CPM and before the residualization and the selection;
  fitting set all abundance-passing features of the block; reference the
  per-feature mean over donors, recomputed at each of `LOESS_ITERATIONS` (3);
  `LOESS_ROBUSTNESS_ITERATIONS` (3, `loessFit`'s default), `delta` =
  `LOESS_DELTA_FRACTION` (0.01) of the range of A; span by `LOESS_SPAN_RULE`,
  limma's adaptive `chooseLowessSpan` defaults, per block; label-free; the
  per-donor tilt accounting over `NORMALIZATION_TILT_BINS` (20) bins written
  as `NORMALIZATION_TREND_FILE`, never an analysis input; the primary export
  byte-identical with `normalization = none`; changes to the numpy cyclic
  loess must be re-checked with `scripts/check_cyclic_loess_limma.R`,
  tolerance 1e-6), the rule that the variant
  analysis keeps the six analysis covariates (double adjustment of the
  library-size proxies accepted and stated), `ADOPTION_CHECK_VARIABLES` (the
  primary export's TMM factors, four technical QC means, RIN, PMI: technical
  variables that are no regressor) and `ADOPTION_RULE` (holds when no
  candidate component reaches `|Spearman| ≥ DIAGNOSTIC_CORRELATION_THRESHOLD`
  against any check variable in both cohorts; p-values, rejections and CV
  AUROC never read; a missing variable is an error, not a skipped cell;
  written to `PREPROCESSING_ADOPTION_TABLE` for the entries with
  `adoption_check` only), `ADOPTION_PRECEDENCE` with
  `ADOPTION_PRECEDENCE_CANDIDATES` (`trend_normalized`,
  `revised_candidate_trend`: if `trend_normalized` meets `ADOPTION_RULE` in
  both cohorts it is the candidate whatever the other gives, because it
  changes the least; if only `revised_candidate_trend` does, that one; if
  neither, the variants stay sensitivity; only the `holds` values are read;
  written as `adoption_precedence` into both candidates' comparison
  manifests, `null` with reason `incomplete` until both have run in both
  cohorts), the rule
  that a variant reruns the whole export at the primary's thresholds and N
  and the analysis at the primary's settings (same `ANALYSIS_SEED` and
  `ANALYSIS_PERMUTATIONS`, CV rank selection), the technical reference check
  (`_check_technical_reference`: pseudobulk manifest hash, gene sets,
  completeness; only for variants that read the table), the output names
  (`PREPROCESSING_*`), `n_jobs = 1` for reported results, and the rule that
  the primary analysis stays primary. Every variant is read as sensitivity
  only: a rejection status that differs between the variant and the primary
  is preprocessing dependence (not a result of the variant), a CV-selected
  rank that differs is reported beside `angle` (not better or worse), the
  non-rejection rules hold in every variant, MOTCO's operating
  characteristics were measured on none of the variant blocks (chrX-kept,
  family-excluded, trend- or mean-selected, residualized, trend-normalized),
  so no variant carries a Type I or power claim of its own, the residualized
  variants carry the stated limitations (degrees of freedom spent before the
  permutation null, stage and sex signal carried by the technical set
  removed, library-size proxies adjusted twice), and the trend-normalized
  variants carry theirs (the cyclic loess assumes most features do not
  change at any abundance level, so a stage or sex effect shaped like a tilt
  is removed with the artefact; the curve is least determined at the rare
  end, where the selection draws most features; the lncRNA-versus-coding
  offset at matched abundance is not corrected; no Type I check). If the
  adoption check holds in both cohorts a later change may propose the
  candidate (between the two trend candidates, the one `ADOPTION_PRECEDENCE`
  names) as a revised primary with the current primary kept as the
  reference; otherwise the variants stay sensitivity analyses. The first latent axis's technical
  correlations per variant (`axis-diagnostics` on a variant root) may be
  attributed to each single setting as a measurement-space finding only.
- Axis-diagnostic constants (flag and ask before changing any of them):
  `DIAGNOSTIC_VARIABLES` (the variable sources; resolved at run time and
  recorded with sources, never chosen or extended from results),
  `DIAGNOSTIC_CORRELATION_THRESHOLD` (0.3, a reporting flag, never an
  adoption rule: adoption rules belong to the change that runs a variant and
  are written in its proposal before its runs), `DIAGNOSTIC_N_BLOCK_PCS` (3),
  `DIAGNOSTIC_QUANTILES`, `DIAGNOSTIC_MEAN_CUTS`, Spearman (not Pearson; the
  technical variables are skewed) with the partial correlation defined as
  the Spearman correlation of OLS residuals on an intercept and the six
  analysis covariates, population-SD z-scoring of the blocks, no p-values.
  The diagnostic reads sex and stage, so it is a diagnostic and never a
  filter: its outputs are an input of no `export`, `analyze`,
  `rank-sensitivity` or `preprocessing-sensitivity` run, and preprocessing
  variants are predeclared in a change proposal, never derived from a
  diagnostic value. Reading rules for the report: an axis correlated with a
  technical variable is a measurement-space finding about the primary
  analysis, not a reinterpretation of its statistics; the primary analysis
  stays primary with its non-rejection rules; a gene-set fraction is a
  residualization candidate only when its between-donor share dominates and
  its within-donor supertype share is small (cell state otherwise).
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
  node. `technical-covariates --omic both` needs the raw objects and
  `results/pseudobulk/` (incl. the per-donor cache): `scripts/technical.sbatch`
  after the pseudobulk job, or nohup on the login node; resubmitting resumes.
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

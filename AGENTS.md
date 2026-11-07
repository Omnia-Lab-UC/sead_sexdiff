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
   `analyze_trajectories` (≥ 0.18.0) with CV rank selection at the rank
   floor `ANALYSIS_RANK_FLOOR` (2), the RRPP null drawn under the residual
   correction `ANALYSIS_RESIDUAL_CORRECTION` (`df`) and the split interaction
   test `ANALYSIS_SPAN_SPLIT_TEST`, writing `results/motco_analysis/<cohort>/`
   (MOTCO's `result.json`, `latent_scores.csv`, `pls_loadings.csv`,
   `stage_means.csv`, plus `analysis_metadata.csv`, `analysis_manifest.json`,
   `summary.md`). Never assemble it from `motco plsr` + `motco de`.
7. `rank-sensitivity` – from `results/motco_input/` and
   `results/motco_analysis/` only, per cohort: MOTCO's dimensionality
   sensitivity analysis, the same `analyze_trajectories` call with the latent
   rank forced to each value of `ANALYSIS_FORCED_RANKS`, everything else
   (inputs, covariates, orientation, correction, split test, permutations,
   seed) equal to `analyze`.
   The CV-selected analysis is read and checked (export hashes, MOTCO
   version, seed, permutations, rank floor), never rerun or modified, and is
   the reference row. Writes `results/motco_rank_sensitivity/<cohort>/`
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
10. `latent-sensitivity` – from `<root>/motco_input/`, `<root>/motco_analysis/`
   and `<root>/pseudobulk/` (`manifest.json`, `donor_technical.csv`,
   `technical_manifest.json`) only, per cell type (`--cell-type`, the root
   from `cell_type_dirs`; the registry, the covariate set, the conditioning
   and every setting the same for every cell type; every manifest read
   checked by `check_manifest_cell_type`, the arm and comparison manifests
   recording the cell type) and per cohort: the two arms of the registry
   `LATENT_SENSITIVITY_ARMS` on the unchanged primary export, both with
   `LATENT_SENSITIVITY_COVARIATES` (the six analysis covariates plus the raw
   `ribosomal_pseudobulk` and `mitochondrial_pseudobulk` fractions):
   `extended_covariates` with the covariate-blind basis, `latent_adjusted`
   with MOTCO's `covariate_latent_adjustment = "stage"` (MOTCO ≥ 0.9.0,
   `LATENT_MOTCO_MIN_VERSION`); everything else equals `analyze`
   (`_analyze_cohort` with `covariates=`, the fractions merged by
   `load_analysis_inputs(extra_covariates=...)`). The primary analysis is the
   reference row, read and checked (`_check_reference`, plus orientation and
   CV settings), never rerun or modified; the technical table is checked by
   `_check_technical_reference`. The predeclared mediator check
   (`mediator_check`: Spearman of each added fraction with stage, all donors
   and per sex, flagged at `DIAGNOSTIC_CORRELATION_THRESHOLD`, no p-value,
   never a gate) is written before the arms. Writes
   `<root>/motco_latent_sensitivity/` (`<arm>/motco_analysis/<cohort>/`,
   mirroring the primary layout with the five export files as
   `inputs.files` so `axis-diagnostics --analysis-dir` can read an arm, and
   `comparison/<cohort>/` with `latent_sensitivity.csv`,
   `latent_sensitivity.md`, `mediator_check.csv`,
   `latent_sensitivity_manifest.json`).
11. `stage-de` – from `results/pseudobulk/` (`rna_counts.h5ad`,
   `atac_counts.h5ad`, `donor_table.csv`, `manifest.json`),
   `results/motco_input/` and `results/motco_analysis/` only, per cohort and
   omic: per-feature DESeq2 (pydeseq2) of the four models of the registry
   `STAGE_DE_MODELS` (`factor_additive`, `factor_interaction`,
   `trend_additive`, `trend_interaction`) with sex, stage and
   `STAGE_DE_COVARIATES` on the raw counts of the cohort donors, on the
   export's abundance-passing features (`features_<omic>.csv` status
   `variance` or `kept`; `exported` marks `kept`), one Wald test per tested
   coefficient (a unit contrast on the design column found by name). The
   counts are bound to the export by the pseudobulk manifest hash and a
   reproduction (`preprocess_block` rerun on the counts: statuses, library
   sizes, TMM factors, tested-feature log2-CPM variances, exported values),
   since the export records no count-file hash; the analysis by
   `_check_reference`. The predeclared
   `STAGE_DE_CONCORDANCE_RULE` compares the exported interaction count of
   `factor_interaction` with the primary `delta` rejection status. Writes
   `results/stage_de/<cohort>/` (`<omic>/<model>__<coefficient>.csv`,
   `<omic>/<model>__fit.json` (the resume record), `stage_de_summary.csv`,
   `stage_de.md`, `stage_de_manifest.json`). It reads sex and stage, so its
   outputs feed no export, analysis or sensitivity command.
12. `stage-de-limma` – from the same inputs as `stage-de` plus its outputs
   (`results/stage_de/<cohort>/`), per cohort and omic: the four models of
   `STAGE_DE_MODELS` fitted a second time with limma-voom through the
   packaged R script `stage_de_limma.R` (`Rscript`, one BLAS thread, exchange
   files; never rpy2): edgeR TMM on the tested counts with the export's
   trims, checked against the export manifest's factors to
   `STAGE_DE_LIMMA_TMM_ATOL` → `voomWithQualityWeights` → `lmFit` → `eBayes`
   (`STAGE_DE_LIMMA_METHOD`, limma's defaults, checked against the installed
   ones). A moderated t per tested coefficient, a joint moderated F per factor
   model (`STAGE_DE_LIMMA_JOINT_TESTS`: `stage_joint`, `interaction_joint`)
   and the union, BH over all tested features. The same donors, counts,
   designs (the formulaic matrix, rebuilt by R's `model.matrix` as a check)
   and features as `stage-de`, by reusing `_stage_de_prepare` unchanged; the
   DESeq2 outputs are checked (complete, made from the inputs now on disk,
   tables hashing as recorded, summary equal to a recomputation) and never
   written. The predeclared `STAGE_DE_LIMMA_CONCORDANCE_RULE`
   (`interaction_joint`) sits beside the DESeq2 cell with `methods_agree`.
   Writes `results/stage_de_limma/<cohort>/` (`<omic>/<model>__<label>.csv`,
   `<omic>/<model>__fit.json` (the resume record), `stage_de_limma_summary.csv`,
   `method_comparison.csv`, `stage_de_limma.md`, `stage_de_limma_manifest.json`).
   It reads sex and stage, so its outputs feed no export, analysis or
   sensitivity command.
13. `positive-control` – from the generated files of both cell types only
   (control `results/microglia/`, reference `results/`), per cohort: the
   predeclared `POSITIVE_CONTROL_EXPECTATIONS` E1 (exported RNA stage counts
   at FDR 0.05, control > reference within each stage-DE method: DESeq2
   `factor_additive` union, limma `stage_joint`; holds only under both, under
   one only `method_dependent`, not holding), E2 (primary `cv_mean_auroc`,
   control > reference) and E3 (control's primary `pls_0`: |Spearman| with
   stage > the largest with a technical variable of
   `POSITIVE_CONTROL_TECHNICAL_SOURCES`, from `axis_correlations.csv`); the
   control holds when all three hold in `POSITIVE_CONTROL_PRIMARY_COHORT`.
   Reads counts, the CV AUROC and raw Spearman correlations only (no p-value,
   rejection status or trajectory statistic; nothing tested). Every manifest
   read must record its root's cell type (`POSITIVE_CONTROL_CELL_TYPES`) and
   the chain must hold (both stage-DE steps and the diagnostic made from the
   analysis read for E2, the limma step from the DESeq2 summary read for E1).
   `--basis` (`POSITIVE_CONTROL_BASES`: `primary`, the default; or
   `latent_adjusted`, E2 and E3 read from each cell type's `latent_adjusted`
   arm of `latent-sensitivity` and the diagnostic made on it, E1 unchanged,
   the arm chain checked: made from the primary export on disk, compared with
   the primary analysis on disk, the registered conditioning, the diagnostic
   made from the arm). Writes `results/microglia/positive_control/<cohort>/`
   (`positive_control.csv`, `positive_control.md`,
   `positive_control_manifest.json`), or `positive_control/<basis>/<cohort>/`
   on a non-primary basis.
14. `donor-flags` – from the raw donor objects' `obs` only (cluster login
   node, under a second per object, modelled on `covariates`): SEA-AD's
   `Severely Affected Donor` coded as `severely_affected` (`DONOR_FLAG_VALUES`:
   `Y` → true, `N` → false, anything else an error naming the donor) and,
   recorded as read and used by nothing, `Cognitive Status`, `Braak`,
   `CERAD score`, `Thal`, `LATE` (`DONOR_FLAG_COLUMNS`) and the `Continuous
   Pseudo-progression Score` where the objects hold it
   (`DONOR_FLAG_ATAC_ONLY_COLUMNS`: the ATAC release only; its absence from
   the RNA objects is a finding). RNA values are used; with `--omic both`
   every column is cross-checked against ATAC (a flag or Reference-status
   disagreement is an error, any other a finding with both values). Reference
   donors (ADNC `Reference`, as the donor table) are left out and listed.
   Writes `<root>/pseudobulk/donor_flags.csv` and `donor_flags_manifest.json`
   (`DONOR_FLAG_OUTPUTS`) and nothing else, never `donor_covariates.csv` or
   `covariates_manifest.json`; a flagged count other than
   `DONOR_FLAG_EXPECTED_FLAGGED` (11, SEA-AD's) is a finding, not corrected.
   The table is an input of `exclusion-sensitivity` and of nothing else.
15. `exclusion-sensitivity` – from `<root>/pseudobulk/` (the export's inputs,
   `donor_technical.csv`, `donor_flags.csv` and their manifests), the
   annotation, `<root>/motco_input/` and `<root>/motco_analysis/` only, per
   cell type (`--cell-type`), per variant of the registry
   `DONOR_EXCLUSION_VARIANTS` (`excluding_severely_affected`: the cohort
   donors whose `severely_affected` is true) and per cohort: the flag check
   (`flag_check.csv`, written and put in place before the variant export:
   the flagged donors per sex × stage cell; for the primary `pls_0`, stage,
   the technical table's numeric columns and the primary export's TMM
   factors and log10 library sizes, the rank-biserial AUC of flagged against
   unflagged donors over the cohort and within stage 2, and among the
   unflagged donors the Spearman with stage and with `pls_0`; no p-value,
   never a gate); the export rerun on the cohort donors minus the excluded
   ones (`_export_cohort` with the primary settings and `donor_exclusion=`,
   the manifest's `donor_exclusion` block: variant, flag, flag-table and
   manifest sha256, excluded donors, cells and `k` before and after; cohort
   thresholds unchanged, `k` recomputed by the export's rule; `metadata.csv`
   and `covariates.csv` checked equal to the primary's minus the excluded
   rows); the analysis rerun unchanged (`_analyze_cohort` with the `analyze`
   parameters); the predeclared axis check (`axis_check.csv`,
   `DONOR_EXCLUSION_AXIS_RULE`: along the variant's `pls_0`, |Spearman| with
   stage strictly above the largest with a technical variable of
   `POSITIVE_CONTROL_TECHNICAL_SOURCES`, from the variant export; the control
   cell type's value is the expectation, the reference's descriptive). The
   primary export and analysis are read and checked (`_check_export_reference`,
   `_check_reference`, orientation and CV settings, the six covariates), the
   technical table by `_check_technical_reference`, the flag table for its
   cell type and a value for every cohort donor; an exclusion that would
   empty a sex × stage cell stops the step; nothing is rerun or modified.
   Writes `<root>/motco_exclusion_sensitivity/<variant>/` (`flag_check/<cohort>/`,
   `motco_input/<cohort>/` and `motco_analysis/<cohort>/` mirroring the
   primary layout so `axis-diagnostics`, `stage-de` and `stage-de-limma` read
   the variant through their directory flags, and `comparison/<cohort>/` with
   `exclusion_sensitivity.csv`, `axis_check.csv`, `exclusion_sensitivity.md`,
   `exclusion_sensitivity_manifest.json`, written last); each directory as
   `<cohort>.tmp` renamed when complete, the cohort's directories removed on
   an error; `--variant all` skips complete variant × cohort pairs without
   `--force` (resumable). On a variant, `stage-de` and `stage-de-limma`
   (`_stage_de_prepare`) remove the export manifest's recorded
   `donor_exclusion` donors from the donor table's cohort before the
   donor-agreement check (they must be cohort donors and the rest must equal
   the export's donors) and record the block in their manifests; without the
   block nothing changes.
16. `figures` – from the generated outputs of both cell types only (every
   cell type's `cell_type_dirs` plus `results/motco_sensitivity/`; never raw
   objects or `results/superseded/`), in `src/sead_sexdiff/figures.py`: the
   registry `FIGURES` (eleven figures, `01_cohorts` … `11_span_split`), every
   family of analyses enumerated from the registries (`ANALYSIS_FORCED_RANKS`,
   `PREPROCESSING_VARIANTS`, `LATENT_SENSITIVITY_ARMS`,
   `DONOR_EXCLUSION_VARIANTS`, `POSITIVE_CONTROL_BASES`, `EXPORT_COHORTS`),
   a missing registered output an error naming the step that writes it.
   Before anything is drawn, a consistency pre-pass (`check_consistency`)
   checks the manifests the registries name and every other `*manifest.json`
   under the registered directories:
   - the cell type of each, against its root's;
   - every recorded `{"path", "sha256"}` entry existing under `results/`
     outside `superseded/` and hashing as recorded;
   - one MOTCO version (the astrocyte primary's) and the registered
     measurement in every analysis drawn (`ANALYSIS_RESIDUAL_CORRECTION`,
     `ANALYSIS_RANK_FLOOR`, the forced rank under `rank-sensitivity` only,
     `ANALYSIS_PERMUTATIONS`, `ANALYSIS_SEED`, the split test
     `ANALYSIS_SPAN_SPLIT_TEST` with its `span_split` block, except the arms
     of `LATENT_SENSITIVITY_SPLIT_UNAVAILABLE`, which carry neither).

   Per figure:
   - the reproductions (adjusted cell means against `stage_means.csv`,
     Spearman against `axis_correlations.csv`, counts against the stage-DE
     summaries and `method_comparison.csv`, cells against the pseudobulk
     manifest);
   - the chain: each diagnostic, stage-DE and rank manifest read records the
     analysis drawn, and figure 10 checks limma's `inputs.deseq2`.

   Every read goes through one recording reader (`FigureReader`), so a
   figure's inputs are the reader's log. Writes `results/figures/`, a
   gitignored working copy:
   - PNG, PDF and SVG per figure, with deterministic bytes and a fingerprint
     of inputs, module and matplotlib in their metadata;
   - `figures.md`, with generated captions in three parts: construction,
     values read at run time, `FIGURE_READING_RULES`;
   - `figures_manifest.json`.

   The working copy is redrawn whole through `<out>.tmp` when any figure is
   stale and left untouched when all are fresh. `--check` lists the stale
   figures (exit 0 / 1 / 2, nothing written). `--report <stem>` runs the
   checks, redraws if needed and freezes the SVGs, `figures.md` and the
   manifest under `reports/figures/<stem>/`, never overwritten. A display
   step: it tests nothing, and its outputs feed no export, analysis,
   sensitivity, diagnostic, stage-DE or positive-control step.

Cell types: steps 2–12 except the rank and preprocessing sensitivity steps,
and steps 14–15 (the eleven commands `pseudobulk`, `covariates`,
`technical-covariates`, `export`, `analyze`, `axis-diagnostics`,
`latent-sensitivity`, `stage-de`, `stage-de-limma`, `donor-flags`,
`exclusion-sensitivity`) take `--cell-type`, a name
of the registry `CELL_TYPES` (default `astrocyte`, results root `results/`;
`microglia`, results root `results/microglia/`). It sets the population (one
selection path) and every default directory (`cell_type_dirs`: `<root>/pseudobulk/`,
`<root>/motco_input/`, …); explicit directory flags win. Every manifest of
these steps records `cell_type` and every reader refuses another cell type's
manifest (`check_manifest_cell_type`, naming the file and both cell types,
nothing written); a manifest written before the registry has no key and
reads as `astrocyte` with a finding. The per-donor caches record it too; the
cell type stays out of the stage-DE fit-record fingerprints. The sensitivity
steps `rank-sensitivity` and `preprocessing-sensitivity` have no flag and run
on astrocytes only; `positive-control` reads both cell types by role and
takes `--basis` instead; `figures` reads both cell types by role and takes
neither (its reader checks every manifest's cell type against its root's).

`README.md` documents CLI flags and the Python API in detail; keep it in sync
when behavior or defaults change.

## Layout

- `src/sead_sexdiff/cli.py` – argparse CLI (`build_parser`, `main`); thin layer
  that calls into `data.py` (and `figures.py` for `figures`) and prints
  `[command]`-prefixed progress messages.
- `src/sead_sexdiff/data.py` – all logic: S3 download, subsetting/QC,
  pseudobulk, differential expression, covariates, technical covariates,
  MOTCO export, analysis, the rank / preprocessing / latent-adjustment
  sensitivity analyses, the axis diagnostics, the pseudobulk stage DE
  (DESeq2, and limma-voom through the R script below), the cell-type
  registry, the positive control, the donor flags and the donor-exclusion
  sensitivity analysis. Plotting is the one exception (`figures.py`).
- `src/sead_sexdiff/figures.py` – the `figures` command (step 16): the
  registry `FIGURES`, the display constants, the recording reader, the
  consistency pre-pass, the eleven figure functions with their caption
  builders, the manifest, the staleness check and the report freeze. It
  imports constants and helpers from `data.py` (`_git_revision`,
  `_sha256_file`, `check_manifest_cell_type`, `cell_type_dirs`, `spearman`,
  the registries) and nothing imports it but `cli.py`. It is the exception to
  "all logic in `data.py`" for plotting (its own dependency and rcParams;
  the user's decision of 2026-11-02). Its sha256 enters every figure's
  fingerprint, so editing it makes every figure stale.
- `src/sead_sexdiff/stage_de_limma.R` – the R side of `stage-de-limma`
  (TMM, voomWithQualityWeights, lmFit, eBayes, the R-side checks and the
  exact-F diagnostic); packaged beside `data.py` and run with `Rscript` on
  exchange files (`--check`, `--versions` modes; rerunnable by hand on
  `run_stage_de_limma(debug_dir=...)` files). Its sha256 enters every fit
  record, so editing it makes the existing records stale.
- `src/sead_sexdiff/DEG.R` – exploratory R script (MAST mixed model); not wired
  into the CLI.
- `scripts/` – `pseudobulk.sbatch` and `technical.sbatch` (cluster jobs;
  `CELL_TYPE`, default `astrocyte`, passed as `--cell-type`),
  `sync_pseudobulk.sh` (rsync of a cell type's pseudobulk directory to a
  workstation, `--cell-type`, per-donor caches excluded), `check_tmm_edger.R`
  (one-off edgeR cross-check of the export's TMM / log2-CPM),
  `check_cyclic_loess_limma.R` (limma cross-check of the variants' cyclic
  loess on `_export_cohort` debug files), `check_stage_de_deseq2.R` (R
  DESeq2 cross-check of one `stage-de` fit on `run_stage_de(debug_dir=...)`
  files) and `move_superseded.py` (moves superseded output directories under
  `results/superseded/<change>/`, hash-checked, standard library only).
- `data/`, `results/` – gitignored; never commit their contents.
- `results/figures/` – the working copy of `figures`, redrawn freely, never
  moved to `results/superseded/`. The record of a figure is the frozen set a
  report carries.
- `results/superseded/<change>/` – superseded outputs, moved (never
  deleted) by the change that re-issued them, with their path relative to
  `results/` kept and `MOVED.json` (old / new path, host, time, sha256 per
  file); never an input of any step or a default directory, not cited as
  results. Now `mtg-span-split-test/` (the MOTCO 0.13.0 `df` outputs
  without the split interaction test, both cell types: astrocyte
  `motco_analysis`, `motco_rank_sensitivity`, `motco_diagnostics`,
  `motco_sensitivity`, `motco_latent_sensitivity`, `stage_de`,
  `stage_de_limma`, `motco_exclusion_sensitivity`; microglia
  `microglia/motco_analysis`, `motco_diagnostics`,
  `motco_latent_sensitivity`, `stage_de`, `stage_de_limma`,
  `motco_exclusion_sensitivity`, `positive_control`),
  `mtg-residual-correction-primary/` (the floored, uncorrected
  outputs of MOTCO 0.10.0, both cell types: astrocyte `motco_analysis`,
  `motco_rank_sensitivity`, `motco_diagnostics`, `motco_sensitivity`,
  `motco_latent_sensitivity`, `stage_de`, `stage_de_limma`,
  `motco_exclusion_sensitivity`; microglia `microglia/motco_analysis`,
  `motco_diagnostics`, `motco_latent_sensitivity`, `stage_de`,
  `stage_de_limma`, `motco_exclusion_sensitivity`, `positive_control`),
  `mtg-rank-floor-primary/` (the rank-1 outputs of MOTCO 0.8.0 / 0.9.0, both
  cell types) and
  `mtg-astrocyte-latent-adjustment-sensitivity/motco_analysis/` (the former
  `results/motco_analysis_0.8.0/`).
- `reports/` – committed case-study reports (e.g. `motco-astrocyte-<date>.md`).
  Their numbers are pasted from generated outputs (`summary.md`, manifests),
  never retyped, and checked against `result.json` before committing. Every
  report committed from `mtg-case-study-figures` on carries the figure set
  frozen by `figures --report <its stem>` after its last output run
  (`reports/figures/<stem>/`: the eleven SVGs, `figures.md`,
  `figures_manifest.json`; committed, never regenerated or edited). It embeds
  `figures/<stem>/<nn>_<name>.svg` and pastes each caption verbatim from the
  frozen `figures.md`. Before the commit, its check script verifies:
  - the pasted captions against `figures.md`;
  - the frozen SVGs against the frozen manifest;
  - every input the frozen manifest records, and every hash recorded by the
    manifests it read, against the files on disk.

  If an output is rerun after freezing, the uncommitted frozen set is
  removed and `--report` is run again. Every report freezes the full set:
  unchanged figures are byte-identical, so git stores them once. Reports
  committed before that change are not edited.
- `openspec/` – OpenSpec change proposals/specs (see below).

## Environment & commands

- Python `>=3.12,<3.14`, managed with **uv** (`uv.lock` is committed).
- Setup: the Nix flake (`flake.nix`) provides Python 3.12, uv, R, and the R
  packages used by `DEG.R`, `check_tmm_edger.R`, `check_cyclic_loess_limma.R`, `check_stage_de_deseq2.R` (DESeq2) and `stage_de_limma.R` (limma, edgeR; `stage-de-limma` needs `Rscript` on `PATH`, so it runs inside the shell or through `nix develop -c`); its shell hook runs `uv sync --locked`. Enter it
  with direnv (`.envrc` → `use flake`) or `nix develop`. Python packages come
  from uv, not Nix; R packages come from Nix (add them to `flake.nix`).
- Run the CLI: `uv run sead-sexdiff <download|pseudobulk|covariates|technical-covariates|export|analyze|rank-sensitivity|preprocessing-sensitivity|axis-diagnostics|latent-sensitivity|stage-de|stage-de-limma|positive-control|donor-flags|exclusion-sensitivity|figures|subset|de> [options]`.
- `matplotlib` is a direct dependency pinned exactly (`matplotlib==3.10.7`,
  the version scanpy already locked) because `figures` promises
  byte-identical files. Move it only with `uv add "matplotlib==X.Y.Z"`
  (no other locked package may move). A move changes every figure's bytes
  and fingerprint, and `figures --check` then lists every figure as stale.
- `motco` is pinned to an exact PyPI version (`motco==0.18.0` in
  `pyproject.toml`; `uv.lock` holds the wheel's and sdist's sha256). Move it
  only with `uv add "motco==X.Y.Z"`, and before any run under the new pin
  check that the wheel's `motco/` tree is byte-identical to `src/motco/` at
  the release tag `vX.Y.Z` of `https://github.com/Omnia-Lab-UC/MOTCO`
  (`git -C ../MOTCO archive vX.Y.Z src/motco`) and that `uv export` changes
  only the motco line (no other locked package moves); record both checks in
  the re-issue report. Never a git or path dependency (`../MOTCO`). It pins
  scikit-learn exactly; do not loosen that pin to resolve a conflict. `run_analysis` refuses a MOTCO
  older than `ANALYSIS_MOTCO_MIN_VERSION` and `run_latent_sensitivity` one
  older than `LATENT_MOTCO_MIN_VERSION`. A pin move reruns `analyze` (the
  reference checks compare MOTCO versions) after checking the rerun against
  the kept outputs. A pin move that changes the measurement re-issues every
  MOTCO-dependent output (spec `motco-release-migration`): a same-machine
  gate first (new pin with the previous measurement against the previous pin
  from the pre-change commit: byte-identical CSVs, equal records, or stop),
  for a new setting MOTCO documents as moving only named outputs, or only
  adding a named record block, a contract gate on the primaries (new pin
  with and without it: everything else equal, or stop), the superseded
  directories moved with `scripts/move_superseded.py` on
  every machine, every MOTCO-dependent step rerun from scratch in dependency
  order for every cell type, a version-consistency check of every manifest
  outside `results/superseded/` (`results/figures/figures_manifest.json`
  included: its per-figure analyses must record the new version) and of
  every analysis record's measurement setting (an analysis where MOTCO does
  not allow the setting passes only when its manifest records it as not
  available), the
  MOTCO-independent outputs hashed before and after, `figures --report
  <re-issue report stem>` as the last rerun step once every other rerun has
  finished, and one re-issue report with the superseded values beside the
  new ones and its reading written before the reruns. `results/figures/` is
  not moved: it is redrawn, and the superseded figures are those frozen with
  the report that presented them. Figure 01 reads no MOTCO-dependent file,
  so its byte identity with the previous report's frozen copy (same module
  and matplotlib) is a check of the same kind as the MOTCO-independent
  hashes.
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
  - `TAXONOMY`: per-release `Subclass` inventory (and omic) keyed by the
    release date in the file name; every `Subclass` value must be listed,
    otherwise the run fails. The population comes from `CELL_TYPES` (below),
    joined by `taxonomy_for(path, cell_type)`; every `Supertype` of the
    population's subclass must be a population or a recorded non-population
    label, otherwise the run fails. RNA calls microglia `Immune`, ATAC
    `Microglia-PVM`; astrocytes are `Astrocyte` in both. Spike notes:
    `openspec/changes/archive/2026-10-16-mtg-astrocyte-donor-pseudobulks/spike.md`,
    `openspec/changes/archive/2026-10-22-mtg-microglia-positive-control/spike.md`.
  - `MATRIX_SOURCE`, `RNA_MIN_GENES` (RNA nuclei ≥ 200 genes; ATAC no extra QC).
  - `Reference` ADNC donors excluded; `ADNC_STAGES` (`Not AD`+`Low` → `0`,
    `Intermediate` → `1`, `High` → `2`).
  - `COHORT_THRESHOLDS` (≥ 30 / ≥ 50 nuclei of the cell type in both omics)
    and `PREREGISTERED_COHORTS`, per cell type from `CELL_TYPES` (astrocyte
    ≥ 30: n = 80, F 11/10/28, M 9/10/12; ≥ 50: n = 74; microglia, the
    metadata-only spike's counts: ≥ 30: n = 76, F 11/10/24, M 9/11/11; ≥ 50:
    n = 63, F 8/10/18, M 8/10/9). Deviations go to `manifest.json`
    `findings`, never silently.
- Covariate coding (`COVARIATE_COLUMNS`, `COVARIATE_CROSS_CHECK`): RNA values
  are used; `N+` ages are capped and flagged; non-numeric values stay missing
  and are listed; nothing is imputed.
- The nucleus selection of the MTG steps is one code path,
  `select_nuclei(..., cell_type=...)` (`select_astrocyte_nuclei` is its
  astrocyte wrapper; the cell type's subclass, then its population
  supertypes, `Reference` excluded, `RNA_MIN_GENES` on the `UMIs` layer; ATAC
  `read_counts=False` reads `obs` only): `pseudobulk_donor` sums what it
  returns and `technical_donor` summarizes the same nuclei. Changing it
  changes both steps; the pseudobulk outputs must then be re-checked against
  the per-donor cache. The selection's `n_astro` key (population nuclei
  before QC) keeps its name so existing caches stay readable.
- Cell-type registry and positive-control constants (flag and ask before
  changing any of them): `CELL_TYPES` (a registry keyed by canonical name,
  entries added by a change proposal only; the name is a flag, the mapping
  is not): `astrocyte` (the default, `DEFAULT_CELL_TYPE`; `Astrocyte`,
  `Astro_1`…`Astro_6-SEAAD` in both releases; results root `results/`, so no
  astrocyte path moves) and `microglia` (RNA `Immune`, ATAC `Microglia-PVM`,
  the six `Micro-PVM_*` supertypes of both releases with identity mapping,
  `Lymphocyte` and `Monocyte` recorded as non-population; no release has a
  PVM-only supertype; results root `results/microglia/`; labels and cells
  from the metadata-only spike). The rule `results/` for the default,
  `results/<name>/` otherwise, is asserted at import (the cluster scripts
  rely on it); `CELL_TYPE_SUBDIRS` (incl. `motco_latent_sensitivity`,
  `motco_exclusion_sensitivity`) and
  `cell_type_dirs` give the default directories. The microglia reuse of the astrocyte `COHORT_THRESHOLDS`,
  `ADNC_STAGES` and per-cohort `k`, and of every export, analysis and
  stage-DE setting, is a decision of the change. `POSITIVE_CONTROL_CELL_TYPES`
  (control `microglia`, reference `astrocyte`), `POSITIVE_CONTROL_EXPECTATIONS`
  (E1–E3, written before any microglia run, never changed from results; E1
  on the exported RNA set and under both stage-DE methods, a method
  disagreement `method_dependent`, not holding),
  `POSITIVE_CONTROL_TECHNICAL_SOURCES` (E3's technical variables: the
  technical table and the export's TMM / library sizes, not the six analysis
  covariates or `male`), `POSITIVE_CONTROL_PRIMARY_COHORT` (`cohort_30`; the
  control holds when E1–E3 all hold there), strict `>` (ties do not hold),
  the output names (`POSITIVE_CONTROL_*`), and `POSITIVE_CONTROL_BASES` with
  `POSITIVE_CONTROL_DEFAULT_BASIS` (`primary`, today's inputs, path and
  bytes; `latent_adjusted`, the only non-primary basis, tied at import to an
  arm of `LATENT_SENSITIVITY_ARMS` with a latent adjustment: E2 and E3 from
  the arm and its diagnostic, E1 basis-independent, outputs under
  `positive_control/latent_adjusted/<cohort>/`; no restated control for
  `extended_covariates`, whose rank and CV AUROC equal the primary's by
  MOTCO's contract, stated and checked in the report; the registry entry's
  `reading` is the predeclared reading of the restated control, written
  before any microglia arm ran, never changed from results). Reading rules
  for the report:
  the control is about the stage signal, not the sex × stage hypothesis; if
  it holds the astrocyte null is not a pipeline failure, if not the latent
  construction is in question for any cell type; the microglia trajectory
  result is a second case-study result under the existing binding rules
  (with both stage-DE concordance cells), with MOTCO's operating
  characteristics measured on neither cell type's blocks; nothing found in
  microglia changes any astrocyte output or reading. The restated control on
  `latent_adjusted` is a diagnosis of the primary-basis control, which stays
  the control's outcome: if it holds in `cohort_30`, the measured quality
  gradient explains the primary failure, without promoting the arm
  (`candidate_fails` stands) or changing the primary-basis outcome or the
  astrocyte reading; if not, the latent construction is in question beyond
  the quality gradient (composition of the `Micro-PVM_*` bundle, rank-1
  selection, microglia ATAC depth); E3 on that basis carries the mediator
  caveat in both directions and a failure for that reason counts;
  `cohort_50` is descriptive. Since `mtg-rank-floor-primary` the control's
  outcome on the floored primary is the primary-basis outcome, with the
  superseded (rank-1) outcome reported beside it per expectation and cohort
  as dependence on the rank rule; E1–E3 and the all-three rule unchanged.
  Since `mtg-residual-correction-primary` the control's outcome on the
  corrected primary (both bases) is the primary-basis outcome, with the
  superseded (uncorrected) outcome beside it; no expectation reads a p-value,
  so E1–E3's values are compared with the superseded ones (a difference
  beyond rounding is a finding); the restated control reads only the arm's
  CV AUROC and first-component correlations, which the `residualized_basis`
  diagnosis does not bear on, so its predeclared reading stands.
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
  `ANALYSIS_MOTCO_MIN_VERSION` (0.18.0), `ANALYSIS_RANK_FLOOR` (2, MOTCO's
  validated `min_components` default: passed explicitly in every
  `TrajectoryAnalysisParams` the package builds, checked against MOTCO's
  default before any input is read (`check_motco_rank_floor`; a mismatch
  stops the run naming both values), recorded as `measurement.min_components`
  and `--min-components` in `cli_equivalent`, printed in `summary.md`, and
  checked by `_check_reference`, so no step compares against an unfloored
  reference), `ANALYSIS_RESIDUAL_CORRECTION` (`df`: RRPP permutes the reduced
  model's residuals multiplied by κ = √(n / (n − r)), r the reduced model's
  rank, 10 with the six covariates and 12 with the eight of the latent arms,
  so only the null moves (p-values, null summaries, the permutation-pooled
  eigengap); MOTCO's validated opt-in by its verdict `adopt_df`
  (`corrected-permutation-2026-10-30.md`), not its default: one of the two
  MOTCO defaults the package overrides, never a flag; passed explicitly in every
  `TrajectoryAnalysisParams` the package builds, checked in the record after
  each run before anything is written (`measurement.residual_correction` of
  the requested kind with κ, or the run stops naming both settings),
  recorded as the manifest's `measurement.residual_correction` (kind, n,
  reduced rank, κ) and `--residual-correction df` in `cli_equivalent`,
  printed in `summary.md` ("RRPP null: residual correction …"; records
  without one render as before), and checked by `_check_reference`
  (`params.residual_correction`, absent meaning none), so no step compares
  against an uncorrected reference), `ANALYSIS_SPAN_SPLIT_TEST` (`True`:
  MOTCO's split interaction test, its verdict `split_validated`
  (`span-split-2026-11-05.md`), not its default, never a flag: the group ×
  stage interaction energy of the analysed features split into the fitted
  stage-PLS span (dimension = the selected rank), its orthogonal complement
  and all of it, each part tested by a trace pseudo-F under RRPP's draws at
  the analysis's seed, so it only adds the record's `span_split` block;
  passed explicitly in every `TrajectoryAnalysisParams` the package builds
  except a latent-adjusted arm, where MOTCO rejects it
  (`LATENT_SENSITIVITY_SPLIT_UNAVAILABLE`, tied at import to the arms with a
  latent adjustment; their manifests quote MOTCO's reason as
  `measurement.span_split.not_available`); `SPAN_SPLIT_PARTS` checked against
  MOTCO's `INTERACTION_SPLIT_PARTS` before any input is read
  (`check_motco_split_parts`); checked in the record after each run before
  anything is written (`_span_split_problem`: the block present exactly when
  requested, the parts, a pseudo-F p each, the analysis's permutations,
  seed and correction, a `span` dimension equal to the selected rank, a
  `complement` dimension equal to the total minus it, `draws_match_rrpp` at
  `n_jobs = 1`; or the run stops naming the cohort and the entry), recorded
  as the manifest's `measurement.span_split` and `--span-split-test` in
  `cli_equivalent`, printed in `summary.md` (the "Split interaction test"
  table; records without a block render as before), and checked by
  `_check_reference` (`params.span_split_test` and the block), so no step
  compares against a reference without it). MOTCO's other defaults (net-displacement
  orientation, double-CV rank at the floor) are not overridden, and there is
  no unadjusted run. The floored, `df`-corrected analysis is the primary of
  both cell types since `mtg-residual-correction-primary`; the earlier
  floored uncorrected analyses (MOTCO 0.10.0, their p-values withheld by
  MOTCO's roadmap from 2026-10-30) and rank-1 analyses are superseded
  records, not a sensitivity arm. The rank-1 and uncorrected code paths stay
  for analyses made by hand (`_analyze_cohort` with `min_components=1` or
  `residual_correction=None`), and so does the path without the split test
  (`span_split_test=False`, the pin-move gate of `mtg-span-split-test`).
  `--n-jobs` changes MOTCO's RRPP draws (one seed stream per
  worker), so report results from `n_jobs = 1`.
- Split-test policy constants (flag and ask before changing any of them;
  predeclared in the proposal of `mtg-span-split-test`, committed before any
  real-label split test ran, never changed from results): `SPAN_SPLIT_PARTS`
  (`total`, `span`, `complement`), `SPAN_SPLIT_CONFIRMATORY_PART`
  (`complement`: the one confirmatory part, its pseudo-F p-value
  `SPAN_SPLIT_STATISTIC` = `p_value_f` read at `ANALYSIS_ALPHA` under
  `ANALYSIS_REJECTION_RULE`, unadjusted), `SPAN_SPLIT_CONCORDANCE_PART`
  (`span`: a concordance cell against "any of `delta`, `angle`, `shape`
  rejected" in the same analysis, a discordance in either direction a
  finding, never a result), `SPAN_SPLIT_NOT_READ` (`total`, `span_share`,
  `p_value_ss`: reported, not read), `SPAN_SPLIT_ROLES`, and
  `SPAN_SPLIT_COLUMNS` (the split columns of every comparison table, after
  the `shape` columns: F and p per part, `span_rejected`, `span_concordant`,
  `complement_rejected`, `span_share`; empty for a row without the test; no
  `total` status). No adjustment across `delta`, `angle`, `shape` and the
  complement, or across cohorts and cell types: each primary carries four
  tests read at α, stated in the report. Reading rules for the report: each
  cell type's `cohort_30` primary gives its split-test result; `cohort_50`
  and every other analysis are sensitivity, a complement status that differs
  from the primary's being rank, preprocessing, covariate-set or cohort
  dependence (not a result of that analysis); a complement rejection is a
  sex × stage difference with interaction energy outside the fitted
  stage-PLS span, not a different molecular program (it may hold leakage
  from an in-span difference, 0.32 of replicates for a 45° in-span rotation
  on the SEA-AD table, and what the fitted span missed of the true stage
  program, about 0.10 kept at SEA-AD's stage signal), with no geometry,
  direction, feature list or attribution, and it starts no attribution; a
  complement non-rejection is uninformative (no power shown at SEA-AD's
  stage signal for any plant, nor for sparse differences: 0.02–0.05 for 20
  features at 2–3 within-stage SDs; the dense-rotation power 0.68 was
  measured at a strong stage signal); the in-span share is reported beside
  the rank and never read against 0. Carried limitations: size measured on
  real noise only for the astrocyte `cohort_30` export at this
  configuration (complement 0.012, span 0.050, total 0.040), inferred for
  `cohort_50`, microglia, the eight-covariate arm, the exclusion variants,
  the forced ranks and the preprocessing variants; out-of-span power shown
  only for a dense rotation at a strong stage signal; the null in-span share
  grows with the fitted rank (0.236 at rank 5 on this export) and no null
  share is computed here; the test is not available on the `latent_adjusted`
  arm. Nothing else is re-read: the stage-DE concordance rules (which read
  `delta`), the positive control (no p-value), the adoption, mediator, flag
  and axis checks and the `latent_adjusted` reading are unchanged.
- Rank-sensitivity constants (flag and ask before changing any of them):
  `ANALYSIS_FORCED_RANKS` (2, 3, 4, 6, 9, 12: 2 = `n_stages − 1`, 3 the next
  rank, 4/6/9/12 the steps of MOTCO's latent-rank ladder with 4 = 2(K − 1)
  the forced rank with a Type I check; a constant, never a flag, never
  chosen from results), the same `ANALYSIS_SEED` and `ANALYSIS_PERMUTATIONS`
  at every rank (so every rank draws the same permutation indices), rank 1
  not forceable (MOTCO accepts `forced_components ≥ 2`), the CV-selected
  analysis (rank ≥ the floor) as the reference row (the forced run at its
  rank, when in the ladder, must equal it in statistics, p-values and split
  parts; otherwise the manifest says the reference row stands in for it),
  the split test at every rank (each rank's `span` part its own rank, a
  complement status that differs from the reference's rank dependence of the
  split test, the in-span share reported beside the rank and not read),
  `n_jobs = 1` for reported results, and the rule that the CV-selected
  analysis stays primary.
- Reporting MOTCO results follows the binding rules in
  `../MOTCO/docs/roadmap.md` (Phase 6, Next changes item 1): CV-selected rank
  beside `angle`; `angle`/`shape` non-rejections are uninformative; a
  `delta` rejection is a difference, not specifically of size; report the
  covariate block and the feature filtering beside the results. The rank
  floor defines `angle` and `shape` in at least two dimensions without making
  them powerful (two dimensions, not a second axis with signal): give each
  analysis's `config_spectrum` eigengap and null summary beside `angle`. The
  forced-rank sweep is read as sensitivity only: the CV-selected analysis
  stays primary, a rejection status that differs between ranks is rank
  dependence (not a result at that rank), the non-rejection rules hold at
  every rank, MOTCO's forced-rank operating characteristics (Type I at
  forced 4 only, unadjusted; `shape` power loss and `angle` leakage at
  forced 2(K − 1); wider `angle` null above `n_stages − 1`) are carried as
  limitations, and the adjusted analysis has no Type I check at any forced
  rank. Under the `df` residual correction the p-values are size-controlled
  at the configuration MOTCO measured (six covariates, both `cohort_30`
  tables; `corrected-permutation-2026-10-30.md`), with `delta` conservative
  there (0.029 / 0.031 pooled), so a `delta` non-rejection stays weak
  evidence; the size is inferred from the mechanism (n and the reduced rank),
  not measured, for `cohort_50`, the eight-covariate latent arms, the
  exclusion variants, the forced ranks and the preprocessing variants, each
  carried as a limitation.
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
- Latent-adjustment-sensitivity constants (flag and ask before changing any
  of them): `LATENT_SENSITIVITY_ADDED_COVARIATES` (`ribosomal_pseudobulk`,
  `mitochondrial_pseudobulk` on the raw fraction scale of the technical
  table, the `residualized` variant's scale; tied at import to
  `RESIDUALIZATION_COVARIATES` minus `log10_lib_*`, so a change to one set
  is a change to both; the export's `log10_lib_<omic>` stay out),
  `LATENT_SENSITIVITY_COVARIATES` (`ANALYSIS_COVARIATES` + the two, in this
  order), `LATENT_SENSITIVITY_CONDITIONING` (`stage`, MOTCO's admissible
  conditioning; `group_stage` is never run), `LATENT_SENSITIVITY_ARMS`
  (`extended_covariates`, `latent_adjusted`; a registry, never a flag, never
  chosen from results, the same for every cell type the step runs on),
  `LATENT_SENSITIVITY_SPLIT_UNAVAILABLE` (`latent_adjusted`: exactly the
  arms with a latent adjustment, tied at import; MOTCO rejects the split
  test there, so the arm requests it off, its record carries no
  `span_split` block, its manifest quotes MOTCO's reason as
  `measurement.span_split.not_available` and its split cells are empty;
  `extended_covariates` carries the test),
  `LATENT_MOTCO_MIN_VERSION` (0.9.0), the mediator
  check (`mediator_check`, the `DIAGNOSTIC_CORRELATION_THRESHOLD` reused,
  not partial on the covariates, written before the arms, never a gate),
  the output names (`LATENT_SENSITIVITY_*`), `n_jobs = 1` for reported
  results, and the rule that the primary analysis stays primary. An arm
  manifest records `measurement.latent_adjustment` and a `cli_equivalent`
  with `--latent-adjustment` only when the adjustment is set (absent for the
  primary, as MOTCO omits the setting from its record), the technical table
  under `inputs.technical`, the fractions' source under
  `covariates.sources` and the cell type; the comparison manifest records
  the cell type and the cell-type findings of the manifests read. Reading
  rules for the report (one per cell type the arms run on): MOTCO's power
  comparison (`candidate_fails`) excludes the latent adjustment as a
  replacement, so the arms are sensitivity only; a rejection status that
  differs is covariate-set dependence (primary vs `extended_covariates`) or
  basis dependence (`extended_covariates` vs `latent_adjusted`), not a
  result of the arm (a complement status that differs between primary and
  `extended_covariates` is covariate-set dependence of the split test); a
  CV-selected rank that differs is reported beside
  `angle`; the non-rejection rules hold in every arm; the `latent_adjusted`
  arm's `angle` and `shape` p-values are reported but not read (MOTCO's
  case-study diagnosis `residualized_basis`: the `stage` latent adjustment is
  liberal on both at this configuration, an excess of the basis the residual
  correction does not address, so neither a rejection nor a non-rejection is
  interpretable; its `delta` is read under the rules above); an arm's
  first-axis technical correlations (`axis-diagnostics` on the arm, whose
  tables label it `primary`) are a measurement-space finding only. Carried
  limitations: the mediator limitation with the mediator check's flags (no
  simulation can tell a confounder from a mediator), MOTCO's Phase 6-only
  operating characteristics (Type I `adopt_stage`; power `candidate_fails`,
  `shape` power lost under a nuisance axis, small-magnitude power lost
  without one), the `residualized_basis` diagnosis (anchor `angle`
  0.084–0.090, `shape` 0.078–0.092 against 0.0695), no operating
  characteristic measured on these blocks or at eight covariates (the
  residual-corrected null's size measured at six covariates only, inferred at
  eight, reduced rank 12), the split test's size inferred, not measured, at
  eight covariates and not available on `latent_adjusted`, and the nuclei
  counts and the fractions both proxying library quality.
- Stage-DE constants (flag and ask before changing any of them):
  `STAGE_DE_MODELS` (a registry, never a flag, never chosen from results:
  `factor_additive` `~ sex + stage` testing `stage[T.1]`, `stage[T.2]`;
  `factor_interaction` `+ sex:stage` testing `sex[T.Male]:stage[T.1]`,
  `sex[T.Male]:stage[T.2]`, the primary interaction model;
  `trend_additive` `~ sex + stage_num` testing `stage_num`;
  `trend_interaction` `+ sex:stage_num` testing `sex[T.Male]:stage_num`, a
  second predeclared model, not a fallback; the two Wald tests of a factor
  model also counted as their `union` (`STAGE_DE_UNION`), with no joint
  error-rate claim, since pydeseq2 0.5.3 has no LRT),
  `STAGE_DE_REFERENCE_LEVELS` (`Female`, stage `0`; treatment coding via the
  categorical order), `STAGE_DE_COVARIATES` (per omic `age_at_death`, `pmi`,
  `rin`, `years_education` and the omic's own `log10_n_<omic>`, continuous,
  recorded scale, never imputed; tied at import to `ANALYSIS_COVARIATES`
  minus the other omic's count), `STAGE_DE_FEATURE_STATUSES` (`variance`,
  `kept`: the export-bound tested set; size factors and dispersions fitted
  on it; the export's TMM not reused) with `STAGE_DE_EXPORTED_STATUS`
  (`kept`) and `STAGE_DE_FEATURE_SETS` (`exported` the primary reading,
  `abundance`; `padj` never recomputed on a subset), `STAGE_DE_ALPHA` (0.05,
  `padj ≤ α`, pydeseq2's BH with independent filtering and Cook's defaults,
  no LFC shrinkage, no permutation null), `STAGE_DE_CONCORDANCE_RULE`
  (exported `factor_interaction` union count > 0 against the primary `delta`
  rejection status; reads a count and a status, never a p-value magnitude;
  a discordance is a finding), the reproduction tolerances
  `STAGE_DE_TMM_RTOL` (1e-10) and `STAGE_DE_VALUES_RTOL` (1e-9), the output
  names (`STAGE_DE_*`), and the rule that the primary analysis stays
  primary. `--n-cpus` does not change the tables, but warnings raised in pydeseq2's worker processes are not
  captured when it is above 1 (recorded as `warnings_scope` and in the `n_cpus` finding). Reading rules for the report: the
  stage counts describe the stage signal in the pseudobulks; the
  interaction counts and the concordance cells cross-check the primary
  analysis and are not a second test of its hypothesis; a discordance in
  either direction is reported as such (per-feature FDR and one latent-space
  permutation test ask related, not identical, questions); a zero
  interaction count is "none detected at this n"; the union counts carry no
  joint error-rate claim; the MOTCO non-rejection rules are untouched; the
  covariate set, the feature set and DESeq2's own normalization are stated
  beside the counts; the FDR of `padj ≤ α` is approximate and Cook's
  handling is inactive with continuous covariates. Changes to the fit must be
  re-checked with `scripts/check_stage_de_deseq2.R` (full DESeq2 fit with
  pydeseq2's trend type plus a plug-in Wald test on pydeseq2's dispersions;
  exit 0 / 2 / 3 / 4, see README; the full fit is not expected to meet 1e-3
  on every feature).
- Stage-DE limma constants (flag and ask before changing any of them):
  `STAGE_DE_LIMMA_METHOD` (edgeR TMM on the tested counts, library size =
  their column sums, `TMM_LOGRATIO_TRIM` / `TMM_SUM_TRIM` (tied at import),
  `doWeighting`, `Acutoff` −1e10, automatic reference; `voomWithQualityWeights`
  with `normalize.method = "none"`, `span` 0.5, `adaptive.span`, no
  `var.design` / `var.group`, `genebygene`, `maxiter` 50, `tol` 1e-5,
  `arrayWeights`' `prior.n` 10; `lmFit` `ls`; `eBayes` `proportion` 0.01,
  `stdev.coef.lim` (0.1, 4), `trend` / `robust` false, `legacy` NULL; BH over
  all tested features, no independent filter: limma 3.68.4 / edgeR 4.10.1
  defaults, which the R script checks against the installed ones; never a
  flag, never changed from results), `STAGE_DE_LIMMA_JOINT_TESTS`
  (`stage_joint` for `factor_additive`, `interaction_joint` for
  `factor_interaction`, tied at import to the factor models' tested labels;
  limma's F is the reported joint test and the exact-F diagnostic only
  measures its approximation), `STAGE_DE_LIMMA_CONCORDANCE_RULE` (exported
  `interaction_joint` count > 0 against the primary `delta` rejection
  status, tied at import to `STAGE_DE_CONCORDANCE_RULE`'s model and set;
  reported beside the DESeq2 cell with `methods_agree`),
  `STAGE_DE_LIMMA_TMM_ATOL` (1e-8 absolute, the `check_tmm_edger.R`
  tolerance), the output names (`STAGE_DE_LIMMA_*`), the R bridge (one
  packaged script through `Rscript`, one BLAS thread; its sha256 and the R /
  limma / edgeR / jsonlite versions in every fingerprint) and the rule that
  the DESeq2 results and the primary analysis stay primary. Reading rules for
  the report: the change was proposed after the DESeq2 results were known,
  and no later change may promote whichever method finds more; a call that
  differs between the methods is method dependence, not a result of either,
  and agreement is not replication (same donors and counts); a zero count is
  "none detected at this n"; the union counts carry no joint error-rate
  claim and the joint F's claim is approximate; the concordance cells
  cross-check the primary analysis and are not a second test of its
  hypothesis, whose non-rejection rules are untouched; voom's log-CPM, the
  BH family and the normalization (the export's TMM, not DESeq2's median of
  ratios) are stated beside the counts. The calibration of the change
  (simulations under DESeq2's own model, rebuilt; realized FDR of the
  additive tables 0.054 / 0.038 against DESeq2's 0.078 / 0.081) is not a
  gate on real data.
- Donor-flag constants (flag and ask before changing any of them):
  `DONOR_FLAG_COLUMNS` (the flag and the five staging labels, coded names),
  `DONOR_FLAG_ATAC_ONLY_COLUMNS` (the pseudo-progression score, read where
  present, its absence a finding), `DONOR_FLAG_SOURCE` with
  `DONOR_FLAG_VALUES` (`Y` / `N`; any other value an error, never worked
  around: the proposal is amended first), `DONOR_FLAG_REFERENCE_LABEL`
  (ADNC `Reference`, the donor table's rule), `DONOR_FLAG_EXPECTED_FLAGGED`
  (11, a finding when different), `DONOR_FLAG_OUTPUTS`, the labels and the
  score recorded and read by nothing (CPS-based staging is a later
  proposal), and the rule that the table stays apart from
  `donor_covariates.csv` / `covariates_manifest.json` and is an input of
  `exclusion-sensitivity` only: never a covariate, a diagnostic variable or
  a filter of any other step.
- Donor-exclusion constants (flag and ask before changing any of them):
  `DONOR_EXCLUSION_VARIANTS` (a registry keyed by variant name with a
  `description`, the flag column and `exclude_when`; now only
  `excluding_severely_affected`, SEA-AD's own handling of expression tests;
  never a flag, never chosen from results, entries added by a change
  proposal only), the cohort thresholds kept with `k` recomputed by the
  export's rule (never a change to `COHORT_THRESHOLDS`, `ADNC_STAGES` or
  `PREREGISTERED_COHORTS`), the flag check (`DONOR_EXCLUSION_CHECK_VARIABLES`,
  `DONOR_EXCLUSION_COMPONENT` `pls_0`, `DONOR_EXCLUSION_FLAG_CHECK_STAGE` `2`;
  rank-biserial AUC and Spearman, no p-value, written before any run, never a
  gate), `DONOR_EXCLUSION_AXIS_RULE` (E3's technical sources, strict `>`,
  the control cell type's value the expectation), `DONOR_EXCLUSION_READING`
  (predeclared in the proposal before any run), the output names
  (`DONOR_EXCLUSION_*`), the stage DE on a variant (`stage-de --omic both`,
  so the limma step's DESeq2 binding holds, and `stage-de-limma --omic rna`;
  the RNA counts are the reading, the ATAC DESeq2 tables are recorded only),
  `n_jobs = 1` for reported results, and the rule that the primary analyses,
  the primary cohorts and the positive control's outcome are unchanged by
  anything found here. Reading rules for the report (both cell types): the
  variant is sensitivity only, never a revised primary; a rejection status
  that differs from the primary's is cohort dependence, not a result of the
  variant; a CV-selected rank that differs is reported beside `angle`; the
  non-rejection rules hold in the variant; the cross-cell-type comparisons
  (variant CV AUROC microglia against astrocytes; each cell type's variant
  exported RNA stage counts against its primary counts) are compared, not
  tested. Carried limitations: the excluded donors carry the most severe
  pathology, so the exclusion truncates the late epoch with the quality tail
  and a stage effect concentrated in them goes with them; the flag was made
  from QC metrics, so the variant conditions on quality by donor selection;
  the stage-2 cells are smaller (less power); MOTCO's operating
  characteristics were measured on neither cohort; the flag's sex split
  bears on the sex × stage comparison.
- Figure constants (flag and ask before changing any of them):
  - `FIGURES`: a registry keyed by name in drawing order, with number,
    title, draw function, caption builder, fixed construction text and
    reading keys. No flag selects, skips or adds a figure. Entries are added
    or changed by a change proposal only, and a change that adds a step adds
    its figure here, so the step's report freezes it.
  - The display choices:
    - `FIGURE_COHORT` (`cohort_30`, tied at import to
      `POSITIVE_CONTROL_PRIMARY_COHORT`; the cohort of figures 03, 07, 08
      and 10, with `cohort_50` wherever a figure draws both);
    - `FIGURE_STAGE_DE` (`rna`, `factor_additive`, `stage_2_vs_0`, checked
      at import against `STAGE_DE_MODELS`; the caption gives the
      `stage_1_vs_0`, union and limma `stage_joint` counts too);
    - `FIGURE_FIRST_AXIS_RULE` (figure 08's variable: E3's rule over
      `POSITIVE_CONTROL_TECHNICAL_SOURCES`, values from the technical table
      or the export manifest);
    - `DIAGNOSTIC_CORRELATION_THRESHOLD` reused for the heatmap labels;
    - `FIGURE_SENSITIVITY_CELL_TYPE` (astrocytes, the cell type
      `rank-sensitivity` and `preprocessing-sensitivity` run on);
    - `FIGURE_PREPROCESSING_ROOT` (tied at import to
      `run_preprocessing_sensitivity`'s default);
    - figure 11 (`11_span_split`, from `mtg-span-split-test`): the four
      primaries' pseudo-F p per part marked by role, the complement p of
      every registered analysis carrying the test (the `latent_adjusted`
      arms listed as not available), the in-span share against the fitted
      rank with no reference line; not added to figure 05, whose
      construction stays frozen.
  - Rendering: `FIGURE_PALETTE` (sex and cell-type pairs validated; the
    aqua is below 3:1 on white, so every cell-type mark has a legend entry
    or a direct label), `FIGURE_DPI` (200), `FIGURE_FORMATS`,
    `FIGURE_RASTER_MIN_POINTS` (500), `FIGURE_SVG_HASHSALT`, `FIGURE_FONT`
    (matplotlib's bundled DejaVu Sans, checked at run time),
    `FIGURE_JITTER_SEEDS`, and the reproduction tolerances
    `FIGURE_STAGE_MEANS_ATOL` (1e-8) and `FIGURE_SPEARMAN_ATOL` (1e-9).
  - `FIGURE_READING_RULES`: the caption texts that restate the rules of this
    file (the binding MOTCO rules, the `df` size statement, sensitivity as
    dependence, the forced-rank limitations, `latent_adjusted` not read,
    measurement-space findings, the diagnostic rule, the stage-DE count
    rules, the positive control's bases, the split test's policy and
    reading, `span_split`). A rule changed here by a later
    change changes its text there in the same change.
  - The output names (`FIGURE_*`).

  Reading rules: the figures are a display step. They test nothing, are an
  input of no step, and change no reading of the reports that made the
  outputs. A caption contains the fixed construction, values read from the
  files and the registered rules, and no other interpretation. A caption
  value that differs from the report that made the output is a finding about
  that report, recorded with both values.
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
  node, as does `donor-flags --omic both` (once per cell type's pseudobulk
  directory, `--cell-type <name>`; then `scripts/sync_pseudobulk.sh`). `technical-covariates --omic both` needs the raw objects and
  `results/pseudobulk/` (incl. the per-donor cache): `scripts/technical.sbatch`
  after the pseudobulk job, or nohup on the login node; resubmitting resumes.
- Another cell type: both sbatch scripts read `CELL_TYPE` (default
  `astrocyte`) and pass `--cell-type`; their `#SBATCH` log paths are the
  astrocyte ones, so submit with `--export=ALL,CELL_TYPE=<name>` and
  `--output` / `--error` under `results/<name>/pseudobulk/logs/` (create it
  first). The scripts derive the pseudobulk directory from the rule
  `results/` for astrocyte, `results/<name>/` otherwise, which `data.py`
  asserts for every `CELL_TYPES` entry at import.
- Workstations pull only a cell type's pseudobulk directory with
  `scripts/sync_pseudobulk.sh [--cell-type <name>]`; never copy raw donor
  objects around.

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

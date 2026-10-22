# Proposal

## Why

This proposal is written **after** the DESeq2 results of `mtg-astrocyte-stage-de` were known, and they were null for the interaction. No feature reached padj ≤ 0.05 on any interaction coefficient in either cohort or omic; the smallest interaction padj was 0.068–0.42. Stage signal was small: 1–2 RNA genes (`CREB5`, `ADAM33`) and 22–36 ATAC peaks at stage 2 vs 0. All four concordance cells were concordant (`reports/motco-astrocyte-stage-de-2026-10-21.md`). Nothing here re-tests that result or looks for a method that finds more. The DESeq2 results stay the reported results of that change.

The DESeq2 cross-check has four known gaps, each stated in its report:
- **No joint interaction test.** pydeseq2 0.5.3 has no likelihood-ratio test. `factor_interaction`'s two coefficients get separate Wald tests, and the concordance rule reads their union, which carries no joint error rate.
- **Approximate FDR.** In the change's simulations under DESeq2's own model, the realized FDR at padj ≤ 0.05 was about 0.06–0.09.
- **No outlier protection.** Cook's handling is inactive with continuous covariates, so no donor is protected against.
- **Strongly varying library quality.** In cohort_30, RNA library sizes span 45× and ATAC 23×. The export's TMM factors run 0.80–1.18 (RNA) and 0.73–1.24 (ATAC). The primary analysis's first latent axis tracks the ribosomal fraction (ρ −0.84) and the RNA TMM factor (−0.73) (`reports/motco-astrocyte-technical-axis-2026-10-19.md`).

limma-voom addresses each gap:
- a moderated F-test of both coefficients jointly;
- the export's own TMM normalization, which is also MOTCO's, instead of DESeq2's median of ratios;
- per-donor quality weights that down-weight noisy donors without removing them.

Why now: `mtg-microglia-positive-control` reads `stage-de` (its expectation E1 and the concordance rule). Merging this change first lets that change predeclare both methods before any microglia run.

## What Changes

- **New `stage-de-limma` step**, per cohort (`cohort_30` primary, `cohort_50` sensitivity) and omic. It sits beside `stage-de`, never replaces it, and leaves `stage-de`'s code path and outputs unchanged. It reuses, unchanged:
  - `STAGE_DE_MODELS` (the four formulas and tested coefficients);
  - `STAGE_DE_REFERENCE_LEVELS`, `STAGE_DE_COVARIATES`, `STAGE_DE_ALPHA`;
  - the export-bound feature set (`features_<omic>.csv` status `variance` / `kept`; `exported` = `kept`) and the two reporting sets;
  - `stage-de`'s reference checks and its reproduction of the export from the counts.

  Each model's design matrix is the one the DESeq2 step fitted, with the same column names.
- **Method, predeclared before any limma run** (constants, never flags):
  - **Normalization:** edgeR TMM with the export's `TMM_LOGRATIO_TRIM` / `TMM_SUM_TRIM`, fitted on the tested set. The factors must equal the export manifest's TMM factors to 1e-8, or the step stops.
  - **Weights:** `voomWithQualityWeights`, which combines voom's observation weights with one weight per donor (decided here; justified in design D4).
  - **Fit:** `lmFit`, then `eBayes` with limma's defaults.
- **Tests:**
  - **Moderated t:** one per tested coefficient.
  - **Joint moderated F:** one per factor model, over its two coefficients: `stage_joint` for `factor_additive` (`stage[T.1]`, `stage[T.2]`) and `interaction_joint` for `factor_interaction` (the two `sex[T.Male]:stage` terms).
  - **Union:** each factor model's two moderated t-tests are also counted as a union, so the counts can be compared like for like with DESeq2. As in `stage-de`, the union carries no joint error rate.
  - **Significance:** BH at `STAGE_DE_ALPHA` (0.05) over all tested features.
- **Limma concordance rule** `STAGE_DE_LIMMA_CONCORDANCE_RULE`, predeclared here. `interaction_detected` is true when at least one exported feature has adj.P ≤ α on `interaction_joint`. It is compared with the primary analysis's `delta` rejection status, read as `stage-de` reads it.
  - The cell sits beside the DESeq2 cell (the union of two Wald tests), which stays the reported cell of `mtg-astrocyte-stage-de`. A `methods_agree` flag compares the two `interaction_detected` values.
  - A disagreement is written as a finding of method dependence. Neither cell replaces the other or changes a MOTCO reading.
- **Method comparison:** `method_comparison.csv` compares the two methods per coefficient, union and joint test, on both feature sets:
  - the counts per method, the overlap and the features only one method calls;
  - whether each method detects anything;
  - Spearman correlations of the log2 fold changes and of the test statistics.

  It reads the DESeq2 outputs after checking them: their manifest must record the same input hashes, their tables must hash as recorded, and their summary must equal a recomputation from the tables. The DESeq2 directory is never written.
- **R through `Rscript`**, not rpy2. A packaged script, `src/sead_sexdiff/stage_de_limma.R`, exchanges files with Python. Its sha256 and the R, limma and edgeR versions enter every fit record. `rPackages.limma` is added to `flake.nix` explicitly (today it comes only as an edgeR dependency).
- **Outputs** under `results/stage_de_limma/<cohort>/`:
  - per tested coefficient and joint test, `<omic>/<model>__<label>.csv`;
  - per fit, `<omic>/<model>__fit.json`, the resume record;
  - per cohort, `stage_de_limma_summary.csv`, `method_comparison.csv`, `stage_de_limma.md` and `stage_de_limma_manifest.json`.

  The step reads sex and stage, so its outputs feed no export, analysis or sensitivity step.
- **Verification before real data:**
  - **Planted effects:** run on the synthetic dataset of `mtg-astrocyte-stage-de` (planted stage effect and female-only interaction; rebuilt if absent).
  - **Calibration:** limma's realized FDR is measured on that change's calibration simulations and set beside DESeq2's.
  - **Stop point:** if limma's FDR exceeds DESeq2's on the homogeneous simulations, the change stops and asks before reading real data.
- **A committed report**, `reports/motco-astrocyte-stage-de-limma-<date>.md`, compares the methods. Its numbers are pasted from the generated files. It carries the reading rules:
  - a call that differs between methods is method dependence, not a result;
  - a zero count means "none detected at this n";
  - the DESeq2 results stay the reported ones;
  - the concordance cells cross-check the primary analysis and are not a second test of its hypothesis;
  - the MOTCO non-rejection rules are untouched.
- README and AGENTS.md document the command, the constants and the reading rules.

Out of scope:
- replacing DESeq2, or changing `stage-de`, its outputs or its report;
- a plain-voom arm (run only inside the synthetic calibration, as the comparison for the quality-weights decision);
- robust eBayes, `treat` and edgeR quasi-likelihood;
- `duplicateCorrelation` / `dream`, since there is one pseudobulk per donor;
- any upstream use of the results;
- running on microglia (that change is amended to add it, see Impact).

## Capabilities

### New Capabilities
- `pseudobulk-stage-de-limma`: the limma-voom fits of the registered stage-DE models on the export-bound feature set, with the predeclared method (TMM → voomWithQualityWeights → lmFit → eBayes), the per-coefficient moderated t and joint moderated F tests, the limma concordance rule beside the DESeq2 cell, the method comparison against the checked `stage-de` outputs, and the committed comparison report.

### Modified Capabilities
<!-- None. `pseudobulk-stage-de`, `motco-input-export`, `motco-trajectory-analysis` and `astrocyte-donor-pseudobulk` are read and checked, never written; `stage-de`'s outputs stay byte-identical. -->

## Impact

- **`src/sead_sexdiff/data.py`:**
  - **Constants:** `STAGE_DE_LIMMA_METHOD` (the settings), `STAGE_DE_LIMMA_JOINT_TESTS`, `STAGE_DE_LIMMA_CONCORDANCE_RULE`, `STAGE_DE_LIMMA_TMM_ATOL`, the `STAGE_DE_LIMMA_*` output names, and import-time ties to `STAGE_DE_MODELS` and `TMM_*`.
  - **Functions:** `stage_de_limma_fit`, `run_stage_de_limma`, `render_stage_de_limma`.
  - **Reuse:** `_stage_de_prepare`, `stage_de_design` and `stage_de_design_check` are reused without changing behaviour. `stage-de` must stay byte-identical on the synthetic data.
- **`src/sead_sexdiff/stage_de_limma.R`** (new, packaged).
- **`src/sead_sexdiff/cli.py`:** new `stage-de-limma` subcommand (`--cohort`, `--omic`, `--pseudobulk-dir`, `--input-dir`, `--analysis-dir`, `--deseq2-dir`, `--out-dir`, `--force`).
- **`flake.nix`:** `rPackages.limma`. Installed versions: R 4.6.1, limma 3.68.4, edgeR 4.10.1. **Dependencies:** no Python dependency added. rpy2 stays declared and unused.
- **Compute:** four fits per cohort × omic, run on the workstation that holds `results/pseudobulk/`, `results/motco_input/`, `results/motco_analysis/` and `results/stage_de/` (nixos-dell-laptop has all four). These are seconds for RNA. For ATAC (~129k peaks), the quality-weight iterations and the exact-F diagnostic are expected to take minutes.
- **Ordering and the microglia change:** `mtg-microglia-positive-control` (branch `origin/mtg-microglia-positive-control`, not on main) needs an amendment before it runs, made after this change merges and outside this change:
  - `--cell-type` and the `cell_type` manifest check on `stage-de-limma` too (eight MTG commands, not seven);
  - `stage-de-limma` in its D10 byte-identity baseline;
  - the limma concordance rule applied to microglia beside the DESeq2 one;
  - a predeclared decision on whether E1 counts DESeq2's `factor_additive` union (as written), limma's `stage_joint`, or both, and how a disagreement enters the all-three rule.
- **Scientific choices introduced here** (flag and ask before changing):
  - the method chain and its settings, including quality weights and eBayes defaults;
  - the joint tests and their labels;
  - the BH family (all tested features, no independent filter);
  - the TMM tolerance;
  - the limma concordance rule and its place beside the DESeq2 cell;
  - the method-dependence reading rule;
  - the rule that the DESeq2 results and the primary MOTCO analysis stay primary.

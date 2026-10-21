# Proposal

## Why

Nothing is rejected in the astrocyte MOTCO analysis or in any of its sensitivity analyses, and the latent-adjustment arm shows why: once the library-quality axis is removed, the single CV-selected component tracks stage (Spearman −0.64 / −0.75) but the female and male stage means along it overlap (cohort_30: F −36.4 → −37.6 → −60.6, M −34.8 → −40.4 → −60.4; `reports/motco-astrocyte-latent-adjustment-2026-10-21.md`), and the CV AUROC stays at 0.63. We do not know whether that is because the astrocyte pseudobulks carry little stage signal at all, or because they carry a stage signal that is the same in both sexes. Both readings predict the same MOTCO output. A per-feature model on the same donors and features separates them: a count of genes and peaks that move with ADNC stage, and a count whose stage effect differs by sex.

The second count also gives MOTCO a predeclared cross-check. MOTCO's `delta` is a latent-space summary of the sex × stage interaction over the exported features, so a cohort whose exported features show a sex × stage interaction at FDR 0.05 is expected to reject `delta`, and a cohort that rejects `delta` is expected to show one. The expectation is written here before either count exists. Discordance is a finding to report, not a result to resolve by choosing one side.

## What Changes

- **New `stage-de` step**, per cohort (`cohort_30` primary, `cohort_50` sensitivity) and per omic (`rna`, `atac`): DESeq2 (pydeseq2) on the raw pseudobulk counts of the cohort's donors, with the covariate-adjusted models of a predeclared registry `STAGE_DE_MODELS`:
  - `factor_additive`: `~ sex + stage + covariates`, stage a three-level factor (reference `0`), sex reference `Female`; Wald tests of the two stage coefficients (`1 vs 0`, `2 vs 0`).
  - `factor_interaction`: `~ sex + stage + sex:stage + covariates`; Wald tests of the two interaction coefficients. This is the primary interaction model, because it matches MOTCO's ordered-categorical stages.
  - `trend_additive`: `~ sex + stage_num + covariates`, stage as 0/1/2; Wald test of the one stage coefficient.
  - `trend_interaction`: `~ sex + stage_num + sex:stage_num + covariates`; Wald test of the one interaction coefficient. A second, predeclared model, not a fallback.
  - pydeseq2 0.5.3 has Wald tests only (no likelihood-ratio test), so a factor model gives one test per coefficient. A feature counts as stage- or interaction-significant under a factor model when either coefficient passes; the union is reported as such and carries no joint error-rate claim.
- **Covariates** `STAGE_DE_COVARIATES` per omic: `age_at_death`, `pmi`, `rin`, `years_education` and the omic's own `log10_n_<omic>` (the analysis covariates minus the other omic's nuclei count), continuous, taken from the primary export's `covariates.csv`. Nothing is imputed; a missing value is an error.
- **Feature set bound to the primary export.** The counts tested are the features that passed the export's chromosome and abundance filters (`features_<omic>.csv` status `variance` or `kept`), so DESeq2's size factors and dispersions are fitted on the same abundance-passing set MOTCO's TMM was computed on. Each result row carries `exported` (status `kept`: one of MOTCO's 5,000 features). Counts are reported on the exported set (the features MOTCO measured, the primary reading) and on the abundance-passing set.
- **Significance:** DESeq2's Benjamini–Hochberg `padj` at `STAGE_DE_ALPHA = 0.05`, with pydeseq2's defaults (independent filtering, Cook's refit). No permutation null.
- **Predeclared concordance rule** `STAGE_DE_CONCORDANCE_RULE`, evaluated per cohort and omic on the exported set from the generated files: `interaction_detected` = at least one exported feature with `padj ≤ α` on an interaction coefficient of `factor_interaction`; compared with the primary analysis's `delta` rejection status read from its `result.json`. Concordant when both are true or both are false; otherwise discordant, written as a finding. The rule reads a count and a status, never a p-value of MOTCO's. The primary analysis stays primary: a DE result changes no MOTCO reading, and a discordance is reported beside both.
- **Reference checks before anything runs:** the primary export manifest records the sha256 of the pseudobulk files on disk; the primary analysis manifest records the sha256 of the export files on disk, the installed MOTCO version, `ANALYSIS_SEED` and `ANALYSIS_PERMUTATIONS` (as `rank-sensitivity` checks them). Any difference stops the step.
- **Outputs** under `results/stage_de/<cohort>/`: `<omic>/<model>__<term>.csv` (one DESeq2 result table per tested coefficient: `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`, `exported`), `stage_de_summary.csv` (counts per cohort × omic × model × term × feature set, with the concordance cells), `stage_de.md`, `stage_de_manifest.json` (inputs with hashes, models, covariates, α, reference levels, pydeseq2 version, git).
- The step reads sex and stage, so it is a check beside the analysis, never a filter: its outputs are an input of no `export`, `analyze` or sensitivity run, and no feature set is derived from it.
- **A committed report**, `reports/motco-astrocyte-stage-de-<date>.md`, numbers pasted from the generated files and checked against them, with the reading rules: the stage counts are a description of the stage signal in the pseudobulks; the interaction counts and the concordance cells are a cross-check of the primary analysis, not a second test of the hypothesis; a discordance in either direction is reported as such (per-feature FDR over thousands of features and a single latent-space permutation test ask related, not identical, questions); the MOTCO non-rejection rules are untouched.
- README and AGENTS.md document the command, the constants and the reading rules.

Out of scope: a likelihood-ratio test (not in pydeseq2 0.5.3); limma-voom or edgeR alternatives; gene-set or peak-annotation enrichment; any change to the export, the analysis or their feature sets; running the step on another cell type (the `mtg-microglia-positive-control` change does that through this step's directory flags).

## Capabilities

### New Capabilities
- `pseudobulk-stage-de`: the covariate-adjusted per-feature DESeq2 models of ADNC stage and sex × stage on the cohort pseudobulks, bound to the primary export's feature set, with the predeclared concordance rule against the primary MOTCO analysis and the committed report.

### Modified Capabilities
<!-- None. `motco-trajectory-analysis`, `motco-input-export` and `astrocyte-donor-pseudobulk` are read, checked and never written. -->

## Impact

- `src/sead_sexdiff/data.py`: constants `STAGE_DE_MODELS`, `STAGE_DE_COVARIATES`, `STAGE_DE_ALPHA`, `STAGE_DE_REFERENCE_LEVELS`, `STAGE_DE_FEATURE_STATUSES`, `STAGE_DE_CONCORDANCE_RULE`, `STAGE_DE_*` output names; `run_stage_de`, a per-fit helper and `render_stage_de`. The reference checks reuse `_check_reference` and the export's hash checks. The legacy PFC `differential_expression` is untouched.
- `src/sead_sexdiff/cli.py`: new `stage-de` subcommand (`--cohort`, `--omic`, `--pseudobulk-dir`, `--input-dir`, `--analysis-dir`, `--out-dir`, `--n-cpus`, `--force`). The models, covariates and α are constants, not flags.
- Dependencies: `pydeseq2 >= 0.5.3` is already pinned; no new dependency.
- `reports/`: one new committed report. `results/` stays gitignored.
- Compute: four fits per cohort × omic; minutes for ~20,000 genes, longer for ~125,000 peaks (pydeseq2 is multi-core through `n_cpus`; the ATAC fits are the only slow part). Runs on the workstation holding `results/pseudobulk/`, `results/motco_input/` and `results/motco_analysis/`; this machine has all three.
- Scientific choices introduced here (flag and ask before changing): the four models and their tested coefficients, the factor model as the primary interaction model, the reference levels, the covariate set per omic, the export-bound feature set and the two reporting sets, α and the BH / independent-filtering defaults, the concordance rule, and the rule that the primary analysis stays primary.

# Pseudobulk stage DE of the MOTCO astrocyte blocks with limma-voom, beside DESeq2 — 2026-10-22

> A second predeclared method for the per-feature stage and sex × stage models of `mtg-astrocyte-stage-de` (OpenSpec change `mtg-astrocyte-stage-de-limma`, spec `pseudobulk-stage-de-limma`). The change was proposed **after** the DESeq2 results were known (`reports/motco-astrocyte-stage-de-2026-10-21.md`: no interaction at padj ≤ 0.05 in any cohort or omic; 1–2 RNA genes and 22–36 ATAC peaks at stage 2 vs 0; every concordance cell concordant). It re-tests nothing and looks for no method that finds more: it fills four stated gaps of the DESeq2 cross-check — no joint interaction test, an approximate FDR, no outlier protection, strongly varying library quality — with limma-voom on the export's own TMM normalization, per-donor quality weights and a joint moderated F, all fixed before any limma run. The same donors, raw counts, design matrices, covariates and abundance-passing features as the DESeq2 step enter the fits, through that step's own preparation code. Every number below is pasted by script from the generated files (`stage_de_limma.md`, `stage_de_limma_summary.csv`, `method_comparison.csv`, `stage_de_limma_manifest.json`, the fit records and tables of both methods, the primary `result.json`, the run logs) or, where named, from the verification outputs of the change (synthetic data and simulations), and was checked against those files before commit. The hand-written text carries only the reading, under the rules below.

**Reading rules (binding).** The DESeq2 results stay the reported ones, and no later change may promote whichever method finds more. A call that differs between the methods is method dependence, not a result of either; the same donors and counts enter both, so agreement is not replication. A zero count is "none detected at this n". The union counts carry no joint error-rate claim, and the joint F's claim is approximate (its exact-F diagnostic stands beside it). The concordance cells cross-check the primary analysis and are not a second test of its hypothesis; the primary analysis stays primary with its non-rejection rules (`angle` and `shape` non-rejections uninformative; a `delta` non-rejection not evidence that the trajectories agree). The step read sex and stage, so it is a check beside the analysis and never a filter. `cohort_30` is primary and `cohort_50` sensitivity (a donor subset of `cohort_30`, not an independent replication).

## 1. Provenance

- Software: R 4.6.1, limma 3.68.4, edgeR 4.10.1, jsonlite 2.0.0 (the Nix shell); `src/sead_sexdiff/stage_de_limma.R` sha256 `0543a17dc2d12ba77ef00d27b53b08a68095b8b1ce00c7b98dc63e31eb5d70a8`; Python 3.12.15, numpy 2.3.5, pandas 2.3.3, scipy 1.15.3, formulaic 1.2.1, formulaic-contrasts 1.0.0; MOTCO 0.9.0 installed for the reference check. sead-sexdiff revision `695f84f304bef3ae7e602b17e72dafcd75ce6450` (dirty: True; the implementation of this change was uncommitted at run time; the rerun of section 4 reproduced the outputs byte for byte).
- Machine: `nixos-home`, 2026-10-22. The change's tasks named nixos-dell-laptop, which was offline. This machine held `results/pseudobulk/` and `results/motco_input/` hashing as on the cluster (18 files, `sha256sum -c` on the cluster); its `results/motco_analysis/` was the MOTCO 0.8.0 analysis, equal file for file to the cluster's `results/motco_analysis_0.8.0/`, and was moved there (the cluster's layout). The MOTCO 0.9.0 primary (`results/motco_analysis/`) and the DESeq2 outputs (`results/stage_de/`) were pulled from the cluster (`/home1/tgonza/sead_sexdiff/results/`), and all 60 pulled files verify against sums computed there.
- Run: `uv run sead-sexdiff stage-de-limma --cohort all --omic both`; 13 min 16 s wall time (fits 11–14 s per RNA model, 69–92 s per ATAC model, the exact-F diagnostic included).
- Pseudobulk inputs: `rna_counts.h5ad` sha256 `1253bed9e5f4a762393ecfbb74ae00c550c38cecd648c8a98ea817170885de40`, `atac_counts.h5ad` sha256 `2cfbe937d8e94cb1f3e0024ba0dcc1df0fd731d658b55be62f4e867fb46cff0d`, `donor_table.csv` sha256 `1b34658b2d65132472a6af2d7be7a26f239df305e1f454f289983c89f0a45e71`, `manifest.json` sha256 `2e93cc0b526856fc043f908a424a11266187146636bc75afed7edd187989a25c`.
- cohort_30: 80 donors, cells Female/stage 0 11, Female/stage 1 10, Female/stage 2 28, Male/stage 0 9, Male/stage 1 10, Male/stage 2 12, k = 9; primary export manifest sha256 `18f2289e0d4c84ed18f4a437ae569052544e9d150ed2726385273694b4a3b0b6`; primary analysis `result.json` sha256 `1d6154b02d7704ef9313ce752035a124cf6cb3fbd5350f0428439a542d5190cd` (CV-selected rank 1, MOTCO 0.9.0, seed 20261019, 999 permutations; `delta` rejected: no); DESeq2 outputs `stage_de_manifest.json` sha256 `b1015b151db2684d495796a463195b5ed36979fe99feae16df06aabce25b9b55`, `stage_de_summary.csv` sha256 `6b4947e8db81edacdf48ac163af8352d747731b1dea090ed5cec5dd50595e903` (pydeseq2 0.5.3; 23 files).
- cohort_50: 74 donors, cells Female/stage 0 10, Female/stage 1 9, Female/stage 2 25, Male/stage 0 9, Male/stage 1 9, Male/stage 2 12, k = 9; primary export manifest sha256 `dd14a6e7d5583fb527e8a99c7db43d12cf237f5451522e8875c5a0a3401c5cb3`; primary analysis `result.json` sha256 `4c6e7e58e39d8ce82c9b93907b70b2110fc2c415e32da40fef8668b74faff9fb` (CV-selected rank 1, MOTCO 0.9.0, seed 20261019, 999 permutations; `delta` rejected: no); DESeq2 outputs `stage_de_manifest.json` sha256 `f38a00d21dfdb7c29b6258b19394c61de487c82a5acea35fe83b1c4977903e12`, `stage_de_summary.csv` sha256 `9f33acf5b59e15bdc3f267c0396e1f1c3f2d9da0f191090ae25a52103c7de4ef` (pydeseq2 0.5.3; 23 files).
- Checks before any fit, per cohort: every check of `stage-de` through its code path (inputs, the pseudobulk manifest hash, the reproduction of the export from the counts, donors, sex, stage, `k`, cells, statuses, the primary analysis against the export); the DESeq2 outputs complete for both omics, made from the inputs now on disk (every recorded input sha256), every table hashing as recorded and the summary and concordance cells equal to a recomputation from the tables; and in R, per omic, the donor order, edgeR's TMM against the export manifest's factors and every design rebuilt by `model.matrix`:
  - cohort_30 rna: 20387 features tested, 5000 exported; edgeR TMM equals the export's to 1.11e-16 (absolute; tolerance 1e-08); the export reproduced from the counts (TMM factors to 0, tested-feature variances to 4.95e-10, exported values to 4.99e-09); designs rebuilt in R equal the formulaic ones in 4 of 4 models (max abs diff 0).
  - cohort_30 atac: 128916 features tested, 5000 exported; edgeR TMM equals the export's to 2.22e-16 (absolute; tolerance 1e-08); the export reproduced from the counts (TMM factors to 0, tested-feature variances to 4.98e-10, exported values to 5e-10); designs rebuilt in R equal the formulaic ones in 4 of 4 models (max abs diff 0).
  - cohort_50 rna: 20238 features tested, 5000 exported; edgeR TMM equals the export's to 1.11e-16 (absolute; tolerance 1e-08); the export reproduced from the counts (TMM factors to 0, tested-feature variances to 4.96e-10, exported values to 4.99e-09); designs rebuilt in R equal the formulaic ones in 4 of 4 models (max abs diff 0).
  - cohort_50 atac: 124818 features tested, 5000 exported; edgeR TMM equals the export's to 1.11e-16 (absolute; tolerance 1e-08); the export reproduced from the counts (TMM factors to 0, tested-feature variances to 4.99e-10, exported values to 5e-10); designs rebuilt in R equal the formulaic ones in 4 of 4 models (max abs diff 0).
  Every one of the 21 predeclared settings equals the installed default in 16 of 16 fits. Every input and every DESeq2 file hashed identically before and after the run (78 files, `sha256sum -c`).
- Method (`STAGE_DE_LIMMA_METHOD`, predeclared): edgeR TMM on the tested counts (library size = their column sums; trims 0.3 / 0.05) → `voomWithQualityWeights` (voom's precision weights times one quality weight per donor; span 0.5, adaptive; `genebygene`; `prior.n = 10`) → `lmFit` → `eBayes` (no trend, not robust, legacy hyperparameters since every feature has the same residual df). A moderated t per tested coefficient; a joint moderated F per factor model (`stage_joint`, `interaction_joint`; df1 = 2); the union of a factor model's two moderated t's for the like-for-like comparison only. Significance: Benjamini–Hochberg adjusted p ≤ 0.05 over all tested features (no independent filter, so limma's BH families are the whole abundance-passing set while DESeq2's exclude the features its baseMean filter removes). voom's log-CPM, `log2((count + 0.5) / (lib.size × factor + 1) × 1e6)`, is not the export's prior-2 log2-CPM; the normalization is the export's TMM, DESeq2's was median of ratios.
- Models, covariates, reference levels and feature sets are those of `stage-de` (`STAGE_DE_MODELS`, `STAGE_DE_COVARIATES`, `Female` / stage `0`, `features_<omic>.csv` status `variance` or `kept`; `exported` = `kept`).

## 2. Primary cohort: cohort_30 (≥ 30 astrocyte nuclei in both omics, n = 80)

### Concordance (pasted from `results/stage_de_limma/cohort_30/stage_de_limma.md`)

Rule (`STAGE_DE_LIMMA_CONCORDANCE_RULE`, predeclared): `interaction_detected` = at least one exported feature (status kept) with adj.P.Val <= STAGE_DE_ALPHA on the joint moderated F of factor_interaction's two interaction coefficients (interaction_joint); `delta_rejected` = the primary analysis rejects delta (p <= ANALYSIS_ALPHA, ANALYSIS_REJECTION_RULE; result.json); concordant when both are true or both are false; discordant otherwise, written as a finding. Reported beside the DESeq2 cell (STAGE_DE_CONCORDANCE_RULE, the union of two Wald tests), which stays the reported cell of stage-de; methods_agree compares the two interaction_detected values and a disagreement is a method-dependence finding; neither cell replaces the other or changes a MOTCO reading. The primary analysis (`results/motco_analysis/cohort_30`, CV-selected rank 1, MOTCO 0.9.0) is read, never rerun or modified.

| omic | DESeq2: exported features (union of two Wald tests) | DESeq2 cell | limma: exported features (interaction_joint) | limma cell | delta rejected | methods agree |
|---|---:|---|---:|---|:---:|:---:|
| rna | 0 | concordant | 0 | concordant | no | yes |
| atac | 0 | concordant | 0 | concordant | no | yes |

### Counts: rna (pasted from `results/stage_de_limma/cohort_30/stage_de_limma.md`)

5000 exported and 20387 abundance-passing features tested; significant / tested at adjusted p ≤ α. The DESeq2 columns are the same coefficient's Wald test, the union for a union row, and DESeq2's union of the model's two Wald tests for a joint-F row (a different test, `joint_f_vs_union`).

| model | test | kind | design column | limma exported | limma abundance | DESeq2 exported | DESeq2 abundance | min adj. p limma (exported / abundance) |
|---|---|---|---|---:|---:|---:|---:|---:|
| factor_additive | stage_1_vs_0 | moderated t | `stage[T.1]` | 0 / 5000 | 0 / 20387 | 0 / 5000 | 0 / 20387 | 0.475 / 0.223 |
| factor_additive | stage_2_vs_0 | moderated t | `stage[T.2]` | 1 / 5000 | 1 / 20387 | 1 / 5000 | 1 / 20387 | 0.000463 / 0.000463 |
| factor_additive | stage_joint | joint F | `stage[T.1];stage[T.2]` | 1 / 5000 | 1 / 20387 | 1 / 5000 | 1 / 20387 | 0.00107 / 0.00107 |
| factor_additive | union | union (no joint error rate) | `stage[T.1];stage[T.2]` | 1 / 5000 | 1 / 20387 | 1 / 5000 | 1 / 20387 | 0.000463 / 0.000463 |
| factor_interaction | male_x_stage_1 | moderated t | `sex[T.Male]:stage[T.1]` | 0 / 5000 | 0 / 20387 | 0 / 5000 | 0 / 20387 | 0.992 / 0.718 |
| factor_interaction | male_x_stage_2 | moderated t | `sex[T.Male]:stage[T.2]` | 0 / 5000 | 0 / 20387 | 0 / 5000 | 0 / 20387 | 0.793 / 0.793 |
| factor_interaction | interaction_joint | joint F | `sex[T.Male]:stage[T.1];sex[T.Male]:stage[T.2]` | 0 / 5000 | 0 / 20387 | 0 / 5000 | 0 / 20387 | 0.898 / 0.55 |
| factor_interaction | union | union (no joint error rate) | `sex[T.Male]:stage[T.1];sex[T.Male]:stage[T.2]` | 0 / 5000 | 0 / 20387 | 0 / 5000 | 0 / 20387 | 0.793 / 0.718 |
| trend_additive | stage_num | moderated t | `stage_num` | 1 / 5000 | 1 / 20387 | 2 / 5000 | 2 / 20387 | 0.000204 / 0.000204 |
| trend_interaction | male_x_stage_num | moderated t | `sex[T.Male]:stage_num` | 0 / 5000 | 0 / 20387 | 0 / 5000 | 0 / 20387 | 0.763 / 0.763 |

### Counts: atac (pasted from `results/stage_de_limma/cohort_30/stage_de_limma.md`)

5000 exported and 128916 abundance-passing features tested; significant / tested at adjusted p ≤ α. The DESeq2 columns are the same coefficient's Wald test, the union for a union row, and DESeq2's union of the model's two Wald tests for a joint-F row (a different test, `joint_f_vs_union`).

| model | test | kind | design column | limma exported | limma abundance | DESeq2 exported | DESeq2 abundance | min adj. p limma (exported / abundance) |
|---|---|---|---|---:|---:|---:|---:|---:|
| factor_additive | stage_1_vs_0 | moderated t | `stage[T.1]` | 0 / 5000 | 0 / 128916 | 0 / 5000 | 0 / 128916 | 0.763 / 0.0973 |
| factor_additive | stage_2_vs_0 | moderated t | `stage[T.2]` | 0 / 5000 | 9 / 128916 | 0 / 5000 | 22 / 128916 | 0.0929 / 0.0457 |
| factor_additive | stage_joint | joint F | `stage[T.1];stage[T.2]` | 0 / 5000 | 0 / 128916 | 0 / 5000 | 22 / 128916 | 0.27 / 0.108 |
| factor_additive | union | union (no joint error rate) | `stage[T.1];stage[T.2]` | 0 / 5000 | 9 / 128916 | 0 / 5000 | 22 / 128916 | 0.0929 / 0.0457 |
| factor_interaction | male_x_stage_1 | moderated t | `sex[T.Male]:stage[T.1]` | 0 / 5000 | 0 / 128916 | 0 / 5000 | 0 / 128916 | 0.731 / 0.0934 |
| factor_interaction | male_x_stage_2 | moderated t | `sex[T.Male]:stage[T.2]` | 0 / 5000 | 0 / 128916 | 0 / 5000 | 0 / 128916 | 0.963 / 0.963 |
| factor_interaction | interaction_joint | joint F | `sex[T.Male]:stage[T.1];sex[T.Male]:stage[T.2]` | 0 / 5000 | 0 / 128916 | 0 / 5000 | 0 / 128916 | 1 / 0.4 |
| factor_interaction | union | union (no joint error rate) | `sex[T.Male]:stage[T.1];sex[T.Male]:stage[T.2]` | 0 / 5000 | 0 / 128916 | 0 / 5000 | 0 / 128916 | 0.731 / 0.0934 |
| trend_additive | stage_num | moderated t | `stage_num` | 0 / 5000 | 0 / 128916 | 0 / 5000 | 23 / 128916 | 0.0977 / 0.0545 |
| trend_interaction | male_x_stage_num | moderated t | `sex[T.Male]:stage_num` | 0 / 5000 | 0 / 128916 | 0 / 5000 | 0 / 128916 | 1 / 1 |

### Method comparison: rna (pasted from `results/stage_de_limma/cohort_30/stage_de_limma.md`)

Features called by both methods and by one only, whether each method detects anything, and, per coefficient, the Spearman correlations of DESeq2's log2 fold change with limma's logFC and of the Wald statistic with the moderated t (over the features with both values finite). Descriptive: no p-value, no gate; a row whose detection differs is a method-dependence finding.

| model | test | kind | set | DESeq2 | limma | both | only DESeq2 | only limma | detection agrees | ρ log2FC | ρ statistic | compared |
|---|---|---|---|---:|---:|---:|---:|---:|:---:|---:|---:|---:|
| factor_additive | stage_1_vs_0 | coefficient | exported | 0 | 0 | 0 | 0 | 0 | yes | 0.758 | 0.780 | 5000 |
| factor_additive | stage_1_vs_0 | coefficient | abundance | 0 | 0 | 0 | 0 | 0 | yes | 0.808 | 0.849 | 20387 |
| factor_additive | stage_2_vs_0 | coefficient | exported | 1 | 1 | 1 | 0 | 0 | yes | 0.870 | 0.879 | 5000 |
| factor_additive | stage_2_vs_0 | coefficient | abundance | 1 | 1 | 1 | 0 | 0 | yes | 0.915 | 0.931 | 20387 |
| factor_additive | union | union | exported | 1 | 1 | 1 | 0 | 0 | yes | n/a | n/a |  |
| factor_additive | union | union | abundance | 1 | 1 | 1 | 0 | 0 | yes | n/a | n/a |  |
| factor_additive | stage_joint | joint_f_vs_union | exported | 1 | 1 | 1 | 0 | 0 | yes | n/a | n/a |  |
| factor_additive | stage_joint | joint_f_vs_union | abundance | 1 | 1 | 1 | 0 | 0 | yes | n/a | n/a |  |
| factor_interaction | male_x_stage_1 | coefficient | exported | 0 | 0 | 0 | 0 | 0 | yes | 0.812 | 0.825 | 5000 |
| factor_interaction | male_x_stage_1 | coefficient | abundance | 0 | 0 | 0 | 0 | 0 | yes | 0.838 | 0.864 | 20387 |
| factor_interaction | male_x_stage_2 | coefficient | exported | 0 | 0 | 0 | 0 | 0 | yes | 0.860 | 0.868 | 5000 |
| factor_interaction | male_x_stage_2 | coefficient | abundance | 0 | 0 | 0 | 0 | 0 | yes | 0.887 | 0.907 | 20387 |
| factor_interaction | union | union | exported | 0 | 0 | 0 | 0 | 0 | yes | n/a | n/a |  |
| factor_interaction | union | union | abundance | 0 | 0 | 0 | 0 | 0 | yes | n/a | n/a |  |
| factor_interaction | interaction_joint | joint_f_vs_union | exported | 0 | 0 | 0 | 0 | 0 | yes | n/a | n/a |  |
| factor_interaction | interaction_joint | joint_f_vs_union | abundance | 0 | 0 | 0 | 0 | 0 | yes | n/a | n/a |  |
| trend_additive | stage_num | coefficient | exported | 2 | 1 | 1 | 1 | 0 | yes | 0.873 | 0.880 | 5000 |
| trend_additive | stage_num | coefficient | abundance | 2 | 1 | 1 | 1 | 0 | yes | 0.915 | 0.931 | 20387 |
| trend_interaction | male_x_stage_num | coefficient | exported | 0 | 0 | 0 | 0 | 0 | yes | 0.859 | 0.868 | 5000 |
| trend_interaction | male_x_stage_num | coefficient | abundance | 0 | 0 | 0 | 0 | 0 | yes | 0.886 | 0.909 | 20387 |

### Method comparison: atac (pasted from `results/stage_de_limma/cohort_30/stage_de_limma.md`)

Features called by both methods and by one only, whether each method detects anything, and, per coefficient, the Spearman correlations of DESeq2's log2 fold change with limma's logFC and of the Wald statistic with the moderated t (over the features with both values finite). Descriptive: no p-value, no gate; a row whose detection differs is a method-dependence finding.

| model | test | kind | set | DESeq2 | limma | both | only DESeq2 | only limma | detection agrees | ρ log2FC | ρ statistic | compared |
|---|---|---|---|---:|---:|---:|---:|---:|:---:|---:|---:|---:|
| factor_additive | stage_1_vs_0 | coefficient | exported | 0 | 0 | 0 | 0 | 0 | yes | 0.789 | 0.794 | 5000 |
| factor_additive | stage_1_vs_0 | coefficient | abundance | 0 | 0 | 0 | 0 | 0 | yes | 0.862 | 0.869 | 128916 |
| factor_additive | stage_2_vs_0 | coefficient | exported | 0 | 0 | 0 | 0 | 0 | yes | 0.826 | 0.834 | 5000 |
| factor_additive | stage_2_vs_0 | coefficient | abundance | 22 | 9 | 3 | 19 | 6 | yes | 0.886 | 0.897 | 128916 |
| factor_additive | union | union | exported | 0 | 0 | 0 | 0 | 0 | yes | n/a | n/a |  |
| factor_additive | union | union | abundance | 22 | 9 | 3 | 19 | 6 | yes | n/a | n/a |  |
| factor_additive | stage_joint | joint_f_vs_union | exported | 0 | 0 | 0 | 0 | 0 | yes | n/a | n/a |  |
| factor_additive | stage_joint | joint_f_vs_union | abundance | 22 | 0 | 0 | 22 | 0 | no | n/a | n/a |  |
| factor_interaction | male_x_stage_1 | coefficient | exported | 0 | 0 | 0 | 0 | 0 | yes | 0.777 | 0.788 | 5000 |
| factor_interaction | male_x_stage_1 | coefficient | abundance | 0 | 0 | 0 | 0 | 0 | yes | 0.852 | 0.875 | 128916 |
| factor_interaction | male_x_stage_2 | coefficient | exported | 0 | 0 | 0 | 0 | 0 | yes | 0.808 | 0.818 | 5000 |
| factor_interaction | male_x_stage_2 | coefficient | abundance | 0 | 0 | 0 | 0 | 0 | yes | 0.860 | 0.885 | 128916 |
| factor_interaction | union | union | exported | 0 | 0 | 0 | 0 | 0 | yes | n/a | n/a |  |
| factor_interaction | union | union | abundance | 0 | 0 | 0 | 0 | 0 | yes | n/a | n/a |  |
| factor_interaction | interaction_joint | joint_f_vs_union | exported | 0 | 0 | 0 | 0 | 0 | yes | n/a | n/a |  |
| factor_interaction | interaction_joint | joint_f_vs_union | abundance | 0 | 0 | 0 | 0 | 0 | yes | n/a | n/a |  |
| trend_additive | stage_num | coefficient | exported | 0 | 0 | 0 | 0 | 0 | yes | 0.830 | 0.838 | 5000 |
| trend_additive | stage_num | coefficient | abundance | 23 | 0 | 0 | 23 | 0 | no | 0.886 | 0.898 | 128916 |
| trend_interaction | male_x_stage_num | coefficient | exported | 0 | 0 | 0 | 0 | 0 | yes | 0.800 | 0.811 | 5000 |
| trend_interaction | male_x_stage_num | coefficient | abundance | 0 | 0 | 0 | 0 | 0 | yes | 0.856 | 0.882 | 128916 |

### Features at adjusted p ≤ 0.05 (generated from the per-test tables; descriptive)

Every feature at adjusted p ≤ 0.05 per limma test, by adjusted p (ties by p-value, then name); a test without one lists its five smallest adjusted p. DESeq2's padj for the same feature and coefficient stands beside it (for a joint F, DESeq2's smaller padj of the model's two coefficients; empty = removed by its independent filter). `exported` = one of MOTCO's 5,000 features.

- rna `factor_additive` `stage_1_vs_0`: none at adjusted p ≤ 0.05; smallest five:

  | feature | AveExpr | logFC | SE | t | adj. p | DESeq2 padj | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `TMEM176A` | 2.00 | 0.665 | 0.149 | 4.48 | 0.223 | 0.613 | no |
  | `GINS2` | 0.85 | -0.614 | 0.138 | -4.44 | 0.223 | 0.997 | no |
  | `STAC3` | 1.41 | -0.732 | 0.166 | -4.42 | 0.223 | 0.613 | no |
  | `SELENOV` | -0.57 | -1.191 | 0.289 | -4.13 | 0.475 | 0.613 | yes |
  | `ICOSLG` | -0.49 | 1.173 | 0.299 | 3.92 | 0.776 | 0.997 | yes |

- rna `factor_additive` `stage_2_vs_0`: 1 at adjusted p ≤ 0.05 (1 exported):

  | feature | AveExpr | logFC | SE | t | adj. p | DESeq2 padj | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `CREB5` | 4.91 | 1.143 | 0.183 | 6.24 | 0.000463 | 7.17e-06 | yes |

- rna `factor_additive` `stage_joint`: 1 at adjusted p ≤ 0.05 (1 exported):

  | feature | AveExpr | logFC_stage_1_vs_0 | logFC_stage_2_vs_0 | F | adj. p | DESeq2 padj (min) | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `CREB5` | 4.91 | 0.349 | 1.143 | 21.07 | 0.00107 | 7.17e-06 | yes |

- rna `factor_interaction` `male_x_stage_1`: none at adjusted p ≤ 0.05; smallest five:

  | feature | AveExpr | logFC | SE | t | adj. p | DESeq2 padj | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `SLC25A24` | 3.82 | 1.671 | 0.379 | 4.40 | 0.718 | 1 | no |
  | `ORAI1` | 0.82 | -1.309 | 0.327 | -4.01 | 0.992 | 1 | no |
  | `MUL1` | 2.87 | -0.741 | 0.194 | -3.82 | 0.992 | 1 | no |
  | `BMP8A` | 0.48 | -1.387 | 0.365 | -3.80 | 0.992 | 1 | no |
  | `BNIP3` | 5.10 | 1.086 | 0.295 | 3.68 | 0.992 | 1 | no |

- rna `factor_interaction` `male_x_stage_2`: none at adjusted p ≤ 0.05; smallest five:

  | feature | AveExpr | logFC | SE | t | adj. p | DESeq2 padj | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `AL109976.1` | -0.41 | -1.925 | 0.463 | -4.16 | 0.793 | 0.684 | yes |
  | `AC133785.1` | 0.98 | 1.523 | 0.375 | 4.06 | 0.793 | 0.135 | yes |
  | `AC091074.2` | -0.09 | 1.887 | 0.472 | 4.00 | 0.793 | 0.461 | yes |
  | `DNMT1` | 4.26 | -0.572 | 0.146 | -3.91 | 0.793 | 0.684 | no |
  | `NAA11` | 1.60 | 2.702 | 0.691 | 3.91 | 0.793 | 0.68 | yes |

- rna `factor_interaction` `interaction_joint`: none at adjusted p ≤ 0.05; smallest five:

  | feature | AveExpr | logFC_male_x_stage_1 | logFC_male_x_stage_2 | F | adj. p | DESeq2 padj (min) | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `SLC25A24` | 3.82 | 1.671 | 0.259 | 12.17 | 0.55 | 0.932 | no |
  | `AC091074.2` | -0.09 | 0.442 | 1.887 | 9.59 | 0.898 | 0.461 | yes |
  | `BMP8A` | 0.48 | -1.387 | -0.203 | 9.21 | 0.898 | 0.928 | no |
  | `ORAI1` | 0.82 | -1.309 | -1.055 | 9.02 | 0.898 | 0.684 | no |
  | `AL109976.1` | -0.41 | -1.555 | -1.925 | 8.89 | 0.898 | 0.684 | yes |

- rna `trend_additive` `stage_num`: 1 at adjusted p ≤ 0.05 (1 exported):

  | feature | AveExpr | logFC | SE | t | adj. p | DESeq2 padj | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `CREB5` | 4.91 | 0.589 | 0.092 | 6.43 | 0.000204 | 4.34e-06 | yes |

- rna `trend_interaction` `male_x_stage_num`: none at adjusted p ≤ 0.05; smallest five:

  | feature | AveExpr | logFC | SE | t | adj. p | DESeq2 padj | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `AC091074.2` | -0.09 | 0.968 | 0.233 | 4.15 | 0.763 | 0.349 | yes |
  | `AL109976.1` | -0.41 | -0.942 | 0.231 | -4.08 | 0.763 | 0.703 | yes |
  | `AC133785.1` | 0.98 | 0.743 | 0.182 | 4.07 | 0.763 | 0.108 | yes |
  | `AL160408.1` | -0.64 | 1.084 | 0.275 | 3.94 | 0.819 | 0.712 | yes |
  | `NAA11` | 1.60 | 1.308 | 0.339 | 3.86 | 0.819 | 0.569 | yes |

- atac `factor_additive` `stage_1_vs_0`: none at adjusted p ≤ 0.05; smallest five:

  | feature | AveExpr | logFC | SE | t | adj. p | DESeq2 padj | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `chr5:10592071-10592905` | -1.34 | 1.900 | 0.361 | 5.27 | 0.0973 | 0.283 | no |
  | `chr12:68291447-68292819` | -0.87 | 1.932 | 0.421 | 4.58 | 0.364 | 0.865 | no |
  | `chr5:107417659-107418312` | -1.18 | 1.675 | 0.371 | 4.51 | 0.364 | 0.664 | no |
  | `chr8:119152756-119153569` | -0.54 | 1.671 | 0.372 | 4.49 | 0.364 | 0.602 | no |
  | `chr1:245589383-245590226` | 2.21 | -0.802 | 0.179 | -4.49 | 0.364 | 0.861 | no |

- atac `factor_additive` `stage_2_vs_0`: 9 at adjusted p ≤ 0.05 (0 exported):

  | feature | AveExpr | logFC | SE | t | adj. p | DESeq2 padj | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `chr3:194880040-194881125` | 4.03 | -0.576 | 0.109 | -5.29 | 0.0457 | 0.0448 | no |
  | `chr12:108750482-108752712` | 2.89 | -0.809 | 0.154 | -5.26 | 0.0457 | 0.124 | no |
  | `chr12:106328829-106330471` | 3.73 | -0.502 | 0.098 | -5.15 | 0.0457 | 0.0448 | no |
  | `chr1:231053444-231054140` | 2.09 | -0.792 | 0.156 | -5.06 | 0.0457 | 0.0504 | no |
  | `chr1:155481563-155482044` | 2.05 | -1.063 | 0.212 | -5.03 | 0.0457 | 0.123 | no |
  | `chr3:143134550-143135514` | 2.80 | -0.759 | 0.153 | -4.97 | 0.0457 | 0.0489 | no |
  | `chr4:123554731-123556429` | 3.74 | -0.569 | 0.115 | -4.96 | 0.0457 | 0.0936 | no |
  | `chr19:15279962-15281620` | 3.49 | -0.543 | 0.110 | -4.95 | 0.0457 | 0.126 | no |
  | `chr6:11605152-11608407` | 3.25 | -0.669 | 0.136 | -4.93 | 0.0457 | 0.109 | no |

- atac `factor_additive` `stage_joint`: none at adjusted p ≤ 0.05; smallest five:

  | feature | AveExpr | logFC_stage_1_vs_0 | logFC_stage_2_vs_0 | F | adj. p | DESeq2 padj (min) | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `chr5:10592071-10592905` | -1.34 | 1.900 | 1.419 | 15.25 | 0.108 | 0.283 | no |
  | `chr12:108750482-108752712` | 2.89 | -0.632 | -0.809 | 14.68 | 0.108 | 0.124 | no |
  | `chr3:143134550-143135514` | 2.80 | -0.687 | -0.759 | 14.28 | 0.108 | 0.0489 | no |
  | `chr3:194880040-194881125` | 4.03 | -0.397 | -0.576 | 14.27 | 0.108 | 0.0448 | no |
  | `chr19:15279962-15281620` | 3.49 | -0.479 | -0.543 | 13.70 | 0.112 | 0.126 | no |

- atac `factor_interaction` `male_x_stage_1`: none at adjusted p ≤ 0.05; smallest five:

  | feature | AveExpr | logFC | SE | t | adj. p | DESeq2 padj | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `chr17:304482-305830` | 1.99 | -2.046 | 0.387 | -5.28 | 0.0934 | 1 | no |
  | `chr5:172556573-172560033` | 3.01 | -2.176 | 0.469 | -4.64 | 0.642 | 1 | no |
  | `chr12:48753307-48754450` | 1.42 | 1.766 | 0.388 | 4.55 | 0.642 | 0.938 | no |
  | `chr7:90932085-90932904` | -0.64 | -2.883 | 0.650 | -4.44 | 0.681 | 1 | no |
  | `chr19:12513516-12514280` | 2.02 | 2.298 | 0.522 | 4.40 | 0.681 | 1 | no |

- atac `factor_interaction` `male_x_stage_2`: none at adjusted p ≤ 0.05; smallest five:

  | feature | AveExpr | logFC | SE | t | adj. p | DESeq2 padj | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `chr17:46259862-46261495` | 1.85 | 1.786 | 0.401 | 4.45 | 0.963 | 0.638 | no |
  | `chr5:179318308-179318574` | 0.31 | 2.207 | 0.511 | 4.32 | 0.963 | 0.704 | no |
  | `chr19:8585367-8586347` | 0.79 | 2.441 | 0.577 | 4.23 | 0.963 | 0.835 | no |
  | `chr10:29589886-29591981` | 1.21 | -1.896 | 0.459 | -4.13 | 0.963 | 0.517 | no |
  | `chr18:44699213-44701847` | 0.34 | -2.168 | 0.525 | -4.13 | 0.963 | 0.801 | no |

- atac `factor_interaction` `interaction_joint`: none at adjusted p ≤ 0.05; smallest five:

  | feature | AveExpr | logFC_male_x_stage_1 | logFC_male_x_stage_2 | F | adj. p | DESeq2 padj (min) | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `chr17:304482-305830` | 1.99 | -2.046 | -0.825 | 14.40 | 0.4 | 0.891 | no |
  | `chr5:172556573-172560033` | 3.01 | -2.176 | -1.753 | 11.89 | 1 | 0.875 | no |
  | `chr17:46259862-46261495` | 1.85 | 1.850 | 1.786 | 11.61 | 1 | 0.638 | no |
  | `chr12:48753307-48754450` | 1.42 | 1.766 | 1.213 | 10.86 | 1 | 0.801 | no |
  | `chr18:44699213-44701847` | 0.34 | -2.318 | -2.168 | 10.46 | 1 | 0.801 | no |

- atac `trend_additive` `stage_num`: none at adjusted p ≤ 0.05; smallest five:

  | feature | AveExpr | logFC | SE | t | adj. p | DESeq2 padj | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `chr3:194880040-194881125` | 4.03 | -0.283 | 0.055 | -5.20 | 0.0545 | 0.0498 | no |
  | `chr1:155481563-155482044` | 2.05 | -0.536 | 0.105 | -5.13 | 0.0545 | 0.109 | no |
  | `chr12:106328829-106330471` | 3.73 | -0.250 | 0.049 | -5.12 | 0.0545 | 0.0498 | no |
  | `chr12:108750482-108752712` | 2.89 | -0.397 | 0.078 | -5.06 | 0.0545 | 0.145 | no |
  | `chr1:231053444-231054140` | 2.09 | -0.393 | 0.078 | -5.02 | 0.0545 | 0.0571 | no |

- atac `trend_interaction` `male_x_stage_num`: none at adjusted p ≤ 0.05; smallest five:

  | feature | AveExpr | logFC | SE | t | adj. p | DESeq2 padj | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `chr1:33725071-33727945` | -0.72 | 1.512 | 0.344 | 4.39 | 1 | 0.845 | no |
  | `chr6:88007915-88009365` | -1.62 | -1.619 | 0.370 | -4.37 | 1 | 0.845 | no |
  | `chr5:179318308-179318574` | 0.31 | 1.022 | 0.240 | 4.26 | 1 | 0.845 | no |
  | `chr7:72077958-72078469` | -0.94 | 1.333 | 0.329 | 4.05 | 1 | 0.922 | no |
  | `chr9:112913366-112914928` | -1.51 | 1.410 | 0.350 | 4.03 | 1 | 0.953 | no |

### Sample weights (pasted from `results/stage_de_limma/cohort_30/stage_de_limma.md`)

One quality weight per donor and model (arrayWeights inside voomWithQualityWeights, shrunk toward 1 with `prior.n = 10`): minimum, median, maximum and the mean per sex × stage cell. Shown without a gate; a cell with a low mean weight counts less in that model's fit.

| omic | model | min | median | max | Female/stage 0 | Female/stage 1 | Female/stage 2 | Male/stage 0 | Male/stage 1 | Male/stage 2 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| rna | factor_additive | 0.366 | 0.971 | 2.222 | 1.176 | 1.243 | 0.981 | 1.079 | 1.148 | 1.104 |
| rna | factor_interaction | 0.364 | 0.998 | 2.389 | 1.171 | 1.261 | 0.999 | 1.072 | 1.123 | 1.081 |
| rna | trend_additive | 0.369 | 0.989 | 2.072 | 1.177 | 1.196 | 0.979 | 1.088 | 1.150 | 1.104 |
| rna | trend_interaction | 0.375 | 1.013 | 2.082 | 1.176 | 1.194 | 0.987 | 1.099 | 1.143 | 1.090 |
| atac | factor_additive | 0.441 | 1.019 | 1.459 | 1.009 | 1.020 | 1.021 | 0.935 | 1.084 | 1.048 |
| atac | factor_interaction | 0.422 | 1.016 | 1.432 | 1.015 | 1.014 | 1.020 | 0.938 | 1.083 | 1.054 |
| atac | trend_additive | 0.414 | 1.019 | 1.458 | 1.010 | 1.028 | 1.021 | 0.933 | 1.088 | 1.041 |
| atac | trend_interaction | 0.395 | 1.020 | 1.457 | 1.017 | 1.025 | 1.020 | 0.935 | 1.085 | 1.044 |

### TMM check (pasted from `results/stage_de_limma/cohort_30/stage_de_limma.md`)

edgeR's TMM factors on the tested counts (library size = their column sums, the export's trims) against the export manifest's, absolute tolerance 1e-08:

- rna: largest absolute difference 1.11e-16 over 80 donors and 20387 features (passed)
- atac: largest absolute difference 2.22e-16 over 80 donors and 128916 features (passed)

### Joint F diagnostic (pasted from `results/stage_de_limma/cohort_30/stage_de_limma.md`)

limma's F is recomputed from the moderated t's and the unweighted design's coefficient correlation (it must equal topTable's), and compared with the exact F from each feature's own weighted coefficient covariance (same variance and degrees of freedom; BH over all tested features). limma's F stays the reported test.

| omic | model | test | df1 / df2 | limma F recomputed (max rel. diff) | exact F: max / median rel. diff | significant: limma / exact | changed calls |
|---|---|---|---|---:|---:|---:|---|
| rna | factor_additive | stage_joint | 2 / 75.93 | 7.6e-16 | 0.141 / 0.015 | 1 / 1 | 0: none |
| rna | factor_interaction | interaction_joint | 2 / 73.96 | 1.1e-15 | 0.312 / 0.0106 | 0 / 0 | 0: none |
| atac | factor_additive | stage_joint | 2 / 103.7 | 1.1e-15 | 0.414 / 0.0487 | 0 / 0 | 0: none |
| atac | factor_interaction | interaction_joint | 2 / 101.9 | 8.9e-16 | 0.531 / 0.0357 | 0 / 0 | 0: none |

### Fits (pasted from `results/stage_de_limma/cohort_30/stage_de_limma.md`)

| omic | model | features × donors | voom span | df prior | s² prior | runtime (s) |
|---|---|---:|---:|---:|---:|---:|
| rna | factor_additive | 20387 × 80 | 0.3944 | 4.933 | 0.7879 | 14.0 |
| rna | factor_interaction | 20387 × 80 | 0.3944 | 4.956 | 0.7874 | 14.2 |
| rna | trend_additive | 20387 × 80 | 0.3944 | 4.945 | 0.7892 | 12.1 |
| rna | trend_interaction | 20387 × 80 | 0.3944 | 4.959 | 0.789 | 12.1 |
| atac | factor_additive | 128916 × 80 | 0.3510 | 32.75 | 0.9156 | 89.5 |
| atac | factor_interaction | 128916 × 80 | 0.3510 | 32.85 | 0.9121 | 91.7 |
| atac | trend_additive | 128916 × 80 | 0.3510 | 32.77 | 0.9176 | 77.5 |
| atac | trend_interaction | 128916 × 80 | 0.3510 | 32.8 | 0.9158 | 77.5 |

### Findings (pasted from `results/stage_de_limma/cohort_30/stage_de_limma.md`)

- cohort_30 atac factor_additive stage_joint (joint_f_vs_union, abundance): method dependence: DESeq2 union detects 22, limma stage_joint detects 0; a result of neither method
- cohort_30 atac trend_additive stage_num (coefficient, abundance): method dependence: DESeq2 stage_num detects 23, limma stage_num detects 0; a result of neither method

## 3. Sensitivity cohort: cohort_50 (≥ 50 astrocyte nuclei in both omics, n = 74)

### Concordance (pasted from `results/stage_de_limma/cohort_50/stage_de_limma.md`)

Rule (`STAGE_DE_LIMMA_CONCORDANCE_RULE`, predeclared): `interaction_detected` = at least one exported feature (status kept) with adj.P.Val <= STAGE_DE_ALPHA on the joint moderated F of factor_interaction's two interaction coefficients (interaction_joint); `delta_rejected` = the primary analysis rejects delta (p <= ANALYSIS_ALPHA, ANALYSIS_REJECTION_RULE; result.json); concordant when both are true or both are false; discordant otherwise, written as a finding. Reported beside the DESeq2 cell (STAGE_DE_CONCORDANCE_RULE, the union of two Wald tests), which stays the reported cell of stage-de; methods_agree compares the two interaction_detected values and a disagreement is a method-dependence finding; neither cell replaces the other or changes a MOTCO reading. The primary analysis (`results/motco_analysis/cohort_50`, CV-selected rank 1, MOTCO 0.9.0) is read, never rerun or modified.

| omic | DESeq2: exported features (union of two Wald tests) | DESeq2 cell | limma: exported features (interaction_joint) | limma cell | delta rejected | methods agree |
|---|---:|---|---:|---|:---:|:---:|
| rna | 0 | concordant | 0 | concordant | no | yes |
| atac | 0 | concordant | 0 | concordant | no | yes |

### Counts: rna (pasted from `results/stage_de_limma/cohort_50/stage_de_limma.md`)

5000 exported and 20238 abundance-passing features tested; significant / tested at adjusted p ≤ α. The DESeq2 columns are the same coefficient's Wald test, the union for a union row, and DESeq2's union of the model's two Wald tests for a joint-F row (a different test, `joint_f_vs_union`).

| model | test | kind | design column | limma exported | limma abundance | DESeq2 exported | DESeq2 abundance | min adj. p limma (exported / abundance) |
|---|---|---|---|---:|---:|---:|---:|---:|
| factor_additive | stage_1_vs_0 | moderated t | `stage[T.1]` | 0 / 5000 | 0 / 20238 | 0 / 5000 | 0 / 20238 | 0.227 / 0.0509 |
| factor_additive | stage_2_vs_0 | moderated t | `stage[T.2]` | 1 / 5000 | 1 / 20238 | 2 / 5000 | 2 / 20238 | 0.0023 / 0.0023 |
| factor_additive | stage_joint | joint F | `stage[T.1];stage[T.2]` | 1 / 5000 | 2 / 20238 | 2 / 5000 | 2 / 20238 | 0.00622 / 0.00622 |
| factor_additive | union | union (no joint error rate) | `stage[T.1];stage[T.2]` | 1 / 5000 | 1 / 20238 | 2 / 5000 | 2 / 20238 | 0.0023 / 0.0023 |
| factor_interaction | male_x_stage_1 | moderated t | `sex[T.Male]:stage[T.1]` | 0 / 5000 | 0 / 20238 | 0 / 5000 | 0 / 20238 | 0.981 / 0.397 |
| factor_interaction | male_x_stage_2 | moderated t | `sex[T.Male]:stage[T.2]` | 0 / 5000 | 0 / 20238 | 0 / 5000 | 0 / 20238 | 0.749 / 0.192 |
| factor_interaction | interaction_joint | joint F | `sex[T.Male]:stage[T.1];sex[T.Male]:stage[T.2]` | 0 / 5000 | 0 / 20238 | 0 / 5000 | 0 / 20238 | 0.882 / 0.337 |
| factor_interaction | union | union (no joint error rate) | `sex[T.Male]:stage[T.1];sex[T.Male]:stage[T.2]` | 0 / 5000 | 0 / 20238 | 0 / 5000 | 0 / 20238 | 0.749 / 0.192 |
| trend_additive | stage_num | moderated t | `stage_num` | 1 / 5000 | 1 / 20238 | 2 / 5000 | 2 / 20238 | 0.00109 / 0.00109 |
| trend_interaction | male_x_stage_num | moderated t | `sex[T.Male]:stage_num` | 0 / 5000 | 0 / 20238 | 0 / 5000 | 0 / 20238 | 0.555 / 0.147 |

### Counts: atac (pasted from `results/stage_de_limma/cohort_50/stage_de_limma.md`)

5000 exported and 124818 abundance-passing features tested; significant / tested at adjusted p ≤ α. The DESeq2 columns are the same coefficient's Wald test, the union for a union row, and DESeq2's union of the model's two Wald tests for a joint-F row (a different test, `joint_f_vs_union`).

| model | test | kind | design column | limma exported | limma abundance | DESeq2 exported | DESeq2 abundance | min adj. p limma (exported / abundance) |
|---|---|---|---|---:|---:|---:|---:|---:|
| factor_additive | stage_1_vs_0 | moderated t | `stage[T.1]` | 0 / 5000 | 0 / 124818 | 0 / 5000 | 0 / 124818 | 0.763 / 0.329 |
| factor_additive | stage_2_vs_0 | moderated t | `stage[T.2]` | 0 / 5000 | 12 / 124818 | 4 / 5000 | 36 / 124818 | 0.0805 / 0.0325 |
| factor_additive | stage_joint | joint F | `stage[T.1];stage[T.2]` | 0 / 5000 | 0 / 124818 | 4 / 5000 | 36 / 124818 | 0.208 / 0.0965 |
| factor_additive | union | union (no joint error rate) | `stage[T.1];stage[T.2]` | 0 / 5000 | 12 / 124818 | 4 / 5000 | 36 / 124818 | 0.0805 / 0.0325 |
| factor_interaction | male_x_stage_1 | moderated t | `sex[T.Male]:stage[T.1]` | 0 / 5000 | 0 / 124818 | 0 / 5000 | 0 / 124818 | 0.719 / 0.0925 |
| factor_interaction | male_x_stage_2 | moderated t | `sex[T.Male]:stage[T.2]` | 0 / 5000 | 0 / 124818 | 0 / 5000 | 0 / 124818 | 0.969 / 0.957 |
| factor_interaction | interaction_joint | joint F | `sex[T.Male]:stage[T.1];sex[T.Male]:stage[T.2]` | 0 / 5000 | 0 / 124818 | 0 / 5000 | 0 / 124818 | 1 / 0.292 |
| factor_interaction | union | union (no joint error rate) | `sex[T.Male]:stage[T.1];sex[T.Male]:stage[T.2]` | 0 / 5000 | 0 / 124818 | 0 / 5000 | 0 / 124818 | 0.719 / 0.0925 |
| trend_additive | stage_num | moderated t | `stage_num` | 0 / 5000 | 7 / 124818 | 4 / 5000 | 30 / 124818 | 0.0934 / 0.0449 |
| trend_interaction | male_x_stage_num | moderated t | `sex[T.Male]:stage_num` | 0 / 5000 | 0 / 124818 | 0 / 5000 | 0 / 124818 | 1 / 1 |

### Method comparison: rna (pasted from `results/stage_de_limma/cohort_50/stage_de_limma.md`)

Features called by both methods and by one only, whether each method detects anything, and, per coefficient, the Spearman correlations of DESeq2's log2 fold change with limma's logFC and of the Wald statistic with the moderated t (over the features with both values finite). Descriptive: no p-value, no gate; a row whose detection differs is a method-dependence finding.

| model | test | kind | set | DESeq2 | limma | both | only DESeq2 | only limma | detection agrees | ρ log2FC | ρ statistic | compared |
|---|---|---|---|---:|---:|---:|---:|---:|:---:|---:|---:|---:|
| factor_additive | stage_1_vs_0 | coefficient | exported | 0 | 0 | 0 | 0 | 0 | yes | 0.806 | 0.827 | 5000 |
| factor_additive | stage_1_vs_0 | coefficient | abundance | 0 | 0 | 0 | 0 | 0 | yes | 0.843 | 0.877 | 20238 |
| factor_additive | stage_2_vs_0 | coefficient | exported | 2 | 1 | 1 | 1 | 0 | yes | 0.869 | 0.876 | 5000 |
| factor_additive | stage_2_vs_0 | coefficient | abundance | 2 | 1 | 1 | 1 | 0 | yes | 0.911 | 0.928 | 20238 |
| factor_additive | union | union | exported | 2 | 1 | 1 | 1 | 0 | yes | n/a | n/a |  |
| factor_additive | union | union | abundance | 2 | 1 | 1 | 1 | 0 | yes | n/a | n/a |  |
| factor_additive | stage_joint | joint_f_vs_union | exported | 2 | 1 | 1 | 1 | 0 | yes | n/a | n/a |  |
| factor_additive | stage_joint | joint_f_vs_union | abundance | 2 | 2 | 1 | 1 | 1 | yes | n/a | n/a |  |
| factor_interaction | male_x_stage_1 | coefficient | exported | 0 | 0 | 0 | 0 | 0 | yes | 0.813 | 0.826 | 5000 |
| factor_interaction | male_x_stage_1 | coefficient | abundance | 0 | 0 | 0 | 0 | 0 | yes | 0.840 | 0.866 | 20238 |
| factor_interaction | male_x_stage_2 | coefficient | exported | 0 | 0 | 0 | 0 | 0 | yes | 0.838 | 0.846 | 5000 |
| factor_interaction | male_x_stage_2 | coefficient | abundance | 0 | 0 | 0 | 0 | 0 | yes | 0.886 | 0.910 | 20238 |
| factor_interaction | union | union | exported | 0 | 0 | 0 | 0 | 0 | yes | n/a | n/a |  |
| factor_interaction | union | union | abundance | 0 | 0 | 0 | 0 | 0 | yes | n/a | n/a |  |
| factor_interaction | interaction_joint | joint_f_vs_union | exported | 0 | 0 | 0 | 0 | 0 | yes | n/a | n/a |  |
| factor_interaction | interaction_joint | joint_f_vs_union | abundance | 0 | 0 | 0 | 0 | 0 | yes | n/a | n/a |  |
| trend_additive | stage_num | coefficient | exported | 2 | 1 | 1 | 1 | 0 | yes | 0.870 | 0.876 | 5000 |
| trend_additive | stage_num | coefficient | abundance | 2 | 1 | 1 | 1 | 0 | yes | 0.909 | 0.927 | 20238 |
| trend_interaction | male_x_stage_num | coefficient | exported | 0 | 0 | 0 | 0 | 0 | yes | 0.830 | 0.841 | 5000 |
| trend_interaction | male_x_stage_num | coefficient | abundance | 0 | 0 | 0 | 0 | 0 | yes | 0.882 | 0.907 | 20238 |

### Method comparison: atac (pasted from `results/stage_de_limma/cohort_50/stage_de_limma.md`)

Features called by both methods and by one only, whether each method detects anything, and, per coefficient, the Spearman correlations of DESeq2's log2 fold change with limma's logFC and of the Wald statistic with the moderated t (over the features with both values finite). Descriptive: no p-value, no gate; a row whose detection differs is a method-dependence finding.

| model | test | kind | set | DESeq2 | limma | both | only DESeq2 | only limma | detection agrees | ρ log2FC | ρ statistic | compared |
|---|---|---|---|---:|---:|---:|---:|---:|:---:|---:|---:|---:|
| factor_additive | stage_1_vs_0 | coefficient | exported | 0 | 0 | 0 | 0 | 0 | yes | 0.815 | 0.821 | 5000 |
| factor_additive | stage_1_vs_0 | coefficient | abundance | 0 | 0 | 0 | 0 | 0 | yes | 0.867 | 0.887 | 124818 |
| factor_additive | stage_2_vs_0 | coefficient | exported | 4 | 0 | 0 | 4 | 0 | no | 0.831 | 0.840 | 5000 |
| factor_additive | stage_2_vs_0 | coefficient | abundance | 36 | 12 | 6 | 30 | 6 | yes | 0.888 | 0.903 | 124818 |
| factor_additive | union | union | exported | 4 | 0 | 0 | 4 | 0 | no | n/a | n/a |  |
| factor_additive | union | union | abundance | 36 | 12 | 6 | 30 | 6 | yes | n/a | n/a |  |
| factor_additive | stage_joint | joint_f_vs_union | exported | 4 | 0 | 0 | 4 | 0 | no | n/a | n/a |  |
| factor_additive | stage_joint | joint_f_vs_union | abundance | 36 | 0 | 0 | 36 | 0 | no | n/a | n/a |  |
| factor_interaction | male_x_stage_1 | coefficient | exported | 0 | 0 | 0 | 0 | 0 | yes | 0.777 | 0.790 | 5000 |
| factor_interaction | male_x_stage_1 | coefficient | abundance | 0 | 0 | 0 | 0 | 0 | yes | 0.855 | 0.875 | 124818 |
| factor_interaction | male_x_stage_2 | coefficient | exported | 0 | 0 | 0 | 0 | 0 | yes | 0.797 | 0.811 | 5000 |
| factor_interaction | male_x_stage_2 | coefficient | abundance | 0 | 0 | 0 | 0 | 0 | yes | 0.865 | 0.889 | 124818 |
| factor_interaction | union | union | exported | 0 | 0 | 0 | 0 | 0 | yes | n/a | n/a |  |
| factor_interaction | union | union | abundance | 0 | 0 | 0 | 0 | 0 | yes | n/a | n/a |  |
| factor_interaction | interaction_joint | joint_f_vs_union | exported | 0 | 0 | 0 | 0 | 0 | yes | n/a | n/a |  |
| factor_interaction | interaction_joint | joint_f_vs_union | abundance | 0 | 0 | 0 | 0 | 0 | yes | n/a | n/a |  |
| trend_additive | stage_num | coefficient | exported | 4 | 0 | 0 | 4 | 0 | no | 0.832 | 0.842 | 5000 |
| trend_additive | stage_num | coefficient | abundance | 30 | 7 | 2 | 28 | 5 | yes | 0.888 | 0.903 | 124818 |
| trend_interaction | male_x_stage_num | coefficient | exported | 0 | 0 | 0 | 0 | 0 | yes | 0.791 | 0.806 | 5000 |
| trend_interaction | male_x_stage_num | coefficient | abundance | 0 | 0 | 0 | 0 | 0 | yes | 0.859 | 0.885 | 124818 |

### Features at adjusted p ≤ 0.05 (generated from the per-test tables; descriptive)

Every feature at adjusted p ≤ 0.05 per limma test, by adjusted p (ties by p-value, then name); a test without one lists its five smallest adjusted p. DESeq2's padj for the same feature and coefficient stands beside it (for a joint F, DESeq2's smaller padj of the model's two coefficients; empty = removed by its independent filter). `exported` = one of MOTCO's 5,000 features.

- rna `factor_additive` `stage_1_vs_0`: none at adjusted p ≤ 0.05; smallest five:

  | feature | AveExpr | logFC | SE | t | adj. p | DESeq2 padj | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `GINS2` | 0.81 | -0.691 | 0.135 | -5.10 | 0.0509 | 0.197 | no |
  | `STAC3` | 1.41 | -0.792 | 0.160 | -4.95 | 0.0509 | 0.12 | no |
  | `NBPF10` | 4.21 | 0.477 | 0.110 | 4.36 | 0.213 | 0.437 | no |
  | `GP6` | 2.41 | -0.612 | 0.142 | -4.31 | 0.213 | 0.698 | no |
  | `TMEM176A` | 2.02 | 0.667 | 0.155 | 4.31 | 0.213 | 0.663 | no |

- rna `factor_additive` `stage_2_vs_0`: 1 at adjusted p ≤ 0.05 (1 exported):

  | feature | AveExpr | logFC | SE | t | adj. p | DESeq2 padj | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `CREB5` | 4.93 | 1.129 | 0.191 | 5.91 | 0.0023 | 3.21e-05 | yes |

- rna `factor_additive` `stage_joint`: 2 at adjusted p ≤ 0.05 (1 exported):

  | feature | AveExpr | logFC_stage_1_vs_0 | logFC_stage_2_vs_0 | F | adj. p | DESeq2 padj (min) | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `CREB5` | 4.93 | 0.366 | 1.129 | 18.73 | 0.00622 | 3.21e-05 | yes |
  | `STAC3` | 1.41 | -0.792 | -0.655 | 14.79 | 0.0446 | 0.0554 | no |

- rna `factor_interaction` `male_x_stage_1`: none at adjusted p ≤ 0.05; smallest five:

  | feature | AveExpr | logFC | SE | t | adj. p | DESeq2 padj | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `SLC25A24` | 3.83 | 1.816 | 0.395 | 4.59 | 0.397 | 0.999 | no |
  | `BMP8A` | 0.48 | -1.443 | 0.367 | -3.93 | 0.981 | 0.999 | no |
  | `ORAI1` | 0.86 | -1.311 | 0.341 | -3.84 | 0.981 | 0.999 | no |
  | `CD36` | 0.43 | 1.513 | 0.401 | 3.77 | 0.981 | 0.999 | no |
  | `PHYH` | 3.73 | -0.627 | 0.167 | -3.76 | 0.981 | 0.999 | no |

- rna `factor_interaction` `male_x_stage_2`: none at adjusted p ≤ 0.05; smallest five:

  | feature | AveExpr | logFC | SE | t | adj. p | DESeq2 padj | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `DNMT1` | 4.27 | -0.687 | 0.147 | -4.69 | 0.192 | 0.569 | no |
  | `AL160408.1` | -0.71 | 2.573 | 0.559 | 4.60 | 0.192 | 0.72 | no |
  | `TAS2R31` | 0.29 | 1.577 | 0.383 | 4.12 | 0.717 | 0.72 | no |
  | `AL109976.1` | -0.32 | -1.902 | 0.479 | -3.97 | 0.749 | 0.72 | yes |
  | `MTX1` | 2.16 | -0.947 | 0.243 | -3.90 | 0.749 | 0.653 | no |

- rna `factor_interaction` `interaction_joint`: none at adjusted p ≤ 0.05; smallest five:

  | feature | AveExpr | logFC_male_x_stage_1 | logFC_male_x_stage_2 | F | adj. p | DESeq2 padj (min) | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `SLC25A24` | 3.83 | 1.816 | 0.340 | 12.99 | 0.337 | 0.91 | no |
  | `DNMT1` | 4.27 | -0.542 | -0.687 | 11.30 | 0.437 | 0.569 | no |
  | `AL160408.1` | -0.71 | 1.155 | 2.573 | 11.15 | 0.437 | 0.72 | no |
  | `BMP8A` | 0.48 | -1.443 | -0.216 | 9.96 | 0.818 | 0.93 | no |
  | `AC091074.2` | -0.13 | 0.456 | 1.962 | 9.15 | 0.882 | 0.653 | yes |

- rna `trend_additive` `stage_num`: 1 at adjusted p ≤ 0.05 (1 exported):

  | feature | AveExpr | logFC | SE | t | adj. p | DESeq2 padj | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `CREB5` | 4.93 | 0.581 | 0.096 | 6.08 | 0.00109 | 2.01e-05 | yes |

- rna `trend_interaction` `male_x_stage_num`: none at adjusted p ≤ 0.05; smallest five:

  | feature | AveExpr | logFC | SE | t | adj. p | DESeq2 padj | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `AL160408.1` | -0.71 | 1.290 | 0.266 | 4.85 | 0.147 | 0.649 | no |
  | `DNMT1` | 4.27 | -0.317 | 0.071 | -4.45 | 0.326 | 0.567 | no |
  | `TAS2R31` | 0.29 | 0.772 | 0.184 | 4.20 | 0.531 | 0.649 | no |
  | `AC091074.2` | -0.13 | 0.995 | 0.249 | 4.00 | 0.555 | 0.447 | yes |
  | `AC133785.1` | 0.95 | 0.750 | 0.188 | 3.99 | 0.555 | 0.117 | yes |

- atac `factor_additive` `stage_1_vs_0`: none at adjusted p ≤ 0.05; smallest five:

  | feature | AveExpr | logFC | SE | t | adj. p | DESeq2 padj | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `chr1:245589383-245590226` | 2.19 | -0.894 | 0.179 | -5.00 | 0.329 | 0.63 | no |
  | `chr18:47983898-47984482` | -0.20 | -1.521 | 0.322 | -4.72 | 0.511 | 1 | no |
  | `chr5:10592071-10592905` | -1.38 | 1.727 | 0.376 | 4.59 | 0.563 | 0.905 | no |
  | `chr8:119152756-119153569` | -0.60 | 1.719 | 0.391 | 4.40 | 0.695 | 0.63 | no |
  | `chr2:178392322-178392845` | -1.13 | 1.815 | 0.418 | 4.34 | 0.695 | 1 | no |

- atac `factor_additive` `stage_2_vs_0`: 12 at adjusted p ≤ 0.05 (0 exported):

  | feature | AveExpr | logFC | SE | t | adj. p | DESeq2 padj | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `chr12:108750482-108752712` | 2.90 | -0.834 | 0.156 | -5.33 | 0.0325 | 0.0841 | no |
  | `chr6:11605152-11608407` | 3.22 | -0.723 | 0.137 | -5.27 | 0.0325 | 0.0383 | no |
  | `chr12:106328829-106330471` | 3.72 | -0.511 | 0.099 | -5.18 | 0.0325 | 0.033 | no |
  | `chr3:194880040-194881125` | 4.02 | -0.567 | 0.111 | -5.10 | 0.0325 | 0.0486 | no |
  | `chr1:231053444-231054140` | 2.04 | -0.822 | 0.162 | -5.08 | 0.0325 | 0.0486 | no |
  | `chr19:15279962-15281620` | 3.48 | -0.557 | 0.110 | -5.04 | 0.0325 | 0.0728 | no |
  | `chr1:245589383-245590226` | 2.19 | -0.830 | 0.165 | -5.04 | 0.0325 | 0.055 | no |
  | `chr4:123554731-123556429` | 3.75 | -0.578 | 0.116 | -4.96 | 0.0325 | 0.0728 | no |
  | `chr1:155481563-155482044` | 2.05 | -1.095 | 0.221 | -4.96 | 0.0325 | 0.106 | no |
  | `chr3:143134550-143135514` | 2.81 | -0.776 | 0.157 | -4.96 | 0.0325 | 0.0354 | no |
  | `chr1:200749095-200750613` | 1.98 | -0.837 | 0.169 | -4.96 | 0.0325 | 0.0341 | no |
  | `chr17:265751-268156` | 4.24 | -0.626 | 0.126 | -4.95 | 0.0325 | 0.0998 | no |

- atac `factor_additive` `stage_joint`: none at adjusted p ≤ 0.05; smallest five:

  | feature | AveExpr | logFC_stage_1_vs_0 | logFC_stage_2_vs_0 | F | adj. p | DESeq2 padj (min) | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `chr1:245589383-245590226` | 2.19 | -0.894 | -0.830 | 16.36 | 0.0965 | 0.055 | no |
  | `chr12:108750482-108752712` | 2.90 | -0.635 | -0.834 | 14.81 | 0.105 | 0.0841 | no |
  | `chr6:11605152-11608407` | 3.22 | -0.402 | -0.723 | 13.90 | 0.105 | 0.0383 | no |
  | `chr19:15279962-15281620` | 3.48 | -0.478 | -0.557 | 13.86 | 0.105 | 0.0728 | no |
  | `chr18:47983898-47984482` | -0.20 | -1.521 | -1.257 | 13.82 | 0.105 | 0.125 | no |

- atac `factor_interaction` `male_x_stage_1`: none at adjusted p ≤ 0.05; smallest five:

  | feature | AveExpr | logFC | SE | t | adj. p | DESeq2 padj | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `chr17:304482-305830` | 1.97 | -2.046 | 0.385 | -5.31 | 0.0925 | 1 | no |
  | `chr12:48753307-48754450` | 1.39 | 1.821 | 0.393 | 4.63 | 0.428 | 0.714 | no |
  | `chr7:90932085-90932904` | -0.68 | -3.090 | 0.672 | -4.60 | 0.428 | 0.714 | no |
  | `chr19:12513516-12514280` | 1.96 | 2.350 | 0.512 | 4.59 | 0.428 | 1 | no |
  | `chr19:18414054-18442660` | 7.48 | 0.595 | 0.132 | 4.51 | 0.467 | 1 | no |

- atac `factor_interaction` `male_x_stage_2`: none at adjusted p ≤ 0.05; smallest five:

  | feature | AveExpr | logFC | SE | t | adj. p | DESeq2 padj | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `chr17:46259862-46261495` | 1.82 | 1.841 | 0.410 | 4.49 | 0.957 | 0.516 | no |
  | `chr10:52778669-52779188` | 0.08 | 2.408 | 0.558 | 4.32 | 0.957 | 0.516 | no |
  | `chr9:112913366-112914928` | -1.55 | 3.151 | 0.732 | 4.31 | 0.957 | 0.789 | no |
  | `chr5:16664584-16665195` | -0.09 | -2.474 | 0.575 | -4.30 | 0.957 | 0.392 | no |
  | `chr5:179318308-179318574` | 0.30 | 2.161 | 0.510 | 4.24 | 0.957 | 0.589 | no |

- atac `factor_interaction` `interaction_joint`: none at adjusted p ≤ 0.05; smallest five:

  | feature | AveExpr | logFC_male_x_stage_1 | logFC_male_x_stage_2 | F | adj. p | DESeq2 padj (min) | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `chr17:304482-305830` | 1.97 | -2.046 | -0.729 | 14.94 | 0.292 | 0.894 | no |
  | `chr17:46259862-46261495` | 1.82 | 1.859 | 1.841 | 11.66 | 0.973 | 0.516 | no |
  | `chr5:16664584-16665195` | -0.09 | -2.560 | -2.474 | 11.50 | 0.973 | 0.392 | no |
  | `chr12:48753307-48754450` | 1.39 | 1.821 | 1.277 | 11.34 | 0.973 | 0.648 | no |
  | `chr19:12513516-12514280` | 1.96 | 2.350 | 1.635 | 10.96 | 0.973 | 0.854 | no |

- atac `trend_additive` `stage_num`: 7 at adjusted p ≤ 0.05 (0 exported):

  | feature | AveExpr | logFC | SE | t | adj. p | DESeq2 padj | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `chr6:11605152-11608407` | 3.22 | -0.358 | 0.068 | -5.28 | 0.0449 | 0.0481 | no |
  | `chr12:106328829-106330471` | 3.72 | -0.256 | 0.049 | -5.25 | 0.0449 | 0.0407 | no |
  | `chr12:108750482-108752712` | 2.90 | -0.406 | 0.079 | -5.11 | 0.0449 | 0.137 | no |
  | `chr3:194880040-194881125` | 4.02 | -0.280 | 0.055 | -5.11 | 0.0449 | 0.0559 | no |
  | `chr1:155481563-155482044` | 2.05 | -0.549 | 0.108 | -5.07 | 0.0449 | 0.116 | no |
  | `chr17:265751-268156` | 4.24 | -0.312 | 0.062 | -5.01 | 0.0449 | 0.122 | no |
  | `chr1:231053444-231054140` | 2.04 | -0.405 | 0.081 | -5.00 | 0.0449 | 0.0784 | no |

- atac `trend_interaction` `male_x_stage_num`: none at adjusted p ≤ 0.05; smallest five:

  | feature | AveExpr | logFC | SE | t | adj. p | DESeq2 padj | exported |
  |---|---:|---:|---:|---:|---:|---:|:---:|
  | `chr1:33725071-33727945` | -0.86 | 1.557 | 0.348 | 4.47 | 1 | 0.701 | no |
  | `chr9:112913366-112914928` | -1.55 | 1.486 | 0.350 | 4.25 | 1 | 0.81 | no |
  | `chr12:75120269-75121753` | -1.16 | -1.424 | 0.342 | -4.16 | 1 | 0.575 | no |
  | `chr5:179318308-179318574` | 0.30 | 0.994 | 0.240 | 4.15 | 1 | 0.662 | no |
  | `chr2:47157842-47158490` | -1.33 | -1.441 | 0.351 | -4.10 | 1 | 0.662 | no |

### Sample weights (pasted from `results/stage_de_limma/cohort_50/stage_de_limma.md`)

One quality weight per donor and model (arrayWeights inside voomWithQualityWeights, shrunk toward 1 with `prior.n = 10`): minimum, median, maximum and the mean per sex × stage cell. Shown without a gate; a cell with a low mean weight counts less in that model's fit.

| omic | model | min | median | max | Female/stage 0 | Female/stage 1 | Female/stage 2 | Male/stage 0 | Male/stage 1 | Male/stage 2 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| rna | factor_additive | 0.361 | 0.991 | 2.309 | 1.244 | 1.312 | 0.919 | 1.114 | 1.176 | 1.096 |
| rna | factor_interaction | 0.347 | 0.978 | 2.564 | 1.235 | 1.337 | 0.938 | 1.109 | 1.139 | 1.075 |
| rna | trend_additive | 0.364 | 1.028 | 2.138 | 1.247 | 1.258 | 0.920 | 1.131 | 1.168 | 1.096 |
| rna | trend_interaction | 0.370 | 1.021 | 2.100 | 1.248 | 1.248 | 0.928 | 1.140 | 1.164 | 1.089 |
| atac | factor_additive | 0.383 | 1.034 | 1.468 | 1.030 | 1.032 | 1.018 | 0.928 | 1.108 | 1.036 |
| atac | factor_interaction | 0.371 | 1.027 | 1.442 | 1.035 | 1.025 | 1.018 | 0.934 | 1.102 | 1.042 |
| atac | trend_additive | 0.373 | 1.040 | 1.462 | 1.033 | 1.039 | 1.017 | 0.925 | 1.109 | 1.032 |
| atac | trend_interaction | 0.357 | 1.042 | 1.431 | 1.040 | 1.034 | 1.016 | 0.930 | 1.105 | 1.034 |

### TMM check (pasted from `results/stage_de_limma/cohort_50/stage_de_limma.md`)

edgeR's TMM factors on the tested counts (library size = their column sums, the export's trims) against the export manifest's, absolute tolerance 1e-08:

- rna: largest absolute difference 1.11e-16 over 74 donors and 20238 features (passed)
- atac: largest absolute difference 1.11e-16 over 74 donors and 124818 features (passed)

### Joint F diagnostic (pasted from `results/stage_de_limma/cohort_50/stage_de_limma.md`)

limma's F is recomputed from the moderated t's and the unweighted design's coefficient correlation (it must equal topTable's), and compared with the exact F from each feature's own weighted coefficient covariance (same variance and degrees of freedom; BH over all tested features). limma's F stays the reported test.

| omic | model | test | df1 / df2 | limma F recomputed (max rel. diff) | exact F: max / median rel. diff | significant: limma / exact | changed calls |
|---|---|---|---|---:|---:|---:|---|
| rna | factor_additive | stage_joint | 2 / 69.9 | 9.7e-16 | 0.17 / 0.0107 | 2 / 2 | 0: none |
| rna | factor_interaction | interaction_joint | 2 / 67.91 | 8.5e-16 | 0.287 / 0.00725 | 0 / 0 | 0: none |
| atac | factor_additive | stage_joint | 2 / 95.67 | 1.3e-15 | 0.427 / 0.0536 | 0 / 0 | 0: none |
| atac | factor_interaction | interaction_joint | 2 / 93.73 | 1.5e-15 | 0.577 / 0.0356 | 0 / 0 | 0: none |

### Fits (pasted from `results/stage_de_limma/cohort_50/stage_de_limma.md`)

| omic | model | features × donors | voom span | df prior | s² prior | runtime (s) |
|---|---|---:|---:|---:|---:|---:|
| rna | factor_additive | 20238 × 74 | 0.3946 | 4.903 | 0.7898 | 12.7 |
| rna | factor_interaction | 20238 × 74 | 0.3946 | 4.912 | 0.789 | 12.9 |
| rna | trend_additive | 20238 × 74 | 0.3946 | 4.917 | 0.7909 | 11.1 |
| rna | trend_interaction | 20238 × 74 | 0.3946 | 4.918 | 0.7907 | 11.6 |
| atac | factor_additive | 124818 × 74 | 0.3516 | 30.67 | 0.9028 | 80.6 |
| atac | factor_interaction | 124818 × 74 | 0.3516 | 30.73 | 0.8987 | 82.3 |
| atac | trend_additive | 124818 × 74 | 0.3516 | 30.7 | 0.905 | 69.0 |
| atac | trend_interaction | 124818 × 74 | 0.3516 | 30.69 | 0.9029 | 69.7 |

### Findings (pasted from `results/stage_de_limma/cohort_50/stage_de_limma.md`)

- cohort_50 atac factor_additive stage_2_vs_0 (coefficient, exported): method dependence: DESeq2 stage_2_vs_0 detects 4, limma stage_2_vs_0 detects 0; a result of neither method
- cohort_50 atac factor_additive union (union, exported): method dependence: DESeq2 union detects 4, limma union detects 0; a result of neither method
- cohort_50 atac factor_additive stage_joint (joint_f_vs_union, exported): method dependence: DESeq2 union detects 4, limma stage_joint detects 0; a result of neither method
- cohort_50 atac factor_additive stage_joint (joint_f_vs_union, abundance): method dependence: DESeq2 union detects 36, limma stage_joint detects 0; a result of neither method
- cohort_50 atac trend_additive stage_num (coefficient, exported): method dependence: DESeq2 stage_num detects 4, limma stage_num detects 0; a result of neither method

## 4. Reproducibility

- Rerun (task 3.2): `stage-de-limma --cohort cohort_30 --omic both --force` into a scratch directory (7 min 3 s) gave `stage_de_limma_summary.csv`, `method_comparison.csv` and all 16 tables byte-identical to the reported run's (`cmp`, 18 of 18 files); the inputs and DESeq2 outputs hashed as before. R runs with one BLAS thread, so a rerun on the same inputs and versions reproduces the tables exactly.

- Cluster copy (task 3.3): `results/stage_de_limma/` (56 files, 358 MB) was copied to `/home1/tgonza/sead_sexdiff/results/` and all 56 files verify there with `sha256sum -c` against sums computed on nixos-home.

## 5. Synthetic verification and calibration (simulations under DESeq2's own model)

Verification on the synthetic dataset of `mtg-astrocyte-stage-de`, rebuilt to its recorded specification (60 donors, 10 per sex × stage cell; 50 planted stage features and 10 female-only interaction features per omic; the earlier scripts were on the offline laptop). Over 16 fits the design columns equalled the DESeq2 step's, the planted stage features filled 47–49 of the 50 smallest adjusted p of `stage_2_vs_0`, `stage_num` and `stage_joint`, null features had p ≤ 0.05 in 0.043–0.056 of cases, the sample weights stayed within 0.90–1.10, and the DESeq2 step's tables and summaries stayed byte-identical with the changed package. With the planted interaction concentrated in one coefficient, `interaction_joint` recalled 0.40–0.80 of the planted interaction features against 0.70–1.00 for the union of the two moderated t's.

Pooled realized FDR and recall at adjusted p ≤ 0.05 over 6 seeds, additive tables (`stage_1_vs_0`, `stage_2_vs_0`, `stage_num`; non-null = any planted effect), pasted from `calibration_additive.csv`:

| condition | method | significant | false | FDR | recall |
|---|---|---:|---:|---:|---:|
| logunif_atacsyn | deseq2 | 1036 | 84 | 0.081 | 0.881 |
| logunif_atacsyn | limma_qw | 990 | 38 | 0.038 | 0.881 |
| logunif_atacsyn | voom | 987 | 35 | 0.035 | 0.881 |
| logunif_atacsyn_het | deseq2 | 1007 | 85 | 0.084 | 0.854 |
| logunif_atacsyn_het | limma_qw | 956 | 43 | 0.045 | 0.845 |
| logunif_atacsyn_het | voom | 911 | 27 | 0.030 | 0.819 |
| trend_rna | deseq2 | 3305 | 257 | 0.078 | 0.847 |
| trend_rna | limma_qw | 3087 | 167 | 0.054 | 0.811 |
| trend_rna | voom | 3097 | 175 | 0.057 | 0.812 |
| trend_rna_het | deseq2 | 3013 | 177 | 0.059 | 0.788 |
| trend_rna_het | limma_qw | 2948 | 148 | 0.050 | 0.778 |
| trend_rna_het | voom | 2726 | 109 | 0.040 | 0.727 |

Joint and interaction tables (pasted from `calibration_pooled.csv`; limma methods only for the joint F):

| condition | method | model | test | significant | false | FDR | recall |
|---|---|---|---|---:|---:|---:|---:|
| logunif_atacsyn | deseq2 | factor_interaction | male_x_stage_2 | 65 | 9 | 0.138 | 0.933 |
| logunif_atacsyn | limma_qw | factor_additive | stage_joint | 362 | 12 | 0.033 | 0.972 |
| logunif_atacsyn | limma_qw | factor_interaction | interaction_joint | 48 | 1 | 0.021 | 0.783 |
| logunif_atacsyn | limma_qw | factor_interaction | male_x_stage_2 | 55 | 1 | 0.018 | 0.900 |
| logunif_atacsyn | voom | factor_additive | stage_joint | 362 | 12 | 0.033 | 0.972 |
| logunif_atacsyn | voom | factor_interaction | interaction_joint | 48 | 0 | 0.000 | 0.800 |
| logunif_atacsyn | voom | factor_interaction | male_x_stage_2 | 56 | 2 | 0.036 | 0.900 |
| logunif_atacsyn_het | deseq2 | factor_interaction | male_x_stage_2 | 56 | 4 | 0.071 | 0.867 |
| logunif_atacsyn_het | limma_qw | factor_additive | stage_joint | 378 | 25 | 0.066 | 0.981 |
| logunif_atacsyn_het | limma_qw | factor_interaction | interaction_joint | 44 | 1 | 0.023 | 0.717 |
| logunif_atacsyn_het | limma_qw | factor_interaction | male_x_stage_2 | 52 | 3 | 0.058 | 0.817 |
| logunif_atacsyn_het | voom | factor_additive | stage_joint | 357 | 10 | 0.028 | 0.964 |
| logunif_atacsyn_het | voom | factor_interaction | interaction_joint | 40 | 1 | 0.025 | 0.650 |
| logunif_atacsyn_het | voom | factor_interaction | male_x_stage_2 | 49 | 1 | 0.020 | 0.800 |

The two homogeneous conditions are rebuilt from the earlier change's records (its simulation scripts and tables were on the offline laptop), so DESeq2 was refitted on the same simulated datasets; its realized FDR there (0.078 / 0.081) is near the 0.06 / 0.09 that change recorded. `trend_rna`: 4,000 features, 200 non-null, dispersions on the cohort_30 RNA DESeq2 trend with log-normal scatter; `logunif_atacsyn`: the synthetic ATAC set's 5,830 peaks and 60 planted features, dispersions log-uniform on [0.05, 0.10]; `_het`: one donor per sex × stage cell at 4× dispersion (their mean quality weights 0.27–0.38, the other donors' 1.11–1.17). Every condition follows DESeq2's own count model with covariates that have no effect, which favours DESeq2; the heterogeneous conditions were built for the quality-weights rationale. The predeclared stop point (limma-voom with quality weights above DESeq2's pooled FDR on the homogeneous additive tables) was not reached. None of these numbers is a gate on, or a property of, the real blocks.

## 6. Reading

- **Concordant in every cell, in both methods.** limma's joint F finds no exported feature with a sex × stage interaction at adjusted p ≤ 0.05 in either cohort or omic, as DESeq2's union of two Wald tests found none, and the primary analysis rejected `delta` in neither cohort. All eight cells (four per method) are `concordant`, and `methods_agree` holds in every cohort and omic. The joint interaction test that pydeseq2 could not run detects nothing the union missed. The `delta` non-rejection stays uninformative, and so does the absence of DE interactions.
- **No sex × stage interaction detected at this n, in either method or model.** On the abundance-passing set the smallest `interaction_joint` adjusted p is 0.55 (cohort_30 rna), 0.4 (cohort_30 atac), 0.337 (cohort_50 rna) and 0.292 (cohort_50 atac); the smallest on a single interaction coefficient (`male_x_stage_1`, `male_x_stage_2`, `male_x_stage_num`) is 0.718, 0.0934, 0.147 and 0.0925. With cells of 9–28 donors this is "none detected at this n". The joint F is not more powerful than the union when an interaction sits on one coefficient: on the synthetic data it recalled 0.40–0.80 of the planted interaction features against 0.70–1.00 for the union.
- **The RNA stage signal is the same gene in both methods.** `CREB5` (exported: yes) is the only RNA gene limma calls at stage 2 vs 0 and on the stage trend in both cohorts (cohort_30 adjusted p 0.000463 and 0.000204; DESeq2 padj 7.17e-06 and 4.34e-06), and it carries `stage_joint` too (0.00107; cohort_50 0.00622). DESeq2's second gene, `ADAM33`, is not called by limma (cohort_30 `stage_num` adjusted p 0.121 against DESeq2's padj 0.0484; cohort_50 `stage_2_vs_0` 0.144 against 0.00973), and limma's `stage_joint` adds `STAC3` in cohort_50 (not exported; adjusted p 0.0446). Both are method dependence at the margin and a result of neither method.
- **The ATAC stage peaks depend on the method.** At stage 2 vs 0 DESeq2 calls 22 (cohort_30) and 36 (cohort_50) abundance-passing peaks and limma 9 and 12, of which 3 and 6 are shared; on the stage trend DESeq2 calls 23 and 30, limma 0 and 7 (0 and 2 shared). Of MOTCO's exported peaks DESeq2 calls 4 in cohort_50 at stage 2 vs 0 and limma none. limma's calls all sit at one adjusted p (0.0457 in cohort_30, 0.0325 in cohort_50 at stage 2 vs 0; 0.0449 on the cohort_50 trend), BH ties at the edge of α, and DESeq2's peaks have limma adjusted p 0.0457–0.301 (cohort_30) and 0.0325–0.278 (cohort_50). Part of the difference is the BH family: DESeq2's independent filter removed 29,993 (cohort_30) and 38,719 (cohort_50) of the 128,916 / 124,818 tested peaks from the stage-2 family, while limma adjusts over all of them. It is not the whole difference: on the cohort_50 trend DESeq2 removed 0 peaks and still calls 30 against limma's 7, so the normalization, the weights and the test statistic contribute too. By the reading rules these peaks are a result of neither method. Both methods agree that the ATAC stage signal is weak, at the edge of α and confined to stage 2 (no call at stage 1 vs 0 in either method or cohort), and limma's `stage_joint` finds none (0 / 0 against 9 / 12 on the union): as on the synthetic data, the 2-df test spends power on a stage-1 coefficient on which no peak is called. The 7 method-dependence findings of the manifests are these ATAC rows.
- **Agreement on the effect estimates, not only on the calls.** Over the abundance-passing features the Spearman correlation of DESeq2's log2 fold change with limma's logFC is 0.81–0.92 (rna) and 0.85–0.89 (atac) in cohort_30 and 0.84–0.91 / 0.85–0.89 in cohort_50, and of the Wald statistic with the moderated t 0.85–0.93 / 0.87–0.90 and 0.87–0.93 / 0.88–0.90. The two methods order the features similarly; they differ in where the few features near α fall.
- **Quality weights.** The per-donor weights range 0.35–2.56 in RNA and 0.36–1.47 in ATAC. In RNA the Female/stage 2 cell has the lowest mean weight in every model and cohort (0.92–1.00, against 1.07–1.34 for the other cells), the case design Risks names (quality weights that down-weight a sex × stage cell, e.g. where low quality concentrates at stage 2). These donors count less in limma's fits, which is shown without a gate; the concordance reading does not change with it. RNA's `df.prior` is 4.9–5.0 (ATAC 30.7–32.9): the RNA residual variances are heterogeneous, and the moderation borrows little strength.
- **The joint F's approximation moves no call.** limma's F, recomputed from the moderated t's and the unweighted design's correlation, equals topTable's to ≤ 1.5e-15. The exact F from each feature's own weighted covariance differs from it by at most 0.14–0.58 (relative; the range of the per-test maxima) and changes 0 calls at α in all eight joint tests.

## 7. Limitations

- This change was proposed after the DESeq2 results were known. Its method, tests, BH family, concordance rule and stop point were fixed in the proposal and design before any limma run, and the DESeq2 results stay the reported ones; a second method chosen after a result is still a choice made after a result, which is why a difference between the methods is read as method dependence only.
- The same donors and counts enter both methods: agreement is not replication, and `cohort_50` is a donor subset of `cohort_30`.
- Power: n = 80 / 74 with sex × stage cells of 9–28 donors; a zero interaction count is "none detected at this n", and the synthetic recalls (dispersions 0.05–0.10, below the real fits') overstate the real power.
- The union counts carry no joint error-rate claim. The joint F's claim is approximate with voom's observation weights (the exact-F diagnostic measures the approximation, and limma's F stays the reported test). limma's FDR on these blocks is not known: the calibration of section 5 followed DESeq2's own count model on rebuilt synthetic conditions (realized FDR of the additive tables 0.054 and 0.038 for limma-voom with quality weights, 0.078 and 0.081 for DESeq2), with covariates that have no effect, and is no property of the real blocks.
- The earlier change's simulation scripts and tables were not available (offline machine); its two conditions were rebuilt from its records, so section 5 compares the methods on the same rebuilt datasets rather than against that change's stored figures.
- RNA's prior degrees of freedom are low (4.9–5.0, beside the residual degrees of freedom of each fit), so the RNA residual variances are moderated little and the moderated t's are close to ordinary t's; robust empirical Bayes was rejected in the design to keep the predeclared canonical pipeline and was not run.
- The quality weights are one per donor and model; they down-weight a noisy donor in every feature, not an outlying observation, and limma has no Cook's-type protection. They read the design, so a cell's mean weight can shift with the model (shown per model).
- voom's log-CPM is not the export's prior-2 log2-CPM; the normalization is the export's TMM (checked to ≤ 2.2e-16) while DESeq2's was median of ratios; limma's BH families are all tested features while DESeq2's exclude its filtered features.
- The R side ran from the Nix shell (R 4.6.1, limma 3.68.4, edgeR 4.10.1); the versions and the script's sha256 are in every fit record, and a record made under other versions is stale.
- The run was made on nixos-home, not the workstation the tasks named; the inputs are bound by hash to the cluster's copies, and the DESeq2 outputs by the step's own reference check.
- The `delta` side of the concordance cells is MOTCO's joint latent-space test with the six analysis covariates; every limitation of the primary report (`reports/motco-astrocyte-2026-10-19.md` and its sensitivity reports) applies to it.

## 8. Next steps

- Archive `mtg-astrocyte-stage-de-limma` and sync the `pseudobulk-stage-de-limma` spec after the merge.
- Then amend `mtg-microglia-positive-control` (branch `origin/mtg-microglia-positive-control`, rebased onto main, via `/opsx:update`) before any microglia run: `--cell-type` and the `cell_type` manifest check on `stage-de-limma` too (eight MTG commands), `stage-de-limma` in its D10 byte-identity baseline, the limma concordance rule for microglia beside the DESeq2 one, and a predeclared choice for its expectation E1 (DESeq2's `factor_additive` union as written, limma's `stage_joint`, or both) with how a method disagreement enters its all-three rule. That choice belongs to the microglia proposal and is made there, before any microglia run.
- Nothing here changes the primary MOTCO analysis, which stays primary, or the DESeq2 results, which stay the reported ones.

## 9. Reproduction

```
scripts/sync_pseudobulk.sh                                    # results/pseudobulk/ from the cluster
rsync -a <cluster>:/home1/tgonza/sead_sexdiff/results/{motco_input,motco_analysis,stage_de} results/
nix develop -c uv run sead-sexdiff stage-de-limma --cohort all --omic both
# one fit by hand: keep the exchange files, then rerun the R script on a model directory
uv run python -c 'from sead_sexdiff.data import run_stage_de_limma; run_stage_de_limma("cohort_30", omic="rna", out_dir="/tmp/sl", debug_dir="/tmp/sl_debug", force=True)'
Rscript --vanilla src/sead_sexdiff/stage_de_limma.R /tmp/sl_debug/cohort_30_rna/factor_interaction
```


# MOTCO astrocyte feature-selection sensitivity: sex × ADNC stage — 2026-10-19

This report reads the five feature-content variants of the registry `PREPROCESSING_VARIANTS` added by the change `mtg-astrocyte-feature-selection-sensitivity` (`exclude_families`, `rna_trend_ratio`, `atac_mean_accessibility`, `residualized` and the composed `revised_candidate`), each the primary export with the registered setting(s) changed and the covariate-adjusted MOTCO analysis rerun unchanged on the variant blocks, beside the primary analysis of `reports/motco-astrocyte-2026-10-19.md`. The chrX variant is in `reports/motco-astrocyte-chrx-sensitivity-2026-10-19.md` and the diagnostic that motivated these variants (the first latent axis of the primary analysis is a ribosomal/mitochondrial composition axis) in `reports/motco-astrocyte-technical-axis-2026-10-19.md`. Every number below is pasted from the generated files (`sensitivity.md`, `sensitivity_manifest.json`, `adoption_check.csv`, `axis_correlations.csv`, `geometry_shares.csv`); the primary analysis stays primary, every variant is sensitivity to it, and the one decision this report carries, the adoption check of the candidate, is a predeclared rule over generated numbers that reads no p-value.

**Reading rules (binding).** The primary analysis remains the primary result; no variant is chosen from the results. `angle` and `shape` non-rejections are uninformative in every row (they are not evidence that the sexes agree in direction or shape); a `delta` rejection is a difference, not specifically of size. A rejection status that differs between a variant and the primary is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics (`../MOTCO/docs/roadmap.md`, Phase 6) were measured on neither family-excluded, trend- or mean-selected nor residualized blocks, so no variant carries a Type I or power claim of its own. A latent axis correlated with a technical variable is a measurement-space finding about the analysis in question, not a reinterpretation of its statistics.

## 1. Provenance

- Software: MOTCO 0.8.0 (git dependency pinned in `uv.lock`), scikit-learn 1.8.0, numpy 2.3.5, pandas 2.3.3, Python 3.12.15; sead-sexdiff revision `b5d108296fa53eac02304385e0f81ae138a5e410` (dirty: True).
- Settings shared by every row: seed 20261019, 999 RRPP permutations (the same permutation indices in every analysis), α = 0.05 (p <= alpha), `n_jobs = 1`; covariates age_at_death, pmi, rin, years_education, log10_n_rna, log10_n_atac entered additively in MOTCO's regression on the latent scores (unchanged in every variant).
- Inputs: pseudobulk manifest sha256 `2e93cc0b526856fc043f908a424a11266187146636bc75afed7edd187989a25c` (created 2026-10-16T18:46:13+00:00), covariates manifest `b8bf9091b507a1098af91c9390f487ccb0c78aa87fba07edc0bfd1d6b83184ae`, annotation GENCODE v32 basic sha256 `a61ab50fc8b5d63f76d50e1f89bc0f5bba0b64fb29d2ac8cfa1b1ef8d4ef0c91`, technical table `eac1b9a15977dfab683738a4bd0b0ffc207b896a102032409b60fccbfb96247b` and its manifest `c2bdb34499644580f22ae7f6c48a7f846c4dabad5b2b04182f8e9d5a98a106e7` (recorded pseudobulk manifest sha256 `2e93cc0b526856fc043f908a424a11266187146636bc75afed7edd187989a25c`, equal to the export's).
- cohort_30: primary export `results/motco_input/cohort_30` (manifest sha256 `18f2289e0d4c84ed18f4a437ae569052544e9d150ed2726385273694b4a3b0b6`, k = 9, 80 donors), primary analysis `results/motco_analysis/cohort_30` (`result.json` sha256 `ca92876a1a61c40791524390e8603125bcbfacf6f081d34e2ae6a9b49fd15de6`, CV-selected rank 1, cv_mean_auroc 0.6313, MOTCO 0.8.0, seed 20261019, 999 permutations); both read and checked, never rerun.
- cohort_50: primary export `results/motco_input/cohort_50` (manifest sha256 `dd14a6e7d5583fb527e8a99c7db43d12cf237f5451522e8875c5a0a3401c5cb3`, k = 9, 74 donors), primary analysis `results/motco_analysis/cohort_50` (`result.json` sha256 `1e7eefe3e63e7b9b51c3fed618c32d4a496ea3d382c7ff025c2a34f5d0ae048f`, CV-selected rank 1, cv_mean_auroc 0.6616, MOTCO 0.8.0, seed 20261019, 999 permutations); both read and checked, never rerun.
- Variants (registry entries; settings as the manifests record them):
  - `exclude_families`: RNA genes of the ribosomal, mitochondrial and hemoglobin sets (TECHNICAL_GENE_SETS) removed before the abundance filter; the ATAC block and every other export setting equal to the primary. Settings: `{"excluded_gene_sets": ["ribosomal", "mitochondrial", "hemoglobin"]}`.
  - `rna_trend_ratio`: RNA features ranked at the selection step by the ratio of their log2-CPM variance to a binned-median mean–variance trend instead of by raw variance; the ATAC block and every other export setting equal to the primary. Settings: `{"selection_rule": {"rna": "trend_ratio", "atac": "variance"}}`.
  - `atac_mean_accessibility`: ATAC peaks ranked at the selection step by mean log2-CPM (total accessibility) instead of by variance; the RNA block and every other export setting equal to the primary. Settings: `{"selection_rule": {"rna": "variance", "atac": "mean"}}`.
  - `residualized`: after log2-CPM and before the selection step, each abundance-passing feature replaced by its mean plus the residual of its regression on an intercept and RESIDUALIZATION_COVARIATES (the export's log10 library sizes and the technical table's gene-set fractions); the selection ranks the residual variance; every other export setting equal to the primary. Settings: `{"residualize": ["log10_lib_rna", "log10_lib_atac", "ribosomal_pseudobulk", "mitochondrial_pseudobulk"]}`.
  - `revised_candidate`: the four settings together in the primary's order: gene sets removed before the abundance filter, residualization after log2-CPM, trend-ratio ranking for RNA and mean-accessibility ranking for ATAC at the selection step; the adoption rule ADOPTION_RULE is evaluated for this variant. Settings: `{"excluded_gene_sets": ["ribosomal", "mitochondrial", "hemoglobin"], "residualize": ["log10_lib_rna", "log10_lib_atac", "ribosomal_pseudobulk", "mitochondrial_pseudobulk"], "selection_rule": {"rna": "trend_ratio", "atac": "mean"}}`.
- Residualization set `RESIDUALIZATION_COVARIATES` = ['log10_lib_rna', 'log10_lib_atac', 'ribosomal_pseudobulk', 'mitochondrial_pseudobulk']; excluded as cell state: `hemoglobin_pseudobulk` (reports/motco-astrocyte-technical-axis-2026-10-19.md, section 2, 'Classification under the spec's rule': median within-donor supertype share 0.00394 above the between-donor share 0.00285 (both negligible; median pseudobulk fraction 6.7e-07 of UMIs); not a residualization candidate).
- Adoption rule `ADOPTION_RULE`: holds when no latent component of the candidate's CV-selected analysis reaches |Spearman| >= DIAGNOSTIC_CORRELATION_THRESHOLD against any variable of ADOPTION_CHECK_VARIABLES, in both cohorts; p-values, rejection statuses and CV AUROC are not read. If it holds, a later change may propose the candidate as a revised primary with the current primary kept as the reference; if not, the variants stay sensitivity analyses. Check variables `ADOPTION_CHECK_VARIABLES` = ['tmm_rna', 'tmm_atac', 'fraction_mitochondrial_umis_mean', 'doublet_score_mean', 'gex_reads_mapped_confidently_to_intronic_regions_mean', 'gex_fraction_of_transcriptomic_reads_in_cells_mean', 'rin', 'pmi']; threshold 0.3.
- Descriptive checks: `axis-diagnostics` run on each variant root into a scratch directory (section 4), reading the variant export, the variant analysis and the technical table; its outputs feed no export or analysis command.

## 2. Primary cohort: cohort_30 (n = 80)

### Overview (generated from `sensitivity.csv` and `sensitivity_manifest.json`)

| analysis | rank | cv_mean_auroc | delta (p) | angle (p) | shape (p) | rejected | rna kept shared / only primary / only variant | atac kept shared / only primary / only variant | findings |
|---|---:|---:|---:|---:|---:|---|---|---|---|
| primary | 1 | 0.6313 | 14.941 (0.3) | 0 (1) | 0.608832 (0.333) | none | — | — | — |
| exclude_families | 1 | 0.6313 | 14.9196 (0.301) | 0 (1) | 0.609158 (0.334) | none | 4990 / 10 / 10 | 5000 / 0 / 0 | none beyond the export's unmapped-genes note |
| rna_trend_ratio | 1 | 0.6326 | 13.101 (0.354) | 0 (1) | 0.610539 (0.315) | none | 2536 / 2464 / 2464 | 5000 / 0 / 0 | none beyond the export's unmapped-genes note |
| atac_mean_accessibility | 1 | 0.6346 | 2.2945 (0.933) | 0 (1) | 0.885871 (0.31) | none | 5000 / 0 / 0 | 0 / 5000 / 5000 | none beyond the export's unmapped-genes note |
| residualized | 1 | 0.4078 | 4.04454 (0.39) | 0 (1) | 0.0442514 (0.795) | none | 4252 / 748 / 748 | 3447 / 1553 / 1553 | none beyond the export's unmapped-genes note |
| revised_candidate | 1 | 0.4822 | 1.45283 (0.829) | 0 (1) | 0.224695 (0.514) | none | 2282 / 2718 / 2718 | 0 / 5000 / 5000 | adoption check holds in cohort_30: no candidate component reaches |Spearman| >= 0.3 against any of 8 check variables (the rule needs both cohorts; a revised primary is a later change) |

### Attribution: the first latent axis per analysis (generated from `axis_correlations.csv`)

Spearman correlation of `pls_0` with the variables named in the task (stage; the astrocyte ribosomal-protein UMI fraction `ribosomal_pseudobulk`; the RNA TMM factor `tmm_rna`; the ATAC log10 library size `log10_lib_atac`), plus the mitochondrial fraction and sex, per analysis. For the residualized exports the library-size regressor was the step-3 library size and the fractions entered the design, so their near-zero correlations are by construction; `tmm_rna` was not a regressor and is a check variable. Descriptive: no p-value is attached to any correlation.

| analysis | rank | stage | ribosomal fraction | RNA TMM factor | ATAC log10 library size | mitochondrial fraction | sex (male) | component SDs | share of the sex-difference geometry along pls_0 |
|---|---:|---:|---:|---:|---:|---:|---:|---|---:|
| primary | 1 | 0.49 | -0.81 | -0.73 | 0.21 | 0.73 | -0.23 | pls_0 35.0 | 1.000 |
| exclude_families | 1 | 0.49 | -0.81 | -0.73 | 0.21 | 0.73 | -0.23 | pls_0 35.0 | 1.000 |
| rna_trend_ratio | 1 | 0.50 | -0.81 | -0.72 | 0.20 | 0.71 | -0.23 | pls_0 34.3 | 1.000 |
| atac_mean_accessibility | 1 | 0.36 | -0.70 | -0.64 | -0.18 | 0.70 | -0.31 | pls_0 54.4 | 1.000 |
| residualized | 1 | 0.75 | -0.05 | -0.08 | 0.01 | 0.09 | 0.01 | pls_0 13.9 | 1.000 |
| revised_candidate | 1 | 0.61 | -0.06 | -0.19 | 0.03 | 0.13 | -0.00 | pls_0 16.6 | 1.000 |

### Adoption check of `revised_candidate` (pasted from `adoption_check.csv`)

Rule: holds when no latent component of the candidate's CV-selected analysis reaches |Spearman| >= DIAGNOSTIC_CORRELATION_THRESHOLD against any variable of ADOPTION_CHECK_VARIABLES, in both cohorts; p-values, rejection statuses and CV AUROC are not read. If it holds, a later change may propose the candidate as a revised primary with the current primary kept as the reference; if not, the variants stay sensitivity analyses. Threshold 0.3. Outcome in cohort_30: **holds = true** (0 of 8 cells flagged).

| component | variable | source | n | spearman | flag | holds |
|---|---|---|---:|---:|:---:|:---:|
| pls_0 | tmm_rna | manifest.json (primary export, TMM factor) | 80 | -0.18335677449601503 | False |  |
| pls_0 | tmm_atac | manifest.json (primary export, TMM factor) | 80 | 0.16530707923112986 | False |  |
| pls_0 | fraction_mitochondrial_umis_mean | donor_technical.csv | 80 | 0.1391701828410689 | False |  |
| pls_0 | doublet_score_mean | donor_technical.csv | 80 | 0.0012658227848101266 | False |  |
| pls_0 | gex_reads_mapped_confidently_to_intronic_regions_mean | donor_technical.csv | 80 | -0.21857840236194812 | False |  |
| pls_0 | gex_fraction_of_transcriptomic_reads_in_cells_mean | donor_technical.csv | 80 | -0.12182372245663385 | False |  |
| pls_0 | rin | covariates.csv (primary export) | 80 | -0.033373971799641314 | False |  |
| pls_0 | pmi | covariates.csv (primary export) | 80 | -0.020000469609411455 | False |  |
| all | all |  | 80 |  | False | True |

### Variant `exclude_families` (pasted from `results/motco_sensitivity/exclude_families/comparison/cohort_30/sensitivity.md`)

#### MOTCO preprocessing sensitivity: cohort_30, variant `exclude_families`

RNA genes of the ribosomal, mitochondrial and hemoglobin sets (TECHNICAL_GENE_SETS) removed before the abundance filter; the ATAC block and every other export setting equal to the primary. Export setting(s) changed: `excluded_gene_sets` = ['ribosomal', 'mitochondrial', 'hemoglobin'] (primary: `excluded_gene_sets` = []); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 1, which stays the primary result; the `exclude_families` row (CV-selected rank 1) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

##### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 14.941 | 0.3 | no | 11.2868 | 8.27982 | 26.5775 |
| exclude_families | cv | 1 | 14.9196 | 0.301 | no | 11.2826 | 8.2719 | 26.5365 |

##### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0 | 1 | no | 3.06306 | 23.2919 | 0 |
| exclude_families | cv | 1 | 0 | 1 | no | 3.06306 | 23.2919 | 0 |

##### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0.608832 | 0.333 | no | 0.489849 | 0.355498 | 1.18287 |
| exclude_families | cv | 1 | 0.609158 | 0.334 | no | 0.49017 | 0.355641 | 1.18419 |

##### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after gene-set removal | after abundance | kept (after selection) | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 35329 | 20387 | 5000 | 4990 | 10 |
| rna | exclude_families | 36601 | 35329 | 35249 | 20317 | 5000 | 4990 | 10 |
| atac | primary | 218882 | 212816 | 212816 | 128916 | 5000 | 5000 | 0 |
| atac | exclude_families | 218882 | 212816 | 212816 | 128916 | 5000 | 5000 | 0 |

##### Gene sets removed

RNA genes matched by the excluded sets of `TECHNICAL_GENE_SETS` are removed at step 2 after the chromosome removal (status `excluded_gene_set` in `features_rna.csv`); a matched gene on an excluded chromosome keeps `sex_chromosome` and is counted under `already on an excluded chromosome`. ATAC names match no set.

| omic | set | matched | removed here | already on an excluded chromosome |
|---|---|---:|---:|---:|
| rna | ribosomal | 81 | 77 | 4 |
| rna | mitochondrial | 13 | 0 | 13 |
| rna | hemoglobin | 3 | 3 | 0 |
| atac | ribosomal | 0 | 0 | 0 |
| atac | mitochondrial | 0 | 0 | 0 |
| atac | hemoglobin | 0 | 0 | 0 |

##### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 1 | 0.6313 | 148.7 | `results/motco_analysis/cohort_30` |
| exclude_families | cv | 1 | 0.6313 | 114.5 | `motco_analysis/cohort_30` |

##### Findings

- exclude_families export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

### Variant `rna_trend_ratio` (pasted from `results/motco_sensitivity/rna_trend_ratio/comparison/cohort_30/sensitivity.md`)

#### MOTCO preprocessing sensitivity: cohort_30, variant `rna_trend_ratio`

RNA features ranked at the selection step by the ratio of their log2-CPM variance to a binned-median mean–variance trend instead of by raw variance; the ATAC block and every other export setting equal to the primary. Export setting(s) changed: `selection_rule` = {'rna': 'trend_ratio', 'atac': 'variance'} (primary: `selection_rule` = {'rna': 'variance', 'atac': 'variance'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 1, which stays the primary result; the `rna_trend_ratio` row (CV-selected rank 1) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

##### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 14.941 | 0.3 | no | 11.2868 | 8.27982 | 26.5775 |
| rna_trend_ratio | cv | 1 | 13.101 | 0.354 | no | 11.0604 | 8.19554 | 26.1634 |

##### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0 | 1 | no | 3.06306 | 23.2919 | 0 |
| rna_trend_ratio | cv | 1 | 0 | 1 | no | 2.7027 | 21.9012 | 0 |

##### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0.608832 | 0.333 | no | 0.489849 | 0.355498 | 1.18287 |
| rna_trend_ratio | cv | 1 | 0.610539 | 0.315 | no | 0.472458 | 0.342789 | 1.14017 |

##### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after abundance | kept (after selection) | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 20387 | 5000 | 2536 | 2464 |
| rna | rna_trend_ratio | 36601 | 35329 | 20387 | 5000 | 2536 | 2464 |
| atac | primary | 218882 | 212816 | 128916 | 5000 | 5000 | 0 |
| atac | rna_trend_ratio | 218882 | 212816 | 128916 | 5000 | 5000 | 0 |

##### Selection rule

`variance`: population variance (ddof 0) of log2-CPM, the primary's rule. `trend_ratio`: variance over a mean–variance trend, the median variance of 20 equal-count bins of mean log2-CPM placed at each bin's median mean and linearly interpolated (flat beyond the end bins; the bins are in the variant export's manifest). `mean`: mean log2-CPM. Ties by name under every rule; the quantiles are of the kept features' mean log2-CPM.

| omic | export | rule | kept mean q05 | q25 | q50 | q75 | q95 |
|---|---|---|---:|---:|---:|---:|---:|
| rna | primary | variance | -0.721 | -0.137 | 0.774 | 2.160 | 4.968 |
| rna | rna_trend_ratio | trend_ratio | 0.084 | 2.086 | 3.819 | 5.362 | 7.521 |

##### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 1 | 0.6313 | 148.7 | `results/motco_analysis/cohort_30` |
| rna_trend_ratio | cv | 1 | 0.6326 | 106.7 | `motco_analysis/cohort_30` |

##### Findings

- rna_trend_ratio export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

### Variant `atac_mean_accessibility` (pasted from `results/motco_sensitivity/atac_mean_accessibility/comparison/cohort_30/sensitivity.md`)

#### MOTCO preprocessing sensitivity: cohort_30, variant `atac_mean_accessibility`

ATAC peaks ranked at the selection step by mean log2-CPM (total accessibility) instead of by variance; the RNA block and every other export setting equal to the primary. Export setting(s) changed: `selection_rule` = {'rna': 'variance', 'atac': 'mean'} (primary: `selection_rule` = {'rna': 'variance', 'atac': 'variance'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 1, which stays the primary result; the `atac_mean_accessibility` row (CV-selected rank 1) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

##### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 14.941 | 0.3 | no | 11.2868 | 8.27982 | 26.5775 |
| atac_mean_accessibility | cv | 1 | 2.2945 | 0.933 | no | 18.2786 | 13.8278 | 46.0121 |

##### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0 | 1 | no | 3.06306 | 23.2919 | 0 |
| atac_mean_accessibility | cv | 1 | 0 | 1 | no | 15.8559 | 51.0417 | 180 |

##### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0.608832 | 0.333 | no | 0.489849 | 0.355498 | 1.18287 |
| atac_mean_accessibility | cv | 1 | 0.885871 | 0.31 | no | 0.652588 | 0.400874 | 1.31615 |

##### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after abundance | kept (after selection) | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 20387 | 5000 | 5000 | 0 |
| rna | atac_mean_accessibility | 36601 | 35329 | 20387 | 5000 | 5000 | 0 |
| atac | primary | 218882 | 212816 | 128916 | 5000 | 0 | 5000 |
| atac | atac_mean_accessibility | 218882 | 212816 | 128916 | 5000 | 0 | 5000 |

##### Selection rule

`variance`: population variance (ddof 0) of log2-CPM, the primary's rule. `trend_ratio`: variance over a mean–variance trend, the median variance of 20 equal-count bins of mean log2-CPM placed at each bin's median mean and linearly interpolated (flat beyond the end bins; the bins are in the variant export's manifest). `mean`: mean log2-CPM. Ties by name under every rule; the quantiles are of the kept features' mean log2-CPM.

| omic | export | rule | kept mean q05 | q25 | q50 | q75 | q95 |
|---|---|---|---:|---:|---:|---:|---:|
| atac | primary | variance | -0.321 | 0.046 | 0.371 | 0.778 | 1.531 |
| atac | atac_mean_accessibility | mean | 5.436 | 5.613 | 5.876 | 6.274 | 7.156 |

##### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 1 | 0.6313 | 148.7 | `results/motco_analysis/cohort_30` |
| atac_mean_accessibility | cv | 1 | 0.6346 | 101.1 | `motco_analysis/cohort_30` |

##### Findings

- atac_mean_accessibility export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

### Variant `residualized` (pasted from `results/motco_sensitivity/residualized/comparison/cohort_30/sensitivity.md`)

#### MOTCO preprocessing sensitivity: cohort_30, variant `residualized`

after log2-CPM and before the selection step, each abundance-passing feature replaced by its mean plus the residual of its regression on an intercept and RESIDUALIZATION_COVARIATES (the export's log10 library sizes and the technical table's gene-set fractions); the selection ranks the residual variance; every other export setting equal to the primary. Export setting(s) changed: `residualize` = ['log10_lib_rna', 'log10_lib_atac', 'ribosomal_pseudobulk', 'mitochondrial_pseudobulk'] (primary: `residualize` = None); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 1, which stays the primary result; the `residualized` row (CV-selected rank 1) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

##### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 14.941 | 0.3 | no | 11.2868 | 8.27982 | 26.5775 |
| residualized | cv | 1 | 4.04454 | 0.39 | no | 3.78699 | 2.8887 | 9.30296 |

##### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0 | 1 | no | 3.06306 | 23.2919 | 0 |
| residualized | cv | 1 | 0 | 1 | no | 0 | 0 | 0 |

##### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0.608832 | 0.333 | no | 0.489849 | 0.355498 | 1.18287 |
| residualized | cv | 1 | 0.0442514 | 0.795 | no | 0.15126 | 0.115125 | 0.364029 |

##### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after abundance | kept (after selection) | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 20387 | 5000 | 4252 | 748 |
| rna | residualized | 36601 | 35329 | 20387 | 5000 | 4252 | 748 |
| atac | primary | 218882 | 212816 | 128916 | 5000 | 3447 | 1553 |
| atac | residualized | 218882 | 212816 | 128916 | 5000 | 3447 | 1553 |

##### Residualization

After log2-CPM and before the selection step, each abundance-passing feature was replaced by its mean plus the residual of its OLS fit on an intercept and 4 covariates: `log10_lib_rna`, `log10_lib_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` (`RESIDUALIZATION_COVARIATES`: the export's own log10 library sizes after step 2 and the technical table's gene-set pseudobulk fractions; excluded as cell state under the technical-axis report's rule: `hemoglobin_pseudobulk` (reports/motco-astrocyte-technical-axis-2026-10-19.md, section 2, 'Classification under the spec's rule': median within-donor supertype share 0.00394 above the between-donor share 0.00285 (both negligible; median pseudobulk fraction 6.7e-07 of UMIs); not a residualization candidate)). The design (`residualization.csv` beside the variant blocks, sha256 `fdb31ed02c18…`) holds no sex, stage or analysis covariate. Limitations: the fit spends 5 degrees of freedom per feature (intercept plus covariates) before the permutation null, which does not see them; it removes whatever part of the stage and sex signal the covariates carry; the analysis covariates were left unchanged, so the library-size proxies are adjusted twice (log10 library sizes here, log10 nuclei counts in the RRPP model); and no Type I check exists for the residualized blocks.

| omic | R² q05 | q25 | q50 | q75 | q95 |
|---|---:|---:|---:|---:|---:|
| rna | 0.039 | 0.105 | 0.186 | 0.304 | 0.528 |
| atac | 0.034 | 0.081 | 0.131 | 0.198 | 0.365 |

##### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 1 | 0.6313 | 148.7 | `results/motco_analysis/cohort_30` |
| residualized | cv | 1 | 0.4078 | 104.9 | `motco_analysis/cohort_30` |

##### Findings

- residualized export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

### Variant `revised_candidate` (pasted from `results/motco_sensitivity/revised_candidate/comparison/cohort_30/sensitivity.md`)

#### MOTCO preprocessing sensitivity: cohort_30, variant `revised_candidate`

the four settings together in the primary's order: gene sets removed before the abundance filter, residualization after log2-CPM, trend-ratio ranking for RNA and mean-accessibility ranking for ATAC at the selection step; the adoption rule ADOPTION_RULE is evaluated for this variant. Export setting(s) changed: `excluded_gene_sets` = ['ribosomal', 'mitochondrial', 'hemoglobin'], `residualize` = ['log10_lib_rna', 'log10_lib_atac', 'ribosomal_pseudobulk', 'mitochondrial_pseudobulk'], `selection_rule` = {'rna': 'trend_ratio', 'atac': 'mean'} (primary: `excluded_gene_sets` = [], `residualize` = None, `selection_rule` = {'rna': 'variance', 'atac': 'variance'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 1, which stays the primary result; the `revised_candidate` row (CV-selected rank 1) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

##### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 14.941 | 0.3 | no | 11.2868 | 8.27982 | 26.5775 |
| revised_candidate | cv | 1 | 1.45283 | 0.829 | no | 5.64583 | 4.22126 | 13.6914 |

##### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0 | 1 | no | 3.06306 | 23.2919 | 0 |
| revised_candidate | cv | 1 | 0 | 1 | no | 0 | 0 | 0 |

##### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0.608832 | 0.333 | no | 0.489849 | 0.355498 | 1.18287 |
| revised_candidate | cv | 1 | 0.224695 | 0.514 | no | 0.269654 | 0.202313 | 0.647269 |

##### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after gene-set removal | after abundance | kept (after selection) | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 35329 | 20387 | 5000 | 2282 | 2718 |
| rna | revised_candidate | 36601 | 35329 | 35249 | 20317 | 5000 | 2282 | 2718 |
| atac | primary | 218882 | 212816 | 212816 | 128916 | 5000 | 0 | 5000 |
| atac | revised_candidate | 218882 | 212816 | 212816 | 128916 | 5000 | 0 | 5000 |

##### Gene sets removed

RNA genes matched by the excluded sets of `TECHNICAL_GENE_SETS` are removed at step 2 after the chromosome removal (status `excluded_gene_set` in `features_rna.csv`); a matched gene on an excluded chromosome keeps `sex_chromosome` and is counted under `already on an excluded chromosome`. ATAC names match no set.

| omic | set | matched | removed here | already on an excluded chromosome |
|---|---|---:|---:|---:|
| rna | ribosomal | 81 | 77 | 4 |
| rna | mitochondrial | 13 | 0 | 13 |
| rna | hemoglobin | 3 | 3 | 0 |
| atac | ribosomal | 0 | 0 | 0 |
| atac | mitochondrial | 0 | 0 | 0 |
| atac | hemoglobin | 0 | 0 | 0 |

##### Selection rule

`variance`: population variance (ddof 0) of log2-CPM, the primary's rule. `trend_ratio`: variance over a mean–variance trend, the median variance of 20 equal-count bins of mean log2-CPM placed at each bin's median mean and linearly interpolated (flat beyond the end bins; the bins are in the variant export's manifest). `mean`: mean log2-CPM. Ties by name under every rule; the quantiles are of the kept features' mean log2-CPM.

| omic | export | rule | kept mean q05 | q25 | q50 | q75 | q95 |
|---|---|---|---:|---:|---:|---:|---:|
| rna | primary | variance | -0.721 | -0.137 | 0.774 | 2.160 | 4.968 |
| rna | revised_candidate | trend_ratio | 0.127 | 2.387 | 3.969 | 5.478 | 7.641 |
| atac | primary | variance | -0.321 | 0.046 | 0.371 | 0.778 | 1.531 |
| atac | revised_candidate | mean | 5.436 | 5.613 | 5.876 | 6.274 | 7.156 |

##### Residualization

After log2-CPM and before the selection step, each abundance-passing feature was replaced by its mean plus the residual of its OLS fit on an intercept and 4 covariates: `log10_lib_rna`, `log10_lib_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` (`RESIDUALIZATION_COVARIATES`: the export's own log10 library sizes after step 2 and the technical table's gene-set pseudobulk fractions; excluded as cell state under the technical-axis report's rule: `hemoglobin_pseudobulk` (reports/motco-astrocyte-technical-axis-2026-10-19.md, section 2, 'Classification under the spec's rule': median within-donor supertype share 0.00394 above the between-donor share 0.00285 (both negligible; median pseudobulk fraction 6.7e-07 of UMIs); not a residualization candidate)). The design (`residualization.csv` beside the variant blocks, sha256 `90fc1f5d567b…`) holds no sex, stage or analysis covariate. Limitations: the fit spends 5 degrees of freedom per feature (intercept plus covariates) before the permutation null, which does not see them; it removes whatever part of the stage and sex signal the covariates carry; the analysis covariates were left unchanged, so the library-size proxies are adjusted twice (log10 library sizes here, log10 nuclei counts in the RRPP model); and no Type I check exists for the residualized blocks.

| omic | R² q05 | q25 | q50 | q75 | q95 |
|---|---:|---:|---:|---:|---:|
| rna | 0.039 | 0.105 | 0.185 | 0.302 | 0.522 |
| atac | 0.034 | 0.081 | 0.131 | 0.198 | 0.365 |

##### Adoption check

`ADOPTION_RULE`: holds when no latent component of the candidate's CV-selected analysis reaches |Spearman| >= DIAGNOSTIC_CORRELATION_THRESHOLD against any variable of ADOPTION_CHECK_VARIABLES, in both cohorts; p-values, rejection statuses and CV AUROC are not read. If it holds, a later change may propose the candidate as a revised primary with the current primary kept as the reference; if not, the variants stay sensitivity analyses. Threshold 0.3 (`DIAGNOSTIC_CORRELATION_THRESHOLD`). In this cohort the rule holds. The rule needs both cohorts; whether the candidate becomes a revised primary is a later change with the current primary kept as the reference, and the candidate's statistics above are sensitivity, not the result.

| component | variable | source | n | Spearman | flagged |
|---|---|---|---:|---:|:---:|
| pls_0 | `tmm_rna` | manifest.json (primary export, TMM factor) | 80 | -0.183 | no |
| pls_0 | `tmm_atac` | manifest.json (primary export, TMM factor) | 80 | 0.165 | no |
| pls_0 | `fraction_mitochondrial_umis_mean` | donor_technical.csv | 80 | 0.139 | no |
| pls_0 | `doublet_score_mean` | donor_technical.csv | 80 | 0.001 | no |
| pls_0 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 80 | -0.219 | no |
| pls_0 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 80 | -0.122 | no |
| pls_0 | `rin` | covariates.csv (primary export) | 80 | -0.033 | no |
| pls_0 | `pmi` | covariates.csv (primary export) | 80 | -0.020 | no |
| all | all | | 80 | | no → holds = true |

##### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 1 | 0.6313 | 148.7 | `results/motco_analysis/cohort_30` |
| revised_candidate | cv | 1 | 0.4822 | 102.8 | `motco_analysis/cohort_30` |

##### Findings

- adoption check holds in cohort_30: no candidate component reaches |Spearman| >= 0.3 against any of 8 check variables (the rule needs both cohorts; a revised primary is a later change)
- revised_candidate export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

## 3. Sensitivity cohort: cohort_50 (n = 74)

### Overview (generated from `sensitivity.csv` and `sensitivity_manifest.json`)

| analysis | rank | cv_mean_auroc | delta (p) | angle (p) | shape (p) | rejected | rna kept shared / only primary / only variant | atac kept shared / only primary / only variant | findings |
|---|---:|---:|---:|---:|---:|---|---|---|---|
| primary | 1 | 0.6616 | 18.7091 (0.213) | 0 (1) | 0.56471 (0.377) | none | — | — | — |
| exclude_families | 1 | 0.6674 | 19.0327 (0.202) | 0 (1) | 0.563235 (0.377) | none | 4977 / 23 / 23 | 5000 / 0 / 0 | none beyond the export's unmapped-genes note |
| rna_trend_ratio | 1 | 0.6716 | 18.0822 (0.211) | 0 (1) | 0.559086 (0.353) | none | 2572 / 2428 / 2428 | 5000 / 0 / 0 | none beyond the export's unmapped-genes note |
| atac_mean_accessibility | 1 | 0.6820 | 9.1382 (0.665) | 0 (1) | 0.645472 (0.465) | none | 5000 / 0 / 0 | 0 / 5000 / 5000 | none beyond the export's unmapped-genes note |
| residualized | 2 | 0.4791 | 8.94403 (0.444) | 30.0936 (0.104) | 0.317151 (0.338) | none | 4218 / 782 / 782 | 3423 / 1577 / 1577 | CV-selected rank differs: 1 in the primary analysis, 2 in the residualized analysis; reported beside angle, not read as better or worse |
| revised_candidate | 2 | 0.5514 | 6.33891 (0.585) | 2.15631 (0.902) | 0.230511 (0.564) | none | 2300 / 2700 / 2700 | 0 / 5000 / 5000 | adoption check does not hold in cohort_50: pls_0 x tmm_atac (0.37), pls_0 x gex_reads_mapped_confidently_to_intronic_regions_mean (-0.33), pls_1 x tmm_atac (-0.37); the variants stay sensitivity analyses; CV-selected rank differs: 1 in the primary analysis, 2 in the revised_candidate analysis; reported beside angle, not read as better or worse |

### Attribution: the first latent axis per analysis (generated from `axis_correlations.csv`)

Spearman correlation of `pls_0` with the variables named in the task (stage; the astrocyte ribosomal-protein UMI fraction `ribosomal_pseudobulk`; the RNA TMM factor `tmm_rna`; the ATAC log10 library size `log10_lib_atac`), plus the mitochondrial fraction and sex, per analysis. For the residualized exports the library-size regressor was the step-3 library size and the fractions entered the design, so their near-zero correlations are by construction; `tmm_rna` was not a regressor and is a check variable. Descriptive: no p-value is attached to any correlation.

| analysis | rank | stage | ribosomal fraction | RNA TMM factor | ATAC log10 library size | mitochondrial fraction | sex (male) | component SDs | share of the sex-difference geometry along pls_0 |
|---|---:|---:|---:|---:|---:|---:|---:|---|---:|
| primary | 1 | 0.50 | -0.82 | -0.72 | 0.23 | 0.77 | -0.19 | pls_0 35.8 | 1.000 |
| exclude_families | 1 | 0.50 | -0.81 | -0.70 | 0.23 | 0.77 | -0.19 | pls_0 35.9 | 1.000 |
| rna_trend_ratio | 1 | 0.51 | -0.82 | -0.72 | 0.22 | 0.74 | -0.18 | pls_0 35.3 | 1.000 |
| atac_mean_accessibility | 1 | 0.36 | -0.73 | -0.62 | -0.11 | 0.73 | -0.28 | pls_0 56.0 | 1.000 |
| residualized | 2 | 0.68 | -0.03 | -0.16 | -0.05 | 0.16 | 0.01 | pls_0 16.0, pls_1 30.0 | 0.128 |
| revised_candidate | 2 | 0.52 | -0.08 | -0.29 | 0.02 | 0.18 | 0.01 | pls_0 21.8, pls_1 30.5 | 0.230 |

### Adoption check of `revised_candidate` (pasted from `adoption_check.csv`)

Rule: holds when no latent component of the candidate's CV-selected analysis reaches |Spearman| >= DIAGNOSTIC_CORRELATION_THRESHOLD against any variable of ADOPTION_CHECK_VARIABLES, in both cohorts; p-values, rejection statuses and CV AUROC are not read. If it holds, a later change may propose the candidate as a revised primary with the current primary kept as the reference; if not, the variants stay sensitivity analyses. Threshold 0.3. Outcome in cohort_50: **holds = false** (3 of 16 cells flagged: pls_0 × `tmm_atac` (0.375), pls_0 × `gex_reads_mapped_confidently_to_intronic_regions_mean` (-0.325), pls_1 × `tmm_atac` (-0.371)).

| component | variable | source | n | spearman | flag | holds |
|---|---|---|---:|---:|:---:|:---:|
| pls_0 | tmm_rna | manifest.json (primary export, TMM factor) | 74 | -0.2872269529803776 | False |  |
| pls_0 | tmm_atac | manifest.json (primary export, TMM factor) | 74 | 0.37469085523880036 | True |  |
| pls_0 | fraction_mitochondrial_umis_mean | donor_technical.csv | 74 | 0.20393928174750087 | False |  |
| pls_0 | doublet_score_mean | donor_technical.csv | 74 | 0.03386893743058126 | False |  |
| pls_0 | gex_reads_mapped_confidently_to_intronic_regions_mean | donor_technical.csv | 74 | -0.32545224326877603 | True |  |
| pls_0 | gex_fraction_of_transcriptomic_reads_in_cells_mean | donor_technical.csv | 74 | -0.19686042206590149 | False |  |
| pls_0 | rin | covariates.csv (primary export) | 74 | -0.11312147045596636 | False |  |
| pls_0 | pmi | covariates.csv (primary export) | 74 | -0.0145176589153464 | False |  |
| pls_1 | tmm_rna | manifest.json (primary export, TMM factor) | 74 | 0.20873750462791554 | False |  |
| pls_1 | tmm_atac | manifest.json (primary export, TMM factor) | 74 | -0.37131432802665676 | True |  |
| pls_1 | fraction_mitochondrial_umis_mean | donor_technical.csv | 74 | -0.029692706405035168 | False |  |
| pls_1 | doublet_score_mean | donor_technical.csv | 74 | -0.1146094039244724 | False |  |
| pls_1 | gex_reads_mapped_confidently_to_intronic_regions_mean | donor_technical.csv | 74 | 0.13405504669043325 | False |  |
| pls_1 | gex_fraction_of_transcriptomic_reads_in_cells_mean | donor_technical.csv | 74 | 0.1605479452054794 | False |  |
| pls_1 | rin | covariates.csv (primary export) | 74 | 0.202069281126135 | False |  |
| pls_1 | pmi | covariates.csv (primary export) | 74 | 0.06529983724372136 | False |  |
| all | all |  | 74 |  | True | False |

### Variant `exclude_families` (pasted from `results/motco_sensitivity/exclude_families/comparison/cohort_50/sensitivity.md`)

#### MOTCO preprocessing sensitivity: cohort_50, variant `exclude_families`

RNA genes of the ribosomal, mitochondrial and hemoglobin sets (TECHNICAL_GENE_SETS) removed before the abundance filter; the ATAC block and every other export setting equal to the primary. Export setting(s) changed: `excluded_gene_sets` = ['ribosomal', 'mitochondrial', 'hemoglobin'] (primary: `excluded_gene_sets` = []); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 1, which stays the primary result; the `exclude_families` row (CV-selected rank 1) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

##### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 18.7091 | 0.213 | no | 11.6918 | 9.218 | 30.0679 |
| exclude_families | cv | 1 | 19.0327 | 0.202 | no | 11.7595 | 9.29099 | 29.9742 |

##### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0 | 1 | no | 1.98198 | 18.7931 | 0 |
| exclude_families | cv | 1 | 0 | 1 | no | 1.98198 | 18.7931 | 0 |

##### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0.56471 | 0.377 | no | 0.495877 | 0.343673 | 1.159 |
| exclude_families | cv | 1 | 0.563235 | 0.377 | no | 0.495428 | 0.34349 | 1.1612 |

##### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after gene-set removal | after abundance | kept (after selection) | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 35329 | 20238 | 5000 | 4977 | 23 |
| rna | exclude_families | 36601 | 35329 | 35249 | 20166 | 5000 | 4977 | 23 |
| atac | primary | 218882 | 212816 | 212816 | 124818 | 5000 | 5000 | 0 |
| atac | exclude_families | 218882 | 212816 | 212816 | 124818 | 5000 | 5000 | 0 |

##### Gene sets removed

RNA genes matched by the excluded sets of `TECHNICAL_GENE_SETS` are removed at step 2 after the chromosome removal (status `excluded_gene_set` in `features_rna.csv`); a matched gene on an excluded chromosome keeps `sex_chromosome` and is counted under `already on an excluded chromosome`. ATAC names match no set.

| omic | set | matched | removed here | already on an excluded chromosome |
|---|---|---:|---:|---:|
| rna | ribosomal | 81 | 77 | 4 |
| rna | mitochondrial | 13 | 0 | 13 |
| rna | hemoglobin | 3 | 3 | 0 |
| atac | ribosomal | 0 | 0 | 0 |
| atac | mitochondrial | 0 | 0 | 0 |
| atac | hemoglobin | 0 | 0 | 0 |

##### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 1 | 0.6616 | 151.6 | `results/motco_analysis/cohort_50` |
| exclude_families | cv | 1 | 0.6674 | 109.2 | `motco_analysis/cohort_50` |

##### Findings

- exclude_families export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

### Variant `rna_trend_ratio` (pasted from `results/motco_sensitivity/rna_trend_ratio/comparison/cohort_50/sensitivity.md`)

#### MOTCO preprocessing sensitivity: cohort_50, variant `rna_trend_ratio`

RNA features ranked at the selection step by the ratio of their log2-CPM variance to a binned-median mean–variance trend instead of by raw variance; the ATAC block and every other export setting equal to the primary. Export setting(s) changed: `selection_rule` = {'rna': 'trend_ratio', 'atac': 'variance'} (primary: `selection_rule` = {'rna': 'variance', 'atac': 'variance'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 1, which stays the primary result; the `rna_trend_ratio` row (CV-selected rank 1) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

##### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 18.7091 | 0.213 | no | 11.6918 | 9.218 | 30.0679 |
| rna_trend_ratio | cv | 1 | 18.0822 | 0.211 | no | 11.5123 | 8.89768 | 28.9267 |

##### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0 | 1 | no | 1.98198 | 18.7931 | 0 |
| rna_trend_ratio | cv | 1 | 0 | 1 | no | 1.44144 | 16.0512 | 0 |

##### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0.56471 | 0.377 | no | 0.495877 | 0.343673 | 1.159 |
| rna_trend_ratio | cv | 1 | 0.559086 | 0.353 | no | 0.480897 | 0.339893 | 1.14808 |

##### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after abundance | kept (after selection) | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 20238 | 5000 | 2572 | 2428 |
| rna | rna_trend_ratio | 36601 | 35329 | 20238 | 5000 | 2572 | 2428 |
| atac | primary | 218882 | 212816 | 124818 | 5000 | 5000 | 0 |
| atac | rna_trend_ratio | 218882 | 212816 | 124818 | 5000 | 5000 | 0 |

##### Selection rule

`variance`: population variance (ddof 0) of log2-CPM, the primary's rule. `trend_ratio`: variance over a mean–variance trend, the median variance of 20 equal-count bins of mean log2-CPM placed at each bin's median mean and linearly interpolated (flat beyond the end bins; the bins are in the variant export's manifest). `mean`: mean log2-CPM. Ties by name under every rule; the quantiles are of the kept features' mean log2-CPM.

| omic | export | rule | kept mean q05 | q25 | q50 | q75 | q95 |
|---|---|---|---:|---:|---:|---:|---:|
| rna | primary | variance | -0.662 | -0.096 | 0.818 | 2.229 | 4.957 |
| rna | rna_trend_ratio | trend_ratio | 0.100 | 2.105 | 3.813 | 5.342 | 7.526 |

##### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 1 | 0.6616 | 151.6 | `results/motco_analysis/cohort_50` |
| rna_trend_ratio | cv | 1 | 0.6716 | 118.4 | `motco_analysis/cohort_50` |

##### Findings

- rna_trend_ratio export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

### Variant `atac_mean_accessibility` (pasted from `results/motco_sensitivity/atac_mean_accessibility/comparison/cohort_50/sensitivity.md`)

#### MOTCO preprocessing sensitivity: cohort_50, variant `atac_mean_accessibility`

ATAC peaks ranked at the selection step by mean log2-CPM (total accessibility) instead of by variance; the RNA block and every other export setting equal to the primary. Export setting(s) changed: `selection_rule` = {'rna': 'variance', 'atac': 'mean'} (primary: `selection_rule` = {'rna': 'variance', 'atac': 'variance'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 1, which stays the primary result; the `atac_mean_accessibility` row (CV-selected rank 1) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

##### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 18.7091 | 0.213 | no | 11.6918 | 9.218 | 30.0679 |
| atac_mean_accessibility | cv | 1 | 9.1382 | 0.665 | no | 18.2862 | 15.3147 | 47.1689 |

##### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0 | 1 | no | 1.98198 | 18.7931 | 0 |
| atac_mean_accessibility | cv | 1 | 0 | 1 | no | 14.5946 | 49.1573 | 180 |

##### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0.56471 | 0.377 | no | 0.495877 | 0.343673 | 1.159 |
| atac_mean_accessibility | cv | 1 | 0.645472 | 0.465 | no | 0.65144 | 0.407746 | 1.33466 |

##### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after abundance | kept (after selection) | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 20238 | 5000 | 5000 | 0 |
| rna | atac_mean_accessibility | 36601 | 35329 | 20238 | 5000 | 5000 | 0 |
| atac | primary | 218882 | 212816 | 124818 | 5000 | 0 | 5000 |
| atac | atac_mean_accessibility | 218882 | 212816 | 124818 | 5000 | 0 | 5000 |

##### Selection rule

`variance`: population variance (ddof 0) of log2-CPM, the primary's rule. `trend_ratio`: variance over a mean–variance trend, the median variance of 20 equal-count bins of mean log2-CPM placed at each bin's median mean and linearly interpolated (flat beyond the end bins; the bins are in the variant export's manifest). `mean`: mean log2-CPM. Ties by name under every rule; the quantiles are of the kept features' mean log2-CPM.

| omic | export | rule | kept mean q05 | q25 | q50 | q75 | q95 |
|---|---|---|---:|---:|---:|---:|---:|
| atac | primary | variance | -0.432 | -0.057 | 0.271 | 0.698 | 1.479 |
| atac | atac_mean_accessibility | mean | 5.443 | 5.621 | 5.883 | 6.281 | 7.163 |

##### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 1 | 0.6616 | 151.6 | `results/motco_analysis/cohort_50` |
| atac_mean_accessibility | cv | 1 | 0.6820 | 106.3 | `motco_analysis/cohort_50` |

##### Findings

- atac_mean_accessibility export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

### Variant `residualized` (pasted from `results/motco_sensitivity/residualized/comparison/cohort_50/sensitivity.md`)

#### MOTCO preprocessing sensitivity: cohort_50, variant `residualized`

after log2-CPM and before the selection step, each abundance-passing feature replaced by its mean plus the residual of its regression on an intercept and RESIDUALIZATION_COVARIATES (the export's log10 library sizes and the technical table's gene-set fractions); the selection ranks the residual variance; every other export setting equal to the primary. Export setting(s) changed: `residualize` = ['log10_lib_rna', 'log10_lib_atac', 'ribosomal_pseudobulk', 'mitochondrial_pseudobulk'] (primary: `residualize` = None); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 1, which stays the primary result; the `residualized` row (CV-selected rank 2) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

##### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 18.7091 | 0.213 | no | 11.6918 | 9.218 | 30.0679 |
| residualized | cv | 2 | 8.94403 | 0.444 | no | 9.36668 | 7.40074 | 23.7553 |

##### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0 | 1 | no | 1.98198 | 18.7931 | 0 |
| residualized | cv | 2 | 30.0936 | 0.104 | no | 14.4779 | 11.3133 | 36.2531 |

##### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0.56471 | 0.377 | no | 0.495877 | 0.343673 | 1.159 |
| residualized | cv | 2 | 0.317151 | 0.338 | no | 0.264388 | 0.142759 | 0.517872 |

##### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after abundance | kept (after selection) | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 20238 | 5000 | 4218 | 782 |
| rna | residualized | 36601 | 35329 | 20238 | 5000 | 4218 | 782 |
| atac | primary | 218882 | 212816 | 124818 | 5000 | 3423 | 1577 |
| atac | residualized | 218882 | 212816 | 124818 | 5000 | 3423 | 1577 |

##### Residualization

After log2-CPM and before the selection step, each abundance-passing feature was replaced by its mean plus the residual of its OLS fit on an intercept and 4 covariates: `log10_lib_rna`, `log10_lib_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` (`RESIDUALIZATION_COVARIATES`: the export's own log10 library sizes after step 2 and the technical table's gene-set pseudobulk fractions; excluded as cell state under the technical-axis report's rule: `hemoglobin_pseudobulk` (reports/motco-astrocyte-technical-axis-2026-10-19.md, section 2, 'Classification under the spec's rule': median within-donor supertype share 0.00394 above the between-donor share 0.00285 (both negligible; median pseudobulk fraction 6.7e-07 of UMIs); not a residualization candidate)). The design (`residualization.csv` beside the variant blocks, sha256 `cbff6f8e8dcd…`) holds no sex, stage or analysis covariate. Limitations: the fit spends 5 degrees of freedom per feature (intercept plus covariates) before the permutation null, which does not see them; it removes whatever part of the stage and sex signal the covariates carry; the analysis covariates were left unchanged, so the library-size proxies are adjusted twice (log10 library sizes here, log10 nuclei counts in the RRPP model); and no Type I check exists for the residualized blocks.

| omic | R² q05 | q25 | q50 | q75 | q95 |
|---|---:|---:|---:|---:|---:|
| rna | 0.043 | 0.115 | 0.201 | 0.324 | 0.550 |
| atac | 0.032 | 0.080 | 0.131 | 0.202 | 0.378 |

##### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 1 | 0.6616 | 151.6 | `results/motco_analysis/cohort_50` |
| residualized | cv | 2 | 0.4791 | 105.9 | `motco_analysis/cohort_50` |

##### Findings

- CV-selected rank differs: 1 in the primary analysis, 2 in the residualized analysis; reported beside angle, not read as better or worse
- residualized export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

### Variant `revised_candidate` (pasted from `results/motco_sensitivity/revised_candidate/comparison/cohort_50/sensitivity.md`)

#### MOTCO preprocessing sensitivity: cohort_50, variant `revised_candidate`

the four settings together in the primary's order: gene sets removed before the abundance filter, residualization after log2-CPM, trend-ratio ranking for RNA and mean-accessibility ranking for ATAC at the selection step; the adoption rule ADOPTION_RULE is evaluated for this variant. Export setting(s) changed: `excluded_gene_sets` = ['ribosomal', 'mitochondrial', 'hemoglobin'], `residualize` = ['log10_lib_rna', 'log10_lib_atac', 'ribosomal_pseudobulk', 'mitochondrial_pseudobulk'], `selection_rule` = {'rna': 'trend_ratio', 'atac': 'mean'} (primary: `excluded_gene_sets` = [], `residualize` = None, `selection_rule` = {'rna': 'variance', 'atac': 'variance'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 1, which stays the primary result; the `revised_candidate` row (CV-selected rank 2) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

##### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 18.7091 | 0.213 | no | 11.6918 | 9.218 | 30.0679 |
| revised_candidate | cv | 2 | 6.33891 | 0.585 | no | 9.67538 | 7.57586 | 25.2956 |

##### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0 | 1 | no | 1.98198 | 18.7931 | 0 |
| revised_candidate | cv | 2 | 2.15631 | 0.902 | no | 12.4302 | 9.75934 | 30.5331 |

##### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0.56471 | 0.377 | no | 0.495877 | 0.343673 | 1.159 |
| revised_candidate | cv | 2 | 0.230511 | 0.564 | no | 0.276132 | 0.155841 | 0.565357 |

##### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after gene-set removal | after abundance | kept (after selection) | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 35329 | 20238 | 5000 | 2300 | 2700 |
| rna | revised_candidate | 36601 | 35329 | 35249 | 20166 | 5000 | 2300 | 2700 |
| atac | primary | 218882 | 212816 | 212816 | 124818 | 5000 | 0 | 5000 |
| atac | revised_candidate | 218882 | 212816 | 212816 | 124818 | 5000 | 0 | 5000 |

##### Gene sets removed

RNA genes matched by the excluded sets of `TECHNICAL_GENE_SETS` are removed at step 2 after the chromosome removal (status `excluded_gene_set` in `features_rna.csv`); a matched gene on an excluded chromosome keeps `sex_chromosome` and is counted under `already on an excluded chromosome`. ATAC names match no set.

| omic | set | matched | removed here | already on an excluded chromosome |
|---|---|---:|---:|---:|
| rna | ribosomal | 81 | 77 | 4 |
| rna | mitochondrial | 13 | 0 | 13 |
| rna | hemoglobin | 3 | 3 | 0 |
| atac | ribosomal | 0 | 0 | 0 |
| atac | mitochondrial | 0 | 0 | 0 |
| atac | hemoglobin | 0 | 0 | 0 |

##### Selection rule

`variance`: population variance (ddof 0) of log2-CPM, the primary's rule. `trend_ratio`: variance over a mean–variance trend, the median variance of 20 equal-count bins of mean log2-CPM placed at each bin's median mean and linearly interpolated (flat beyond the end bins; the bins are in the variant export's manifest). `mean`: mean log2-CPM. Ties by name under every rule; the quantiles are of the kept features' mean log2-CPM.

| omic | export | rule | kept mean q05 | q25 | q50 | q75 | q95 |
|---|---|---|---:|---:|---:|---:|---:|
| rna | primary | variance | -0.662 | -0.096 | 0.818 | 2.229 | 4.957 |
| rna | revised_candidate | trend_ratio | 0.148 | 2.413 | 3.980 | 5.501 | 7.651 |
| atac | primary | variance | -0.432 | -0.057 | 0.271 | 0.698 | 1.479 |
| atac | revised_candidate | mean | 5.443 | 5.621 | 5.883 | 6.281 | 7.163 |

##### Residualization

After log2-CPM and before the selection step, each abundance-passing feature was replaced by its mean plus the residual of its OLS fit on an intercept and 4 covariates: `log10_lib_rna`, `log10_lib_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` (`RESIDUALIZATION_COVARIATES`: the export's own log10 library sizes after step 2 and the technical table's gene-set pseudobulk fractions; excluded as cell state under the technical-axis report's rule: `hemoglobin_pseudobulk` (reports/motco-astrocyte-technical-axis-2026-10-19.md, section 2, 'Classification under the spec's rule': median within-donor supertype share 0.00394 above the between-donor share 0.00285 (both negligible; median pseudobulk fraction 6.7e-07 of UMIs); not a residualization candidate)). The design (`residualization.csv` beside the variant blocks, sha256 `ca8f17114698…`) holds no sex, stage or analysis covariate. Limitations: the fit spends 5 degrees of freedom per feature (intercept plus covariates) before the permutation null, which does not see them; it removes whatever part of the stage and sex signal the covariates carry; the analysis covariates were left unchanged, so the library-size proxies are adjusted twice (log10 library sizes here, log10 nuclei counts in the RRPP model); and no Type I check exists for the residualized blocks.

| omic | R² q05 | q25 | q50 | q75 | q95 |
|---|---:|---:|---:|---:|---:|
| rna | 0.043 | 0.114 | 0.199 | 0.321 | 0.544 |
| atac | 0.032 | 0.080 | 0.131 | 0.202 | 0.378 |

##### Adoption check

`ADOPTION_RULE`: holds when no latent component of the candidate's CV-selected analysis reaches |Spearman| >= DIAGNOSTIC_CORRELATION_THRESHOLD against any variable of ADOPTION_CHECK_VARIABLES, in both cohorts; p-values, rejection statuses and CV AUROC are not read. If it holds, a later change may propose the candidate as a revised primary with the current primary kept as the reference; if not, the variants stay sensitivity analyses. Threshold 0.3 (`DIAGNOSTIC_CORRELATION_THRESHOLD`). In this cohort the rule does not hold: flagged pls_0 × `tmm_atac` (0.37), pls_0 × `gex_reads_mapped_confidently_to_intronic_regions_mean` (-0.33), pls_1 × `tmm_atac` (-0.37). The rule needs both cohorts; whether the candidate becomes a revised primary is a later change with the current primary kept as the reference, and the candidate's statistics above are sensitivity, not the result.

| component | variable | source | n | Spearman | flagged |
|---|---|---|---:|---:|:---:|
| pls_0 | `tmm_rna` | manifest.json (primary export, TMM factor) | 74 | -0.287 | no |
| pls_0 | `tmm_atac` | manifest.json (primary export, TMM factor) | 74 | 0.375 | yes |
| pls_0 | `fraction_mitochondrial_umis_mean` | donor_technical.csv | 74 | 0.204 | no |
| pls_0 | `doublet_score_mean` | donor_technical.csv | 74 | 0.034 | no |
| pls_0 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 74 | -0.325 | yes |
| pls_0 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 74 | -0.197 | no |
| pls_0 | `rin` | covariates.csv (primary export) | 74 | -0.113 | no |
| pls_0 | `pmi` | covariates.csv (primary export) | 74 | -0.015 | no |
| pls_1 | `tmm_rna` | manifest.json (primary export, TMM factor) | 74 | 0.209 | no |
| pls_1 | `tmm_atac` | manifest.json (primary export, TMM factor) | 74 | -0.371 | yes |
| pls_1 | `fraction_mitochondrial_umis_mean` | donor_technical.csv | 74 | -0.030 | no |
| pls_1 | `doublet_score_mean` | donor_technical.csv | 74 | -0.115 | no |
| pls_1 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 74 | 0.134 | no |
| pls_1 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 74 | 0.161 | no |
| pls_1 | `rin` | covariates.csv (primary export) | 74 | 0.202 | no |
| pls_1 | `pmi` | covariates.csv (primary export) | 74 | 0.065 | no |
| all | all | | 74 | | yes → holds = false |

##### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 1 | 0.6616 | 151.6 | `results/motco_analysis/cohort_50` |
| revised_candidate | cv | 2 | 0.5514 | 100.6 | `motco_analysis/cohort_50` |

##### Findings

- adoption check does not hold in cohort_50: pls_0 x tmm_atac (0.37), pls_0 x gex_reads_mapped_confidently_to_intronic_regions_mean (-0.33), pls_1 x tmm_atac (-0.37); the variants stay sensitivity analyses
- CV-selected rank differs: 1 in the primary analysis, 2 in the revised_candidate analysis; reported beside angle, not read as better or worse
- revised_candidate export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

## 4. Reading

- **Statistics.** No statistic changed rejection status in any variant or cohort: nothing is rejected in the primary analysis and nothing is rejected in any variant analysis, so the non-rejection rules above apply to every row (uninformative `angle` and `shape`, no evidence about `delta`).
- **Ranks.** The CV-selected rank differs in cohort_50 `residualized` (1 in the primary analysis, 2 in the variant analysis); cohort_50 `revised_candidate` (1 in the primary analysis, 2 in the variant analysis); reported beside `angle`, not read as better or worse.
- **Attribution of the first axis (ribosomal fraction / RNA TMM factor / stage).** cohort_30: primary `pls_0` rho(ribosomal) -0.81, rho(tmm_rna) -0.73, rho(stage) 0.49; `exclude_families` -0.81 / -0.73 / 0.49; `rna_trend_ratio` -0.81 / -0.72 / 0.50; `atac_mean_accessibility` -0.70 / -0.64 / 0.36; `residualized` -0.05 / -0.08 / 0.75; `revised_candidate` -0.06 / -0.19 / 0.61 — cohort_50: primary `pls_0` rho(ribosomal) -0.82, rho(tmm_rna) -0.72, rho(stage) 0.50; `exclude_families` -0.81 / -0.70 / 0.50; `rna_trend_ratio` -0.82 / -0.72 / 0.51; `atac_mean_accessibility` -0.73 / -0.62 / 0.36; `residualized` -0.03 / -0.16 / 0.68; `revised_candidate` -0.08 / -0.29 / 0.52. The single-setting variants share every other setting with the primary and the candidate, so the change of each correlation is attributable to the one setting; this is a measurement-space description, and the statistics of each row are read under the rules above regardless of it.

- **What each setting did to the first axis (cohort_30; cohort_50 in parentheses).** `exclude_families` removes 77 autosomal ribosomal genes and 3 hemoglobin genes (0 mitochondrial beyond the primary's chrM exclusion; the primary had kept 9 of the 81 ribosomal genes) and leaves `pls_0` where it was: rho(ribosomal fraction) -0.81 → -0.81 (-0.81), rho(RNA TMM) -0.73 → -0.73, rho(stage) 0.49 → 0.49, with 4990 of 5000 RNA features shared. The axis is therefore not carried by the ribosomal genes themselves but by the composition effect across the block. `rna_trend_ratio` replaces 2464 of 5000 RNA features, raises the kept median mean log2-CPM from 0.77 to 3.82 and keeps 58 of the 81 ribosomal genes (their variance exceeds the trend at their abundance), yet the axis is again unchanged (ribosomal -0.81, TMM -0.72, stage 0.50). `atac_mean_accessibility` replaces all 5000 ATAC peaks (kept median mean 0.37 → 5.88), rescales the axis (SD 35.0 → 54.4) and flips its ATAC library-size correlation (0.21 → -0.18), but it stays the composition axis (ribosomal -0.70, TMM -0.64; stage 0.36). `residualized` is the one setting that removes it: ribosomal -0.05 (-0.03), TMM -0.08 (-0.16), mitochondrial 0.09; the first axis's stage correlation rises to 0.75 (0.68), its SD falls from 35.0 to 13.9, 748 RNA and 1553 ATAC features change, and no ribosomal gene stays in the kept set (0 of 81). `revised_candidate` inherits the residualized picture (ribosomal -0.06 (-0.08), TMM -0.19 (-0.29), stage 0.61 (0.52)).
- **What the variant analyses show, under the rules.** Nothing is rejected in any row, so every `angle` and `shape` row is uninformative and no `delta` row is evidence about a size difference. In cohort_50 the residualized blocks select rank 2 (`residualized` and `revised_candidate`; angle 30.0936° at p 0.104 and 2.15631° at p 0.902, both non-rejections and therefore uninformative), reported beside `angle` and not read as better or worse. The held-out `cv_mean_auroc` falls from 0.6313 / 0.6616 in the primary to 0.4078 / 0.4791 with residualization and 0.4822 / 0.5514 for the candidate: once the technical set is removed, the stage discrimination of the latent space in held-out folds is at or below chance, so the higher in-sample stage correlation of the residualized first axis does not generalize across folds. The technical-axis report anticipated this: stage and the fractions are not independent in these donors, and the residualization removes the stage-related variance that travels with them together with the preparation effect. This is a measurement-space observation reported beside the adoption check (design D7 rejected a CV-AUROC floor as a rule), and it says nothing about the primary analysis's statistics.
- **Adoption outcome.** The check holds in cohort_30 (largest |rho| 0.22, below 0.3) and fails in cohort_50 on the ATAC TMM factor (`pls_0` 0.37, `pls_1` -0.37) and the intronic-read fraction (`pls_0` -0.33). The rule needs both cohorts, so it is not met: the five variants stay sensitivity analyses, no revised-primary change follows from this report, and the primary preprocessing stays the primary. The ATAC TMM factor was kept out of the residualization set by design (D3: a normalization by-product, used as a check variable), and the candidate's cohort_50 axes still track it although the ATAC log10 library size was a regressor; whether a different technical set or a MOTCO-side residualization (the follow-up proposed in the MOTCO repository) would meet the rule is a question for a new proposal, written before any run.


## 5. Limitations

- **Residualization (`residualized`, `revised_candidate`).** Covariates ['log10_lib_rna', 'log10_lib_atac', 'ribosomal_pseudobulk', 'mitochondrial_pseudobulk']: the fit spends 5 degrees of freedom per feature (intercept plus 4 covariates) before the permutation null, which does not see them; it removes whatever part of the stage and sex signal the technical set carries (the technical-axis report found `pls_0` correlated with stage as well as with the fractions, and this report cannot separate preparation noise from stage-related composition); the analysis covariates were left unchanged, so the library-size proxies are adjusted twice (log10 library sizes in the design, log10 nuclei counts in the RRPP model); and there is no Type I check for the residualized blocks.
- **Gene-set exclusion (`exclude_families`, `revised_candidate`).** Removing the ribosomal, mitochondrial and hemoglobin genes removes the measured proxies of the composition effect from the block, not the effect on the remaining genes; the `MT-` genes were already removed by the chrM exclusion of the primary, and the four chrX ribosomal genes by the chrX exclusion.
- **Selection rules (`rna_trend_ratio`, `atac_mean_accessibility`, `revised_candidate`).** The trend is a 20-bin median fit with flat extrapolation beyond the end bins, so the ratio can favour the high-mean tail where the trend is flat; mean accessibility ranks peaks by abundance alone, so a low-abundance peak cannot enter the block whatever its stage variation. Both rules change which features enter the blocks; MOTCO z-scores features without block weights, so N = 5,000 per block stays equal and the comparison with the primary is one of feature content only.
- **No Type I or power claim.** MOTCO's forced-rank and adjusted operating characteristics were measured on the primary-style blocks only; no variant has a Type I check, and the adoption rule is a technical-correlation criterion, not a validation of the candidate's inference.
- **Descriptive correlations at n = 80 and 74.** No p-value is attached to any correlation of the attribution or adoption tables; the two cohorts share 74 donors and are not independent replicates; `cohort_50` is a sensitivity cohort.
- **Double role of the technical table.** The gene-set fractions are regressors of the residualization and the QC means are check variables; both come from the same astrocyte nuclei and the same libraries, so the check variables are not independent of the regressors.

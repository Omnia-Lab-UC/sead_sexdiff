# MOTCO astrocyte trend-normalization sensitivity: sex × ADNC stage — 2026-10-20

This report reads the two trend-normalization variants of the registry `PREPROCESSING_VARIANTS` added by the change `mtg-astrocyte-trend-normalization-sensitivity`: `trend_normalized` (limma's fast cyclic loess of each donor's log2-CPM against the per-feature mean over donors, after log2-CPM and before the selection, in both omics; every other setting the primary's) and `revised_candidate_trend` (the archived candidate's four settings plus the same normalization). Each is the primary export with the registered settings changed and the covariate-adjusted MOTCO analysis rerun unchanged on the variant blocks, beside the primary analysis of `reports/motco-astrocyte-2026-10-19.md`. The variants were declared, with the adoption check on both and the precedence rule between them, in the change proposal before any run; they follow from a label-free check that located a per-donor abundance tilt in the export's log2-CPM (proposal). The feature-content variants and the adoption rule are in `reports/motco-astrocyte-feature-selection-2026-10-19.md`, the diagnostic of the first latent axis in `reports/motco-astrocyte-technical-axis-2026-10-19.md`. Every number below is pasted from the generated files (`sensitivity.md`, `sensitivity.csv`, `sensitivity_manifest.json`, `adoption_check.csv`, `normalization_trend.csv`, `latent_scores.csv`, and `axis_correlations.csv` / `geometry_shares.csv` of scratch `axis-diagnostics` runs) or, for the six archived variants, from the feature-selection report's attribution tables. The primary analysis stays primary, both variants are sensitivity to it, and the two decisions this report carries, the adoption check of each variant and the precedence between them, are predeclared rules over generated numbers that read no p-value.

**Reading rules (binding).** The primary analysis remains the primary result; no variant is chosen from the results. `angle` and `shape` non-rejections are uninformative in every row (they are not evidence that the sexes agree in direction or shape); a `delta` rejection is a difference, not specifically of size. A rejection status that differs between a variant and the primary is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics (`../MOTCO/docs/roadmap.md`, Phase 6) were measured on neither residualized nor family-excluded nor trend-selected nor trend-normalized blocks, so neither variant carries a Type I or power claim of its own. A latent axis correlated with a technical variable is a measurement-space finding about the analysis in question, not a reinterpretation of its statistics.

## 1. Provenance

- Software: MOTCO 0.8.0 (git dependency pinned in `uv.lock`), scikit-learn 1.8.0, numpy 2.3.5, scipy 1.15.3, pandas 2.3.3, Python 3.12.15; lowess: statsmodels.nonparametric.smoothers_lowess.lowess (statsmodels 0.14.5); frac = span, it = 3, delta = 0.01 x range(A), as R's lowess in limma's loessFit; sead-sexdiff revision `18a6bc612231660cfb94f7a4dca19e7e2f8810de` (dirty: True; the implementation of this change was uncommitted at run time; a rerun of `trend_normalized` for cohort_50 reproduced its export, analysis and comparison files byte for byte, task 4.4).
- Machine: `nixos-dell-laptop`, by decision of 2026-10-20 (`nixos-home`, where the archived variant trees live, was not reachable). Its kept pseudobulk, covariates, annotation, primary export and primary analysis files have the sha256 values recorded in the feature-selection report (below), and the technical table was synced from the cluster with `scripts/sync_pseudobulk.sh`. The two variants ran one call each (`--variant trend_normalized`, then `--variant revised_candidate_trend`), since `--variant all` would have rerun the six archived variants, which are not on this machine; their numbers below come from the feature-selection report.
- Settings shared by every row: seed 20261019, 999 RRPP permutations (the same permutation indices in every analysis), α = 0.05 (p <= alpha), `n_jobs = 1`; covariates age_at_death, pmi, rin, years_education, log10_n_rna, log10_n_atac entered additively in MOTCO's regression on the latent scores (unchanged in both variants).
- Inputs: pseudobulk manifest sha256 `2e93cc0b526856fc043f908a424a11266187146636bc75afed7edd187989a25c` (created 2026-10-16T18:46:13+00:00), covariates manifest `b8bf9091b507a1098af91c9390f487ccb0c78aa87fba07edc0bfd1d6b83184ae`, annotation GENCODE v32 basic sha256 `a61ab50fc8b5d63f76d50e1f89bc0f5bba0b64fb29d2ac8cfa1b1ef8d4ef0c91`, technical table `eac1b9a15977dfab683738a4bd0b0ffc207b896a102032409b60fccbfb96247b` and its manifest `c2bdb34499644580f22ae7f6c48a7f846c4dabad5b2b04182f8e9d5a98a106e7` (recorded pseudobulk manifest sha256 `2e93cc0b526856fc043f908a424a11266187146636bc75afed7edd187989a25c`, equal to the export's).
- cohort_30: primary export `results/motco_input/cohort_30` (manifest sha256 `18f2289e0d4c84ed18f4a437ae569052544e9d150ed2726385273694b4a3b0b6`, k = 9, 80 donors), primary analysis `results/motco_analysis/cohort_30` (`result.json` sha256 `ca92876a1a61c40791524390e8603125bcbfacf6f081d34e2ae6a9b49fd15de6`, CV-selected rank 1, cv_mean_auroc 0.6313, MOTCO 0.8.0, seed 20261019, 999 permutations); both read and checked, never rerun.
- cohort_50: primary export `results/motco_input/cohort_50` (manifest sha256 `dd14a6e7d5583fb527e8a99c7db43d12cf237f5451522e8875c5a0a3401c5cb3`, k = 9, 74 donors), primary analysis `results/motco_analysis/cohort_50` (`result.json` sha256 `1e7eefe3e63e7b9b51c3fed618c32d4a496ea3d382c7ff025c2a34f5d0ae048f`, CV-selected rank 1, cv_mean_auroc 0.6616, MOTCO 0.8.0, seed 20261019, 999 permutations); both read and checked, never rerun.
- Variants (registry entries; settings as the manifests record them):
  - `trend_normalized`: after log2-CPM and before the selection step, each donor's log2-CPM over the abundance-passing features corrected by limma's fast cyclic loess (a lowess curve of its deviation from the per-feature mean over donors against that mean, three iterations, adaptive span) in both omics; every other export setting equal to the primary; the adoption rule ADOPTION_RULE is evaluated for this variant. Settings: `{"normalization": {"rna": "cyclic_loess", "atac": "cyclic_loess"}}`.
  - `revised_candidate_trend`: the revised candidate's four settings plus the cyclic-loess trend normalization of both omics, in the primary's order: gene sets removed before the abundance filter, cyclic loess after log2-CPM, residualization after the normalization, trend-ratio ranking for RNA and mean-accessibility ranking for ATAC at the selection step; the adoption rule ADOPTION_RULE is evaluated for this variant. Settings: `{"excluded_gene_sets": ["ribosomal", "mitochondrial", "hemoglobin"], "normalization": {"rna": "cyclic_loess", "atac": "cyclic_loess"}, "residualize": ["log10_lib_rna", "log10_lib_atac", "ribosomal_pseudobulk", "mitochondrial_pseudobulk"], "selection_rule": {"rna": "trend_ratio", "atac": "mean"}}`.
- Trend normalization (`NORMALIZATION_RULES` = ['none', 'cyclic_loess']): limma's `normalizeCyclicLoess(method = "fast")`, `LOESS_ITERATIONS` = 3 cyclic iterations against the per-feature mean over donors, `LOESS_ROBUSTNESS_ITERATIONS` = 3, `LOESS_DELTA_FRACTION` = 0.01, span by `LOESS_SPAN_RULE` = {'small_n': 50, 'min_span': 0.3, 'power': 0.3333333333333333} (limma's `chooseLowessSpan` defaults) over the block's abundance-passing features; tilt accounting over `NORMALIZATION_TILT_BINS` = 20 equal-count bins of A.
- Cross-check against limma (3.68.4; `scripts/check_cyclic_loess_limma.R` on the debug files of a scratch `_export_cohort` call with the `trend_normalized` settings, whose `normalization_trend.csv` equals the kept one): tolerance 1e-06 log2-CPM; every recorded span equals limma's adaptive span; largest absolute differences cohort_30 rna 1.48e-08 (median 1.18e-10, limma corrections up to 1.22), cohort_30 atac 3.48e-09 (median 8.16e-11, limma corrections up to 1.13), cohort_50 rna 1.77e-08 (median 1.33e-10, limma corrections up to 1.24), cohort_50 atac 2.85e-09 (median 1.46e-10, limma corrections up to 1.13): PASS in all four blocks.
- Adoption rule `ADOPTION_RULE`: holds when no latent component of the candidate's CV-selected analysis reaches |Spearman| >= DIAGNOSTIC_CORRELATION_THRESHOLD against any variable of ADOPTION_CHECK_VARIABLES, in both cohorts; p-values, rejection statuses and CV AUROC are not read. If it holds, a later change may propose the candidate as a revised primary with the current primary kept as the reference; if not, the variants stay sensitivity analyses. Check variables `ADOPTION_CHECK_VARIABLES` = ['tmm_rna', 'tmm_atac', 'fraction_mitochondrial_umis_mean', 'doublet_score_mean', 'gex_reads_mapped_confidently_to_intronic_regions_mean', 'gex_fraction_of_transcriptomic_reads_in_cells_mean', 'rin', 'pmi']; threshold 0.3. Evaluated for both variants.
- Precedence rule `ADOPTION_PRECEDENCE` (candidates ['trend_normalized', 'revised_candidate_trend']): if trend_normalized meets ADOPTION_RULE in both cohorts it is the candidate a later revised-primary change proposes, whatever revised_candidate_trend gives (no gene removed, no degrees of freedom spent, no technical covariate read); if only revised_candidate_trend meets it in both cohorts, that one is; if neither does, the variants stay sensitivity analyses. Only the holds values of the adoption tables are read.
- Descriptive checks: `axis-diagnostics` run on the primary root and on each variant root into a scratch directory (sections 2–3), reading the export, the analysis and the technical table; the tilt statistics of `normalization_trend.csv` correlated with the latent scores; its outputs feed no export or analysis command.

## 2. Primary cohort: cohort_30 (n = 80)

### Overview (generated from `sensitivity.csv` and `sensitivity_manifest.json`)

| analysis | rank | cv_mean_auroc | delta (p) | angle (p) | shape (p) | rejected | rna kept shared / only primary / only variant | atac kept shared / only primary / only variant | adoption check |
|---|---:|---:|---:|---:|---:|---|---|---|---|
| primary | 1 | 0.6313 | 14.941 (0.3) | 0 (1) | 0.608832 (0.333) | none | — | — | — |
| trend_normalized | 1 | 0.6561 | 8.67659 (0.504) | 0 (1) | 0.591169 (0.293) | none | 4329 / 671 / 671 | 3360 / 1640 / 1640 | holds = false |
| revised_candidate_trend | 3 | 0.4946 | 14.3836 (0.166) | 16.7741 (0.435) | 0.0792388 (0.852) | none | 2238 / 2762 / 2762 | 0 / 5000 / 5000 | holds = false |

### Normalization: tilt before and after (generated from `sensitivity_manifest.json`, quantiles over donors of `normalization_trend.csv`)

Per donor, A is the per-feature mean over donors and M the donor's value minus A, of the log2-CPM before or after the cyclic loess; slope = OLS slope of M on A; the bins are 20 equal-count bins of A. The same rows are in each variant's pasted `sensitivity.md` below.

| variant | omic | features | span | statistic | before q05 | q25 | q50 | q75 | q95 | after q05 | q25 | q50 | q75 | q95 |
|---|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| trend_normalized | rna | 20387 | 0.394400 | slope of M on A | -0.0639 | -0.0342 | 0.0012 | 0.0294 | 0.0655 | -0.0140 | -0.0066 | -0.0026 | 0.0043 | 0.0188 |
| trend_normalized | rna | 20387 | 0.394400 | median M, lowest bin | -0.3761 | -0.1549 | 0.0267 | 0.1926 | 0.4588 | -0.1090 | -0.0214 | 0.0102 | 0.0419 | 0.0651 |
| trend_normalized | rna | 20387 | 0.394400 | median M, highest bin | -0.2846 | -0.1030 | 0.0048 | 0.1100 | 0.2457 | -0.0145 | -0.0086 | -0.0031 | 0.0022 | 0.0164 |
| trend_normalized | rna | 20387 | 0.394400 | largest absolute correction (per donor) | 0.1311 | 0.2392 | 0.3396 | 0.5065 | 0.9425 | | | | | |
| trend_normalized | atac | 128916 | 0.351049 | slope of M on A | -0.1610 | -0.0483 | 0.0200 | 0.0550 | 0.0826 | -0.0247 | -0.0043 | 0.0035 | 0.0072 | 0.0154 |
| trend_normalized | atac | 128916 | 0.351049 | median M, lowest bin | -0.9403 | -0.8906 | -0.0687 | 0.1940 | 0.6540 | -0.6540 | -0.4934 | -0.0842 | -0.0499 | 0.1628 |
| trend_normalized | atac | 128916 | 0.351049 | median M, highest bin | -0.4164 | -0.1617 | 0.0252 | 0.1658 | 0.3549 | -0.0118 | -0.0007 | 0.0085 | 0.0160 | 0.0297 |
| trend_normalized | atac | 128916 | 0.351049 | largest absolute correction (per donor) | 0.3428 | 0.4949 | 0.6157 | 0.7275 | 0.8875 | | | | | |
| revised_candidate_trend | rna | 20317 | 0.394508 | slope of M on A | -0.0639 | -0.0342 | 0.0013 | 0.0296 | 0.0660 | -0.0139 | -0.0065 | -0.0025 | 0.0043 | 0.0189 |
| revised_candidate_trend | rna | 20317 | 0.394508 | median M, lowest bin | -0.3774 | -0.1554 | 0.0277 | 0.1933 | 0.4569 | -0.0996 | -0.0237 | 0.0094 | 0.0385 | 0.0657 |
| revised_candidate_trend | rna | 20317 | 0.394508 | median M, highest bin | -0.2835 | -0.1029 | 0.0028 | 0.1104 | 0.2436 | -0.0141 | -0.0087 | -0.0033 | 0.0022 | 0.0161 |
| revised_candidate_trend | rna | 20317 | 0.394508 | largest absolute correction (per donor) | 0.1312 | 0.2384 | 0.3403 | 0.5049 | 0.9387 | | | | | |
| revised_candidate_trend | atac | 128916 | 0.351049 | slope of M on A | -0.1610 | -0.0483 | 0.0200 | 0.0550 | 0.0826 | -0.0247 | -0.0043 | 0.0035 | 0.0072 | 0.0154 |
| revised_candidate_trend | atac | 128916 | 0.351049 | median M, lowest bin | -0.9403 | -0.8906 | -0.0687 | 0.1940 | 0.6540 | -0.6540 | -0.4934 | -0.0842 | -0.0499 | 0.1628 |
| revised_candidate_trend | atac | 128916 | 0.351049 | median M, highest bin | -0.4164 | -0.1617 | 0.0252 | 0.1658 | 0.3549 | -0.0118 | -0.0007 | 0.0085 | 0.0160 | 0.0297 |
| revised_candidate_trend | atac | 128916 | 0.351049 | largest absolute correction (per donor) | 0.3428 | 0.4949 | 0.6157 | 0.7275 | 0.8875 | | | | | |

### Tilt and the first latent axis (generated from `normalization_trend.csv` and `latent_scores.csv`)

Spearman correlation over donors of each tilt statistic with `pls_0` of the primary analysis and of the variant's analysis. For `trend_normalized` the `before` values are computed on the primary pipeline's log2-CPM, so its primary-`pls_0` column of the `before` rows reproduces the label-free check of the proposal from a step output; for `revised_candidate_trend` they are computed after its gene-set removal. Descriptive: no p-value is attached.

| variant | omic | statistic | before × primary pls_0 | before × variant pls_0 | after × primary pls_0 | after × variant pls_0 |
|---|---|---|---:|---:|---:|---:|
| trend_normalized | rna | slope of M on A | 0.90 | 0.84 | -0.29 | -0.25 |
| trend_normalized | rna | median M, lowest bin | -0.89 | -0.84 | -0.21 | -0.21 |
| trend_normalized | rna | median M, highest bin | 0.78 | 0.75 | -0.03 | -0.05 |
| trend_normalized | rna | largest absolute correction | -0.18 | -0.19 | | |
| trend_normalized | atac | slope of M on A | -0.77 | -0.70 | 0.52 | 0.51 |
| trend_normalized | atac | median M, lowest bin | 0.65 | 0.57 | 0.28 | 0.23 |
| trend_normalized | atac | median M, highest bin | -0.59 | -0.60 | -0.58 | -0.59 |
| trend_normalized | atac | largest absolute correction | -0.27 | -0.18 | | |
| revised_candidate_trend | rna | slope of M on A | 0.90 | 0.14 | -0.29 | -0.20 |
| revised_candidate_trend | rna | median M, lowest bin | -0.89 | -0.20 | -0.19 | -0.02 |
| revised_candidate_trend | rna | median M, highest bin | 0.77 | -0.02 | -0.06 | -0.14 |
| revised_candidate_trend | rna | largest absolute correction | -0.17 | 0.04 | | |
| revised_candidate_trend | atac | slope of M on A | -0.77 | -0.22 | 0.52 | 0.08 |
| revised_candidate_trend | atac | median M, lowest bin | 0.65 | 0.24 | 0.28 | 0.21 |
| revised_candidate_trend | atac | median M, highest bin | -0.59 | -0.14 | -0.58 | -0.15 |
| revised_candidate_trend | atac | largest absolute correction | -0.27 | 0.10 | | |

Per-donor mean z-score over the selected block (each feature z-scored over donors, population SD) against `pls_0`:

| export | omic | × primary pls_0 | × pls_0 of its own analysis |
|---|---|---:|---:|
| primary | rna | -0.94 | -0.94 |
| primary | atac | 0.51 | 0.51 |
| trend_normalized | rna | -0.72 | -0.72 |
| trend_normalized | atac | -0.31 | -0.43 |
| revised_candidate_trend | rna | -0.31 | -0.32 |
| revised_candidate_trend | atac | -0.02 | -0.24 |

### Attribution: the first latent axis per analysis (generated from `axis_correlations.csv` and `geometry_shares.csv`; archived rows from the feature-selection report)

Spearman correlation of each latent component with stage, the astrocyte ribosomal-protein UMI fraction `ribosomal_pseudobulk`, the RNA TMM factor `tmm_rna`, the ATAC TSS enrichment `atac_tss_enrichment_score_mean`, the ATAC log10 library size `log10_lib_atac`, the mitochondrial fraction `mitochondrial_pseudobulk` and sex, per analysis; the component's score SD; and the share of the sex-difference geometry (Σ‖sex difference by stage‖²) along it. The archived rows (`pls_0` only) are pasted from `reports/motco-astrocyte-feature-selection-2026-10-19.md`, which did not tabulate the TSS enrichment (`n/r`); the primary row computed here equals its archived row. Descriptive: no p-value is attached.

| analysis | rank | component | stage | ribosomal fraction | RNA TMM factor | ATAC TSS enrichment | ATAC log10 library size | mitochondrial fraction | sex (male) | component SD | sex-difference share |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| primary | 1 | pls_0 | 0.49 | -0.81 | -0.73 | -0.77 | 0.21 | 0.73 | -0.23 | 35.0 | 1.000 |
| exclude_families (archived) | 1 | pls_0 | 0.49 | -0.81 | -0.73 | n/r | 0.21 | 0.73 | -0.23 | pls_0 35.0 | 1.000 |
| rna_trend_ratio (archived) | 1 | pls_0 | 0.50 | -0.81 | -0.72 | n/r | 0.20 | 0.71 | -0.23 | pls_0 34.3 | 1.000 |
| atac_mean_accessibility (archived) | 1 | pls_0 | 0.36 | -0.70 | -0.64 | n/r | -0.18 | 0.70 | -0.31 | pls_0 54.4 | 1.000 |
| residualized (archived) | 1 | pls_0 | 0.75 | -0.05 | -0.08 | n/r | 0.01 | 0.09 | 0.01 | pls_0 13.9 | 1.000 |
| revised_candidate (archived) | 1 | pls_0 | 0.61 | -0.06 | -0.19 | n/r | 0.03 | 0.13 | -0.00 | pls_0 16.6 | 1.000 |
| trend_normalized | 1 | pls_0 | 0.53 | -0.81 | -0.67 | -0.77 | 0.13 | 0.71 | -0.26 | 32.2 | 1.000 |
| revised_candidate_trend | 3 | pls_0 | 0.64 | -0.07 | -0.16 | -0.14 | 0.03 | 0.15 | -0.04 | 17.8 | 0.051 |
| revised_candidate_trend | 3 | pls_1 | -0.44 | -0.03 | -0.32 | -0.28 | 0.02 | 0.08 | -0.02 | 25.1 | 0.699 |
| revised_candidate_trend | 3 | pls_2 | -0.05 | -0.03 | -0.23 | -0.13 | 0.06 | 0.17 | -0.01 | 14.5 | 0.250 |

### Adoption check of `trend_normalized` (pasted from `adoption_check.csv`)

Rule: holds when no latent component of the candidate's CV-selected analysis reaches |Spearman| >= DIAGNOSTIC_CORRELATION_THRESHOLD against any variable of ADOPTION_CHECK_VARIABLES, in both cohorts; p-values, rejection statuses and CV AUROC are not read. If it holds, a later change may propose the candidate as a revised primary with the current primary kept as the reference; if not, the variants stay sensitivity analyses. Threshold 0.3. Outcome in cohort_30: **holds = false** (6 of 8 cells flagged).

| component | variable | source | n | spearman | flag | holds |
|---|---|---|---:|---:|:---:|:---:|
| pls_0 | tmm_rna | manifest.json (primary export, TMM factor) | 80 | -0.673230192217534 | True |  |
| pls_0 | tmm_atac | manifest.json (primary export, TMM factor) | 80 | 0.5054149085794655 | True |  |
| pls_0 | fraction_mitochondrial_umis_mean | donor_technical.csv | 80 | 0.7298406000937647 | True |  |
| pls_0 | doublet_score_mean | donor_technical.csv | 80 | 0.01809657759024848 | False |  |
| pls_0 | gex_reads_mapped_confidently_to_intronic_regions_mean | donor_technical.csv | 80 | -0.5649118900231207 | True |  |
| pls_0 | gex_fraction_of_transcriptomic_reads_in_cells_mean | donor_technical.csv | 80 | -0.6321143928738866 | True |  |
| pls_0 | rin | covariates.csv (primary export) | 80 | -0.334618905873116 | True |  |
| pls_0 | pmi | covariates.csv (primary export) | 80 | 0.031630285466701116 | False |  |
| all | all |  | 80 |  | True | False |

### Adoption check of `revised_candidate_trend` (pasted from `adoption_check.csv`)

Rule: holds when no latent component of the candidate's CV-selected analysis reaches |Spearman| >= DIAGNOSTIC_CORRELATION_THRESHOLD against any variable of ADOPTION_CHECK_VARIABLES, in both cohorts; p-values, rejection statuses and CV AUROC are not read. If it holds, a later change may propose the candidate as a revised primary with the current primary kept as the reference; if not, the variants stay sensitivity analyses. Threshold 0.3. Outcome in cohort_30: **holds = false** (1 of 24 cells flagged).

| component | variable | source | n | spearman | flag | holds |
|---|---|---|---:|---:|:---:|:---:|
| pls_0 | tmm_rna | manifest.json (primary export, TMM factor) | 80 | -0.15850914205344585 | False |  |
| pls_0 | tmm_atac | manifest.json (primary export, TMM factor) | 80 | 0.14174871073605252 | False |  |
| pls_0 | fraction_mitochondrial_umis_mean | donor_technical.csv | 80 | 0.15860290670417254 | False |  |
| pls_0 | doublet_score_mean | donor_technical.csv | 80 | 0.015822784810126583 | False |  |
| pls_0 | gex_reads_mapped_confidently_to_intronic_regions_mean | donor_technical.csv | 80 | -0.2125188262977556 | False |  |
| pls_0 | gex_fraction_of_transcriptomic_reads_in_cells_mean | donor_technical.csv | 80 | -0.12116736990154713 | False |  |
| pls_0 | rin | covariates.csv (primary export) | 80 | -0.020385436234484106 | False |  |
| pls_0 | pmi | covariates.csv (primary export) | 80 | -0.0018875003558706002 | False |  |
| pls_1 | tmm_rna | manifest.json (primary export, TMM factor) | 80 | -0.3225035161744023 | True |  |
| pls_1 | tmm_atac | manifest.json (primary export, TMM factor) | 80 | 0.1456399437412096 | False |  |
| pls_1 | fraction_mitochondrial_umis_mean | donor_technical.csv | 80 | 0.0890295358649789 | False |  |
| pls_1 | doublet_score_mean | donor_technical.csv | 80 | 0.12894983591186124 | False |  |
| pls_1 | gex_reads_mapped_confidently_to_intronic_regions_mean | donor_technical.csv | 80 | -0.11486237026902739 | False |  |
| pls_1 | gex_fraction_of_transcriptomic_reads_in_cells_mean | donor_technical.csv | 80 | -0.03900609470229723 | False |  |
| pls_1 | rin | covariates.csv (primary export) | 80 | -0.04762853790725067 | False |  |
| pls_1 | pmi | covariates.csv (primary export) | 80 | 0.08449201903577276 | False |  |
| pls_2 | tmm_rna | manifest.json (primary export, TMM factor) | 80 | -0.2342240975152368 | False |  |
| pls_2 | tmm_atac | manifest.json (primary export, TMM factor) | 80 | 0.1721050164088139 | False |  |
| pls_2 | fraction_mitochondrial_umis_mean | donor_technical.csv | 80 | 0.16727613689639007 | False |  |
| pls_2 | doublet_score_mean | donor_technical.csv | 80 | 0.017229254571026725 | False |  |
| pls_2 | gex_reads_mapped_confidently_to_intronic_regions_mean | donor_technical.csv | 80 | -0.06265859504675718 | False |  |
| pls_2 | gex_fraction_of_transcriptomic_reads_in_cells_mean | donor_technical.csv | 80 | 0.04631973745897796 | False |  |
| pls_2 | rin | covariates.csv (primary export) | 80 | 0.056408694169256764 | False |  |
| pls_2 | pmi | covariates.csv (primary export) | 80 | 0.03531149734088353 | False |  |
| all | all |  | 80 |  | True | False |

### Variant `trend_normalized` (pasted from `results/motco_sensitivity/trend_normalized/comparison/cohort_30/sensitivity.md`)

#### MOTCO preprocessing sensitivity: cohort_30, variant `trend_normalized`

after log2-CPM and before the selection step, each donor's log2-CPM over the abundance-passing features corrected by limma's fast cyclic loess (a lowess curve of its deviation from the per-feature mean over donors against that mean, three iterations, adaptive span) in both omics; every other export setting equal to the primary; the adoption rule ADOPTION_RULE is evaluated for this variant. Export setting(s) changed: `normalization` = {'rna': 'cyclic_loess', 'atac': 'cyclic_loess'} (primary: `normalization` = {'rna': 'none', 'atac': 'none'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 1, which stays the primary result; the `trend_normalized` row (CV-selected rank 1) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected, residualized nor trend-normalized), so the variant carries no Type I or power claim of its own.

##### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 14.941 | 0.3 | no | 11.2868 | 8.27982 | 26.5775 |
| trend_normalized | cv | 1 | 8.67659 | 0.504 | no | 10.2371 | 7.61107 | 24.4672 |

##### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0 | 1 | no | 3.06306 | 23.2919 | 0 |
| trend_normalized | cv | 1 | 0 | 1 | no | 0.900901 | 12.7087 | 0 |

##### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0.608832 | 0.333 | no | 0.489849 | 0.355498 | 1.18287 |
| trend_normalized | cv | 1 | 0.591169 | 0.293 | no | 0.429787 | 0.319469 | 1.06008 |

##### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after abundance | kept (after selection) | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 20387 | 5000 | 4329 | 671 |
| rna | trend_normalized | 36601 | 35329 | 20387 | 5000 | 4329 | 671 |
| atac | primary | 218882 | 212816 | 128916 | 5000 | 3360 | 1640 |
| atac | trend_normalized | 218882 | 212816 | 128916 | 5000 | 3360 | 1640 |

##### Normalization

After log2-CPM and before the selection step, each donor's log2-CPM over the abundance-passing features was corrected by limma's fast cyclic loess (`normalizeCyclicLoess(method = "fast")`): at each of 3 iterations the reference is the per-feature mean over donors (A), and each donor's deviation from it (M) is fitted by a lowess curve on A (3 robustifying iterations, delta = 0.01 of the range of A) and subtracted. The span is limma's adaptive default from the number of features (`LOESS_SPAN_RULE`: 0.3 + (1 − 0.3) · (50 / n)^(1/3)). Only the block's values are read (no label, covariate or technical variable). The tilt statistics are per donor, with A and M taken against the per-feature mean of the values before or after the normalization: the OLS slope of M on A and the median M of the lowest and highest of 20 equal-count bins of A; the quantiles are over donors and the per-donor values are in `normalization_trend.csv` beside the variant blocks (sha256 `f7192b22689e…`). Limitations: the normalization assumes that most features do not change at any abundance level, so a stage or sex effect shaped like a tilt is removed with the artefact; the curve is least determined at the rare end, where the selection draws most features; the lncRNA-versus-coding offset at matched abundance is not corrected; and no Type I check exists for the normalized blocks.

| omic | primary | variant | features | span | iterations | largest absolute correction |
|---|---|---|---:|---:|---:|---:|
| rna | none | cyclic_loess | 20387 | 0.394400 | 3 | 1.2155 |
| atac | none | cyclic_loess | 128916 | 0.351049 | 3 | 1.1286 |

| omic | statistic | values | q05 | q25 | q50 | q75 | q95 |
|---|---|---|---:|---:|---:|---:|---:|
| rna | slope of M on A | before | -0.0639 | -0.0342 | 0.0012 | 0.0294 | 0.0655 |
| rna | slope of M on A | after | -0.0140 | -0.0066 | -0.0026 | 0.0043 | 0.0188 |
| rna | median M, lowest bin | before | -0.3761 | -0.1549 | 0.0267 | 0.1926 | 0.4588 |
| rna | median M, lowest bin | after | -0.1090 | -0.0214 | 0.0102 | 0.0419 | 0.0651 |
| rna | median M, highest bin | before | -0.2846 | -0.1030 | 0.0048 | 0.1100 | 0.2457 |
| rna | median M, highest bin | after | -0.0145 | -0.0086 | -0.0031 | 0.0022 | 0.0164 |
| rna | largest absolute correction | per donor | 0.1311 | 0.2392 | 0.3396 | 0.5065 | 0.9425 |
| atac | slope of M on A | before | -0.1610 | -0.0483 | 0.0200 | 0.0550 | 0.0826 |
| atac | slope of M on A | after | -0.0247 | -0.0043 | 0.0035 | 0.0072 | 0.0154 |
| atac | median M, lowest bin | before | -0.9403 | -0.8906 | -0.0687 | 0.1940 | 0.6540 |
| atac | median M, lowest bin | after | -0.6540 | -0.4934 | -0.0842 | -0.0499 | 0.1628 |
| atac | median M, highest bin | before | -0.4164 | -0.1617 | 0.0252 | 0.1658 | 0.3549 |
| atac | median M, highest bin | after | -0.0118 | -0.0007 | 0.0085 | 0.0160 | 0.0297 |
| atac | largest absolute correction | per donor | 0.3428 | 0.4949 | 0.6157 | 0.7275 | 0.8875 |

##### Adoption check

`ADOPTION_RULE`: holds when no latent component of the candidate's CV-selected analysis reaches |Spearman| >= DIAGNOSTIC_CORRELATION_THRESHOLD against any variable of ADOPTION_CHECK_VARIABLES, in both cohorts; p-values, rejection statuses and CV AUROC are not read. If it holds, a later change may propose the candidate as a revised primary with the current primary kept as the reference; if not, the variants stay sensitivity analyses. Threshold 0.3 (`DIAGNOSTIC_CORRELATION_THRESHOLD`). In this cohort the rule does not hold: flagged pls_0 × `tmm_rna` (-0.67), pls_0 × `tmm_atac` (0.51), pls_0 × `fraction_mitochondrial_umis_mean` (0.73), pls_0 × `gex_reads_mapped_confidently_to_intronic_regions_mean` (-0.56), pls_0 × `gex_fraction_of_transcriptomic_reads_in_cells_mean` (-0.63), pls_0 × `rin` (-0.33). The rule needs both cohorts; whether the candidate becomes a revised primary is a later change with the current primary kept as the reference, and the candidate's statistics above are sensitivity, not the result. Between `trend_normalized` and `revised_candidate_trend`, `ADOPTION_PRECEDENCE` decides from both cohorts' `holds` values: if trend_normalized meets ADOPTION_RULE in both cohorts it is the candidate a later revised-primary change proposes, whatever revised_candidate_trend gives (no gene removed, no degrees of freedom spent, no technical covariate read); if only revised_candidate_trend meets it in both cohorts, that one is; if neither does, the variants stay sensitivity analyses. Only the holds values of the adoption tables are read; its outcome is recorded as `adoption_precedence` in the comparison manifests.

| component | variable | source | n | Spearman | flagged |
|---|---|---|---:|---:|:---:|
| pls_0 | `tmm_rna` | manifest.json (primary export, TMM factor) | 80 | -0.673 | yes |
| pls_0 | `tmm_atac` | manifest.json (primary export, TMM factor) | 80 | 0.505 | yes |
| pls_0 | `fraction_mitochondrial_umis_mean` | donor_technical.csv | 80 | 0.730 | yes |
| pls_0 | `doublet_score_mean` | donor_technical.csv | 80 | 0.018 | no |
| pls_0 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 80 | -0.565 | yes |
| pls_0 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 80 | -0.632 | yes |
| pls_0 | `rin` | covariates.csv (primary export) | 80 | -0.335 | yes |
| pls_0 | `pmi` | covariates.csv (primary export) | 80 | 0.032 | no |
| all | all | | 80 | | yes → holds = false |

##### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 1 | 0.6313 | 148.7 | `results/motco_analysis/cohort_30` |
| trend_normalized | cv | 1 | 0.6561 | 193.8 | `motco_analysis/cohort_30` |

##### Findings

- adoption check does not hold in cohort_30: pls_0 x tmm_rna (-0.67), pls_0 x tmm_atac (0.51), pls_0 x fraction_mitochondrial_umis_mean (0.73), pls_0 x gex_reads_mapped_confidently_to_intronic_regions_mean (-0.56), pls_0 x gex_fraction_of_transcriptomic_reads_in_cells_mean (-0.63), pls_0 x rin (-0.33); the variants stay sensitivity analyses
- trend_normalized export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

### Variant `revised_candidate_trend` (pasted from `results/motco_sensitivity/revised_candidate_trend/comparison/cohort_30/sensitivity.md`)

#### MOTCO preprocessing sensitivity: cohort_30, variant `revised_candidate_trend`

the revised candidate's four settings plus the cyclic-loess trend normalization of both omics, in the primary's order: gene sets removed before the abundance filter, cyclic loess after log2-CPM, residualization after the normalization, trend-ratio ranking for RNA and mean-accessibility ranking for ATAC at the selection step; the adoption rule ADOPTION_RULE is evaluated for this variant. Export setting(s) changed: `excluded_gene_sets` = ['ribosomal', 'mitochondrial', 'hemoglobin'], `normalization` = {'rna': 'cyclic_loess', 'atac': 'cyclic_loess'}, `residualize` = ['log10_lib_rna', 'log10_lib_atac', 'ribosomal_pseudobulk', 'mitochondrial_pseudobulk'], `selection_rule` = {'rna': 'trend_ratio', 'atac': 'mean'} (primary: `excluded_gene_sets` = [], `normalization` = {'rna': 'none', 'atac': 'none'}, `residualize` = None, `selection_rule` = {'rna': 'variance', 'atac': 'variance'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 1, which stays the primary result; the `revised_candidate_trend` row (CV-selected rank 3) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected, residualized nor trend-normalized), so the variant carries no Type I or power claim of its own.

##### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 14.941 | 0.3 | no | 11.2868 | 8.27982 | 26.5775 |
| revised_candidate_trend | cv | 3 | 14.3836 | 0.166 | no | 8.49356 | 6.54021 | 22.1517 |

##### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0 | 1 | no | 3.06306 | 23.2919 | 0 |
| revised_candidate_trend | cv | 3 | 16.7741 | 0.435 | no | 16.4397 | 9.28967 | 34.0151 |

##### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0.608832 | 0.333 | no | 0.489849 | 0.355498 | 1.18287 |
| revised_candidate_trend | cv | 3 | 0.0792388 | 0.852 | no | 0.188868 | 0.104185 | 0.3839 |

##### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after gene-set removal | after abundance | kept (after selection) | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 35329 | 20387 | 5000 | 2238 | 2762 |
| rna | revised_candidate_trend | 36601 | 35329 | 35249 | 20317 | 5000 | 2238 | 2762 |
| atac | primary | 218882 | 212816 | 212816 | 128916 | 5000 | 0 | 5000 |
| atac | revised_candidate_trend | 218882 | 212816 | 212816 | 128916 | 5000 | 0 | 5000 |

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
| rna | revised_candidate_trend | trend_ratio | -0.031 | 2.348 | 4.000 | 5.532 | 7.735 |
| atac | primary | variance | -0.321 | 0.046 | 0.371 | 0.778 | 1.531 |
| atac | revised_candidate_trend | mean | 5.433 | 5.614 | 5.881 | 6.285 | 7.179 |

##### Normalization

After log2-CPM and before the residualization and the selection step, each donor's log2-CPM over the abundance-passing features was corrected by limma's fast cyclic loess (`normalizeCyclicLoess(method = "fast")`): at each of 3 iterations the reference is the per-feature mean over donors (A), and each donor's deviation from it (M) is fitted by a lowess curve on A (3 robustifying iterations, delta = 0.01 of the range of A) and subtracted. The span is limma's adaptive default from the number of features (`LOESS_SPAN_RULE`: 0.3 + (1 − 0.3) · (50 / n)^(1/3)). Only the block's values are read (no label, covariate or technical variable). The tilt statistics are per donor, with A and M taken against the per-feature mean of the values before or after the normalization: the OLS slope of M on A and the median M of the lowest and highest of 20 equal-count bins of A; the quantiles are over donors and the per-donor values are in `normalization_trend.csv` beside the variant blocks (sha256 `e60d1a63f939…`). Limitations: the normalization assumes that most features do not change at any abundance level, so a stage or sex effect shaped like a tilt is removed with the artefact; the curve is least determined at the rare end, where the selection draws most features; the lncRNA-versus-coding offset at matched abundance is not corrected; and no Type I check exists for the normalized blocks.

| omic | primary | variant | features | span | iterations | largest absolute correction |
|---|---|---|---:|---:|---:|---:|
| rna | none | cyclic_loess | 20317 | 0.394508 | 3 | 1.2130 |
| atac | none | cyclic_loess | 128916 | 0.351049 | 3 | 1.1286 |

| omic | statistic | values | q05 | q25 | q50 | q75 | q95 |
|---|---|---|---:|---:|---:|---:|---:|
| rna | slope of M on A | before | -0.0639 | -0.0342 | 0.0013 | 0.0296 | 0.0660 |
| rna | slope of M on A | after | -0.0139 | -0.0065 | -0.0025 | 0.0043 | 0.0189 |
| rna | median M, lowest bin | before | -0.3774 | -0.1554 | 0.0277 | 0.1933 | 0.4569 |
| rna | median M, lowest bin | after | -0.0996 | -0.0237 | 0.0094 | 0.0385 | 0.0657 |
| rna | median M, highest bin | before | -0.2835 | -0.1029 | 0.0028 | 0.1104 | 0.2436 |
| rna | median M, highest bin | after | -0.0141 | -0.0087 | -0.0033 | 0.0022 | 0.0161 |
| rna | largest absolute correction | per donor | 0.1312 | 0.2384 | 0.3403 | 0.5049 | 0.9387 |
| atac | slope of M on A | before | -0.1610 | -0.0483 | 0.0200 | 0.0550 | 0.0826 |
| atac | slope of M on A | after | -0.0247 | -0.0043 | 0.0035 | 0.0072 | 0.0154 |
| atac | median M, lowest bin | before | -0.9403 | -0.8906 | -0.0687 | 0.1940 | 0.6540 |
| atac | median M, lowest bin | after | -0.6540 | -0.4934 | -0.0842 | -0.0499 | 0.1628 |
| atac | median M, highest bin | before | -0.4164 | -0.1617 | 0.0252 | 0.1658 | 0.3549 |
| atac | median M, highest bin | after | -0.0118 | -0.0007 | 0.0085 | 0.0160 | 0.0297 |
| atac | largest absolute correction | per donor | 0.3428 | 0.4949 | 0.6157 | 0.7275 | 0.8875 |

##### Residualization

After log2-CPM, the trend normalization and before the selection step, each abundance-passing feature was replaced by its mean plus the residual of its OLS fit on an intercept and 4 covariates: `log10_lib_rna`, `log10_lib_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` (`RESIDUALIZATION_COVARIATES`: the export's own log10 library sizes after step 2 and the technical table's gene-set pseudobulk fractions; excluded as cell state under the technical-axis report's rule: `hemoglobin_pseudobulk` (reports/motco-astrocyte-technical-axis-2026-10-19.md, section 2, 'Classification under the spec's rule': median within-donor supertype share 0.00394 above the between-donor share 0.00285 (both negligible; median pseudobulk fraction 6.7e-07 of UMIs); not a residualization candidate)). The design (`residualization.csv` beside the variant blocks, sha256 `e95605f9069f…`) holds no sex, stage or analysis covariate. Limitations: the fit spends 5 degrees of freedom per feature (intercept plus covariates) before the permutation null, which does not see them; it removes whatever part of the stage and sex signal the covariates carry; the analysis covariates were left unchanged, so the library-size proxies are adjusted twice (log10 library sizes here, log10 nuclei counts in the RRPP model); and no Type I check exists for the residualized blocks.

| omic | R² q05 | q25 | q50 | q75 | q95 |
|---|---:|---:|---:|---:|---:|
| rna | 0.037 | 0.096 | 0.169 | 0.267 | 0.478 |
| atac | 0.015 | 0.041 | 0.074 | 0.124 | 0.251 |

##### Adoption check

`ADOPTION_RULE`: holds when no latent component of the candidate's CV-selected analysis reaches |Spearman| >= DIAGNOSTIC_CORRELATION_THRESHOLD against any variable of ADOPTION_CHECK_VARIABLES, in both cohorts; p-values, rejection statuses and CV AUROC are not read. If it holds, a later change may propose the candidate as a revised primary with the current primary kept as the reference; if not, the variants stay sensitivity analyses. Threshold 0.3 (`DIAGNOSTIC_CORRELATION_THRESHOLD`). In this cohort the rule does not hold: flagged pls_1 × `tmm_rna` (-0.32). The rule needs both cohorts; whether the candidate becomes a revised primary is a later change with the current primary kept as the reference, and the candidate's statistics above are sensitivity, not the result. Between `trend_normalized` and `revised_candidate_trend`, `ADOPTION_PRECEDENCE` decides from both cohorts' `holds` values: if trend_normalized meets ADOPTION_RULE in both cohorts it is the candidate a later revised-primary change proposes, whatever revised_candidate_trend gives (no gene removed, no degrees of freedom spent, no technical covariate read); if only revised_candidate_trend meets it in both cohorts, that one is; if neither does, the variants stay sensitivity analyses. Only the holds values of the adoption tables are read; its outcome is recorded as `adoption_precedence` in the comparison manifests.

| component | variable | source | n | Spearman | flagged |
|---|---|---|---:|---:|:---:|
| pls_0 | `tmm_rna` | manifest.json (primary export, TMM factor) | 80 | -0.159 | no |
| pls_0 | `tmm_atac` | manifest.json (primary export, TMM factor) | 80 | 0.142 | no |
| pls_0 | `fraction_mitochondrial_umis_mean` | donor_technical.csv | 80 | 0.159 | no |
| pls_0 | `doublet_score_mean` | donor_technical.csv | 80 | 0.016 | no |
| pls_0 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 80 | -0.213 | no |
| pls_0 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 80 | -0.121 | no |
| pls_0 | `rin` | covariates.csv (primary export) | 80 | -0.020 | no |
| pls_0 | `pmi` | covariates.csv (primary export) | 80 | -0.002 | no |
| pls_1 | `tmm_rna` | manifest.json (primary export, TMM factor) | 80 | -0.323 | yes |
| pls_1 | `tmm_atac` | manifest.json (primary export, TMM factor) | 80 | 0.146 | no |
| pls_1 | `fraction_mitochondrial_umis_mean` | donor_technical.csv | 80 | 0.089 | no |
| pls_1 | `doublet_score_mean` | donor_technical.csv | 80 | 0.129 | no |
| pls_1 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 80 | -0.115 | no |
| pls_1 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 80 | -0.039 | no |
| pls_1 | `rin` | covariates.csv (primary export) | 80 | -0.048 | no |
| pls_1 | `pmi` | covariates.csv (primary export) | 80 | 0.084 | no |
| pls_2 | `tmm_rna` | manifest.json (primary export, TMM factor) | 80 | -0.234 | no |
| pls_2 | `tmm_atac` | manifest.json (primary export, TMM factor) | 80 | 0.172 | no |
| pls_2 | `fraction_mitochondrial_umis_mean` | donor_technical.csv | 80 | 0.167 | no |
| pls_2 | `doublet_score_mean` | donor_technical.csv | 80 | 0.017 | no |
| pls_2 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 80 | -0.063 | no |
| pls_2 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 80 | 0.046 | no |
| pls_2 | `rin` | covariates.csv (primary export) | 80 | 0.056 | no |
| pls_2 | `pmi` | covariates.csv (primary export) | 80 | 0.035 | no |
| all | all | | 80 | | yes → holds = false |

##### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 1 | 0.6313 | 148.7 | `results/motco_analysis/cohort_30` |
| revised_candidate_trend | cv | 3 | 0.4946 | 186.4 | `motco_analysis/cohort_30` |

##### Findings

- adoption check does not hold in cohort_30: pls_1 x tmm_rna (-0.32); the variants stay sensitivity analyses
- CV-selected rank differs: 1 in the primary analysis, 3 in the revised_candidate_trend analysis; reported beside angle, not read as better or worse
- revised_candidate_trend export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

## 3. Sensitivity cohort: cohort_50 (n = 74)

### Overview (generated from `sensitivity.csv` and `sensitivity_manifest.json`)

| analysis | rank | cv_mean_auroc | delta (p) | angle (p) | shape (p) | rejected | rna kept shared / only primary / only variant | atac kept shared / only primary / only variant | adoption check |
|---|---:|---:|---:|---:|---:|---|---|---|---|
| primary | 1 | 0.6616 | 18.7091 (0.213) | 0 (1) | 0.56471 (0.377) | none | — | — | — |
| trend_normalized | 1 | 0.6529 | 11.9963 (0.362) | 0 (1) | 0.504277 (0.363) | none | 4342 / 658 / 658 | 3257 / 1743 / 1743 | holds = false |
| revised_candidate_trend | 1 | 0.5461 | 5.95142 (0.495) | 0 (1) | 0.104564 (0.768) | none | 2272 / 2728 / 2728 | 0 / 5000 / 5000 | holds = true |

### Normalization: tilt before and after (generated from `sensitivity_manifest.json`, quantiles over donors of `normalization_trend.csv`)

Per donor, A is the per-feature mean over donors and M the donor's value minus A, of the log2-CPM before or after the cyclic loess; slope = OLS slope of M on A; the bins are 20 equal-count bins of A. The same rows are in each variant's pasted `sensitivity.md` below.

| variant | omic | features | span | statistic | before q05 | q25 | q50 | q75 | q95 | after q05 | q25 | q50 | q75 | q95 |
|---|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| trend_normalized | rna | 20238 | 0.394631 | slope of M on A | -0.0626 | -0.0337 | 0.0001 | 0.0299 | 0.0673 | -0.0138 | -0.0061 | -0.0029 | 0.0049 | 0.0176 |
| trend_normalized | rna | 20238 | 0.394631 | median M, lowest bin | -0.3766 | -0.1769 | 0.0506 | 0.2020 | 0.4238 | -0.1472 | -0.0196 | 0.0179 | 0.0407 | 0.0767 |
| trend_normalized | rna | 20238 | 0.394631 | median M, highest bin | -0.2576 | -0.1153 | 0.0109 | 0.1169 | 0.2688 | -0.0146 | -0.0080 | -0.0035 | 0.0025 | 0.0143 |
| trend_normalized | rna | 20238 | 0.394631 | largest absolute correction (per donor) | 0.1478 | 0.2433 | 0.3612 | 0.5583 | 0.9439 | | | | | |
| trend_normalized | atac | 124818 | 0.351601 | slope of M on A | -0.1629 | -0.0510 | 0.0186 | 0.0605 | 0.0913 | -0.0235 | -0.0040 | 0.0019 | 0.0069 | 0.0169 |
| trend_normalized | atac | 124818 | 0.351601 | median M, lowest bin | -1.0531 | -0.8624 | -0.0226 | 0.2940 | 0.6626 | -0.6301 | -0.3593 | -0.0674 | -0.0205 | 0.2357 |
| trend_normalized | atac | 124818 | 0.351601 | median M, highest bin | -0.4462 | -0.1349 | 0.0489 | 0.1695 | 0.3353 | -0.0134 | -0.0016 | 0.0081 | 0.0161 | 0.0298 |
| trend_normalized | atac | 124818 | 0.351601 | largest absolute correction (per donor) | 0.3368 | 0.4902 | 0.6114 | 0.8138 | 0.9805 | | | | | |
| revised_candidate_trend | rna | 20166 | 0.394743 | slope of M on A | -0.0627 | -0.0337 | 0.0009 | 0.0301 | 0.0674 | -0.0135 | -0.0061 | -0.0028 | 0.0049 | 0.0175 |
| revised_candidate_trend | rna | 20166 | 0.394743 | median M, lowest bin | -0.3886 | -0.1829 | 0.0441 | 0.2033 | 0.4315 | -0.1469 | -0.0235 | 0.0148 | 0.0415 | 0.0769 |
| revised_candidate_trend | rna | 20166 | 0.394743 | median M, highest bin | -0.2553 | -0.1044 | 0.0040 | 0.1160 | 0.2627 | -0.0149 | -0.0085 | -0.0033 | 0.0027 | 0.0142 |
| revised_candidate_trend | rna | 20166 | 0.394743 | largest absolute correction (per donor) | 0.1561 | 0.2340 | 0.3503 | 0.5367 | 0.9357 | | | | | |
| revised_candidate_trend | atac | 124818 | 0.351601 | slope of M on A | -0.1629 | -0.0510 | 0.0186 | 0.0605 | 0.0913 | -0.0235 | -0.0040 | 0.0019 | 0.0069 | 0.0169 |
| revised_candidate_trend | atac | 124818 | 0.351601 | median M, lowest bin | -1.0531 | -0.8624 | -0.0226 | 0.2940 | 0.6626 | -0.6301 | -0.3593 | -0.0674 | -0.0205 | 0.2357 |
| revised_candidate_trend | atac | 124818 | 0.351601 | median M, highest bin | -0.4462 | -0.1349 | 0.0489 | 0.1695 | 0.3353 | -0.0134 | -0.0016 | 0.0081 | 0.0161 | 0.0298 |
| revised_candidate_trend | atac | 124818 | 0.351601 | largest absolute correction (per donor) | 0.3368 | 0.4902 | 0.6114 | 0.8138 | 0.9805 | | | | | |

### Tilt and the first latent axis (generated from `normalization_trend.csv` and `latent_scores.csv`)

Spearman correlation over donors of each tilt statistic with `pls_0` of the primary analysis and of the variant's analysis. For `trend_normalized` the `before` values are computed on the primary pipeline's log2-CPM, so its primary-`pls_0` column of the `before` rows reproduces the label-free check of the proposal from a step output; for `revised_candidate_trend` they are computed after its gene-set removal. Descriptive: no p-value is attached.

| variant | omic | statistic | before × primary pls_0 | before × variant pls_0 | after × primary pls_0 | after × variant pls_0 |
|---|---|---|---:|---:|---:|---:|
| trend_normalized | rna | slope of M on A | 0.89 | 0.82 | -0.29 | -0.25 |
| trend_normalized | rna | median M, lowest bin | -0.91 | -0.84 | -0.35 | -0.37 |
| trend_normalized | rna | median M, highest bin | 0.74 | 0.72 | -0.09 | -0.10 |
| trend_normalized | rna | largest absolute correction | -0.17 | -0.18 | | |
| trend_normalized | atac | slope of M on A | -0.81 | -0.73 | 0.52 | 0.54 |
| trend_normalized | atac | median M, lowest bin | 0.72 | 0.62 | 0.25 | 0.18 |
| trend_normalized | atac | median M, highest bin | -0.64 | -0.64 | -0.61 | -0.62 |
| trend_normalized | atac | largest absolute correction | -0.30 | -0.20 | | |
| revised_candidate_trend | rna | slope of M on A | 0.89 | 0.18 | -0.30 | -0.26 |
| revised_candidate_trend | rna | median M, lowest bin | -0.91 | -0.23 | -0.35 | 0.02 |
| revised_candidate_trend | rna | median M, highest bin | 0.74 | -0.01 | -0.08 | -0.18 |
| revised_candidate_trend | rna | largest absolute correction | -0.17 | 0.03 | | |
| revised_candidate_trend | atac | slope of M on A | -0.81 | -0.29 | 0.52 | 0.07 |
| revised_candidate_trend | atac | median M, lowest bin | 0.72 | 0.28 | 0.25 | 0.16 |
| revised_candidate_trend | atac | median M, highest bin | -0.64 | -0.21 | -0.61 | -0.12 |
| revised_candidate_trend | atac | largest absolute correction | -0.30 | 0.05 | | |

Per-donor mean z-score over the selected block (each feature z-scored over donors, population SD) against `pls_0`:

| export | omic | × primary pls_0 | × pls_0 of its own analysis |
|---|---|---:|---:|
| primary | rna | -0.93 | -0.93 |
| primary | atac | 0.59 | 0.59 |
| trend_normalized | rna | -0.71 | -0.70 |
| trend_normalized | atac | -0.36 | -0.49 |
| revised_candidate_trend | rna | -0.31 | -0.36 |
| revised_candidate_trend | atac | -0.04 | -0.21 |

### Attribution: the first latent axis per analysis (generated from `axis_correlations.csv` and `geometry_shares.csv`; archived rows from the feature-selection report)

Spearman correlation of each latent component with stage, the astrocyte ribosomal-protein UMI fraction `ribosomal_pseudobulk`, the RNA TMM factor `tmm_rna`, the ATAC TSS enrichment `atac_tss_enrichment_score_mean`, the ATAC log10 library size `log10_lib_atac`, the mitochondrial fraction `mitochondrial_pseudobulk` and sex, per analysis; the component's score SD; and the share of the sex-difference geometry (Σ‖sex difference by stage‖²) along it. The archived rows (`pls_0` only) are pasted from `reports/motco-astrocyte-feature-selection-2026-10-19.md`, which did not tabulate the TSS enrichment (`n/r`); the primary row computed here equals its archived row. Descriptive: no p-value is attached.

| analysis | rank | component | stage | ribosomal fraction | RNA TMM factor | ATAC TSS enrichment | ATAC log10 library size | mitochondrial fraction | sex (male) | component SD | sex-difference share |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| primary | 1 | pls_0 | 0.50 | -0.82 | -0.72 | -0.77 | 0.23 | 0.77 | -0.19 | 35.8 | 1.000 |
| exclude_families (archived) | 1 | pls_0 | 0.50 | -0.81 | -0.70 | n/r | 0.23 | 0.77 | -0.19 | pls_0 35.9 | 1.000 |
| rna_trend_ratio (archived) | 1 | pls_0 | 0.51 | -0.82 | -0.72 | n/r | 0.22 | 0.74 | -0.18 | pls_0 35.3 | 1.000 |
| atac_mean_accessibility (archived) | 1 | pls_0 | 0.36 | -0.73 | -0.62 | n/r | -0.11 | 0.73 | -0.28 | pls_0 56.0 | 1.000 |
| residualized (archived) | 2 | pls_0 | 0.68 | -0.03 | -0.16 | n/r | -0.05 | 0.16 | 0.01 | pls_0 16.0, pls_1 30.0 | 0.128 |
| revised_candidate (archived) | 2 | pls_0 | 0.52 | -0.08 | -0.29 | n/r | 0.02 | 0.18 | 0.01 | pls_0 21.8, pls_1 30.5 | 0.230 |
| trend_normalized | 1 | pls_0 | 0.55 | -0.81 | -0.66 | -0.78 | 0.14 | 0.76 | -0.23 | 33.0 | 1.000 |
| revised_candidate_trend | 1 | pls_0 | 0.63 | -0.07 | -0.21 | -0.19 | 0.00 | 0.19 | -0.04 | 19.8 | 1.000 |

### Adoption check of `trend_normalized` (pasted from `adoption_check.csv`)

Rule: holds when no latent component of the candidate's CV-selected analysis reaches |Spearman| >= DIAGNOSTIC_CORRELATION_THRESHOLD against any variable of ADOPTION_CHECK_VARIABLES, in both cohorts; p-values, rejection statuses and CV AUROC are not read. If it holds, a later change may propose the candidate as a revised primary with the current primary kept as the reference; if not, the variants stay sensitivity analyses. Threshold 0.3. Outcome in cohort_50: **holds = false** (6 of 8 cells flagged).

| component | variable | source | n | spearman | flag | holds |
|---|---|---|---:|---:|:---:|:---:|
| pls_0 | tmm_rna | manifest.json (primary export, TMM factor) | 74 | -0.6610144390966308 | True |  |
| pls_0 | tmm_atac | manifest.json (primary export, TMM factor) | 74 | 0.5777267678637541 | True |  |
| pls_0 | fraction_mitochondrial_umis_mean | donor_technical.csv | 74 | 0.7825101814142909 | True |  |
| pls_0 | doublet_score_mean | donor_technical.csv | 74 | 0.010914476119955569 | False |  |
| pls_0 | gex_reads_mapped_confidently_to_intronic_regions_mean | donor_technical.csv | 74 | -0.5957837525802175 | True |  |
| pls_0 | gex_fraction_of_transcriptomic_reads_in_cells_mean | donor_technical.csv | 74 | -0.6339726027397259 | True |  |
| pls_0 | rin | covariates.csv (primary export) | 74 | -0.3407567667722281 | True |  |
| pls_0 | pmi | covariates.csv (primary export) | 74 | 0.07235127157403248 | False |  |
| all | all |  | 74 |  | True | False |

### Adoption check of `revised_candidate_trend` (pasted from `adoption_check.csv`)

Rule: holds when no latent component of the candidate's CV-selected analysis reaches |Spearman| >= DIAGNOSTIC_CORRELATION_THRESHOLD against any variable of ADOPTION_CHECK_VARIABLES, in both cohorts; p-values, rejection statuses and CV AUROC are not read. If it holds, a later change may propose the candidate as a revised primary with the current primary kept as the reference; if not, the variants stay sensitivity analyses. Threshold 0.3. Outcome in cohort_50: **holds = true** (0 of 8 cells flagged).

| component | variable | source | n | spearman | flag | holds |
|---|---|---|---:|---:|:---:|:---:|
| pls_0 | tmm_rna | manifest.json (primary export, TMM factor) | 74 | -0.2094187338022954 | False |  |
| pls_0 | tmm_atac | manifest.json (primary export, TMM factor) | 74 | 0.2332913735653461 | False |  |
| pls_0 | fraction_mitochondrial_umis_mean | donor_technical.csv | 74 | 0.20592373195112915 | False |  |
| pls_0 | doublet_score_mean | donor_technical.csv | 74 | 0.026019992595335054 | False |  |
| pls_0 | gex_reads_mapped_confidently_to_intronic_regions_mean | donor_technical.csv | 74 | -0.2753667187540782 | False |  |
| pls_0 | gex_fraction_of_transcriptomic_reads_in_cells_mean | donor_technical.csv | 74 | -0.1572602739726027 | False |  |
| pls_0 | rin | covariates.csv (primary export) | 74 | -0.03738622383538812 | False |  |
| pls_0 | pmi | covariates.csv (primary export) | 74 | 0.035153474087874496 | False |  |
| all | all |  | 74 |  | False | True |

### Variant `trend_normalized` (pasted from `results/motco_sensitivity/trend_normalized/comparison/cohort_50/sensitivity.md`)

#### MOTCO preprocessing sensitivity: cohort_50, variant `trend_normalized`

after log2-CPM and before the selection step, each donor's log2-CPM over the abundance-passing features corrected by limma's fast cyclic loess (a lowess curve of its deviation from the per-feature mean over donors against that mean, three iterations, adaptive span) in both omics; every other export setting equal to the primary; the adoption rule ADOPTION_RULE is evaluated for this variant. Export setting(s) changed: `normalization` = {'rna': 'cyclic_loess', 'atac': 'cyclic_loess'} (primary: `normalization` = {'rna': 'none', 'atac': 'none'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 1, which stays the primary result; the `trend_normalized` row (CV-selected rank 1) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected, residualized nor trend-normalized), so the variant carries no Type I or power claim of its own.

##### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 18.7091 | 0.213 | no | 11.6918 | 9.218 | 30.0679 |
| trend_normalized | cv | 1 | 11.9963 | 0.362 | no | 10.5288 | 8.17474 | 26.3734 |

##### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0 | 1 | no | 1.98198 | 18.7931 | 0 |
| trend_normalized | cv | 1 | 0 | 1 | no | 0.18018 | 5.69495 | 0 |

##### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0.56471 | 0.377 | no | 0.495877 | 0.343673 | 1.159 |
| trend_normalized | cv | 1 | 0.504277 | 0.363 | no | 0.432218 | 0.30948 | 1.01479 |

##### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after abundance | kept (after selection) | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 20238 | 5000 | 4342 | 658 |
| rna | trend_normalized | 36601 | 35329 | 20238 | 5000 | 4342 | 658 |
| atac | primary | 218882 | 212816 | 124818 | 5000 | 3257 | 1743 |
| atac | trend_normalized | 218882 | 212816 | 124818 | 5000 | 3257 | 1743 |

##### Normalization

After log2-CPM and before the selection step, each donor's log2-CPM over the abundance-passing features was corrected by limma's fast cyclic loess (`normalizeCyclicLoess(method = "fast")`): at each of 3 iterations the reference is the per-feature mean over donors (A), and each donor's deviation from it (M) is fitted by a lowess curve on A (3 robustifying iterations, delta = 0.01 of the range of A) and subtracted. The span is limma's adaptive default from the number of features (`LOESS_SPAN_RULE`: 0.3 + (1 − 0.3) · (50 / n)^(1/3)). Only the block's values are read (no label, covariate or technical variable). The tilt statistics are per donor, with A and M taken against the per-feature mean of the values before or after the normalization: the OLS slope of M on A and the median M of the lowest and highest of 20 equal-count bins of A; the quantiles are over donors and the per-donor values are in `normalization_trend.csv` beside the variant blocks (sha256 `8a60e38aab48…`). Limitations: the normalization assumes that most features do not change at any abundance level, so a stage or sex effect shaped like a tilt is removed with the artefact; the curve is least determined at the rare end, where the selection draws most features; the lncRNA-versus-coding offset at matched abundance is not corrected; and no Type I check exists for the normalized blocks.

| omic | primary | variant | features | span | iterations | largest absolute correction |
|---|---|---|---:|---:|---:|---:|
| rna | none | cyclic_loess | 20238 | 0.394631 | 3 | 1.2402 |
| atac | none | cyclic_loess | 124818 | 0.351601 | 3 | 1.1331 |

| omic | statistic | values | q05 | q25 | q50 | q75 | q95 |
|---|---|---|---:|---:|---:|---:|---:|
| rna | slope of M on A | before | -0.0626 | -0.0337 | 0.0001 | 0.0299 | 0.0673 |
| rna | slope of M on A | after | -0.0138 | -0.0061 | -0.0029 | 0.0049 | 0.0176 |
| rna | median M, lowest bin | before | -0.3766 | -0.1769 | 0.0506 | 0.2020 | 0.4238 |
| rna | median M, lowest bin | after | -0.1472 | -0.0196 | 0.0179 | 0.0407 | 0.0767 |
| rna | median M, highest bin | before | -0.2576 | -0.1153 | 0.0109 | 0.1169 | 0.2688 |
| rna | median M, highest bin | after | -0.0146 | -0.0080 | -0.0035 | 0.0025 | 0.0143 |
| rna | largest absolute correction | per donor | 0.1478 | 0.2433 | 0.3612 | 0.5583 | 0.9439 |
| atac | slope of M on A | before | -0.1629 | -0.0510 | 0.0186 | 0.0605 | 0.0913 |
| atac | slope of M on A | after | -0.0235 | -0.0040 | 0.0019 | 0.0069 | 0.0169 |
| atac | median M, lowest bin | before | -1.0531 | -0.8624 | -0.0226 | 0.2940 | 0.6626 |
| atac | median M, lowest bin | after | -0.6301 | -0.3593 | -0.0674 | -0.0205 | 0.2357 |
| atac | median M, highest bin | before | -0.4462 | -0.1349 | 0.0489 | 0.1695 | 0.3353 |
| atac | median M, highest bin | after | -0.0134 | -0.0016 | 0.0081 | 0.0161 | 0.0298 |
| atac | largest absolute correction | per donor | 0.3368 | 0.4902 | 0.6114 | 0.8138 | 0.9805 |

##### Adoption check

`ADOPTION_RULE`: holds when no latent component of the candidate's CV-selected analysis reaches |Spearman| >= DIAGNOSTIC_CORRELATION_THRESHOLD against any variable of ADOPTION_CHECK_VARIABLES, in both cohorts; p-values, rejection statuses and CV AUROC are not read. If it holds, a later change may propose the candidate as a revised primary with the current primary kept as the reference; if not, the variants stay sensitivity analyses. Threshold 0.3 (`DIAGNOSTIC_CORRELATION_THRESHOLD`). In this cohort the rule does not hold: flagged pls_0 × `tmm_rna` (-0.66), pls_0 × `tmm_atac` (0.58), pls_0 × `fraction_mitochondrial_umis_mean` (0.78), pls_0 × `gex_reads_mapped_confidently_to_intronic_regions_mean` (-0.60), pls_0 × `gex_fraction_of_transcriptomic_reads_in_cells_mean` (-0.63), pls_0 × `rin` (-0.34). The rule needs both cohorts; whether the candidate becomes a revised primary is a later change with the current primary kept as the reference, and the candidate's statistics above are sensitivity, not the result. Between `trend_normalized` and `revised_candidate_trend`, `ADOPTION_PRECEDENCE` decides from both cohorts' `holds` values: if trend_normalized meets ADOPTION_RULE in both cohorts it is the candidate a later revised-primary change proposes, whatever revised_candidate_trend gives (no gene removed, no degrees of freedom spent, no technical covariate read); if only revised_candidate_trend meets it in both cohorts, that one is; if neither does, the variants stay sensitivity analyses. Only the holds values of the adoption tables are read; its outcome is recorded as `adoption_precedence` in the comparison manifests.

| component | variable | source | n | Spearman | flagged |
|---|---|---|---:|---:|:---:|
| pls_0 | `tmm_rna` | manifest.json (primary export, TMM factor) | 74 | -0.661 | yes |
| pls_0 | `tmm_atac` | manifest.json (primary export, TMM factor) | 74 | 0.578 | yes |
| pls_0 | `fraction_mitochondrial_umis_mean` | donor_technical.csv | 74 | 0.783 | yes |
| pls_0 | `doublet_score_mean` | donor_technical.csv | 74 | 0.011 | no |
| pls_0 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 74 | -0.596 | yes |
| pls_0 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 74 | -0.634 | yes |
| pls_0 | `rin` | covariates.csv (primary export) | 74 | -0.341 | yes |
| pls_0 | `pmi` | covariates.csv (primary export) | 74 | 0.072 | no |
| all | all | | 74 | | yes → holds = false |

##### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 1 | 0.6616 | 151.6 | `results/motco_analysis/cohort_50` |
| trend_normalized | cv | 1 | 0.6529 | 201.9 | `motco_analysis/cohort_50` |

##### Findings

- adoption check does not hold in cohort_50: pls_0 x tmm_rna (-0.66), pls_0 x tmm_atac (0.58), pls_0 x fraction_mitochondrial_umis_mean (0.78), pls_0 x gex_reads_mapped_confidently_to_intronic_regions_mean (-0.60), pls_0 x gex_fraction_of_transcriptomic_reads_in_cells_mean (-0.63), pls_0 x rin (-0.34); the variants stay sensitivity analyses
- trend_normalized export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

### Variant `revised_candidate_trend` (pasted from `results/motco_sensitivity/revised_candidate_trend/comparison/cohort_50/sensitivity.md`)

#### MOTCO preprocessing sensitivity: cohort_50, variant `revised_candidate_trend`

the revised candidate's four settings plus the cyclic-loess trend normalization of both omics, in the primary's order: gene sets removed before the abundance filter, cyclic loess after log2-CPM, residualization after the normalization, trend-ratio ranking for RNA and mean-accessibility ranking for ATAC at the selection step; the adoption rule ADOPTION_RULE is evaluated for this variant. Export setting(s) changed: `excluded_gene_sets` = ['ribosomal', 'mitochondrial', 'hemoglobin'], `normalization` = {'rna': 'cyclic_loess', 'atac': 'cyclic_loess'}, `residualize` = ['log10_lib_rna', 'log10_lib_atac', 'ribosomal_pseudobulk', 'mitochondrial_pseudobulk'], `selection_rule` = {'rna': 'trend_ratio', 'atac': 'mean'} (primary: `excluded_gene_sets` = [], `normalization` = {'rna': 'none', 'atac': 'none'}, `residualize` = None, `selection_rule` = {'rna': 'variance', 'atac': 'variance'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 1, which stays the primary result; the `revised_candidate_trend` row (CV-selected rank 1) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected, residualized nor trend-normalized), so the variant carries no Type I or power claim of its own.

##### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 18.7091 | 0.213 | no | 11.6918 | 9.218 | 30.0679 |
| revised_candidate_trend | cv | 1 | 5.95142 | 0.495 | no | 6.89976 | 5.24076 | 17.0778 |

##### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0 | 1 | no | 1.98198 | 18.7931 | 0 |
| revised_candidate_trend | cv | 1 | 0 | 1 | no | 0 | 0 | 0 |

##### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0.56471 | 0.377 | no | 0.495877 | 0.343673 | 1.159 |
| revised_candidate_trend | cv | 1 | 0.104564 | 0.768 | no | 0.278634 | 0.199673 | 0.646356 |

##### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after gene-set removal | after abundance | kept (after selection) | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 35329 | 20238 | 5000 | 2272 | 2728 |
| rna | revised_candidate_trend | 36601 | 35329 | 35249 | 20166 | 5000 | 2272 | 2728 |
| atac | primary | 218882 | 212816 | 212816 | 124818 | 5000 | 0 | 5000 |
| atac | revised_candidate_trend | 218882 | 212816 | 212816 | 124818 | 5000 | 0 | 5000 |

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
| rna | revised_candidate_trend | trend_ratio | -0.027 | 2.354 | 4.000 | 5.524 | 7.718 |
| atac | primary | variance | -0.432 | -0.057 | 0.271 | 0.698 | 1.479 |
| atac | revised_candidate_trend | mean | 5.440 | 5.620 | 5.886 | 6.290 | 7.183 |

##### Normalization

After log2-CPM and before the residualization and the selection step, each donor's log2-CPM over the abundance-passing features was corrected by limma's fast cyclic loess (`normalizeCyclicLoess(method = "fast")`): at each of 3 iterations the reference is the per-feature mean over donors (A), and each donor's deviation from it (M) is fitted by a lowess curve on A (3 robustifying iterations, delta = 0.01 of the range of A) and subtracted. The span is limma's adaptive default from the number of features (`LOESS_SPAN_RULE`: 0.3 + (1 − 0.3) · (50 / n)^(1/3)). Only the block's values are read (no label, covariate or technical variable). The tilt statistics are per donor, with A and M taken against the per-feature mean of the values before or after the normalization: the OLS slope of M on A and the median M of the lowest and highest of 20 equal-count bins of A; the quantiles are over donors and the per-donor values are in `normalization_trend.csv` beside the variant blocks (sha256 `a88681aba998…`). Limitations: the normalization assumes that most features do not change at any abundance level, so a stage or sex effect shaped like a tilt is removed with the artefact; the curve is least determined at the rare end, where the selection draws most features; the lncRNA-versus-coding offset at matched abundance is not corrected; and no Type I check exists for the normalized blocks.

| omic | primary | variant | features | span | iterations | largest absolute correction |
|---|---|---|---:|---:|---:|---:|
| rna | none | cyclic_loess | 20166 | 0.394743 | 3 | 1.2390 |
| atac | none | cyclic_loess | 124818 | 0.351601 | 3 | 1.1331 |

| omic | statistic | values | q05 | q25 | q50 | q75 | q95 |
|---|---|---|---:|---:|---:|---:|---:|
| rna | slope of M on A | before | -0.0627 | -0.0337 | 0.0009 | 0.0301 | 0.0674 |
| rna | slope of M on A | after | -0.0135 | -0.0061 | -0.0028 | 0.0049 | 0.0175 |
| rna | median M, lowest bin | before | -0.3886 | -0.1829 | 0.0441 | 0.2033 | 0.4315 |
| rna | median M, lowest bin | after | -0.1469 | -0.0235 | 0.0148 | 0.0415 | 0.0769 |
| rna | median M, highest bin | before | -0.2553 | -0.1044 | 0.0040 | 0.1160 | 0.2627 |
| rna | median M, highest bin | after | -0.0149 | -0.0085 | -0.0033 | 0.0027 | 0.0142 |
| rna | largest absolute correction | per donor | 0.1561 | 0.2340 | 0.3503 | 0.5367 | 0.9357 |
| atac | slope of M on A | before | -0.1629 | -0.0510 | 0.0186 | 0.0605 | 0.0913 |
| atac | slope of M on A | after | -0.0235 | -0.0040 | 0.0019 | 0.0069 | 0.0169 |
| atac | median M, lowest bin | before | -1.0531 | -0.8624 | -0.0226 | 0.2940 | 0.6626 |
| atac | median M, lowest bin | after | -0.6301 | -0.3593 | -0.0674 | -0.0205 | 0.2357 |
| atac | median M, highest bin | before | -0.4462 | -0.1349 | 0.0489 | 0.1695 | 0.3353 |
| atac | median M, highest bin | after | -0.0134 | -0.0016 | 0.0081 | 0.0161 | 0.0298 |
| atac | largest absolute correction | per donor | 0.3368 | 0.4902 | 0.6114 | 0.8138 | 0.9805 |

##### Residualization

After log2-CPM, the trend normalization and before the selection step, each abundance-passing feature was replaced by its mean plus the residual of its OLS fit on an intercept and 4 covariates: `log10_lib_rna`, `log10_lib_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` (`RESIDUALIZATION_COVARIATES`: the export's own log10 library sizes after step 2 and the technical table's gene-set pseudobulk fractions; excluded as cell state under the technical-axis report's rule: `hemoglobin_pseudobulk` (reports/motco-astrocyte-technical-axis-2026-10-19.md, section 2, 'Classification under the spec's rule': median within-donor supertype share 0.00394 above the between-donor share 0.00285 (both negligible; median pseudobulk fraction 6.7e-07 of UMIs); not a residualization candidate)). The design (`residualization.csv` beside the variant blocks, sha256 `443e90258051…`) holds no sex, stage or analysis covariate. Limitations: the fit spends 5 degrees of freedom per feature (intercept plus covariates) before the permutation null, which does not see them; it removes whatever part of the stage and sex signal the covariates carry; the analysis covariates were left unchanged, so the library-size proxies are adjusted twice (log10 library sizes here, log10 nuclei counts in the RRPP model); and no Type I check exists for the residualized blocks.

| omic | R² q05 | q25 | q50 | q75 | q95 |
|---|---:|---:|---:|---:|---:|
| rna | 0.039 | 0.104 | 0.182 | 0.286 | 0.501 |
| atac | 0.016 | 0.044 | 0.079 | 0.135 | 0.277 |

##### Adoption check

`ADOPTION_RULE`: holds when no latent component of the candidate's CV-selected analysis reaches |Spearman| >= DIAGNOSTIC_CORRELATION_THRESHOLD against any variable of ADOPTION_CHECK_VARIABLES, in both cohorts; p-values, rejection statuses and CV AUROC are not read. If it holds, a later change may propose the candidate as a revised primary with the current primary kept as the reference; if not, the variants stay sensitivity analyses. Threshold 0.3 (`DIAGNOSTIC_CORRELATION_THRESHOLD`). In this cohort the rule holds. The rule needs both cohorts; whether the candidate becomes a revised primary is a later change with the current primary kept as the reference, and the candidate's statistics above are sensitivity, not the result. Between `trend_normalized` and `revised_candidate_trend`, `ADOPTION_PRECEDENCE` decides from both cohorts' `holds` values: if trend_normalized meets ADOPTION_RULE in both cohorts it is the candidate a later revised-primary change proposes, whatever revised_candidate_trend gives (no gene removed, no degrees of freedom spent, no technical covariate read); if only revised_candidate_trend meets it in both cohorts, that one is; if neither does, the variants stay sensitivity analyses. Only the holds values of the adoption tables are read; its outcome is recorded as `adoption_precedence` in the comparison manifests.

| component | variable | source | n | Spearman | flagged |
|---|---|---|---:|---:|:---:|
| pls_0 | `tmm_rna` | manifest.json (primary export, TMM factor) | 74 | -0.209 | no |
| pls_0 | `tmm_atac` | manifest.json (primary export, TMM factor) | 74 | 0.233 | no |
| pls_0 | `fraction_mitochondrial_umis_mean` | donor_technical.csv | 74 | 0.206 | no |
| pls_0 | `doublet_score_mean` | donor_technical.csv | 74 | 0.026 | no |
| pls_0 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 74 | -0.275 | no |
| pls_0 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 74 | -0.157 | no |
| pls_0 | `rin` | covariates.csv (primary export) | 74 | -0.037 | no |
| pls_0 | `pmi` | covariates.csv (primary export) | 74 | 0.035 | no |
| all | all | | 74 | | no → holds = true |

##### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 1 | 0.6616 | 151.6 | `results/motco_analysis/cohort_50` |
| revised_candidate_trend | cv | 1 | 0.5461 | 183.4 | `motco_analysis/cohort_50` |

##### Findings

- adoption check holds in cohort_50: no candidate component reaches |Spearman| >= 0.3 against any of 8 check variables (the rule needs both cohorts; a revised primary is a later change)
- revised_candidate_trend export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

## 4. Adoption and precedence

| variant | cohort_30 holds | cohort_30 flagged cells | cohort_50 holds | cohort_50 flagged cells | meets ADOPTION_RULE (both cohorts) |
|---|:---:|---|:---:|---|:---:|
| trend_normalized | false | `pls_0` × `tmm_rna` (-0.67), `pls_0` × `tmm_atac` (0.51), `pls_0` × `fraction_mitochondrial_umis_mean` (0.73), `pls_0` × `gex_reads_mapped_confidently_to_intronic_regions_mean` (-0.56), `pls_0` × `gex_fraction_of_transcriptomic_reads_in_cells_mean` (-0.63), `pls_0` × `rin` (-0.33) | false | `pls_0` × `tmm_rna` (-0.66), `pls_0` × `tmm_atac` (0.58), `pls_0` × `fraction_mitochondrial_umis_mean` (0.78), `pls_0` × `gex_reads_mapped_confidently_to_intronic_regions_mean` (-0.60), `pls_0` × `gex_fraction_of_transcriptomic_reads_in_cells_mean` (-0.63), `pls_0` × `rin` (-0.34) | no |
| revised_candidate_trend | false | `pls_1` × `tmm_rna` (-0.32) | true | none | no |

`adoption_precedence` (identical in the four comparison manifests): candidate **null**; reason: "neither trend_normalized nor revised_candidate_trend meets ADOPTION_RULE in both cohorts; no revised-primary change follows and the variants stay sensitivity analyses". The rule read only the `holds` values; the flagged cells are listed beside it.

## 5. Reading

- **Statistics.** No statistic changed rejection status in either variant or cohort: nothing is rejected in the primary analysis and nothing is rejected in either variant analysis, so the non-rejection rules apply to every row. In particular the `trend_normalized` rows do not reject `angle` (0° at rank 1 in both cohorts, by construction) or `shape`, which is uninformative about orientation and shape, not evidence that the directions or shapes agree; its `delta` (8.67659, p 0.504; cohort_50 11.9963, p 0.362) is no evidence about a size difference.
- **Ranks.** The CV-selected rank differs in cohort_30 `revised_candidate_trend` (1 in the primary analysis, 3 in the variant analysis; angle 16.7741° at p 0.435, a non-rejection and therefore uninformative), reported beside `angle` and not read as better or worse; every other row selects rank 1.
- **The tilt is removed.** In `trend_normalized` the median over donors of the per-donor slope of M on A stays near zero while its spread shrinks: rna q05–q95 -0.0639 to 0.0655 before, -0.0140 to 0.0188 after (cohort_30), atac -0.1610 to 0.0826 before, -0.0247 to 0.0154 after; the highest-bin medians shrink alike. The slope's correlation with the primary `pls_0` falls from 0.90 to -0.29 for RNA (cohort_50 0.89 to -0.29). The lowest-abundance bin of ATAC is the exception: its median M after normalization keeps a wide spread (q05 -0.6540, q95 0.1628 in cohort_30), where the curve is least determined. Over twenty equal-count bins, every normalized bin median is closer to zero than the largest un-normalized one in all four blocks (largest |median| after / before: rna 1.14 / 2.48, atac 1.05 / 1.98 in cohort_30), and the donor with the largest slope loses most of it (task 4.5).
- **Attribution of the first axis: the tilt does not carry it.** With the tilt removed, `trend_normalized`'s `pls_0` is still the composition axis: rho(ribosomal fraction) -0.81 → -0.81, rho(RNA TMM) -0.73 → -0.67, rho(ATAC TSS enrichment) -0.77 → -0.77, rho(mitochondrial fraction) 0.73 → 0.71, rho(stage) 0.49 → 0.53 in cohort_30 (cohort_50: -0.82 → -0.81, -0.72 → -0.66, -0.77 → -0.78, 0.77 → 0.76, 0.50 → 0.55), with SD 35.0 → 32.2. The per-donor mean z over the selected RNA block, which tracked the primary `pls_0` at -0.94, still tracks it at -0.72 after normalization. The proposal's check showed that the tilt predicted most of the selected set's shift; this run shows that removing the tilt leaves the axis in place, so the shift that remains, and the composition axis with it, is carried by something the abundance trend does not describe (the proposal recorded one candidate, the lncRNA-versus-coding offset at matched abundance; it is not tested here). Because `trend_normalized` changes one setting, this attribution is to the normalization alone; it is a measurement-space description, and the statistics of each row are read under the rules above regardless of it.
- **The composed variant.** `revised_candidate_trend` behaves like the archived `revised_candidate`: the residualization removes the composition correlations of the first axis (cohort_30 `pls_0` ribosomal -0.07, TMM -0.16, TSS -0.14, stage 0.64; cohort_50 -0.07, -0.21, -0.19, 0.63), against the archived candidate's -0.06, -0.19, —, 0.61 (cohort_30). In cohort_30 it selects rank 3, and its second component carries the RNA TMM factor (`pls_1` -0.32, stage -0.44) and 0.699 of the sex-difference geometry. The held-out `cv_mean_auroc` is 0.4946 / 0.5461 against 0.6313 / 0.6616 in the primary and 0.6561 / 0.6529 for `trend_normalized`: once the technical set is regressed out, held-out stage discrimination is again near chance, as in the archived residualized variants, while the normalization alone leaves it where it was. This is a measurement-space observation reported beside the adoption check; the rules read no CV AUROC.
- **Adoption and precedence.** `trend_normalized` does not meet the adoption rule: it fails in both cohorts on six of eight check variables (`pls_0` × `tmm_rna` (-0.67), `pls_0` × `tmm_atac` (0.51), `pls_0` × `fraction_mitochondrial_umis_mean` (0.73), `pls_0` × `gex_reads_mapped_confidently_to_intronic_regions_mean` (-0.56), `pls_0` × `gex_fraction_of_transcriptomic_reads_in_cells_mean` (-0.63), `pls_0` × `rin` (-0.33) in cohort_30). `revised_candidate_trend` fails in cohort_30 (`pls_1` × `tmm_rna` (-0.32)) and holds in cohort_50, so it does not meet the rule either. Under `ADOPTION_PRECEDENCE` neither candidate is named: no revised-primary change follows from this report, both variants stay sensitivity analyses, and the primary preprocessing stays the primary. The archived `revised_candidate` failed in the other cohort (it held in cohort_30 and failed in cohort_50), so each residualized candidate fails by one cohort.
- **What follows.** Since no candidate is named, the variants stay sensitivity analyses. Under the change's plan the next predeclared variant to consider on the sead_sexdiff side is a per-gene-class normalization of the lncRNA-versus-coding offset at matched abundance, in a new proposal written before any run. This run is a reason for caution there: the abundance tilt, which the label-free check found to explain most of the selected set's shift, does not carry the first axis once it is removed, so such a proposal should first state why the gene-class offset would. The MOTCO-side covariate-aware latent construction, proposed in the MOTCO repository, is unaffected by this result.

## 6. Limitations

- **Trend normalization (`trend_normalized`, `revised_candidate_trend`).** The cyclic loess assumes that most features do not change at any abundance level, a stronger assumption than TMM's, so a stage or sex effect shaped like a tilt would be removed with the artefact; the curve is least determined at the rare end, where the selection draws most features (the ATAC lowest-bin medians above keep their spread); the lncRNA-versus-coding offset at matched abundance is not corrected; and there is no Type I check for the trend-normalized blocks.
- **Residualization (`revised_candidate_trend`).** Covariates ['log10_lib_rna', 'log10_lib_atac', 'ribosomal_pseudobulk', 'mitochondrial_pseudobulk']: the fit spends 5 degrees of freedom per feature before the permutation null, which does not see them; it removes whatever part of the stage and sex signal the technical set carries; the library-size proxies are adjusted twice (log10 library sizes in the design, log10 nuclei counts in the RRPP model), after the loess, so the composed blocks pass through three adjustments in sequence; and there is no Type I check.
- **Gene-set exclusion and selection rules (`revised_candidate_trend`).** As in the feature-selection report: the excluded sets remove measured proxies of the composition effect, not the effect on the remaining genes; the trend ratio can favour the high-mean tail and mean accessibility ranks peaks by abundance alone.
- **No Type I or power claim.** MOTCO's operating characteristics were measured on the primary-style blocks only; neither variant has a Type I check, and the adoption rule is a technical-correlation criterion, not a validation of a variant's inference.
- **Descriptive correlations at n = 80 and 74.** No p-value is attached to any correlation of the attribution, tilt or adoption tables; the two cohorts share 74 donors and are not independent replicates; `cohort_50` is a sensitivity cohort.
- **Double role of the technical table.** The gene-set fractions are regressors of the residualization and the QC means are check variables; both come from the same astrocyte nuclei and libraries.
- **Machine and archived rows.** The two variants ran on `nixos-dell-laptop`, whose inputs are hash-identical to those of the archived runs; the six archived variants were not rerun, and their rows are pasted from the feature-selection report, which has no TSS-enrichment column.


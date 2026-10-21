# MOTCO astrocyte latent-adjustment sensitivity: sex × ADNC stage — 2026-10-21

> Covariate-set and latent-basis sensitivity of the covariate-adjusted comparison of female and male astrocyte trajectories across ADNC stages in the SEA-AD MTG (OpenSpec change `mtg-astrocyte-latent-adjustment-sensitivity`, spec `motco-latent-adjustment-sensitivity`). In both cohorts the primary analysis's single CV-selected latent component is a donor-level library-quality axis (`reports/motco-astrocyte-technical-axis-2026-10-19.md`); the analysis covariates adjust the stage means within that basis but never touch which basis PLS-DA chooses. MOTCO 0.9.0's `covariate_latent_adjustment` fits the basis on features residualized on the declared covariates conditional on stage. Per cohort, two predeclared arms ran on the unchanged primary export beside the primary analysis: `extended_covariates` (the six analysis covariates plus the ribosomal and mitochondrial pseudobulk fractions, covariate-blind basis) and `latent_adjusted` (the same eight covariates, basis conditioned on `stage`). Every number below is pasted by script from the generated files (`latent_sensitivity.md`, `latent_sensitivity.csv`, `mediator_check.csv`, the comparison and analysis manifests, `result.json`, and `axis_correlations.csv` / `geometry_shares.csv` of `axis-diagnostics` on the primary and on each arm) and was checked against the six `result.json` records before commit; MOTCO's operating characteristics are cited from MOTCO's own reports. The hand-written text carries only the reading, under the rules below.

**Reading rules (binding).** The primary analysis remains the primary result, and MOTCO's power comparison (`candidate_fails`) excludes the latent adjustment as its replacement; both arms are sensitivity to it. `cohort_30` is primary and `cohort_50` sensitivity. A rejection status that differs between `primary` and `extended_covariates` is covariate-set dependence, one that differs between `extended_covariates` and `latent_adjusted` is basis dependence; neither is a result of the arm. A CV-selected rank that differs is reported beside `angle`, not read as better or worse. `angle` and `shape` non-rejections are uninformative and a `delta` rejection is a difference, not specifically of size, in every row. A technical correlation of an arm's latent axis is a measurement-space finding only, never a reinterpretation of any arm's statistics and never grounds to adopt an arm.

## 1. Provenance

- Software: MOTCO 0.9.0 (git dependency pinned to `v0.9.0`, commit `0150adc27b48528f1999a31a146a3df17ae3047f` in `uv.lock`), scikit-learn 1.8.0, numpy 2.3.5, scipy 1.15.3, pandas 2.3.3, Python 3.12.15 (the Nix shell); sead-sexdiff revision `ef388f7c13feb79c2ebea503cbfe235b9ca753c1` (dirty: True; the implementation of this change was uncommitted at run time; the reproducibility rerun below reproduced the outputs byte for byte).
- Machine: `nixos-office` (NixOS workstation), 2026-10-21. The pseudobulk outputs (with the technical table), the primary export and the primary analysis were pulled from the cluster (`/home1/tgonza/sead_sexdiff/results/`; `scripts/sync_pseudobulk.sh` and `rsync`); all 39 files hashed identically to the cluster's before the run. The cluster venv runs Python 3.13, against the 3.12 the kept primary recorded, so the run was made on the workstation.
- Settings shared by every row: seed 20261019, 999 RRPP permutations (the same permutation indices in every row), α = 0.05 (p <= alpha), `n_jobs = 1`; blocks `rna` and `atac` in identity units, stage order 0, 1, 2, net-displacement orientation, CV rank selection at MOTCO's defaults, missing-covariate policy `error` (nothing imputed).
- Covariates: primary `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`; both arms `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` (`LATENT_SENSITIVITY_COVARIATES`), the two added fractions from `donor_technical.csv` on their raw fraction scale (as in the `residualized` variant). Conditioning of `latent_adjusted`: `stage` (`group_stage` is not run).
- Inputs: pseudobulk manifest sha256 `2e93cc0b526856fc043f908a424a11266187146636bc75afed7edd187989a25c`, technical table `eac1b9a15977dfab683738a4bd0b0ffc207b896a102032409b60fccbfb96247b` and its manifest `c2bdb34499644580f22ae7f6c48a7f846c4dabad5b2b04182f8e9d5a98a106e7` (checked against the pseudobulk manifest, `TECHNICAL_GENE_SETS` and completeness for every cohort donor); the same values as in `reports/motco-astrocyte-trend-normalization-2026-10-20.md`.
- cohort_30: primary export `results/motco_input/cohort_30` (manifest sha256 `18f2289e0d4c84ed18f4a437ae569052544e9d150ed2726385273694b4a3b0b6`, 80 donors); reference `results/motco_analysis/cohort_30` (`result.json` sha256 `1d6154b02d7704ef9313ce752035a124cf6cb3fbd5350f0428439a542d5190cd`, CV-selected rank 1, cv_mean_auroc 0.6313, MOTCO 0.9.0, seed 20261019, 999 permutations); read and checked, never rerun or modified after the pin-bump check.
- cohort_50: primary export `results/motco_input/cohort_50` (manifest sha256 `dd14a6e7d5583fb527e8a99c7db43d12cf237f5451522e8875c5a0a3401c5cb3`, 74 donors); reference `results/motco_analysis/cohort_50` (`result.json` sha256 `4c6e7e58e39d8ce82c9b93907b70b2110fc2c415e32da40fef8668b74faff9fb`, CV-selected rank 1, cv_mean_auroc 0.6616, MOTCO 0.9.0, seed 20261019, 999 permutations); read and checked, never rerun or modified after the pin-bump check.
- Pin-bump check (design D7, amended 2026-10-21 by decision of the user): `analyze --cohort all --force` under MOTCO 0.9.0 with this change's code, and the pre-change code (`git archive` of `ef388f7`, its own `uv sync --locked`, MOTCO 0.8.0) on the same machine into a scratch directory. The gate passed in both cohorts:
  - cohort_30: `latent_scores.csv`, `pls_loadings.csv`, `stage_means.csv` and `analysis_metadata.csv` byte-identical between the two reruns; `result.json` and `analysis_manifest.json` equal apart from the runtime, timestamps, git, the MOTCO version and the `cli_equivalent` output path.
  - cohort_50: `latent_scores.csv`, `pls_loadings.csv`, `stage_means.csv` and `analysis_metadata.csv` byte-identical between the two reruns; `result.json` and `analysis_manifest.json` equal apart from the runtime, timestamps, git, the MOTCO version and the `cli_equivalent` output path.
- The 0.9.0 rerun was also compared with the kept 0.8.0 primary, which was made on the laptop workstation at `c2f25e3` (Python 3.12.15). This comparison is recorded, not a gate:
  - cohort_30: `latent_scores.csv` max abs difference 2.63e-13 (largest |value| 87.2), `pls_loadings.csv` max abs difference 1.01e-16 (largest |value| 0.0327), `stage_means.csv` max abs difference 4.26e-14 (largest |value| 24.1), `analysis_metadata.csv` byte-identical; statistics equal to 4.6e-15 relative and null summaries to 1e-14 relative; p-values and the CV-selected rank identical; keys only in the rerun's manifest `preprocessing.excluded_gene_sets`, `preprocessing.normalization`, `preprocessing.normalization_fit`, `preprocessing.residualization`, `preprocessing.selection_rule` (added by later archived sead-sexdiff changes). `result.json` sha256 `ca92876a1a61c40791524390e8603125bcbfacf6f081d34e2ae6a9b49fd15de6` (kept) → `1d6154b02d7704ef9313ce752035a124cf6cb3fbd5350f0428439a542d5190cd` (reference of the arms); `analysis_manifest.json` `17b2150ecba0285ecb8cadf846cda1b7426e232c615003561f8f1934d24b7f85` → `dc126b8a3b33a55bbab6bf4e4d46c5815ae1666122078bdfb1e883adbe2ac86d`.
  - cohort_50: `latent_scores.csv` max abs difference 8.17e-14 (largest |value| 84.6), `pls_loadings.csv` max abs difference 1.01e-16 (largest |value| 0.0324), `stage_means.csv` max abs difference 4.26e-14 (largest |value| 23), `analysis_metadata.csv` byte-identical; statistics equal to 1.2e-15 relative and null summaries to 1.3e-14 relative; p-values and the CV-selected rank identical; keys only in the rerun's manifest `preprocessing.excluded_gene_sets`, `preprocessing.normalization`, `preprocessing.normalization_fit`, `preprocessing.residualization`, `preprocessing.selection_rule` (added by later archived sead-sexdiff changes). `result.json` sha256 `1e7eefe3e63e7b9b51c3fed618c32d4a496ea3d382c7ff025c2a34f5d0ae048f` (kept) → `4c6e7e58e39d8ce82c9b93907b70b2110fc2c415e32da40fef8668b74faff9fb` (reference of the arms); `analysis_manifest.json` `43fe406690fa112a0c674c67c9e9440f17a5736618ffb38a8099fbd0b18679d3` → `7cfc9c5d820f3caa566f67f300d2fbad711a64aa6aaaae9c1fa71c0748081f9a`.
  - The rounding differences come from the machine, not from MOTCO (the same-machine gate above is byte-identical). The kept directories remain as `results/motco_analysis_0.8.0/` on the cluster, where the earlier reports' hashes stay resolvable; the cluster's `results/motco_analysis/` now holds the 0.9.0 rerun.
- Reproducibility: `latent-sensitivity --cohort cohort_50` rerun into a scratch directory with the same inputs, MOTCO release, seed and `n_jobs` gave byte-identical `latent_sensitivity.csv`, `mediator_check.csv` and arm CSVs (all 10 files compared).
- Diagnostics: `sead-sexdiff axis-diagnostics --analysis-dir <root> --rank-dir <absent> --out-dir results/motco_latent_sensitivity/diagnostics/<name>` on the primary and on each arm, with the primary export as input; every analysis manifest's export hashes matched the files on disk. The diagnostic labels the analysis it reads `primary`; the tables below name the analysis. Its partial correlations condition on the six analysis covariates (its definition), also for the arms. No rank sweep exists for the 0.9.0 reference or the arms (recorded by the diagnostic as a finding).

## 2. Primary cohort: cohort_30 (≥ 30 astrocyte nuclei in both omics, n = 80)

### Overview (generated from `latent_sensitivity.csv`)

| analysis | covariates | conditioning | rank | cv_mean_auroc | delta (p) | angle (p) | shape (p) | rejected |
|---|---:|---|---:|---:|---:|---:|---:|---|
| primary | 6 | none | 1 | 0.6313 | 14.941 (0.3) | 0 (1) | 0.608832 (0.333) | none |
| extended_covariates | 8 | none | 1 | 0.6313 | 3.93908 (0.641) | 0 (1) | 1.01735 (0.127) | none |
| latent_adjusted | 8 | stage | 1 | 0.6299 | 1.3709 (0.717) | 0 (1) | 0.16446 (0.327) | none |

### Mediator check (pasted from `mediator_check.csv`, written before the arms ran)

| covariate | subset | n | Spearman | flagged |
|---|---|---:|---:|:---:|
| `ribosomal_pseudobulk` | all | 80 | -0.399 | yes |
| `ribosomal_pseudobulk` | F | 49 | -0.432 | yes |
| `ribosomal_pseudobulk` | M | 31 | -0.285 | no |
| `mitochondrial_pseudobulk` | all | 80 | 0.312 | yes |
| `mitochondrial_pseudobulk` | F | 49 | 0.342 | yes |
| `mitochondrial_pseudobulk` | M | 31 | 0.189 | no |

### First latent axis per analysis (generated from `axis_correlations.csv`)

Spearman correlation of `pls_0` with the technical-axis report's headline variables, stage and sex, raw and (in parentheses) partial given the six analysis covariates; n/a where a variable is one of them. No p-value is attached.

| variable | primary | extended_covariates | latent_adjusted |
|---|---:|---:|---:|
| `ribosomal_pseudobulk` | -0.81 (-0.79) | -0.81 (-0.79) | -0.15 (-0.05) |
| `ribosomal_mean` | -0.84 (-0.82) | -0.84 (-0.82) | -0.15 (-0.04) |
| `mitochondrial_pseudobulk` | 0.73 (0.56) | 0.73 (0.56) | 0.12 (-0.03) |
| `fraction_mitochondrial_umis_mean` | 0.75 (0.59) | 0.75 (0.59) | 0.12 (-0.03) |
| `atac_tss_enrichment_score_mean` | -0.77 (-0.67) | -0.77 (-0.67) | -0.08 (-0.07) |
| `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | -0.77 (-0.64) | -0.77 (-0.64) | -0.06 (0.01) |
| `tmm_rna` | -0.73 (-0.67) | -0.73 (-0.67) | -0.24 (-0.18) |
| `tmm_atac` | 0.49 (0.50) | 0.49 (0.50) | -0.07 (0.01) |
| `log10_lib_rna` | 0.27 (0.08) | 0.27 (0.08) | 0.13 (0.25) |
| `log10_lib_atac` | 0.21 (-0.00) | 0.21 (-0.00) | 0.38 (0.30) |
| `rin` | -0.31 (n/a) | -0.31 (n/a) | 0.04 (n/a) |
| `pmi` | 0.05 (n/a) | 0.05 (n/a) | 0.06 (n/a) |
| `stage` | 0.49 (0.41) | 0.49 (0.41) | -0.64 (-0.77) |
| `male` | -0.23 (-0.33) | -0.23 (-0.33) | 0.06 (0.03) |

- `primary` `pls_0` (score SD 34.99): 20 of 47 variables flagged at |ρ| ≥ 0.3 (raw or partial): `ribosomal_mean` -0.84, `ribosomal_pseudobulk` -0.81, `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` -0.77, `atac_tss_enrichment_score_mean` -0.77, `fraction_mitochondrial_umis_mean` 0.75, `mitochondrial_mean` 0.75, `tmm_rna` -0.73, `mitochondrial_pseudobulk` 0.73, `ribosomal_sd` -0.70, `gex_sequencing_saturation_mean` 0.67, `fraction_mitochondrial_umis_sd` 0.65, `mitochondrial_sd` 0.65, `gex_fraction_of_transcriptomic_reads_in_cells_mean` -0.61, `atac_fraction_of_high_quality_fragments_in_cells_mean` -0.55, `gex_reads_mapped_confidently_to_intronic_regions_mean` -0.54, `stage` 0.49, `tmm_atac` 0.49, `gex_median_umi_counts_per_cell_mean` -0.47, `rin` -0.31, `male` -0.23. Geometry share along it: net_displacement_Female 1.000, net_displacement_Male 1.000, sex_difference_by_stage 1.000, sex_difference_centred 1.000.
- `extended_covariates` `pls_0` (score SD 34.99): 20 of 47 variables flagged at |ρ| ≥ 0.3 (raw or partial): `ribosomal_mean` -0.84, `ribosomal_pseudobulk` -0.81, `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` -0.77, `atac_tss_enrichment_score_mean` -0.77, `fraction_mitochondrial_umis_mean` 0.75, `mitochondrial_mean` 0.75, `tmm_rna` -0.73, `mitochondrial_pseudobulk` 0.73, `ribosomal_sd` -0.70, `gex_sequencing_saturation_mean` 0.67, `fraction_mitochondrial_umis_sd` 0.65, `mitochondrial_sd` 0.65, `gex_fraction_of_transcriptomic_reads_in_cells_mean` -0.61, `atac_fraction_of_high_quality_fragments_in_cells_mean` -0.55, `gex_reads_mapped_confidently_to_intronic_regions_mean` -0.54, `stage` 0.49, `tmm_atac` 0.49, `gex_median_umi_counts_per_cell_mean` -0.47, `rin` -0.31, `male` -0.23. Geometry share along it: net_displacement_Female 1.000, net_displacement_Male 1.000, sex_difference_by_stage 1.000, sex_difference_centred 1.000.
- `latent_adjusted` `pls_0` (score SD 13.53): 3 of 47 variables flagged at |ρ| ≥ 0.3 (raw or partial): `stage` -0.64, `log10_lib_atac` 0.38, `log10_n_atac` 0.31. Geometry share along it: net_displacement_Female 1.000, net_displacement_Male 1.000, sex_difference_by_stage 1.000, sex_difference_centred 1.000.

### Comparison (pasted from `results/motco_latent_sensitivity/comparison/cohort_30/latent_sensitivity.md`)

#### MOTCO latent-adjustment sensitivity: cohort_30

Covariate-set and latent-basis sensitivity of the covariate-adjusted trajectory analysis. The `primary` row is the primary analysis (CV-selected rank 1, the six analysis covariates, covariate-blind basis), read from its `result.json` and never rerun; it stays the primary result, and MOTCO's power comparison (`candidate_fails`) excludes the latent adjustment as its replacement. `extended_covariates` adds `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` (raw fractions from `donor_technical.csv`) with the covariate-blind basis; `latent_adjusted` enters the same eight covariates and fits the basis on features residualized on them conditional on `stage`, with the trajectory model, LS-means and RRPP of `extended_covariates`, so those two rows differ through the basis and the CV-selected rank only. Both arms are read as sensitivity, never as the result: a rejection status that differs between `primary` and `extended_covariates` is covariate-set dependence, one that differs between `extended_covariates` and `latent_adjusted` is basis dependence; a CV-selected rank that differs is reported beside `angle`, not read as better or worse; `angle` and `shape` non-rejections are uninformative and a `delta` rejection is a difference, not specifically of size, in every row. 999 RRPP permutations with seed 20261019 in every row (the same permutation indices); α = 0.05 (rejected when p ≤ α). MOTCO measured the latent adjustment's operating characteristics at its Phase 6 design only (Type I `adopt_stage`; power `candidate_fails`, with `shape` power lost under a nuisance axis and small-magnitude power lost without one), on neither these blocks nor this eight-covariate set, so neither arm carries a Type I or power claim of its own.

#### delta

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 1 | 14.941 | 0.3 | no | 11.2868 | 8.27982 | 26.5775 |
| extended_covariates | 8 | none | 1 | 3.93908 | 0.641 | no | 6.98845 | 5.36768 | 17.2685 |
| latent_adjusted | 8 | stage | 1 | 1.3709 | 0.717 | no | 2.98333 | 2.23983 | 7.22261 |

#### angle (degrees)

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 1 | 0 | 1 | no | 3.06306 | 23.2919 | 0 |
| extended_covariates | 8 | none | 1 | 0 | 1 | no | 4.32432 | 27.5761 | 0 |
| latent_adjusted | 8 | stage | 1 | 0 | 1 | no | 0 | 0 | 0 |

#### shape

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 1 | 0.608832 | 0.333 | no | 0.489849 | 0.355498 | 1.18287 |
| extended_covariates | 8 | none | 1 | 1.01735 | 0.127 | no | 0.54722 | 0.36351 | 1.2249 |
| latent_adjusted | 8 | stage | 1 | 0.16446 | 0.327 | no | 0.132841 | 0.098463 | 0.328221 |

#### Runs

| arm | selection | rank | cv_mean_auroc | covariates (in order) | conditioning | runtime (s) | directory |
|---|---|---:|---:|---|---|---:|---|
| primary | cv | 1 | 0.6313 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac` | none | 61.9 | `results/motco_analysis/cohort_30` |
| extended_covariates | cv | 1 | 0.6313 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` | none | 61.1 | `extended_covariates/motco_analysis/cohort_30` |
| latent_adjusted | cv | 1 | 0.6299 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` | stage | 61.3 | `latent_adjusted/motco_analysis/cohort_30` |

#### Latent adjustment

MOTCO's `latent_adjustment` block: the residualization's degrees of freedom on the full sample (`n − 1 − k − q`) and the smallest over the CV training sets, and the quantiles over features of the partial R² of the covariates given the intercept and the conditioning indicators (pooled-standardized features).

| arm | conditioning | residual df | smallest CV training df | R² q05 | R² q25 | R² q50 | R² q75 | R² q95 | R² mean | R² fraction > 0.5 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| latent_adjusted | stage | 69 | 29 | 0.0862 | 0.1547 | 0.2230 | 0.3192 | 0.5743 | 0.2578 | 0.0907 |

#### Mediator check

Spearman correlation of each added covariate with stage (ordinal 0, 1, 2) over all cohort donors and within each sex, written before the arms ran; |ρ| ≥ 0.3 (`DIAGNOSTIC_CORRELATION_THRESHOLD`) flags a mediator candidate (flagged rows first). No p-value is computed, and a flag removes no covariate and changes no arm: no simulation can tell a confounder from a mediator, so a flagged fraction may carry stage signal that the `latent_adjusted` basis removed (and that the covariate adjustment of both arms absorbed).

| covariate | subset | n | Spearman | flagged |
|---|---|---:|---:|:---:|
| `ribosomal_pseudobulk` | all | 80 | -0.399 | yes |
| `ribosomal_pseudobulk` | F | 49 | -0.432 | yes |
| `mitochondrial_pseudobulk` | all | 80 | 0.312 | yes |
| `mitochondrial_pseudobulk` | F | 49 | 0.342 | yes |
| `ribosomal_pseudobulk` | M | 31 | -0.285 | no |
| `mitochondrial_pseudobulk` | M | 31 | 0.189 | no |

#### Findings

- mediator candidate: ribosomal_pseudobulk has Spearman -0.399 with stage (all, n = 80); the latent_adjusted basis may have removed stage signal carried by it (no simulation can tell a confounder from a mediator; the flag changes no arm)
- mediator candidate: ribosomal_pseudobulk has Spearman -0.432 with stage (F, n = 49); the latent_adjusted basis may have removed stage signal carried by it (no simulation can tell a confounder from a mediator; the flag changes no arm)
- mediator candidate: mitochondrial_pseudobulk has Spearman 0.312 with stage (all, n = 80); the latent_adjusted basis may have removed stage signal carried by it (no simulation can tell a confounder from a mediator; the flag changes no arm)
- mediator candidate: mitochondrial_pseudobulk has Spearman 0.342 with stage (F, n = 49); the latent_adjusted basis may have removed stage signal carried by it (no simulation can tell a confounder from a mediator; the flag changes no arm)

## 3. Sensitivity cohort: cohort_50 (≥ 50 astrocyte nuclei in both omics, n = 74)

### Overview (generated from `latent_sensitivity.csv`)

| analysis | covariates | conditioning | rank | cv_mean_auroc | delta (p) | angle (p) | shape (p) | rejected |
|---|---:|---|---:|---:|---:|---:|---:|---|
| primary | 6 | none | 1 | 0.6616 | 18.7091 (0.213) | 0 (1) | 0.56471 (0.377) | none |
| extended_covariates | 8 | none | 1 | 0.6616 | 5.61916 (0.585) | 0 (1) | 0.697892 (0.201) | none |
| latent_adjusted | 8 | stage | 3 | 0.6645 | 4.96136 (0.494) | 24.7286 (0.154) | 0.132168 (0.337) | none |

### Mediator check (pasted from `mediator_check.csv`, written before the arms ran)

| covariate | subset | n | Spearman | flagged |
|---|---|---:|---:|:---:|
| `ribosomal_pseudobulk` | all | 74 | -0.352 | yes |
| `ribosomal_pseudobulk` | F | 44 | -0.360 | yes |
| `ribosomal_pseudobulk` | M | 30 | -0.286 | no |
| `mitochondrial_pseudobulk` | all | 74 | 0.356 | yes |
| `mitochondrial_pseudobulk` | F | 44 | 0.422 | yes |
| `mitochondrial_pseudobulk` | M | 30 | 0.189 | no |

### First latent axis per analysis (generated from `axis_correlations.csv`)

Spearman correlation of `pls_0` with the technical-axis report's headline variables, stage and sex, raw and (in parentheses) partial given the six analysis covariates; n/a where a variable is one of them. No p-value is attached.

| variable | primary | extended_covariates | latent_adjusted |
|---|---:|---:|---:|
| `ribosomal_pseudobulk` | -0.82 (-0.77) | -0.82 (-0.77) | -0.06 (0.02) |
| `ribosomal_mean` | -0.84 (-0.80) | -0.84 (-0.80) | -0.05 (0.04) |
| `mitochondrial_pseudobulk` | 0.77 (0.57) | 0.77 (0.57) | -0.03 (-0.09) |
| `fraction_mitochondrial_umis_mean` | 0.79 (0.60) | 0.79 (0.60) | -0.03 (-0.09) |
| `atac_tss_enrichment_score_mean` | -0.77 (-0.69) | -0.77 (-0.69) | 0.03 (0.02) |
| `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | -0.78 (-0.68) | -0.78 (-0.68) | 0.07 (0.10) |
| `tmm_rna` | -0.72 (-0.62) | -0.72 (-0.62) | -0.13 (-0.05) |
| `tmm_atac` | 0.56 (0.58) | 0.56 (0.58) | -0.14 (-0.11) |
| `log10_lib_rna` | 0.28 (0.02) | 0.28 (0.02) | 0.11 (0.23) |
| `log10_lib_atac` | 0.23 (-0.04) | 0.23 (-0.04) | 0.33 (0.34) |
| `rin` | -0.30 (n/a) | -0.30 (n/a) | 0.09 (n/a) |
| `pmi` | 0.10 (n/a) | 0.10 (n/a) | -0.00 (n/a) |
| `stage` | 0.50 (0.43) | 0.50 (0.43) | -0.75 (-0.80) |
| `male` | -0.19 (-0.29) | -0.19 (-0.29) | 0.07 (0.05) |

- `primary` `pls_0` (score SD 35.76): 21 of 47 variables flagged at |ρ| ≥ 0.3 (raw or partial): `ribosomal_mean` -0.84, `ribosomal_pseudobulk` -0.82, `fraction_mitochondrial_umis_mean` 0.79, `mitochondrial_mean` 0.79, `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` -0.78, `atac_tss_enrichment_score_mean` -0.77, `mitochondrial_pseudobulk` 0.77, `gex_sequencing_saturation_mean` 0.73, `tmm_rna` -0.72, `ribosomal_sd` -0.69, `fraction_mitochondrial_umis_sd` 0.68, `mitochondrial_sd` 0.68, `gex_fraction_of_transcriptomic_reads_in_cells_mean` -0.61, `atac_fraction_of_high_quality_fragments_in_cells_mean` -0.60, `gex_reads_mapped_confidently_to_intronic_regions_mean` -0.56, `tmm_atac` 0.56, `gex_median_umi_counts_per_cell_mean` -0.51, `stage` 0.50, `gex_reads_mapped_confidently_to_intronic_regions_sd` 0.33, `log10_n_atac` 0.32, `atac_percent_duplicates_mean` -0.30. Geometry share along it: net_displacement_Female 1.000, net_displacement_Male 1.000, sex_difference_by_stage 1.000, sex_difference_centred 1.000.
- `extended_covariates` `pls_0` (score SD 35.76): 21 of 47 variables flagged at |ρ| ≥ 0.3 (raw or partial): `ribosomal_mean` -0.84, `ribosomal_pseudobulk` -0.82, `fraction_mitochondrial_umis_mean` 0.79, `mitochondrial_mean` 0.79, `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` -0.78, `atac_tss_enrichment_score_mean` -0.77, `mitochondrial_pseudobulk` 0.77, `gex_sequencing_saturation_mean` 0.73, `tmm_rna` -0.72, `ribosomal_sd` -0.69, `fraction_mitochondrial_umis_sd` 0.68, `mitochondrial_sd` 0.68, `gex_fraction_of_transcriptomic_reads_in_cells_mean` -0.61, `atac_fraction_of_high_quality_fragments_in_cells_mean` -0.60, `gex_reads_mapped_confidently_to_intronic_regions_mean` -0.56, `tmm_atac` 0.56, `gex_median_umi_counts_per_cell_mean` -0.51, `stage` 0.50, `gex_reads_mapped_confidently_to_intronic_regions_sd` 0.33, `log10_n_atac` 0.32, `atac_percent_duplicates_mean` -0.30. Geometry share along it: net_displacement_Female 1.000, net_displacement_Male 1.000, sex_difference_by_stage 1.000, sex_difference_centred 1.000.
- `latent_adjusted` `pls_0` (score SD 13.85): 2 of 47 variables flagged at |ρ| ≥ 0.3 (raw or partial): `stage` -0.75, `log10_lib_atac` 0.33. Geometry share along it: net_displacement_Female 0.800, net_displacement_Male 0.613, sex_difference_by_stage 0.028, sex_difference_centred 0.022.
- `latent_adjusted` `pls_1` (score SD 17.78): 19 of 47 variables flagged at |ρ| ≥ 0.3 (raw or partial): `atac_tss_enrichment_score_mean` 0.67, `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` 0.65, `tmm_rna` 0.61, `ribosomal_mean` 0.61, `fraction_mitochondrial_umis_mean` -0.60, `mitochondrial_mean` -0.60, `ribosomal_pseudobulk` 0.60, `mitochondrial_pseudobulk` -0.59, `atac_fraction_of_high_quality_fragments_in_cells_mean` 0.56, `gex_sequencing_saturation_mean` -0.53, `fraction_mitochondrial_umis_sd` -0.53, `mitochondrial_sd` -0.53, `ribosomal_sd` 0.51, `tmm_atac` -0.46, `gex_median_umi_counts_per_cell_mean` 0.45, `gex_reads_mapped_confidently_to_intronic_regions_mean` 0.44, `atac_percent_duplicates_mean` 0.44, `gex_fraction_of_transcriptomic_reads_in_cells_mean` 0.34, `n_nuclei_rna` -0.17. Geometry share along it: net_displacement_Female 0.104, net_displacement_Male 0.002, sex_difference_by_stage 0.186, sex_difference_centred 0.234.
- `latent_adjusted` `pls_2` (score SD 18.01): 17 of 47 variables flagged at |ρ| ≥ 0.3 (raw or partial): `atac_tss_enrichment_score_mean` -0.50, `tmm_rna` -0.50, `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` -0.45, `ribosomal_mean` -0.42, `ribosomal_pseudobulk` -0.41, `atac_fraction_of_high_quality_fragments_in_cells_mean` -0.40, `ribosomal_sd` -0.36, `gex_fraction_of_transcriptomic_reads_in_cells_sd` 0.36, `gex_reads_mapped_confidently_to_intronic_regions_mean` -0.35, `gex_sequencing_saturation_mean` 0.35, `gex_fraction_of_transcriptomic_reads_in_cells_mean` -0.35, `tmm_atac` 0.35, `fraction_mitochondrial_umis_sd` 0.35, `mitochondrial_sd` 0.35, `fraction_mitochondrial_umis_mean` 0.35, `mitochondrial_mean` 0.35, `mitochondrial_pseudobulk` 0.32. Geometry share along it: net_displacement_Female 0.096, net_displacement_Male 0.384, sex_difference_by_stage 0.786, sex_difference_centred 0.744.

### Comparison (pasted from `results/motco_latent_sensitivity/comparison/cohort_50/latent_sensitivity.md`)

#### MOTCO latent-adjustment sensitivity: cohort_50

Covariate-set and latent-basis sensitivity of the covariate-adjusted trajectory analysis. The `primary` row is the primary analysis (CV-selected rank 1, the six analysis covariates, covariate-blind basis), read from its `result.json` and never rerun; it stays the primary result, and MOTCO's power comparison (`candidate_fails`) excludes the latent adjustment as its replacement. `extended_covariates` adds `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` (raw fractions from `donor_technical.csv`) with the covariate-blind basis; `latent_adjusted` enters the same eight covariates and fits the basis on features residualized on them conditional on `stage`, with the trajectory model, LS-means and RRPP of `extended_covariates`, so those two rows differ through the basis and the CV-selected rank only. Both arms are read as sensitivity, never as the result: a rejection status that differs between `primary` and `extended_covariates` is covariate-set dependence, one that differs between `extended_covariates` and `latent_adjusted` is basis dependence; a CV-selected rank that differs is reported beside `angle`, not read as better or worse; `angle` and `shape` non-rejections are uninformative and a `delta` rejection is a difference, not specifically of size, in every row. 999 RRPP permutations with seed 20261019 in every row (the same permutation indices); α = 0.05 (rejected when p ≤ α). MOTCO measured the latent adjustment's operating characteristics at its Phase 6 design only (Type I `adopt_stage`; power `candidate_fails`, with `shape` power lost under a nuisance axis and small-magnitude power lost without one), on neither these blocks nor this eight-covariate set, so neither arm carries a Type I or power claim of its own.

#### delta

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 1 | 18.7091 | 0.213 | no | 11.6918 | 9.218 | 30.0679 |
| extended_covariates | 8 | none | 1 | 5.61916 | 0.585 | no | 7.90176 | 6.05969 | 19.7526 |
| latent_adjusted | 8 | stage | 3 | 4.96136 | 0.494 | no | 5.99914 | 4.52218 | 14.382 |

#### angle (degrees)

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 1 | 0 | 1 | no | 1.98198 | 18.7931 | 0 |
| extended_covariates | 8 | none | 1 | 0 | 1 | no | 0.720721 | 11.3728 | 0 |
| latent_adjusted | 8 | stage | 3 | 24.7286 | 0.154 | no | 15.5719 | 9.05433 | 32.1693 |

#### shape

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 1 | 0.56471 | 0.377 | no | 0.495877 | 0.343673 | 1.159 |
| extended_covariates | 8 | none | 1 | 0.697892 | 0.201 | no | 0.420275 | 0.319695 | 1.03213 |
| latent_adjusted | 8 | stage | 3 | 0.132168 | 0.337 | no | 0.110548 | 0.0582071 | 0.219488 |

#### Runs

| arm | selection | rank | cv_mean_auroc | covariates (in order) | conditioning | runtime (s) | directory |
|---|---|---:|---:|---|---|---:|---|
| primary | cv | 1 | 0.6616 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac` | none | 62.4 | `results/motco_analysis/cohort_50` |
| extended_covariates | cv | 1 | 0.6616 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` | none | 62.6 | `extended_covariates/motco_analysis/cohort_50` |
| latent_adjusted | cv | 3 | 0.6645 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` | stage | 61.4 | `latent_adjusted/motco_analysis/cohort_50` |

#### Latent adjustment

MOTCO's `latent_adjustment` block: the residualization's degrees of freedom on the full sample (`n − 1 − k − q`) and the smallest over the CV training sets, and the quantiles over features of the partial R² of the covariates given the intercept and the conditioning indicators (pooled-standardized features).

| arm | conditioning | residual df | smallest CV training df | R² q05 | R² q25 | R² q50 | R² q75 | R² q95 | R² mean | R² fraction > 0.5 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| latent_adjusted | stage | 63 | 25 | 0.0906 | 0.1613 | 0.2350 | 0.3406 | 0.6031 | 0.2721 | 0.1082 |

#### Mediator check

Spearman correlation of each added covariate with stage (ordinal 0, 1, 2) over all cohort donors and within each sex, written before the arms ran; |ρ| ≥ 0.3 (`DIAGNOSTIC_CORRELATION_THRESHOLD`) flags a mediator candidate (flagged rows first). No p-value is computed, and a flag removes no covariate and changes no arm: no simulation can tell a confounder from a mediator, so a flagged fraction may carry stage signal that the `latent_adjusted` basis removed (and that the covariate adjustment of both arms absorbed).

| covariate | subset | n | Spearman | flagged |
|---|---|---:|---:|:---:|
| `ribosomal_pseudobulk` | all | 74 | -0.352 | yes |
| `ribosomal_pseudobulk` | F | 44 | -0.360 | yes |
| `mitochondrial_pseudobulk` | all | 74 | 0.356 | yes |
| `mitochondrial_pseudobulk` | F | 44 | 0.422 | yes |
| `ribosomal_pseudobulk` | M | 30 | -0.286 | no |
| `mitochondrial_pseudobulk` | M | 30 | 0.189 | no |

#### Findings

- CV-selected rank differs: 1 in extended_covariates, 3 in latent_adjusted; reported beside angle, not read as better or worse
- mediator candidate: ribosomal_pseudobulk has Spearman -0.352 with stage (all, n = 74); the latent_adjusted basis may have removed stage signal carried by it (no simulation can tell a confounder from a mediator; the flag changes no arm)
- mediator candidate: ribosomal_pseudobulk has Spearman -0.360 with stage (F, n = 44); the latent_adjusted basis may have removed stage signal carried by it (no simulation can tell a confounder from a mediator; the flag changes no arm)
- mediator candidate: mitochondrial_pseudobulk has Spearman 0.356 with stage (all, n = 74); the latent_adjusted basis may have removed stage signal carried by it (no simulation can tell a confounder from a mediator; the flag changes no arm)
- mediator candidate: mitochondrial_pseudobulk has Spearman 0.422 with stage (F, n = 44); the latent_adjusted basis may have removed stage signal carried by it (no simulation can tell a confounder from a mediator; the flag changes no arm)

## 4. Reading

- **No rejection in any row.** No statistic is rejected at α = 0.05 in the primary analysis or in either arm, in either cohort, so no rejection status differs between `primary` and `extended_covariates` (no covariate-set dependence of a rejection) or between `extended_covariates` and `latent_adjusted` (no basis dependence of a rejection). The non-rejections are uninformative: `angle` and `shape` non-rejections are not evidence that the sexes agree in direction or shape, and a `delta` non-rejection is not evidence against a size difference, in any row.
- **`extended_covariates` keeps the primary's basis.** The basis is fitted without covariates in both, so the CV-selected rank (1 and 1), the CV AUROC (0.6313, 0.6616) and the first axis's correlations (`ribosomal_mean` -0.84 and -0.84) equal the primary's. What changes is the trajectory model: with the two fractions entered, `delta` is 3.93908 against 14.941 (cohort_30) and 5.61916 against 18.7091 (cohort_50), and `shape` 1.01735 against 0.608832 and 0.697892 against 0.56471, with the null summaries moving with them. These are covariate-set dependence of the magnitudes, not of a rejection status, and they cannot say whether the fractions removed a confounded or a mediated part of the trajectories (section 5).
- **`latent_adjusted` moves the first axis off the library-quality variables.** Its `pls_0` correlates at most 0.38 (`log10_lib_atac`) with a technical variable in cohort_30 and 0.33 (`log10_lib_atac`) in cohort_50, against `ribosomal_mean` -0.84 / -0.84 for the primary, and -0.64 / -0.75 with stage, against 0.49 / 0.50. This is a measurement-space finding about the arm: it does not reinterpret the primary's statistics, and it is not grounds to adopt the arm.
- **cohort_50 selects rank 3 under the latent adjustment.** The CV-selected rank differs from `extended_covariates` (1), so `angle` is no longer 0° by construction: 24.7286° (p 0.154), reported beside rank 3, not read as better or worse, and its non-rejection is uninformative. Its second and third components carry the library-quality variables (`ribosomal_mean` 0.61 on `pls_1`, `tmm_rna` 0.61 on `pls_1` and -0.50 on `pls_2`): MOTCO measures the raw features in the adjusted basis, so a component can still track a technical variable, which the trajectory model's eight covariates then adjust for. In cohort_30 the rank stays 1 and `angle` is 0° by construction in all three rows.
- **Mediator candidates in both cohorts.** The mediator check flags `ribosomal_pseudobulk` all -0.399 (n = 80), `ribosomal_pseudobulk` F -0.432 (n = 49), `mitochondrial_pseudobulk` all 0.312 (n = 80), `mitochondrial_pseudobulk` F 0.342 (n = 49) in cohort_30 and `ribosomal_pseudobulk` all -0.352 (n = 74), `ribosomal_pseudobulk` F -0.360 (n = 44), `mitochondrial_pseudobulk` all 0.356 (n = 74), `mitochondrial_pseudobulk` F 0.422 (n = 44) in cohort_50; the male cells are not flagged (`ribosomal_pseudobulk` M -0.285, `mitochondrial_pseudobulk` M 0.189; `ribosomal_pseudobulk` M -0.286, `mitochondrial_pseudobulk` M 0.189). Both fractions move with stage, mostly in the female donors, so the `latent_adjusted` basis, and the covariate adjustment of both arms, may have removed stage signal carried by them. No simulation can tell a confounder from a mediator, and the flags changed no arm.

## 5. Limitations

- **Mediators.** A residualized variable that mediates the stage effect takes real signal with it. The ribosomal and mitochondrial fractions are flagged in both cohorts (section 4); a pathology-driven change in astrocyte RNA quality would make them mediators rather than confounders, and the data cannot tell the two apart (MOTCO `docs/reports/latent-adjustment-type-i-2026-10-20.md` §4, `latent-adjustment-power-2026-10-20.md` §6).
- **MOTCO's operating characteristics were measured at its Phase 6 design only.** Type I: `adopt_stage`; with `stage` conditioning `delta` 0.031, `angle` 0.067 and `shape` 0.058 against the bound 0.0706 (R = 448), `group_stage` exceeded on `angle` (Type I report, summary). Power: `candidate_fails` for `stage` (power report, summary): it gains orientation power (`angle` 0.922 against 0.784 with a covariate-aligned nuisance axis, 0.931 against 0.642 without), loses `shape` power with the axis (0.243 against 0.428), loses small-magnitude power without it (0.400 against 0.520 at magnitude e = 0.1), and its translation control rejects `angle` at 0.085–0.092, above the bound (the covariate-blind basis at 0.076–0.086). So a `delta` or `shape` non-rejection under the latent adjustment is weaker than under the covariate-blind basis, and an `angle` rejection seen only under it would have to be read against the translation control (none occurred here).
- **No operating characteristic on these blocks or at this covariate set.** MOTCO's checks used semisynthetic blocks of 498 features with two injected and six pure-noise covariates (eight declared); here 10,000 exported features and eight real covariates, two of them stage-associated. Neither arm carries a Type I or power claim of its own, and the `extended_covariates` adjustment with eight covariates has no Type I check beyond MOTCO's.
- **Degrees of freedom (cohort_30).** The residualization spends 11 columns (intercept, stage indicators, eight covariates) for 80 donors: residual df 69, smallest CV training df 29; the covariates' partial R² has median 0.2230 and exceeds 0.5 for a fraction 0.0907 of the features.
- **Degrees of freedom (cohort_50).** The residualization spends 11 columns (intercept, stage indicators, eight covariates) for 74 donors: residual df 63, smallest CV training df 25; the covariates' partial R² has median 0.2350 and exceeds 0.5 for a fraction 0.1082 of the features.
- **Library quality enters twice.** The nuclei counts (`log10_n_rna`, `log10_n_atac`) and the two fractions both proxy library quality; the adjustment is joint, so the collinear proxies do not change what is adjusted for, but the arms adjust for library quality through more than one proxy.
- **Rank above `n_stages − 1`.** cohort_50's `latent_adjusted` rank exceeds 2; MOTCO's latent-rank ladder found a wider `angle` null above `n_stages − 1` (unadjusted, at a larger n), cited as what it is.
- **Reference rerun under the new pin.** The arms' reference is the primary rerun under MOTCO 0.9.0 on this machine, byte-identical to MOTCO 0.8.0 on the same machine and equal at rounding level to the kept laptop primary (section 1); the earlier reports' numbers and reading are unchanged at the precision they print.
- **Uncommitted implementation at run time** (`dirty: True`); the cohort_50 rerun reproduced the outputs byte for byte, and a post-commit rerun will carry a clean `git` field.
- **Everything the primary report lists** (additive linear adjustment, structural analogue for the operating characteristics, fixed preprocessing choices, Monte Carlo error and `n_jobs`) applies unchanged.

## 6. Next steps

- None follows from these results for the primary analysis: it stays primary. A revised primary, `group_stage` conditioning, other covariate sets, latent adjustment of a preprocessing variant, and MOTCO's block weighting and rank-1 `angle` degeneracy (roadmap items 4–5) are out of scope and would each need their own predeclared change.

## 7. Reproduction

```
scripts/sync_pseudobulk.sh                                    # results/pseudobulk/ from the cluster
rsync -a <cluster>:/home1/tgonza/sead_sexdiff/results/{motco_input,motco_analysis} results/
uv run sead-sexdiff latent-sensitivity --cohort all
for a in extended_covariates latent_adjusted; do
  uv run sead-sexdiff axis-diagnostics --analysis-dir results/motco_latent_sensitivity/$a/motco_analysis \
      --out-dir results/motco_latent_sensitivity/diagnostics/$a
done
```

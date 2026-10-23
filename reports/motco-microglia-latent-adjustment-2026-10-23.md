# MOTCO microglia latent-adjustment sensitivity and the restated positive control: sex × ADNC stage — 2026-10-23

> The two predeclared arms of `latent-sensitivity` (`extended_covariates`: the six analysis covariates plus the ribosomal and mitochondrial pseudobulk fractions with the covariate-blind basis; `latent_adjusted`: the same eight covariates with MOTCO 0.9.0's latent basis conditioned on `stage`) run on the unchanged microglia primary export beside the microglia primary analysis, and the positive control restated on the `latent_adjusted` basis (`positive-control --basis latent_adjusted`: E2 and E3 read from each cell type's `latent_adjusted` arm and the diagnostic made on it, E1 unchanged), OpenSpec change `mtg-microglia-latent-sensitivity` (specs `motco-latent-adjustment-sensitivity`, `cell-type-positive-control`). The primary-basis control did not hold (`reports/motco-microglia-positive-control-2026-10-22.md`: E3 failed in `cohort_30`, the first latent component following the mitochondrial fraction more closely than stage); this change asks whether that failure is a failure of the latent construction or a consequence of the measured quality gradient, with the reading of both outcomes written in the proposal before any microglia arm ran. Every number below is pasted by script from the generated files (`latent_sensitivity.md`, `latent_sensitivity.csv`, `mediator_check.csv`, the comparison and analysis manifests, `result.json`, `axis_correlations.csv` / `geometry_shares.csv` of `axis-diagnostics` on the primary and on each arm, and `positive_control.md` of both bases) and was checked against the six microglia `result.json` records and both positive-control tables before commit. The hand-written text carries only the reading, under the rules below.

**Reading rules (binding).** The microglia primary analysis remains the primary result, and MOTCO's power comparison (`candidate_fails`) excludes the latent adjustment as its replacement; both arms are sensitivity to it. `cohort_30` is primary and `cohort_50` sensitivity. A rejection status that differs between `primary` and `extended_covariates` is covariate-set dependence, one that differs between `extended_covariates` and `latent_adjusted` is basis dependence; neither is a result of the arm. A CV-selected rank that differs is reported beside `angle`, not read as better or worse. `angle` and `shape` non-rejections stay uninformative and a `delta` rejection is a difference, not specifically of size, in every arm. An arm's first-axis technical correlation is a measurement-space finding only. The restated control is a diagnosis of the primary-basis control, which stays the control's outcome (section 4); nothing here changes any astrocyte output or reading.

## 1. Provenance

- Software: MOTCO 0.9.0 (git dependency pinned to `v0.9.0`, commit `0150adc27b48528f1999a31a146a3df17ae3047f` in `uv.lock`), scikit-learn 1.8.0, numpy 2.3.5, scipy 1.15.3, pandas 2.3.3, Python 3.12.15 (the Nix shell); sead-sexdiff revision `d2a14e04890c567e4bf54cf62b0f1b8aa551913f` (dirty: True; the implementation of this change was uncommitted at run time, as in the earlier sensitivity changes; the reproducibility rerun below reproduced the outputs byte for byte).
- Machine: `nixos-home` (NixOS workstation), 2026-10-23, `n_jobs = 1`. The microglia chain (`results/microglia/`: the pseudobulk with the technical table, the export, the MOTCO 0.9.0 primary analysis, its diagnostic, both stage-DE steps and the primary-basis control) was made on this machine on 2026-10-22 for the positive-control change, copied to the cluster and verified; nothing of it was rerun. The refactor of this change (`--cell-type` on `latent-sensitivity`, `--basis` on `positive-control`) reproduced the astrocyte `latent-sensitivity` outputs and the primary-basis control on this machine byte for byte before the microglia runs (every CSV and arm Markdown; `latent_sensitivity.md` apart from its runtime column; manifests apart from runtime, timestamps, git, the scratch output path and the added `cell_type` / `basis` keys and cell-type findings), and every pre-existing file under `results/` (718) hashed unchanged after the microglia run.
- cohort_30: primary export `results/microglia/motco_input/cohort_30` (manifest sha256 `fb6afd9c24d8e6eee88bad4d8cb5a3ead2ddbf282b5ecff4b196bb893b8a9936`, 76 donors); reference `results/microglia/motco_analysis/cohort_30` (`result.json` sha256 `c72967f920b2c766b31cd1a6957761f5e1d2c94c64f3cb25bda0c1e60ec2ceb4`, CV-selected rank 1, cv_mean_auroc 0.6921, MOTCO 0.9.0, seed 20261019, 999 permutations); read and checked, never rerun or modified; technical table `results/microglia/pseudobulk/donor_technical.csv` sha256 `ac54eef6a253135db7260eec0365b837c103495f1d3e1041c9b9e06ff288b5cd`, pseudobulk manifest `49b937223a485a5233fcc7953fc51f9cd343ba65ccbc5c2dc2c3800ce8d8cb1c`. Every manifest read records `microglia`.
- cohort_50: primary export `results/microglia/motco_input/cohort_50` (manifest sha256 `8f0fd84ab291f93c60ea40f8f864ce8aad830005331918a8f2c65d97f2389bad`, 63 donors); reference `results/microglia/motco_analysis/cohort_50` (`result.json` sha256 `bf56d2c959dbffbf31fb8bd44698b191b52ceb2b3c2eeb7065df0e24eb56d77e`, CV-selected rank 1, cv_mean_auroc 0.6414, MOTCO 0.9.0, seed 20261019, 999 permutations); read and checked, never rerun or modified; technical table `results/microglia/pseudobulk/donor_technical.csv` sha256 `ac54eef6a253135db7260eec0365b837c103495f1d3e1041c9b9e06ff288b5cd`, pseudobulk manifest `49b937223a485a5233fcc7953fc51f9cd343ba65ccbc5c2dc2c3800ce8d8cb1c`. Every manifest read records `microglia`.
- Settings shared by every row: seed 20261019, 999 RRPP permutations (the same permutation indices in every row), α = 0.05 (p ≤ α), `n_jobs = 1`; blocks `rna` and `atac` in identity units, stage order 0, 1, 2, net-displacement orientation, CV rank selection at MOTCO's defaults, missing-covariate policy `error` (nothing imputed). The registry of arms, the covariate set (`LATENT_SENSITIVITY_COVARIATES`: `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk`), the `stage` conditioning and the mediator check are the astrocyte change's, unchanged for microglia.
- Astrocyte arms (the reference cell type of the restated control): `results/motco_latent_sensitivity/` (both arms, `comparison/`, `diagnostics/`; 72 files, made on `nixos-office` on 2026-10-21) was pulled from the cluster's copy with rsync on 2026-10-23 and hashed identically to the cluster's; all 13 sha256 values the astrocyte report cites are found on disk, both comparison manifests' reference, arm and export hashes equal the files now on disk, and each arm diagnostic records its arm's `result.json`. The astrocyte arm manifests predate the cell-type registry and read as astrocyte with a finding; no astrocyte file was regenerated.
- Diagnostics: `sead-sexdiff axis-diagnostics --cell-type microglia --analysis-dir <arm root> --out-dir results/microglia/motco_latent_sensitivity/diagnostics/<arm>` on each arm, with the microglia primary export as input; the input-hash and cell-type checks passed; no rank sweep exists for microglia (recorded by the diagnostic as a finding). The diagnostic labels the analysis it reads `primary`; the tables below name the analysis. Its partial correlations condition on the six analysis covariates (its definition), also for the arms.
- Reproducibility: `latent-sensitivity --cell-type microglia --cohort cohort_30 --force` rerun into a scratch directory with the same inputs, MOTCO release, seed and `n_jobs` gave byte-identical `latent_sensitivity.csv`, `mediator_check.csv` and both arms' `latent_scores.csv`, `pls_loadings.csv`, `stage_means.csv` and `analysis_metadata.csv` (all 10 files compared). `results/microglia/motco_latent_sensitivity/` and `results/microglia/positive_control/latent_adjusted/` were copied to the cluster (`/home1/tgonza/sead_sexdiff/results/microglia/`) and all 66 files hashed identically there.
- Positive control: `positive-control --basis latent_adjusted --cohort all` on 2026-10-23 into `results/microglia/positive_control/latent_adjusted/<cohort>/`; the primary-basis outputs (`results/microglia/positive_control/<cohort>/`, made 2026-10-22) hashed unchanged. Before anything was written the chain was checked per cell type: the arm's manifest records the cell type and the sha256 of the five primary export files now on disk, the arm is CV-selected with `stage` conditioning, the arm's comparison manifest names the primary `result.json` now on disk (the analysis E1's stage-DE manifests are checked against) and this arm directory, and the arm diagnostic records the arm's `result.json`. E1 is read unchanged from both stage-DE summaries (`e1_basis_independent`). The step was also checked to stop with nothing written on a copy of the microglia root with a comparison manifest naming another primary, with the arm removed and with the arm diagnostic removed.

## 2. Primary cohort: cohort_30 (≥ 30 microglia nuclei in both omics, n = 76)

### Overview (generated from `latent_sensitivity.csv`)

| analysis | covariates | conditioning | rank | cv_mean_auroc | delta (p) | angle (p) | shape (p) | rejected |
|---|---:|---|---:|---:|---:|---:|---:|---|
| primary | 6 | none | 1 | 0.6921 | 8.0642 (0.391) | 0 (1.0) | 0.154684 (0.637) | none |
| extended_covariates | 8 | none | 1 | 0.6921 | 3.00017 (0.587) | 0 (1.0) | 0.190112 (0.611) | none |
| latent_adjusted | 8 | stage | 1 | 0.6415 | 4.19879 (0.35) | 0 (1.0) | 0.170999 (0.274) | none |

### Mediator check (pasted from `mediator_check.csv`, written before the arms ran)

| covariate | subset | n | Spearman | flagged |
|---|---|---:|---:|:---:|
| `ribosomal_pseudobulk` | all | 76 | 0.053 | no |
| `ribosomal_pseudobulk` | F | 45 | -0.030 | no |
| `ribosomal_pseudobulk` | M | 31 | 0.127 | no |
| `mitochondrial_pseudobulk` | all | 76 | 0.498 | yes |
| `mitochondrial_pseudobulk` | F | 45 | 0.577 | yes |
| `mitochondrial_pseudobulk` | M | 31 | 0.329 | yes |

### First latent axis per analysis (generated from `axis_correlations.csv`)

Spearman correlation of `pls_0` with the technical-axis report's headline variables, stage and sex, raw and (in parentheses) partial given the six analysis covariates; n/a where a variable is one of them. No p-value is attached.

| variable | primary | extended_covariates | latent_adjusted |
|---|---:|---:|---:|
| `ribosomal_pseudobulk` | 0.01 (-0.11) | 0.01 (-0.11) | -0.23 (-0.20) |
| `ribosomal_mean` | -0.12 (-0.22) | -0.12 (-0.22) | -0.30 (-0.27) |
| `mitochondrial_pseudobulk` | 0.84 (0.80) | 0.84 (0.80) | 0.11 (0.09) |
| `fraction_mitochondrial_umis_mean` | 0.84 (0.80) | 0.84 (0.80) | 0.13 (0.11) |
| `atac_tss_enrichment_score_mean` | -0.72 (-0.70) | -0.72 (-0.70) | -0.24 (-0.22) |
| `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | -0.73 (-0.69) | -0.73 (-0.69) | -0.31 (-0.25) |
| `tmm_rna` | -0.34 (-0.38) | -0.34 (-0.38) | -0.42 (-0.34) |
| `tmm_atac` | 0.46 (0.56) | 0.46 (0.56) | -0.27 (-0.10) |
| `log10_lib_rna` | 0.11 (0.41) | 0.11 (0.41) | 0.05 (0.05) |
| `log10_lib_atac` | -0.35 (-0.46) | -0.35 (-0.46) | 0.31 (0.13) |
| `rin` | -0.30 (n/a) | -0.30 (n/a) | 0.03 (n/a) |
| `pmi` | 0.12 (n/a) | 0.12 (n/a) | 0.08 (n/a) |
| `stage` | 0.69 (0.65) | 0.69 (0.65) | -0.53 (-0.65) |
| `male` | -0.28 (-0.22) | -0.28 (-0.22) | 0.05 (-0.00) |

- `primary` `pls_0` (score SD 24.24; CV-selected rank 1): 28 of 47 variables flagged at |ρ| ≥ 0.3 (raw or partial): `mitochondrial_pseudobulk` 0.84, `fraction_mitochondrial_umis_mean` 0.84, `mitochondrial_mean` 0.84, `fraction_mitochondrial_umis_sd` 0.83, `mitochondrial_sd` 0.83, `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` -0.73, `gex_sequencing_saturation_mean` 0.72, `atac_tss_enrichment_score_mean` -0.72…. Geometry share along it: net_displacement_Female 1.000, net_displacement_Male 1.000, sex_difference_by_stage 1.000, sex_difference_centred 1.000.
- `extended_covariates` `pls_0` (score SD 24.24; CV-selected rank 1): 28 of 47 variables flagged at |ρ| ≥ 0.3 (raw or partial): `mitochondrial_pseudobulk` 0.84, `fraction_mitochondrial_umis_mean` 0.84, `mitochondrial_mean` 0.84, `fraction_mitochondrial_umis_sd` 0.83, `mitochondrial_sd` 0.83, `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` -0.73, `gex_sequencing_saturation_mean` 0.72, `atac_tss_enrichment_score_mean` -0.72…. Geometry share along it: net_displacement_Female 1.000, net_displacement_Male 1.000, sex_difference_by_stage 1.000, sex_difference_centred 1.000.
- `latent_adjusted` `pls_0` (score SD 15.46; CV-selected rank 1): 10 of 47 variables flagged at |ρ| ≥ 0.3 (raw or partial): `stage` -0.53, `tmm_rna` -0.42, `log10_n_atac` 0.38, `gex_reads_mapped_confidently_to_intronic_regions_sd` 0.34, `doublet_score_mean` 0.34, `n_libraries_rna` 0.32, `log10_lib_atac` 0.31, `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` -0.31…. Geometry share along it: net_displacement_Female 1.000, net_displacement_Male 1.000, sex_difference_by_stage 1.000, sex_difference_centred 1.000.

### `latent_adjusted` feature_r2 summary (from `result.json`)

- conditioning `stage`, 8 covariates, 76 donors, residual df 65, smallest CV training df 27; per-feature partial R² of the covariates: q05 0.0772, q25 0.1414, q50 0.2005, q75 0.2697, q95 0.4172, mean 0.2158, fraction above 0.5 0.0173.

### Comparison (pasted from `results/microglia/motco_latent_sensitivity/comparison/cohort_30/latent_sensitivity.md`)

### MOTCO latent-adjustment sensitivity: cohort_30

Covariate-set and latent-basis sensitivity of the covariate-adjusted trajectory analysis. The `primary` row is the primary analysis (CV-selected rank 1, the six analysis covariates, covariate-blind basis), read from its `result.json` and never rerun; it stays the primary result, and MOTCO's power comparison (`candidate_fails`) excludes the latent adjustment as its replacement. `extended_covariates` adds `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` (raw fractions from `donor_technical.csv`) with the covariate-blind basis; `latent_adjusted` enters the same eight covariates and fits the basis on features residualized on them conditional on `stage`, with the trajectory model, LS-means and RRPP of `extended_covariates`, so those two rows differ through the basis and the CV-selected rank only. Both arms are read as sensitivity, never as the result: a rejection status that differs between `primary` and `extended_covariates` is covariate-set dependence, one that differs between `extended_covariates` and `latent_adjusted` is basis dependence; a CV-selected rank that differs is reported beside `angle`, not read as better or worse; `angle` and `shape` non-rejections are uninformative and a `delta` rejection is a difference, not specifically of size, in every row. 999 RRPP permutations with seed 20261019 in every row (the same permutation indices); α = 0.05 (rejected when p ≤ α). MOTCO measured the latent adjustment's operating characteristics at its Phase 6 design only (Type I `adopt_stage`; power `candidate_fails`, with `shape` power lost under a nuisance axis and small-magnitude power lost without one), on neither these blocks nor this eight-covariate set, so neither arm carries a Type I or power claim of its own.

#### delta

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 1 | 8.0642 | 0.391 | no | 7.23504 | 5.27064 | 17.4616 |
| extended_covariates | 8 | none | 1 | 3.00017 | 0.587 | no | 4.36875 | 3.25874 | 10.7028 |
| latent_adjusted | 8 | stage | 1 | 4.19879 | 0.35 | no | 3.6585 | 2.82658 | 8.95446 |

#### angle (degrees)

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 1 | 0 | 1 | no | 0 | 0 | 0 |
| extended_covariates | 8 | none | 1 | 0 | 1 | no | 0 | 0 | 0 |
| latent_adjusted | 8 | stage | 1 | 0 | 1 | no | 0 | 0 | 0 |

#### shape

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 1 | 0.154684 | 0.637 | no | 0.25893 | 0.192607 | 0.638734 |
| extended_covariates | 8 | none | 1 | 0.190112 | 0.611 | no | 0.27856 | 0.198728 | 0.657107 |
| latent_adjusted | 8 | stage | 1 | 0.170999 | 0.274 | no | 0.124711 | 0.0935683 | 0.302866 |

#### Runs

| arm | selection | rank | cv_mean_auroc | covariates (in order) | conditioning | runtime (s) | directory |
|---|---|---:|---:|---|---|---:|---|
| primary | cv | 1 | 0.6921 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac` | none | 88.9 | `results/microglia/motco_analysis/cohort_30` |
| extended_covariates | cv | 1 | 0.6921 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` | none | 109.0 | `extended_covariates/motco_analysis/cohort_30` |
| latent_adjusted | cv | 1 | 0.6415 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` | stage | 118.1 | `latent_adjusted/motco_analysis/cohort_30` |

#### Latent adjustment

MOTCO's `latent_adjustment` block: the residualization's degrees of freedom on the full sample (`n − 1 − k − q`) and the smallest over the CV training sets, and the quantiles over features of the partial R² of the covariates given the intercept and the conditioning indicators (pooled-standardized features).

| arm | conditioning | residual df | smallest CV training df | R² q05 | R² q25 | R² q50 | R² q75 | R² q95 | R² mean | R² fraction > 0.5 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| latent_adjusted | stage | 65 | 27 | 0.0772 | 0.1414 | 0.2005 | 0.2697 | 0.4172 | 0.2158 | 0.0173 |

#### Mediator check

Spearman correlation of each added covariate with stage (ordinal 0, 1, 2) over all cohort donors and within each sex, written before the arms ran; |ρ| ≥ 0.3 (`DIAGNOSTIC_CORRELATION_THRESHOLD`) flags a mediator candidate (flagged rows first). No p-value is computed, and a flag removes no covariate and changes no arm: no simulation can tell a confounder from a mediator, so a flagged fraction may carry stage signal that the `latent_adjusted` basis removed (and that the covariate adjustment of both arms absorbed).

| covariate | subset | n | Spearman | flagged |
|---|---|---:|---:|:---:|
| `mitochondrial_pseudobulk` | all | 76 | 0.498 | yes |
| `mitochondrial_pseudobulk` | F | 45 | 0.577 | yes |
| `mitochondrial_pseudobulk` | M | 31 | 0.329 | yes |
| `ribosomal_pseudobulk` | all | 76 | 0.053 | no |
| `ribosomal_pseudobulk` | F | 45 | -0.030 | no |
| `ribosomal_pseudobulk` | M | 31 | 0.127 | no |

#### Findings

- mediator candidate: mitochondrial_pseudobulk has Spearman 0.498 with stage (all, n = 76); the latent_adjusted basis may have removed stage signal carried by it (no simulation can tell a confounder from a mediator; the flag changes no arm)
- mediator candidate: mitochondrial_pseudobulk has Spearman 0.577 with stage (F, n = 45); the latent_adjusted basis may have removed stage signal carried by it (no simulation can tell a confounder from a mediator; the flag changes no arm)
- mediator candidate: mitochondrial_pseudobulk has Spearman 0.329 with stage (M, n = 31); the latent_adjusted basis may have removed stage signal carried by it (no simulation can tell a confounder from a mediator; the flag changes no arm)

## 3. Sensitivity cohort: cohort_50 (≥ 50 microglia nuclei in both omics, n = 63)

### Overview (generated from `latent_sensitivity.csv`)

| analysis | covariates | conditioning | rank | cv_mean_auroc | delta (p) | angle (p) | shape (p) | rejected |
|---|---:|---|---:|---:|---:|---:|---:|---|
| primary | 6 | none | 1 | 0.6414 | 12.4641 (0.227) | 0 (1.0) | 0.269568 (0.399) | none |
| extended_covariates | 8 | none | 1 | 0.6414 | 2.46906 (0.699) | 0 (1.0) | 0.270584 (0.351) | none |
| latent_adjusted | 8 | stage | 3 | 0.6205 | 1.11008 (0.875) | 12.1741 (0.611) | 0.176571 (0.114) | none |

### Mediator check (pasted from `mediator_check.csv`, written before the arms ran)

| covariate | subset | n | Spearman | flagged |
|---|---|---:|---:|:---:|
| `ribosomal_pseudobulk` | all | 63 | -0.029 | no |
| `ribosomal_pseudobulk` | F | 36 | -0.226 | no |
| `ribosomal_pseudobulk` | M | 27 | 0.158 | no |
| `mitochondrial_pseudobulk` | all | 63 | 0.440 | yes |
| `mitochondrial_pseudobulk` | F | 36 | 0.541 | yes |
| `mitochondrial_pseudobulk` | M | 27 | 0.254 | no |

### First latent axis per analysis (generated from `axis_correlations.csv`)

Spearman correlation of `pls_0` with the technical-axis report's headline variables, stage and sex, raw and (in parentheses) partial given the six analysis covariates; n/a where a variable is one of them. No p-value is attached.

| variable | primary | extended_covariates | latent_adjusted |
|---|---:|---:|---:|
| `ribosomal_pseudobulk` | -0.09 (-0.09) | -0.09 (-0.09) | -0.18 (-0.11) |
| `ribosomal_mean` | -0.23 (-0.26) | -0.23 (-0.26) | -0.26 (-0.18) |
| `mitochondrial_pseudobulk` | 0.78 (0.82) | 0.78 (0.82) | 0.14 (0.11) |
| `fraction_mitochondrial_umis_mean` | 0.78 (0.81) | 0.78 (0.81) | 0.15 (0.13) |
| `atac_tss_enrichment_score_mean` | -0.62 (-0.65) | -0.62 (-0.65) | -0.29 (-0.28) |
| `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | -0.69 (-0.70) | -0.69 (-0.70) | -0.33 (-0.27) |
| `tmm_rna` | -0.41 (-0.48) | -0.41 (-0.48) | -0.31 (-0.21) |
| `tmm_atac` | 0.47 (0.67) | 0.47 (0.67) | -0.25 (-0.07) |
| `log10_lib_rna` | 0.13 (0.49) | 0.13 (0.49) | 0.13 (0.02) |
| `log10_lib_atac` | -0.32 (-0.50) | -0.32 (-0.50) | 0.34 (0.23) |
| `rin` | -0.19 (n/a) | -0.19 (n/a) | 0.02 (n/a) |
| `pmi` | 0.04 (n/a) | 0.04 (n/a) | 0.08 (n/a) |
| `stage` | 0.73 (0.67) | 0.73 (0.67) | -0.53 (-0.67) |
| `male` | -0.27 (-0.28) | -0.27 (-0.28) | 0.04 (0.01) |

- `primary` `pls_0` (score SD 23.64; CV-selected rank 1): 25 of 47 variables flagged at |ρ| ≥ 0.3 (raw or partial): `mitochondrial_pseudobulk` 0.78, `mitochondrial_mean` 0.78, `fraction_mitochondrial_umis_mean` 0.78, `fraction_mitochondrial_umis_sd` 0.77, `mitochondrial_sd` 0.77, `stage` 0.73, `gex_sequencing_saturation_mean` 0.72, `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` -0.69…. Geometry share along it: net_displacement_Female 1.000, net_displacement_Male 1.000, sex_difference_by_stage 1.000, sex_difference_centred 1.000.
- `extended_covariates` `pls_0` (score SD 23.64; CV-selected rank 1): 25 of 47 variables flagged at |ρ| ≥ 0.3 (raw or partial): `mitochondrial_pseudobulk` 0.78, `mitochondrial_mean` 0.78, `fraction_mitochondrial_umis_mean` 0.78, `fraction_mitochondrial_umis_sd` 0.77, `mitochondrial_sd` 0.77, `stage` 0.73, `gex_sequencing_saturation_mean` 0.72, `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` -0.69…. Geometry share along it: net_displacement_Female 1.000, net_displacement_Male 1.000, sex_difference_by_stage 1.000, sex_difference_centred 1.000.
- `latent_adjusted` `pls_0` (score SD 15.72; CV-selected rank 3): 8 of 47 variables flagged at |ρ| ≥ 0.3 (raw or partial): `stage` -0.53, `gex_reads_mapped_confidently_to_intronic_regions_sd` 0.41, `n_libraries_rna` 0.39, `log10_n_atac` 0.35, `log10_lib_atac` 0.34, `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` -0.33, `tmm_rna` -0.31, `atac_percent_duplicates_mean` -0.24. Geometry share along it: net_displacement_Female 0.737, net_displacement_Male 0.654, sex_difference_by_stage 0.122, sex_difference_centred 0.064.

### `latent_adjusted` feature_r2 summary (from `result.json`)

- conditioning `stage`, 8 covariates, 63 donors, residual df 52, smallest CV training df 20; per-feature partial R² of the covariates: q05 0.0813, q25 0.1480, q50 0.2078, q75 0.2794, q95 0.4109, mean 0.2217, fraction above 0.5 0.0141.

### Comparison (pasted from `results/microglia/motco_latent_sensitivity/comparison/cohort_50/latent_sensitivity.md`)

### MOTCO latent-adjustment sensitivity: cohort_50

Covariate-set and latent-basis sensitivity of the covariate-adjusted trajectory analysis. The `primary` row is the primary analysis (CV-selected rank 1, the six analysis covariates, covariate-blind basis), read from its `result.json` and never rerun; it stays the primary result, and MOTCO's power comparison (`candidate_fails`) excludes the latent adjustment as its replacement. `extended_covariates` adds `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` (raw fractions from `donor_technical.csv`) with the covariate-blind basis; `latent_adjusted` enters the same eight covariates and fits the basis on features residualized on them conditional on `stage`, with the trajectory model, LS-means and RRPP of `extended_covariates`, so those two rows differ through the basis and the CV-selected rank only. Both arms are read as sensitivity, never as the result: a rejection status that differs between `primary` and `extended_covariates` is covariate-set dependence, one that differs between `extended_covariates` and `latent_adjusted` is basis dependence; a CV-selected rank that differs is reported beside `angle`, not read as better or worse; `angle` and `shape` non-rejections are uninformative and a `delta` rejection is a difference, not specifically of size, in every row. 999 RRPP permutations with seed 20261019 in every row (the same permutation indices); α = 0.05 (rejected when p ≤ α). MOTCO measured the latent adjustment's operating characteristics at its Phase 6 design only (Type I `adopt_stage`; power `candidate_fails`, with `shape` power lost under a nuisance axis and small-magnitude power lost without one), on neither these blocks nor this eight-covariate set, so neither arm carries a Type I or power claim of its own.

#### delta

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 1 | 12.4641 | 0.227 | no | 8.45222 | 5.99525 | 19.5118 |
| extended_covariates | 8 | none | 1 | 2.46906 | 0.699 | no | 5.02046 | 3.79773 | 12.4959 |
| latent_adjusted | 8 | stage | 3 | 1.11008 | 0.875 | no | 5.86056 | 4.59823 | 14.6532 |

#### angle (degrees)

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 1 | 0 | 1 | no | 0 | 0 | 0 |
| extended_covariates | 8 | none | 1 | 0 | 1 | no | 0 | 0 | 0 |
| latent_adjusted | 8 | stage | 3 | 12.1741 | 0.611 | no | 16.8349 | 10.5021 | 36.7317 |

#### shape

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 1 | 0.269568 | 0.399 | no | 0.255498 | 0.194219 | 0.624989 |
| extended_covariates | 8 | none | 1 | 0.270584 | 0.351 | no | 0.230102 | 0.170008 | 0.5502 |
| latent_adjusted | 8 | stage | 3 | 0.176571 | 0.114 | no | 0.108847 | 0.0564956 | 0.216514 |

#### Runs

| arm | selection | rank | cv_mean_auroc | covariates (in order) | conditioning | runtime (s) | directory |
|---|---|---:|---:|---|---|---:|---|
| primary | cv | 1 | 0.6414 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac` | none | 84.1 | `results/microglia/motco_analysis/cohort_50` |
| extended_covariates | cv | 1 | 0.6414 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` | none | 102.7 | `extended_covariates/motco_analysis/cohort_50` |
| latent_adjusted | cv | 3 | 0.6205 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` | stage | 112.3 | `latent_adjusted/motco_analysis/cohort_50` |

#### Latent adjustment

MOTCO's `latent_adjustment` block: the residualization's degrees of freedom on the full sample (`n − 1 − k − q`) and the smallest over the CV training sets, and the quantiles over features of the partial R² of the covariates given the intercept and the conditioning indicators (pooled-standardized features).

| arm | conditioning | residual df | smallest CV training df | R² q05 | R² q25 | R² q50 | R² q75 | R² q95 | R² mean | R² fraction > 0.5 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| latent_adjusted | stage | 52 | 20 | 0.0813 | 0.1480 | 0.2078 | 0.2794 | 0.4109 | 0.2217 | 0.0141 |

#### Mediator check

Spearman correlation of each added covariate with stage (ordinal 0, 1, 2) over all cohort donors and within each sex, written before the arms ran; |ρ| ≥ 0.3 (`DIAGNOSTIC_CORRELATION_THRESHOLD`) flags a mediator candidate (flagged rows first). No p-value is computed, and a flag removes no covariate and changes no arm: no simulation can tell a confounder from a mediator, so a flagged fraction may carry stage signal that the `latent_adjusted` basis removed (and that the covariate adjustment of both arms absorbed).

| covariate | subset | n | Spearman | flagged |
|---|---|---:|---:|:---:|
| `mitochondrial_pseudobulk` | all | 63 | 0.440 | yes |
| `mitochondrial_pseudobulk` | F | 36 | 0.541 | yes |
| `ribosomal_pseudobulk` | all | 63 | -0.029 | no |
| `ribosomal_pseudobulk` | F | 36 | -0.226 | no |
| `ribosomal_pseudobulk` | M | 27 | 0.158 | no |
| `mitochondrial_pseudobulk` | M | 27 | 0.254 | no |

#### Findings

- CV-selected rank differs: 1 in extended_covariates, 3 in latent_adjusted; reported beside angle, not read as better or worse
- mediator candidate: mitochondrial_pseudobulk has Spearman 0.440 with stage (all, n = 63); the latent_adjusted basis may have removed stage signal carried by it (no simulation can tell a confounder from a mediator; the flag changes no arm)
- mediator candidate: mitochondrial_pseudobulk has Spearman 0.541 with stage (F, n = 36); the latent_adjusted basis may have removed stage signal carried by it (no simulation can tell a confounder from a mediator; the flag changes no arm)

## 4. The positive control restated on the `latent_adjusted` basis

`positive-control --basis latent_adjusted` reads E2 and E3 from each cell type's `latent_adjusted` arm and the diagnostic made on it; E1 is basis-independent (the stage-DE summaries are read unchanged). The expectations, their comparisons, E3's technical variables (the two adjusted fractions stay in the compared set) and the all-three rule are the primary basis's; the chain (arm made from the primary export on disk, compared with the primary analysis on disk, the registered conditioning, the diagnostic made from the arm, every manifest's cell type) was checked before anything was written. The primary-basis tables are repeated from `reports/motco-microglia-positive-control-2026-10-22.md` for the comparison.

### cohort_30, basis `primary` (primary cohort; pasted from `results/microglia/positive_control/cohort_30/positive_control.md`)

| Expectation | Method | Quantity | Control | Value | Compared with | Value | Holds | Outcome | Detail |
|---|---|---|---|---|---|---|---|---|---|
| E1 | deseq2 | exported RNA features at FDR 0.05: factor_additive union | microglia | 61 | astrocyte | 1 | yes | holds | tested 5000 (control), 5000 (reference) |
| E1 | limma | exported RNA features at FDR 0.05: factor_additive stage_joint | microglia | 2 | astrocyte | 1 | yes | holds | tested 5000 (control), 5000 (reference) |
| E1 | both | E1 under both stage-DE methods |  |  |  |  | yes | holds | holds under both methods |
| E2 |  | CV mean AUROC of the primary analysis (latent.cv_mean_auroc) | microglia | 0.692109 | astrocyte | 0.631278 | yes | holds |  |
| E3 |  | \|Spearman\| of pls_0 (primary) with stage vs the largest with a technical variable | stage | 0.688175 | mitochondrial_pseudobulk | 0.843828 | no | fails | 39 technical variables compared (donor_technical.csv, manifest.json (export, tmm)) |
| overall |  | E1, E2 and E3 all hold |  |  |  |  | no | does_not_hold | the primary cohort: this row is the control's outcome; not holding: E3 |

Overall: **does not hold** (E3 not holding); E3 technical maximum `mitochondrial_pseudobulk` (donor_technical.csv) |ρ| 0.8438.

### cohort_30, basis `latent_adjusted` (primary cohort; pasted from `results/microglia/positive_control/latent_adjusted/cohort_30/positive_control.md`)

| Expectation | Method | Quantity | Control | Value | Compared with | Value | Holds | Outcome | Detail |
|---|---|---|---|---|---|---|---|---|---|
| E1 | deseq2 | exported RNA features at FDR 0.05: factor_additive union | microglia | 61 | astrocyte | 1 | yes | holds | tested 5000 (control), 5000 (reference) |
| E1 | limma | exported RNA features at FDR 0.05: factor_additive stage_joint | microglia | 2 | astrocyte | 1 | yes | holds | tested 5000 (control), 5000 (reference) |
| E1 | both | E1 under both stage-DE methods |  |  |  |  | yes | holds | holds under both methods |
| E2 |  | CV mean AUROC of the latent_adjusted arm analysis (basis latent_adjusted) (latent.cv_mean_auroc) | microglia | 0.641503 | astrocyte | 0.629944 | yes | holds |  |
| E3 |  | \|Spearman\| of pls_0 (latent_adjusted arm, labelled primary by its diagnostic) with stage vs the largest with a technical variable | stage | 0.527877 | tmm_rna | 0.415256 | yes | holds | 39 technical variables compared (donor_technical.csv, manifest.json (export, tmm)) |
| overall |  | E1, E2 and E3 all hold |  |  |  |  | yes | holds | the primary cohort: this row is the control's outcome |

Overall: **holds**; E3 technical maximum `tmm_rna` (manifest.json (export, tmm)) |ρ| 0.4153.

### cohort_50, basis `primary` (descriptive; pasted from `results/microglia/positive_control/cohort_50/positive_control.md`)

| Expectation | Method | Quantity | Control | Value | Compared with | Value | Holds | Outcome | Detail |
|---|---|---|---|---|---|---|---|---|---|
| E1 | deseq2 | exported RNA features at FDR 0.05: factor_additive union | microglia | 35 | astrocyte | 2 | yes | holds | tested 5000 (control), 5000 (reference) |
| E1 | limma | exported RNA features at FDR 0.05: factor_additive stage_joint | microglia | 0 | astrocyte | 1 | no | fails | tested 5000 (control), 5000 (reference) |
| E1 | both | E1 under both stage-DE methods |  |  |  |  | no | method_dependent | held under deseq2 only; counts as not holding |
| E2 |  | CV mean AUROC of the primary analysis (latent.cv_mean_auroc) | microglia | 0.641431 | astrocyte | 0.661645 | no | fails |  |
| E3 |  | \|Spearman\| of pls_0 (primary) with stage vs the largest with a technical variable | stage | 0.730773 | mitochondrial_pseudobulk | 0.784562 | no | fails | 39 technical variables compared (donor_technical.csv, manifest.json (export, tmm)) |
| overall |  | E1, E2 and E3 all hold |  |  |  |  | no | does_not_hold | not the primary cohort (cohort_30); descriptive; not holding: E1, E2, E3 |

Overall: **does not hold** (E1, E2, E3 not holding); E3 technical maximum `mitochondrial_pseudobulk` (donor_technical.csv) |ρ| 0.7846.

### cohort_50, basis `latent_adjusted` (descriptive; pasted from `results/microglia/positive_control/latent_adjusted/cohort_50/positive_control.md`)

| Expectation | Method | Quantity | Control | Value | Compared with | Value | Holds | Outcome | Detail |
|---|---|---|---|---|---|---|---|---|---|
| E1 | deseq2 | exported RNA features at FDR 0.05: factor_additive union | microglia | 35 | astrocyte | 2 | yes | holds | tested 5000 (control), 5000 (reference) |
| E1 | limma | exported RNA features at FDR 0.05: factor_additive stage_joint | microglia | 0 | astrocyte | 1 | no | fails | tested 5000 (control), 5000 (reference) |
| E1 | both | E1 under both stage-DE methods |  |  |  |  | no | method_dependent | held under deseq2 only; counts as not holding |
| E2 |  | CV mean AUROC of the latent_adjusted arm analysis (basis latent_adjusted) (latent.cv_mean_auroc) | microglia | 0.620498 | astrocyte | 0.664479 | no | fails |  |
| E3 |  | \|Spearman\| of pls_0 (latent_adjusted arm, labelled primary by its diagnostic) with stage vs the largest with a technical variable | stage | 0.532681 | gex_reads_mapped_confidently_to_intronic_regions_sd | 0.406702 | yes | holds | 39 technical variables compared (donor_technical.csv, manifest.json (export, tmm)) |
| overall |  | E1, E2 and E3 all hold |  |  |  |  | no | does_not_hold | not the primary cohort (cohort_30); descriptive; not holding: E1, E2 |

Overall: **does not hold** (E1, E2 not holding); E3 technical maximum `gex_reads_mapped_confidently_to_intronic_regions_sd` (donor_technical.csv) |ρ| 0.4067.

### `extended_covariates` equals the primary (the covariate-blind basis)

By MOTCO's contract the covariate-blind basis, its CV folds, selected rank and CV AUROC do not depend on the covariate set, so `extended_covariates` repeats the primary's `pls_0` and `cv_mean_auroc` and gets no restated control. Checked from the records:

cohort_30:

| cell type | analysis | rank | cv_mean_auroc | equal to the primary |
|---|---|---:|---:|:---:|
| microglia | primary | 1 | 0.6921 | |
| microglia | extended_covariates | 1 | 0.6921 | yes |
| astrocyte | primary | 1 | 0.6313 | |
| astrocyte | extended_covariates | 1 | 0.6313 | yes |

cohort_50:

| cell type | analysis | rank | cv_mean_auroc | equal to the primary |
|---|---|---:|---:|:---:|
| microglia | primary | 1 | 0.6414 | |
| microglia | extended_covariates | 1 | 0.6414 | yes |
| astrocyte | primary | 1 | 0.6616 | |
| astrocyte | extended_covariates | 1 | 0.6616 | yes |

Equal in every cell type and cohort.

## 5. Reading

- **The restated control holds in the primary cohort.** On the `latent_adjusted` basis in `cohort_30`, E1 holds under both methods (61 against 1 exported RNA stage genes under DESeq2, 2 against 1 under limma; unchanged from the primary basis, as E1 is basis-independent), E2 holds (microglia CV AUROC 0.6415 against the astrocytes' 0.6299) and E3 holds (|ρ| of the microglia `pls_0` with stage 0.5279 against 0.4153 for the largest technical variable, `tmm_rna`). By the reading written in the proposal before any microglia arm ran, the measured quality gradient explains the primary control's failure: once the library-quality direction is conditioned out of the basis, the pipeline resolves stage in the cell type with the stronger stage signal. The first microglia axis moves off the mitochondrial fraction (0.84 on the primary basis, 0.11 on the adjusted one; 28 of 47 diagnostic variables flagged against 10) and the largest technical correlation falls below the stage correlation, in both cohorts. This does not promote the arm (MOTCO's power comparison, `candidate_fails`, still excludes the latent adjustment as a replacement of the primary), does not change the primary-basis outcome (the control does not hold on the primary basis, `reports/motco-microglia-positive-control-2026-10-22.md`) and does not change the astrocyte reading. It is the diagnosis a later MOTCO-side validation of the latent adjustment on blocks with a dominant nuisance axis, or a later revised-primary proposal, would cite.
- **What the hold rests on, and the mediator caveat.** E3 holds with a smaller stage correlation than on the primary basis (0.5279 against 0.6882): the adjustment removed the quality direction and, with it, stage signal the mitochondrial fraction carries (the mediator check flags it in every subset of `cohort_30`: all 0.498, F 0.577, M 0.329; `ribosomal_pseudobulk` is not flagged in microglia). The caveat cuts both ways, as predeclared: the stage |ρ| fell rather than rose, so the quality gradient and the stage signal are entangled in these donors, and the adjusted basis reads the part of stage that is not carried by the fraction. E2 holds by 0.0116 (the microglia CV AUROC fell from 0.6921 to 0.6415 under the adjustment, the astrocytes' from 0.6313 to 0.6299); the expectation is strict `>` and ties do not hold, and nothing is tested, so the margin is reported, not interpreted.
- **`cohort_50` is descriptive.** E3 holds there too on the adjusted basis (0.5327 against 0.4067, `gex_reads_mapped_confidently_to_intronic_regions_sd`), E2 fails (0.6205 against 0.6645) and E1 stays method-dependent (35 against 2 under DESeq2, 0 against 1 under limma), as on the primary basis; the overall row does not hold, and the control's outcome is `cohort_30`'s.
- **`extended_covariates` keeps the primary's basis.** In both cell types and both cohorts the arm's CV-selected rank and CV AUROC equal the primary's (section 4), as MOTCO's contract implies; the arm changes only the trajectory model. With the two fractions entered, `delta` is 3.00017 against 8.0642 (`cohort_30`) and 2.46906 against 12.4641 (`cohort_50`), `shape` 0.190112 against 0.154684 and 0.270584 against 0.269568; none rejected in either.
- **No rejection in any row.** No statistic is rejected at α = 0.05 in the microglia primary analysis or in either arm, in either cohort, so no rejection status differs between `primary` and `extended_covariates` (no covariate-set dependence of a rejection) or between `extended_covariates` and `latent_adjusted` (no basis dependence of a rejection). The non-rejections are uninformative: `angle` and `shape` non-rejections carry no information under MOTCO's rules, and `delta` non-rejections are read with MOTCO's Phase 6-only operating characteristics. The microglia sex × stage result stays the rank-1, nothing-rejected result of the positive-control report.
- **`cohort_50` selects rank 3 under the latent adjustment.** The CV-selected rank differs from `extended_covariates` (1), so `angle` is no longer 0° by construction: 12.1741° (p 0.611), reported beside rank 3, not read as better or worse, and its non-rejection is uninformative. `cohort_30` keeps rank 1 in every row.
- **The adjusted first axis in microglia.** `pls_0` of `latent_adjusted` correlates at most 0.42 (`tmm_rna`, `cohort_30`) and 0.41 (`gex_reads_mapped_confidently_to_intronic_regions_sd`, `cohort_50`) with a technical variable, against 0.84 and 0.78 (`mitochondrial_pseudobulk`) for the primary, and -0.53 / -0.53 with stage (partial -0.65 / -0.67), against 0.69 / 0.73. The sign of a PLS component is arbitrary; the expectation reads |ρ|. This is a measurement-space finding about the arm, not a reinterpretation of any arm's statistics and not grounds to adopt it.
- **The two cell types under the same adjustment.** In astrocytes the adjustment moved the first axis onto stage (|ρ| 0.64–0.75 against 0.49–0.50, `reports/motco-astrocyte-latent-adjustment-2026-10-21.md`); in microglia it lowered the stage correlation (0.53 against 0.69–0.73) while removing the mitochondrial one. The difference is consistent with the mediator check: the mitochondrial fraction tracks stage more closely in microglia (0.498) than in astrocytes (0.312), so more of the microglia stage signal travels with the quality direction. Compared, not tested.

## 6. Limitations

- **Mediators.** The mitochondrial fraction is flagged in every subset of `cohort_30` (all 0.498, F 0.577, M 0.329) and in `cohort_50` (all 0.440, F 0.541); a pathology-driven change in microglial RNA quality would make it a mediator rather than a confounder, and the data cannot tell the two apart (MOTCO `docs/reports/latent-adjustment-type-i-2026-10-20.md` §4). The adjusted basis therefore reads the stage signal that is not carried by the fraction; the stage |ρ| along the adjusted first axis (0.53) is lower than along the primary one (0.69 / 0.73), and E3 holds on the adjusted basis with that smaller correlation. A failure for that reason would have counted; a hold for that reason is a hold with the caveat stated.
- **MOTCO's operating characteristics were measured at its Phase 6 design only.** Type I: `adopt_stage`; with `stage` conditioning `delta` 0.031, `angle` 0.067 and `shape` 0.058 against the bound 0.0706 (R = 448), `group_stage` exceeded on `angle` (Type I report, summary). Power: `candidate_fails` for `stage` (power report, summary), with `shape` power lost under a nuisance axis and small-magnitude power lost without one. The latent adjustment is a sensitivity analysis in MOTCO's binding conditions and in this report.
- **No operating characteristic on these blocks or at this covariate set.** MOTCO's checks used semisynthetic blocks of 498 features with two injected and six pure-noise covariates (eight declared); here 10,000 exported features of a second cell type and eight real covariates, one of them stage-associated. Neither arm carries a Type I or power claim of its own, and the restated control compares values, it tests nothing.
- **Degrees of freedom (cohort_30).** The residualization spends 11 columns (intercept, stage indicators, eight covariates) for 76 donors: residual df 65, smallest CV training df 27; the covariates' partial R² has median 0.2005 and exceeds 0.5 for a fraction 0.0173 of the features.
- **Degrees of freedom (cohort_50).** The same 11 columns for 63 donors: residual df 52, smallest CV training df 20; the covariates' partial R² has median 0.2078 and exceeds 0.5 for a fraction 0.0141 of the features. The smallest sex × stage cells (9 and 10 donors) make eight covariates a heavy adjustment at this n; no design was singular.
- **E2 compares CV AUROCs of adjusted bases.** On the adjusted basis each cell type's AUROC is the cross-validated separation of its own stage-conditioned basis; the comparison between cell types is as on the primary basis (the same covariate set, conditioning and CV settings in both), and its margin in `cohort_30` (0.0116) is reported without a test, as the expectation prescribes.
- **Library quality enters twice.** The nuclei counts (`log10_n_rna`, `log10_n_atac`) and the two fractions both proxy library quality; the adjustment is joint, so the collinear proxies do not change what is adjusted for, but the arms adjust for library quality through more than one proxy. `log10_n_atac` (0.38) and `log10_lib_atac` (0.31) remain flagged along the adjusted first axis in `cohort_30`; the former is an analysis covariate and not a technical variable of E3, the latter is and lies below the stage correlation.
- **The microglia population and depth.** The `Micro-PVM_*` bundle includes perivascular macrophages, and the microglia ATAC nuclei per donor are few (the positive-control report's cohort tables); the candidates named in the proposal for a failure remain limitations of a hold, since the adjusted basis is built on the same blocks.
- **Rank above `n_stages − 1`.** `cohort_50`'s `latent_adjusted` rank is 3; MOTCO's latent-rank ladder found a wider `angle` null above `n_stages − 1` (unadjusted, at a larger n), cited as what it is.
- **Uncommitted implementation at run time** (`dirty: True`); the `cohort_30` rerun reproduced the outputs byte for byte, and a post-commit rerun would carry a clean `git` field.
- **Everything the primary and positive-control reports list** (additive linear adjustment, structural analogue for the operating characteristics, fixed preprocessing choices, Monte Carlo error and `n_jobs`, the stage × quality entanglement at the donor level in SEA-AD) applies unchanged.

## 7. Next steps

- None follows for the primary analyses or the primary-basis control: both stay as they are. The hold is the diagnosis that a later MOTCO-side validation of the latent adjustment on blocks with a dominant nuisance axis, a later revised-primary proposal, or the predeclared donor-exclusion sensitivity (SEA-AD's severely-affected-donor flag, a donor-selection route to the same end) would cite; each needs its own predeclared change. `group_stage` conditioning, other covariate sets, and `rank-sensitivity` or `preprocessing-sensitivity` on microglia are out of scope.

## 8. Reproduction

```
scripts/sync_pseudobulk.sh --cell-type microglia            # results/microglia/pseudobulk/ from the cluster
rsync -a <cluster>:/home1/tgonza/sead_sexdiff/results/microglia/{motco_input,motco_analysis,motco_diagnostics,stage_de,stage_de_limma} results/microglia/
rsync -a <cluster>:/home1/tgonza/sead_sexdiff/results/motco_latent_sensitivity results/   # the astrocyte arms
uv run sead-sexdiff latent-sensitivity --cell-type microglia --cohort all
for a in extended_covariates latent_adjusted; do
  uv run sead-sexdiff axis-diagnostics --cell-type microglia \
      --analysis-dir results/microglia/motco_latent_sensitivity/$a/motco_analysis \
      --out-dir results/microglia/motco_latent_sensitivity/diagnostics/$a
done
uv run sead-sexdiff positive-control --basis latent_adjusted --cohort all
```

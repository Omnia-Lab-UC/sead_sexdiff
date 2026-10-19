# MOTCO astrocyte rank sensitivity: sex × ADNC stage — 2026-10-19

> Dimensionality sensitivity analysis of the covariate-adjusted comparison of female and male astrocyte trajectories
> across ADNC stages in the SEA-AD MTG (OpenSpec change `mtg-astrocyte-rank-sensitivity`, spec
> `motco-rank-sensitivity`). The primary analysis is the CV-selected one in `reports/motco-astrocyte-2026-10-19.md`;
> this report reruns it with the latent rank forced to 2, 3, 4, 6, 9, 12 and reads the result as
> sensitivity to the primary analysis, never as its replacement. Every statistic, p-value, null summary and rank
> below is pasted by script from the generated `rank_sensitivity.md` of each cohort
> (`results/motco_rank_sensitivity/<cohort>/`), which is rendered from the `result.json` of the reference analysis
> and of each forced-rank run, and was checked against those 14 records before commit. The stage means come from
> `stage_means.csv`.

**Run:** `sead-sexdiff rank-sensitivity --cohort all` on the NixOS workstation `nixos-home`, 2026-10-19, at commit
`9cd07a1` with the change's implementation uncommitted (the manifests record `dirty: true`),
`n_jobs = 1`: 6.0 s for the twelve forced runs (0.4–0.7 s each, no cross-validation) beside
reference analyses of 149 s and 152 s. The reference analyses and the export were produced on
the laptop workstation and pulled from the cluster (`/home1/tgonza/sead_sexdiff/results/`); their hashes matched
the sweep's inputs and settings. **Result:** nothing is rejected at α = 0.05 at any rank in either cohort. The
CV-selected analysis (rank 1 in both cohorts, nothing rejected) remains the primary result. The sweep gives no
evidence of a difference at any rank, and that is not evidence of agreement in size, direction or shape.

**In one paragraph.** For each cohort, MOTCO 0.8.0's `analyze_trajectories` was called six more times with
`forced_components` set to 2, 3, 4, 6, 9, 12, with the same RNA and ATAC blocks (5,000 + 5,000 log2-CPM
features), the same six covariates, the same net-displacement orientation, the same 999 RRPP permutations and the
same seed as the primary analysis, so that every rank permutes the same residual indices and the ranks differ only
in the latent space the trajectories are measured in. The first PLS component is the same at every rank, and the
female path along it stays the longer one. In the primary cohort (n = 80) the forced-rank `delta` p-values range
from 0.052 (rank 4) to 0.866; the closest to α is `delta` at rank 4
(p = 0.052), which is not a rejection and, had it crossed α, would have been rank dependence of the sensitivity
analysis rather than a result at rank 4. `angle`, degenerate at rank 1, is 22.1423° at rank 4 and between
22.0° and 35.8° across the forced ranks in the primary cohort, with p-values
0.152–0.455: non-rejections that are uninformative about orientation, not evidence that the
sexes' directions agree. `shape` is not rejected at any rank, and its non-rejections are uninformative too. The
sensitivity cohort (n = 74) agrees at every rank. MOTCO's forced-rank operating characteristics were measured
without covariates and, at this design, only at rank 4, so the forced-rank p-values are carried as a sensitivity
description, not as calibrated tests.

## 1. Provenance

| item | value |
|---|---|
| MOTCO | 0.8.0, tag `v0.8.0`, commit `979657c3d0e6c9ad4d4bff0dc595a18df854ca61` (pinned in `uv.lock`) |
| entry point | `motco.analysis.analyze_trajectories` with `TrajectoryAnalysisParams(forced_components=k)` + `write_analysis_outputs`, through `sead_sexdiff.data.run_rank_sensitivity` (the same `_analyze_cohort` helper as the `analyze` step) |
| sead_sexdiff | `9cd07a1` (the change's implementation uncommitted at run time; both sweep manifests record `dirty: true`) |
| machine | `nixos-home` (NixOS workstation); the reference analyses were produced on the laptop workstation and pulled from the cluster |
| reference analyses | `results/motco_analysis/<cohort>/`, the CV-selected analysis of `reports/motco-astrocyte-2026-10-19.md` (sead_sexdiff `c2f25e3`, `dirty: false`): rank 1 in both cohorts; `result.json` sha256 `ca92876a1a61…` (cohort_30), `1e7eefe3e63e…` (cohort_50); `analysis_manifest.json` sha256 `17b2150ecba0…`, `43fe406690fa…`; read and checked, never rerun or modified |
| inputs | the export of the primary analysis: export manifest sha256 `18f2289e0d4c…` (cohort_30), `dd14a6e7d558…` (cohort_50); the sweep checked all five files of each cohort against the hashes in the reference `analysis_manifest.json` before running |
| forced ranks | `ANALYSIS_FORCED_RANKS = (2, 3, 4, 6, 9, 12)`, a constant: 2 = `n_stages − 1`, the lowest rank at which `angle` carries orientation information and the rank MOTCO's stage-supervised rule selects at three stages in its operating studies; 3 the next rank; 4, 6, 9, 12 the steps of MOTCO's latent-rank ladder, with 4 = 2(K − 1) the forced rank with a Type I check at this design |
| rank 1 | rank 1 cannot be forced (MOTCO accepts forced_components >= 2 only); the CV-selected analysis is the reference row |
| measurement | as the primary analysis: blocks `rna`, `atac` in identity units; group `sex`; stages 0, 1, 2; covariates `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac` (additive, `missing_covariates = error`); net-displacement orientation; no attribution; cross-validation skipped (`component_selection = forced`) |
| permutations, seed | 999 RRPP permutations, seed 20261019, at every rank and equal to the primary analysis's; `n_jobs` 1 (both cohorts) |
| α | 0.05, rejected when p ≤ α |
| software | python 3.12.15, numpy 2.3.5, scipy 1.15.3, scikit-learn 1.8.0, pandas 2.3.3 |

Every forced-rank directory (`rank_<kk>/`) holds the same seven files as a CV-selected analysis directory, with
`measurement.cv = null` in its manifest and a `cli_equivalent` that reproduces it with `motco analyze
--forced-components k` (section 5).

## 2. Primary cohort: cohort_30 (≥ 30 astrocyte nuclei in both omics, n = 80)

Reference: `results/motco_analysis/cohort_30` (CV-selected rank 1, CV mean AUROC
0.6313; `result.json` sha256 `ca92876a1a61…`,
`analysis_manifest.json` sha256 `17b2150ecba0…`; export
manifest sha256 `18f2289e0d4c…`). Nothing is rejected at α = 0.05 at
the CV rank or at any forced rank. The smallest forced-rank p-value is `delta` at rank 4
(p = 0.052). Findings recorded by the sweep: CV-selected rank 1 is not in the ladder (2, 3, 4, 6, 9, 12): rank 1 cannot be forced (MOTCO accepts forced_components >= 2 only); the CV-selected analysis is the reference row.

Forced-rank sweep of the covariate-adjusted trajectory analysis. The reference row (`cv`) is the CV-selected analysis at rank 1, which stays the primary result; the forced ranks 2, 3, 4, 6, 9, 12 are read as sensitivity to it, and no rank is chosen from the results. 999 RRPP permutations with seed 20261019 at every rank (the same permutation indices throughout); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative at every rank; a `delta` rejection is a difference, not specifically of size.

### delta (pasted from `rank_sensitivity.md`)

| selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|:---:|---:|---:|---:|
| cv | 1 | 14.941 | 0.3 | no | 11.2868 | 8.27982 | 26.5775 |
| forced | 2 | 3.91847 | 0.744 | no | 9.63089 | 7.18465 | 23.2427 |
| forced | 3 | 10.2517 | 0.338 | no | 8.48407 | 6.74633 | 21.5709 |
| forced | 4 | 19.247 | 0.052 | no | 7.78873 | 6.01757 | 19.3458 |
| forced | 6 | 6.08738 | 0.369 | no | 5.69703 | 4.5972 | 14.8882 |
| forced | 9 | 1.79797 | 0.778 | no | 5.26197 | 4.29954 | 13.432 |
| forced | 12 | 1.15591 | 0.866 | no | 5.32251 | 4.29312 | 13.705 |

### angle (degrees) (pasted from `rank_sensitivity.md`)

| selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|:---:|---:|---:|---:|
| cv | 1 | 0 | 1 | no | 3.06306 | 23.2919 | 0 |
| forced | 2 | 22.832 | 0.241 | no | 15.7454 | 12.6226 | 40.5824 |
| forced | 3 | 23.0297 | 0.316 | no | 18.8414 | 10.5607 | 39.0162 |
| forced | 4 | 22.1423 | 0.356 | no | 19.6796 | 9.68581 | 37.4196 |
| forced | 6 | 21.9664 | 0.455 | no | 22.0223 | 8.9536 | 38.8732 |
| forced | 9 | 31.5192 | 0.21 | no | 25.4414 | 8.08815 | 40.3043 |
| forced | 12 | 35.8438 | 0.152 | no | 28.4101 | 7.42662 | 41.5921 |

### shape (pasted from `rank_sensitivity.md`)

| selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|:---:|---:|---:|---:|
| cv | 1 | 0.608832 | 0.333 | no | 0.489849 | 0.355498 | 1.18287 |
| forced | 2 | 0.295865 | 0.445 | no | 0.290846 | 0.158807 | 0.589015 |
| forced | 3 | 0.121943 | 0.789 | no | 0.22084 | 0.117102 | 0.427894 |
| forced | 4 | 0.0296112 | 0.975 | no | 0.162616 | 0.08214 | 0.307316 |
| forced | 6 | 0.153228 | 0.224 | no | 0.111023 | 0.0578495 | 0.215647 |
| forced | 9 | 0.147921 | 0.161 | no | 0.0957448 | 0.0508231 | 0.187395 |
| forced | 12 | 0.131817 | 0.208 | no | 0.0925515 | 0.0492176 | 0.181653 |

### Runs (pasted from `rank_sensitivity.md`)

| selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---:|---:|---:|---|
| cv | 1 | 0.6313 | 148.7 | `results/motco_analysis/cohort_30` |
| forced | 2 | n/a | 0.4 | `rank_02` |
| forced | 3 | n/a | 0.5 | `rank_03` |
| forced | 4 | n/a | 0.7 | `rank_04` |
| forced | 6 | n/a | 0.5 | `rank_06` |
| forced | 9 | n/a | 0.5 | `rank_09` |
| forced | 12 | n/a | 0.5 | `rank_12` |

### Descriptive checks

These checks describe the geometry behind the tables and carry no inference.

- **The first component does not move.** `pls_0` in every forced rank's `latent_scores.csv` equals the CV run's
  single component (`numpy.allclose`, same sign, largest absolute difference 2.8e-14 over the six ranks). PLS
  components are fitted sequentially, so forcing a higher rank appends components to the CV-selected space rather
  than replacing it.
- **Stage means along `pls_0` are therefore identical at every rank** (LS-means at the pooled covariate means,
  `stage_means.csv`), and the female path stays the longer one: net displacement (stage 2 − stage 0) along
  `pls_0` 32.63 for females and 17.69 for males, both monotone from stage 0 to 2, at the CV rank and
  at every forced rank.

| group | stage 0 | stage 1 | stage 2 | net displacement (2 − 0) along `pls_0` |
|---|---:|---:|---:|---:|
| Female | -11.75 | -5.33 | 20.89 | 32.63 |
| Male | -24.11 | -10.83 | -6.42 | 17.69 |

- **Null width across ranks.** The RRPP null of every statistic narrows as the rank grows, while the `angle`
  null's 95 % quantile stays between 37.4° and 41.6° at the forced ranks. At the CV rank the
  `angle` null is degenerate (0° or 180° draws: q95 0, sd 23.2919). The observed
  `angle` at the forced ranks (22.0° to 35.8°) lies below that quantile at
  every rank.

| selection | rank | `delta` null sd | `angle` null sd | `angle` null q95 | `shape` null sd |
|---|---:|---:|---:|---:|---:|
| cv | 1 | 8.27982 | 23.2919 | 0 | 0.355498 |
| forced | 2 | 7.18465 | 12.6226 | 40.5824 | 0.158807 |
| forced | 3 | 6.74633 | 10.5607 | 39.0162 | 0.117102 |
| forced | 4 | 6.01757 | 9.68581 | 37.4196 | 0.08214 |
| forced | 6 | 4.5972 | 8.9536 | 38.8732 | 0.0578495 |
| forced | 9 | 4.29954 | 8.08815 | 40.3043 | 0.0508231 |
| forced | 12 | 4.29312 | 7.42662 | 41.5921 | 0.0492176 |

## 3. Sensitivity cohort: cohort_50 (≥ 50 astrocyte nuclei in both omics, n = 74)

Reference: `results/motco_analysis/cohort_50` (CV-selected rank 1, CV mean AUROC
0.6616; `result.json` sha256 `1e7eefe3e63e…`,
`analysis_manifest.json` sha256 `43fe406690fa…`; export
manifest sha256 `dd14a6e7d558…`). Nothing is rejected at α = 0.05 at
the CV rank or at any forced rank. The smallest forced-rank p-value is `delta` at rank 6
(p = 0.208). Findings recorded by the sweep: CV-selected rank 1 is not in the ladder (2, 3, 4, 6, 9, 12): rank 1 cannot be forced (MOTCO accepts forced_components >= 2 only); the CV-selected analysis is the reference row.

Forced-rank sweep of the covariate-adjusted trajectory analysis. The reference row (`cv`) is the CV-selected analysis at rank 1, which stays the primary result; the forced ranks 2, 3, 4, 6, 9, 12 are read as sensitivity to it, and no rank is chosen from the results. 999 RRPP permutations with seed 20261019 at every rank (the same permutation indices throughout); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative at every rank; a `delta` rejection is a difference, not specifically of size.

### delta (pasted from `rank_sensitivity.md`)

| selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|:---:|---:|---:|---:|
| cv | 1 | 18.7091 | 0.213 | no | 11.6918 | 9.218 | 30.0679 |
| forced | 2 | 6.01569 | 0.596 | no | 9.78355 | 7.88773 | 25.288 |
| forced | 3 | 9.36583 | 0.496 | no | 11.1486 | 8.59682 | 27.3245 |
| forced | 4 | 14.1822 | 0.217 | no | 9.01492 | 6.56113 | 21.4611 |
| forced | 6 | 10.7046 | 0.208 | no | 6.81048 | 5.21583 | 16.7231 |
| forced | 9 | 6.51697 | 0.386 | no | 6.1817 | 4.92541 | 15.8524 |
| forced | 12 | 5.88663 | 0.436 | no | 6.20365 | 4.94692 | 15.8513 |

### angle (degrees) (pasted from `rank_sensitivity.md`)

| selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|:---:|---:|---:|---:|
| cv | 1 | 0 | 1 | no | 1.98198 | 18.7931 | 0 |
| forced | 2 | 23.1596 | 0.222 | no | 15.1566 | 12.4608 | 40.0109 |
| forced | 3 | 23.4732 | 0.299 | no | 19.2555 | 11.5529 | 42.4812 |
| forced | 4 | 24.5502 | 0.292 | no | 20.4756 | 10.6534 | 42.0583 |
| forced | 6 | 25.2911 | 0.33 | no | 22.5825 | 9.44256 | 40.8348 |
| forced | 9 | 30.5012 | 0.251 | no | 25.7816 | 8.57882 | 42.7693 |
| forced | 12 | 32.5938 | 0.279 | no | 28.7713 | 7.98762 | 43.9592 |

### shape (pasted from `rank_sensitivity.md`)

| selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|:---:|---:|---:|---:|
| cv | 1 | 0.56471 | 0.377 | no | 0.495877 | 0.343673 | 1.159 |
| forced | 2 | 0.284932 | 0.484 | no | 0.296547 | 0.162566 | 0.599664 |
| forced | 3 | 0.16487 | 0.739 | no | 0.251509 | 0.126233 | 0.483667 |
| forced | 4 | 0.110409 | 0.719 | no | 0.168021 | 0.0888066 | 0.33709 |
| forced | 6 | 0.0845199 | 0.652 | no | 0.116727 | 0.062907 | 0.231105 |
| forced | 9 | 0.0826203 | 0.593 | no | 0.103863 | 0.0566925 | 0.208077 |
| forced | 12 | 0.0785119 | 0.598 | no | 0.100357 | 0.0553857 | 0.200242 |

### Runs (pasted from `rank_sensitivity.md`)

| selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---:|---:|---:|---|
| cv | 1 | 0.6616 | 151.6 | `results/motco_analysis/cohort_50` |
| forced | 2 | n/a | 0.4 | `rank_02` |
| forced | 3 | n/a | 0.4 | `rank_03` |
| forced | 4 | n/a | 0.4 | `rank_04` |
| forced | 6 | n/a | 0.5 | `rank_06` |
| forced | 9 | n/a | 0.5 | `rank_09` |
| forced | 12 | n/a | 0.5 | `rank_12` |

### Descriptive checks

These checks describe the geometry behind the tables and carry no inference.

- **The first component does not move.** `pls_0` in every forced rank's `latent_scores.csv` equals the CV run's
  single component (`numpy.allclose`, same sign, largest absolute difference 2.8e-14 over the six ranks). PLS
  components are fitted sequentially, so forcing a higher rank appends components to the CV-selected space rather
  than replacing it.
- **Stage means along `pls_0` are therefore identical at every rank** (LS-means at the pooled covariate means,
  `stage_means.csv`), and the female path stays the longer one: net displacement (stage 2 − stage 0) along
  `pls_0` 37.15 for females and 18.44 for males, both monotone from stage 0 to 2, at the CV rank and
  at every forced rank.

| group | stage 0 | stage 1 | stage 2 | net displacement (2 − 0) along `pls_0` |
|---|---:|---:|---:|---:|
| Female | -16.43 | -3.39 | 20.72 | 37.15 |
| Male | -22.95 | -6.94 | -4.51 | 18.44 |

- **Null width across ranks.** The RRPP null of every statistic narrows as the rank grows, while the `angle`
  null's 95 % quantile stays between 40.0° and 44.0° at the forced ranks. At the CV rank the
  `angle` null is degenerate (0° or 180° draws: q95 0, sd 18.7931). The observed
  `angle` at the forced ranks (23.2° to 32.6°) lies below that quantile at
  every rank.

| selection | rank | `delta` null sd | `angle` null sd | `angle` null q95 | `shape` null sd |
|---|---:|---:|---:|---:|---:|
| cv | 1 | 9.218 | 18.7931 | 0 | 0.343673 |
| forced | 2 | 7.88773 | 12.4608 | 40.0109 | 0.162566 |
| forced | 3 | 8.59682 | 11.5529 | 42.4812 | 0.126233 |
| forced | 4 | 6.56113 | 10.6534 | 42.0583 | 0.0888066 |
| forced | 6 | 5.21583 | 9.44256 | 40.8348 | 0.062907 |
| forced | 9 | 4.92541 | 8.57882 | 42.7693 | 0.0566925 |
| forced | 12 | 4.94692 | 7.98762 | 43.9592 | 0.0553857 |

## 4. Covariates and feature filtering

Unchanged from the primary analysis, by hash: the sweep read the same five export files and entered the same six
covariates in the same order, and the reference `analysis_manifest.json` hashes matched. For the covariate block
(design columns, reference row, coefficients, dropped donors: none) see sections 2 and 3 of
`reports/motco-astrocyte-2026-10-19.md`, and for the feature filtering (chrX/Y/M removal, CPM ≥ 1 in ≥ k donors,
TMM, log2-CPM, top 5,000 by variance per block) its section 4. The latent space is fitted without the covariates
at every rank; only the trajectory regression and the RRPP reduced model carry them.

## 5. Reproducibility checks

- **Refactor.** The `analyze` step's per-cohort body moved into the helper the sweep shares. On this machine the
  pre-refactor code (`git archive HEAD`) and the refactored code give byte-identical `latent_scores.csv`,
  `pls_loadings.csv`, `stage_means.csv`, `analysis_metadata.csv` and `summary.md` for `cohort_50`, and
  `result.json` / `analysis_manifest.json` equal apart from `runtime_seconds`, `created_at` and `git`.
- **Across machines.** Rerunning `analyze --cohort cohort_50` here reproduces the laptop's `summary.md`, p-values
  and rank exactly, but not its CSVs byte for byte: `analysis_metadata.csv` differs by one ulp in one
  `log10_n_atac` value (`numpy.log10` differs between the two machines at identical package versions), and the
  statistics and null summaries follow at ≤ 1e-13 relative. The forced runs reported here were computed on this
  machine against the laptop's reference records; no decision at α = 0.05 depends on that.
- **Rerun.** `cohort_50` swept again with `--force` into a scratch directory: `rank_sensitivity.csv` and, at every
  rank, `latent_scores.csv`, `pls_loadings.csv` and `stage_means.csv` byte-identical; `statistics`, `p_values`
  and `null_summary` equal.
- **MOTCO's CLI.** `cohort_30/rank_04` reproduced with its manifest's `cli_equivalent` (`motco analyze` on
  `analysis_metadata.csv` with the six `--covariate` flags and `--forced-components 4`): `statistics`,
  `p_values`, `null_summary`, `selected_lv = 4` and `component_selection = forced` equal.
- **Reference row.** In both cohorts the `cv` row of `rank_sensitivity.csv` equals the kept `result.json` cell
  for cell, and the reference directories are byte-identical before and after the sweep.
- **Monte Carlo spread.** Not re-measured here. The primary report's 30-seed check on `cohort_50` put the SD of
  p-values near p = 0.2–0.35 at about 0.015; the forced-rank p-values differ between ranks by more than that
  where they differ at all, and no decision depends on them.

## 6. Interpretation

The binding rules come from `../MOTCO/docs/roadmap.md` (Phase 6; Next changes item 1) and from the spec of this
change. Applied to these results:

- **The CV-selected analysis is primary.** It selected rank 1 in both cohorts and rejected nothing. The forced
  ranks are its sensitivity analysis; no rank is chosen after seeing the results, and no forced-rank row replaces
  the primary row.
- **No rank dependence of any decision.** No statistic is rejected at any rank in either cohort, so the rejection
  status is the same at the CV rank and at every forced rank. The observed statistics differ between ranks (they
  are measured in different latent spaces and are compared only through their p-values), and the p-values move
  with the rank; `delta` at rank 4 in `cohort_30` (p = 0.052) comes closest to α. Had it been rejected, the
  report would have recorded rank dependence of the sensitivity analysis and kept the CV-selected non-rejection
  as the primary result.
- **No evidence of a difference at any rank, not evidence of agreement.** The sweep does not show that the sexes'
  trajectories agree in size, direction or shape at any rank: `angle` and `shape` non-rejections stay
  uninformative at every rank (section 7), and the adjusted analysis carries no power claim.
- **`angle` is defined at the forced ranks but still uninformative.** At rank 1 the angle is 0° by construction. At
  the forced ranks it is 22.0°–35.8° (`cohort_30`) and 23.2°–32.6° (`cohort_50`), below the null's 95 %
  quantile at every rank. These non-rejections say nothing about whether the directions agree.
- **A `delta` rejection would have been a difference, not specifically of size**, because direction and shape
  differences also move `delta`.

## 7. Limitations

- **Type I evidence at a forced rank covers rank 4 only, without covariates.** MOTCO's net-displacement Type I
  study (`docs/reports/net-displacement-type-i-2026-10-13.md` §2, at the v0.8.0 tag) holds size at
  `pls_forced4 × nd` at this design (n = 80, K = 3, cells 11/10/28 vs 9/10/12: anchor `delta` 0.006, `angle`
  0.048, `shape` 0.006), with a narrower and less dispersed `angle` null than CV (q95 32.8°, IQR 5.9, §4). There
  is no Type I measurement at ranks 2, 3, 6, 9 or 12 at this design, and none at any forced rank with
  covariates: the covariate Type I check (`docs/reports/covariate-type-i-2026-10-16.md`) ran under CV only.
- **Forced 2(K − 1) trades `shape` for `angle`.** In MOTCO's orientation probe
  (`docs/reports/orientation-measurement-probe-2026-10-12.md` §3, §6) `pls_forced` at rank 4 recovers
  orientation (`angle` power 0.99 at e = 1) but drops `shape` power from 0.19 to 0.01 and leaks `angle` into the
  magnitude cell (0.15). A `shape` non-rejection at rank 4 is therefore even less informative than at the CV rank.
- **Ranks above `n_stages − 1` widen the `angle` null.** MOTCO's latent-rank ladder
  (`docs/reports/latent-rank-ladder-2026-09-08.md`, n = 1200, K = 4) found no `angle` gain above `n_stages − 1`,
  `shape` and `delta` power losses, a wider `angle` null, and no anchor inflation at any rank. The last point is
  cited as what it is, an unadjusted result at a larger n, not as a guarantee here. In these data the `angle`
  null's q95 at the forced ranks (37.4°–41.6° in `cohort_30`) is wider than the 32.8° of MOTCO's
  forced-4 study at this design, a descriptive comparison across different data.
- **Twelve model columns at cells of 9–28 donors.** Each forced run fits an OLS with 80 (or 74) rows, six sex ×
  stage columns and six covariates on up to twelve latent columns; MOTCO did not report a rank-deficient design
  at any rank, and the per-statistic null sd is in the tables.
- **Same seed at every rank.** Chosen so that between-rank differences are measurement differences; it also
  means the ranks share one permutation draw, so the sweep says nothing about Monte Carlo spread.
- **Uncommitted implementation at run time.** The sweep ran at `9cd07a1` with the change's code uncommitted; the
  `--force` rerun check above shows the outputs are reproducible, and a post-commit rerun will carry a clean
  `git` field.
- **Everything the primary report lists** (additive linear adjustment, covariate-blind basis, six covariate
  columns against three validated, no power claim, structural analogue, fixed preprocessing choices, Monte Carlo
  error and `n_jobs`) applies unchanged.

## 8. Next steps

- Preprocessing sensitivity exports (other N, keeping chrX, residualized features).
- Pathway/module enrichment, exploratory without an orientation result.
- MOTCO side: whether CV may select rank 1, where `angle` is degenerate and `forced_components` cannot reproduce
  the fit; whether parallel RRPP should reproduce the serial draws.

## 9. Reproduction

```
uv sync --locked
uv run sead-sexdiff rank-sensitivity --cohort all          # needs results/motco_input/ and results/motco_analysis/
```

Outputs go to `results/motco_rank_sensitivity/<cohort>/`. A single rank is reproduced with MOTCO's own CLI by
running the `cli_equivalent.command` of `rank_<kk>/analysis_manifest.json` from the repository root with `OUT_DIR`
replaced by a new directory. The tables above are `rank_sensitivity.md`; the same values with full precision are in
`rank_sensitivity.csv` (shortest round-trip floats; read with `float_precision="round_trip"`).

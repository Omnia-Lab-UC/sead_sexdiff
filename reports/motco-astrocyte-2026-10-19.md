# MOTCO astrocyte trajectory analysis: sex × ADNC stage — 2026-10-19

> Covariate-adjusted comparison of female and male astrocyte trajectories across ADNC stages in the SEA-AD MTG
> (OpenSpec change `mtg-astrocyte-motco-analysis`, spec `motco-trajectory-analysis`). Every statistic, p-value,
> null summary, rank and cell size below is pasted by script from the generated `summary.md` of each cohort, which
> is rendered from MOTCO's `result.json` in `results/motco_analysis/<cohort>/`. The latent stage means come from
> `stage_means.csv` and the feature counts from the export manifests (`results/motco_input/<cohort>/manifest.json`).

**Run:** `sead-sexdiff analyze --cohort all --force` on the laptop workstation, 2026-10-19, at commit
`c2f25e3`, `n_jobs = 1` (149 s for cohort_30, 152 s for cohort_50). **Result:** nothing is rejected at α = 0.05 in either cohort. Cross-validation selects a
one-dimensional latent space (rank 1) in both cohorts, so `angle` is 0° by construction there and says nothing
about orientation.

**In one paragraph.** For each cohort, MOTCO 0.8.0's `analyze_trajectories` fitted one stage-supervised PLS-DA
latent space on the RNA and ATAC blocks (5,000 + 5,000 log2-CPM features). It then compared the female and male
trajectories across stages 0 (Not AD + Low), 1 (Intermediate) and 2 (High), adjusting for age at death, PMI, RIN,
years of education and log10 astrocyte nuclei per omic, with p-values from 999 RRPP permutations. In the primary
cohort (≥ 30 nuclei, n = 80) and in the sensitivity cohort (≥ 50 nuclei, n = 74), none of `delta`, `angle` or
`shape` is rejected. The latent space separates the stages only weakly (CV mean AUROC 0.63 and 0.66) and has one
component. In one dimension the two sexes' net displacements are either parallel (0°) or antiparallel (180°), and
here they are parallel. These results are not evidence that the sexes' trajectories are the same. Under the binding
interpretation rules the `angle` and `shape` non-rejections are uninformative at this n, and the adjusted analysis
carries no power claim. Rank 1 is also below the ranks (2 or more) at which MOTCO's operating characteristics were
measured.

## 1. Provenance

| item | value |
|---|---|
| MOTCO | 0.8.0, tag `v0.8.0`, commit `979657c3d0e6c9ad4d4bff0dc595a18df854ca61` (pinned in `uv.lock`) |
| entry point | `motco.analysis.analyze_trajectories` + `write_analysis_outputs` (not `motco plsr` + `motco de`) |
| sead_sexdiff | `c2f25e3` (branch `mtg-astrocyte-motco-analysis`; both analysis manifests record it with `dirty: false`) |
| inputs | export of 2026-10-16 (`ab3423a`); export manifest sha256 `18f2289e0d4c…` (cohort_30), `dd14a6e7d558…` (cohort_50); the analysis manifests hold the sha256 of all five input files |
| measurement | net-displacement orientation; double-CV rank (`cv2_splits` 4, `cv1_splits` 3, `n_repeats` 5, `max_components` 20, `random_state` 1203); no attribution |
| permutations, seed | 999 RRPP permutations, seed 20261019, `n_jobs` 1 (both cohorts) |
| α | 0.05, rejected when p ≤ α (the rule of MOTCO's covariate Type I check) |
| software | python 3.12.15, numpy 2.3.5, scipy 1.15.3, scikit-learn 1.8.0, pandas 2.3.3 |

The model: blocks `rna` and `atac` in identity units (log2-CPM as exported), group `sex` (Female, Male), stage
`stage` in the declared order 0, 1, 2. Covariates enter additively in this order: `age_at_death`, `pmi`, `rin`,
`years_education`, `log10_n_rna`, `log10_n_atac`. The last two are log10 of the astrocyte nuclei counts per omic and
are derived at analysis time. `education_level` is not entered (categorical, by decision) and neither is
`age_capped` (constant: no capped ages in this release). The missing-covariate policy is `error`. No value is
missing and nothing is imputed. There is no unadjusted run.

## 2. Primary cohort: cohort_30 (≥ 30 astrocyte nuclei in both omics, n = 80)

MOTCO 0.8.0; n = 80; groups Female vs Male; stages 0 → 1 → 2; orientation `net_displacement`; 999 RRPP permutations, seed 20261019; α = 0.05 (rejected when p ≤ α).

### Statistics

| statistic | observed | p (RRPP) | rejected | null mean | null sd | null q95 | note |
|---|---:|---:|:---:|---:|---:|---:|---|
| delta | 14.941 | 0.3 | no | 11.2868 | 8.27982 | 26.5775 |  |
| angle | 0 | 1 | no | 3.06306 | 23.2919 | 0 | degrees; CV-selected rank 1 (at rank 1 the angle is 0° or 180° by construction) |
| shape | 0.608832 | 0.333 | no | 0.489849 | 0.355498 | 1.18287 |  |

### Cell sizes

| group | stage 0 | stage 1 | stage 2 | total |
|---|---:|---:|---:|---:|
| Female | 11 | 10 | 28 | 49 |
| Male | 9 | 10 | 12 | 31 |

### Latent space

- PLS-DA rank: 1 (`cv`)
- cv_mean_auroc: 0.6313
- CV: cv2_splits 4, cv1_splits 3, n_repeats 5 (5 completed), max_components 20, random_state 1203
- blocks: rna (5000 features, identity), atac (5000 features, identity)

### Covariates

- entered (in order): `age_at_death` (numeric), `pmi` (numeric), `rin` (numeric), `years_education` (numeric), `log10_n_rna` (numeric), `log10_n_atac` (numeric)
- design columns: 6 (age_at_death, pmi, rin, years_education, log10_n_rna, log10_n_atac)
- missing policy: `error`; dropped donors: 0
- LS-means evaluated at the reference row (`pooled_mean`):

| design column | reference value |
|---|---:|
| age_at_death | 89.025 |
| pmi | 7.0935 |
| rin | 8.15087 |
| years_education | 16.1375 |
| log10_n_rna | 2.82437 |
| log10_n_atac | 2.22423 |

### Latent stage means (descriptive)

LS-means on the single PLS component, evaluated at the pooled covariate means (`stage_means.csv`). They show the
geometry behind the statistics and carry no inference.

| group | stage 0 | stage 1 | stage 2 | net displacement (2 − 0) |
|---|---:|---:|---:|---:|
| Female | -11.75 | -5.33 | 20.89 | 32.63 |
| Male | -24.11 | -10.83 | -6.42 | 17.69 |

## 3. Sensitivity cohort: cohort_50 (≥ 50 astrocyte nuclei in both omics, n = 74)

MOTCO 0.8.0; n = 74; groups Female vs Male; stages 0 → 1 → 2; orientation `net_displacement`; 999 RRPP permutations, seed 20261019; α = 0.05 (rejected when p ≤ α).

### Statistics

| statistic | observed | p (RRPP) | rejected | null mean | null sd | null q95 | note |
|---|---:|---:|:---:|---:|---:|---:|---|
| delta | 18.7091 | 0.213 | no | 11.6918 | 9.218 | 30.0679 |  |
| angle | 0 | 1 | no | 1.98198 | 18.7931 | 0 | degrees; CV-selected rank 1 (at rank 1 the angle is 0° or 180° by construction) |
| shape | 0.56471 | 0.377 | no | 0.495877 | 0.343673 | 1.159 |  |

### Cell sizes

| group | stage 0 | stage 1 | stage 2 | total |
|---|---:|---:|---:|---:|
| Female | 10 | 9 | 25 | 44 |
| Male | 9 | 9 | 12 | 30 |

### Latent space

- PLS-DA rank: 1 (`cv`)
- cv_mean_auroc: 0.6616
- CV: cv2_splits 4, cv1_splits 3, n_repeats 5 (5 completed), max_components 20, random_state 1203
- blocks: rna (5000 features, identity), atac (5000 features, identity)

### Covariates

- entered (in order): `age_at_death` (numeric), `pmi` (numeric), `rin` (numeric), `years_education` (numeric), `log10_n_rna` (numeric), `log10_n_atac` (numeric)
- design columns: 6 (age_at_death, pmi, rin, years_education, log10_n_rna, log10_n_atac)
- missing policy: `error`; dropped donors: 0
- LS-means evaluated at the reference row (`pooled_mean`):

| design column | reference value |
|---|---:|
| age_at_death | 89.3108 |
| pmi | 7.23847 |
| rin | 8.14784 |
| years_education | 16 |
| log10_n_rna | 2.832 |
| log10_n_atac | 2.2729 |

### Latent stage means (descriptive)

| group | stage 0 | stage 1 | stage 2 | net displacement (2 − 0) |
|---|---:|---:|---:|---:|
| Female | -16.43 | -3.39 | 20.72 | 37.15 |
| Male | -22.95 | -6.94 | -4.51 | 18.44 |

## 4. Feature filtering

The blocks MOTCO analysed are the export's TMM log2-CPM features (archived change
`2026-10-16-mtg-astrocyte-motco-export`). Per cohort and omic, without using sex, stage or covariates:

1. drop chrX, chrY and chrM features;
2. keep features with CPM ≥ 1 in at least `k` donors, `k` = the smallest sex × stage cell;
3. compute TMM factors and log2-CPM with prior count 2;
4. keep the 5,000 most variable features (equal N per block, because MOTCO z-scores features and concatenates the
   blocks without weights).

| cohort | omic | input | after chrX/Y/M | after CPM ≥ 1 in ≥ k | kept (top N) | k (smallest cell) |
|---|---|---:|---:|---:|---:|---|
| cohort_30 | rna | 36,601 | 35,329 | 20,387 | 5,000 | 9 (Male/stage 0) |
| cohort_30 | atac | 218,882 | 212,816 | 128,916 | 5,000 | 9 (Male/stage 0) |
| cohort_50 | rna | 36,601 | 35,329 | 20,238 | 5,000 | 9 (Female/stage 1) |
| cohort_50 | atac | 218,882 | 212,816 | 124,818 | 5,000 | 9 (Female/stage 1) |

Features removed by chromosome (identical in both cohorts): RNA: chrX 1,148, chrY 111, chrM 13; ATAC: chrX 5,886, chrY 180, chrM 0.

29 of the 36,601 RNA genes are absent from GENCODE v32 basic (24 on unplaced scaffolds, 5 retired ids). They are
kept with chromosome `unmapped`, which the export manifests record as a finding.

## 5. Reproducibility checks

- **Commit rerun.** The outputs reported here come from a rerun of both cohorts at the committed revision. Its
  CSVs and `summary.md` are byte-identical to the run before the commit, and its records are equal apart from
  `runtime_seconds`.
- **Rerun.** `cohort_50` rerun with `--force`: `latent_scores.csv`, `pls_loadings.csv`, `stage_means.csv` and
  `analysis_metadata.csv` byte-identical; `statistics`, `p_values`, `null_summary` and `selected_lv` equal.
- **MOTCO's CLI.** `cohort_30` reproduced with the manifest's `cli_equivalent` (`motco analyze` on
  `analysis_metadata.csv` with the six `--covariate` flags): statistics, p-values, null summaries, rank, covariate
  block and the three MOTCO CSVs equal.
- **`n_jobs`.** `cohort_50` at `n_jobs = 4` gives the same statistics, rank, latent scores and loadings, but a
  different permutation null: p(`delta`) 0.199 and p(`shape`) 0.319 instead of the reported 0.213 and 0.377, and
  p(`angle`) 1 either way. MOTCO 0.8.0's parallel RRPP seeds one stream per worker, while the serial path uses
  one, so the two runs use different permutation draws of the same null. The results are reproducible at a fixed
  `n_jobs` only, and this report uses `n_jobs = 1`.
- **Monte Carlo spread.** Rerunning only the RRPP step on the saved `cohort_50` latent scores with MOTCO's own
  design matrices reproduces both runs exactly. Over 30 further seeds (`n_jobs = 1`), p(`delta`) ranges
  0.184–0.239 (SD 0.016) and p(`shape`) 0.317–0.378 (SD 0.013). Both `n_jobs` results lie within that spread, so
  the difference is Monte Carlo variation. The reported p(`shape`) 0.377 sits at the top of the range. No decision
  at α = 0.05 depends on the draws.

## 6. Interpretation

The binding rules come from `../MOTCO/docs/roadmap.md` (Phase 6; Next changes item 1), the Phase 6 small-n report
§8, the net-displacement power report §7 and the covariate Type I report §5. Applied to these results:

- **`delta` is not rejected** in either cohort. Had it been rejected, it would be evidence that the sexes'
  trajectories differ, not specifically in size, because direction and shape differences also reject `delta`.
  The non-rejection is not evidence against a size difference. MOTCO's Phase 6 statement that a `delta`
  non-rejection argues against size differences of about 25 % or more was measured without covariates and in
  latent spaces of rank 2 or more. The covariate Type I report makes no power claim for the adjusted analysis, and
  this latent space has rank 1.
  Descriptively, the female path along `pls_0` is longer than the male path in both cohorts (stage means above),
  but the permutation null is wide enough that this is not distinguishable from chance.
- **`angle` is uninformative about orientation.** `angle` is the angle between the sexes' net displacement vectors
  (stage 2 mean − stage 0 mean), so it reads only the endpoint stages. The CV-selected rank is 1 in both cohorts.
  In one dimension the two net displacements can only be parallel (0°) or antiparallel (180°). The observed 0°
  says that both sexes move in the same direction along the single component. It does not say that their
  directions of change agree in the 10,000-feature space, and its non-rejection (p = 1) is uninformative, not
  evidence of agreement. Even at higher ranks MOTCO finds `angle` non-rejections uninformative at n = 80: power is
  0.084–0.178 at e ≤ 0.75 and 0.280 at the construction maximum when CV retains three or fewer components.
- **`shape` is not rejected**, and a `shape` non-rejection is uninformative at n = 80 (power at most 0.228 at
  the construction maximum, measured without covariates).
- **Cohort roles.** `cohort_30` is the primary analysis and `cohort_50` the sensitivity analysis. The two agree:
  rank 1, nothing rejected, stage means increasing from stage 0 to 2 in both sexes, with the female path the
  longer one. MOTCO's operating study expects the smaller
  cohort to leave the `delta` statements unchanged and to make `angle` and `shape` weaker still.

**What this analysis supports:** after adjusting for the six covariates, it finds no evidence that female and
male MTG astrocyte trajectories across ADNC stages differ in the RNA + ATAC latent space. **What it does not
support:** any claim that they are the same, in size, direction or shape.

## 7. Limitations

- **Rank 1 is outside the simulated range.** In MOTCO's Phase 6 operating studies CV selected rank 2 or 3 in
  essentially every replicate, and more (up to 9) only at the strongest orientation construction. The covariate
  Type I check selected rank 2 or 3 in all 500 replicates. The reported operating characteristics therefore
  describe latent spaces of rank 2 or more. The weak stage separation (CV mean AUROC 0.63 / 0.66) is consistent
  with CV retaining a single component. MOTCO is also inconsistent here: its CV tries ranks 1 to
  `max_components` − 1 (1–19), but `forced_components` accepts only 2 and up, so the rank-1 fit cannot be
  reproduced as a forced-rank diagnostic, and the forced-rank sweep starts at 2.
- **Additive, linear adjustment.** The covariate Type I check covers covariate effects that are additive, linear
  and correctly specified. Non-linear effects, or covariates that interact with sex or stage, are outside it.
- **Covariate-blind basis.** The PLS latent space is fitted without the covariates. Variation aligned with a
  covariate that also tracks stage can take up latent capacity, here the single component. The covariates enter
  only the trajectory regression and the RRPP reduced model.
- **More covariate columns than validated.** Six numeric covariates give six design columns. The Type I check used
  three at n = 80 and does not bound the cost of more. Cells of 9–28 donors carry 12 model columns (six for the
  sex × stage cells, six covariates).
- **No power claim** for the adjusted analysis. Adjustment spends degrees of freedom and, under confounding,
  removes confounded signal. The operating studies' power figures were measured without covariates.
- **Structural analogue.** MOTCO's operating characteristics come from InterSIM blocks (498 measured features,
  Gaussian / M-value) at this design (n, cell sizes, three stages, two blocks). They transfer as statements about
  the design, not about RNA / ATAC pseudobulk distributions or about 10,000 log2-CPM features.
- **Preprocessing choices.** chrX/Y/M removal, the CPM ≥ 1 in ≥ k filter, TMM, prior count 2 and the top 5,000 by
  variance per block are fixed choices (section 4). Their sensitivity analysis is a follow-up.
- **Monte Carlo error and `n_jobs`.** p-values from 999 permutations vary by about 0.015 (SD) between seeds near
  p = 0.2–0.35, and they depend on `n_jobs` through MOTCO's per-worker seeding (section 5).
- **No attribution.** MOTCO refuses orientation attribution with covariates, and nothing here calls for it.

## 8. Next steps

- Forced-rank sensitivity (`forced_components` 2, 3, …): the first priority given rank 1. It is MOTCO's documented
  dimensionality sensitivity analysis, and its results are to be read as sensitivity, not as a replacement for the
  CV-selected analysis.
- MOTCO side: decide whether CV may select rank 1, where `angle` is degenerate, given that the forced path rejects
  it; and whether parallel RRPP should reproduce the serial draws, so that `n_jobs` stops changing p-values.
- Preprocessing sensitivity exports: other N, keeping chrX, residualized features.
- Pathway/module enrichment. Under the roadmap it follows a significant orientation result, so without one it is
  exploratory.

## 9. Reproduction

```
uv sync --locked                                   # MOTCO 0.8.0 from the lock
uv run sead-sexdiff analyze --cohort all --force   # writes results/motco_analysis/<cohort>/
# one cohort through MOTCO's own CLI: the manifest's cli_equivalent, OUT_DIR replaced
python -c "import json; print(json.load(open('results/motco_analysis/cohort_30/analysis_manifest.json'))['cli_equivalent']['command'])"
```

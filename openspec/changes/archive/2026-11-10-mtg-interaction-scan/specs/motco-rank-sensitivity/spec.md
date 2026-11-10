## MODIFIED Requirements

### Requirement: Every setting other than the rank equals the analysis
Each forced-rank run SHALL use the same settings as the CV-selected analysis. These are the blocks and units, the stage order, the six covariates in the same order, the missing-covariate policy, the orientation estimator, the residual correction, the split test, the interaction scan, the permutation count and the seed. Every rank then draws the same permutation indices, and differences between ranks are measurement differences. Cross-validation SHALL NOT run.

#### Scenario: Record carries the settings
- **WHEN** the rank-4 run of `cohort_30` finishes
- **THEN** its `result.json` records `forced_components = 4`, `component_selection = forced`, `residual_correction = df` with its κ, `span_split_test = true` with a `span` dimension of 4, `interaction_scan = true`, `permutations = 999`, the analysis seed, `orientation_estimator = net_displacement`, the six covariate names in order and no cross-validation AUROC

#### Scenario: CV-selected rank inside the ladder
- **WHEN** the reference analysis's CV-selected rank is one of the ladder's ranks
- **THEN** the forced run at that rank gives the same statistics, p-values and split parts as the reference analysis, and any difference is recorded as a finding

#### Scenario: Scan equal at every rank
- **WHEN** any forced-rank run finishes
- **THEN** its `interaction_scan.csv` is byte-identical to the reference analysis's and its `interaction_scan` block equals the reference's; otherwise the sweep stops naming the rank and the first differing entry, and the cross-rank table is not written

### Requirement: Per-rank outputs
Each forced-rank run SHALL be written to `results/motco_rank_sensitivity/<cohort>/rank_<kk>/` (`kk` zero-padded to two digits) holding the same files as a CV-selected analysis directory: MOTCO's artifact (`result.json`, `latent_scores.csv`, `pls_loadings.csv`, `stage_means.csv`, `interaction_scan.csv`), `analysis_metadata.csv`, `analysis_manifest.json` and `summary.md`. The sweep SHALL refuse to overwrite an existing cohort directory unless forced.

#### Scenario: Rerun without force
- **WHEN** `results/motco_rank_sensitivity/cohort_30/` exists and the sweep runs again without `--force`
- **THEN** it stops with an error naming the directory and writes nothing

#### Scenario: Forced rank beside angle
- **WHEN** `rank_06/summary.md` is written
- **THEN** its `angle` row names forced rank 6 and says the rank was forced, and its latent block says cross-validation was skipped

### Requirement: Cross-rank table columns
Each row SHALL give, per statistic (`delta`, `angle`, `shape`), the observed value, the RRPP p-value, whether p ≤ α and the null mean, sd and q95. It SHALL also give the split columns of spec `motco-span-split-test`, the scan columns of spec `motco-interaction-scan`, the reference row's CV AUROC and the row's directory. The CSV SHALL hold only values that reproduce across reruns; each run's runtime goes to the Markdown and the manifest.

#### Scenario: Table values equal the records
- **WHEN** a cell of `rank_sensitivity.csv` is compared with the `result.json` of its row's directory
- **THEN** they are equal

#### Scenario: Floating-point cells round-trip
- **WHEN** a `delta` null sd whose shortest representation needs 17 significant digits is written to `rank_sensitivity.csv`
- **THEN** the cell is written in its shortest round-trip representation and parses back to the record's value exactly

#### Scenario: Runtime and directories
- **WHEN** the table of `cohort_30` is written
- **THEN** it has no runtime column, its forced rows name their directory as `rank_<kk>` relative to the cohort directory and its reference row names the reference analysis directory, while `rank_sensitivity.md` and the manifest carry each run's runtime

#### Scenario: Split share across ranks
- **WHEN** the table of `cohort_30` is written
- **THEN** each row carries its in-span share beside its rank (the `span` dimension), so the share's dependence on the rank can be read from the table, and no row reads the share as a test

#### Scenario: Scan columns equal across rows
- **WHEN** the table of `cohort_30` is written
- **THEN** every row's scan cells equal the reference row's, and `rank_sensitivity.md` states that the scan does not depend on the latent rank

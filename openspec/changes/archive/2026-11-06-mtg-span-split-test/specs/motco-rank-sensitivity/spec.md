# Spec Delta

## MODIFIED Requirements

### Requirement: Reference analysis must match the sweep's inputs and settings
Before any rank runs, the sweep SHALL check the reference analysis manifest against the sweep. Each of the five export files' sha256 SHALL equal the file now on disk. The MOTCO version, seed, permutation count, latent-rank floor, residual correction and split test SHALL equal the sweep's. On any difference it SHALL stop with an error naming the file or setting.

#### Scenario: Export regenerated after the analysis
- **WHEN** `rna.csv` was rewritten after `analyze` ran and its sha256 differs from the one in `analysis_manifest.json`
- **THEN** the sweep stops with an error naming `rna.csv` and the `analyze` command, and writes nothing

#### Scenario: Matching reference
- **WHEN** every hash and setting matches
- **THEN** the reference row of the sweep carries the analysis record's statistics, p-values, null summaries, split parts and CV-selected rank unchanged

#### Scenario: Unfloored reference
- **WHEN** the reference analysis was made with the rank floor 1
- **THEN** the sweep stops with an error naming the reference `result.json`, both floors and the `analyze` command, and writes nothing

#### Scenario: Uncorrected reference
- **WHEN** the reference analysis was made without the residual correction
- **THEN** the sweep stops with an error naming the reference `result.json`, the recorded setting, `df` and the `analyze` command, and writes nothing

#### Scenario: Reference without the split test
- **WHEN** the reference analysis carries no `span_split` block
- **THEN** the sweep stops with an error naming the reference `result.json` and the `analyze` command, and writes nothing

### Requirement: Every setting other than the rank equals the analysis
Each forced-rank run SHALL use the same settings as the CV-selected analysis. These are the blocks and units, the stage order, the six covariates in the same order, the missing-covariate policy, the orientation estimator, the residual correction, the split test, the permutation count and the seed. Every rank then draws the same permutation indices, and differences between ranks are measurement differences. Cross-validation SHALL NOT run.

#### Scenario: Record carries the settings
- **WHEN** the rank-4 run of `cohort_30` finishes
- **THEN** its `result.json` records `forced_components = 4`, `component_selection = forced`, `residual_correction = df` with its κ, `span_split_test = true` with a `span` dimension of 4, `permutations = 999`, the analysis seed, `orientation_estimator = net_displacement`, the six covariate names in order and no cross-validation AUROC

#### Scenario: CV-selected rank inside the ladder
- **WHEN** the reference analysis's CV-selected rank is one of the ladder's ranks
- **THEN** the forced run at that rank gives the same statistics, p-values and split parts as the reference analysis, and any difference is recorded as a finding

### Requirement: Cross-rank table columns
Each row SHALL give, per statistic (`delta`, `angle`, `shape`), the observed value, the RRPP p-value, whether p ≤ α and the null mean, sd and q95. It SHALL also give the split columns of spec `motco-span-split-test`, the reference row's CV AUROC and the row's directory. The CSV SHALL hold only values that reproduce across reruns; each run's runtime goes to the Markdown and the manifest.

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

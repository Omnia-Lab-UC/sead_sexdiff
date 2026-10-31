## MODIFIED Requirements

### Requirement: Reference analysis must match the sweep's inputs and settings
Before any rank runs, the sweep SHALL check that the reference analysis manifest records the same sha256 for each of the five export files as the files now on disk, and the same MOTCO version, seed, permutation count, latent-rank floor and residual correction as the sweep uses. On any difference it SHALL stop with an error naming the file or setting.

#### Scenario: Export regenerated after the analysis
- **WHEN** `rna.csv` was rewritten after `analyze` ran and its sha256 differs from the one in `analysis_manifest.json`
- **THEN** the sweep stops with an error naming `rna.csv` and the `analyze` command, and writes nothing

#### Scenario: Matching reference
- **WHEN** every hash and setting matches
- **THEN** the reference row of the sweep carries the analysis record's statistics, p-values, null summaries and CV-selected rank unchanged

#### Scenario: Unfloored reference
- **WHEN** the reference analysis was made with the rank floor 1
- **THEN** the sweep stops with an error naming the reference `result.json`, both floors and the `analyze` command, and writes nothing

#### Scenario: Uncorrected reference
- **WHEN** the reference analysis was made without the residual correction
- **THEN** the sweep stops with an error naming the reference `result.json`, the recorded setting, `df` and the `analyze` command, and writes nothing

### Requirement: Every setting other than the rank equals the analysis
Each forced-rank run SHALL use the same blocks and units, stage order, six covariates in the same order, missing-covariate policy, orientation estimator, residual correction, permutation count and seed as the CV-selected analysis, so that every rank draws the same permutation indices and differences between ranks are measurement differences. Cross-validation SHALL NOT run.

#### Scenario: Record carries the settings
- **WHEN** the rank-4 run of `cohort_30` finishes
- **THEN** its `result.json` records `forced_components = 4`, `component_selection = forced`, `residual_correction = df` with its κ, `permutations = 999`, the analysis seed, `orientation_estimator = net_displacement`, the six covariate names in order and no cross-validation AUROC

#### Scenario: CV-selected rank inside the ladder
- **WHEN** the reference analysis's CV-selected rank is one of the ladder's ranks
- **THEN** the forced run at that rank gives the same statistics and p-values as the reference analysis, and any difference is recorded as a finding

### Requirement: Sweep manifest records provenance
`rank_sensitivity_manifest.json` SHALL record the cohort, the ladder, the seed, the permutation count, the rank floor, the residual correction, α, `n_jobs`, the reference analysis directory with the sha256 of its `result.json` and `analysis_manifest.json`, the sha256 of each export input, each per-rank directory, the MOTCO, Python and library versions, the git revision and findings, including that the CV-selected rank is not in the ladder when it is not.

#### Scenario: Manifest after a cohort sweep
- **WHEN** `cohort_30` finishes with a CV-selected rank of 2
- **THEN** its manifest is valid JSON listing the six ranks, `permutations = 999`, the seed, the rank floor 2, the residual correction `df`, the MOTCO version, hashes for the five export files and the two reference files, and the comparison of the forced rank-2 run with the reference analysis

#### Scenario: CV-selected rank outside the ladder
- **WHEN** a cohort's CV-selected rank is 5
- **THEN** the manifest carries a finding that rank 5 is not in the ladder and the reference row stands in for it

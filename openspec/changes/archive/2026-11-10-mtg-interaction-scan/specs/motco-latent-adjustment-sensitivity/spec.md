## MODIFIED Requirements

### Requirement: Predeclared arms
The step SHALL run exactly two arms per cohort, each once:
- `extended_covariates`, with the covariate-blind basis (no latent adjustment), the split test and the interaction scan;
- `latent_adjusted`, with MOTCO's covariate latent adjustment conditioned on `stage`, with the interaction scan and without the split test, which MOTCO does not allow with a latent adjustment.

The arms SHALL be a fixed registry, not a runtime option. `group_stage` conditioning SHALL NOT be run. Apart from the covariate set, the latent-adjustment setting and the split test it excludes, both arms SHALL use the primary analysis's settings: the same export, blocks and units, group and stage columns, stage order, orientation estimator, CV rank selection, residual correction, interaction scan, seed and permutation count.

#### Scenario: Arms complete
- **WHEN** the step finishes for `cohort_30`
- **THEN** two arm records exist, both with `component_selection = cv`, `residual_correction = df` with a reduced rank of 12 (eight covariates), the primary's seed, 999 permutations and an `interaction_scan` block with `df_residual` n − 14:
  - `extended_covariates`: no `latent_adjustment` block, and a `span_split` block;
  - `latent_adjusted`: a `latent_adjustment` block whose `conditioning` is `stage`, and no `span_split` block, its manifest recording the split test as not available

#### Scenario: Latent adjustment changes only the basis
- **WHEN** both arms of a cohort have finished
- **THEN** the two records carry the same covariate block, and any difference between their statistics arises through the latent basis and the CV-selected rank only, as MOTCO's latent-adjustment contract guarantees

#### Scenario: Scan equal across the arms
- **WHEN** both arms of a cohort have finished
- **THEN** their `interaction_scan.csv` files are byte-identical and their `interaction_scan` blocks are equal; otherwise the step stops naming the cohort and the first differing entry, and the comparison directory is not written

### Requirement: Comparison table generated from the records
Per cohort the step SHALL write `comparison/<cohort>/latent_sensitivity.csv`. It holds the reference row and then one row per arm, each with:
- arm name, covariate names, conditioning (or none);
- selected rank, CV mean AUROC;
- for `delta`, `angle` and `shape`: observed value, p-value, rejection under the primary's α and rejection rule, and null mean, SD and q95;
- the split columns of spec `motco-span-split-test`, empty for `latent_adjusted`;
- the scan columns of spec `motco-interaction-scan`, filled for every row;
- for `latent_adjusted`: the latent-adjustment block's `residual_df`, `cv_min_training_df` and `feature_r2` summary.

The step SHALL also write its Markdown rendering, `latent_sensitivity.md`, and `latent_sensitivity_manifest.json`. The manifest records:
- the covariate set, the arms, seed, permutations and `n_jobs`;
- every input path and sha256, the reference analysis and the arm directories;
- the mediator flags, software, git revision and findings.

Every value SHALL be read from the records, not computed again.

#### Scenario: Values traceable
- **WHEN** `latent_sensitivity.csv` is written
- **THEN** every statistic, p-value, split cell, scan cell and rank in it equals the corresponding value in the reference or arm `result.json`

#### Scenario: Split cells of the latent-adjusted arm
- **WHEN** `latent_sensitivity.csv` is written
- **THEN** the `latent_adjusted` row's split cells are empty and `latent_sensitivity.md` states that the split test is not available with a latent adjustment

#### Scenario: Scan status differs from the primary
- **WHEN** the `extended_covariates` arm's scan is rejected in a cohort and the primary's is not
- **THEN** the comparison manifest records a finding naming the cohort and both p-values as covariate-set dependence of the scan, and `latent_sensitivity.md` states that the `latent_adjusted` row's scan equals the `extended_covariates` row's by construction

## MODIFIED Requirements

### Requirement: Analysis through MOTCO's validated entry point
The analysis SHALL run through MOTCO's trajectory-analysis entry point (`motco analyze` / `analyze_trajectories`) of a pinned MOTCO release ≥ 0.13.0, the first release with the residual-corrected RRPP null. It SHALL NOT assemble the result from `motco plsr`, `motco de` or the statistics primitives.

#### Scenario: Older MOTCO installed
- **WHEN** the installed MOTCO is older than 0.13.0
- **THEN** the analysis stops with an error naming the installed and required versions before reading any input

#### Scenario: Record names the entry point's release
- **WHEN** a cohort analysis finishes
- **THEN** `result.json` records the MOTCO version under `software`

### Requirement: Fixed measurement settings
The analysis SHALL declare blocks `rna` and `atac` in identity units, group column `sex`, stage column `stage` with the explicit order `0`, `1`, `2`, the net-displacement orientation estimator, cross-validated rank selection at the MOTCO defaults with the latent-rank floor `ANALYSIS_RANK_FLOOR = 2` passed explicitly, the residual correction `ANALYSIS_RESIDUAL_CORRECTION = "df"` passed explicitly, 999 permutations and a fixed recorded seed. The floor SHALL equal MOTCO's default; if the installed MOTCO's default differs, the analysis SHALL stop naming both values before reading any input. The residual correction is MOTCO's validated opt-in (`adopt_df`) and not its default; the analysis SHALL check that the record states the correction it applied, with κ, and SHALL write nothing otherwise.

#### Scenario: Record carries the measurement settings
- **WHEN** a cohort analysis finishes
- **THEN** `result.json` records `orientation_estimator = net_displacement`, `component_selection = cv`, `min_components = 2`, `residual_correction = df` under `params` and the applied correction (kind `df`, κ, n, reduced rank) under `measurement`, a selected rank of at least 2, `permutations = 999`, the seed, the ordered stage levels `0, 1, 2` and the group levels, and `analysis_manifest.json` and its `cli_equivalent` carry the floor and `--residual-correction df`

#### Scenario: MOTCO default diverges
- **WHEN** the installed MOTCO's default `min_components` is not 2
- **THEN** the analysis stops with an error naming MOTCO's default and `ANALYSIS_RANK_FLOOR`, and writes nothing

#### Scenario: Correction absent from the record
- **WHEN** MOTCO returns a record whose `measurement` carries no residual correction, or one of another kind than `df`
- **THEN** the analysis stops with an error naming the recorded and the expected correction, and writes nothing

#### Scenario: κ follows the design
- **WHEN** `cohort_30` (n = 80) is analysed with the six analysis covariates
- **THEN** the record's reduced rank is 10 and its κ is √(80 / 70) ≈ 1.069

### Requirement: Summary table generated from the record
`summary.md` SHALL be generated from `result.json`, not typed by hand, and SHALL show per statistic (`delta`, `angle`, `shape`) the observed value, the RRPP p-value and the null mean, sd and 95 % quantile; the CV-selected rank, the rank floor and `cv_mean_auroc`; the residual correction the null was drawn under (kind, κ, n, reduced rank); the sex × stage cell sizes; and the covariate block (design columns, reference row, missing policy, dropped donors).

#### Scenario: Rank beside angle
- **WHEN** `summary.md` is written
- **THEN** the `angle` row names the CV-selected rank from `result.json` and the summary states the rank floor it was selected under

#### Scenario: Correction stated
- **WHEN** `summary.md` is written for a corrected analysis
- **THEN** it states that the RRPP null used the `df` residual correction, with κ, n and the reduced rank from the record

#### Scenario: Record without a correction
- **WHEN** `summary.md` is rendered from a record that carries no residual correction (a superseded analysis)
- **THEN** the output is the one rendered before this change, with no correction line

### Requirement: Report applies the binding interpretation
The report SHALL state that a `delta` rejection is evidence of a difference and not specifically of size; that `angle` and `shape` non-rejections are uninformative at this n; that `angle` reads only the endpoint stages; that the adjustment is additive-linear on a covariate-blind basis; that the operating characteristics are a structural analogue (InterSIM blocks, 498 features); that the rank floor defines `angle` and `shape` in at least two dimensions without making them powerful, giving each analysis's `config_spectrum` eigengap and null summary beside `angle`; and that under the `df` correction the p-values are size-controlled at the configuration MOTCO measured (six covariates, both `cohort_30` tables), with `delta` conservative there so that a `delta` non-rejection stays weak evidence, and inferred, not measured, elsewhere.

#### Scenario: Non-rejection wording
- **WHEN** `angle` is not rejected at α = 0.05 in a cohort
- **THEN** the report states that the result is uninformative about orientation, not that the directions agree

#### Scenario: Delta rejection wording
- **WHEN** `delta` is rejected at α = 0.05 in a cohort
- **THEN** the report states that the sexes' trajectories differ, without claiming the difference is in size

#### Scenario: Floored angle
- **WHEN** a floored analysis selects rank 2 and does not reject `angle`
- **THEN** the report gives the rank, the eigengap and the `angle` null summary beside the statistic, and states that the floor guarantees two dimensions, not a second axis with signal

#### Scenario: Corrected delta non-rejection
- **WHEN** the corrected primary does not reject `delta` in a cohort
- **THEN** the report states that MOTCO measured the corrected `delta` as conservative at this configuration (0.029 / 0.031 pooled), so the non-rejection is weak evidence against a size difference

#### Scenario: Sensitivity cohort
- **WHEN** the report presents the corrected `cohort_50` analysis
- **THEN** it states that the corrected null's size there is inferred from MOTCO's mechanism (n and the reduced rank), not measured

## ADDED Requirements

### Requirement: Residual correction checked by every reader of the primary
Every step that reads a primary or arm analysis as its reference (rank, preprocessing, latent-adjustment and donor-exclusion sensitivity, both stage-DE steps) SHALL check that the reference record was made with `ANALYSIS_RESIDUAL_CORRECTION`, beside its MOTCO-version, seed, permutation and rank-floor checks, and SHALL stop naming the file and both settings otherwise.

#### Scenario: Uncorrected reference
- **WHEN** `stage-de` reads a primary `result.json` whose parameters carry no `residual_correction`
- **THEN** it stops with an error naming the file, the recorded setting (none), `df` and the `analyze` command, and writes nothing

#### Scenario: Corrected reference
- **WHEN** `rank-sensitivity` reads a primary made under MOTCO 0.13.0 at floor 2 with `residual_correction = df`
- **THEN** the reference check passes and the sweep runs

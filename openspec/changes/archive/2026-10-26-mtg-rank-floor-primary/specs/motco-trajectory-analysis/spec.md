# Spec Delta

## MODIFIED Requirements

### Requirement: Analysis through MOTCO's validated entry point
The analysis SHALL run through MOTCO's trajectory-analysis entry point (`motco analyze` / `analyze_trajectories`) of a pinned MOTCO release ≥ 0.10.0, the first release with the latent-rank floor. It SHALL NOT assemble the result from `motco plsr`, `motco de` or the statistics primitives.

#### Scenario: Older MOTCO installed
- **WHEN** the installed MOTCO is older than 0.10.0
- **THEN** the analysis stops with an error naming the installed and required versions before reading any input

#### Scenario: Record names the entry point's release
- **WHEN** a cohort analysis finishes
- **THEN** `result.json` records the MOTCO version under `software`

### Requirement: Fixed measurement settings
The analysis SHALL declare blocks `rna` and `atac` in identity units, group column `sex`, stage column `stage` with the explicit order `0`, `1`, `2`, the net-displacement orientation estimator, cross-validated rank selection at the MOTCO defaults with the latent-rank floor `ANALYSIS_RANK_FLOOR = 2` passed explicitly, 999 permutations and a fixed recorded seed. The floor SHALL equal MOTCO's default; if the installed MOTCO's default differs, the analysis SHALL stop naming both values before reading any input.

#### Scenario: Record carries the measurement settings
- **WHEN** a cohort analysis finishes
- **THEN** `result.json` records `orientation_estimator = net_displacement`, `component_selection = cv`, `min_components = 2`, a selected rank of at least 2, `permutations = 999`, the seed, the ordered stage levels `0, 1, 2` and the group levels, and `analysis_manifest.json` and its `cli_equivalent` carry the floor

#### Scenario: MOTCO default diverges
- **WHEN** the installed MOTCO's default `min_components` is not 2
- **THEN** the analysis stops with an error naming MOTCO's default and `ANALYSIS_RANK_FLOOR`, and writes nothing

### Requirement: Summary table generated from the record
`summary.md` SHALL be generated from `result.json`, not typed by hand, and SHALL show per statistic (`delta`, `angle`, `shape`) the observed value, the RRPP p-value and the null mean, sd and 95 % quantile; the CV-selected rank, the rank floor and `cv_mean_auroc`; the sex × stage cell sizes; and the covariate block (design columns, reference row, missing policy, dropped donors).

#### Scenario: Rank beside angle
- **WHEN** `summary.md` is written
- **THEN** the `angle` row names the CV-selected rank from `result.json` and the summary states the rank floor it was selected under

### Requirement: Report applies the binding interpretation
The report SHALL state that a `delta` rejection is evidence of a difference and not specifically of size; that `angle` and `shape` non-rejections are uninformative at this n; that `angle` reads only the endpoint stages; that the adjustment is additive-linear on a covariate-blind basis; that the operating characteristics are a structural analogue (InterSIM blocks, 498 features); and that the rank floor defines `angle` and `shape` in at least two dimensions without making them powerful, giving each analysis's `config_spectrum` eigengap and null summary beside `angle`.

#### Scenario: Non-rejection wording
- **WHEN** `angle` is not rejected at α = 0.05 in a cohort
- **THEN** the report states that the result is uninformative about orientation, not that the directions agree

#### Scenario: Delta rejection wording
- **WHEN** `delta` is rejected at α = 0.05 in a cohort
- **THEN** the report states that the sexes' trajectories differ, without claiming the difference is in size

#### Scenario: Floored angle
- **WHEN** a floored analysis selects rank 2 and does not reject `angle`
- **THEN** the report gives the rank, the eigengap and the `angle` null summary beside the statistic, and states that the floor guarantees two dimensions, not a second axis with signal

## ADDED Requirements

### Requirement: Rank floor checked by every reader of the primary
Every step that reads a primary or arm analysis as its reference (rank, preprocessing, latent-adjustment and donor-exclusion sensitivity, both stage-DE steps) SHALL check that the reference record was made at `ANALYSIS_RANK_FLOOR`, beside its existing MOTCO-version, seed and permutation checks, and SHALL stop naming the file and both floors otherwise.

#### Scenario: Unfloored reference
- **WHEN** `latent-sensitivity` reads a primary `result.json` whose parameters carry no `min_components` (floor 1)
- **THEN** it stops with an error naming the file, floor 1 and floor 2, and the `analyze` command, and writes nothing

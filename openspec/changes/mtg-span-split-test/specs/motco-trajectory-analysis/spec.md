# Spec Delta

## MODIFIED Requirements

### Requirement: Analysis through MOTCO's validated entry point
The analysis SHALL run through MOTCO's trajectory-analysis entry point (`motco analyze` / `analyze_trajectories`) of a pinned MOTCO release ≥ 0.18.0, the first release with the split interaction test (and ≥ 0.13.0's residual-corrected RRPP null). It SHALL NOT assemble the result from `motco plsr`, `motco de` or the statistics primitives, including `motco.stats.interaction`.

#### Scenario: Older MOTCO installed
- **WHEN** the installed MOTCO is older than 0.18.0
- **THEN** the analysis stops with an error naming the installed and required versions before reading any input

#### Scenario: Record names the entry point's release
- **WHEN** a cohort analysis finishes
- **THEN** `result.json` records the MOTCO version under `software`

### Requirement: Fixed measurement settings
The analysis SHALL declare:
- blocks `rna` and `atac` in identity units;
- group column `sex`, and stage column `stage` with the explicit order `0`, `1`, `2`;
- the net-displacement orientation estimator;
- cross-validated rank selection at the MOTCO defaults, with the latent-rank floor `ANALYSIS_RANK_FLOOR = 2` passed explicitly;
- the residual correction `ANALYSIS_RESIDUAL_CORRECTION = "df"` and the split test `ANALYSIS_SPAN_SPLIT_TEST = True`, both passed explicitly;
- 999 permutations and a fixed recorded seed.

The floor SHALL equal MOTCO's default; if the installed MOTCO's default differs, the analysis SHALL stop naming both values before reading any input. The correction (MOTCO's validated opt-in `adopt_df`) and the split test (`split_validated`) are not MOTCO defaults. The analysis SHALL check that the record states the correction with κ and carries the split block, and SHALL write nothing otherwise.

#### Scenario: Record carries the measurement settings
- **WHEN** a cohort analysis finishes
- **THEN** `result.json` records under `params` `orientation_estimator = net_displacement`, `component_selection = cv`, `min_components = 2`, `residual_correction = df` and `span_split_test = true`. It records the applied correction (kind `df`, κ, n, reduced rank) under `measurement`, a `span_split` block, a selected rank of at least 2, `permutations = 999`, the seed, the ordered stage levels `0, 1, 2` and the group levels. `analysis_manifest.json` and its `cli_equivalent` carry the floor, `--residual-correction df` and `--span-split-test`

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
`summary.md` SHALL be generated from `result.json`, not typed by hand. It SHALL show:
- per statistic (`delta`, `angle`, `shape`), the observed value, the RRPP p-value and the null mean, sd and 95 % quantile;
- the CV-selected rank, the rank floor and `cv_mean_auroc`;
- the residual correction the null was drawn under (kind, κ, n, reduced rank);
- when the record carries a `span_split` block, the split table of spec `motco-span-split-test`;
- the sex × stage cell sizes;
- the covariate block (design columns, reference row, missing policy, dropped donors).

#### Scenario: Rank beside angle
- **WHEN** `summary.md` is written
- **THEN** the `angle` row names the CV-selected rank from `result.json` and the summary states the rank floor it was selected under

#### Scenario: Correction stated
- **WHEN** `summary.md` is written for a corrected analysis
- **THEN** it states that the RRPP null used the `df` residual correction, with κ, n and the reduced rank from the record

#### Scenario: Split table shown
- **WHEN** `summary.md` is written for an analysis whose record carries a `span_split` block
- **THEN** it shows the three parts with dimension, pseudo-F, pseudo-F p-value and role, and the in-span share labelled as not read

#### Scenario: Record without a correction
- **WHEN** `summary.md` is rendered from a record that carries no residual correction (a superseded analysis)
- **THEN** the output is the one rendered before the correction was introduced, with no correction line and no split table

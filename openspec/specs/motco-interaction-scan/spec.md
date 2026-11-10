# motco-interaction-scan Specification

## Purpose
The per-feature interaction scan that every MOTCO analysis of the case study carries: how it is set, checked, recorded and tabulated, and the predeclared multiplicity policy, power statement and reading under which its global p-value and hits are reported.

## Requirements

### Requirement: Scan in every analysis
Every MOTCO analysis the package runs SHALL request MOTCO's per-feature interaction scan through the registered constant `ANALYSIS_INTERACTION_SCAN = True`, passed explicitly, never as a flag. This covers the primaries, the forced ranks, the preprocessing variants, both latent arms and the exclusion variants. The scan is available with a latent adjustment, so no analysis is exempt.

#### Scenario: Primary and sensitivity analyses
- **WHEN** `analyze`, `rank-sensitivity`, `preprocessing-sensitivity`, `latent-sensitivity` or `exclusion-sensitivity` runs an analysis
- **THEN** its `result.json` holds an `interaction_scan` block and `params.interaction_scan = true`, and its directory holds `interaction_scan.csv`

#### Scenario: Latent-adjusted arm
- **WHEN** the `latent_adjusted` arm is run
- **THEN** its record carries the scan block although it carries no split block

### Requirement: Scan record checked before anything is written
After each analysis and before anything is written, the `interaction_scan` block SHALL be present exactly when requested and SHALL agree with the analysis. That covers the features per block, the design's degrees of freedom, a global p-value, the hit level `ANALYSIS_ALPHA`, the analysis's permutations, seed and correction, and `draws_match_rrpp = true` at `n_jobs = 1`. Otherwise the analysis SHALL stop, naming the cohort and the entry, with nothing written.

#### Scenario: Block missing
- **WHEN** MOTCO returns a record without `interaction_scan` for an analysis that requested the scan
- **THEN** the analysis stops naming the cohort and the missing block, and no directory is written

#### Scenario: Degrees of freedom follow the design
- **WHEN** the astrocyte `cohort_30` primary (80 donors, six covariates, 5,000 RNA and 5,000 ATAC features) is analysed
- **THEN** the block records 10,000 features (5,000 per block), `df_interaction` 2 (stages − 1), `df_residual` 68 (donors − 6 cell means − 6 covariates), 999 permutations, the seed 20261019 and the `df` correction

### Requirement: Scan table checked before anything is written
The scan table SHALL hold one row per analysed feature in the record's block order. Its columns are block, feature, interaction SS, residual SS, F, adjusted p-value, degeneracy, and the difference at each stage `0`, `1`, `2`. Its hits SHALL be exactly the rows with adjusted p ≤ `ANALYSIS_ALPHA` and SHALL equal the record's hit list and per-block counts. The global p-value SHALL be ≤ α exactly when there is at least one hit. Otherwise the analysis SHALL stop with nothing written.

#### Scenario: Rejection without a hit
- **WHEN** a record's global p-value is ≤ 0.05 but no row of its table has adjusted p ≤ 0.05
- **THEN** the analysis stops naming the cohort and both values, and no directory is written

#### Scenario: Hit-level constant moved
- **WHEN** the installed MOTCO's hit level (`INTERACTION_SCAN_ALPHA`) is not `ANALYSIS_ALPHA`
- **THEN** the analysis stops before reading any input, naming both values

### Requirement: Scan recorded in the manifest, the CLI equivalent and the summary
The analysis manifest SHALL record the scan under `measurement.interaction_scan`: requested, the confirmatory statistic, the hit statistic and level, and the values not read. `cli_equivalent` SHALL carry `--interaction-scan`. `summary.md` SHALL show, generated from the record: features and degrees of freedom, the maximum F with its feature and block, the global p-value and status, hits per block, the null-maximum shares labelled descriptive, and one row per hit.

#### Scenario: Summary of a primary without hits
- **WHEN** `summary.md` is written for a primary whose scan has no hit
- **THEN** it has a scan table with the global p-value marked confirmatory, zero hits per block, the null-maximum shares marked descriptive and the line "no hit", and no per-feature F beyond the maximum

#### Scenario: Record without the scan
- **WHEN** `summary.md` is rendered from a record without an `interaction_scan` block (a superseded analysis)
- **THEN** it has no scan table and is otherwise unchanged

### Requirement: Scan checked by every reader of the primary
Every step that reads a primary analysis as its reference SHALL check that the reference record carries `params.interaction_scan = true` and an `interaction_scan` block, beside its version, seed, permutation, rank-floor, correction and split checks. These steps are rank, preprocessing, latent-adjustment and donor-exclusion sensitivity and both stage-DE steps. Otherwise the step SHALL stop, naming the file and the `analyze` command.

#### Scenario: Reference without the scan
- **WHEN** `stage-de` reads a primary `result.json` made under MOTCO 0.18.0
- **THEN** it stops naming the file and the `analyze` command, and writes nothing

### Requirement: Scan identical where only the latent space differs
The scan reads the features and trajectory models only. So every forced-rank analysis's scan block and table SHALL equal its CV-selected reference's, and the `latent_adjusted` arm's SHALL equal the `extended_covariates` arm's. Only the provenance entries may differ. Otherwise the step SHALL stop, naming the analysis and the first differing entry, with the comparison outputs unwritten. These equalities are checks, never sensitivity readings.

#### Scenario: Forced rank
- **WHEN** `rank-sensitivity` finishes forced rank 9 of astrocyte `cohort_30`
- **THEN** its `interaction_scan.csv` is byte-identical to the primary's, and its block equals the primary's

#### Scenario: Latent arms differ
- **WHEN** the `latent_adjusted` arm's global p-value differs from the `extended_covariates` arm's in a cohort
- **THEN** `latent-sensitivity` stops naming the cohort, both arms and the value, and the comparison directory is not written

### Requirement: Scan columns in every comparison table
Every comparison table generated from analysis records SHALL give, per row, after the split columns:
- the maximum F and its block;
- the global p-value and status (p ≤ α);
- hits per block.

The tables are `rank_sensitivity.csv`, `sensitivity.csv`, `latent_sensitivity.csv` and `exclusion_sensitivity.csv`. These cells SHALL be empty for a row whose analysis carries no scan. A scan status that differs from the reference row's SHALL be written as a dependence finding.

#### Scenario: Table values equal the records
- **WHEN** a scan cell of `exclusion_sensitivity.csv` is compared with the `result.json` of its row's directory
- **THEN** they are equal

#### Scenario: Variant status differs
- **WHEN** the `residualized` variant's scan is rejected in `cohort_30` and the primary's is not
- **THEN** the comparison manifest records a finding naming the variant and both p-values as preprocessing dependence, not a result of the variant

### Requirement: Global p-value is the scan's confirmatory reading
The scan's global p-value SHALL be its only confirmatory reading. It is read at `ANALYSIS_ALPHA` under `ANALYSIS_REJECTION_RULE` (p ≤ 0.05), unadjusted. No adjustment SHALL be made across `delta`, `angle`, `shape`, the complement and the scan, or across cohorts and cell types. The report SHALL state that each primary carries five tests read at α. This policy is predeclared and SHALL NOT be changed from results.

#### Scenario: Scan rejected in a primary
- **WHEN** the `cohort_30` primary of a cell type has a global scan p-value ≤ 0.05
- **THEN** the report gives it as that cell type's scan result under the reading requirement, whatever `delta`, `angle`, `shape` and the complement read

### Requirement: Hits are the scan's familywise detections
A hit SHALL be a feature whose step-down max-T adjusted p-value is ≤ `ANALYSIS_ALPHA`. The familywise error rate over all scanned features is then controlled at α within the scan. Hits SHALL be reported with block, feature, F, adjusted p-value and the per-stage differences. A `cohort_30` primary's hits SHALL be part of its result. Any other analysis's hits SHALL be reported as sensitivity.

#### Scenario: Hits in a primary
- **WHEN** the microglia `cohort_30` primary's scan has three RNA hits and one ATAC hit
- **THEN** the report lists the four features with their adjusted p-values and differences at stages `0`, `1` and `2`, the ATAC hit named by its peak coordinates

### Requirement: Descriptive scan values not read
The per-feature F and interaction SS, the maximum F, the per-stage differences and the per-block null-maximum shares SHALL be reported where tabulated and SHALL NOT be read as results. No status SHALL be derived from them, and no per-stage difference SHALL be tested.

#### Scenario: Large maximum F without a rejection
- **WHEN** a primary's maximum F is the largest of all analyses but its global p-value is 0.31
- **THEN** the report states the non-rejection under the reading requirement and makes no claim about the feature with the maximum F

### Requirement: Power statement for the astrocyte cohort_30 primary
When the astrocyte `cohort_30` primary's scan has no hit, the report SHALL state: no sparse ramp of the plasmode's form was detected. At this export's design, such a ramp (20 features, 10 per block, in `Male`, stage weights 0 / 0.5 / 1, random signs) is detected in 1.000 of replicates at 3 within-stage SDs and in 0.692 at 2 SDs. This SHALL require the five export files to hash as MOTCO's plasmode settings record. Every other analysis's non-rejection SHALL be uninformative.

#### Scenario: Export unchanged
- **WHEN** the astrocyte `cohort_30` primary has no hit and its five export files hash as `results/interaction-scan-2026-11-09/plasmode/settings.json` records
- **THEN** the report gives the power statement, and states that it excludes neither smaller effects, nor effects spread thinly over many features, nor other shapes

#### Scenario: Export hashes differ
- **WHEN** any of the five files hashes differently from the plasmode's record
- **THEN** the non-rejection is uninformative, and the report records a finding naming the file and both hashes

#### Scenario: Microglia without hits
- **WHEN** the microglia `cohort_30` primary has no hit
- **THEN** the report reads its non-rejection as uninformative, with no power statement

### Requirement: Reading of the scan in the report
The report SHALL read the scan as follows:
- each cell type's `cohort_30` primary gives the result, and `cohort_50` and every other analysis are sensitivity;
- a status differing from the primary's is preprocessing, covariate-set or cohort dependence;
- a rejection names features whose covariate-adjusted sex difference changes across stages, at familywise level α. It is not a driver claim, carries no geometry, starts no attribution and changes no other test's reading.

#### Scenario: Rejection read
- **WHEN** a `cohort_30` primary's scan is rejected
- **THEN** the report names the hits as detections, states that correlated neighbours can be hits too, and gives no trajectory, direction-of-program or driver statement

#### Scenario: Forced rank or latent-adjusted arm
- **WHEN** the report tabulates the forced ranks and the `latent_adjusted` arm
- **THEN** their scan rows are stated as identical to their reference's by construction, and no dependence is read from them

### Requirement: Scan limitations carried
The report SHALL carry the scan's limitations, each with MOTCO's numbers and source:
- size measured on real noise for one export only;
- the null maxima's block imbalance;
- the plain, unmoderated per-feature F;
- the stylized plant behind the power statement;
- the scan's weakness for differences spread thinly over many features.

#### Scenario: Limitations in the report
- **WHEN** the report is committed
- **THEN** its limitations section states, beside the split-test and earlier limitations:
  - real-noise size only for the astrocyte `cohort_30` export at six covariates (global 0.042; familywise 0.042 under the plasmode's partial null), inferred for every other analysis;
  - ATAC giving 0.618 of the null maxima with half the features, so the threshold leans on ATAC;
  - the plain F at n = 80;
  - the 20-feature ramp in one sex as the only plant with measured power;
  - magnitude power 0.01 against `delta`'s 0.21 on the simulated SEA-AD table

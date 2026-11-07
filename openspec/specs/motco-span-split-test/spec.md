# motco-span-split-test Specification

## Purpose

The split interaction test that every MOTCO analysis of the case study carries where MOTCO allows it: how it is set, checked, recorded and tabulated, and the predeclared multiplicity policy and reading under which its three parts are reported.

## Requirements

### Requirement: Split test in every analysis that MOTCO allows
Every MOTCO analysis the package runs SHALL request MOTCO's split interaction test through the registered constant `ANALYSIS_SPAN_SPLIT_TEST = True`, passed explicitly, never as a flag. The one exception is an analysis with a latent adjustment, which MOTCO rejects with the test. That analysis SHALL request it off, and its manifest SHALL record the test as not available, with MOTCO's reason.

#### Scenario: Primary and sensitivity analyses
- **WHEN** `analyze`, `rank-sensitivity`, `preprocessing-sensitivity`, the `extended_covariates` arm of `latent-sensitivity` or `exclusion-sensitivity` runs an analysis
- **THEN** its `result.json` holds a `span_split` block and `params.span_split_test = true`

#### Scenario: Latent-adjusted arm
- **WHEN** the `latent_adjusted` arm is run
- **THEN** its `result.json` has no `span_split` block and no `params.span_split_test`, and its manifest's `measurement.span_split` states that the test is not available with a latent adjustment

### Requirement: Split record checked before anything is written
After each analysis and before anything is written, the record's `span_split` block SHALL be present exactly when the test was requested. The block SHALL hold:
- the parts `total`, `span` and `complement`, each with a pseudo-F p-value;
- the analysis's permutations, seed and residual correction;
- a `span` dimension equal to the selected rank, and a `complement` dimension equal to the `total` dimension minus that rank;
- `draws_match_rrpp = true` when `n_jobs` is 1.

Otherwise the analysis SHALL stop naming the cohort and the differing entry, with nothing written.

#### Scenario: Block missing
- **WHEN** MOTCO returns a record without `span_split` for an analysis that requested the test
- **THEN** the analysis stops naming the cohort and the missing block, and no directory is written

#### Scenario: Dimensions follow the rank
- **WHEN** the astrocyte `cohort_30` primary selects rank 5 on 5,000 RNA and 5,000 ATAC features
- **THEN** the block records `span` dimension 5, `complement` dimension 9,995 and `total` dimension 10,000, with 999 permutations, the seed 20261019 and the `df` correction

### Requirement: Split test recorded in the manifest, the CLI equivalent and the summary
The analysis manifest SHALL record the test under `measurement.span_split`: requested or not available, the confirmatory part, the statistic read, and the parts' roles. `cli_equivalent` SHALL carry `--span-split-test` when the test was requested. `summary.md` SHALL show, generated from the record, per part:
- the dimension;
- the pseudo-F and its p-value;
- its role (confirmatory, concordance, descriptive).

It SHALL also show the in-span share, labelled as not read. The interaction-SS p-values SHALL stay in `result.json` only.

#### Scenario: Summary of a primary
- **WHEN** `summary.md` is written for a primary
- **THEN** it has a split table with the `complement` row marked confirmatory, the `span` row marked as the concordance cell and the `total` row marked descriptive, and no interaction-SS p-value

#### Scenario: Record without the test
- **WHEN** `summary.md` is rendered from a record without a `span_split` block (a superseded analysis or the `latent_adjusted` arm)
- **THEN** it has no split table and is otherwise unchanged

### Requirement: Split test checked by every reader of the primary
Every step that reads a primary analysis as its reference SHALL check that the reference record carries `params.span_split_test = true` and a `span_split` block, beside its MOTCO-version, seed, permutation, rank-floor and residual-correction checks. These steps are rank, preprocessing, latent-adjustment and donor-exclusion sensitivity and both stage-DE steps. Otherwise the step SHALL stop naming the file and the `analyze` command.

#### Scenario: Reference without the test
- **WHEN** `rank-sensitivity` reads a primary `result.json` made under MOTCO 0.18.0 without the split test
- **THEN** it stops naming the file and the `analyze` command, and writes nothing

### Requirement: Split columns in every comparison table
Every comparison table generated from analysis records SHALL give, per row:
- the pseudo-F and its p-value for each part;
- the `complement` status (p ≤ α);
- the `span` status and the span concordance cell;
- the in-span share.

The tables are `rank_sensitivity.csv`, `sensitivity.csv`, `latent_sensitivity.csv` and `exclusion_sensitivity.csv`. These cells SHALL be empty for a row whose analysis carries no split test. No `total` status SHALL be written.

#### Scenario: Table values equal the records
- **WHEN** a split cell of `rank_sensitivity.csv` is compared with the `result.json` of its row's directory
- **THEN** they are equal

#### Scenario: Latent-adjusted row
- **WHEN** `latent_sensitivity.csv` is written
- **THEN** the `latent_adjusted` row's split cells are empty and the `extended_covariates` row's are filled

### Requirement: Complement is the one confirmatory part
The `complement` part's pseudo-F p-value SHALL be the split test's only confirmatory reading. It is read at `ANALYSIS_ALPHA` under `ANALYSIS_REJECTION_RULE` (p ≤ 0.05), unadjusted. No adjustment SHALL be made across `delta`, `angle`, `shape` and the complement, or across cohorts and cell types. The report SHALL state that each primary carries four tests read at α. This policy is predeclared and SHALL NOT be changed from results.

#### Scenario: Complement rejected in a primary
- **WHEN** the `cohort_30` primary of a cell type has a complement pseudo-F p-value ≤ 0.05
- **THEN** the report gives it as that cell type's split-test result under the reading requirement, whatever `total` and `span` read

### Requirement: Span part is a concordance cell
The `span` part's status (pseudo-F p ≤ α) SHALL be compared with whether any of `delta`, `angle` and `shape` is rejected in the same analysis. The two agree or the cell is discordant. A discordance in either direction SHALL be reported as a finding, never as a result, and SHALL change no reading.

#### Scenario: Span discordant
- **WHEN** the span part rejects and none of `delta`, `angle`, `shape` does, or the reverse
- **THEN** the concordance cell reads discordant, the report records it as a finding with the documented power ordering (the span part trails `delta` and `shape`), and no reading changes

### Requirement: Total, share and interaction-SS p-values not read
`total`'s pseudo-F p-value, the in-span share and every interaction-SS p-value SHALL be reported and SHALL NOT be read: no rejection, non-rejection or magnitude of them is interpreted.

#### Scenario: Total rejected alone
- **WHEN** `total` has p ≤ 0.05 and the complement does not
- **THEN** the report gives the value and states that `total` is not read, and no result is claimed

### Requirement: Reading of the split test in the report
The report SHALL read the split test as follows:
- **Results.** Each cell type's `cohort_30` primary gives its result. `cohort_50` and every other analysis are sensitivity: a complement status that differs from the primary's is rank, preprocessing, covariate-set or cohort dependence.
- **A complement rejection** is a sex × stage difference with interaction energy outside the fitted stage-PLS span, not a different molecular program. It may hold in-span leakage and the fitted span's misses of the true stage program, and it carries no geometry, direction, feature list or attribution.
- **A complement non-rejection** is uninformative.
- **The in-span share** is reported beside the rank and is never read against 0.

#### Scenario: Complement non-rejection
- **WHEN** the complement is not rejected in a primary
- **THEN** the report states that the result is uninformative, because no power was shown at SEA-AD's stage signal or for sparse differences (MOTCO's plasmode: 0.02–0.05 for 20 features at 2–3 within-stage SDs), and does not state that the sexes agree outside the span

#### Scenario: Complement rejection
- **WHEN** the complement is rejected in a primary
- **THEN** the report states that the sexes' stage trajectories differ with energy outside the fitted span, names the two other sources such energy can have (in-span leakage, rejected in 0.32 of replicates for a 45° in-span rotation on the SEA-AD table; the fitted span keeping about 0.10 of the population stage program at SEA-AD's signal), and makes no claim about direction, program or features

#### Scenario: Status differs in a sensitivity analysis
- **WHEN** a forced rank rejects the complement and the CV-selected primary does not
- **THEN** the report gives it as rank dependence of the split test, not a result at that rank

### Requirement: Split-test limitations carried
The report SHALL carry these limitations:
- **Size.** On real noise it was measured only for the astrocyte `cohort_30` export at this configuration (complement 0.012, span 0.050, total 0.040 at the plasmode anchor). For `cohort_50`, microglia, the eight-covariate arm, the exclusion variants, the forced ranks and the preprocessing variants it is inferred, not measured.
- **Out-of-span power** was shown only for a dense rotation at a strong stage signal.
- **The null in-span share grows with the fitted rank** (0.236 at rank 5 on this export), and no null share is computed here.
- **The test is not available** on the `latent_adjusted` arm.

#### Scenario: Limitations section
- **WHEN** the report is written
- **THEN** its limitations list each of these with MOTCO's numbers and the report they come from

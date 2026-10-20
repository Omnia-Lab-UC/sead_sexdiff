# Spec Delta

## MODIFIED Requirements

### Requirement: Predeclared variant registry
The variants SHALL be a fixed registry, not a runtime choice of settings. The registry SHALL hold `keep_chrx`, `exclude_families`, `rna_trend_ratio`, `atac_mean_accessibility`, `residualized`, `revised_candidate`, `trend_normalized` and `revised_candidate_trend`, each with a description and the export settings it changes; every setting not listed for a variant SHALL equal the primary export's. The step SHALL run only registered variants and SHALL NOT add, drop or alter a variant from the results.

#### Scenario: keep_chrx settings
- **WHEN** the `keep_chrx` export of a cohort finishes
- **THEN** its manifest records excluded chromosomes `chrY` and `chrM` only, the variant name, and the same `min_cpm`, `k`, `prior_count`, `N` and TMM settings as the primary export

#### Scenario: Unregistered variant
- **WHEN** a variant name not in the registry is requested
- **THEN** the step stops with an error naming the registered variants, and writes nothing

#### Scenario: Registry listed in the help
- **WHEN** `preprocessing-sensitivity --help` is shown
- **THEN** `--variant` lists the eight names and `all`

### Requirement: Variant export is the primary pipeline with one setting changed
The variant export SHALL apply the primary export's steps in the primary order (restrict to the cohort, remove excluded chromosomes and excluded gene sets, CPM ≥ `min_cpm` in ≥ `k` donors, TMM, log2-CPM, trend normalization when declared, residualization when declared, top `N` by the selection rule) with the primary's `k`, `min_cpm`, `prior_count`, `N` and TMM settings, differing only in the settings the registry lists for the variant, computed on the cohort's donors without sex, stage or analysis covariates. It SHALL write the same seven files as the primary export, plus `residualization.csv` when it residualizes and `normalization_trend.csv` when it trend-normalizes.

#### Scenario: chrX features pass through the filters
- **WHEN** `keep_chrx` exports the RNA block of `cohort_30`
- **THEN** no chrX gene in `features_rna.csv` has status `sex_chromosome`, every chrY and chrM gene does, chrX genes carry `abundance`, `variance` or `kept`, the status counts sum to 36,601, and the `kept` count equals the number of columns of `rna.csv`

#### Scenario: Whole pipeline recomputed
- **WHEN** `keep_chrx` exports a cohort
- **THEN** the abundance filter, TMM factors, log2-CPM and variance ranking are computed with chrX features present, and the block holds `N` features (or all that pass, with a finding), not the primary's features plus chrX

#### Scenario: Donor tables unchanged
- **WHEN** any variant exports a cohort
- **THEN** its `metadata.csv` and `covariates.csv` are byte-identical to the primary export's

#### Scenario: Primary export unchanged by the new settings
- **WHEN** `export` runs with its defaults after this change
- **THEN** its six CSVs are byte-identical to the kept primary export and its manifest has no `excluded_gene_sets`, `selection_rule`, `residualization` or `normalization` key beyond the primary's recorded defaults

#### Scenario: Normalization precedes residualization and selection
- **WHEN** `revised_candidate_trend` exports a cohort
- **THEN** the manifest's preprocessing order lists the trend normalization after log2-CPM and before the residualization, and the residualization and the selection are computed on the normalized values

### Requirement: Comparison carries the feature accounting
`sensitivity.md` and `sensitivity_manifest.json` SHALL give, per omic and for both the primary and the variant export, the feature counts after each export step and the numbers of kept features shared with the primary, kept only by the primary and kept only by the variant; for `keep_chrx` the number and names of kept chrX features; for a variant that excludes gene sets the number of genes removed per set; for a variant with a non-variance selection rule the rule per omic and the quantiles (5, 25, 50, 75, 95) of the kept features' mean log2-CPM in both exports; for a residualized variant the covariates and the quantiles of the per-feature R²; for a trend-normalized variant the normalization rule per omic, the span and iterations, the quantiles over donors of the tilt statistics before and after normalization and the largest absolute correction.

#### Scenario: Feature accounting of a cohort
- **WHEN** the comparison of `cohort_30` is written
- **THEN** for `rna` it gives input 36,601 for both exports, the primary's counts after chromosome removal, abundance and variance, the variant's counts, and shared / primary-only / variant-only counts that sum to each export's kept count

#### Scenario: Accounting of the trend-ratio variant
- **WHEN** the comparison of `rna_trend_ratio` is written for a cohort
- **THEN** it records `selection_rule = trend_ratio` for `rna` and `variance` for `atac`, and the kept-mean quantiles of both exports' RNA blocks

#### Scenario: Accounting of the trend-normalized variant
- **WHEN** the comparison of `trend_normalized` is written for a cohort
- **THEN** it records `normalization = cyclic_loess` for both omics with the span and iterations per omic, and for each omic the quantiles over donors of the M-on-A slope and of the median M of the lowest and highest abundance twentieth before and after normalization

### Requirement: Composed candidate
`revised_candidate` SHALL apply the four settings together in the primary's order: gene sets removed before the abundance filter, residualization after log2-CPM, `trend_ratio` for RNA and `mean` for ATAC at the selection step. `revised_candidate_trend` SHALL apply the same four settings plus `cyclic_loess` normalization for both omics, in the primary's order (normalization after log2-CPM and before the residualization). Each manifest SHALL record every setting it applies.

#### Scenario: Candidate settings
- **WHEN** `revised_candidate` exports a cohort
- **THEN** its manifest records the three excluded gene sets, the residualization covariates, `selection_rule = {rna: trend_ratio, atac: mean}` and the primary's `k`, `min_cpm`, `prior_count`, `N` and TMM settings

#### Scenario: Trend candidate settings
- **WHEN** `revised_candidate_trend` exports a cohort
- **THEN** its manifest records the three excluded gene sets, `normalization = {rna: cyclic_loess, atac: cyclic_loess}`, the residualization covariates, `selection_rule = {rna: trend_ratio, atac: mean}` and the primary's `k`, `min_cpm`, `prior_count`, `N` and TMM settings, and `residualization.csv` and `normalization_trend.csv` are both written

### Requirement: Adoption check generated from the records
For every variant the registry marks for the adoption check (`revised_candidate`, `trend_normalized`, `revised_candidate_trend`) the step SHALL write `comparison/<cohort>/adoption_check.csv` with, per latent component of the variant's analysis and per variable of `ADOPTION_CHECK_VARIABLES`, the Spearman correlation and whether its absolute value reaches `DIAGNOSTIC_CORRELATION_THRESHOLD`, and an overall `holds` value that is true only when no cell reaches it. The check SHALL read no p-value, rejection status or CV AUROC, and the rule and variables SHALL be constants.

#### Scenario: Check holds
- **WHEN** every correlation of every candidate component with every check variable is below the threshold in absolute value
- **THEN** `adoption_check.csv` records `holds = true` for that cohort

#### Scenario: Check fails on one variable
- **WHEN** one component correlates at 0.35 with the RNA TMM factor
- **THEN** the row is flagged and `holds = false`

#### Scenario: Check written for the single-setting variant
- **WHEN** `trend_normalized` finishes a cohort
- **THEN** `comparison/<cohort>/adoption_check.csv` exists with one row per component × variable and the `holds` row, computed with the same variables and threshold as `revised_candidate`'s

### Requirement: Report attributes the axis change and states the adoption outcome
The report SHALL, per cohort, present each variant beside the primary with the comparison tables, attribute the change of the first latent axis's technical correlations to each single setting using the variants' analyses, state the adoption check's outcome for every checked variant in both cohorts, apply the precedence rule `ADOPTION_PRECEDENCE` to name which candidate, if any, a later revised-primary change proposes, and state that whether a candidate becomes a revised primary is a later change with the current primary kept as the reference.

#### Scenario: Adoption check holds in both cohorts
- **WHEN** `adoption_check.csv` records `holds = true` in `cohort_30` and `cohort_50`
- **THEN** the report states that the rule is met and that a revised-primary change may be proposed, without presenting the candidate's statistics as the result

#### Scenario: Adoption check fails
- **WHEN** `holds` is false in either cohort
- **THEN** the report states that the rule is not met, names the flagged variables, and keeps every variant as sensitivity

#### Scenario: Precedence between candidates
- **WHEN** `trend_normalized` and `revised_candidate_trend` both meet the rule in both cohorts
- **THEN** the report names `trend_normalized` as the candidate under `ADOPTION_PRECEDENCE` and states the reason the rule gives (no gene removed, no degrees of freedom spent, no technical covariate read)

### Requirement: Report carries the reading rules of the new variants
For every new variant the report SHALL state that `angle` and `shape` non-rejections stay uninformative, that a `delta` rejection is a difference and not specifically of size, and that MOTCO's operating characteristics were measured on neither residualized nor family-excluded nor trend-selected nor trend-normalized blocks, so no variant carries a Type I or power claim of its own.

#### Scenario: Non-rejection wording
- **WHEN** `shape` is not rejected in the `rna_trend_ratio` analysis of a cohort
- **THEN** the report states that the result is uninformative about shape, not that the trajectories agree in shape

#### Scenario: Non-rejection wording in the trend-normalized variant
- **WHEN** `angle` is not rejected in the `trend_normalized` analysis of a cohort
- **THEN** the report states that the result is uninformative about orientation, not that the directions agree

## ADDED Requirements

### Requirement: Trend normalization of a block
When a variant declares `cyclic_loess` for an omic, the export SHALL, after log2-CPM and before any residualization and the selection step, correct each donor's log2-CPM over the abundance-passing features by a smooth curve of its deviation from the per-feature mean over donors (M) against that mean (A), iterated against the updated mean, using no label and no donor covariate. `none` SHALL leave the values unchanged.

#### Scenario: Tilt removed
- **WHEN** a block is normalized with `cyclic_loess`
- **THEN** for every donor the median of M over each of twenty equal-count abundance bins of the normalized values is closer to zero in absolute value than the largest of the un-normalized medians over all donors and bins (a donor without tilt has noise-level medians before and after, so the comparison is not made donor by donor), and the per-donor slope of M on A is reduced in absolute value for the donor with the largest un-normalized slope

#### Scenario: Label-free
- **WHEN** the normalization of a cohort is computed
- **THEN** it reads only the block's log2-CPM values; sex, stage, the analysis covariates and the technical table do not enter

#### Scenario: Primary unchanged
- **WHEN** an omic is normalized with `none`
- **THEN** its log2-CPM values equal the un-normalized values exactly

### Requirement: Trend normalization recorded and cross-checked
A trend-normalizing export SHALL record per omic the rule, span, iterations and reference in its manifest, write `normalization_trend.csv` (per donor and omic: the M-on-A slope and the median M of the lowest and highest abundance twentieth, before and after), and the numpy normalization SHALL be cross-checked against limma's cyclic loess on the same log2-CPM by a committed script whose tolerance and outcome are recorded in the report.

#### Scenario: Trend table written
- **WHEN** `trend_normalized` exports `cohort_30`
- **THEN** `normalization_trend.csv` has one row per donor and omic (160 rows) whose donors equal `metadata.csv`'s, with the four tilt statistics before and after

#### Scenario: Cross-check against limma
- **WHEN** the cross-check script is run on the export's debug files of a cohort
- **THEN** it reports per omic the largest absolute difference between the numpy and the limma normalized values and exits non-zero when it exceeds the recorded tolerance

### Requirement: Precedence rule between adoption candidates
`ADOPTION_PRECEDENCE` SHALL be a constant written before any run: when `trend_normalized` meets `ADOPTION_RULE` in both cohorts it is the candidate a later revised-primary change proposes regardless of `revised_candidate_trend`; when only `revised_candidate_trend` meets it, that one is; when neither does, the variants stay sensitivity analyses. The rule SHALL read only the `holds` values of the adoption tables.

#### Scenario: Single setting preferred
- **WHEN** both checked variants hold in both cohorts
- **THEN** the step's manifests and the report name `trend_normalized` as the candidate

#### Scenario: Neither holds
- **WHEN** neither checked variant holds in both cohorts
- **THEN** the report states that no revised-primary change follows and names the flagged variables per variant and cohort

### Requirement: Report carries the trend-normalization limitations
When discussing `trend_normalized` or `revised_candidate_trend` the report SHALL state that the normalization assumes most features do not change at any abundance level, so a stage or sex effect shaped like a tilt would be removed with the artefact; that the curve is least determined at the rare end, where the selection draws most features; that the lncRNA-versus-coding offset at matched abundance is not corrected; and that there is no Type I check for the variant.

#### Scenario: Limitations stated
- **WHEN** the report discusses `trend_normalized`
- **THEN** it names the assumption, the rare-end uncertainty, the uncorrected gene-class offset and the absence of a Type I check

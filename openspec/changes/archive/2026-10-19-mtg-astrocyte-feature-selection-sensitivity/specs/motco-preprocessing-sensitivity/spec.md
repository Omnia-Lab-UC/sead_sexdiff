# Spec Delta

## MODIFIED Requirements

### Requirement: Inputs are the pseudobulk outputs, the annotation, the primary export and the primary analysis
For each requested variant and cohort the step SHALL read the pseudobulk directory (the files the export reads), the gene annotation, the primary export under `results/motco_input/<cohort>/` and the primary analysis under `results/motco_analysis/<cohort>/`; for a variant that excludes gene sets, residualizes or composes them it SHALL also read `results/pseudobulk/donor_technical.csv` and `technical_manifest.json`; and nothing else. When one is missing it SHALL stop with an error naming the file or directory and the command that produces it, and write nothing.

#### Scenario: Pseudobulk outputs not synced
- **WHEN** `results/pseudobulk/donor_covariates.csv` is absent
- **THEN** the step stops with an error naming the file and the `covariates` command, and writes nothing

#### Scenario: Primary analysis not yet run
- **WHEN** `results/motco_analysis/cohort_30/` does not exist
- **THEN** the step stops with an error naming the directory and the `analyze` command, and writes nothing

#### Scenario: Primary export not yet run
- **WHEN** `results/motco_input/cohort_50/` does not exist
- **THEN** the step stops with an error naming the directory and the `export` command, and writes nothing

#### Scenario: Technical table absent for a residualized variant
- **WHEN** `residualized` is requested and `results/pseudobulk/donor_technical.csv` is absent
- **THEN** the step stops with an error naming the file and the `technical-covariates` command, and writes nothing

#### Scenario: Technical table not needed
- **WHEN** only `keep_chrx` is requested and the technical table is absent
- **THEN** the step runs

### Requirement: Predeclared variant registry
The variants SHALL be a fixed registry, not a runtime choice of settings. The registry SHALL hold `keep_chrx`, `exclude_families`, `rna_trend_ratio`, `atac_mean_accessibility`, `residualized` and `revised_candidate`, each with a description and the export settings it changes; every setting not listed for a variant SHALL equal the primary export's. The step SHALL run only registered variants and SHALL NOT add, drop or alter a variant from the results.

#### Scenario: keep_chrx settings
- **WHEN** the `keep_chrx` export of a cohort finishes
- **THEN** its manifest records excluded chromosomes `chrY` and `chrM` only, the variant name, and the same `min_cpm`, `k`, `prior_count`, `N` and TMM settings as the primary export

#### Scenario: Unregistered variant
- **WHEN** a variant name not in the registry is requested
- **THEN** the step stops with an error naming the registered variants, and writes nothing

#### Scenario: Registry listed in the help
- **WHEN** `preprocessing-sensitivity --help` is shown
- **THEN** `--variant` lists the six names and `all`

### Requirement: Variant export is the primary pipeline with one setting changed
The variant export SHALL apply the primary export's steps in the primary order (restrict to the cohort, remove excluded chromosomes and excluded gene sets, CPM ≥ `min_cpm` in ≥ `k` donors, TMM, log2-CPM, residualization when declared, top `N` by the selection rule) with the primary's `k`, `min_cpm`, `prior_count`, `N` and TMM settings, differing only in the settings the registry lists for the variant, computed on the cohort's donors without sex, stage or analysis covariates. It SHALL write the same seven files as the primary export, plus `residualization.csv` when it residualizes.

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
- **THEN** its six CSVs are byte-identical to the kept primary export and its manifest has no `excluded_gene_sets`, `selection_rule` or `residualization` key beyond the primary's recorded defaults

### Requirement: Comparison carries the feature accounting
`sensitivity.md` and `sensitivity_manifest.json` SHALL give, per omic and for both the primary and the variant export, the feature counts after each export step and the numbers of kept features shared with the primary, kept only by the primary and kept only by the variant; for `keep_chrx` the number and names of kept chrX features; for a variant that excludes gene sets the number of genes removed per set; for a variant with a non-variance selection rule the rule per omic and the quantiles (5, 25, 50, 75, 95) of the kept features' mean log2-CPM in both exports; for a residualized variant the covariates and the quantiles of the per-feature R².

#### Scenario: Feature accounting of a cohort
- **WHEN** the comparison of `cohort_30` is written
- **THEN** for `rna` it gives input 36,601 for both exports, the primary's counts after chromosome removal, abundance and variance, the variant's counts, and shared / primary-only / variant-only counts that sum to each export's kept count

#### Scenario: Accounting of the trend-ratio variant
- **WHEN** the comparison of `rna_trend_ratio` is written for a cohort
- **THEN** it records `selection_rule = trend_ratio` for `rna` and `variance` for `atac`, and the kept-mean quantiles of both exports' RNA blocks

## ADDED Requirements

### Requirement: Technical table must match the pseudobulk
Before a variant that reads the technical table is run, the step SHALL check that `technical_manifest.json` records the sha256 of the pseudobulk manifest now on disk and that every cohort donor has a non-missing value for every residualization covariate, and SHALL stop naming the donor and column otherwise.

#### Scenario: Donor lacks a technical value
- **WHEN** a `cohort_30` donor has an empty ribosomal pseudobulk fraction
- **THEN** the `residualized` run stops naming the donor and the column, and writes nothing

### Requirement: Gene-family exclusion variant
`exclude_families` SHALL remove, before the abundance filter, every RNA gene matched by the gene sets `ribosomal`, `mitochondrial` and `hemoglobin` of `TECHNICAL_GENE_SETS`, with status `excluded_gene_set` in `features_rna.csv`, and SHALL leave the ATAC block's settings equal to the primary's. The manifest SHALL record the matched genes per set.

#### Scenario: Ribosomal genes removed
- **WHEN** `exclude_families` exports the RNA block of a cohort
- **THEN** no kept gene matches the ribosomal pattern, every matched gene carries `excluded_gene_set`, the status counts sum to 36,601, and `features_atac.csv` equals the primary's

### Requirement: Selection rules
The selection step SHALL support three rules: `variance` (the primary's), `trend_ratio` (variance divided by the median variance of twenty equal-count bins of mean log2-CPM, linearly interpolated at each feature's mean) and `mean` (mean log2-CPM). Each SHALL keep the top `N` features, ties broken by name, and mark dropped features `selection` (`variance` keeps the primary's `variance` status).

#### Scenario: Trend ratio values
- **WHEN** a block is ranked by `trend_ratio`
- **THEN** `features_<omic>.csv` holds each abundance-passing feature's mean, variance, trend and ratio, and the kept features are the `N` highest ratios

#### Scenario: Primary status unchanged
- **WHEN** a block is ranked by `variance`
- **THEN** dropped features carry `variance`, as in the primary export

### Requirement: Selection-rule variants
`rna_trend_ratio` SHALL apply `trend_ratio` to RNA and leave ATAC at `variance`; `atac_mean_accessibility` SHALL apply `mean` to ATAC and leave RNA at `variance`. The manifest SHALL record the rule per omic and, for `trend_ratio`, the twenty bin edges and medians.

#### Scenario: Trend ratio recorded
- **WHEN** `rna_trend_ratio` exports a cohort
- **THEN** the manifest records `selection_rule = {rna: trend_ratio, atac: variance}` with twenty bin edges and medians, and `features_atac.csv` equals the primary's

#### Scenario: Mean accessibility
- **WHEN** `atac_mean_accessibility` exports a cohort
- **THEN** the kept peaks are the 5,000 with the highest mean log2-CPM and `features_rna.csv` equals the primary's

### Requirement: Residualized variant
`residualized` SHALL, after log2-CPM and before the selection step, regress each abundance-passing feature on an intercept and `RESIDUALIZATION_COVARIATES`, replace its values by residual plus feature mean, and rank by the residuals' variance. `RESIDUALIZATION_COVARIATES` SHALL be the per-omic log10 library size and the gene-set pseudobulk fractions of the technical table, minus any fraction the technical-axis report classifies as cell state.

#### Scenario: Residuals orthogonal to the design
- **WHEN** `residualized` exports a cohort
- **THEN** every kept feature's values have zero Pearson correlation (within floating-point tolerance) with every design column

#### Scenario: Units kept
- **WHEN** a feature is residualized
- **THEN** its exported mean equals its log2-CPM mean before residualization

### Requirement: Residualization design recorded and label-free
A residualizing export SHALL write the design matrix as `residualization.csv` (the cohort's donors in export order) and record the covariates and the per-feature R² quantiles in the manifest. The design SHALL contain no sex, stage or analysis covariate.

#### Scenario: Design excludes labels
- **WHEN** the residualization design is built
- **THEN** its columns are exactly `RESIDUALIZATION_COVARIATES` and none is `sex`, `stage` or a column of `ANALYSIS_COVARIATES`

#### Scenario: Design file written
- **WHEN** `residualized` exports `cohort_30`
- **THEN** `residualization.csv` has 80 rows whose `donor_id` column equals `metadata.csv`'s

### Requirement: Composed candidate
`revised_candidate` SHALL apply the four settings together in the primary's order: gene sets removed before the abundance filter, residualization after log2-CPM, `trend_ratio` for RNA and `mean` for ATAC at the selection step. Its manifest SHALL record all four settings.

#### Scenario: Candidate settings
- **WHEN** `revised_candidate` exports a cohort
- **THEN** its manifest records the three excluded gene sets, the residualization covariates, `selection_rule = {rna: trend_ratio, atac: mean}` and the primary's `k`, `min_cpm`, `prior_count`, `N` and TMM settings

### Requirement: Adoption check generated from the records
For `revised_candidate` the step SHALL write `comparison/<cohort>/adoption_check.csv` with, per latent component of the candidate's analysis and per variable of `ADOPTION_CHECK_VARIABLES`, the Spearman correlation and whether its absolute value reaches `DIAGNOSTIC_CORRELATION_THRESHOLD`, and an overall `holds` value that is true only when no cell reaches it. The check SHALL read no p-value, rejection status or CV AUROC, and the rule and variables SHALL be constants.

#### Scenario: Check holds
- **WHEN** every correlation of every candidate component with every check variable is below the threshold in absolute value
- **THEN** `adoption_check.csv` records `holds = true` for that cohort

#### Scenario: Check fails on one variable
- **WHEN** one component correlates at 0.35 with the RNA TMM factor
- **THEN** the row is flagged and `holds = false`

### Requirement: Report attributes the axis change and states the adoption outcome
The report SHALL, per cohort, present each variant beside the primary with the comparison tables, attribute the change of the first latent axis's technical correlations to each single setting using the variants' analyses, state the adoption check's outcome for the candidate in both cohorts, and state that whether the candidate becomes a revised primary is a later change with the current primary kept as the reference.

#### Scenario: Adoption check holds in both cohorts
- **WHEN** `adoption_check.csv` records `holds = true` in `cohort_30` and `cohort_50`
- **THEN** the report states that the rule is met and that a revised-primary change may be proposed, without presenting the candidate's statistics as the result

#### Scenario: Adoption check fails
- **WHEN** `holds` is false in either cohort
- **THEN** the report states that the rule is not met, names the flagged variables, and keeps every variant as sensitivity

### Requirement: Report carries the reading rules of the new variants
For every new variant the report SHALL state that `angle` and `shape` non-rejections stay uninformative, that a `delta` rejection is a difference and not specifically of size, and that MOTCO's operating characteristics were measured on neither residualized nor family-excluded nor trend-selected blocks, so no variant carries a Type I or power claim of its own.

#### Scenario: Non-rejection wording
- **WHEN** `shape` is not rejected in the `rna_trend_ratio` analysis of a cohort
- **THEN** the report states that the result is uninformative about shape, not that the trajectories agree in shape

### Requirement: Report carries the residualization limitations
When discussing `residualized` or `revised_candidate` the report SHALL name the residualization covariates, state the degrees of freedom the residualization spends that the permutation null does not see, state that it removes whatever part of the stage and sex signal the technical set carries, and state that the analysis covariates were left unchanged so library-size proxies are adjusted twice.

#### Scenario: Residualized limitations stated
- **WHEN** the report discusses `residualized`
- **THEN** it names the covariates, the degrees of freedom spent, the double adjustment and the absence of a Type I check for the variant

# motco-preprocessing-sensitivity Specification

## Purpose

Measures how the covariate-adjusted sex × ADNC-stage trajectory comparison of the astrocyte cohorts depends on a preprocessing choice of the export, by rerunning the export with one predeclared setting changed and the analysis unchanged beside the primary analysis, and records the comparison so that it is read as sensitivity to the primary analysis, never as its replacement.

## Requirements

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

### Requirement: Primary export must match the inputs and the constants
Before any variant runs, the step SHALL check that the primary export manifest records the same sha256 for the pseudobulk manifest, the covariates manifest and the annotation as the files now on disk, the same `min_cpm`, `prior_count`, `N`, TMM settings and excluded chromosomes as the constants, and the same `k`, cells and donor count as the donor table now gives. On any difference it SHALL stop with an error naming the file or setting and the `export` command.

#### Scenario: Pseudobulk regenerated after the primary export
- **WHEN** `results/pseudobulk/manifest.json` has a sha256 different from the one in the primary export manifest
- **THEN** the step stops with an error naming `manifest.json` and the `export` command, and writes nothing

#### Scenario: Annotation differs
- **WHEN** the GTF under `<data root>/annotation/` has a sha256 different from the one in the primary export manifest
- **THEN** the step stops with an error naming the annotation file, and writes nothing

#### Scenario: Matching primary export
- **WHEN** every hash and setting matches
- **THEN** the variant export proceeds with the same `k` and cells as the primary export

### Requirement: Primary analysis must match the primary export and the settings
Before any variant runs, the step SHALL check that the primary analysis manifest records the same sha256 for each of the five primary export files as the files now on disk, the same MOTCO version as the one installed, the same seed and permutation count as the step uses, and CV rank selection. On any difference it SHALL stop with an error naming the file or setting and the `analyze` command.

#### Scenario: Primary export regenerated after the primary analysis
- **WHEN** `results/motco_input/cohort_30/rna.csv` has a sha256 different from the one in `results/motco_analysis/cohort_30/analysis_manifest.json`
- **THEN** the step stops with an error naming `rna.csv` and the `analyze` command, and writes nothing

#### Scenario: Matching primary analysis
- **WHEN** every hash and setting matches
- **THEN** the primary row of the comparison carries the primary record's statistics, p-values, null summaries and CV-selected rank unchanged

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

### Requirement: Variant analysis uses the primary analysis settings
The variant analysis SHALL run through MOTCO's trajectory-analysis entry point with the same blocks and units, stage order, six covariates in order, missing-covariate policy, orientation estimator, CV rank selection at the same settings, permutation count and seed as the primary analysis, so that the only difference is the exported blocks. It SHALL write the same seven files as the primary analysis.

#### Scenario: Record carries the settings
- **WHEN** the `keep_chrx` analysis of `cohort_30` finishes
- **THEN** its `result.json` records `component_selection = cv`, `permutations = 999`, the analysis seed, `orientation_estimator = net_displacement` and the six covariate names in order, and its manifest's preprocessing block records excluded chromosomes `chrY` and `chrM`

#### Scenario: Rank selected on the variant
- **WHEN** CV selects a rank different from the primary analysis's
- **THEN** the record carries that rank, the comparison names it beside `angle`, and the manifest records a finding that the ranks differ

### Requirement: Output layout mirrors the primary pipeline
Each variant SHALL be written under `results/motco_sensitivity/<variant>/`: the variant export to `motco_input/<cohort>/`, the variant analysis to `motco_analysis/<cohort>/` and the comparison to `comparison/<cohort>/` (`sensitivity.csv`, `sensitivity.md`, `sensitivity_manifest.json`). A cohort SHALL count as complete when its comparison directory exists. The step SHALL refuse to overwrite any of a cohort's three directories unless forced.

#### Scenario: Rerun without force
- **WHEN** `results/motco_sensitivity/keep_chrx/motco_analysis/cohort_30/` exists and the step runs again for that cohort without `--force`
- **THEN** it stops with an error naming the directory and writes nothing

#### Scenario: Variant roots usable by the rank sweep
- **WHEN** the `rank-sensitivity` command is given `results/motco_sensitivity/keep_chrx/motco_input` and `results/motco_sensitivity/keep_chrx/motco_analysis` as its export and analysis roots
- **THEN** it reads the variant export and analysis of a cohort as it reads the primary ones, and its reference check passes

### Requirement: The primary export and analysis are never modified
The step SHALL NOT write under `results/motco_input/` or `results/motco_analysis/`.

#### Scenario: Primary directories after a run
- **WHEN** the step finishes for both cohorts
- **THEN** every file under `results/motco_input/` and `results/motco_analysis/` is byte-identical to before the run

### Requirement: Comparison table generated from the records
The step SHALL write `comparison/<cohort>/sensitivity.csv` and its Markdown rendering `sensitivity.md`, generated from the records and never typed by hand, with the primary row first (`variant = primary`) and the variant row second. Each row SHALL give the rank selection and rank, and per statistic (`delta`, `angle`, `shape`) the observed value, the RRPP p-value, whether p ≤ α and the null mean, sd and q95, plus the CV AUROC and the row's directory.

#### Scenario: Table rows
- **WHEN** the comparison of `cohort_30` is written
- **THEN** `sensitivity.csv` has two rows, the first with `variant = primary` and the primary's rank, the second with `variant = keep_chrx` and the variant's rank

#### Scenario: Table values equal the records
- **WHEN** a cell of `sensitivity.csv` is compared with the `result.json` of its row's directory
- **THEN** they are equal, floating-point cells being written in their shortest round-trip representation

#### Scenario: Runtime and directories
- **WHEN** the table of `cohort_30` is written
- **THEN** it has no runtime column, its variant row names `motco_analysis/cohort_30` relative to the variant root and its primary row names the primary analysis directory, while `sensitivity.md` and the manifest carry each run's runtime

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

### Requirement: Comparison manifest records provenance
`sensitivity_manifest.json` SHALL record the variant name and settings, the cohort, the seed, permutation count, α and `n_jobs`, the sha256 of the pseudobulk manifest and the annotation, of the primary export's files and manifest, of the primary analysis's `result.json` and `analysis_manifest.json`, of the variant export's files and of the variant analysis's two files, every directory, the MOTCO, Python and library versions, the git revision and findings.

#### Scenario: Manifest after a cohort run
- **WHEN** the `keep_chrx` run of `cohort_30` finishes
- **THEN** its manifest is valid JSON naming the variant, listing excluded chromosomes `chrY` and `chrM`, `permutations = 999`, the seed, the MOTCO version, hashes for the pseudobulk manifest, the annotation, the five primary export files, the two primary analysis files, the five variant export files and the two variant analysis files

#### Scenario: Findings
- **WHEN** a statistic's rejection status at α differs between the variant and the primary, or the CV-selected ranks differ, or `n_jobs` is not 1
- **THEN** the manifest records a finding for each, naming the statistic or rank and calling a rejection difference preprocessing dependence of the sensitivity analysis

### Requirement: Reproducible run
Two runs of a variant on the same pseudobulk outputs, annotation, primary export, primary analysis, MOTCO release, seed and `n_jobs` SHALL give identical statistics, p-values, rank and null summaries, and byte-identical variant export CSVs, `latent_scores.csv`, `pls_loadings.csv`, `stage_means.csv`, `analysis_metadata.csv` and `sensitivity.csv`.

#### Scenario: Rerun with force
- **WHEN** `keep_chrx` is run twice for `cohort_50` with `--force` and the same `n_jobs`
- **THEN** the two `result.json` files have equal statistics, p-values, rank and null summaries, and the CSVs are byte-identical

### Requirement: Parallel runs are recorded
Because MOTCO's RRPP draws depend on `n_jobs`, the comparison manifest SHALL record `n_jobs`, a run with `n_jobs` other than 1 SHALL carry a finding that its permutation draws differ from an `n_jobs = 1` run, and the report SHALL use `n_jobs = 1`.

#### Scenario: Run with several workers
- **WHEN** `keep_chrx` is run with `--n-jobs 4`
- **THEN** its comparison manifest records `n_jobs = 4` and a finding that the p-values and null summaries differ from an `n_jobs = 1` run

### Requirement: Committed preprocessing-sensitivity report
A Markdown report under `reports/` SHALL be committed, giving for each cohort the comparison table (primary row and variant row) and the feature accounting pasted from the generated files, the MOTCO version, seed, permutation count and `n_jobs`, and a pointer to the primary report for the covariate block. Every number SHALL be checked against the matching `result.json` or manifest before commit. `cohort_30` SHALL be presented as primary and `cohort_50` as sensitivity.

#### Scenario: Report matches the records
- **WHEN** the report is compared with the primary and variant `result.json` files and the comparison manifests
- **THEN** every statistic, p-value, null summary, rank and feature count in the report equals the record's value

### Requirement: Report keeps the primary analysis primary
The report SHALL state that the primary analysis remains the primary result and that no variant replaces it or is chosen after seeing the results. A statistic whose rejection status differs between the variant and the primary SHALL be reported as preprocessing dependence of the sensitivity analysis, not as a result of the variant. A CV-selected rank that differs SHALL be reported beside `angle`, not read as better or worse.

#### Scenario: Rejection in the variant only
- **WHEN** `delta` is rejected at α = 0.05 in the `keep_chrx` analysis but not in the primary analysis
- **THEN** the report records the rejection as preprocessing dependence of the sensitivity analysis, keeps the primary non-rejection as the primary result, and does not present the variant as the analysis

#### Scenario: No rejection in either
- **WHEN** no statistic is rejected in the primary or the variant analysis of a cohort
- **THEN** the report states that keeping chrX gives no evidence of a difference either and that this is not evidence of agreement in size, direction or shape

### Requirement: Report carries the reading rules and the chrX limitations
The report SHALL state that `angle` and `shape` non-rejections stay uninformative in the variant, that a `delta` rejection is a difference and not specifically of size, that with chrX kept features that differ between the sexes by X dosage enter the blocks, that MOTCO's operating characteristics were not measured with such features, and that the variant carries no Type I or power claim of its own.

#### Scenario: Non-rejection wording in the variant
- **WHEN** `angle` is not rejected in the `keep_chrx` analysis of a cohort
- **THEN** the report states that the result is uninformative about orientation, not that the directions agree

#### Scenario: chrX limitation stated
- **WHEN** the report discusses the `keep_chrx` variant
- **THEN** it names the X-dosage features (such as `XIST`) that the primary export removed by design, says whether they were kept by the variant export, and states that MOTCO's Type I and power evidence covers neither such features nor this variant

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

### Requirement: Report carries the residualization limitations
When discussing `residualized` or `revised_candidate` the report SHALL name the residualization covariates, state the degrees of freedom the residualization spends that the permutation null does not see, state that it removes whatever part of the stage and sex signal the technical set carries, and state that the analysis covariates were left unchanged so library-size proxies are adjusted twice.

#### Scenario: Residualized limitations stated
- **WHEN** the report discusses `residualized`
- **THEN** it names the covariates, the degrees of freedom spent, the double adjustment and the absence of a Type I check for the variant

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

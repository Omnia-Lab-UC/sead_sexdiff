# motco-preprocessing-sensitivity Specification

## Purpose

Measures how the covariate-adjusted sex × ADNC-stage trajectory comparison of the astrocyte cohorts depends on a preprocessing choice of the export, by rerunning the export with one predeclared setting changed and the analysis unchanged beside the primary analysis, and records the comparison so that it is read as sensitivity to the primary analysis, never as its replacement.

## Requirements

### Requirement: Inputs are the pseudobulk outputs, the annotation, the primary export and the primary analysis
For each requested variant and cohort the step SHALL read the pseudobulk directory (the files the export reads) and the gene annotation, the primary export under `results/motco_input/<cohort>/` and the primary analysis under `results/motco_analysis/<cohort>/`, and nothing else. When one is missing it SHALL stop with an error naming the file or directory and the command that produces it, and write nothing.

#### Scenario: Pseudobulk outputs not synced
- **WHEN** `results/pseudobulk/donor_covariates.csv` is absent
- **THEN** the step stops with an error naming the file and the `covariates` command, and writes nothing

#### Scenario: Primary analysis not yet run
- **WHEN** `results/motco_analysis/cohort_30/` does not exist
- **THEN** the step stops with an error naming the directory and the `analyze` command, and writes nothing

#### Scenario: Primary export not yet run
- **WHEN** `results/motco_input/cohort_50/` does not exist
- **THEN** the step stops with an error naming the directory and the `export` command, and writes nothing

### Requirement: Predeclared variant registry
The variants SHALL be a fixed registry, not a runtime choice of settings. The registry SHALL hold `keep_chrx`, which removes features on chrY and chrM and keeps chrX, with every other export setting equal to the primary export's. The step SHALL run only registered variants and SHALL NOT add, drop or alter a variant from the results.

#### Scenario: keep_chrx settings
- **WHEN** the `keep_chrx` export of a cohort finishes
- **THEN** its manifest records excluded chromosomes `chrY` and `chrM` only, the variant name, and the same `min_cpm`, `k`, `prior_count`, `N` and TMM settings as the primary export

#### Scenario: Unregistered variant
- **WHEN** a variant name not in the registry is requested
- **THEN** the step stops with an error naming the registered variants, and writes nothing

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
The variant export SHALL apply the primary export's steps in the primary order (restrict to the cohort, remove excluded chromosomes, CPM ≥ `min_cpm` in ≥ `k` donors, TMM, log2-CPM, top `N` by variance) with the primary's `k`, `min_cpm`, `prior_count`, `N` and TMM settings, differing only in the excluded chromosomes, computed on the cohort's donors without sex, stage or covariates. It SHALL write the same seven files as the primary export.

#### Scenario: chrX features pass through the filters
- **WHEN** `keep_chrx` exports the RNA block of `cohort_30`
- **THEN** no chrX gene in `features_rna.csv` has status `sex_chromosome`, every chrY and chrM gene does, chrX genes carry `abundance`, `variance` or `kept`, the status counts sum to 36,601, and the `kept` count equals the number of columns of `rna.csv`

#### Scenario: Whole pipeline recomputed
- **WHEN** `keep_chrx` exports a cohort
- **THEN** the abundance filter, TMM factors, log2-CPM and variance ranking are computed with chrX features present, and the block holds `N` features (or all that pass, with a finding), not the primary's features plus chrX

#### Scenario: Donor tables unchanged
- **WHEN** `keep_chrx` exports a cohort
- **THEN** its `metadata.csv` and `covariates.csv` are byte-identical to the primary export's

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
`sensitivity.md` and `sensitivity_manifest.json` SHALL give, per omic and for both the primary and the variant export, the feature counts after each export step, and for the variant the number of kept features on chrX and the numbers of kept features shared with the primary, kept only by the primary and kept only by the variant. The manifest SHALL list the kept chrX feature names.

#### Scenario: Feature accounting of a cohort
- **WHEN** the comparison of `cohort_30` is written
- **THEN** for `rna` it gives input 36,601 for both exports, the primary's counts after chromosome removal, abundance and variance, the variant's counts, the variant's kept chrX genes, and shared / primary-only / variant-only counts that sum to each export's kept count

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

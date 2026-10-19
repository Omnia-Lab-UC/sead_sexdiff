# latent-axis-diagnostic Specification

## Purpose

Measures, reproducibly and from the kept outputs only, what the MOTCO latent axes and the exported blocks' principal axes correlate with (technical, covariate and design variables), how much of the trajectory geometry lies along each axis, and what the feature selection kept, so that preprocessing decisions can be grounded in recorded diagnostics rather than in scratch checks.

## Requirements

### Requirement: Inputs are kept outputs only
Per cohort the diagnostic SHALL read the primary export (`results/motco_input/<cohort>/`), the primary analysis (`results/motco_analysis/<cohort>/`), the technical table (`results/pseudobulk/donor_technical.csv`) and, when present, the forced-rank analyses under `results/motco_rank_sensitivity/<cohort>/rank_*/`, and nothing else. A missing required input SHALL stop the run with an error naming the file and the command that produces it. It SHALL write nothing into the directories it reads.

#### Scenario: Technical table absent
- **WHEN** `results/pseudobulk/donor_technical.csv` is missing
- **THEN** the diagnostic stops with an error naming the file and the `technical-covariates` command, and writes nothing

#### Scenario: Rank sweep absent
- **WHEN** `results/motco_rank_sensitivity/cohort_30/` does not exist
- **THEN** the diagnostic of `cohort_30` runs on the primary analysis alone and records a finding that no forced-rank analyses were read

### Requirement: Analysis and export must match
Before writing, the diagnostic SHALL verify that each analysis manifest's recorded input hashes equal the sha256 of the export files it reads, and SHALL stop with an error naming the file and the recorded and computed hashes when they differ.

#### Scenario: Export regenerated after the analysis
- **WHEN** `rna.csv` of a cohort differs from the hash in the primary analysis manifest
- **THEN** the diagnostic stops naming `rna.csv` and both hashes, and writes nothing

### Requirement: Predeclared variable set
The variables SHALL be the constant `DIAGNOSTIC_VARIABLES`: the technical table's donor columns, the export manifest's TMM factor and log10 library size per omic, the six analysis covariates (nuclei counts on the log10 scale), stage (ordinal) and sex (indicator). The set SHALL NOT be chosen or extended from results.

#### Scenario: Variables recorded
- **WHEN** the diagnostic of a cohort is written
- **THEN** its manifest lists every variable with its source file or derivation

### Requirement: Axis correlations
For every latent component of every analysis read, the diagnostic SHALL write the Spearman correlation with each variable, the partial Spearman correlation given the six analysis covariates, and the component's score standard deviation; for the first three principal components of each z-scored block it SHALL write the same correlations and the variance share. Each correlation whose absolute value is at least `DIAGNOSTIC_CORRELATION_THRESHOLD` SHALL be flagged.

#### Scenario: Primary rank-1 analysis
- **WHEN** the primary analysis of a cohort selected rank 1
- **THEN** `axis_correlations.csv` holds one row per variable for `pls_0`, with raw and partial correlations and the flag

#### Scenario: Nested components across forced ranks
- **WHEN** the forced-rank analyses are read
- **THEN** the table holds the correlations of every component at every rank, and a finding states whether `pls_0` is identical across ranks

### Requirement: Geometry shares along components
From each analysis's `stage_means.csv` the diagnostic SHALL write, per component, the share of the squared length of each sex's net displacement, of the per-stage sex-difference vectors and of the stage-centred sex-difference vectors that lies along that component.

#### Scenario: Rank-1 shares
- **WHEN** the analysis has one component
- **THEN** every share equals 1 and the component's absolute displacements per sex are recorded

### Requirement: Feature-selection accounting
From the export's `features_<omic>.csv` and block files the diagnostic SHALL write, per omic: the quantiles (5, 25, 50, 75, 95) of the kept features' mean log2-CPM, the share of kept features with mean below 0 and below 1, the share of the abundance-passing variance held by the kept features, the variance at the cut (kept minimum, dropped median), and the share of kept features whose absolute Spearman correlation with the first latent component is at least `DIAGNOSTIC_CORRELATION_THRESHOLD`.

#### Scenario: RNA accounting
- **WHEN** the diagnostic of `cohort_30` is written
- **THEN** `feature_selection.csv` has an `rna` row with the five quantiles, both shares, the variance share, the two cut values and the correlated share

### Requirement: Diagnostic outputs and manifest
Per cohort the diagnostic SHALL write `results/motco_diagnostics/<cohort>/` with `axis_correlations.csv`, `block_pcs.csv`, `geometry_shares.csv`, `feature_selection.csv`, `diagnostics.md` (tables rendered from the CSVs, flagged correlations listed first) and `diagnostics_manifest.json` (every input path and sha256, the variable set, the threshold, the analyses read, findings, software versions, git revision). Reruns on the same inputs SHALL produce byte-identical CSVs.

#### Scenario: Reproducible run
- **WHEN** the diagnostic of a cohort is run twice into two directories
- **THEN** the four CSVs are byte-identical

### Requirement: Labels are read but feed no choice
The diagnostic reads stage and sex. Its outputs SHALL NOT be an input of `export`, `analyze`, `rank-sensitivity` or `preprocessing-sensitivity`, and the documentation SHALL state that preprocessing variants are predeclared in a change proposal, never derived from a diagnostic value.

#### Scenario: No downstream read
- **WHEN** any export or analysis command runs
- **THEN** it reads nothing under `results/motco_diagnostics/`

### Requirement: Committed technical-axis report
A report `reports/motco-astrocyte-technical-axis-<date>.md` SHALL paste both cohorts' `diagnostics.md` and the donor-level technical summary, with provenance (repository revision, MOTCO version, the primary export and analysis hashes, the technical manifest hash), and every number in it SHALL be checkable against the generated files.

#### Scenario: Numbers match the records
- **WHEN** the report is committed
- **THEN** a check finds every correlation, share, quantile and count of its tables in the matching cohort's CSVs

### Requirement: Report reads axis correlations as measurement-space findings
The report SHALL state that an axis correlated with a technical variable is a finding about the measurement space of the primary analysis, not a reinterpretation of its statistics or p-values; that the primary analysis stays primary and its non-rejections keep their reading rules; and that no preprocessing variant is chosen in this report.

#### Scenario: Composition axis found
- **WHEN** the first latent component is flagged against the ribosomal fraction and the TMM factor
- **THEN** the report describes the measurement space and names the residualization candidate without restating `delta`, `angle` or `shape` as different

### Requirement: Report classifies gene-set fractions by within-donor share
The report SHALL read a gene-set fraction whose between-donor variance dominates and whose within-donor supertype share is small as a preparation property and a residualization candidate, and one with a large within-donor supertype share as cell state that is not a candidate, using the technical table's recorded shares.

#### Scenario: Within-donor biology
- **WHEN** a gene-set fraction's within-donor supertype share exceeds its between-donor share
- **THEN** the report reads it as cell state and does not list it as a residualization candidate

#### Scenario: Donor-level preparation property
- **WHEN** a gene-set fraction's between-donor share dominates and its supertype share is small in most donors
- **THEN** the report lists it as a residualization candidate for the next change

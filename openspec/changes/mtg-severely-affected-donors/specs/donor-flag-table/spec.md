# Spec Delta

## Purpose

Records SEA-AD's donor-level flags and staging labels (the severely-affected flag, cognitive status, Braak, CERAD, Thal, LATE and the continuous pseudo-progression score) from the donor objects' `obs` into a table kept apart from the covariate table, so the donor-exclusion sensitivity analysis can select donors by the data producers' own flag without anyone re-reading the raw objects.

## ADDED Requirements

### Requirement: Flags read from donor object metadata only
The step SHALL read, for every donor object of each requested omic, the donor-level `obs` columns `Donor ID`, `Severely Affected Donor`, `Cognitive Status`, `Braak`, `CERAD score`, `Thal` and `LATE`, plus `Continuous Pseudo-progression Score` from the objects that hold it, without loading any count or expression matrix. Each column SHALL hold one value per donor object, otherwise the run SHALL fail naming the file and the column. A missing column SHALL be an error naming the file, except the pseudo-progression score, whose absence is recorded.

#### Scenario: Matrices are not loaded
- **WHEN** the step runs over all RNA and ATAC donor objects
- **THEN** only `obs` columns are read and the step finishes without loading `X` or any layer

#### Scenario: Score absent in a release
- **WHEN** the RNA objects lack `Continuous Pseudo-progression Score`
- **THEN** the step records a finding naming the column and the omic, takes the score from the ATAC objects, and continues

### Requirement: Flag coded, labels recorded, Reference donors excluded
The step SHALL write one row per non-Reference donor with `donor_id`, `severely_affected` (true for `Y`, false for `N`; any other value SHALL stop the run naming the donor and the value), and the labels `cognitive_status`, `braak`, `cerad_score`, `thal`, `late` and `cps` as read, not coded. Donors whose staging labels read `Reference` SHALL be left out of the table and listed in the manifest. Nothing SHALL be imputed.

#### Scenario: Flagged donor
- **WHEN** a donor's `Severely Affected Donor` is `Y`
- **THEN** its row has `severely_affected = true`

#### Scenario: Unexpected flag value
- **WHEN** a donor's `Severely Affected Donor` is `Unknown`
- **THEN** the run stops with an error naming the donor and the value, and writes nothing

### Requirement: Cross-omic agreement
For donors with both an RNA and an ATAC object the step SHALL compare every column between the two objects. A disagreement on `Severely Affected Donor` SHALL stop the run naming the donor and both values. Any other disagreement SHALL NOT stop the run; it SHALL be written to the manifest's findings with the donor, the column and both values, and the RNA value SHALL be used.

#### Scenario: Flag disagrees between omics
- **WHEN** a donor's RNA object has `Y` and its ATAC object has `N`
- **THEN** the run stops with an error naming the donor and both values, and writes nothing

#### Scenario: Label disagrees between omics
- **WHEN** a donor's `Braak` differs between the two objects
- **THEN** the table holds the RNA value and the findings record the donor, `Braak` and both values

### Requirement: Flag table and manifest apart from the covariate table
The step SHALL write `donor_flags.csv` and `donor_flags_manifest.json` into the cell type's pseudobulk directory and SHALL NOT write `donor_covariates.csv`, `covariates_manifest.json` or any other existing file. The manifest SHALL record the cell type, the source columns and their coded names, the objects read per omic, the excluded Reference donors, the findings, software versions and the git revision.

#### Scenario: Existing files untouched
- **WHEN** the step completes in `results/microglia/pseudobulk/`
- **THEN** `donor_flags.csv` and `donor_flags_manifest.json` exist, and every file that existed before has the same sha256

### Requirement: Flag table is an input of the exclusion step only
No export, analysis, diagnostic or sensitivity step other than the donor-exclusion sensitivity SHALL read the flag table, and the flag SHALL enter no model as a covariate or a variable. The exclusion step SHALL record the table's sha256 in every manifest it writes.

#### Scenario: Primary export unchanged by the table
- **WHEN** `export` runs after `donor-flags` has written the table
- **THEN** the export reads no flag file and its outputs equal the outputs without the table

# donor-covariate-table Specification

## Purpose

Records donor-level covariates (age, PMI, RIN, education) from the SEA-AD MTG donor objects as a numerically coded table. The covariate-adjusted trajectory analysis can then use them without anyone re-reading the raw objects.

## Requirements

### Requirement: Covariates read from donor object metadata only
The system SHALL read, for every donor object of each requested omic, the donor-level `obs` columns `Donor ID`, `Age at Death`, `PMI`, `Years of education` and `Highest level of education`, plus `RIN` from the RNA objects only. It SHALL do so without loading any count or expression matrix. Each column SHALL hold exactly one value per donor object, and otherwise the run SHALL fail with an error naming the file and the column. A missing column SHALL be an error naming the file, except `RIN` in ATAC objects, which is not read.

#### Scenario: Matrices are not loaded
- **WHEN** the covariate step runs over all 84 RNA and 86 ATAC donor objects
- **THEN** only `obs` columns are read and the step finishes without loading `X` or any layer

#### Scenario: Value varies within a donor object
- **WHEN** a donor object has two distinct `PMI` values across its nuclei
- **THEN** the run stops with an error naming the file and `PMI`

### Requirement: Numeric coding with recorded caps and missing values
The system SHALL write `results/pseudobulk/donor_covariates.csv` with one row per donor that has an RNA object and the columns `donor_id`, `age_at_death`, `age_capped`, `pmi`, `rin`, `years_education`, `education_level`. Numeric columns SHALL be parsed as floats. An age reported as `90+` (or any `N+` form) SHALL be coded as `N`, with `age_capped` true. Values that are empty, `nan`, `Unknown` or otherwise non-numeric in a numeric column SHALL be coded as missing, and the manifest SHALL list each such donor and column. The table SHALL never contain an imputed value.

#### Scenario: Capped age
- **WHEN** a donor's `Age at Death` is `90+`
- **THEN** its row has `age_at_death = 90` and `age_capped = true`

#### Scenario: Missing PMI is recorded, not imputed
- **WHEN** a donor's `PMI` is `nan`
- **THEN** its `pmi` cell is empty and the manifest lists the donor under missing `pmi`

### Requirement: Cross-omic agreement
For donors with both an RNA and an ATAC object, the system SHALL compare `Age at Death`, `PMI`, `Years of education` and `Highest level of education` between the two objects. Disagreements SHALL NOT stop the run. Each one SHALL be written to the covariate manifest's findings with the donor, the column and both values, and the RNA value SHALL be used.

#### Scenario: PMI disagrees between omics
- **WHEN** a donor's RNA object has `PMI = 8.13` and its ATAC object has `PMI = 8.5`
- **THEN** the table holds 8.13 and the findings record the donor, `PMI`, 8.13 and 8.5

#### Scenario: RNA-only run
- **WHEN** the step runs with only the RNA omic
- **THEN** the table is written, and the manifest records that no cross-omic check was done

### Requirement: Covariate manifest
The system SHALL write `results/pseudobulk/covariates_manifest.json` recording the source columns and their coded names, the coding rules, the donor objects read per omic, the missing-value list, the capped-age count, the findings, software versions and the git revision.

#### Scenario: Manifest after a full run
- **WHEN** the step completes over all donor objects
- **THEN** the manifest is valid JSON listing 84 RNA and 86 ATAC objects read, and the findings list (possibly empty)

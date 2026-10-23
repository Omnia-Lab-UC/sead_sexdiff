# Spec Delta

## MODIFIED Requirements

### Requirement: References must match
Before any fit the step SHALL check that the export manifest records the sha256 of the pseudobulk manifest now on disk, that donors, sex and stage agree between the export, the donor table and the count files, and that the analysis manifest records the export files' sha256 now on disk, the installed MOTCO version, the seed and the permutation count. When the export manifest records a donor exclusion (an exclusion-cohort export), the step SHALL remove the recorded excluded donors from the donor table's cohort before the donor-agreement check and SHALL record the exclusion in its manifests. On any difference it SHALL stop naming the file or setting and write nothing. The manifest SHALL record the sha256 of every input file read.

#### Scenario: Pseudobulk manifest changed since the export
- **WHEN** the pseudobulk `manifest.json` has a sha256 different from the one the export manifest records
- **THEN** the step stops with an error naming the file and the `export` command, and writes nothing

#### Scenario: Exclusion-cohort export
- **WHEN** the export manifest records `donor_exclusion` with three excluded donors and the export's `metadata.csv` lists the cohort donors minus those three
- **THEN** the donor-agreement check passes, the fits run on the export's donors, and the fit manifests record the exclusion

#### Scenario: Exclusion recorded but donors differ
- **WHEN** the export manifest records an exclusion and `metadata.csv` still differs from the donor table's cohort minus the excluded donors
- **THEN** the step stops with an error naming both files and the differing donors, and writes nothing

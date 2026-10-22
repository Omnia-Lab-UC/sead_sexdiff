# Spec Delta

## MODIFIED Requirements

### Requirement: Cohort summary
The system SHALL write `cohort_summary.csv` in the cell type's pseudobulk directory with, for each cohort, the count of donors per sex × stage cell and the cohort total, and SHALL print the same table at the end of the run. The manifest SHALL record whether every sex × stage cell is non-empty in each cohort, and the cell type's preregistered cohort cells from the registry, which for `microglia` are the counts of the metadata-only spike recorded before the pseudobulk ran.

#### Scenario: Summary matches the preregistered cohort
- **WHEN** the table is built from the complete MTG releases for `astrocyte`
- **THEN** the ≥30 cohort summary reports n = 80 with F 11/10/28 and M 9/10/12 across stages 0/1/2, and the ≥50 cohort reports n = 74; any deviation is written to the manifest as a finding rather than silently accepted

#### Scenario: Microglia summary matches the spike
- **WHEN** the table is built from the complete MTG releases for `microglia`
- **THEN** the cohort cells equal the registry's preregistered microglia cells (the spike's counts); any deviation is written to the manifest as a finding rather than silently accepted

#### Scenario: Empty cell is flagged
- **WHEN** a sex × stage cell has zero donors in a cohort
- **THEN** the manifest marks that cohort as not analysable and the summary prints a warning

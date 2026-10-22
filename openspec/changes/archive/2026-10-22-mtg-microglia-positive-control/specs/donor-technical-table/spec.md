# Spec Delta

## MODIFIED Requirements

### Requirement: Same nuclei as the pseudobulk
The technical step SHALL select, in every RNA donor object, exactly the nuclei the pseudobulk step of the same cell type summed: the registry's mapping of the object's release for that cell type, `Reference` donors excluded, and the `RNA_MIN_GENES` nucleus filter applied to the `UMIs` layer. The selection SHALL be made by the code path the pseudobulk uses, with the same cell type, and a donor's nucleus count SHALL equal the pseudobulk's `n_nuclei` for that donor, otherwise the run SHALL fail naming the donor and both counts. The step SHALL refuse a pseudobulk directory whose manifest records another cell type, naming both.

#### Scenario: Nucleus counts agree
- **WHEN** the technical step processes a donor whose pseudobulk cache records `n_nuclei = 812`
- **THEN** it reads 812 nuclei of the cell type from that object and records 812

#### Scenario: Excluded donor
- **WHEN** a donor object has ADNC `Reference`
- **THEN** the step records the donor as excluded with no technical row, as the pseudobulk does

#### Scenario: Cell-type mismatch
- **WHEN** `technical-covariates --cell-type microglia --pseudobulk-dir results/pseudobulk` is run against the astrocyte pseudobulk
- **THEN** the step stops with an error naming `astrocyte` and `microglia`, and writes nothing

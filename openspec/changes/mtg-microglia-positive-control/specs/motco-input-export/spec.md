# Spec Delta

## MODIFIED Requirements

### Requirement: Inputs are the derived pseudobulk outputs only
The export SHALL read `rna_counts.h5ad`, `atac_counts.h5ad`, `donor_table.csv`, `manifest.json` and `donor_covariates.csv` from the cell type's pseudobulk directory (default `results/pseudobulk/` for `astrocyte`, `results/microglia/pseudobulk/` for `microglia`), plus the gene annotation, and SHALL NOT read donor objects. It SHALL refuse to run when the pseudobulk manifest is marked partial, when its cell type differs from the requested one, when a required file is missing, or when a cohort is marked not analysable. Each refusal SHALL be an error naming the file, cohort or cell types. The export manifest SHALL record the cell type.

#### Scenario: Partial pseudobulk input
- **WHEN** `manifest.json` has `partial: true`
- **THEN** the export stops with an error saying the pseudobulk run is partial

#### Scenario: Covariates not yet synced
- **WHEN** `donor_covariates.csv` is absent
- **THEN** the export stops with an error naming the file and the `covariates` command that produces it

#### Scenario: Cell-type mismatch
- **WHEN** `export --cell-type microglia --pseudobulk-dir results/pseudobulk` is run against the astrocyte pseudobulk
- **THEN** the export stops with an error naming `astrocyte` and `microglia`, and writes nothing

#### Scenario: Astrocyte manifest without a cell type
- **WHEN** the astrocyte pseudobulk manifest was written before the registry existed and has no `cell_type` key
- **THEN** the export treats it as `astrocyte`, records that as a finding, and its outputs are unchanged

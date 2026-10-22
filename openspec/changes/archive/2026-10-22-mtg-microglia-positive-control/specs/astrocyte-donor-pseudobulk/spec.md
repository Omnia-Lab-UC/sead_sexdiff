# Spec Delta

## MODIFIED Requirements

### Requirement: Astrocyte identification through a recorded taxonomy mapping
The system SHALL identify the nuclei of the requested cell type (default `astrocyte`) through an explicit per-release mapping, taken from the cell-type registry, from the release's `Subclass` label(s) to the canonical population and from the release's `Supertype` labels of that population to canonical supertype names. For `astrocyte` the mapping SHALL be the current one (`Astrocyte` in both releases; `Astro_1`…`Astro_6-SEAAD`), so the astrocyte outputs are unchanged. The mapping and the cell type SHALL be written to the run manifest. If a donor object lacks the `Subclass` or `Supertype` column, or its labels are not covered by the mapping for that release (population or recorded non-population supertypes of the selected subclass), the run SHALL fail with an error naming the file and the unmapped label rather than silently yielding zero nuclei.

#### Scenario: Both releases map to the same population
- **WHEN** the RNA (2026-06-22) and ATAC (2024-12-06) donor objects are processed for a registered cell type
- **THEN** that population is selected in both and the manifest records the cell type and the raw-to-canonical label mapping used for each release

#### Scenario: Unmapped label is an error
- **WHEN** a donor object contains a `Subclass` value absent from the release mapping, or a `Supertype` of the selected subclass absent from both the population and the recorded non-population lists
- **THEN** the run stops with an error naming the file and the label

#### Scenario: Astrocyte default unchanged
- **WHEN** the pseudobulk runs without a cell type
- **THEN** it selects the same nuclei as before this change, and a two-donor dry run reproduces the cached per-donor counts and compositions of those donors

### Requirement: Pseudobulk outputs and run manifest
The system SHALL write, under the cell type's pseudobulk directory (`results/pseudobulk/` for `astrocyte`, `results/microglia/pseudobulk/` for `microglia`), one donors × features integer count matrix per omic in a compact loadable format that preserves feature annotations, and a run manifest (JSON) recording the cell type, the release file list and sizes, the taxonomy mapping, QC parameters, the matrix source per omic, excluded donors, software versions and the pipeline git revision. Outputs SHALL be small enough to transfer between machines (on the order of 100 MB in total).

#### Scenario: Outputs present after a full run
- **WHEN** the pseudobulk step completes for both omics
- **THEN** the cell type's pseudobulk directory contains the RNA matrix, the ATAC matrix, `supertype_composition.csv` and `manifest.json` with `cell_type`, and the two matrices contain only non-negative integers

#### Scenario: Microglia directory
- **WHEN** `pseudobulk --cell-type microglia --omic both` completes
- **THEN** the outputs are under `results/microglia/pseudobulk/` and nothing under `results/pseudobulk/` changes

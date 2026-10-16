# astrocyte-donor-pseudobulk Specification

## Purpose

Turns the SEA-AD MTG per-donor RNA and ATAC objects into per-donor raw-count astrocyte pseudobulks with reconciled taxonomy labels, recorded nucleus QC, and per-donor supertype composition, as the frozen raw input of the MOTCO case study.

## Requirements

### Requirement: One donor file at a time
The system SHALL process donor objects sequentially, loading at most one donor object into memory at a time and freeing it before the next, and SHALL never concatenate nucleus-level data across donors.

#### Scenario: Peak memory is bounded by one donor
- **WHEN** the RNA pseudobulk runs over all 84 donor files
- **THEN** peak resident memory stays within a small multiple of the largest single donor object and does not grow with the number of donors processed

### Requirement: Astrocyte identification through a recorded taxonomy mapping
The system SHALL identify astrocyte nuclei through an explicit per-release mapping from the release's `Subclass` label(s) to the canonical population `Astrocyte`, and from the release's astrocyte `Supertype` labels to canonical supertype names. The mapping SHALL be written to the run manifest. If a donor object lacks the `Subclass` or `Supertype` column, or its labels are not covered by the mapping for that release, the run SHALL fail with an error naming the file and the unmapped label rather than silently yielding zero nuclei.

#### Scenario: Both releases map to the same population
- **WHEN** the RNA (2026-06-22) and ATAC (2024-12-06) donor objects are processed
- **THEN** astrocytes are selected in both and the manifest records the raw-to-canonical label mapping used for each release

#### Scenario: Unmapped label is an error
- **WHEN** a donor object contains a `Subclass` value absent from the release mapping
- **THEN** the run stops with an error naming the file and the label

### Requirement: Reference donors excluded
The system SHALL exclude nuclei whose `Overall AD neuropathological Change` value is `Reference` and SHALL record the number of excluded donors.

#### Scenario: Reference donor present
- **WHEN** a donor object carries `Reference` as its ADNC value
- **THEN** that donor contributes no pseudobulk row and is listed as excluded in the manifest

### Requirement: RNA pseudobulk from raw UMIs with nucleus QC
The system SHALL build the RNA pseudobulk by summing the `UMIs` layer over astrocyte nuclei that have at least 200 detected genes, per donor, across all astrocyte supertypes. No gene filtering and no normalization SHALL be applied. The QC threshold SHALL be recorded in the manifest.

#### Scenario: Sum over QC-passing astrocyte nuclei
- **WHEN** a donor has 120 astrocyte nuclei of which 115 pass the 200-gene filter
- **THEN** its pseudobulk row equals the column sum of `UMIs` over those 115 nuclei and its RNA nucleus count is 115

#### Scenario: Missing UMIs layer is an error
- **WHEN** a donor object has no `UMIs` layer
- **THEN** the run stops with an error naming the file

### Requirement: ATAC pseudobulk from raw peak counts without additional nucleus QC
The system SHALL build the ATAC pseudobulk by summing the raw peak counts over all astrocyte nuclei per donor, across all astrocyte supertypes, applying no nucleus QC beyond the release's own `final-nuclei` filtering. The matrix used (`X` or a named layer) SHALL be verified to hold non-negative integers and SHALL be recorded in the manifest.

#### Scenario: Sum over all astrocyte nuclei
- **WHEN** a donor has 80 astrocyte nuclei in the ATAC object
- **THEN** its pseudobulk row is the column sum of raw peak counts over all 80 and its ATAC nucleus count is 80

#### Scenario: Non-integer counts are an error
- **WHEN** the selected matrix contains non-integer or negative values
- **THEN** the run stops with an error identifying the matrix

### Requirement: Consistent feature set across donors
The system SHALL require every donor object of an omic to carry the same ordered feature set (genes for RNA, peaks for ATAC). A mismatch SHALL be an error that names the offending file and the number of differing features.

#### Scenario: Identical peak set
- **WHEN** all 86 ATAC donor objects carry the same 218,882 peaks in the same order
- **THEN** the ATAC pseudobulk has 218,882 columns in that order

### Requirement: Per-donor nucleus counts and supertype composition
The system SHALL write, for each omic, a long-format table with one row per donor × supertype giving the number of astrocyte nuclei that contributed to the pseudobulk, and a per-donor total. Donors with zero astrocyte nuclei SHALL appear with a zero total.

#### Scenario: Composition table
- **WHEN** the RNA pseudobulk completes
- **THEN** `supertype_composition.csv` has columns `donor_id`, `omic`, `supertype`, `n_nuclei` and a total per donor equal to the RNA nucleus count in the donor table

### Requirement: Pseudobulk outputs and run manifest
The system SHALL write, under `results/pseudobulk/`, one donors × features integer count matrix per omic in a compact loadable format that preserves feature annotations, and a run manifest (JSON) recording the release file list and sizes, the taxonomy mapping, QC parameters, the matrix source per omic, excluded donors, software versions and the pipeline git revision. Outputs SHALL be small enough to transfer between machines (on the order of 100 MB in total).

#### Scenario: Outputs present after a full run
- **WHEN** the pseudobulk step completes for both omics
- **THEN** `results/pseudobulk/` contains the RNA matrix, the ATAC matrix, `supertype_composition.csv` and `manifest.json`, and the two matrices contain only non-negative integers

### Requirement: Resumable runs and donor subsets
The system SHALL cache each donor's completed pseudobulk row and composition so that a rerun skips completed donors, and SHALL accept an explicit donor list or count limit for dry runs on a subset.

#### Scenario: Rerun after interruption
- **WHEN** the run is interrupted after 40 of 84 RNA donors and started again
- **THEN** the first 40 donors are not re-read and the run completes from donor 41

#### Scenario: Two-donor dry run
- **WHEN** the step runs with a donor limit of 2
- **THEN** outputs are produced for exactly those 2 donors and the manifest marks the run as partial

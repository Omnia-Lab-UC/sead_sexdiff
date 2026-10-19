# Spec Delta

## Purpose

Records, per donor, technical quality measures of the astrocyte nuclei that entered the pseudobulks (gene-set UMI fractions, nucleus and library QC, their within-donor dispersion and supertype share), from the raw donor objects and the pseudobulk, as a table kept apart from the donor covariate table so that the primary pipeline's inputs and reference checks are unchanged.

## ADDED Requirements

### Requirement: Same nuclei as the pseudobulk
The technical step SHALL select, in every RNA donor object, exactly the nuclei the pseudobulk step summed: the taxonomy mapping of the object's release, `Reference` donors excluded, and the `RNA_MIN_GENES` nucleus filter applied to the `UMIs` layer. The selection SHALL be made by the code path the pseudobulk uses, and a donor's nucleus count SHALL equal the pseudobulk's `n_nuclei` for that donor, otherwise the run SHALL fail naming the donor and both counts.

#### Scenario: Nucleus counts agree
- **WHEN** the technical step processes a donor whose pseudobulk cache records `n_nuclei = 812`
- **THEN** it reads 812 astrocyte nuclei from that object and records 812

#### Scenario: Excluded donor
- **WHEN** a donor object has ADNC `Reference`
- **THEN** the step records the donor as excluded with no technical row, as the pseudobulk does

### Requirement: Technical reads touch no normalized matrix
The step SHALL read only `obs` columns, `var` and the selected nuclei's rows of the `UMIs` layer from RNA objects, and only `obs` columns from ATAC objects. It SHALL never load the RNA `X`.

#### Scenario: Matrices are not loaded
- **WHEN** the step runs over all RNA and ATAC donor objects
- **THEN** the RNA `X` is never read and ATAC count matrices are never read

### Requirement: Predeclared gene sets and QC columns
The gene sets SHALL be the constant `TECHNICAL_GENE_SETS`: `ribosomal` (gene names matching `^RP[LS]\d+[A-Z]?$`, `^RPLP[0-2]$` or `RPSA`), `mitochondrial` (names starting `MT-`) and `hemoglobin` (`HBB`, `HBA1`, `HBA2`). The QC columns SHALL be the constant `TECHNICAL_OBS_COLUMNS` per omic. The manifest SHALL list the matched genes per set and count, and a set that matches no gene SHALL be an error.

#### Scenario: Ribosomal set excludes pseudogenes
- **WHEN** the RNA `var` holds `RPL13A`, `RPL13AP5`, `RPLP0`, `RPSA` and `RPS6KA1`
- **THEN** `ribosomal` matches `RPL13A`, `RPLP0` and `RPSA` and matches neither `RPL13AP5` nor `RPS6KA1`

#### Scenario: Thirteen MT genes
- **WHEN** the step runs on the RNA objects
- **THEN** the manifest lists 13 genes for `mitochondrial`

### Requirement: Per-nucleus fractions, per-donor summaries
For every selected RNA nucleus the step SHALL compute, from the `UMIs` row, the UMI fraction in each gene set, and SHALL read the QC columns. Per donor it SHALL write the nucleus-weighted mean and the within-donor standard deviation of each fraction and QC column, the UMI-weighted pseudobulk fraction of each gene set, the share of each fraction's within-donor variance explained by `Supertype`, the number of nuclei and the number of distinct libraries (`load_name`).

#### Scenario: Two summaries of one fraction
- **WHEN** a donor's astrocyte nuclei have ribosomal fractions averaging 0.012 while the summed UMIs give 0.015
- **THEN** its row holds `ribosomal_mean = 0.012` and `ribosomal_pseudobulk = 0.015`

#### Scenario: Supertype share
- **WHEN** a donor's ribosomal fraction differs between its astrocyte supertypes
- **THEN** its row holds the share of that fraction's variance explained by `Supertype` in `[0, 1]`, computed within that donor only

### Requirement: Pseudobulk fractions checked against the pseudobulk file
The UMI-weighted fraction per donor and gene set SHALL also be computed from `results/pseudobulk/rna_counts.h5ad`, and the two SHALL agree to `COVARIATE_RTOL`, otherwise the run SHALL fail naming the donor and the set.

#### Scenario: Agreement with the pseudobulk file
- **WHEN** a donor's summed `UMIs` rows give a ribosomal fraction of 0.0151
- **THEN** the same fraction computed from the donor's row of `rna_counts.h5ad` equals 0.0151 within tolerance

### Requirement: ATAC library metrics per donor
From every ATAC donor object the step SHALL read, for the astrocyte nuclei, the ATAC QC columns of `TECHNICAL_OBS_COLUMNS["atac"]` and write their nucleus-weighted means, recording whether each column is constant within the object.

#### Scenario: Single-library ATAC object
- **WHEN** an ATAC object holds one value of `ATAC_TSS_enrichment_score` across its nuclei
- **THEN** the donor's row holds that value and the manifest records the column as constant for that donor

### Requirement: Outputs kept apart from the covariate table
The step SHALL write `results/pseudobulk/donor_technical.csv` (one row per donor with a technical row, `donor_id` first), `results/pseudobulk/nucleus_technical_summary.csv` (per donor and gene set: within-donor SD, supertype share, nucleus count) and `results/pseudobulk/technical_manifest.json`. It SHALL NOT write or modify `donor_covariates.csv` or `covariates_manifest.json`.

#### Scenario: Covariate files untouched
- **WHEN** the technical step completes
- **THEN** the sha256 of `donor_covariates.csv` and `covariates_manifest.json` are the same as before the run

#### Scenario: Donors with no ATAC object
- **WHEN** a donor has an RNA object but no ATAC object
- **THEN** its ATAC columns are empty and the manifest lists the donor under missing ATAC

### Requirement: Technical manifest
`technical_manifest.json` SHALL record the gene sets with matched genes, the QC columns per omic, the donor objects read per omic, the excluded donors, the per-donor nucleus counts beside the pseudobulk's, the pseudobulk manifest hash the fractions were checked against, the findings, software versions and the git revision.

#### Scenario: Manifest after a full run
- **WHEN** the step completes over all donor objects
- **THEN** the manifest is valid JSON listing 84 RNA and 86 ATAC objects read, the pseudobulk manifest sha256, and a findings list

### Requirement: Resumable runs
The step SHALL cache each donor's result under `results/pseudobulk/per_donor_technical/` and SHALL skip cached donors on rerun unless forced, so a stopped run resumes.

#### Scenario: Resume after interruption
- **WHEN** the step is rerun after stopping with 30 donors cached
- **THEN** those 30 are read from the cache and only the remaining donor objects are opened

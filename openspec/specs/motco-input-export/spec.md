# motco-input-export Specification

## Purpose

Turns the raw-count astrocyte donor pseudobulks into MOTCO-ready RNA and ATAC blocks for each cohort, using a fixed, recorded preprocessing order. Every feature can be traced from the raw matrix to the exported block.

## Requirements

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

### Requirement: Per-cohort preprocessing in a fixed order
For each cohort (`cohort_30`, `cohort_50`) and each omic, the system SHALL apply these steps in order, each computed only on that cohort's donors:
(1) restrict to the cohort's donors;
(2) remove features on chrX, chrY and chrM;
(3) keep features with CPM ≥ 1 in at least `k` donors, where `k` is the size of the cohort's smallest sex × stage cell, and CPM uses each donor's total count over the features left after step 2;
(4) compute TMM normalization factors on the filtered counts;
(5) transform to log2-CPM with prior count 2, using effective library sizes (total count × TMM factor) and a prior scaled by each donor's library size relative to the mean (edgeR `cpm(log = TRUE)` semantics);
(6) keep the `N = 5000` features with the highest variance of log2-CPM across the cohort's donors, ties broken by feature name.
No step SHALL use sex, stage or covariates.

#### Scenario: Same thresholds, separate per cohort
- **WHEN** both cohorts are exported
- **THEN** each cohort's filter, TMM factors and feature selection are computed from its own donors, and the two cohorts may keep different features

#### Scenario: Abundance threshold follows the smallest cell
- **WHEN** the `cohort_30` smallest sex × stage cell has 9 donors
- **THEN** a feature passes step 3 only if at least 9 cohort donors have CPM ≥ 1

#### Scenario: Fewer features than N survive
- **WHEN** fewer than 5000 features of an omic pass step 3
- **THEN** all of them are kept and the manifest records a finding with the number kept

### Requirement: Sex-chromosome removal through recorded annotation
The system SHALL determine each RNA gene's chromosome from a pinned gene annotation matched on the Ensembl gene id with the version suffix ignored. For ATAC it SHALL use the chromosome prefix of the peak name (`chrN:start-end`). The annotation's source URL, release and checksum SHALL be recorded. RNA genes absent from the annotation SHALL be kept and counted in the manifest findings. A peak name that does not parse SHALL be an error.

#### Scenario: Y-linked gene removed
- **WHEN** the RNA block is exported
- **THEN** no exported gene is annotated on chrX, chrY or chrM, and the manifest records how many genes were removed per chromosome

#### Scenario: chrY peaks removed
- **WHEN** the ATAC block is exported
- **THEN** none of the 180 chrY or 5,886 chrX peaks is in the exported block

### Requirement: Export contract per cohort
The system SHALL write `results/motco_input/<cohort>/` containing:
- `rna.csv` and `atac.csv`: donors × features, first column `donor_id`, string feature names, finite values only, in log2-CPM units.
- `metadata.csv`: columns `donor_id`, `sex`, `stage`, `stage_label`. `stage` holds `0`/`1`/`2` and `stage_label` holds `Not AD + Low` / `Intermediate` / `High`.
- `covariates.csv`: `donor_id`, the coded covariates, and `n_rna`, `n_atac`.
- `features_rna.csv` and `features_atac.csv`: every input feature with its chromosome, the step at which it was removed (`sex_chromosome`, `abundance`, `variance`) or `kept`, and its log2-CPM variance where computed.
All donor-indexed files SHALL list the same donors in the same order, sorted by `donor_id`.

#### Scenario: Rows align across files
- **WHEN** a cohort is exported
- **THEN** `rna.csv`, `atac.csv`, `metadata.csv` and `covariates.csv` have identical `donor_id` columns, and `cohort_30` has 80 rows

#### Scenario: Feature accounting closes
- **WHEN** `features_atac.csv` is tallied by status
- **THEN** the counts sum to 218,882 and the `kept` count equals the number of columns in `atac.csv`

### Requirement: Preprocessing manifest
The system SHALL write `results/motco_input/<cohort>/manifest.json` recording:
- the declared units per block (`log2-CPM`, MOTCO transform `identity`);
- the stage column and the explicit stage order `["0", "1", "2"]`, with labels; the group column and its labels;
- every threshold (`min_cpm`, `k` and its source cell, `prior_count`, `N`, excluded chromosomes);
- feature counts after each step per omic;
- TMM factors and library sizes per donor;
- the annotation source and checksum;
- the hash of the input pseudobulk manifest;
- covariate missingness within the cohort;
- software versions, the git revision and findings.
It SHALL also record that covariates are exported but not applied.

#### Scenario: Manifest declares the MOTCO inputs
- **WHEN** a cohort export finishes
- **THEN** its manifest names the blocks `rna` and `atac` with transform `identity`, the stage order `0, 1, 2`, and the group column `sex`

### Requirement: Reruns are deterministic
The system SHALL produce byte-identical CSV outputs when rerun on the same inputs and annotation. It SHALL refuse to overwrite an existing cohort directory unless forced.

#### Scenario: Rerun without force
- **WHEN** `results/motco_input/cohort_30/` exists and the export runs again without `--force`
- **THEN** it stops with an error naming the directory

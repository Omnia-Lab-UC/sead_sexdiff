# Spec Delta

## ADDED Requirements

### Requirement: Arms run per registered cell type
The step SHALL take a registered cell type (default `astrocyte`) and SHALL run the registered arms, with the covariate set, conditioning, mediator check, seed, permutation count and CV rank selection unchanged, on that cell type's primary export against that cell type's primary analysis. The registry of arms and the covariate set SHALL NOT vary with the cell type. The default SHALL reproduce the astrocyte outputs byte for byte.

#### Scenario: Microglia arms complete
- **WHEN** the step runs with the cell type `microglia` for `cohort_30`
- **THEN** two arm records exist under the microglia results root, `extended_covariates` and `latent_adjusted`, both listing the eight covariates in the registered order, with the microglia primary analysis as the reference row, the primary's seed and 999 permutations

#### Scenario: Default unchanged
- **WHEN** the step runs without a cell type on the astrocyte inputs
- **THEN** every output equals the output of the step before this change, apart from the recorded cell type in its manifests

## MODIFIED Requirements

### Requirement: Inputs are the primary export, the primary analysis and the technical table
For each requested cohort the step SHALL read only, under the cell type's results root (`results/` for `astrocyte`, `results/<name>/` otherwise):
- the primary export (`rna.csv`, `atac.csv`, `metadata.csv`, `covariates.csv`, `manifest.json` under `<root>/motco_input/<cohort>/`);
- the primary analysis (`result.json`, `analysis_manifest.json` under `<root>/motco_analysis/<cohort>/`);
- the technical table with its manifest and the pseudobulk manifest (`donor_technical.csv`, `technical_manifest.json`, `manifest.json` under `<root>/pseudobulk/`).

Every manifest read SHALL record the cell type the step runs for; a manifest recording another cell type SHALL stop the step with an error naming the file and both cell types, and a manifest without a cell type SHALL read as `astrocyte` with a finding. It SHALL stop with an error naming the missing file and the command that produces it, and write nothing.

#### Scenario: Technical table absent
- **WHEN** `results/pseudobulk/donor_technical.csv` is missing
- **THEN** the step stops with an error naming the file and the `technical-covariates` command, and writes nothing

#### Scenario: Primary analysis absent
- **WHEN** `results/motco_analysis/cohort_30/` does not exist
- **THEN** the step stops with an error naming the directory and the `analyze` command, and writes nothing

#### Scenario: Cell-type mismatch
- **WHEN** the step runs for `microglia` and the export manifest under the given input directory records `astrocyte`
- **THEN** the step stops with an error naming the manifest and both cell types, and writes nothing

#### Scenario: Microglia defaults
- **WHEN** the step runs for `microglia` without directory flags
- **THEN** it reads `results/microglia/motco_input/`, `results/microglia/motco_analysis/` and `results/microglia/pseudobulk/`, each manifest recording `microglia`

### Requirement: Output layout mirrors the primary analysis
Each arm SHALL be written as a complete analysis directory, `<root>/motco_latent_sensitivity/<arm>/motco_analysis/<cohort>/` under the cell type's results root (`results/motco_latent_sensitivity/` for `astrocyte`), with the files `analyze` writes. Its manifest SHALL record the five export hashes as the analysis inputs, the technical table's path and sha256 separately, and the cell type. Pointing `axis-diagnostics` (the same cell type, `--analysis-dir` at the arm's `motco_analysis/`) at an arm SHALL then pass that command's input-hash and cell-type checks. The comparison manifest SHALL record the cell type. A cohort SHALL be written under a temporary name and renamed when both arms and the comparison are complete. An existing cohort output SHALL be an error without `--force`.

#### Scenario: Diagnostics on an arm
- **WHEN** `axis-diagnostics --analysis-dir results/motco_latent_sensitivity/latent_adjusted/motco_analysis` runs with the primary export as input
- **THEN** its input-hash check passes, and it reports the arm's latent components against the predeclared variables

#### Scenario: Diagnostics on a microglia arm
- **WHEN** `axis-diagnostics --cell-type microglia --analysis-dir results/microglia/motco_latent_sensitivity/latent_adjusted/motco_analysis` runs with the microglia export as input
- **THEN** its input-hash and cell-type checks pass, and it reports the arm's latent components against the predeclared variables

#### Scenario: Existing output
- **WHEN** `results/motco_latent_sensitivity/comparison/cohort_30/` exists and `--force` is not given
- **THEN** the step stops with an error naming the directory, and writes nothing

### Requirement: Committed latent-adjustment report
For each cell type the arms are run on, the change that runs them SHALL commit `reports/motco-<cell type>-latent-adjustment-<date>.md`. It SHALL contain:
- per cohort, the comparison table and the mediator check, pasted from the generated files;
- the first latent axis's technical correlations for the primary and both arms, from `axis-diagnostics`;
- each `latent_adjusted` arm's `feature_r2` summary;
- provenance: the MOTCO commit, this repository's revision, seed, permutations, `n_jobs` and, for the astrocyte report, the pin-bump check.

Every statistic, p-value and rank in the report SHALL be checked against the records before committing.

#### Scenario: Report checked
- **WHEN** the report is committed
- **THEN** a check finds every statistic, p-value and rank of the two reference and four arm records in the report text

#### Scenario: Microglia report
- **WHEN** the arms have run on microglia in both cohorts
- **THEN** `reports/motco-microglia-latent-adjustment-<date>.md` is committed with the two comparison tables, both mediator checks, the first-axis correlations of the primary and both arms, the `feature_r2` summaries and the provenance, every value checked against the six microglia records

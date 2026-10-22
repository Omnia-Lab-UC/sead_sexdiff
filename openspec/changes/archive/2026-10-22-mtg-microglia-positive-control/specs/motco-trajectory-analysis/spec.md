# Spec Delta

## MODIFIED Requirements

### Requirement: Inputs are the exported MOTCO blocks only
For each requested cohort the analysis SHALL read `rna.csv`, `atac.csv`, `metadata.csv`, `covariates.csv` and `manifest.json` from the cell type's export root (default `results/motco_input/<cohort>/` for `astrocyte`, `results/microglia/motco_input/<cohort>/` for `microglia`) and nothing else. It SHALL stop with an error naming the file when one is missing, when the donor columns of the four tables differ, when the export manifest's donor count or block feature counts disagree with the tables, or when the export manifest's cell type differs from the requested one (naming both). The analysis manifest SHALL record the cell type.

#### Scenario: Export not yet run
- **WHEN** `results/motco_input/cohort_50/` does not exist
- **THEN** the analysis stops with an error naming the directory and the `export` command that produces it

#### Scenario: Tables disagree
- **WHEN** `covariates.csv` lists a donor that `metadata.csv` does not
- **THEN** the analysis stops with an error naming both files

#### Scenario: Cell-type mismatch
- **WHEN** `analyze --cell-type microglia --input-dir results/motco_input` is run against the astrocyte export
- **THEN** the analysis stops with an error naming `astrocyte` and `microglia`, and writes nothing

#### Scenario: Microglia defaults
- **WHEN** `analyze --cell-type microglia` runs without directory flags
- **THEN** it reads `results/microglia/motco_input/<cohort>/` and writes `results/microglia/motco_analysis/<cohort>/` with the primary analysis's settings, and nothing under `results/motco_analysis/` changes

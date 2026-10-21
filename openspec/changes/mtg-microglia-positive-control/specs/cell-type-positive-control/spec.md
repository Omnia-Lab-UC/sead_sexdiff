# Spec Delta

## Purpose

Runs the unchanged MTG pipeline on a second, registered cell type with a documented ADNC-stage response (microglia) and evaluates predeclared expectations about the stage signal, so that the astrocyte non-rejections can be told apart from a pipeline that does not resolve stage structure on real SEA-AD data.

## ADDED Requirements

### Requirement: Cell-type registry
The system SHALL hold a fixed registry of cell types, each with a canonical name, the per-release `Subclass` label(s) and `Supertype` labels that define the population, the results root and the preregistered cohort cells. It SHALL contain `astrocyte` (the current mapping, results root `results/`, the current preregistered cells) and `microglia` (RNA 2026-06-22 `Subclass` `Immune` and ATAC 2024-12-06 `Subclass` `Microglia-PVM`, restricted to the microglia / PVM `Supertype` labels confirmed by the spike; results root `results/microglia/`). Every MTG step SHALL take the cell type by name and SHALL refuse a name not in the registry. The mapping SHALL NOT be a runtime option.

#### Scenario: Default is astrocyte
- **WHEN** any MTG step runs without a cell type
- **THEN** it selects astrocytes with the current mapping and reads and writes the current directories, and its outputs are unchanged

#### Scenario: Unregistered cell type
- **WHEN** `pseudobulk --cell-type oligodendrocyte` is run
- **THEN** the step stops with an error listing the registered names, and writes nothing

#### Scenario: Microglia selection excludes other immune supertypes
- **WHEN** the RNA pseudobulk runs for `microglia` on an object whose `Immune` subclass holds lymphocyte or monocyte nuclei
- **THEN** only nuclei with a registered microglia / PVM `Supertype` are summed, and any `Immune` supertype absent from the registry's microglia list and from its recorded non-microglia list is an error naming the file and the label

### Requirement: Spike before any pseudobulk
Before the microglia pseudobulk runs, the change SHALL record in `spike.md`, from the donor objects' `obs` columns only: the `Subclass` × `Supertype` labels of the immune / microglia population per release, the microglia nuclei per donor and omic over all donors, and the resulting sex × stage cells at the ≥ 30 and ≥ 50 thresholds. The preregistered microglia cells of the registry SHALL be those counts.

#### Scenario: Spike recorded
- **WHEN** the microglia registry entry is filled in
- **THEN** `spike.md` lists the labels per release, the per-donor nucleus counts and the cohort cells they imply, and the registry's preregistered cells equal them

### Requirement: Cell type recorded and checked through the chain
Every manifest of the MTG steps SHALL record `cell_type`. `export` SHALL refuse a pseudobulk manifest whose cell type differs from the requested one; `analyze`, `axis-diagnostics` and `stage-de` SHALL refuse an export manifest whose cell type differs; `positive-control` SHALL refuse a control root whose manifests are not `microglia` or a reference root whose manifests are not `astrocyte`. Each refusal SHALL be an error naming both cell types and the file.

#### Scenario: Astrocyte export analysed as microglia
- **WHEN** `analyze --cell-type microglia --input-dir results/motco_input` is run
- **THEN** the step stops with an error naming `astrocyte` (the export manifest's) and `microglia` (the requested), and writes nothing

### Requirement: Directory defaults follow the cell type
The MTG steps SHALL default their directories to the cell type's results root (`<root>/pseudobulk/`, `<root>/motco_input/`, `<root>/motco_analysis/`, `<root>/motco_diagnostics/`, `<root>/stage_de/`), with explicit directory flags taking precedence, so that the microglia outputs mirror the astrocyte layout.

#### Scenario: Microglia defaults
- **WHEN** `export --cell-type microglia` runs without directory flags
- **THEN** it reads `results/microglia/pseudobulk/` and writes `results/microglia/motco_input/<cohort>/`

### Requirement: Predeclared expectations evaluated from the generated files
The `positive-control` step SHALL evaluate, per cohort, from the generated files of both cell types and nothing else:
- E1: the number of exported RNA features with a stage effect at FDR 0.05 (`stage-de`, `factor_additive`, either coefficient, exported set) is larger for microglia than for astrocytes;
- E2: the microglia primary analysis's CV mean AUROC exceeds the astrocyte primary analysis's;
- E3: along the first latent component of the microglia primary analysis, |Spearman| with stage exceeds |Spearman| with every technical variable of the diagnostic's variable set.

The control SHALL hold when E1, E2 and E3 all hold in the primary cohort. The step SHALL read no p-value, rejection status or trajectory statistic, SHALL compute no test, and SHALL write the values compared, each expectation's outcome and the overall outcome to `positive_control.csv`, `positive_control.md` and `positive_control_manifest.json` with the sha256 of every file read. The expectations SHALL be constants, written before the microglia runs, and SHALL NOT be changed from results.

#### Scenario: Control holds
- **WHEN** microglia `cohort_30` has 1,840 exported stage-significant genes against the astrocytes' 310, CV AUROC 0.81 against 0.6313, and first-axis |ρ| 0.72 with stage against a largest technical |ρ| of 0.41
- **THEN** `positive_control.csv` records E1, E2 and E3 as holding and the overall outcome as holding, with all six compared values

#### Scenario: One expectation fails
- **WHEN** E1 and E2 hold but the largest technical |ρ| exceeds the stage |ρ| along the first microglia axis
- **THEN** the overall outcome is recorded as not holding, with E3 named, and nothing else changes

#### Scenario: Missing input
- **WHEN** `results/microglia/stage_de/cohort_30/stage_de_summary.csv` is absent
- **THEN** the step stops with an error naming the file and the `stage-de` command, and writes nothing

### Requirement: Microglia trajectory analysis under the existing rules
The microglia sex × stage analysis SHALL run through the existing `analyze` step on the microglia export with the primary analysis's settings (six analysis covariates, seed, permutations, orientation, CV rank selection, `n_jobs = 1` for reported results) and SHALL be read under the existing rules: a `delta` rejection is a difference, not specifically of size; `angle` and `shape` non-rejections are uninformative; the CV-selected rank beside `angle`; no Type I or power claim of its own. Nothing found in microglia SHALL change any astrocyte output or reading.

#### Scenario: Astrocyte outputs untouched
- **WHEN** the whole microglia chain has run
- **THEN** the sha256 of every file under `results/pseudobulk/`, `results/motco_input/`, `results/motco_analysis/` and `results/stage_de/` equals its value before

### Requirement: Committed report
The change SHALL commit `reports/motco-microglia-positive-control-<date>.md` with numbers pasted from the generated files and checked against them: provenance, the spike, the cohort tables, the primary analysis summary, the stage-DE counts of both cell types side by side, the E1–E3 table, the first-axis correlations, and the reading rules: the control is about the stage signal, not the sex × stage hypothesis; if it holds, the astrocyte null is not a pipeline failure; if not, the latent construction is in question for any cell type; the microglia trajectory result is read under the binding rules; MOTCO's operating characteristics were measured on neither cell type's blocks.

#### Scenario: Report numbers checked
- **WHEN** the report is committed
- **THEN** every value in its E1–E3 table equals `positive_control.csv`, and every statistic of the microglia analysis equals its `result.json`

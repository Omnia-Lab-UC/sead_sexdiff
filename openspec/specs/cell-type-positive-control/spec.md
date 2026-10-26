# cell-type-positive-control Specification

## Purpose

Runs the unchanged MTG pipeline on a second, registered cell type with a documented ADNC-stage response (microglia) and evaluates predeclared expectations about the stage signal, so that the astrocyte non-rejections can be told apart from a pipeline that does not resolve stage structure on real SEA-AD data.

## Requirements

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
Every manifest of the MTG steps SHALL record `cell_type`. `export` SHALL refuse a pseudobulk manifest whose cell type differs from the requested one; `analyze`, `axis-diagnostics`, `stage-de` and `stage-de-limma` SHALL refuse an export manifest whose cell type differs; `stage-de-limma` SHALL also refuse a DESeq2 stage-DE manifest whose cell type differs; `positive-control` SHALL refuse a control root whose manifests are not `microglia` or a reference root whose manifests are not `astrocyte`. Each refusal SHALL be an error naming both cell types and the file.

#### Scenario: Astrocyte export analysed as microglia
- **WHEN** `analyze --cell-type microglia --input-dir results/motco_input` is run
- **THEN** the step stops with an error naming `astrocyte` (the export manifest's) and `microglia` (the requested), and writes nothing

#### Scenario: Astrocyte DESeq2 outputs read for microglia
- **WHEN** `stage-de-limma --cell-type microglia --deseq2-dir results/stage_de` is run
- **THEN** the step stops with an error naming `astrocyte` (the DESeq2 manifest's) and `microglia` (the requested), and writes nothing

### Requirement: Directory defaults follow the cell type
The MTG steps SHALL default their directories to the cell type's results root (`<root>/pseudobulk/`, `<root>/motco_input/`, `<root>/motco_analysis/`, `<root>/motco_diagnostics/`, `<root>/stage_de/`, `<root>/stage_de_limma/`; `stage-de-limma` reads the same cell type's `<root>/stage_de/`), with explicit directory flags taking precedence, so that the microglia outputs mirror the astrocyte layout.

#### Scenario: Microglia defaults
- **WHEN** `export --cell-type microglia` runs without directory flags
- **THEN** it reads `results/microglia/pseudobulk/` and writes `results/microglia/motco_input/<cohort>/`

### Requirement: Predeclared expectations evaluated from the generated files
The `positive-control` step SHALL evaluate, per cohort, from the generated files of both cell types and nothing else:
- E1: under each stage-DE method, the number of exported RNA features with a stage effect at FDR 0.05 is larger for microglia than for astrocytes: `stage-de` (`factor_additive`, the union of its two coefficients, exported set) and `stage-de-limma` (`factor_additive`, the joint moderated F `stage_joint`, exported set). E1 holds only when it holds under both methods; when it holds under one only, it is recorded as `method_dependent` and counts as not holding;
- E2: the microglia primary analysis's CV mean AUROC exceeds the astrocyte primary analysis's;
- E3: along the first latent component of the microglia primary analysis, |Spearman| with stage exceeds |Spearman| with every technical variable of the diagnostic's variable set.

The control SHALL hold when E1, E2 and E3 all hold in the primary cohort. The step SHALL read no p-value, rejection status or trajectory statistic, SHALL compute no test, and SHALL write the values compared, each expectation's outcome and the overall outcome to `positive_control.csv`, `positive_control.md` and `positive_control_manifest.json` with the sha256 of every file read. The expectations SHALL be constants, written before the microglia runs, and SHALL NOT be changed from results.

#### Scenario: Control holds
- **WHEN** microglia `cohort_30` has 1,840 (DESeq2 union) and 1,610 (limma `stage_joint`) exported stage-significant genes against the astrocytes' 1 and 1, CV AUROC 0.81 against 0.6313, and first-axis |ρ| 0.72 with stage against a largest technical |ρ| of 0.41
- **THEN** `positive_control.csv` records E1 (under both methods), E2 and E3 as holding and the overall outcome as holding, with all eight compared values

#### Scenario: One expectation fails
- **WHEN** E1 and E2 hold but the largest technical |ρ| exceeds the stage |ρ| along the first microglia axis
- **THEN** the overall outcome is recorded as not holding, with E3 named, and nothing else changes

#### Scenario: Methods disagree on E1
- **WHEN** microglia `cohort_30` has more exported stage-significant genes than the astrocytes under DESeq2's union but not under limma's `stage_joint`
- **THEN** E1 is recorded as `method_dependent`, naming DESeq2 as the method under which it held, the overall outcome as not holding, and nothing else changes

#### Scenario: Missing input
- **WHEN** `results/microglia/stage_de/cohort_30/stage_de_summary.csv` or `results/microglia/stage_de_limma/cohort_30/stage_de_limma_summary.csv` is absent
- **THEN** the step stops with an error naming the file and the `stage-de` or `stage-de-limma` command, and writes nothing

### Requirement: Microglia trajectory analysis under the existing rules
The microglia sex × stage analysis SHALL run through the existing `analyze` step on the microglia export with the primary analysis's settings (six analysis covariates, seed, permutations, orientation, CV rank selection, `n_jobs = 1` for reported results) and SHALL be read under the existing rules: a `delta` rejection is a difference, not specifically of size; `angle` and `shape` non-rejections are uninformative; the CV-selected rank beside `angle`; no Type I or power claim of its own. Nothing found in microglia SHALL change any astrocyte output or reading.

#### Scenario: Astrocyte outputs untouched
- **WHEN** the whole microglia chain has run
- **THEN** the sha256 of every file under `results/pseudobulk/`, `results/motco_input/`, `results/motco_analysis/`, `results/stage_de/` and `results/stage_de_limma/` equals its value before

### Requirement: Both stage-DE methods on microglia
The change SHALL run `stage-de` and then `stage-de-limma` on the microglia export with their registered models, constants and settings unchanged, and SHALL report both concordance cells per cohort and omic: the DESeq2 cell (`STAGE_DE_CONCORDANCE_RULE`), the reported one, and the limma cell (`STAGE_DE_LIMMA_CONCORDANCE_RULE`) beside it with `methods_agree`. A disagreement SHALL be read as method dependence, a result of neither method, and SHALL change no MOTCO reading.

#### Scenario: Both cells reported
- **WHEN** both stage-DE steps have run on microglia `cohort_30`
- **THEN** the report lists, per omic, the DESeq2 cell, the limma cell, the `delta` status and `methods_agree`, each equal to the two steps' summaries

### Requirement: Committed report
The change SHALL commit `reports/motco-microglia-positive-control-<date>.md` with numbers pasted from the generated files and checked against them: provenance, the spike, the cohort tables, the primary analysis summary, the stage-DE counts of both cell types and both methods side by side, both concordance cells, the E1–E3 table, the first-axis correlations, and the reading rules: the control is about the stage signal, not the sex × stage hypothesis; if it holds, the astrocyte null is not a pipeline failure; if not, the latent construction is in question for any cell type; the microglia trajectory result is read under the binding rules; MOTCO's operating characteristics were measured on neither cell type's blocks.

#### Scenario: Report numbers checked
- **WHEN** the report is committed
- **THEN** every value in its E1–E3 table equals `positive_control.csv`, every stage-DE count and concordance cell equals the two methods' summaries, and every statistic of the microglia analysis equals its `result.json`

### Requirement: Expectations evaluated on a registered basis
The `positive-control` step SHALL evaluate E1–E3 on a basis of a fixed registry, `primary` (the default) or `latent_adjusted`. On `latent_adjusted`, E2 and E3 SHALL read each cell type's `latent_adjusted` arm analysis and the diagnostic made from it instead of the primary analysis and its diagnostic, and E1 SHALL be read unchanged. The expectations, their comparisons, the technical variables of E3 and the all-three rule SHALL NOT change with the basis. The basis SHALL be a registry, never extended from results.

#### Scenario: Primary basis unchanged
- **WHEN** the step runs without a basis
- **THEN** it reads the primary analyses and diagnostics, writes to `results/microglia/positive_control/<cohort>/`, and its three outputs equal the outputs of the step before this change

#### Scenario: Adjusted basis
- **WHEN** the step runs with the basis `latent_adjusted` for `cohort_30`
- **THEN** E2 compares the two cell types' `latent_adjusted` arms' CV mean AUROC, E3 reads the first component of the microglia `latent_adjusted` arm from the diagnostic made on that arm against the same technical variables, E1 equals the primary-basis E1, and the outputs go to `results/microglia/positive_control/latent_adjusted/<cohort>/`

#### Scenario: Arm absent
- **WHEN** the basis is `latent_adjusted` and `results/microglia/motco_latent_sensitivity/latent_adjusted/motco_analysis/cohort_30/` or the diagnostic on it is missing
- **THEN** the step stops with an error naming the directory and the `latent-sensitivity --cell-type microglia` or `axis-diagnostics` command, and writes nothing

### Requirement: Arm chain checked on the adjusted basis
On the `latent_adjusted` basis the step SHALL check, for each cell type, that the arm's manifest records the cell type, the hashes of the primary export files now on disk and the primary analysis it was compared with (the analysis read by E1's stage-DE steps), and that the diagnostic read for E3 records the arm analysis it was made from. On any difference the step SHALL stop with an error naming the file, and write nothing.

#### Scenario: Arm made from another analysis
- **WHEN** the microglia `latent_adjusted` comparison manifest records a primary analysis whose `result.json` sha256 differs from the one on disk
- **THEN** the step stops with an error naming both files, and writes nothing

### Requirement: Restated control is read as a diagnosis of the primary control
The restated control on `latent_adjusted` SHALL be reported beside the primary-basis control, with the predeclared reading: if it holds in the primary cohort, the measured quality gradient explains the primary control's failure, without promoting the arm, changing the primary-basis outcome or changing the astrocyte reading; if it does not hold, the latent construction is in question beyond the quality gradient. E3 on the adjusted basis SHALL carry the mediator caveat, and a failure for that reason SHALL count as a failure.

#### Scenario: Restated control holds
- **WHEN** E1, E2 and E3 all hold on `latent_adjusted` in `cohort_30` and the primary-basis control did not hold
- **THEN** the report states that the quality gradient explains the primary failure, that MOTCO's power comparison still excludes the arm as a replacement, and that the primary-basis outcome and the astrocyte reading are unchanged

#### Scenario: Restated control fails
- **WHEN** E3 does not hold on `latent_adjusted` in `cohort_30`
- **THEN** the report states that the latent construction is in question beyond the quality gradient, names the candidates (composition, rank-1 selection, microglia ATAC depth) and the mediator caveat, and the expectation is not re-tuned

### Requirement: Extended-covariates basis equals the primary's
No restated control SHALL be evaluated for `extended_covariates`. The report SHALL state and check that, for each cell type and cohort, the `extended_covariates` arm's CV-selected rank and CV mean AUROC equal the primary analysis's, as the covariate-blind basis implies; a difference SHALL be recorded as a finding, not resolved.

#### Scenario: Equality checked
- **WHEN** the microglia report is written
- **THEN** it lists, per cohort, the primary and `extended_covariates` rank and CV AUROC of both cell types and states whether they are equal

### Requirement: Control outcome re-issued on the floored primary
When the primary analyses are re-issued under the rank floor, the control SHALL be evaluated again on both bases with E1–E3, their comparisons, technical variables and the all-three rule unchanged. The floored outcome SHALL be the primary-basis outcome from then on; the report SHALL give the superseded outcome beside it, per expectation and cohort, with the values compared.

#### Scenario: Outcome unchanged
- **WHEN** E1–E3 give the same outcomes on the floored primary as on the superseded one in `cohort_30`
- **THEN** the report states that the floor leaves the control's outcome as it was, with both sets of values

#### Scenario: Outcome changes through E2
- **WHEN** the floored microglia CV mean AUROC no longer exceeds the astrocytes' and the superseded one did
- **THEN** E2 is recorded as not holding on the floored primary, the overall outcome follows the all-three rule, the superseded outcome is reported beside it as dependence on the rank rule, and no expectation is changed

#### Scenario: First component reproduces
- **WHEN** E3 is evaluated on the floored primary
- **THEN** the report states whether the first-component correlations equal the superseded ones (PLS components nest across ranks), and a difference beyond rounding is recorded as a finding

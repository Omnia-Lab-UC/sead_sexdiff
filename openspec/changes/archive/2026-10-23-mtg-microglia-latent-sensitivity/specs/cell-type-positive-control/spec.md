# Spec Delta

## ADDED Requirements

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

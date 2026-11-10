## ADDED Requirements

### Requirement: Scan concordance cell beside the limma cell
For each cohort and omic the step SHALL compare `scan_detected` (at least one hit of the analysis read in that omic's block) with limma's `interaction_detected` on `interaction_joint`. The cell SHALL be `concordant` when both agree and `discordant` otherwise, a discordance written as a finding, never a result. It SHALL count beside it the hits whose `interaction_joint` adjusted p ≤ α. The rule SHALL be a constant, read no p-value magnitude and change no MOTCO output.

#### Scenario: Both negative
- **WHEN** no exported RNA feature of `cohort_30` has an `interaction_joint` adjusted p ≤ 0.05 and the primary's scan has no RNA hit
- **THEN** the limma summary's rule row for `cohort_30`, `rna` has `scan_detected = false` and `scan_concordance = concordant`, beside the unchanged `delta` cell

#### Scenario: Hit outside the tables
- **WHEN** a hit names a feature absent from the exported set of the limma tables
- **THEN** the step stops naming the feature and the analysis, and writes nothing

### Requirement: Scan cells in the method comparison
The `interaction_joint` row of `method_comparison.csv` SHALL carry `scan_detected` and the scan concordance cell of each method: the DESeq2 cell read from the checked DESeq2 summary, and the limma cell. A scan cell that differs between the methods SHALL be written as a method-dependence finding. Neither scan cell SHALL replace the other or the `delta` cells.

#### Scenario: Methods disagree on the scan cell
- **WHEN** the scan has one RNA hit in `cohort_30`, DESeq2's exported union count there is 1 and limma's `interaction_joint` count is 0
- **THEN** the row has `scan_concordance_deseq2 = concordant` and `scan_concordance_limma = discordant`, and the manifest records both a scan discordance finding and a method-dependence finding

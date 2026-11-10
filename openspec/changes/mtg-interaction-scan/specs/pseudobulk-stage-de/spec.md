## ADDED Requirements

### Requirement: Predeclared scan concordance cell against the primary analysis
For each cohort and omic the step SHALL evaluate, from the generated files, `scan_detected`: at least one hit of the analysis read in that omic's block, from its `interaction_scan` block. It SHALL compare that with the existing `interaction_detected`. The cell SHALL be `concordant` when both agree and `discordant` otherwise. A discordance SHALL be written as a finding, never a result. The rule SHALL be a constant, read no p-value magnitude and change no MOTCO output.

#### Scenario: Both negative
- **WHEN** no exported RNA feature of `cohort_30` has an interaction `padj ≤ 0.05` and the primary's scan has no RNA hit
- **THEN** the summary's rule row for `cohort_30`, `rna` has `scan_hits = 0`, `scan_detected = false` and `scan_concordance = concordant`, beside the unchanged `delta` cell

#### Scenario: Scan detects, DESeq2 does not
- **WHEN** the primary's scan has two ATAC hits in `cohort_30` and DESeq2's exported union count there is 0
- **THEN** the cell is `discordant`, and the manifest records a finding naming the cohort, the omic, both counts and the hits, with the hit overlap 0

### Requirement: Scan hit overlap counted beside the cell
Beside the scan cell the step SHALL count `scan_hits_significant`: the analysis's hits in that omic's block whose `padj ≤ α` on either interaction coefficient of `factor_interaction`. Every hit SHALL be found among the exported features of the table, otherwise the step SHALL stop naming the feature. The manifest SHALL record the hits read, the analysis `result.json` hash and the rule.

#### Scenario: Hits are exported features
- **WHEN** the primary's scan names an RNA hit
- **THEN** that feature is in the exported set of the `factor_interaction` tables, and the overlap count says whether DESeq2 also calls it

#### Scenario: Hit outside the tables
- **WHEN** a hit names a feature absent from the exported set of the tables
- **THEN** the step stops naming the feature and the analysis, and writes nothing

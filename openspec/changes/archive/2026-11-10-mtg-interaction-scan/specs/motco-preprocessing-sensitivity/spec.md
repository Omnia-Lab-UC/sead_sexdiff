## MODIFIED Requirements

### Requirement: Comparison table generated from the records
The step SHALL write `comparison/<cohort>/sensitivity.csv` and its Markdown rendering `sensitivity.md`, generated from the records and never typed by hand. The primary row comes first (`variant = primary`) and the variant row second. Each row SHALL give:
- the rank selection and the rank;
- per statistic (`delta`, `angle`, `shape`), the observed value, the RRPP p-value, whether p ≤ α and the null mean, sd and q95;
- the split columns of spec `motco-span-split-test`;
- the scan columns of spec `motco-interaction-scan`;
- the CV AUROC and the row's directory.

#### Scenario: Table rows
- **WHEN** the comparison of `cohort_30` is written
- **THEN** `sensitivity.csv` has two rows, the first with `variant = primary` and the primary's rank, the second with `variant = keep_chrx` and the variant's rank

#### Scenario: Table values equal the records
- **WHEN** a cell of `sensitivity.csv` is compared with the `result.json` of its row's directory
- **THEN** they are equal, floating-point cells being written in their shortest round-trip representation

#### Scenario: Runtime and directories
- **WHEN** the table of `cohort_30` is written
- **THEN** it has no runtime column, its variant row names `motco_analysis/cohort_30` relative to the variant root and its primary row names the primary analysis directory, while `sensitivity.md` and the manifest carry each run's runtime

#### Scenario: Variant carries the split test
- **WHEN** a variant analysis finishes
- **THEN** its record carries a `span_split` block made with the primary's permutations, seed and correction, and its row's split cells are filled

#### Scenario: Variant carries the scan
- **WHEN** a variant analysis finishes
- **THEN** its record carries an `interaction_scan` block over the variant's own exported features, made with the primary's permutations, seed and correction, and its row's scan cells are filled; a scan status that differs from the primary's is written as a preprocessing-dependence finding

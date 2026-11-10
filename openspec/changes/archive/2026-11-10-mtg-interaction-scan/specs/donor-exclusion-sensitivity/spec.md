## MODIFIED Requirements

### Requirement: Comparison table generated from the records
Per cohort the step SHALL write `exclusion_sensitivity.csv` with the primary row and the variant row. Each row gives:
- donors, cells and `k`;
- the selected rank and CV mean AUROC;
- for `delta`, `angle` and `shape`, the observed value, p-value, rejection and null mean, SD and q95;
- the split columns of spec `motco-span-split-test`;
- the scan columns of spec `motco-interaction-scan`.

The step SHALL also write the table's Markdown rendering and a manifest. The manifest records the variant, the excluded donors, every input path and sha256, seed, permutations, `n_jobs`, software, git revision and findings. Every value SHALL be read from the records, not computed again.

#### Scenario: Values traceable
- **WHEN** `exclusion_sensitivity.csv` is written
- **THEN** every statistic, p-value, split cell, scan cell and rank equals the corresponding value in the primary or variant `result.json`

#### Scenario: Scan status differs from the primary
- **WHEN** the variant's scan is rejected in a cohort and the primary's is not, or the reverse
- **THEN** the manifest records a finding naming the variant, the cohort and both p-values as cohort dependence of the scan, not a result of the variant

# Spec Delta

## MODIFIED Requirements

### Requirement: Comparison table generated from the records
Per cohort the step SHALL write `exclusion_sensitivity.csv` with the primary row and the variant row. Each row gives:
- donors, cells and `k`;
- the selected rank and CV mean AUROC;
- for `delta`, `angle` and `shape`, the observed value, p-value, rejection and null mean, SD and q95;
- the split columns of spec `motco-span-split-test`.

The step SHALL also write the table's Markdown rendering and a manifest. The manifest records the variant, the excluded donors, every input path and sha256, seed, permutations, `n_jobs`, software, git revision and findings. Every value SHALL be read from the records, not computed again.

#### Scenario: Values traceable
- **WHEN** `exclusion_sensitivity.csv` is written
- **THEN** every statistic, p-value, split cell and rank equals the corresponding value in the primary or variant `result.json`

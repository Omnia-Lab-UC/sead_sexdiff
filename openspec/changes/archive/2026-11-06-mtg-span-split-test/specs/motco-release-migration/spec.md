# Spec Delta

## MODIFIED Requirements

### Requirement: Measurement contract checked on the primaries
When MOTCO documents that the new setting changes only named outputs, or only adds a named record block, the primary of every cell type and cohort SHALL be run on the same machine under the new pin with and without it. Every output outside the named ones SHALL be equal; otherwise the change SHALL stop before any superseded output is moved.

#### Scenario: Correction moves only the null
- **WHEN** the four primaries are run under MOTCO 0.13.0 with and without `residual_correction = df`
- **THEN**:
  - `latent_scores.csv`, `pls_loadings.csv`, `stage_means.csv` and `analysis_metadata.csv` are byte-identical per cohort;
  - `result.json` is equal after dropping `p_values`, `null_summary`, `config_spectrum.permutation_pooled_eigengap` (drawn from the null), the correction's entries under `params` and `measurement`, and the runtime;
  - each κ equals √(n / (n − r)), with r the recorded reduced rank (10 for six covariates)

#### Scenario: Split test only adds a block
- **WHEN** the four primaries are run under MOTCO 0.18.0 with the `df` correction, with and without `span_split_test`
- **THEN**:
  - the four CSVs are byte-identical per cohort;
  - `result.json` is equal after dropping the `span_split` block, `params.span_split_test` and the runtime, including every p-value and null summary, because the split test draws from its own generator at the analysis's seed;
  - each block's permutations, seed and correction equal the analysis's

#### Scenario: Contract broken
- **WHEN** a latent score, the selected rank, the CV AUROC, an observed statistic or the observed spectrum differs between the two runs
- **THEN** nothing is moved or re-run, and the change stops naming the cohort and the differing value

### Requirement: Version consistency checked
After the re-issue, every manifest under `results/` outside `results/superseded/` that records a MOTCO version SHALL record the new pin's version. Every analysis record SHALL carry the new measurement setting, except where MOTCO does not allow the setting in that analysis and its manifest records it as not available. The figures manifest SHALL record the new pin's version for every analysis it drew. A manifest that does not SHALL be listed, and the report SHALL NOT be committed until it is resolved.

#### Scenario: Stale manifest found
- **WHEN** a diagnostic manifest under `results/motco_sensitivity/` still records MOTCO 0.8.0
- **THEN** the check lists it and the step that made it is rerun before the report is written

#### Scenario: Figures drawn before the last rerun
- **WHEN** `results/figures/figures_manifest.json` records MOTCO 0.13.0 for an analysis after the pin moved
- **THEN** the check lists it, the uncommitted frozen set is removed, and `figures --report` is run again before the report is written

#### Scenario: Setting not available in an analysis
- **WHEN** a `latent_adjusted` arm's record carries no `span_split` block after the move to 0.18.0
- **THEN** the check passes it only if its manifest records the split test as not available, and lists it otherwise; any other record without the block is listed

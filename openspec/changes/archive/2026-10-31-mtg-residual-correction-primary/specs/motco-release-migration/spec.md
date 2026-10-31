## ADDED Requirements

### Requirement: MOTCO pinned to a published release
The MOTCO dependency SHALL be an exact version from PyPI, with the wheel's sha256 in `uv.lock`. Before any run under a new pin, the wheel's `motco/` package SHALL be byte-identical to `src/motco/` at that version's release tag, and no other locked package's version SHALL change. Both checks SHALL be recorded in the re-issue report.

#### Scenario: Wheel matches the tag
- **WHEN** the pin moves to `motco==0.13.0`
- **THEN** the downloaded wheel's `motco/` tree equals `src/motco/` at tag `v0.13.0` file by file, `uv.lock` records motco 0.13.0 from the PyPI registry with the wheel's sha256, and the report names the wheel's sha256 and the tag's commit

#### Scenario: Wheel differs from the tag
- **WHEN** any file of the wheel's `motco/` differs from the tag's `src/motco/`
- **THEN** no analysis runs under the new pin and the change stops naming the file

#### Scenario: Another package moves
- **WHEN** re-locking with the new pin changes the version of any package other than motco
- **THEN** the pin move is not committed until the change is reverted or proposed separately

### Requirement: Measurement contract checked on the primaries
When MOTCO documents that the new setting changes only named outputs, the primary of every cell type and cohort SHALL be run on the same machine under the new pin with and without it. Every output outside the named ones SHALL be equal; otherwise the change SHALL stop before any superseded output is moved.

#### Scenario: Correction moves only the null
- **WHEN** the four primaries are run under MOTCO 0.13.0 with and without `residual_correction = df`
- **THEN** `latent_scores.csv`, `pls_loadings.csv`, `stage_means.csv` and `analysis_metadata.csv` are byte-identical per cohort, `result.json` is equal after dropping `p_values`, `null_summary`, `config_spectrum.permutation_pooled_eigengap` (drawn from the null), the correction's entries under `params` and `measurement`, and the runtime, and each κ equals √(n / (n − r)) with r the recorded reduced rank (10 for six covariates)

#### Scenario: Contract broken
- **WHEN** a latent score, the selected rank, the CV AUROC, an observed statistic or the observed spectrum differs between the two runs
- **THEN** nothing is moved or re-run, and the change stops naming the cohort and the differing value

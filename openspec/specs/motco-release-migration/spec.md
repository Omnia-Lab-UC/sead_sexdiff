# motco-release-migration Specification

## Purpose

Governs how a MOTCO pin move that changes the measurement re-issues the case study, so that every result rests on one MOTCO release, the earlier results stay resolvable as superseded records, and the new primary is read under a reading written before it was run.

## Requirements

### Requirement: Same-machine gate before the move
Before any run under the new pin's changed measurement, the primary analysis of every cell type and cohort SHALL be run under the new pin with the previous measurement restored and under the previous pin from the pre-change commit, on the same machine and inputs. The CSVs SHALL be byte-identical and the JSON records equal apart from runtimes, times, revisions, versions and output paths; otherwise the change SHALL stop.

#### Scenario: Gate passes
- **WHEN** MOTCO 0.10.0 with the rank floor set to 1 and MOTCO 0.9.0 from the pre-change commit analyse the same four exports on `nixos-home`
- **THEN** `latent_scores.csv`, `pls_loadings.csv`, `stage_means.csv` and `analysis_metadata.csv` are byte-identical per cohort, the records are equal after dropping the listed keys, and the floored reruns may start

#### Scenario: Gate fails
- **WHEN** any compared file differs
- **THEN** no floored analysis is run, the superseded outputs stay in place, and the change stops with the differing file named

### Requirement: Comparison with the kept outputs recorded
The gate's new-pin rerun SHALL also be compared with the kept primary directories. The outcome SHALL be recorded in the report (byte-identical, or the largest differences and whether p-values, ranks and metadata are equal) and SHALL NOT be a gate, because analyses made on another machine differ at rounding level.

#### Scenario: Kept primary from another machine
- **WHEN** the kept astrocyte primary was made on the cluster and differs from the workstation rerun only at rounding level
- **THEN** the report records the largest absolute differences and the equal p-values and ranks, and the change continues

### Requirement: Superseded outputs moved, not deleted
Before the reruns, every MOTCO-dependent output directory of every cell type SHALL be moved to `results/superseded/<change>/` with its path relative to `results/` kept, on the workstation that runs the re-issue and on the cluster (the canonical copy). Every moved file's sha256 SHALL equal its value before the move. Nothing SHALL be deleted.

#### Scenario: Astrocyte and microglia analyses moved
- **WHEN** the move runs for this change
- **THEN** `results/motco_analysis/` is at `results/superseded/mtg-rank-floor-primary/motco_analysis/` and `results/microglia/motco_analysis/` at `results/superseded/mtg-rank-floor-primary/microglia/motco_analysis/`, each file hashing as before

#### Scenario: Hash cited by an earlier report
- **WHEN** an earlier report cites the sha256 of a superseded `result.json`
- **THEN** a file with that sha256 exists under `results/superseded/mtg-rank-floor-primary/`

### Requirement: MOTCO-independent outputs untouched
The pseudobulk, covariate, technical and flag tables, the exports and the annotation SHALL NOT be moved, rewritten or rerun. Their sha256 SHALL be recorded before the move and SHALL be equal after the last rerun.

#### Scenario: Exports after the re-issue
- **WHEN** every rerun has finished
- **THEN** every file under `results/pseudobulk/`, `results/motco_input/`, `results/microglia/pseudobulk/` and `results/microglia/motco_input/` hashes as before the move

### Requirement: Every MOTCO-dependent output re-issued
Every step whose superseded output was moved SHALL be rerun from scratch under the new pin for every cell type it ran on, in dependency order, with the commands of the reports that made it and `n_jobs = 1`, so that no output depends on the superseded ones. `figures` SHALL be the last step, run with `--report` for the re-issue report after every other rerun has finished. `results/figures/` is a working copy and SHALL NOT be moved.

#### Scenario: Microglia chain re-issued
- **WHEN** the re-issue finishes
- **THEN** the microglia `analyze`, `axis-diagnostics`, `latent-sensitivity`, both stage-DE steps, `exclusion-sensitivity` with its diagnostics and stage DE, and `positive-control` on both bases each have outputs made under the new pin

#### Scenario: Step added since the superseded run
- **WHEN** a superseded output came from an astrocyte-only step (`rank-sensitivity`, `preprocessing-sensitivity`)
- **THEN** it is rerun for astrocytes only, with its registry unchanged

#### Scenario: Figures redrawn last
- **WHEN** every other rerun has finished
- **THEN** `figures --report <re-issue report stem>` redraws the set from the re-issued outputs with its consistency check passing, and freezes it with the re-issue report. Figure 01, which reads no MOTCO-dependent file, is byte-identical to its copy in the previous report's frozen set when the module and matplotlib are unchanged since then

#### Scenario: Superseded figures
- **WHEN** a reader wants the figures of analyses superseded by a re-issue made after `mtg-case-study-figures`
- **THEN** they are the frozen set of the report that presented those analyses, under `reports/figures/<that report's stem>/`, and nothing under `results/superseded/` holds figures

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

#### Scenario: Setting available everywhere
- **WHEN** the check runs after the move to 0.20.0
- **THEN** every analysis record, the four `latent_adjusted` arms included, must carry an `interaction_scan` block and an `interaction_scan.csv` listed in its manifest's `outputs`, and any record without them is listed

### Requirement: Re-issue report with superseded values beside
The change SHALL commit one report that gives the new primary analyses in full, each re-issued step's regenerated tables, and per step an outcome table of superseded against new values. Numbers SHALL be pasted from generated files and checked by a script. Earlier reports SHALL NOT be edited.

#### Scenario: Outcome table
- **WHEN** the report presents the latent-adjustment sensitivity
- **THEN** it lists per cell type, cohort and arm the superseded and the new selected rank, CV AUROC, statistics, p-values and rejections, each equal to the corresponding record

### Requirement: Reading predeclared before the reruns
The reading of the new primary and of every re-issued sensitivity SHALL be written in the change's proposal before any rerun. A rejection status that differs between a new analysis and its superseded one SHALL be reported as dependence on the measurement change, with the new status as the result under the binding rules. No registry, expectation, adoption rule or threshold SHALL change.

#### Scenario: Status differs after the move
- **WHEN** the superseded primary did not reject `delta` and the new primary does
- **THEN** the report gives the new rejection as the primary result, read as a difference not specifically of size, and states the superseded status beside it as dependence on the rank rule

#### Scenario: Expected reproduction fails
- **WHEN** a value the construction says should reproduce (the first latent component, a stage-DE table) differs from its superseded value beyond rounding
- **THEN** the report records it as a finding with both values, and no reading is changed because of it

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
When MOTCO documents that the new setting changes only named outputs, or only adds a named record block (with or without named output files), the primary of every cell type and cohort SHALL be run on the same machine under the new pin with and without it. Every output outside the named ones SHALL be equal; otherwise the change SHALL stop before any superseded output is moved.

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

#### Scenario: Scan only adds a block and a table
- **WHEN** the four primaries are run under MOTCO 0.20.0 with the `df` correction and the split test, with and without `interaction_scan`
- **THEN**:
  - the four CSVs are byte-identical per cohort, and only the run with the scan writes `interaction_scan.csv`;
  - `result.json` is equal after dropping the `interaction_scan` block, `params.interaction_scan` and the runtime, including every p-value, null summary and the `span_split` block;
  - each scan block's permutations, seed and correction equal the analysis's

#### Scenario: Contract broken
- **WHEN** a latent score, the selected rank, the CV AUROC, an observed statistic or the observed spectrum differs between the two runs
- **THEN** nothing is moved or re-run, and the change stops naming the cohort and the differing value

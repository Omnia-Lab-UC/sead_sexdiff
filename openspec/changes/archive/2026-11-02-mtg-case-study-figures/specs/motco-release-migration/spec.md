# Spec Delta

## MODIFIED Requirements

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
After the re-issue, every manifest under `results/` outside `results/superseded/` that records a MOTCO version SHALL record the new pin's version, and every analysis record SHALL carry the new measurement setting. The figures manifest SHALL record the new pin's version for every analysis it drew. A manifest that does not SHALL be listed and the report SHALL NOT be committed until it is resolved.

#### Scenario: Stale manifest found
- **WHEN** a diagnostic manifest under `results/motco_sensitivity/` still records MOTCO 0.8.0
- **THEN** the check lists it and the step that made it is rerun before the report is written

#### Scenario: Figures drawn before the last rerun
- **WHEN** `results/figures/figures_manifest.json` records MOTCO 0.13.0 for an analysis after the pin moved
- **THEN** the check lists it, the uncommitted frozen set is removed, and `figures --report` is run again before the report is written

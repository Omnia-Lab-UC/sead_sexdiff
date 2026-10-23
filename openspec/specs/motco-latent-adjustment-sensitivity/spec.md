# motco-latent-adjustment-sensitivity Specification

## Purpose

Measures how the covariate-adjusted sex × ADNC-stage trajectory comparison of the astrocyte cohorts depends on the covariate set and on how the latent basis is chosen. It reruns MOTCO's analysis entry point on the primary export twice beside the primary analysis: with a predeclared extended covariate set, and with the same set plus MOTCO's `stage`-conditioned latent adjustment. The arms are recorded so that they are read as sensitivity to the primary analysis, never as its replacement.

## Requirements

### Requirement: Inputs are the primary export, the primary analysis and the technical table
For each requested cohort the step SHALL read only, under the cell type's results root (`results/` for `astrocyte`, `results/<name>/` otherwise):
- the primary export (`rna.csv`, `atac.csv`, `metadata.csv`, `covariates.csv`, `manifest.json` under `<root>/motco_input/<cohort>/`);
- the primary analysis (`result.json`, `analysis_manifest.json` under `<root>/motco_analysis/<cohort>/`);
- the technical table with its manifest and the pseudobulk manifest (`donor_technical.csv`, `technical_manifest.json`, `manifest.json` under `<root>/pseudobulk/`).

Every manifest read SHALL record the cell type the step runs for; a manifest recording another cell type SHALL stop the step with an error naming the file and both cell types, and a manifest without a cell type SHALL read as `astrocyte` with a finding. It SHALL stop with an error naming the missing file and the command that produces it, and write nothing.

#### Scenario: Technical table absent
- **WHEN** `results/pseudobulk/donor_technical.csv` is missing
- **THEN** the step stops with an error naming the file and the `technical-covariates` command, and writes nothing

#### Scenario: Primary analysis absent
- **WHEN** `results/motco_analysis/cohort_30/` does not exist
- **THEN** the step stops with an error naming the directory and the `analyze` command, and writes nothing

#### Scenario: Cell-type mismatch
- **WHEN** the step runs for `microglia` and the export manifest under the given input directory records `astrocyte`
- **THEN** the step stops with an error naming the manifest and both cell types, and writes nothing

#### Scenario: Microglia defaults
- **WHEN** the step runs for `microglia` without directory flags
- **THEN** it reads `results/microglia/motco_input/`, `results/microglia/motco_analysis/` and `results/microglia/pseudobulk/`, each manifest recording `microglia`

### Requirement: References must match
Before any arm runs, the step SHALL check that:
- the primary analysis manifest records, for each of the five export files, the sha256 of the file now on disk;
- it records the installed MOTCO version and the step's seed and permutation count;
- the technical manifest records the sha256 of the pseudobulk manifest now on disk and gene sets equal to `TECHNICAL_GENE_SETS`;
- the technical table has a finite value of every added covariate for every cohort donor.

On any difference the step SHALL stop with an error naming the file, setting, donor or column, and write nothing.

#### Scenario: Primary analysis made under another MOTCO
- **WHEN** `analysis_manifest.json` records MOTCO 0.8.0 and 0.9.0 is installed
- **THEN** the step stops with an error naming both versions and the `analyze` command, and writes nothing

#### Scenario: Fraction missing for a donor
- **WHEN** a cohort donor has an empty `mitochondrial_pseudobulk` cell
- **THEN** the step stops with an error naming the donor and the column, and writes nothing

### Requirement: MOTCO release with latent adjustment
The step SHALL run through MOTCO's trajectory-analysis entry point of a pinned release ≥ 0.9.0. It SHALL stop with an error naming the installed and required versions before reading any input when the installed MOTCO is older.

#### Scenario: Older MOTCO installed
- **WHEN** MOTCO 0.8.0 is installed
- **THEN** the step stops naming 0.8.0 and 0.9.0, before reading any input

### Requirement: Predeclared covariate set
Both arms SHALL declare, additively and in this order, the six analysis covariates of the primary analysis followed by `ribosomal_pseudobulk` and `mitochondrial_pseudobulk`. The fractions SHALL be taken from the technical table on the raw fraction scale, joined to the export by `donor_id`. The missing-covariate policy SHALL be `error`, and no value SHALL be imputed. The set SHALL be a fixed constant, not a runtime option, and SHALL NOT be chosen or changed from results.

#### Scenario: Record lists the eight covariates
- **WHEN** an arm finishes
- **THEN** its `result.json` covariate block lists `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` in that order

### Requirement: Predeclared arms
The step SHALL run exactly two arms per cohort, each once:
- `extended_covariates`, with the covariate-blind basis (no latent adjustment);
- `latent_adjusted`, with MOTCO's covariate latent adjustment conditioned on `stage`.

The arms SHALL be a fixed registry, not a runtime option. `group_stage` conditioning SHALL NOT be run. Apart from the covariate set and the latent-adjustment setting, both arms SHALL use the primary analysis's settings: the same export, blocks and units, group and stage columns, stage order, orientation estimator, CV rank selection, seed and permutation count.

#### Scenario: Arms complete
- **WHEN** the step finishes for `cohort_30`
- **THEN** two arm records exist: `extended_covariates` without a `latent_adjustment` block, and `latent_adjusted` with one whose `conditioning` is `stage`; both have `component_selection = cv`, the primary's seed and 999 permutations

#### Scenario: Latent adjustment changes only the basis
- **WHEN** both arms of a cohort have finished
- **THEN** the two records carry the same covariate block, and any difference between their statistics arises through the latent basis and the CV-selected rank only, as MOTCO's latent-adjustment contract guarantees

### Requirement: Primary analysis is the reference and is never modified
The primary analysis SHALL be the reference row of the comparison, taken from its `result.json` unchanged. The step SHALL NOT rerun it, and SHALL NOT write under `results/motco_analysis/`, `results/motco_input/` or `results/pseudobulk/`.

#### Scenario: Reference row
- **WHEN** the comparison of a cohort is written
- **THEN** its reference row carries the primary record's statistics, p-values, null summaries and CV-selected rank, and the primary files' sha256 are unchanged

### Requirement: Predeclared mediator check
For each cohort, before any arm runs, the step SHALL write `mediator_check.csv`. It gives, for each of the two added covariates, the Spearman correlation with stage (ordinal 0/1/2) over all cohort donors and within each sex, with the number of donors. A cell with `|ρ| ≥ DIAGNOSTIC_CORRELATION_THRESHOLD` SHALL be flagged as a mediator candidate. The check SHALL compute no p-value. A flag SHALL NOT remove a covariate, change an arm or stop a run.

#### Scenario: Fraction associated with stage
- **WHEN** `mitochondrial_pseudobulk` has Spearman 0.35 with stage in `cohort_30`
- **THEN** that row is flagged as a mediator candidate, both arms still run with the eight covariates, and the comparison manifest lists the flag

### Requirement: Output layout mirrors the primary analysis
Each arm SHALL be written as a complete analysis directory, `<root>/motco_latent_sensitivity/<arm>/motco_analysis/<cohort>/` under the cell type's results root (`results/motco_latent_sensitivity/` for `astrocyte`), with the files `analyze` writes. Its manifest SHALL record the five export hashes as the analysis inputs, the technical table's path and sha256 separately, and the cell type. Pointing `axis-diagnostics` (the same cell type, `--analysis-dir` at the arm's `motco_analysis/`) at an arm SHALL then pass that command's input-hash and cell-type checks. The comparison manifest SHALL record the cell type. A cohort SHALL be written under a temporary name and renamed when both arms and the comparison are complete. An existing cohort output SHALL be an error without `--force`.

#### Scenario: Diagnostics on an arm
- **WHEN** `axis-diagnostics --analysis-dir results/motco_latent_sensitivity/latent_adjusted/motco_analysis` runs with the primary export as input
- **THEN** its input-hash check passes, and it reports the arm's latent components against the predeclared variables

#### Scenario: Diagnostics on a microglia arm
- **WHEN** `axis-diagnostics --cell-type microglia --analysis-dir results/microglia/motco_latent_sensitivity/latent_adjusted/motco_analysis` runs with the microglia export as input
- **THEN** its input-hash and cell-type checks pass, and it reports the arm's latent components against the predeclared variables

#### Scenario: Existing output
- **WHEN** `results/motco_latent_sensitivity/comparison/cohort_30/` exists and `--force` is not given
- **THEN** the step stops with an error naming the directory, and writes nothing

### Requirement: Comparison table generated from the records
Per cohort the step SHALL write `comparison/<cohort>/latent_sensitivity.csv`. It holds the reference row and then one row per arm, each with:
- arm name, covariate names, conditioning (or none);
- selected rank, CV mean AUROC;
- for `delta`, `angle` and `shape`: observed value, p-value, rejection under the primary's α and rejection rule, and null mean, SD and q95;
- for `latent_adjusted`: the latent-adjustment block's `residual_df`, `cv_min_training_df` and `feature_r2` summary.

The step SHALL also write its Markdown rendering, `latent_sensitivity.md`, and `latent_sensitivity_manifest.json`. The manifest records the covariate set, the arms, seed, permutations, `n_jobs`, every input path and sha256, the reference analysis, the arm directories, the mediator flags, software, git revision and findings. Every value SHALL be read from the records, not computed again.

#### Scenario: Values traceable
- **WHEN** `latent_sensitivity.csv` is written
- **THEN** every statistic, p-value and rank in it equals the corresponding value in the reference or arm `result.json`

### Requirement: Reproducible run
Two runs on the same inputs, MOTCO release, seed and `n_jobs` SHALL give identical arm statistics, p-values, selected ranks and null summaries, a byte-identical `latent_sensitivity.csv` and `mediator_check.csv`, and byte-identical arm CSVs. The manifest SHALL record `n_jobs`. A run with `n_jobs` other than 1 SHALL carry a finding, and the report SHALL use `n_jobs = 1`.

#### Scenario: Rerun with force
- **WHEN** a cohort is run twice with `--force` and `n_jobs = 1`
- **THEN** both `latent_sensitivity.csv` files are byte-identical

### Requirement: Primary reproduced under the new MOTCO pin
Before the arms are run on real data, the primary analysis SHALL be rerun under the new MOTCO pin and compared with a rerun under the previous pin (MOTCO 0.8.0, the pre-change code in its own locked environment) on the same machine and inputs. Between the two reruns `latent_scores.csv`, `pls_loadings.csv`, `stage_means.csv` and `analysis_metadata.csv` SHALL be byte-identical, and `result.json` SHALL be equal apart from the runtime and the MOTCO version. On any other difference the kept 0.8.0 directories SHALL be restored and the change SHALL stop.

The new-pin rerun SHALL also be compared with the kept 0.8.0 directories. The outcome (byte-identical, or the largest differences found and whether p-values and ranks are equal) SHALL be recorded in the report's provenance and SHALL NOT stop the change. When the same-machine comparison passes, the new-pin rerun SHALL be the reference of the arms, and the kept directories SHALL be kept until the report is committed.

#### Scenario: Pin bump verified
- **WHEN** `analyze --force` is rerun for both cohorts under MOTCO 0.9.0, and the pre-change code under MOTCO 0.8.0 on the same machine
- **THEN** the same-machine comparison passes, the 0.9.0 rerun is the reference of the arms, and both comparisons are recorded in the report's provenance

#### Scenario: Rerun on another machine
- **WHEN** the kept 0.8.0 directories were made on another machine and the 0.9.0 rerun differs from them at rounding level
- **THEN** the change continues on the same-machine comparison, and the report records the largest differences and that p-values and ranks are equal

### Requirement: Committed latent-adjustment report
For each cell type the arms are run on, the change that runs them SHALL commit `reports/motco-<cell type>-latent-adjustment-<date>.md`. It SHALL contain:
- per cohort, the comparison table and the mediator check, pasted from the generated files;
- the first latent axis's technical correlations for the primary and both arms, from `axis-diagnostics`;
- each `latent_adjusted` arm's `feature_r2` summary;
- provenance: the MOTCO commit, this repository's revision, seed, permutations, `n_jobs` and, for the astrocyte report, the pin-bump check.

Every statistic, p-value and rank in the report SHALL be checked against the records before committing.

#### Scenario: Report checked
- **WHEN** the report is committed
- **THEN** a check finds every statistic, p-value and rank of the two reference and four arm records in the report text

#### Scenario: Microglia report
- **WHEN** the arms have run on microglia in both cohorts
- **THEN** `reports/motco-microglia-latent-adjustment-<date>.md` is committed with the two comparison tables, both mediator checks, the first-axis correlations of the primary and both arms, the `feature_r2` summaries and the provenance, every value checked against the six microglia records

### Requirement: Report keeps the primary analysis primary
The report SHALL state that:
- the primary analysis stays primary, and that MOTCO's power comparison (`candidate_fails`) excludes latent adjustment as its replacement;
- a rejection status that differs between an arm and the primary is covariate-set dependence (primary vs `extended_covariates`) or basis dependence (`extended_covariates` vs `latent_adjusted`), not a result of the arm;
- a CV-selected rank that differs is reported beside `angle`, not as better or worse;
- `angle` and `shape` non-rejections stay uninformative and a `delta` rejection is a difference, not specifically of size, in every arm.

A technical correlation of an arm's first axis SHALL be described as a measurement-space finding only, never as a reinterpretation of any arm's statistics, and never as grounds to adopt an arm.

#### Scenario: Arm rejects where the primary does not
- **WHEN** `latent_adjusted` rejects `angle` in `cohort_30` and the primary does not
- **THEN** the report calls this basis dependence, reports the arm's selected rank beside `angle`, and keeps the primary's non-rejection as the result

### Requirement: Report carries the limitations
The report SHALL carry:
- the mediator limitation, with the mediator check's flags: a flagged fraction may carry stage signal that the `latent_adjusted` basis removed, and no simulation can tell a confounder from a mediator;
- MOTCO's latent-adjustment operating characteristics, measured at the Phase 6 design only: Type I `adopt_stage`; power `candidate_fails`, with `shape` power lost under a nuisance axis and small-magnitude power lost without one;
- the facts that no operating characteristic was measured on the SEA-AD blocks or at this eight-covariate set, that the two arms carry no Type I or power claim of their own, and that the nuclei counts and the fractions both proxy library quality, which the report states.

#### Scenario: Mediator flag reported
- **WHEN** a fraction is flagged in `mediator_check.csv`
- **THEN** the report names the fraction, the cohort and ρ, and states that the `latent_adjusted` arm may have removed stage signal carried by it

### Requirement: Arms run per registered cell type
The step SHALL take a registered cell type (default `astrocyte`) and SHALL run the registered arms, with the covariate set, conditioning, mediator check, seed, permutation count and CV rank selection unchanged, on that cell type's primary export against that cell type's primary analysis. The registry of arms and the covariate set SHALL NOT vary with the cell type. The default SHALL reproduce the astrocyte outputs byte for byte.

#### Scenario: Microglia arms complete
- **WHEN** the step runs with the cell type `microglia` for `cohort_30`
- **THEN** two arm records exist under the microglia results root, `extended_covariates` and `latent_adjusted`, both listing the eight covariates in the registered order, with the microglia primary analysis as the reference row, the primary's seed and 999 permutations

#### Scenario: Default unchanged
- **WHEN** the step runs without a cell type on the astrocyte inputs
- **THEN** every output equals the output of the step before this change, apart from the recorded cell type in its manifests

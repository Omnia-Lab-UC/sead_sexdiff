# motco-rank-sensitivity Specification

## Purpose

Measures how the covariate-adjusted sex × ADNC-stage trajectory comparison of the astrocyte cohorts depends on the latent rank, by rerunning MOTCO's analysis entry point at a predeclared ladder of forced ranks beside the CV-selected analysis, and records the sweep so that it is read as sensitivity to the primary analysis, never as its replacement.

## Requirements

### Requirement: Inputs are the export and the CV-selected analysis
For each requested cohort the sweep SHALL read the export (`rna.csv`, `atac.csv`, `metadata.csv`, `covariates.csv`, `manifest.json` under `results/motco_input/<cohort>/`) and the CV-selected analysis (`result.json`, `analysis_manifest.json` under `results/motco_analysis/<cohort>/`), and nothing else. It SHALL stop with an error naming the directory and the command that produces it when either is missing.

#### Scenario: Analysis not yet run
- **WHEN** `results/motco_analysis/cohort_30/` does not exist
- **THEN** the sweep stops with an error naming the directory and the `analyze` command, and writes nothing

#### Scenario: Export not yet run
- **WHEN** `results/motco_input/cohort_50/` does not exist
- **THEN** the sweep stops with an error naming the directory and the `export` command, and writes nothing

### Requirement: Reference analysis must match the sweep's inputs and settings
Before any rank runs, the sweep SHALL check that the reference analysis manifest records the same sha256 for each of the five export files as the files now on disk, and the same MOTCO version, seed, permutation count and latent-rank floor as the sweep uses. On any difference it SHALL stop with an error naming the file or setting.

#### Scenario: Export regenerated after the analysis
- **WHEN** `rna.csv` was rewritten after `analyze` ran and its sha256 differs from the one in `analysis_manifest.json`
- **THEN** the sweep stops with an error naming `rna.csv` and the `analyze` command, and writes nothing

#### Scenario: Matching reference
- **WHEN** every hash and setting matches
- **THEN** the reference row of the sweep carries the analysis record's statistics, p-values, null summaries and CV-selected rank unchanged

#### Scenario: Unfloored reference
- **WHEN** the reference analysis was made with the rank floor 1
- **THEN** the sweep stops with an error naming the reference `result.json`, both floors and the `analyze` command, and writes nothing

### Requirement: Predeclared rank ladder
The sweep SHALL run the ranks 2, 3, 4, 6, 9 and 12, each once, through MOTCO's trajectory-analysis entry point with the rank forced. The ladder SHALL be a fixed constant, not a runtime option, and the sweep SHALL NOT add, drop or choose ranks from the results. Before any run it SHALL stop with an error naming the rank when a rank exceeds what the cohort can support.

#### Scenario: Ladder completes
- **WHEN** the sweep finishes for `cohort_30`
- **THEN** six per-rank records exist, each with `component_selection = forced` and a selected rank equal to its forced rank, and no other rank was run

#### Scenario: Infeasible rank
- **WHEN** a cohort has fewer donors than the largest rank of the ladder
- **THEN** the sweep stops with an error naming that rank before writing anything

### Requirement: Every setting other than the rank equals the analysis
Each forced-rank run SHALL use the same blocks and units, stage order, six covariates in the same order, missing-covariate policy, orientation estimator, permutation count and seed as the CV-selected analysis, so that every rank draws the same permutation indices and differences between ranks are measurement differences. Cross-validation SHALL NOT run.

#### Scenario: Record carries the settings
- **WHEN** the rank-4 run of `cohort_30` finishes
- **THEN** its `result.json` records `forced_components = 4`, `component_selection = forced`, `permutations = 999`, the analysis seed, `orientation_estimator = net_displacement`, the six covariate names in order and no cross-validation AUROC

#### Scenario: CV-selected rank inside the ladder
- **WHEN** the reference analysis's CV-selected rank is one of the ladder's ranks
- **THEN** the forced run at that rank gives the same statistics and p-values as the reference analysis, and any difference is recorded as a finding

### Requirement: Per-rank outputs
Each forced-rank run SHALL be written to `results/motco_rank_sensitivity/<cohort>/rank_<kk>/` (`kk` zero-padded to two digits) holding the same files as a CV-selected analysis directory: MOTCO's artifact (`result.json`, `latent_scores.csv`, `pls_loadings.csv`, `stage_means.csv`), `analysis_metadata.csv`, `analysis_manifest.json` and `summary.md`. The sweep SHALL refuse to overwrite an existing cohort directory unless forced.

#### Scenario: Rerun without force
- **WHEN** `results/motco_rank_sensitivity/cohort_30/` exists and the sweep runs again without `--force`
- **THEN** it stops with an error naming the directory and writes nothing

#### Scenario: Forced rank beside angle
- **WHEN** `rank_06/summary.md` is written
- **THEN** its `angle` row names forced rank 6 and says the rank was forced, and its latent block says cross-validation was skipped

### Requirement: Cross-rank table generated from the records
The sweep SHALL write `results/motco_rank_sensitivity/<cohort>/rank_sensitivity.csv` and its Markdown rendering `rank_sensitivity.md`, generated from the records and never typed by hand, with the reference row first (`cv`, the CV-selected rank) followed by one row per forced rank in ascending order.

#### Scenario: Table rows
- **WHEN** the sweep finishes for `cohort_30`
- **THEN** `rank_sensitivity.csv` has seven rows, the first with `component_selection = cv` and the reference rank, the rest with ranks 2, 3, 4, 6, 9, 12 in that order

### Requirement: Cross-rank table columns
Each row SHALL give, per statistic (`delta`, `angle`, `shape`), the observed value, the RRPP p-value, whether p ≤ α and the null mean, sd and q95, plus the reference row's CV AUROC and the row's directory. The CSV SHALL hold only values that reproduce across reruns; each run's runtime goes to the Markdown and the manifest.

#### Scenario: Table values equal the records
- **WHEN** a cell of `rank_sensitivity.csv` is compared with the `result.json` of its row's directory
- **THEN** they are equal

#### Scenario: Floating-point cells round-trip
- **WHEN** a `delta` null sd whose shortest representation needs 17 significant digits is written to `rank_sensitivity.csv`
- **THEN** the cell is written in its shortest round-trip representation and parses back to the record's value exactly

#### Scenario: Runtime and directories
- **WHEN** the table of `cohort_30` is written
- **THEN** it has no runtime column, its forced rows name their directory as `rank_<kk>` relative to the cohort directory and its reference row names the reference analysis directory, while `rank_sensitivity.md` and the manifest carry each run's runtime

### Requirement: Sweep manifest records provenance
`rank_sensitivity_manifest.json` SHALL record the cohort, the ladder, the seed, the permutation count, the rank floor, α, `n_jobs`, the reference analysis directory with the sha256 of its `result.json` and `analysis_manifest.json`, the sha256 of each export input, each per-rank directory, the MOTCO, Python and library versions, the git revision and findings, including that the CV-selected rank is not in the ladder when it is not.

#### Scenario: Manifest after a cohort sweep
- **WHEN** `cohort_30` finishes with a CV-selected rank of 2
- **THEN** its manifest is valid JSON listing the six ranks, `permutations = 999`, the seed, the rank floor 2, the MOTCO version, hashes for the five export files and the two reference files, and the comparison of the forced rank-2 run with the reference analysis

#### Scenario: CV-selected rank outside the ladder
- **WHEN** a cohort's CV-selected rank is 5
- **THEN** the manifest carries a finding that rank 5 is not in the ladder and the reference row stands in for it

### Requirement: The CV-selected analysis is never modified
The sweep SHALL NOT write into `results/motco_analysis/`.

#### Scenario: Analysis directory after a sweep
- **WHEN** the sweep finishes for both cohorts
- **THEN** every file under `results/motco_analysis/` is byte-identical to before the sweep

### Requirement: Reproducible sweep
Two sweeps on the same inputs, reference, MOTCO release, seed and `n_jobs` SHALL give identical statistics, p-values and null summaries at every rank, and byte-identical `latent_scores.csv`, `pls_loadings.csv`, `stage_means.csv` and `rank_sensitivity.csv`. A run with `n_jobs` other than 1 SHALL carry a finding that its permutation draws differ from an `n_jobs = 1` run, and the report SHALL use `n_jobs = 1`.

#### Scenario: Rerun with force
- **WHEN** `cohort_50` is swept twice with `--force` and the same `n_jobs`
- **THEN** every per-rank `result.json` pair has equal statistics, p-values and null summaries, and the CSVs are byte-identical

### Requirement: Committed rank-sensitivity report
A Markdown report under `reports/` SHALL be committed, giving for each cohort the cross-rank table (reference row and every forced rank) pasted from the generated files, the MOTCO version, seed, permutation count and `n_jobs`, and a pointer to the primary report for the covariate block and the feature filtering. Every number SHALL be checked against the matching `result.json` before commit. `cohort_30` SHALL be presented as primary and `cohort_50` as sensitivity.

#### Scenario: Report matches the records
- **WHEN** the report is compared with the per-rank and reference `result.json` files
- **THEN** every statistic, p-value, null summary and rank in the report equals the record's value

### Requirement: Report keeps the CV-selected analysis primary
The report SHALL state that the CV-selected analysis remains the primary result and that no rank is chosen after seeing the results. A statistic whose rejection status differs between ranks SHALL be reported as rank dependence of the sensitivity analysis, not as a result at the rank where it is rejected.

#### Scenario: Rejection at a forced rank only
- **WHEN** `delta` is rejected at α = 0.05 at rank 4 but not in the CV-selected analysis
- **THEN** the report records the rejection as rank dependence of the sensitivity analysis, keeps the CV-selected non-rejection as the primary result, and does not present rank 4 as the analysis

#### Scenario: No rejection at any rank
- **WHEN** no statistic is rejected at any rank in a cohort
- **THEN** the report states that the sweep gives no evidence of a difference at any rank and that this is not evidence of agreement in size, direction or shape

### Requirement: Report carries the forced-rank reading rules
The report SHALL state that `angle` and `shape` non-rejections stay uninformative at every rank, that a `delta` rejection is a difference and not specifically of size, that MOTCO's forced-rank operating characteristics are carried as limitations, and that the covariate-adjusted analysis has no Type I check at any forced rank.

#### Scenario: Forced-rank limitations
- **WHEN** the report discusses a forced rank
- **THEN** it states that MOTCO's Type I evidence at a forced rank covers rank 4 at this design without covariates, that forced 2(K − 1) lost `shape` power and leaked `angle` into the magnitude cell in MOTCO's probe, and that ranks above `n_stages − 1` widened the `angle` null in MOTCO's rank ladder

#### Scenario: Non-rejection wording at a forced rank
- **WHEN** `angle` is not rejected at rank 6 in a cohort
- **THEN** the report states that the result is uninformative about orientation at that rank, not that the directions agree

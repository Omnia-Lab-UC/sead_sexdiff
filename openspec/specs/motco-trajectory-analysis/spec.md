# motco-trajectory-analysis Specification

## Purpose

Runs the covariate-adjusted comparison of female and male astrocyte trajectories across ADNC stages on the exported cohorts through MOTCO's validated analysis entry point, and records the result so that every number in the report can be traced to a pinned MOTCO release, fixed inputs and a fixed seed.

## Requirements

### Requirement: Inputs are the exported MOTCO blocks only
For each requested cohort the analysis SHALL read `rna.csv`, `atac.csv`, `metadata.csv`, `covariates.csv` and `manifest.json` from `results/motco_input/<cohort>/` and nothing else. It SHALL stop with an error naming the file when one is missing, when the donor columns of the four tables differ, or when the export manifest's donor count or block feature counts disagree with the tables.

#### Scenario: Export not yet run
- **WHEN** `results/motco_input/cohort_50/` does not exist
- **THEN** the analysis stops with an error naming the directory and the `export` command that produces it

#### Scenario: Tables disagree
- **WHEN** `covariates.csv` lists a donor that `metadata.csv` does not
- **THEN** the analysis stops with an error naming both files

### Requirement: Fixed, recorded covariate set
The analysis SHALL enter exactly these covariates, additively and in this order: `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`, where `log10_n_rna` and `log10_n_atac` are the base-10 logarithms of `n_rna` and `n_atac`. `education_level` and `age_capped` SHALL NOT enter the model. The missing-covariate policy SHALL be `error`, and no value SHALL be imputed.

#### Scenario: Derived nuclei covariates
- **WHEN** a donor has `n_rna = 1000` and `n_atac = 100`
- **THEN** the model sees `log10_n_rna = 3` and `log10_n_atac = 2`, and the record's covariate block lists the six covariate names in the order above

#### Scenario: Missing covariate
- **WHEN** a cohort donor has an empty `rin` cell
- **THEN** the analysis stops with an error naming `rin` and the donor, and no result is written

### Requirement: Analysis through MOTCO's validated entry point
The analysis SHALL run through MOTCO's trajectory-analysis entry point (`motco analyze` / `analyze_trajectories`) of a pinned MOTCO release ≥ 0.8.0. It SHALL NOT assemble the result from `motco plsr`, `motco de` or the statistics primitives.

#### Scenario: Older MOTCO installed
- **WHEN** the installed MOTCO is older than 0.8.0
- **THEN** the analysis stops with an error naming the installed and required versions before reading any input

#### Scenario: Record names the entry point's release
- **WHEN** a cohort analysis finishes
- **THEN** `result.json` records the MOTCO version under `software`

### Requirement: Fixed measurement settings
The analysis SHALL declare blocks `rna` and `atac` in identity units, group column `sex`, stage column `stage` with the explicit order `0`, `1`, `2`, the net-displacement orientation estimator, cross-validated rank selection at the MOTCO defaults, 999 permutations and a fixed recorded seed.

#### Scenario: Record carries the measurement settings
- **WHEN** a cohort analysis finishes
- **THEN** `result.json` records `orientation_estimator = net_displacement`, `component_selection = cv`, `permutations = 999`, the seed, the ordered stage levels `0, 1, 2` and the group levels

### Requirement: Analysis outputs per cohort
The analysis SHALL write `results/motco_analysis/<cohort>/` containing MOTCO's artifact (`result.json`, `latent_scores.csv`, `pls_loadings.csv`, `stage_means.csv`), `analysis_metadata.csv` (the merged metadata and covariate table the model saw, one row per donor sorted by `donor_id`), `analysis_manifest.json` and `summary.md`. It SHALL refuse to overwrite an existing cohort directory unless forced.

#### Scenario: Rerun without force
- **WHEN** `results/motco_analysis/cohort_30/` exists and the analysis runs again without `--force`
- **THEN** it stops with an error naming the directory and writes nothing

#### Scenario: Metadata table is self-sufficient
- **WHEN** `analysis_metadata.csv` is passed to `motco analyze` with the covariate flags listed in the manifest
- **THEN** MOTCO accepts it without further reshaping

### Requirement: Analysis manifest records provenance
`analysis_manifest.json` SHALL record the cohort, the sha256 of every input file and of the export manifest, the measurement settings (estimator, rank selection, permutations, seed, α), the MOTCO, Python and library versions, the git revision and findings.

#### Scenario: Manifest after a cohort run
- **WHEN** `cohort_30` finishes
- **THEN** its manifest is valid JSON holding a sha256 for each of the five input files, `permutations = 999`, the seed and the MOTCO version

### Requirement: Analysis manifest records the model inputs
`analysis_manifest.json` SHALL record the covariate names in model order, the derivation of the log10 columns, the exported columns that were not entered with the reason, and the preprocessing summary copied from the export manifest (excluded chromosomes, `min_cpm`, `k`, `prior_count`, `N`, feature counts per step).

#### Scenario: Covariates and filtering are traceable
- **WHEN** `cohort_30` finishes
- **THEN** its manifest lists the six covariate names in order, names `education_level` and `age_capped` as not entered, and holds `k = 9` and `N = 5000`

### Requirement: Summary table generated from the record
`summary.md` SHALL be generated from `result.json`, not typed by hand, and SHALL show per statistic (`delta`, `angle`, `shape`) the observed value, the RRPP p-value and the null mean, sd and 95 % quantile; the CV-selected rank and `cv_mean_auroc`; the sex × stage cell sizes; and the covariate block (design columns, reference row, missing policy, dropped donors).

#### Scenario: Rank beside angle
- **WHEN** `summary.md` is written
- **THEN** the `angle` row names the CV-selected rank from `result.json`

### Requirement: Reproducible results
Two runs on the same inputs, MOTCO release, seed and `n_jobs` SHALL give identical statistics, p-values, selected rank and null summaries, and byte-identical `latent_scores.csv`, `pls_loadings.csv`, `stage_means.csv` and `analysis_metadata.csv`. Because MOTCO's RRPP draws depend on `n_jobs`, the manifest SHALL record `n_jobs`, a run with `n_jobs` other than 1 SHALL carry a finding saying so, and the report SHALL use `n_jobs = 1`.

#### Scenario: Rerun with force
- **WHEN** a cohort is analysed twice with `--force` and the same `n_jobs`
- **THEN** the statistics, p-values and selected rank in the two `result.json` files are equal and the CSVs are byte-identical

### Requirement: Committed case-study report
A Markdown report under `reports/` SHALL be committed, giving for each cohort the statistics, p-values and null summaries, the CV-selected rank beside `angle`, the covariate block, the cell sizes, the feature filtering beside the results, and the MOTCO version and seed. `cohort_30` SHALL be presented as primary and `cohort_50` as sensitivity.

#### Scenario: Report matches the records
- **WHEN** the report is compared with the two `result.json` files
- **THEN** every statistic, p-value and selected rank in the report equals the record's value

#### Scenario: Cohort roles
- **WHEN** the report presents both cohorts
- **THEN** `cohort_30` is labelled primary and `cohort_50` sensitivity

### Requirement: Report applies the binding interpretation
The report SHALL state that a `delta` rejection is evidence of a difference and not specifically of size; that `angle` and `shape` non-rejections are uninformative at this n; that `angle` reads only the endpoint stages; that the adjustment is additive-linear on a covariate-blind basis; and that the operating characteristics are a structural analogue (InterSIM blocks, 498 features).

#### Scenario: Non-rejection wording
- **WHEN** `angle` is not rejected at α = 0.05 in a cohort
- **THEN** the report states that the result is uninformative about orientation, not that the directions agree

#### Scenario: Delta rejection wording
- **WHEN** `delta` is rejected at α = 0.05 in a cohort
- **THEN** the report states that the sexes' trajectories differ, without claiming the difference is in size

# Spec Delta

## Purpose

Measures how the covariate-adjusted sex × ADNC-stage trajectory comparison and the pseudobulk stage DE depend on the donors SEA-AD flagged as severely affected, by rerunning the whole export and the primary analysis on the cohort donors without them, per registered cell type, with a predeclared check of whether the quality gradient of the primary analyses sits in those donors. The variant is recorded so that it is read as sensitivity to the primary analysis, never as its replacement.

## ADDED Requirements

### Requirement: Inputs are the pseudobulk outputs, the flag table, the primary export and analysis
For each requested cohort the step SHALL read only, under the cell type's results root: the pseudobulk outputs and the covariate table, the gene annotation, the flag table with its manifest, the technical table with its manifest, the primary export and the primary analysis. Every manifest read SHALL record the cell type; a mismatch SHALL stop the step naming the file and both cell types. A missing file SHALL stop the step naming it and the command that produces it. Nothing SHALL be written before every check passes.

#### Scenario: Flag table absent
- **WHEN** `results/microglia/pseudobulk/donor_flags.csv` is missing
- **THEN** the step stops with an error naming the file and the `donor-flags --cell-type microglia` command, and writes nothing

#### Scenario: Cell-type mismatch
- **WHEN** the step runs for `microglia` and the flag table's manifest records `astrocyte`
- **THEN** the step stops with an error naming the manifest and both cell types, and writes nothing

### Requirement: Primary export and analysis must match
Before any run the step SHALL check that the primary export records the pseudobulk, covariate and annotation hashes now on disk, the thresholds, `k` and the cells of the cohort, and that the primary analysis records the export files' sha256, the installed MOTCO version, the seed, the permutation count and CV rank selection. The flag table SHALL hold a value of the flag for every cohort donor. On any difference the step SHALL stop naming the file, setting or donor, and write nothing.

#### Scenario: Flag missing for a cohort donor
- **WHEN** a cohort donor has no row in the flag table
- **THEN** the step stops with an error naming the donor, and writes nothing

### Requirement: Predeclared exclusion registry
The step SHALL run only variants of a fixed registry. The registry SHALL hold `excluding_severely_affected`: exclude the cohort donors whose `severely_affected` is true. A variant SHALL NOT be a runtime option beyond its name, SHALL NOT be derived from any result, and entries SHALL be added by a change proposal only.

#### Scenario: Unknown variant
- **WHEN** the step is asked for a variant not in the registry
- **THEN** it stops with an error listing the registered variants, and writes nothing

### Requirement: Flag check written before any run
Per cohort, before the variant export runs, the step SHALL write `flag_check.csv`: the flagged donors per sex × stage cell; for each numeric column of the technical table, each TMM factor and log10 library size of the primary export and the primary analysis's first latent component, the rank-biserial AUC of flagged against unflagged donors, overall and within stage 2; and the first component's Spearman correlation with stage and with each of those variables among the unflagged donors. It SHALL compute no p-value, and its values SHALL change no run.

#### Scenario: Flagged donors carry the tail
- **WHEN** the flagged donors of microglia `cohort_30` have AUC 0.95 for `mitochondrial_pseudobulk` and 0.9 for the primary first component
- **THEN** `flag_check.csv` records both values and the counts per cell, and the variant runs unchanged

### Requirement: Variant export is the primary pipeline on the unflagged cohort donors
The variant export SHALL run the export's pipeline unchanged (the primary's excluded chromosomes, abundance threshold, TMM, prior count, selection and normalization) on the cohort donors minus the excluded ones, with `k` recomputed from the remaining sex × stage cells by the export's rule and the cohort thresholds unchanged. Its manifest SHALL record the exclusion: the variant, the flag, the flag table's sha256, the excluded donors, and the cells and `k` before and after.

#### Scenario: Smallest cell recomputed
- **WHEN** excluding the flagged donors reduces the Male stage-2 cell of microglia `cohort_30` from 11 to 8
- **THEN** the variant export uses `k = 8`, and its manifest records the cells and `k` before and after and the excluded donors

#### Scenario: No donor excluded
- **WHEN** no cohort donor is flagged
- **THEN** the variant export equals the primary export's values, the manifest records an empty exclusion, and the comparison carries a finding

### Requirement: Variant analysis uses the primary analysis settings
The variant analysis SHALL run on the variant export with the primary analysis's covariates, orientation, CV rank selection, seed and permutation count, so the variant and the primary differ in the donors only.

#### Scenario: Settings recorded
- **WHEN** the variant analysis of `cohort_30` finishes
- **THEN** its record lists the six analysis covariates, `component_selection = cv`, the primary's seed and 999 permutations

### Requirement: Predeclared axis check
After the variant analysis, the step SHALL write `axis_check.csv`: along the variant analysis's first latent component, the |Spearman| with stage and the largest |Spearman| with a technical variable (the technical table's numeric columns for the variant donors and the variant export's TMM factors and log10 library sizes), and whether the first exceeds the second. The rule SHALL be a constant written before any run. For the control cell type the outcome is the change's expectation; for the reference cell type it is descriptive. No p-value SHALL be computed.

#### Scenario: Axis check holds
- **WHEN** the microglia variant first component has |ρ| 0.71 with stage and a largest technical |ρ| of 0.52
- **THEN** `axis_check.csv` records both values, the variable, and `holds = true`

#### Scenario: Axis check fails
- **WHEN** the largest technical |ρ| exceeds the stage |ρ|
- **THEN** `axis_check.csv` records `holds = false` naming the variable, and nothing else changes

### Requirement: Output layout mirrors the primary pipeline
The step SHALL write, under `<root>/motco_exclusion_sensitivity/<variant>/`: `flag_check/<cohort>/`, `motco_input/<cohort>/` (the variant export, the files `export` writes), `motco_analysis/<cohort>/` (the files `analyze` writes) and `comparison/<cohort>/`. Pointing `axis-diagnostics`, `stage-de` or `stage-de-limma` at the variant's `motco_input/` and `motco_analysis/` SHALL pass their input checks. A cohort SHALL be written under temporary names and renamed when complete; an existing cohort output SHALL be an error without `--force`.

#### Scenario: Diagnostics on the variant
- **WHEN** `axis-diagnostics --cell-type microglia` runs with the variant's `motco_input/` and `motco_analysis/`
- **THEN** its input-hash and cell-type checks pass and it reports the variant's components against the predeclared variables

### Requirement: Comparison table generated from the records
Per cohort the step SHALL write `exclusion_sensitivity.csv` with the primary row and the variant row (donors, cells, `k`, selected rank, CV mean AUROC, and for `delta`, `angle` and `shape` the observed value, p-value, rejection and null mean, SD and q95), its Markdown rendering and a manifest recording the variant, the excluded donors, every input path and sha256, seed, permutations, `n_jobs`, software, git revision and findings. Every value SHALL be read from the records, not computed again.

#### Scenario: Values traceable
- **WHEN** `exclusion_sensitivity.csv` is written
- **THEN** every statistic, p-value and rank equals the corresponding value in the primary or variant `result.json`

### Requirement: The primary export and analysis are never modified
The step SHALL NOT write under the cell type's `pseudobulk/`, `motco_input/` or `motco_analysis/`, and the primary row SHALL be taken from the primary record unchanged.

#### Scenario: Primary files unchanged
- **WHEN** the step finishes for both cohorts
- **THEN** the sha256 of every file under the primary directories equals its value before

### Requirement: Reproducible run
Two runs on the same inputs, MOTCO release, seed and `n_jobs` SHALL give byte-identical variant export CSVs, byte-identical `flag_check.csv`, `axis_check.csv` and `exclusion_sensitivity.csv`, and identical variant statistics, p-values and ranks. A run with `n_jobs` other than 1 SHALL carry a finding, and the report SHALL use `n_jobs = 1`.

#### Scenario: Rerun with force
- **WHEN** a cohort is run twice with `--force` and `n_jobs = 1`
- **THEN** both `exclusion_sensitivity.csv` files are byte-identical

### Requirement: Stage DE on the exclusion cohort
The change SHALL run the DESeq2 and the limma-voom stage-DE steps, RNA only, with their registered models and settings unchanged, on the variant export and analysis of both cell types and both cohorts, writing under the variant root. The report SHALL list their exported RNA stage counts beside the primary's and both concordance cells against the variant analysis.

#### Scenario: Counts beside the primary's
- **WHEN** both stage-DE steps have run on the microglia exclusion cohort of `cohort_30`
- **THEN** the report lists the DESeq2 `factor_additive` union and the limma `stage_joint` exported counts of the variant beside the primary's 61 and 2, with the two concordance cells

### Requirement: Committed report with the reading rules
The change SHALL commit `reports/motco-donor-exclusion-sensitivity-<date>.md`, both cell types, with numbers pasted from the generated files and checked against them: provenance, the flag table summary, the flag checks, the comparison tables, the axis checks, the stage-DE counts and concordance cells, the variant's first-axis correlations, and the predeclared reading: the primary analyses and the primary-basis positive control stay as they are; a rejection status that differs is cohort dependence, not a result of the variant; a CV-selected rank that differs is reported beside `angle`; the non-rejection rules hold in the variant; the cross-cell-type comparisons are compared, not tested.

#### Scenario: Report numbers checked
- **WHEN** the report is committed
- **THEN** a check finds every statistic, p-value, rank and AUROC of the four primary and four variant records, every flag-check and axis-check value, and every stage-DE count in the text

### Requirement: Report carries the limitations
The report SHALL state that the excluded donors carry the most severe pathology, so the exclusion truncates the late epoch along with the quality tail and a stage effect concentrated in them is removed with them; that the flag was made by SEA-AD from QC metrics, so the variant conditions on quality by donor selection; that the stage-2 cells are smaller, so the variant has less power than the primary; that MOTCO's operating characteristics were measured on neither cohort; and that the flag's sex split, reported, bears on the sex × stage comparison.

#### Scenario: Limitations present
- **WHEN** the report is committed
- **THEN** it names the late-epoch truncation, the QC origin of the flag, the smaller cells, the absent operating characteristics and the flag's sex split

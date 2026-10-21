# Spec Delta

## Purpose

Fits covariate-adjusted per-feature DESeq2 models of ADNC stage and of the sex × stage interaction on the cohort donors' raw pseudobulk counts, on the feature set of the primary MOTCO export, and evaluates a predeclared concordance rule between the interaction count and the primary analysis's `delta` rejection. It describes the stage signal the latent analysis works on and cross-checks that analysis; it never replaces it.

## ADDED Requirements

### Requirement: Inputs are the pseudobulk, the primary export and the primary analysis
For each requested cohort the step SHALL read only:
- the pseudobulk outputs (`rna_counts.h5ad`, `atac_counts.h5ad`, `donor_table.csv`, `manifest.json` under the pseudobulk directory);
- the primary export (`metadata.csv`, `covariates.csv`, `features_rna.csv`, `features_atac.csv`, `manifest.json` under `<export root>/<cohort>/`, and `rna.csv`, `atac.csv` for the analysis hash check and the reproduction check only);
- the primary analysis (`result.json`, `analysis_manifest.json` under `<analysis root>/<cohort>/`).

It SHALL NOT read donor objects, and it SHALL stop with an error naming the missing file and the command that produces it, writing nothing.

#### Scenario: Export absent
- **WHEN** `results/motco_input/cohort_50/features_atac.csv` is missing
- **THEN** the step stops with an error naming the file and the `export` command, and writes nothing

#### Scenario: Analysis absent
- **WHEN** `results/motco_analysis/cohort_30/result.json` is missing
- **THEN** the step stops with an error naming the file and the `analyze` command, and writes nothing

### Requirement: References must match
Before any fit the step SHALL check that the export manifest records the sha256 of the pseudobulk manifest now on disk, that donors, sex and stage agree between the export, the donor table and the count files, and that the analysis manifest records the export files' sha256 now on disk, the installed MOTCO version, the seed and the permutation count. On any difference it SHALL stop naming the file or setting and write nothing. The manifest SHALL record the sha256 of every input file read.

#### Scenario: Pseudobulk manifest changed since the export
- **WHEN** the pseudobulk `manifest.json` has a sha256 different from the one the export manifest records
- **THEN** the step stops with an error naming the file and the `export` command, and writes nothing

### Requirement: Counts reproduce the export
Because the export records no count-file hash, the step SHALL rerun the export's preprocessing on the cohort donors' counts (chromosomes from `features_<omic>.csv`, recorded thresholds) before any fit, and SHALL require the recorded statuses and per-donor library sizes exactly and the TMM factors, tested-feature log2-CPM variances and exported values to rounding; otherwise it SHALL stop naming the count file and the `export` command, writing nothing.

#### Scenario: Pseudobulk changed since the export
- **WHEN** `rna_counts.h5ad` no longer reproduces the feature statuses, library sizes, TMM factors, tested-feature variances or exported values the export recorded for the cohort donors
- **THEN** the step stops with an error naming the file and the `export` command, and writes nothing

### Requirement: Feature set bound to the primary export
For each omic the step SHALL test exactly the features that passed the export's chromosome and abundance filters, as recorded in `features_<omic>.csv` (status `variance` or `kept`), in the export's feature order. Size factors and dispersions SHALL be fitted on that set. Every result row SHALL carry `exported`, true for the features of MOTCO's block (status `kept`). The step SHALL NOT select features from sex, stage, covariates or its own results.

#### Scenario: Tested set equals the abundance-passing set
- **WHEN** the RNA model of `cohort_30` is fitted
- **THEN** the number of features in every result table equals the export manifest's post-abundance RNA count (20,387), of which 5,000 are marked `exported`

#### Scenario: Sex-chromosome features are not tested
- **WHEN** a gene has status `sex_chromosome` in `features_rna.csv`
- **THEN** it appears in no result table

### Requirement: Predeclared models and tested coefficients
The step SHALL fit, per cohort and omic, exactly the models of a fixed registry, each with sex (reference `Female`), stage and the omic's covariates:
- `factor_additive`: sex + stage (three-level factor, reference `0`) + covariates; tested coefficients `stage 1 vs 0`, `stage 2 vs 0`.
- `factor_interaction`: sex + stage + sex × stage + covariates; tested coefficients `Male × stage 1`, `Male × stage 2`.
- `trend_additive`: sex + stage as 0/1/2 + covariates; tested coefficient `stage_num`.
- `trend_interaction`: sex + stage_num + sex × stage_num + covariates; tested coefficient `Male × stage_num`.

Each tested coefficient SHALL get one Wald test with Benjamini–Hochberg adjustment at α = 0.05 and the DESeq2 defaults (independent filtering, Cook's outlier refit). The registry, the reference levels and α SHALL be constants, not runtime options, and SHALL NOT be changed from results.

#### Scenario: Result tables per coefficient
- **WHEN** the step finishes `cohort_30`, omic `rna`
- **THEN** `results/stage_de/cohort_30/rna/` holds six tables, one per tested coefficient, each with `feature`, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`, `exported`, and the manifest records the four formulas, the tested coefficients, the reference levels and α

#### Scenario: Union under a factor model
- **WHEN** the summary counts stage-significant features under `factor_additive`
- **THEN** it reports the count per coefficient and the count of features significant on either coefficient, labelled as a union without a joint error-rate claim

### Requirement: Predeclared covariates, never imputed
Each omic's model SHALL enter, as continuous covariates, `age_at_death`, `pmi`, `rin`, `years_education` and the omic's own log10 nuclei count (`log10_n_rna` for RNA, `log10_n_atac` for ATAC), read from the primary export's `covariates.csv` and the metadata. A missing value SHALL be an error naming the donor and the column; nothing SHALL be imputed.

#### Scenario: Covariate columns recorded
- **WHEN** the ATAC models of a cohort are fitted
- **THEN** the manifest lists `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_atac` as that omic's covariates, and `log10_n_rna` is not in the ATAC design

#### Scenario: Missing covariate
- **WHEN** a cohort donor has an empty `rin` cell
- **THEN** the step stops with an error naming the donor and the column, and writes nothing for that cohort

### Requirement: Summary counts on two feature sets
The step SHALL write `stage_de_summary.csv` with, per cohort × omic × model × tested coefficient (and the union per factor model), the number of tested features and the number significant at α, each on the exported set and on the abundance-passing set, plus the smallest `padj`. The exported-set counts SHALL be the primary reading.

#### Scenario: Summary rows
- **WHEN** both omics of a cohort have run
- **THEN** that cohort's `stage_de_summary.csv` has a row for every omic × coefficient and union × `feature_set` (`exported`, `abundance`), with `n_tested` and `n_significant`, and the exported `n_tested` is 5,000 in every row

### Requirement: Predeclared concordance rule against the primary analysis
For each cohort and omic the step SHALL evaluate, from the generated files: `interaction_detected` (at least one exported feature with `padj ≤ α` on either interaction coefficient of `factor_interaction`) against `delta_rejected` read from the primary analysis's `result.json` under the analysis rejection rule. The cell SHALL be `concordant` when both are true or both are false, `discordant` otherwise, and a discordance SHALL be written as a finding. The rule SHALL read only the binary `delta` rejection status (derived from `result.json` under the analysis rejection rule, as the comparison tables derive it) and no p-value magnitude, SHALL change no MOTCO output, and SHALL be a constant.

#### Scenario: Both negative
- **WHEN** no exported RNA feature has an interaction `padj ≤ 0.05` in `cohort_30` and the primary analysis did not reject `delta`
- **THEN** the summary's concordance cell for `cohort_30`, `rna` is `concordant` with `interaction_detected = false`, `delta_rejected = false`

#### Scenario: Interaction without a delta rejection
- **WHEN** 12 exported ATAC features have an interaction `padj ≤ 0.05` in `cohort_50` and the primary analysis did not reject `delta`
- **THEN** the cell is `discordant`, the manifest records a finding naming the cohort, omic and count, and the primary analysis's files are unchanged

### Requirement: Outputs, manifest and reproducibility
The step SHALL write, per cohort, the result tables, a fit record per omic × model (`<omic>/<model>__fit.json`: the fingerprint of the inputs and settings, the fit, the reproduction of the export, each table's sha256), `stage_de_summary.csv`, `stage_de.md` (rendered from the summary and the manifest) and `stage_de_manifest.json` recording: the input files with sha256, the pseudobulk and export manifest hashes, the primary analysis's `result.json` hash and `delta` status, the models, covariates, reference levels and α, the feature counts per set, the concordance cells, findings, pydeseq2 and dependency versions, git revision and timestamps. A rerun on the same inputs and versions SHALL give equal result tables (equal `pvalue` and `padj`, equal counts) and a byte-identical summary. The step SHALL write nothing under the pseudobulk, export or analysis directories, and no `export`, `analyze` or sensitivity step SHALL read its outputs.

#### Scenario: Resume
- **WHEN** a model's fit record matches the inputs and settings now on disk and its tables hash as recorded, and the step runs without `--force`
- **THEN** the model is skipped, and a record made from other inputs or settings stops the step with an error naming it, writing nothing

#### Scenario: Rerun
- **WHEN** `stage-de --cohort cohort_30` is run twice with `--force`
- **THEN** the two `stage_de_summary.csv` files are byte-identical and the result tables agree on every `padj`

#### Scenario: Inputs untouched
- **WHEN** the step finishes
- **THEN** the sha256 of every pseudobulk, export and analysis file equals its value before the run

### Requirement: Committed report with the reading rules
The change SHALL commit `reports/motco-astrocyte-stage-de-<date>.md` whose numbers are pasted from the generated files and checked against them before commit. It SHALL state: the stage counts describe the stage signal in the pseudobulks; the interaction counts and concordance cells cross-check the primary analysis and are not a second test of its hypothesis; a discordance in either direction is reported as such, since per-feature FDR over thousands of features and one latent-space permutation test ask related, not identical, questions; the union counts carry no joint error-rate claim; the primary analysis stays primary with its non-rejection rules; the covariate set and the feature set beside the counts.

#### Scenario: Report numbers checked
- **WHEN** the report is committed
- **THEN** every count and concordance cell in it equals the value in `stage_de_summary.csv`, and the `delta` statuses equal those in the primary `result.json` records

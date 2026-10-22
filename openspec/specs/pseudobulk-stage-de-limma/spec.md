# pseudobulk-stage-de-limma Specification

## Purpose

Fits the registered stage-DE models a second time with limma-voom: TMM normalization, quality weights, moderated t and joint moderated F tests. It runs on the same donors, counts, designs and feature set as the DESeq2 step. Its concordance cell sits beside DESeq2's and it compares the two methods. It is a second predeclared method beside DESeq2, never a replacement, and it changes no MOTCO reading.

## Requirements

### Requirement: Inputs include the DESeq2 stage-DE outputs
For each requested cohort the step SHALL read only the pseudobulk, primary export and primary analysis files the DESeq2 stage-DE step reads, and that step's outputs for the cohort (manifest, summary, per-coefficient tables). A missing file SHALL stop the step with an error naming the file and the command that produces it. Nothing SHALL be written.

#### Scenario: DESeq2 outputs absent
- **WHEN** `results/stage_de/cohort_30/stage_de_manifest.json` is missing
- **THEN** the step stops with an error naming the file and the `stage-de` command, and writes nothing

#### Scenario: Analysis absent
- **WHEN** `results/motco_analysis/cohort_50/result.json` is missing
- **THEN** the step stops with an error naming the file and the `analyze` command, and writes nothing

### Requirement: References must match, as in the DESeq2 step
Before any fit the step SHALL apply every reference check of the DESeq2 stage-DE step and its reproduction of the export from the counts: the pseudobulk manifest hash, donors, sex and stage, k and cells, feature statuses, the analysis's export hashes, MOTCO version, seed and permutation count. On any difference it SHALL stop naming the file or setting and the command to rerun, writing nothing.

#### Scenario: Pseudobulk changed since the export
- **WHEN** `atac_counts.h5ad` no longer reproduces the feature statuses, library sizes, TMM factors, tested-feature variances or exported values the export recorded
- **THEN** the step stops with an error naming the file and the `export` command, and writes nothing

### Requirement: DESeq2 outputs bound to the same inputs
Before any fit the step SHALL check:
- the DESeq2 manifest is complete for both omics;
- it records the same sha256 for every input file as the files now on disk;
- every DESeq2 result table hashes as that manifest records;
- the DESeq2 summary's counts and concordance cells equal a recomputation from those tables.

Otherwise it SHALL stop naming the file and the `stage-de` command, writing nothing. It SHALL never write under the DESeq2 directory.

#### Scenario: DESeq2 run made from other inputs
- **WHEN** the DESeq2 manifest of `cohort_30` records a `rna_counts.h5ad` sha256 other than the file's now on disk
- **THEN** the step stops with an error naming the manifest and `stage-de --force`, and writes nothing

#### Scenario: DESeq2 table edited
- **WHEN** a DESeq2 per-coefficient table no longer hashes as its manifest records
- **THEN** the step stops with an error naming the table, and writes nothing

### Requirement: Same models, designs and feature set as the DESeq2 step
The step SHALL fit exactly the registered stage-DE models:
- the same formulas, tested coefficients, reference levels (`Female`, stage `0`) and per-omic covariates;
- each model's design matrix equal, column by column under the same names, to the one the DESeq2 step fits;
- the features with export status `variance` or `kept`, in the export's order, each row carrying `exported` (status `kept`).

It SHALL select no feature from sex, stage, covariates or results.

#### Scenario: Tested set equals the DESeq2 tested set
- **WHEN** the ATAC models of `cohort_30` are fitted
- **THEN** every result table lists the 128,916 abundance-passing peaks of the DESeq2 tables in the same order, 5,000 of them `exported`

#### Scenario: Treatment coding agrees
- **WHEN** a model's design is rebuilt from the design table and formula by the R side's own model-matrix construction
- **THEN** it equals the design matrix passed in, column by column under the name map (`sex[T.Male]` ↔ `sexMale`, `stage[T.1]` ↔ `stage1`), or the step stops naming the model

### Requirement: Predeclared normalization bound to the export
Each omic's counts SHALL be TMM-normalized on the tested features, with library sizes equal to their column sums and the export's TMM trims. Before any model is fitted, the factors SHALL equal the export manifest's per-donor TMM factors to an absolute 1e-8. Otherwise the step SHALL stop naming the omic, the largest difference and the `export` command, writing nothing for that cohort.

#### Scenario: TMM factors differ
- **WHEN** a donor's TMM factor computed by the step differs from the export manifest's by 3e-6
- **THEN** the step stops with an error naming the omic, the difference and the `export` command

### Requirement: Predeclared limma-voom fit with quality weights
Each model SHALL be fitted by voom with per-donor quality weights on the TMM-normalized counts and that model's design, then a linear model fit and empirical-Bayes moderation, all with the method's default settings. The settings, the R, limma and edgeR versions, and each donor's quality weight per model SHALL be recorded. None SHALL be a runtime option or be changed from results.

#### Scenario: Quality weights recorded
- **WHEN** the RNA `factor_interaction` model of `cohort_30` is fitted
- **THEN** its fit record lists one quality weight per donor, with the weights' minimum, median and maximum overall and their mean per sex × stage cell

### Requirement: Moderated t per coefficient and joint moderated F per factor model
Per tested coefficient the step SHALL report a moderated t-test. Per factor model it SHALL report one moderated F-test of its two coefficients: `stage_joint` for the additive factor model, `interaction_joint` for the interaction factor model. Adjusted p-values SHALL be Benjamini–Hochberg over all tested features, with no independent filtering. A feature SHALL count as significant when its adjusted p-value is ≤ 0.05.

#### Scenario: Tables per cohort and omic
- **WHEN** the step finishes `cohort_30`, omic `rna`
- **THEN** `results/stage_de_limma/cohort_30/rna/` holds six per-coefficient tables (`feature`, `AveExpr`, `logFC`, `SE`, `t`, `P.Value`, `adj.P.Val`, `B`, `exported`) and two joint tables (`feature`, `AveExpr`, the two coefficients' `logFC`, `F`, `P.Value`, `adj.P.Val`, `exported`)

#### Scenario: No missing adjusted p-values
- **WHEN** any table is written
- **THEN** every tested feature has an adjusted p-value

### Requirement: Joint F approximation measured
Per joint test the step SHALL also compute, for each feature, the exact moderated F from the feature's own weighted coefficient covariance, with the same variance and degrees of freedom. It SHALL record the largest relative difference from the reported F. Features whose call at α differs when the exact F's p-values are BH-adjusted SHALL be counted and listed as findings. The reported F SHALL remain limma's.

#### Scenario: Calls unchanged
- **WHEN** the exact and reported F give the same calls at α for `interaction_joint`
- **THEN** the fit record gives the largest relative difference and a changed-call count of 0, and no finding is written

### Requirement: Summary on two feature sets with union rows
The step SHALL write `stage_de_limma_summary.csv` with a row per cohort × omic × model × test × feature set (`exported`, `abundance`). Each row SHALL give the number tested, the number significant and the smallest adjusted p-value. The tests SHALL include each moderated t, each joint F and, per factor model, the union of its two moderated t calls, labelled as a union with no joint error-rate claim. The exported-set counts SHALL be the primary reading.

#### Scenario: Summary rows
- **WHEN** both omics of a cohort have run
- **THEN** the summary has, per omic and feature set, six moderated-t rows, two joint-F rows and two union rows, and the exported `n_tested` is 5,000 in every row

### Requirement: Limma concordance rule beside the DESeq2 cell
For each cohort and omic the step SHALL set `interaction_detected` (at least one exported feature with adjusted p ≤ α on `interaction_joint`). It SHALL compare that with the primary analysis's `delta` rejection status, read as the DESeq2 step reads it. The cell SHALL be `concordant` when both agree and `discordant` otherwise, a discordance written as a finding. The rule SHALL be a constant, read no p-value magnitude of MOTCO's and change no MOTCO output.

#### Scenario: Both negative
- **WHEN** no exported RNA feature of `cohort_30` has an `interaction_joint` adjusted p ≤ 0.05 and the primary analysis did not reject `delta`
- **THEN** the limma cell for `cohort_30`, `rna` is `concordant` with `interaction_detected = false`, `delta_rejected = false`

### Requirement: Method agreement of the concordance cells
Per cohort and omic, the step SHALL report the DESeq2 cell (read from the checked DESeq2 summary) beside the limma cell, with `methods_agree` (true when both `interaction_detected` values are equal). A disagreement SHALL be written as a method-dependence finding naming both counts. Neither cell SHALL replace the other, and the DESeq2 cell SHALL stay the reported cell of the DESeq2 step.

#### Scenario: Methods disagree
- **WHEN** limma finds 3 exported features on `interaction_joint` in `cohort_50` ATAC and DESeq2's union count there is 0
- **THEN** `methods_agree` is false, the manifest records a method-dependence finding naming 3 and 0, and the DESeq2 outputs and the primary analysis are unchanged

### Requirement: Method comparison table
The step SHALL write `method_comparison.csv` per cohort × omic × model × test × feature set. Each row SHALL give:
- the significant count per method;
- the features called by both methods and by one only;
- whether each method detects anything, and whether that agrees.

Per-coefficient rows SHALL add the Spearman correlations of the two methods' log2 fold changes and test statistics, with the number of features compared. The table SHALL carry no p-value of its own.

#### Scenario: Per-coefficient comparison
- **WHEN** both methods have run for `cohort_30` RNA
- **THEN** `method_comparison.csv` has a row for `factor_additive` / `stage_2_vs_0` on each feature set with both counts, the overlap, the one-method counts and the two Spearman correlations

### Requirement: Joint tests compared with the DESeq2 union, disagreements as findings
Each joint-F row of the comparison SHALL compare limma's joint-F set with DESeq2's union set of the same model, labelled as a comparison of different tests. Every row whose detection differs between the methods SHALL be written as a method-dependence finding.

#### Scenario: Joint against union
- **WHEN** the `interaction_joint` row is written
- **THEN** it compares limma's joint-F set with DESeq2's union set and is labelled `joint_f_vs_union`

#### Scenario: Detection differs
- **WHEN** DESeq2 detects 2 exported stage features on `factor_additive` / `stage_2_vs_0` in `cohort_50` RNA and limma detects 0
- **THEN** the row has `detection_agrees = false` and the manifest lists a method-dependence finding naming the row

### Requirement: R execution isolated and recorded
The R side SHALL run as a separate process from one packaged script with one BLAS thread. The script's sha256 and the R, limma and edgeR versions SHALL be recorded in every fit record and in the manifest. If R is unavailable, a required package is missing or the process fails, the step SHALL stop naming the cause and the model. It SHALL write no table or record for that model.

#### Scenario: R missing
- **WHEN** the step runs outside the Nix shell, where `Rscript` is not on the path
- **THEN** it stops with an error naming `Rscript` and the Nix shell, and writes nothing

#### Scenario: R process fails
- **WHEN** the R process exits non-zero during the ATAC `trend_additive` fit
- **THEN** the step stops with the end of R's error output, and no `trend_additive` table or fit record is written for ATAC

### Requirement: Outputs, manifest and reproducibility
Per cohort the step SHALL write the tables, a fit record per omic × model (fingerprint of inputs, settings, script hash and versions; the fit; the tables' sha256), the summary, `method_comparison.csv`, `stage_de_limma.md` (rendered from them) and `stage_de_limma_manifest.json`. A rerun on the same inputs and versions SHALL give byte-identical tables, summary and comparison. The step SHALL write nothing under its input directories. No export, analysis or sensitivity step SHALL read its outputs.

#### Scenario: Resume
- **WHEN** a model's fit record matches the inputs, settings, script and versions now present and its tables hash as recorded, and the step runs without `--force`
- **THEN** the model is skipped, and a record made from other inputs, settings, script or versions stops the step with an error naming it, writing nothing

#### Scenario: Rerun
- **WHEN** `stage-de-limma --cohort cohort_30` is run twice with `--force`
- **THEN** the two `stage_de_limma_summary.csv` files, the two `method_comparison.csv` files and every table are byte-identical

#### Scenario: Inputs and DESeq2 outputs untouched
- **WHEN** the step finishes
- **THEN** the sha256 of every pseudobulk, export, analysis and DESeq2 stage-DE file equals its value before the run

### Requirement: Committed comparison report
The change SHALL commit `reports/motco-astrocyte-stage-de-limma-<date>.md`. Its numbers SHALL be pasted from the generated files of both methods and checked against them before commit.

#### Scenario: Report numbers checked
- **WHEN** the report is committed
- **THEN** every count, concordance cell and comparison value in it equals `stage_de_limma_summary.csv` or `method_comparison.csv`, every DESeq2 number equals `stage_de_summary.csv`, and the `delta` statuses equal the primary `result.json` records

### Requirement: Report reading rules
The report SHALL state:
- the proposal was written after the DESeq2 results were known;
- the DESeq2 results stay the reported ones;
- a call that differs between methods is method dependence, not a result;
- a zero count means "none detected at this n";
- the union counts carry no joint error-rate claim, and the joint F's claim is approximate;
- the concordance cells cross-check the primary analysis, whose non-rejection rules are untouched.

#### Scenario: Method dependence stated
- **WHEN** a feature is significant under limma and not under DESeq2
- **THEN** the report lists it as method dependence and reads it as a result of neither method

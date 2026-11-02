# case-study-figures Specification

## Purpose

Draws the case study's registered figures from the generated outputs of both cell types, with generated captions and a manifest that ties every figure to the files it was made from, so that a figure shows exactly the current results and visibly goes stale when they are re-issued.

## Requirements

### Requirement: Registered figure set
The `figures` command SHALL draw exactly the figures of the registry `FIGURES`, in registry order, each written as `<nn>_<name>` (`01_cohorts` … `10_stage_de`). No flag SHALL select, skip or add a figure. Entries SHALL be added or changed by a change proposal only.

#### Scenario: Default run
- **WHEN** `sead-sexdiff figures` runs with every input present
- **THEN** ten figures `01_cohorts` to `10_stage_de` are written, each as PNG, PDF and SVG, with one caption each in `figures.md`

### Requirement: Generated outputs only
The figures SHALL read only files written by the package's steps under the default directories of the registered cell types (`cell_type_dirs`, plus `results/motco_sensitivity/` for the astrocyte preprocessing variants). They SHALL NOT read raw donor objects, any path under `results/superseded/`, or any file a step did not write.

#### Scenario: Superseded directory present
- **WHEN** `results/superseded/mtg-residual-correction-primary/motco_analysis/` exists beside `results/motco_analysis/`
- **THEN** only `results/motco_analysis/` is read, and the manifest lists no path under `results/superseded/`

### Requirement: Every registered analysis drawn
Figures that draw a family of analyses SHALL enumerate them from the registries: forced ranks from `ANALYSIS_FORCED_RANKS`, preprocessing variants from `PREPROCESSING_VARIANTS`, arms from `LATENT_SENSITIVITY_ARMS`, exclusion variants from `DONOR_EXCLUSION_VARIANTS`, bases from `POSITIVE_CONTROL_BASES`, cohorts from `EXPORT_COHORTS`. A registered analysis whose output is missing SHALL stop the run naming it, with nothing written.

#### Scenario: Overview of all analyses
- **WHEN** the p-value overview is drawn from today's outputs
- **THEN** it holds 44 analyses (astrocyte: primary, six forced ranks, eight preprocessing variants, two arms, one exclusion variant, per cohort; microglia: primary, two arms, one exclusion variant, per cohort) and 132 p-values

#### Scenario: Variant output missing
- **WHEN** `results/motco_sensitivity/keep_chrx/motco_analysis/cohort_50/result.json` does not exist
- **THEN** the run stops naming that path, and no figure, caption or manifest is written

### Requirement: Cell type of every manifest checked
Every manifest the figures read SHALL record the cell type of the root it was read from, checked as every other reader checks it. A manifest of another cell type SHALL stop the run naming the file and both cell types.

#### Scenario: Microglia manifest under the astrocyte root
- **WHEN** `results/motco_analysis/cohort_30/analysis_manifest.json` records `microglia`
- **THEN** the run stops naming the file, `astrocyte` and `microglia`, with nothing written

### Requirement: Recorded hashes hold
Before anything is drawn, every file that a manifest read by the figures records with a path and a sha256 (its inputs, references and outputs) SHALL exist and hash as recorded. A missing or differing file SHALL stop the run naming the manifest, the file and both hashes, with nothing written.

#### Scenario: Today's outputs
- **WHEN** the check runs on the outputs of `mtg-residual-correction-primary`
- **THEN** every recorded entry hashes as recorded and drawing starts

#### Scenario: Diagnostic made from an earlier analysis
- **WHEN** `results/motco_analysis/cohort_30/result.json` was re-made after `results/motco_diagnostics/cohort_30/diagnostics_manifest.json` recorded its sha256
- **THEN** the run stops naming the diagnostic manifest and the analysis file, so no figure pairs a diagnostic with an analysis it was not made from

### Requirement: One measurement across the analyses drawn
Every MOTCO analysis record drawn SHALL carry the same MOTCO version and the registered measurement: the residual correction `ANALYSIS_RESIDUAL_CORRECTION`, the rank floor `ANALYSIS_RANK_FLOOR` (or a forced rank of `ANALYSIS_FORCED_RANKS` in `rank-sensitivity`), `ANALYSIS_PERMUTATIONS` and `ANALYSIS_SEED`. A record that differs SHALL stop the run naming it and the differing setting.

#### Scenario: Re-issue half finished
- **WHEN** the microglia primary records MOTCO 0.14.0 and the astrocyte primary MOTCO 0.13.0
- **THEN** the run stops naming both records and versions, and no figure mixes the two

#### Scenario: Uncorrected analysis
- **WHEN** an arm's `result.json` has no `params.residual_correction`
- **THEN** the run stops naming the file and the correction expected (`df`)

### Requirement: Recomputed values reproduce their source
Every value a figure computes that an existing output already holds SHALL be checked against that output before it is drawn. This covers the covariate-adjusted donor scores' cell means against `stage_means.csv` (to 1e-8), the Spearman ρ against `axis_correlations.csv` (to 1e-9), and the DESeq2/limma overlap counts against `method_comparison.csv` (exact). A failed check SHALL stop the run naming the value and both numbers.

#### Scenario: Features dropped by independent filtering
- **WHEN** 2,159 exported microglia features have no DESeq2 `padj` and the figure would count limma-only calls among the rest
- **THEN** the count (8) differs from `method_comparison.csv` (10) and the run stops; with missing `padj` read as not significant the counts are 43 / 10 / 18 and drawing continues

### Requirement: Captions in three parts
Each figure SHALL have a caption made of a fixed description of its construction, sentences whose values are read from the files at run time, and the registered reading rules that apply to it. A caption SHALL contain no other interpretation of a result. Captions SHALL be written to `figures.md` and to the manifest.

#### Scenario: First-axis caption
- **WHEN** figure 08 is drawn from today's outputs
- **THEN** its caption gives the astrocyte ρ with stage (0.49) and with the ribosomal mean fraction (−0.84), and the microglia ρ (0.69, mitochondrial pseudobulk fraction 0.84), read from `axis_correlations.csv`, and states that they are measurement-space findings that do not reinterpret the primary analysis

#### Scenario: Value changes after a re-issue
- **WHEN** a re-issued primary has a different `delta` p-value
- **THEN** the regenerated caption gives the new value with no edit to its text

### Requirement: Binding reporting context drawn with MOTCO results
Every figure that draws a MOTCO statistic SHALL show, in the figure itself, the covariate block, the feature filtering, the RRPP permutations and residual correction, and the non-rejection rules. Every `angle` SHALL be shown with its CV-selected rank, and in figure 04 with its eigengap and null summary. The `latent_adjusted` arm's `angle` and `shape` p-values SHALL be marked as reported, not read.

#### Scenario: Overview figure
- **WHEN** figure 05 is drawn
- **THEN** each row shows the selected or forced rank per cohort, the `latent_adjusted` `angle` and `shape` markers are drawn as not read, and the footnote states that sensitivity analyses are read as dependence and that `angle` / `shape` non-rejections are uninformative

### Requirement: Figures manifest
The run SHALL write `figures_manifest.json`. Per figure it records every file read with its sha256, the manifests read, the MOTCO version and measurement of the analyses drawn, the caption, and the sha256 of each image. For the run it records the module's sha256, the package version, the git revision, the matplotlib version and the display constants.

#### Scenario: Tracing a figure to its inputs
- **WHEN** a reader asks which analysis figure 04 drew for microglia `cohort_50`
- **THEN** the manifest's entry for `04_statistics_vs_null` lists `results/microglia/motco_analysis/cohort_50/result.json` with the sha256 that file has

### Requirement: Deterministic files
Running the command twice on the same inputs with the same code and matplotlib version SHALL produce byte-identical PNG, PDF and SVG files. No creation date or random identifier SHALL enter a file. SVG text SHALL be stored as outlines and the font SHALL be the one bundled with matplotlib.

#### Scenario: Rerun on the same machine
- **WHEN** the command is run twice on unchanged outputs into two directories
- **THEN** every image file in one is byte-identical to its namesake in the other

### Requirement: Dense layers rasterized
In the PDF and SVG files, a scatter layer with more than `FIGURE_RASTER_MIN_POINTS` (500) points SHALL be embedded as an image at the PNG resolution. Axes, text, lines and smaller layers SHALL stay vector.

#### Scenario: Stage-DE figure
- **WHEN** figure 10 is written as SVG
- **THEN** the volcano and method-comparison point clouds are embedded images, the highlighted features and every label are vector, and the file is under 1 MB

### Requirement: Fingerprint embedded in every file
Each image file SHALL carry in its metadata the figure name and a fingerprint: the sha256 over the figure's input entries, the module's sha256 and the matplotlib version. The git revision SHALL NOT enter the file.

#### Scenario: Figure pasted into a slide deck
- **WHEN** a PNG is found outside the repository
- **THEN** its metadata gives the figure name and fingerprint, which match one entry of a `figures_manifest.json`

### Requirement: Staleness check
`figures --check` SHALL compare the existing manifest with the files on disk. A figure is stale when a file it read, or a file recorded with a sha256 by a manifest it read, hashes differently or is missing, or when the module's sha256 or the matplotlib version differs. It SHALL list the stale figures with reasons, exit 0 when all are fresh, 1 when any is stale and 2 without a manifest, and write nothing.

#### Scenario: Primary re-issued
- **WHEN** `results/motco_analysis/` was moved to `superseded/` and re-made since the figures were drawn
- **THEN** `--check` lists figures 02–10 with the changed or missing files (figure 10 through the stage-DE manifests' recorded reference), leaves figure 01 fresh, and exits 1

### Requirement: Working copy redrawn when stale
The output directory (default `results/figures/`) SHALL be a working copy. A run SHALL do nothing when every figure is fresh, and SHALL redraw the whole set when any figure is stale or there is no manifest. The set SHALL be written to `<out>.tmp/` and renamed into place only after every file and the manifest are written. On an error the temporary directory SHALL be removed and the previous set left unchanged.

#### Scenario: Fresh set
- **WHEN** a run finds every figure fresh
- **THEN** it prints that the set is fresh and writes nothing

#### Scenario: Stale set
- **WHEN** a run finds figures 02–10 stale after a re-issue
- **THEN** it lists them with their reasons, passes the consistency check on the re-issued outputs, and replaces the whole set, figure 01 byte-identical to its previous file

#### Scenario: Error while drawing
- **WHEN** a check fails while figure 07 is being drawn
- **THEN** `<out>.tmp/` is removed and the previous `<out>/` is unchanged

### Requirement: Frozen set for a report
`figures --report <report-stem>` SHALL run the checks and redraw the working copy if it is stale or absent. It SHALL then copy every SVG, `figures.md` and `figures_manifest.json` to `reports/figures/<report-stem>/`, check each copy's sha256 against the manifest, and never overwrite an existing destination. A committed frozen set SHALL never be regenerated or edited.

#### Scenario: Freezing for the figures report
- **WHEN** `figures --report motco-case-study-figures-2026-11-03` runs
- **THEN** `reports/figures/motco-case-study-figures-2026-11-03/` holds ten SVGs, `figures.md` and the manifest, each SVG hashing as the manifest records

#### Scenario: Destination exists
- **WHEN** that directory already exists
- **THEN** nothing is copied and the command stops naming it

#### Scenario: Unchanged figure in a later report
- **WHEN** a later report freezes a set in which figure 01's inputs, the module and matplotlib are unchanged
- **THEN** its `01_cohorts.svg` is byte-identical to the earlier report's copy

### Requirement: Every report carries a frozen set
Every report committed under `reports/` from this change on SHALL carry the set frozen by `figures --report <its stem>` from the outputs on disk when it is committed. Its check script SHALL verify, before commit, that the frozen SVGs hash as the frozen manifest records and that every input the frozen manifest records hashes as recorded on disk. Reports committed before this change SHALL NOT be edited.

#### Scenario: Report of a later change
- **WHEN** a change adds a preprocessing variant and commits its report
- **THEN** the report's frozen set includes the variant in figure 05, made by `figures --report` after the variant's runs, and its check script passed before the commit

#### Scenario: Outputs changed after freezing
- **WHEN** an output is rerun after `figures --report` and before the commit
- **THEN** the check script finds a recorded input that no longer hashes as recorded, and the report is not committed until its uncommitted frozen set is removed and `figures --report` is run again

### Requirement: Fixed display choices
The display choices SHALL be constants recorded in the manifest, never flags. They are `cohort_30` for figures 03, 07, 08 and 10, `factor_additive` `stage_2_vs_0` on RNA for figure 10, and E3's rule over `POSITIVE_CONTROL_TECHNICAL_SOURCES` for figure 08's variable. They also cover `DIAGNOSTIC_CORRELATION_THRESHOLD` for the heatmap labels, the palette, the resolution and `FIGURE_RASTER_MIN_POINTS`. Every count of a figure 10 coefficient not drawn SHALL be given in its caption.

#### Scenario: Stage-DE caption
- **WHEN** figure 10 is drawn
- **THEN** its caption gives, per cell type, the exported RNA counts at FDR 0.05 for `stage_1_vs_0`, `stage_2_vs_0` and their union from `stage_de_summary.csv`, and limma's `stage_joint` count

### Requirement: Display only
The figures SHALL test nothing and SHALL NOT be an input of any export, analysis, sensitivity, diagnostic, stage-DE or positive-control step. Nothing drawn or captioned SHALL change a reading of the reports that made the outputs.

#### Scenario: Caption differs from a report
- **WHEN** a caption value differs from the value pasted in the report that made the output
- **THEN** the difference is recorded as a finding about that report with both values, and neither the output nor its reading changes

# Proposal

## Why

The case study has fourteen committed reports and no figures. Every result is a table pasted from a generated file. The trajectories, the nulls, the 44 MOTCO analyses, the technical axis and the positive control are hard to take in from text alone. Ten prototype figures were drawn on 2026-11-02 from the re-issued outputs of `mtg-residual-correction-primary`. They read only generated files and checked their recomputed numbers against the pipeline's own tables. That check caught a plotting bug before anyone saw the plot: pydeseq2 leaves `padj` empty for 2,159 exported microglia features, and dropping those rows miscounted the limma-only calls (8 instead of 10).

The figures need the same guarantees as every other output:
- they are made from the files on disk and say which ones;
- they become stale, visibly, when an upstream output is re-issued;
- they cannot mix superseded and current results;
- the copy a report cites is frozen with the report, as its pasted numbers are.

## What Changes

- **A new `figures` command** (`sead-sexdiff figures`), in a new module `src/sead_sexdiff/figures.py`. It draws the registry `FIGURES` of ten figures from the generated outputs of both cell types only:

  | # | Name | Shows |
  |---|---|---|
  | 01 | `cohorts` | Nuclei per donor against the cohort thresholds; donors per sex × stage cell, both cohorts and cell types |
  | 02 | `trajectories` | LS-mean trajectories per sex in `pls_0` × `pls_1`, with covariate-adjusted donors; the four primaries |
  | 03 | `components_by_stage` | LS-means per component by stage, with the geometry shares; `cohort_30` |
  | 04 | `statistics_vs_null` | Observed `delta`, `angle`, `shape` against their RRPP null quantiles; the four primaries |
  | 05 | `pvalue_overview` | Every MOTCO analysis outside `results/superseded/` × 3 statistics |
  | 06 | `rank_sensitivity` | Astrocyte forced ranks against their nulls, with the CV-selected rank |
  | 07 | `axis_correlations` | Spearman ρ of every component with every diagnostic variable; `cohort_30` |
  | 08 | `first_axis` | `pls_0` against stage and against E3's largest technical variable; `cohort_30` |
  | 09 | `positive_control` | E1–E3, control against reference, both bases and cohorts |
  | 10 | `stage_de` | Volcano (`factor_additive` `stage_2_vs_0`, RNA) and DESeq2 against limma; `cohort_30` |

  Analyses are enumerated from the registries (`ANALYSIS_FORCED_RANKS`, `PREPROCESSING_VARIANTS`, `LATENT_SENSITIVITY_ARMS`, `DONOR_EXCLUSION_VARIANTS`, `POSITIVE_CONTROL_BASES`). A missing input is an error, never a skipped panel.
- **A consistency check before anything is drawn.** It covers:
  - the cell type recorded by every manifest read;
  - every file a manifest records with a sha256 (inputs, references, outputs), which must hash as recorded (1,498 entries hash as recorded today);
  - one MOTCO version and the registered measurement (`ANALYSIS_RESIDUAL_CORRECTION`, `ANALYSIS_RANK_FLOOR` or a forced rank, `ANALYSIS_PERMUTATIONS`, `ANALYSIS_SEED`) across every analysis drawn;
  - reproduction of every number the figures recompute against the file that already holds it: adjusted donor means against `stage_means.csv`, Spearman ρ against `axis_correlations.csv`, overlap counts against `method_comparison.csv`.

  A failure names the file and writes nothing. So a half-finished re-issue cannot produce a figure that mixes two measurements.
- **Generated captions.** Each figure gets a caption in three parts:
  - a fixed description of the construction;
  - values filled in from the files at run time, never typed;
  - the registered reading rules.

  No caption interprets a result beyond those rules. The captions are written to `figures.md` and into the manifest. The figures carry the binding reporting context (covariate block, feature filtering, non-rejection rules) as an in-figure footnote where MOTCO results are drawn.
- **A manifest, `figures_manifest.json`.** Per figure it records:
  - every file read with its sha256, and the manifests read;
  - the MOTCO version and measurement of the analyses drawn;
  - the module's sha256, the package version and the git revision;
  - the matplotlib version and the display constants;
  - the caption and the sha256 of each image.
- **Deterministic files.** Each figure is written as PNG (200 dpi), PDF and SVG:
  - SVG text is stored as outlines;
  - the SVG id salt is fixed;
  - creation dates are left out;
  - the font is matplotlib's bundled DejaVu Sans;
  - scatter layers with more than `FIGURE_RASTER_MIN_POINTS` (500) points are rasterized inside the vector files.

  Two runs on the same inputs are byte-identical. The prototype measured 20 of 20 SVG and PDF files. Rasterizing takes the stage-DE SVG from 4.6 MB to 420 KB. A fingerprint of the inputs and the code is embedded in each file's metadata, so a copied figure stays traceable. matplotlib becomes a direct dependency, pinned exactly at the locked 3.10.7.
- **`results/figures/` is a working copy.** A plain run does nothing when every figure is fresh. Otherwise it redraws the whole set through `figures.tmp/`, renamed into place when complete. The record of a figure is the frozen copy a report carries (below), so the working copy is redrawn freely, needs no `--force`, and is gitignored.
- **`figures --check`.** It lists the figures whose inputs, module or matplotlib version changed since they were drawn, or whose inputs are gone (moved to `superseded/`). It exits non-zero and writes nothing.
- **Figures generated with every report.** `figures --report <report-stem>` is one command that:
  - runs the checks;
  - redraws the set if it is stale or absent;
  - copies the SVGs, `figures.md` and the manifest to `reports/figures/<report-stem>/`, never overwriting.

  It is a rule (AGENTS.md) that every change that commits a report ends with it. The report's check script verifies the frozen set against the outputs on disk at commit time. A frozen set is never regenerated, like a report's pasted numbers. Every report freezes the full set. Unchanged figures are byte-identical, so git stores them once and only changed figures add size. A change that adds a step adds its figure to `FIGURES` in the same proposal, so its report freezes it.
- **Migration.** Spec `motco-release-migration` gains:
  - `figures` as the last rerun step, its set frozen with the re-issue report;
  - the figures manifest in the version-consistency check.

  `results/figures/` is not moved to `superseded/`. It is redrawn, and the superseded figures are the ones frozen with earlier reports.
- **One committed report**, `reports/motco-case-study-figures-<date>.md`. It presents the first frozen set (`figures --report`) with its captions. Earlier reports are not edited. It states that the figures and captions show outputs already reported, that the readings are those of the reports that made each output, and that no figure is a new analysis.
- **Predeclared display choices.** They are fixed in this proposal, not taken from results:
  - `cohort_30` (the primary cohort) for figures 03, 07, 08 and 10, with `cohort_50` shown wherever a figure draws both cohorts;
  - `factor_additive` `stage_2_vs_0` on RNA for figure 10, the contrast between the extreme stages; its caption also gives the `stage_1_vs_0` and union counts from `stage_de_summary.csv`, so no count is hidden;
  - figure 08's technical variable chosen by E3's rule over `POSITIVE_CONTROL_TECHNICAL_SOURCES`;
  - the reporting flag `DIAGNOSTIC_CORRELATION_THRESHOLD` reused for the heatmap labels;
  - astrocytes only for figure 06, because `rank-sensitivity` runs on astrocytes only.

  The prototype drew the current results. None of these choices selects which results are shown: every analysis, both cell types and both bases are drawn.
- **Reading, predeclared.** The figures are a display step. They test nothing, feed no export, analysis, sensitivity or diagnostic step, and change no reading. A caption value that differs from the report that made the output is a finding about that report, recorded with both values. Every binding rule holds in the captions as in the reports:
  - CV rank and eigengap beside `angle`;
  - `angle` and `shape` non-rejections uninformative;
  - `delta` conservative under the `df` correction;
  - the `latent_adjusted` arm's `angle` and `shape` reported but not read;
  - sensitivity analyses read as dependence, never as results.

## Capabilities

### New Capabilities
- `case-study-figures`: the `figures` command. Covers the registered figure set drawn from generated outputs only, the consistency check before drawing, generated captions, the manifest, deterministic files, the working copy and its staleness check, and the frozen set every report carries.

### Modified Capabilities
- `motco-release-migration`: `figures` is the last step of the rerun order, and its set is frozen with the re-issue report. Its manifest is included in the version-consistency check.

## Impact

- **Code:**
  - new `src/sead_sexdiff/figures.py`, with the registry, constants, consistency check, ten figure functions, captions, manifest, check and report freeze;
  - `cli.py` gains the `figures` command;
  - `data.py` is unchanged apart from any helper made importable (`_git_revision`, `_sha256_file`, `check_manifest_cell_type` are reused, not copied).
- **Dependencies:** `matplotlib==3.10.7` as a direct dependency via `uv add`. It is locked today through scanpy, and no other locked package may move.
- **Outputs:** `results/figures/` (gitignored), with PNG, PDF and SVG per figure, `figures.md` and `figures_manifest.json`. Committed: `reports/figures/<report-stem>/` (about 2.1 MB of SVG per set) and the figures report.
- **Docs:**
  - AGENTS.md:
    - pipeline step 16;
    - the layout: the new module, an exception to "all logic in `data.py`" for plotting, approved by the user on 2026-11-02;
    - figure constants under "flag and ask";
    - the migration rerun order;
    - the rule that every committed report carries a frozen set made by `figures --report`, checked by its check script;
  - README: a "Case-study figures" section.
- **Not affected:** every existing step, output and manifest. The figures read them and write nothing outside `results/figures/` and `reports/figures/`.

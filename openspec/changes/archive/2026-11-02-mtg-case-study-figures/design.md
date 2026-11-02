# Design

## Context

- **The prototype.** A scratch script drew the ten figures on 2026-11-02 on `nixos-dell-laptop`, which holds both cell types' re-issued outputs (MOTCO 0.13.0, `df`, floor 2). The script is `scratchpad/figures.py`, a session directory, so it is ported, not imported.
  - It imports the registries and constants from `sead_sexdiff.data` and reads directories through `cell_type_dirs`.
  - It already makes three reproduction checks: adjusted cell means against `stage_means.csv`, ρ against `axis_correlations.csv`, overlap counts against `method_comparison.csv`. The third caught the `padj` bug.
  - Rendering was inspected figure by figure and label collisions were fixed: arrowheads over the stage markers, p-value labels, heatmap group labels, the E3 variable name.
- **What the manifests record.** Every manifest written by the MTG steps records `cell_type` (pre-registry manifests excepted, read as astrocyte), `software` and `git`. Its inputs, references and outputs are recorded as `{"path", "sha256"}` objects, nested at different depths per step. Walking every manifest under `results/` outside `superseded/` finds 1,498 such entries. All are repository-relative paths under `results/`, and all hash as recorded (checked 2026-11-02).
- **What MOTCO records.** `latent_scores.csv` holds raw PLS scores. `stage_means.csv` holds LS-means: `ls_means @ betas`, evaluated at `covariates.reference_row` (`pooled_mean`). `result.json` holds per-component `covariates.coefficients`, so adjusted donor scores are raw − (x − x_ref)·β. Their cell means equal `stage_means.csv` to 1.3e-13 for both cell types' `cohort_30`. `null_summary` holds the mean, sd and q50/q90/q95/q99 (no null draws).
- **Helpers in `data.py`.** `_git_revision`, `_sha256_file`, `check_manifest_cell_type`, `cell_type_dirs`, the registries and the constants. matplotlib 3.10.7 is locked as a transitive dependency of scanpy.
- **Patterns the steps share.** Output through `<dir>.tmp` renamed when complete, and removed on error. `[command]`-prefixed prints. Constants under "flag and ask". Readers that refuse rather than skip.

## Goals / Non-Goals

**Goals:**
- A figure's inputs are known by construction, not by a hand-written list, and are checked before the figure is trusted.
- The same inputs, code and matplotlib version produce the same bytes, so a hash comparison means something.
- A re-issue makes the old figures visibly stale and leaves the committed copies alone.

**Non-Goals:**
- Interactive or HTML figures, and journal styling (fonts, sizes, panel letters). Those can come from the same functions later.
- Figures of the legacy PFC pipeline or of superseded outputs.
- Byte-identical figures across machines. The inputs themselves differ across machines at rounding level (the DESeq2 tables do), so staleness is judged on input hashes, never on image hashes.
- A per-figure flag, a cell-type flag or directory flags for inputs. The figure set is the registry, and both cell types are read by role, as `positive-control` does.

## Decisions

**D1. A new module, `src/sead_sexdiff/figures.py`, and a `figures` command.**
- The module imports constants and helpers from `data.py` and never the other way round. `cli.py` wires `figures` with `--out` (default `results/figures`), `--check` and `--report <report-stem>`.
- *Alternative: inside `data.py`, per AGENTS.md's "all logic in `data.py`".* Rejected:
  - plotting is a separate concern with its own dependency and rcParams;
  - `data.py` is past 13,000 lines;
  - nothing in the pipeline imports the figures.

  The user approved the separate module on 2026-11-02. AGENTS.md's layout records it as the exception for plotting.

**D2. A registry with predeclared display constants.**
- `FIGURES` is a dict keyed by name, in drawing order. Each entry holds its number, title, function and caption builder.
- The display constants:
  - `FIGURE_COHORT = "cohort_30"`, tied at import to `POSITIVE_CONTROL_PRIMARY_COHORT`;
  - `FIGURE_STAGE_DE = {omic: "rna", model: "factor_additive", coefficient: "stage_2_vs_0"}`, checked at import against `STAGE_DE_MODELS`;
  - `FIGURE_RASTER_MIN_POINTS = 500`, `FIGURE_DPI = 200`, `FIGURE_SVG_HASHSALT = "sead-sexdiff"`;
  - `FIGURE_PALETTE`: sex blue `#2a78d6` / orange `#eb6834`; cell type aqua `#1baf7a` / violet `#6250d6`; ink, grid and status colours.
- Figure 08's variable rule reuses E3's: the largest |ρ| with `pls_0` over the variables whose source is in `POSITIVE_CONTROL_TECHNICAL_SOURCES`. Values come from the technical table, or from the export manifest for the TMM factors and library sizes.
- The prototype used only `donor_technical.csv`. Both cell types' choice is the same today (ribosomal mean, mitochondrial pseudobulk fraction), but the rule is E3's from now on.
- *Alternative: flags for cohort and coefficient.* Rejected. Every other display of results in this repository is fixed before it is read.

**D3. Reads go through a recording reader.**
- Every file a figure function reads goes through one reader object. It hashes the file once (cached), records `(path, sha256)` under the figure being drawn, and refuses any path that is absolute, outside `results/`, or under `results/superseded/`.
- Manifests are read through the same object, which also records them as manifests read.
- The manifest's input list is the reader's log, so it cannot omit a file the figure used.
- *Alternative: declare inputs per figure.* Rejected. A declaration drifts from the code, while the log cannot.

**D4. The consistency check: a global pre-pass, then per-figure reproductions.**
- *Pre-pass, before any drawing.* The set of manifests is derived from the registries: per cell type and cohort, the analysis, diagnostics, stage-DE, limma and pseudobulk manifests; the rank, preprocessing, latent and exclusion manifests and their analysis manifests; and the positive-control manifests per basis. Each is checked as follows:
  - `check_manifest_cell_type` against its root's cell type;
  - every nested `{"path", "sha256"}` entry exists under `results/` outside `superseded/` and hashes as recorded;
  - every `result.json` named under an analysis directory has `software.motco` equal to the astrocyte primary's;
  - `params.residual_correction == ANALYSIS_RESIDUAL_CORRECTION`;
  - `params.min_components == ANALYSIS_RANK_FLOOR`;
  - `params.forced_components` is `None`, or one of `ANALYSIS_FORCED_RANKS` under `rank-sensitivity`;
  - `params.permutations` and `params.seed` equal the constants.

  Any failure raises naming the manifest, file and values.
- *Per figure.* The three reproduction checks of the prototype, plus Spearman of figure 08's technical variable against `axis_correlations.csv`.
- A failure removes `<out>.tmp/`.
- *Alternative: check only the files the figures read.* Rejected. The chain (the diagnostic made from this analysis, the stage DE made from this analysis) is recorded only in the upstream manifests, so walking them is the only way to refuse a mixed set.

**D5. Captions are built, not written.**
- Each figure has a caption builder that receives the values the figure function loaded. It returns three parts:
  - the construction: fixed text, describing how the figure was made, not what it shows;
  - the values: f-strings over loaded numbers;
  - the reading: keys into `FIGURE_READING_RULES`, a dict whose texts restate the rules already in AGENTS.md (binding MOTCO rules, the `df` size statement, the `latent_adjusted` note, the sensitivity-as-dependence rule, measurement-space findings, stage-DE counts not being a second test, the positive control's basis rule).
- `figures.md` gets, per figure, `## Figure N. <title>`, an `![](NN_name.svg)` link, and the three parts as paragraphs.
- The captions are generated whenever the figures are drawn, so a re-issue regenerates them.
- *Alternative: hand-written captions in the report.* Rejected. They would carry retyped numbers, which AGENTS.md forbids.

**D6. Determinism.**
- rcParams: `font.family = "DejaVu Sans"` (bundled), `svg.fonttype = "path"`, `svg.hashsalt = FIGURE_SVG_HASHSALT`.
- `savefig` metadata:
  - PNG: `Title`, `Description` (the fingerprint);
  - PDF: `Title`, `Subject` (the fingerprint), `CreationDate=None`;
  - SVG: `Title`, `Description`, `Date=None`.
- Before each save, every `PathCollection` with more than `FIGURE_RASTER_MIN_POINTS` offsets gets `set_rasterized(True)`, rendered at `FIGURE_DPI`. This is done by walking `fig.axes`, not by patching `scatter`.
- Jitter uses `numpy.random.default_rng` with fixed seeds per figure.
- matplotlib is pinned with `uv add "matplotlib==3.10.7"`, the version already locked.
- The prototype showed SVG and PDF byte-identical across two runs. PNG is checked in the tasks.

**D7. The fingerprint.**
- It is the sha256 of a canonical JSON (sorted keys, no whitespace) of `{figure, inputs: [[path, sha256], …] sorted, module_sha256, matplotlib}`.
- The git revision is left out, so a README commit does not change any image. The manifest still records it.
- The module's sha256 stands in for code: a change to the drawing code changes every fingerprint, as intended. A change to `data.py` constants that matters reaches the figures through the outputs it changes.
- `--check` recomputes the same parts and compares them with the manifest.

**D8. Run modes.**
- **Plain run.** With a manifest at `--out` and every figure fresh, it prints and exits. Otherwise it redraws the whole set:
  - the pre-pass runs;
  - all ten figures are drawn into `<out>.tmp/` and the manifest is written last;
  - the old directory becomes `<out>.old/`, the new one moves into place, and `<out>.old/` is removed only after the rename succeeded.

  Redrawing the whole set rather than the stale figures keeps one manifest from one run. Unchanged figures come out byte-identical anyway.
- **No `--force`.** `results/figures/` is a working copy, and the record is the frozen set a report carries (D9). Replacing the working copy loses nothing that a report cites.
- **`--check`.** Runs only the staleness comparison, which also re-hashes the entries of the manifests the figures read (spec: a figure is stale through its manifests' recorded hashes). Exit codes 0 / 1 / 2, nothing written.
- **`--report <stem>`.** The one command a report needs:
  - the stem must match `^[a-z0-9][a-z0-9.-]*$`, and `reports/figures/<stem>/` must not exist (checked first, before any drawing);
  - it runs the plain run, which redraws if stale or absent;
  - it copies the ten SVGs, `figures.md` and `figures_manifest.json` to `reports/figures/<stem>/` through a temporary directory, verifies each SVG's sha256 against the manifest, and renames into place.
- *Alternative: draw at the end of each pipeline step.* Rejected. The steps run one cell type and one step at a time, and the set is consistent only after the last of them (D4), so a hook in each step would fail the check mid-chain or draw a mixed state. A report is the point where the outputs are known to be consistent.
- *Alternative: a git pre-commit hook refusing a report without a frozen set.* Not adopted (the user's decision). The AGENTS.md rule and the report's check script are the enforcement.

**D9. Every report carries a frozen set.**
- AGENTS.md's rule: every change that commits a report runs `figures --report <report stem>` after its last output run. The report embeds `figures/<stem>/NN_name.svg` and pastes each caption from the frozen `figures.md`, verbatim. The report's check script compares:
  - the pasted captions with `figures.md`;
  - the frozen SVG hashes with the frozen manifest;
  - the frozen manifest's input hashes, and the recorded hashes of the manifests it read, with the files on disk at commit time.
- Every report freezes all ten figures, not only the ones it discusses. That keeps one rule and one check, and a reader of any report sees the whole current state. Deterministic bytes make this cheap: git stores identical content once, so a report whose figures did not change adds only the manifest and `figures.md`.
- A change that adds a step adds its figure (when it has one) to `FIGURES` in its own proposal, so the step's report freezes it.
- Committed reports before this change are not edited. Their outputs are mostly superseded, and the command does not read `superseded/`, so it could not draw them anyway.

**D10. Migration.**
- AGENTS.md's re-issue text lists `figures --report <re-issue stem>` as the last step. `results/figures/` is not among the moved directories: it is redrawn, and the superseded figures are those frozen with the report that presented them.
- After a re-issue, the version-consistency check reads `figures_manifest.json`'s per-figure MOTCO versions.
- Figure 01 reads only the donor tables. Its byte identity with the previous report's frozen copy (same module and matplotlib) is a check of the same kind as the MOTCO-independent hashes.

## Risks / Trade-offs

- **Unseen data breaks the layout.** A re-issue may bring a longer variable name, a new registry entry or a count of 0 that collides with a label. → Every change looks at the redrawn figures before `--report` freezes them. The checks guarantee the numbers, not the layout.
- **A working-copy figure that was never frozen is lost on the next redraw.** → No report cited it. Anything worth keeping is frozen with a report.
- **The rule is not enforced by tooling.** A report could be committed without a frozen set. → The report's check script fails without one, and AGENTS.md states the rule beside the existing "numbers pasted, never retyped" rule, which is enforced the same way (a rule plus each report's check script).
- **A matplotlib or font change moves every byte.** → The exact pin. The fingerprint includes the matplotlib version, so `--check` reports every figure stale after an upgrade instead of hiding it.
- **The pre-pass rejects a legitimate partial state.** For example, a cell type whose latent arm has not been run yet. → That is intended: the registry says the arm exists, so the figures wait. The error names the missing output and the step that makes it.
- **Captions misstate a rule.** → The reading texts are written once in `FIGURE_READING_RULES` and reviewed in this change against AGENTS.md. The report's check pastes them verbatim. A rule changed in AGENTS.md by a later change must change them too, so the tasks list `FIGURE_READING_RULES` as a "flag and ask" constant.
- **Rasterized layers look soft when zoomed.** → 200 dpi, the PNG resolution. The marks that carry meaning (significant features, LS-means, labels) stay vector.
- **The figures read sex and stage.** → They are a display step and an input of nothing, as the diagnostics are (spec: Display only).
- **The aqua cell-type colour is below 3:1 contrast on white.** The validator returned WARN. → Every cell-type mark has a legend entry or direct label, never colour alone.

## Migration Plan

1. Pin matplotlib. Port the prototype into `figures.py` with D2–D8. Wire the CLI.
2. On `nixos-dell-laptop` (the holder of the re-issued outputs), run the following, which writes nothing in `results/` except `results/figures/`:
   - `figures`;
   - the determinism check (two runs into scratch directories);
   - the negative checks on a scratch copy of `results/`;
   - `--check`.
3. `figures --report motco-case-study-figures-<date>`. Write and check the report.
4. Update the docs.

Rollback: delete the module, the command and `results/figures/`. Nothing else reads them. The committed report and its frozen copies stay as a dated record.

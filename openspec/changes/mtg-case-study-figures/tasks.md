# Tasks

## 1. Dependency

- [x] 1.1 Pin matplotlib directly (design D6): run `uv add "matplotlib==3.10.7"`.
  - Verify that `pyproject.toml` lists it.
  - Verify that `uv.lock` still records matplotlib 3.10.7.
  - Verify that `uv export --no-hashes` before and after differs at most in the `# via` comments, with no package version changed.

## 2. Module core: registry, reader, consistency check

- [x] 2.1 Create `src/sead_sexdiff/figures.py` with the constants and the registry (design D2):
  - `FIGURES`, `FIGURE_COHORT`, `FIGURE_STAGE_DE`, `FIGURE_RASTER_MIN_POINTS`, `FIGURE_DPI`, `FIGURE_SVG_HASHSALT`, `FIGURE_PALETTE`, `FIGURE_READING_RULES`;
  - the import-time ties: `FIGURE_COHORT == POSITIVE_CONTROL_PRIMARY_COHORT`, and `FIGURE_STAGE_DE`'s model and coefficient present in `STAGE_DE_MODELS`.

  Verify that `nix develop -c uv run python -c "import sead_sexdiff.figures"` succeeds, and that changing either tie in a scratch copy of the module raises at import naming the constant.
- [x] 2.2 Add the recording reader (design D3). Verify in a scratch script that:
  - it logs `(path, sha256)` per figure;
  - it hashes each file once;
  - it refuses an absolute path, a path outside `results/` and a path under `results/superseded/`, naming the path.
- [x] 2.3 Add the pre-pass consistency check (design D4: manifests from the registries, cell types, recorded hashes, one measurement). Verify:
  - on the real outputs it passes, and reports the number of manifests and entries checked (at least the 1,498 entries hashing today);
  - on a full scratch copy of `results/` outside `superseded/` (`cp -r`, never hard links; run with the scratch directory as working directory), each of the following stops the run naming the manifest, file and values, with no `results/figures*` written:
    - one byte appended to `results/motco_analysis/cohort_30/latent_scores.csv`;
    - `params.residual_correction` removed from a latent arm's `result.json`;
    - `software.motco` changed in `results/microglia/motco_analysis/cohort_50/result.json`;
    - `cell_type` changed in an astrocyte diagnostics manifest;
    - `results/motco_sensitivity/keep_chrx/motco_analysis/cohort_50/` removed.

## 3. Figures and captions

- [x] 3.1 Port the ten figure functions from the prototype through the reader, with its layout fixes. Figure 08 uses E3's rule over `POSITIVE_CONTROL_TECHNICAL_SOURCES`, with TMM factors and library sizes read from the export manifest (design D2). Each reproduction check raises naming the value and both numbers (spec: Recomputed values reproduce their source). Verify on the real outputs that:
  - all ten figures render;
  - figure 05 holds 44 analyses and 132 p-values;
  - figure 08 picks `ribosomal_mean` (astrocyte, ρ −0.841) and `mitochondrial_pseudobulk` (microglia, ρ 0.844);
  - figure 10's microglia overlap counts are 43 / 10 / 18;
  - a scratch copy of `method_comparison.csv` with `n_only_limma` changed stops the run.

  Look at every figure, and record that no label collides with data or with another label.

  (Done 2026-11-02: all ten render; 44 analyses / 132 p-values; `ribosomal_mean` −0.841 and `mitochondrial_pseudobulk` 0.844; microglia 43 / 10 / 18; the changed `n_only_limma` stops the run naming both numbers, nothing written. Looked at every figure. Collisions inherited from the prototype were fixed: 03's component titles ran into each other (now three lines), 09's value labels overlapped for close pairs (they now diverge from the pair) and the long E3 name left its panel, and 10's volcano labels overlapped each other and dots (now a deterministic placer with leader lines). The one-line titles of 07 and 08 widened the page (now two lines). After the fixes no label collides with data or with another label. Figure 10 also checks the limma step's `inputs.deseq2` record against the DESeq2 files drawn, a binding the `{path, sha256}` walk does not cover.)
- [x] 3.2 Add the caption builders and `FIGURE_READING_RULES` (design D5), and write `figures.md`.
  - Verify with a scratch script that every number in each caption equals its source value in the file named in the manifest (p-values, ρ, counts, n, ranks).
  - Verify that figure 10's caption gives the `stage_1_vs_0`, `stage_2_vs_0`, union and limma `stage_joint` exported RNA counts per cell type from the summaries.
  - Verify that each reading text matches the rule it restates in AGENTS.md: the binding MOTCO rules, the `df` size statement, `latent_adjusted` not read, sensitivity as dependence, measurement-space findings, stage-DE counts not a second test, the positive control's bases.

  (Done 2026-11-02: an independent rebuild from the input files matches all 382 caption numbers; figure 10 names all four counts. Reading texts compared with AGENTS.md: `latent_adjusted_not_read` gained "neither a rejection nor a non-rejection is interpretable". `stage_de_counts` now says what AGENTS.md says: the stage counts describe the stage signal and the MOTCO non-rejection rules are untouched, since "not a second test" is AGENTS.md's rule for the interaction counts, which figure 10 does not draw. It also gained the Cook's note. `positive_control` gained the registered holds / does-not-hold reading.)
- [x] 3.3 Draw the binding reporting context in every figure that shows a MOTCO statistic (02, 04, 05, 06): covariate block, feature filtering, permutations and correction, non-rejection rules, the rank beside `angle`, and the eigengap and null summary beside `angle` in 04. Verify that the text is read from the records (covariate columns, blocks, permutations, correction) and not typed, by checking it against `result.json`.

  (Done 2026-11-02: in 02, 04, 05 (the primaries' context) and 06 every record's covariate columns, blocks, permutations, correction and κ appear in the drawn footnote; a reader returning altered records changes all of them.)

## 4. Files, manifest and run modes

- [x] 4.1 Make the saves deterministic and rasterize dense layers (design D6), and embed the fingerprint (design D7). Verify:
  - two runs into scratch directories give byte-identical files (all 30);
  - `10_stage_de.svg` is under 1 MB;
  - every PNG's text chunks and every PDF's and SVG's metadata carry the figure name and the fingerprint recorded in the manifest;
  - the git revision appears in no image file.
- [x] 4.2 Write `figures_manifest.json` (spec: Figures manifest). Verify that for every figure the input list equals the reader's log, with no path under `superseded/`. Verify that each image's sha256 matches, and that the run-level module sha256, git revision, matplotlib version and display constants are present. Verify that `04_statistics_vs_null` lists `results/microglia/motco_analysis/cohort_50/result.json` with its current sha256.
- [x] 4.3 Add the working-copy run (design D8): fresh means nothing written; otherwise a whole-set redraw through `<out>.tmp`, with rename and cleanup. Verify on the scratch copy:
  - a second run prints that the set is fresh and writes nothing (directory mtimes unchanged);
  - after the working copy's recorded module sha256 is changed in its manifest (stale through the code hash, with consistent outputs), the run lists every figure as stale, redraws the set and leaves no `.tmp` or `.old` directory;
  - every image is byte-identical before and after that redraw;
  - a failure injected in figure 07 leaves the previous set byte-identical and no `.tmp`.
- [x] 4.4 Add `--check` (spec: Staleness check). Verify on the scratch copy:
  - exit 0 on a fresh set;
  - exit 2 without a manifest;
  - after `results/motco_analysis/` is moved aside and re-copied with one changed byte in `result.json`: exit 1, figures 02–10 listed with reasons, figure 01 fresh;
  - after an edit to the module, every figure is listed as stale through the code hash;
  - nothing is written in any case.
- [x] 4.5 Add `--report <stem>` (spec: Frozen set for a report). Verify on the scratch copy, with `reports/` pointed at a scratch directory:
  - with a fresh set, the set is copied to `reports/figures/<stem>/` (ten SVGs, `figures.md`, manifest), each SVG hashing as recorded;
  - with a stale set, it redraws first and then copies;
  - a second `--report` to the same stem refuses before drawing;
  - an invalid stem (`../x`) refuses;
  - a second stem frozen from unchanged inputs gives byte-identical SVGs.

## 5. CLI and documentation

- [x] 5.1 Wire `figures` in `cli.py` (`--out`, `--check`, `--report`, `[figures]` prints, exit codes passed through). Verify that `sead-sexdiff figures --help` lists the three options and the defaults, and that `sead-sexdiff --help` lists the command.
- [x] 5.2 Document the command.
  - **README:** a "Case-study figures" section covering the ten figures, the inputs, the checks, the files, the working copy, `--check`, `--report` and the report embedding.
  - **AGENTS.md:**
    - pipeline step 16;
    - the layout entry for `figures.py`, recording the exception to "all logic in `data.py`" for plotting (user's decision);
    - the figure constants, incl. `FIGURE_READING_RULES`, under "flag and ask";
    - beside the `reports/` rule on pasted numbers: every committed report carries the set frozen by `figures --report <its stem>` after its last output run, checked by its check script;
    - a change that adds a step adds its figure to `FIGURES`;
    - `figures --report` as the last rerun step of a re-issue, and `results/figures/` as a working copy that is not moved;
    - the "display only" reading.

  Verify that every constant and flag named in both files exists in the code, and that the README commands run as written.

## 6. Real run, frozen set and report

- [x] 6.1 On `nixos-dell-laptop` (the holder of the re-issued outputs), run `sead-sexdiff figures`.
  - Verify that the consistency check passes and that `figures --check` exits 0.
  - Run the dataviz palette validator on the sex and cell-type pairs, and record the result (all pass; aqua below 3:1, relieved by legends and direct labels).
  - Look at all ten figures once more.

  (Done 2026-11-02 on nixos-dell-laptop: 156 manifests and 1,498 recorded entries hash as recorded, and the 44 analyses are on MOTCO 0.13.0 / df / floor 2; `--check` exits 0; the working copy is byte-identical to the scratch draws. Validator, light mode on #ffffff: the sex pair passes all checks (CVD ΔE 24.7 protan, normal 33.6). The cell-type pair passes (CVD 27.1 deutan, normal 32.8) with the expected WARN: #1baf7a at 2.82:1, relieved by the legend entries and direct labels of 04, 05, 09 and 10 (06 draws one cell type, named in its title). All ten figures looked at once more.)
- [x] 6.2 Run `figures --report motco-case-study-figures-<date>`. Verify that `reports/figures/motco-case-study-figures-<date>/` holds ten SVGs, `figures.md` and the manifest, and that the SVGs total about 2 MB.
- [x] 6.3 Write `reports/motco-case-study-figures-<date>.md`:
  - provenance: machine, git revision, MOTCO version and measurement of every analysis drawn, manifest sha256;
  - the ten frozen figures, each with its caption pasted from the frozen `figures.md`;
  - the reading: a display step, readings as in the reports that made each output (the latest being `motco-residual-correction-primary-2026-10-31.md`), no figure a new analysis;
  - the display choices and why each is fixed;
  - the prototype's `padj` miscount as the reason for the reproduction checks;
  - findings.

  (Done 2026-11-02: `reports/motco-case-study-figures-2026-11-02.md`. Its tables, counts and captions were filled by a scratch script from the frozen set and the files its manifest records, so no number is retyped.)
- [x] 6.4 Check the report with a scratch script. It verifies that:
  - every pasted caption equals the frozen `figures.md` text;
  - every frozen SVG hashes as the frozen manifest records;
  - every input the frozen manifest records hashes as recorded on disk;
  - every value shared by a caption and `motco-residual-correction-primary-2026-10-31.md` (primary statistics, p-values, ranks, CV AUROC, E1–E3 values, stage-DE counts) is equal.

  Any difference is recorded in the report's findings with both values. Commit only after the script passes, or after each difference is recorded.

  (Done 2026-11-02: 30 of 30 caption paragraphs pasted verbatim; 10 of 10 frozen SVGs hash as recorded; 212 recorded inputs and 498 entries of the manifests read hash as recorded on disk; 257 caption values equal the 2026-10-31 report, no rounding ties, no difference. PASS.)
- [x] 6.5 Copy `results/figures/` to the cluster's `results/figures/` (the canonical copy). Verify every file's sha256 against the laptop's.

  (Done 2026-11-02: the cluster had no `results/figures/`; rsync of the 32 files (30 images, `figures.md`, `figures_manifest.json`); every sha256 equal on both machines.)

## Workflow follow-up

- Open a PR with the code, docs, frozen figures and report; merge after review.
- Archive the change and sync the `case-study-figures` and `motco-release-migration` specs.
- Update the project memory (`motco-case-study.md`) with the merged change.

# Design

## Context

- **MOTCO 0.18.0** (tag `v0.18.0`, commit `bceead1`; PyPI wheel sha256 `b142a095c9fd68331bee69c70d64dbcf2c58e53cf37cf1be50f2bd1c8af951ed`, sdist `e583b40a…`, uploaded 2026-11-05).
  - **Between `v0.13.0` and `v0.18.0`** the package changed `analysis.py`, `cli.py`, `stats/attribution.py`, the new `stats/interaction.py` and `simulations/*`. `stats/permutation.py`, `stats/pls.py`, `stats/trajectory.py` and `stats/design.py` are unchanged, and `pyproject.toml`'s dependency pins are unchanged (scikit-learn `==1.8.0`).
  - **On the analysis path:**
    - attribution became available with covariates (0.14.0); it stays off by MOTCO's default, `AttributionDiagnosticSettings.enabled = False`, which the package never sets;
    - `TrajectoryAnalysisResult.latent_artifacts` (in memory only);
    - `span_split_test`.
  - **MOTCO's contract:** without the setting, every record and table equals 0.17.0's, and each release since 0.13.0 kept absent-setting records byte-identical.
- **The split test as `analyze_trajectories` runs it.**
  - After the fit, `_span_split_record` calls `interaction_split_test` on the standardized features in the PLS input coordinates against `qr(x_rotations_)`. It uses the trajectory models with the covariates, `permutations`, `seed` and the residual correction.
  - Its draws come from its own `np.random.default_rng(seed)`, RRPP's serial sequence, so RRPP's p-values do not move.
  - The record gains `span_split`:
    - `parts.{total,span,complement}`, each with `dimension`, `ss_interaction`, `rss`, `df_interaction` (2), `df_residual` (n − 12 with six covariates), `f`, `p_value_f` and `p_value_ss`;
    - `span_share`, `permutations`, `seed`, `residual_correction` (the string `"df"`) and `draws_match_rrpp` (`n_jobs in (None, 1)`).
  - `params.span_split_test` is present only when true.
  - With `covariate_latent_adjustment` set, the request raises before fitting.
- **Where the package builds analyses.** Five `TrajectoryAnalysisParams` constructions feed `_analyze_cohort`: `run_analysis`, `run_rank_sensitivity`, `run_preprocessing_sensitivity`, the arm registry in `run_latent_sensitivity`, and `run_exclusion_sensitivity`.
  - `_analyze_cohort` checks the residual correction between the MOTCO call and the first write. The split check goes in the same place.
  - `_sensitivity_row` is the base of every comparison row: `_variant_row` (preprocessing, exclusion), `_latent_row` and the rank sweep. The stage-DE concordance reads only its `delta_rejected`.
  - `_check_reference` serves the rank, preprocessing, latent and exclusion steps and `_stage_de_prepare`.
- **Current outputs.** The 2026-10-31 re-issue (MOTCO 0.13.0, `df`, floor 2) was made on `nixos-dell-laptop` and rsynced to the cluster and to `nixos-home` (sha256-verified 2026-11-02). `figures --report motco-case-study-figures-2026-11-02` froze the ten-figure set on `nixos-home`. The astrocyte `cohort_30` export hashes exactly as the input of MOTCO's plasmode and fidelity gate (`results/span-split-2026-11-05/plasmode-fidelity/fidelity.json`).
- **Primaries today:** astrocyte rank 5 / 2 (κ 1.0690 / 1.0753), microglia rank 2 / 4 (1.0731 / 1.0903), for `cohort_30` / `cohort_50`. So the `span` dimensions will be 5, 2, 2, 4, and the complement 10,000 minus those.

## Goals / Non-Goals

**Goals:**
- The test is set in one constant, passed explicitly, verified in every record and recorded wherever it acts. Every reader checks it, and the one exception is recorded rather than silent.
- The multiplicity policy and reading are fixed in code constants and in the specs before any real-label split test is run.
- After the re-issue every output outside `results/superseded/` rests on MOTCO 0.18.0, and every pre-existing value is shown to reproduce.

**Non-Goals:**
- A flag for the test, an arm without it, or a split test run outside `analyze_trajectories`.
- A null in-span share, `p_value_ss` readings, or MOTCO's null distributions (`return_null`).
- A stage-DE or positive-control rule that reads the split test.
- Bringing `nixos-dell-laptop` up to date; it re-syncs from the cluster when next used.

## Decisions

**D1. Pin `motco==0.18.0` from PyPI, bound to the tag.**
- **The move:** `uv add "motco==0.18.0"`. Then check:
  - `uv.lock` records the PyPI wheel with the sha256 above;
  - `uv export --no-hashes` differs from main's only in the motco line;
  - the downloaded wheel's `motco/` tree equals `git -C ../MOTCO archive v0.18.0 src/motco` file by file.
- **Why PyPI and the exact tag:** the existing pin rule (spec `motco-release-migration`).
- **Alternative rejected:** pinning 0.14.0–0.17.0 first. None of them carries the test, and each would need its own re-issue.

**D2. One constant, passed explicitly, checked in the record.**
- **Constants in `data.py`, beside the analysis constants:**
  - `ANALYSIS_MOTCO_MIN_VERSION = "0.18.0"` and `ANALYSIS_SPAN_SPLIT_TEST = True`;
  - the policy: `SPAN_SPLIT_PARTS = ("total", "span", "complement")`, `SPAN_SPLIT_CONFIRMATORY_PART = "complement"`, `SPAN_SPLIT_CONCORDANCE_PART = "span"`, `SPAN_SPLIT_STATISTIC = "p_value_f"` and `SPAN_SPLIT_NOT_READ = ("total", "span_share", "p_value_ss")`.
  - `SPAN_SPLIT_PARTS` is checked against MOTCO's `INTERACTION_SPLIT_PARTS` at the first analysis, beside `check_motco_rank_floor`. A MOTCO that renamed a part stops the run.
- **Where it is passed:** `span_split_test=ANALYSIS_SPAN_SPLIT_TEST` goes into all five constructions. The latent arm construction passes `span_split_test=ANALYSIS_SPAN_SPLIT_TEST and spec["latent_adjustment"] is None`, so `latent_adjusted` asks for it off. An import-time assertion ties the exception to `LATENT_SENSITIVITY_ARMS`: exactly the arms with a latent adjustment.
- **The record check in `_analyze_cohort`:** after the correction check and before any write. When requested, the record must have:
  - a `span_split` block with exactly `SPAN_SPLIT_PARTS`, each with a non-null `p_value_f`;
  - `permutations`, `seed` and `residual_correction` equal to `params`;
  - `parts.span.dimension == latent.selected_lv`, and `parts.total.dimension` equal to the summed block widths;
  - `parts.complement.dimension` equal to the total minus the span;
  - `df_interaction == len(STAGES) − 1`;
  - `draws_match_rrpp == (n_jobs == 1)`.

  When not requested, there must be no block. Any failure raises naming the cohort and the entry, with nothing written.
- **The manifest:** `measurement.span_split` is either `{"requested": true, "confirmatory": "complement", "concordance": "span", "statistic": "p_value_f", "not_read": [...]}` or `{"requested": false, "not_available": "<MOTCO's reason, quoted>"}`.
- **Alternative rejected:** letting `_analyze_cohort` drop the setting whenever a latent adjustment is present. That would hide the exception inside the helper instead of the registry that defines the arms.

**D3. Summary, CLI equivalent and comparison columns.**
- **`render_analysis_summary`** adds, after the statistics table and only when the record has a block, a "Split interaction test" table with the columns part, role, dimension, pseudo-F and p (pseudo-F). It is followed by the in-span share "(descriptive, not read)" and one line on the draws ("RRPP's permutation sequence at the analysis's seed" when `draws_match_rrpp`). A record without the block renders byte-identically to today, so superseded records and the `latent_adjusted` arm are unaffected.
- **`_cli_equivalent`** appends `--span-split-test` when `params.span_split_test`.
- **`_sensitivity_row`** gains `total_f`, `total_p`, `span_f`, `span_p`, `span_rejected`, `span_concordant`, `complement_f`, `complement_p`, `complement_rejected` and `span_share`, after the `shape` columns.
  - `span_concordant` = `span_rejected == any(delta_rejected, angle_rejected, shape_rejected)`.
  - Every cell is empty (`None`) when the record has no block.
  - There is no `total_rejected`, so no table invites reading `total`.
  - Every comparison table (rank, preprocessing, latent, exclusion) inherits the columns, and their Markdown renderings print them.

**D4. The reference check.** `_check_reference` adds `_record_span_split(record)`: `params.span_split_test` true and a block present. Otherwise it raises naming the file and `analyze`, like the floor and correction checks. All five readers inherit it.

**D5. Figure 11 and the figures checks.**
- **Registry:** `FIGURES` gains `span_split` (number 11, title "Split interaction test"), with reading keys `display_only`, `span_split`, `sensitivity_dependence` and `forced_rank`.
- **Panels.** All three are read through `FigureReader.result`:
  - (a) the four primaries (cell type × cohort) on a log p axis with α: complement filled in the cell-type colour, span hollow with its concordance label, total grey;
  - (b) the complement p of every registered analysis carrying the test, enumerated as figure 05 enumerates (`registered_analyses()`), `cohort_30` filled and `cohort_50` hollow, with the `latent_adjusted` rows labelled "not available";
  - (c) `span_share` against the rank for the same analyses.
- **Caption builder:** gives the primaries' three p-values, ranks, dimensions and shares, and the count of analyses drawn (40).
- **`FIGURE_READING_RULES["span_split"]`:** the policy, the complement reading (rejection and non-rejection), the fitted-span caveat, the share not read, and the size statement (measured on real noise for the astrocyte `cohort_30` export only).
- **Consistency checks:** `_measurement` gains `span_split_test` and `span_split_present`. `check_consistency` requires both for every analysis drawn except the `latent_adjusted` arms, which must have neither.
- **Fingerprint:** editing `figures.py` changes every figure's fingerprint. That is harmless here because the re-issue redraws all eleven.
- **Alternative rejected:** adding the complement to figure 05. That would change a frozen construction text and mix a fourth p-value into a figure read under the production statistics' rules.

**D6. Two gates on `nixos-home`, before the move.**
1. *Pin-move gate.*
   - Create a scratch worktree of the pre-change commit (main `a5e9008`) and run `uv sync --locked` there (0.13.0). Run `analyze --cell-type {astrocyte,microglia} --cohort all --out-dir <scratch>/v0130/<cell type>` on the real exports.
   - From the change branch, a scratch script runs `_analyze_cohort` with the analyze parameters and `span_split_test=False` into `<scratch>/v0180_none/<cell type>`.
   - Pass condition: the four CSVs are byte-identical per cohort. `result.json` and the manifest are equal after dropping `runtime_seconds`, `created_at`, `git`, `software.motco`, `cli_equivalent` and the new `measurement.span_split`.
   - This checks that nothing on the default path moved across 0.14.0–0.18.0.
2. *Kept-output comparison* (recorded, not a gate). Compare `v0180_none` with the kept primaries. They were made on the laptop, so rounding-level differences (one ulp in `np.log10`) are possible; the report records the largest and whether ranks, p-values and `analysis_metadata.csv` are equal.
3. *Contract gate.*
   - Run the change branch's `analyze --cell-type … --cohort all --out-dir <scratch>/v0180_split/<cell type>` and compare it with `v0180_none`.
   - Pass condition: the four CSVs are byte-identical, and `result.json` is equal after dropping `span_split`, `params.span_split_test` and `measurement.runtime_seconds`. That includes every p-value and null summary.
   - A failure stops the change and is reported to MOTCO.
   - `v0180_split` also serves as the reproducibility rerun of D8 step 1.

- *Why the gates run before the move:* a failed gate leaves `results/` untouched.

**D7. Superseded layout.**
- **The move:** `python3 scripts/move_superseded.py --change mtg-span-split-test` with the 2026-10-31 directory list.
  - Astrocyte: `motco_analysis`, `motco_rank_sensitivity`, `motco_diagnostics`, `motco_sensitivity`, `motco_latent_sensitivity`, `stage_de`, `stage_de_limma`, `motco_exclusion_sensitivity`.
  - Microglia, under `microglia/`: `motco_analysis`, `motco_diagnostics`, `motco_latent_sensitivity`, `stage_de`, `stage_de_limma`, `motco_exclusion_sensitivity`, `positive_control`.
  - It runs on `nixos-home`, then on the cluster. `results/figures/` is not moved.
- **Before the move:**
  - the sha256 of every file under both cell types' `pseudobulk/` and `motco_input/` and under `data/annotation/` is recorded;
  - so is every input hash recorded by `reports/figures/motco-case-study-figures-2026-11-02/figures_manifest.json`.
- **After the move:**
  - every input of that frozen set must be found by hash, either under `results/superseded/mtg-span-split-test/` or at its unchanged MOTCO-independent path;
  - so must every `result.json` hash cited by `reports/motco-residual-correction-primary-2026-10-31.md`.

  Those reports therefore stay resolvable without being edited.

**D8. Rerun order and commands.** These are the 2026-10-31 re-issue's commands (its design D6), through `nix develop -c uv run sead-sexdiff …` with `--n-jobs 1` and `--n-cpus 10` for `stage-de`:
1. `analyze`;
2. `rank-sensitivity` (astrocytes);
3. `axis-diagnostics`;
4. `preprocessing-sensitivity --variant all` with `axis-diagnostics` per variant root;
5. `latent-sensitivity` with its diagnostics;
6. `stage-de --omic both` and `stage-de-limma --omic both`;
7. `exclusion-sensitivity`, then on its variant root `axis-diagnostics`, `stage-de --omic both` and `stage-de-limma --omic rna`;
8. `positive-control` on both bases;
9. after every other step, `figures --report motco-span-split-<date>`.

The two cell types' chains run side by side.

**D9. Comparison script (superseded → new).** It runs over the 44 analysis pairs. For each pair:
- the D6.3 contract comparison against the superseded record: CSVs, and records outside `span_split`, `params.span_split_test`, `software` and runtimes. Byte-identity is expected only on the machine that made the superseded outputs, so on `nixos-home` it checks equality at rounding level with equal ranks, CV AUROCs, statistics, p-values and statuses;
- the new split block's values;
- the complement status and the span concordance cell.

The script also checks every p-free output listed in the 2026-10-31 design D7, at the same expectation:
- the axis-diagnostic CSVs, `adoption_check.csv` and the precedence;
- `mediator_check.csv`, `flag_check.csv` and `axis_check.csv`;
- E1–E3 on both bases;
- the stage-DE and limma tables and counts. These are beyond rounding across machines, with equal calls (finding 11 of 2026-10-26).

It produces the report's tables. Outside the D6 gates every difference is a finding, never a gate.

**D10. Version and measurement consistency check.**
- **Manifests:** every `*manifest*.json` outside `results/superseded/` that records a MOTCO version records `0.18.0`, `results/figures/figures_manifest.json` included.
- **Records:** every `result.json` carries `df`, floor 2 or a forced rank, and a `span_split` block. The exception is a `latent_adjusted` arm, which has no block and whose manifest records `not_available`.
- **When:** it runs before the report is assembled and again before the commit.

**D11. Report.** `reports/motco-span-split-<date>.md` is assembled by a scratch script from generated files. It contains:
- **Provenance:** the machine, the pin and wheel check, the gates, the kept-output comparison and the move.
- **MOTCO's validation** with its numbers:
  - C1–C3, and the real-noise size on this exact export (with its hashes);
  - the localization and power tables relevant to the reading;
  - the plasmode null share 0.236.
- **The predeclared policy and reading,** quoted from the proposal.
- **The split results:** the four primaries in full (rank, dimensions, F and p per part with roles, complement status, span concordance, share), then per family one table with the complement status against the primary's.
- **One section per re-issued step,** with the D9 outcome table.
- **The reading in the outcome that occurred,** the findings and the limitations.
- **The frozen set:** the eleven SVGs, with captions pasted verbatim from the frozen `figures.md`.

A check script verifies that:
- every split p-value, F, dimension and share, and every statistic, p-value, rank, CV AUROC and κ of every `result.json` is in the text;
- the captions match `figures.md`, and the frozen SVGs and inputs match the frozen manifest.

**D12. Machine.** The re-issue runs on `nixos-home`, which has the GENCODE annotation that `preprocessing-sensitivity` needs and the current outputs (sha256-equal to the cluster's). The gates compare runs made on the same machine. Every cross-machine comparison is recorded, not gated.

## Risks / Trade-offs

- **[A version-only move re-runs everything]** → Unavoidable: `_check_reference` compares MOTCO versions, and `figures` requires one version throughout. The compute is that of 2026-10-31 (a few hours). The split test adds 999 O(n²) draws per analysis (seconds).
- **[Forty-four analyses × three parts invite selective reading]** → Only the two `cohort_30` primary complements are results. Every other status is dependence. `total`, the share and the SS p-values are never read, and no table carries a `total` status.
- **[A complement rejection is over-read as a different program]** → The reading names the two other sources of out-of-span energy and forbids direction, feature and program claims. The split test has no geometry.
- **[A complement non-rejection is read as agreement]** → It is uninformative by rule. No power was shown at SEA-AD's stage signal, and none for sparse differences at any signal.
- **[Span concordance is trivially "agree" while nothing rejects]** → That is the expected cell and is stated as such. It carries no evidence either way.
- **[MOTCO's intermediate releases moved the default path]** → The D6.1 gate stops the change before anything is moved.
- **[Cross-machine superseded comparison]** → Expected at rounding level and recorded. The gates are same-machine.
- **[Earlier reports cite moved hashes]** → D7 checks that every cited hash and every frozen-set input is found after the move.
- **[`draws_match_rrpp` false]** → It can only happen at `n_jobs` ≠ 1. The record check requires it at `n_jobs = 1`, which every reported run uses.
- **[No split test on `latent_adjusted`]** → This is MOTCO's design (defined on the covariate-blind basis only). It is recorded per arm, and the arm's other readings are unchanged.

## Migration Plan

1. Branch `mtg-span-split-test`: pin, code, docs (D1–D5), verified on a scratch fixture.
2. Run the gates on `nixos-home` (D6). Stop on failure.
3. Move the superseded outputs (D7), then rerun everything (D8).
4. Run the comparisons and the consistency check (D9, D10), then write the report and its check (D11).
5. On the cluster: pull, `uv sync --locked`, the same move, then rsync the re-issued directories and verify their sha256.
6. PR; after the merge, archive and sync the specs.

**Rollback.** Before step 3 nothing under `results/` has changed. After it, `MOVED.json` lists every move, and moving the directories back with their hashes checked restores the 0.13.0 state on the main commit.

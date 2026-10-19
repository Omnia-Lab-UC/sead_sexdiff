# Proposal

## Why

The covariate-adjusted analysis (`reports/motco-astrocyte-2026-10-19.md`, archived change `2026-10-19-mtg-astrocyte-motco-analysis`) found that MOTCO's double CV selects a rank-1 latent space in both cohorts. At rank 1 `angle` is 0° or 180° by construction and carries no orientation information, and rank 1 lies outside the range MOTCO's operating studies cover (rank 2 or more in essentially every replicate). The report names the forced-rank sweep as the first follow-up: MOTCO's documented dimensionality sensitivity analysis (`docs/api/analysis.md`, "re-run with `--forced-components 2`, `3`, … into separate output directories"), to be read as sensitivity to the CV-selected analysis, never as its replacement. Nothing in this repo runs it yet.

## What Changes

- **New `rank-sensitivity` step.** For each cohort (`cohort_30` primary, `cohort_50` sensitivity), rerun the covariate-adjusted trajectory analysis with the latent rank fixed at each value of a predeclared ladder, `ANALYSIS_FORCED_RANKS = (2, 3, 4, 6, 9, 12)`, through the same MOTCO entry point (`analyze_trajectories` with `forced_components = k`). Everything else is identical to the `analyze` step: the same export inputs, blocks, stage order, six covariates, net-displacement orientation, 999 RRPP permutations and the same seed `ANALYSIS_SEED` for every rank, so that every rank draws the same permutation indices and differences between ranks are measurement differences, not draw differences. The ladder: 2 is `n_stages − 1`, the rank MOTCO's stage-supervised rule is expected to select with three stages and the lowest rank at which `angle` can carry orientation information; 3 is the next rank; 4, 6, 9 and 12 are the steps of MOTCO's own latent-rank ladder (`docs/reports/latent-rank-ladder-2026-09-08.md`), and 4 = 2(K − 1) is the forced rank with a Type I check at this design (`pls_forced4 × nd`, `docs/reports/net-displacement-type-i-2026-10-13.md` §2). Rank 1, the CV choice, cannot be forced (MOTCO accepts `forced_components ≥ 2` only); the CV-selected analysis is the reference row instead.
- **The CV-selected analysis is the reference, read, not rerun.** The sweep requires `results/motco_analysis/<cohort>/` from the `analyze` step, checks that its input hashes, MOTCO version, seed and permutation count equal the sweep's, and places its record as the `cv` row of the sweep table. The primary analysis is never modified.
- **Outputs per cohort** (`results/motco_rank_sensitivity/<cohort>/`): one MOTCO artifact per forced rank in `rank_<kk>/` (the same seven files the `analyze` step writes, with `component_selection = forced`), a cross-rank table `rank_sensitivity.csv` (reference row plus one row per rank: `delta`, `angle`, `shape`, their p-values and null summaries, runtime), its Markdown rendering `rank_sensitivity.md`, and `rank_sensitivity_manifest.json` (ladder, seed, permutations, reference analysis hashes, per-rank directories, software, git revision, findings).
- **A committed report** `reports/motco-astrocyte-rank-sensitivity-<date>.md` with the sweep table per cohort pasted from the generated files, checked against every `result.json` before commit, and the binding reading: the CV-selected analysis stays primary; rank-dependence is reported, never resolved by choosing a rank after seeing the results; `angle` and `shape` non-rejections stay uninformative at every rank; MOTCO's forced-rank operating characteristics (Type I holds at forced 4 at the Phase 6 design, unadjusted; forced 2(K − 1) loses `shape` power and leaks `angle` into the magnitude cell; ranks above `n_stages − 1` widen the `angle` null at the Phase 5 design) are carried as limitations, and the covariate adjustment has no Type I check at any forced rank.
- The `summary.md` renderer handles forced-rank records (forced rank beside `angle`, no CV block) without changing its output for CV records. README and AGENTS.md document the command, the new constants and the reading rules.

Out of scope, as later changes: preprocessing sensitivity exports (other N, keeping chrX, residualized features); pathway/module enrichment; any MOTCO-side change (whether CV may select rank 1; per-worker RRPP seeding).

## Capabilities

### New Capabilities
- `motco-rank-sensitivity`: the forced-rank dimensionality sensitivity sweep of the covariate-adjusted astrocyte trajectory analysis, with a predeclared rank ladder, the CV-selected analysis as its reference, reproducible per-rank outputs, a generated cross-rank table and a committed report that reads the sweep as sensitivity only.

### Modified Capabilities
<!-- None. `motco-trajectory-analysis` is unchanged: the CV-selected analysis, its outputs and its report keep every requirement; the sweep reads them. -->

## Impact

- `src/sead_sexdiff/data.py`: new constants `ANALYSIS_FORCED_RANKS` and `RANK_SENSITIVITY_*` (output names), a `run_rank_sensitivity` function, and a refactor that moves the per-cohort body of `run_analysis` (dataset construction, MOTCO call, artifact, manifest, summary) into a shared helper so the forced runs write the same artifact as the CV run. `run_analysis` output must stay byte-identical. `render_analysis_summary` gains a forced-rank branch.
- `src/sead_sexdiff/cli.py`: new `rank-sensitivity` subcommand (`--cohort`, `--input-dir`, `--analysis-dir`, `--out-dir`, `--n-jobs`, `--force`). The ladder is a constant, not a flag.
- `reports/`: one new committed report. `results/` stays gitignored.
- No dependency change: MOTCO 0.8.0 already provides `forced_components`.
- Compute: forced runs skip cross-validation, so each is one PLS fit plus 999 RRPP refits in a ≤ 12-dimensional space; seconds per rank, a few minutes for both cohorts, on any workstation holding `results/motco_input/` and `results/motco_analysis/`. The machine this change is proposed on holds neither; the real run (tasks group 3) needs the workstation that ran `analyze`, or those two directories synced to it.
- Scientific choices introduced here (flag and ask before changing): the rank ladder `(2, 3, 4, 6, 9, 12)`, the same seed for every rank, `n_jobs = 1` for the reported results, and the rule that the CV-selected analysis stays primary.

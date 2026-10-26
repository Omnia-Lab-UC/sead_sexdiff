# Design

## Context

- **MOTCO 0.10.0** (`v0.10.0`, commit `098b170`) adds `TrajectoryAnalysisParams.min_components` (default 2). The double CV then selects from the candidates `min_components … max_components − 1`, and every outer-fold selection, the modal rank and the pooled fit are at or above it. `1` normalizes to `None`, which is left out of `to_dict()` and of the record's `latent` block, so a floor-1 analysis reproduces every 0.9.0 artifact apart from `software`. A floor above 1 is recorded in `params` and in `latent.integration_params`. `forced_components` skips CV, so the floor does not act on a forced run (its `params` still carries `min_components = 2`). The CLI gains `motco analyze --min-components` (default 2).
- **`_analyze_cohort`** is the one place an analysis is run and recorded (`analyze`, `rank-sensitivity`, `preprocessing-sensitivity`, `latent-sensitivity`, `exclusion-sensitivity`). The five callers build `TrajectoryAnalysisParams(permutations=ANALYSIS_PERMUTATIONS, seed=ANALYSIS_SEED, n_jobs=n_jobs[, forced_components=k | covariate_latent_adjustment=…])`. `_cli_equivalent` and `render_analysis_summary` derive from `params` and the record.
- **`_check_reference(cohort_dir, analysis_dir, inputs, motco_version)`** is the reference check shared by the rank, preprocessing, latent and exclusion sensitivity steps and `_stage_de_prepare`, which `stage-de-limma` reuses. It already compares the export hashes, the MOTCO version, the seed, the permutations and `component_selection = cv`. A pin move therefore already makes every downstream step refuse the superseded primary.
- **What is on disk** (`nixos-home`): astrocyte `motco_analysis/` (0.9.0, made on the cluster), `motco_analysis_0.8.0/` (the record kept by the latent-adjustment change), `motco_rank_sensitivity/`, `motco_diagnostics/` and six of the eight `motco_sensitivity/` variants (all made against the 0.8.0 analysis), `motco_latent_sensitivity/` (0.9.0, cluster), `stage_de/` (cluster), `stage_de_limma/`, `motco_exclusion_sensitivity/`; microglia `motco_analysis/`, `motco_diagnostics/`, `motco_latent_sensitivity/`, `stage_de/`, `stage_de_limma/`, `motco_exclusion_sensitivity/`, `positive_control/` (all made here). The cluster is the canonical copy. The two trend variants and the variant diagnostics of the feature-selection and trend-normalization reports were made elsewhere or in scratch directories.
- **PLS nesting**: the first latent component is identical at every rank. The 0.8.0 diagnostic recorded `pls_0` as `np.allclose`-identical across the primary and all six forced ranks.

## Goals / Non-Goals

**Goals:**
- The floor is set in one constant, passed explicitly, recorded everywhere it acts, and checked by every reader.
- A same-machine gate proves that, at floor 1, the pin move changes nothing before anything floored is run.
- After the re-issue, every output outside `results/superseded/` rests on MOTCO 0.10.0 at floor 2, and every superseded byte is still resolvable by its hash.

**Non-Goals:**
- A CLI flag for the floor, or a floor-1 "unfloored arm": the superseded outputs are that record.
- Rerunning MOTCO-independent steps, or moving their outputs.
- Bringing the laptop or `nixos-office` copies up to date (they re-sync from the cluster when next used).

## Decisions

**D1. One constant, passed explicitly, tied to MOTCO's default.** `ANALYSIS_RANK_FLOOR = 2` next to the other `ANALYSIS_*` constants. The five callers build `TrajectoryAnalysisParams(..., min_components=ANALYSIS_RANK_FLOOR)` (the forced-rank sweep too, so its `params` stays uniform). `_analyze_cohort` checks `TrajectoryAnalysisParams().min_components == ANALYSIS_RANK_FLOOR` before reading inputs; the error names both values. It records `measurement.min_components` (`params.resolved_min_components`) beside `component_selection`, and `_cli_equivalent` adds `--min-components <floor>` for CV runs only (MOTCO ignores it under `--forced-components`). `render_analysis_summary` prints "CV-selected rank k (rank floor f)". `ANALYSIS_MOTCO_MIN_VERSION = "0.10.0"`; `LATENT_MOTCO_MIN_VERSION` stays `0.9.0`, now implied by the analysis minimum.
- *Alternative: rely on MOTCO's default without passing it.* Rejected: AGENTS.md keeps MOTCO's defaults, and passing the constant does keep them, but a later MOTCO that moved its default would silently change the primary. Passing it and checking equality makes such a move an error instead.
- *Alternative: a `--rank-floor` flag.* Rejected: measurement settings here are constants, never flags.

**D2. The reference check reads the floor.** `_check_reference` reads `record["params"].get("min_components", 1)` (absent means floor 1, MOTCO's omitted spelling) and stops naming the file, both floors and `analyze` when it differs from `ANALYSIS_RANK_FLOOR`. The MOTCO-version check already catches a 0.9.0 reference; the floor check also catches a 0.10.0 reference run at floor 1 by hand (`motco analyze --min-components 1`). The rank-1 code paths stay (the summary note, the `_analyze_cohort` finding, `RANK_1_NOTE`) because a floor-1 analysis can still be read. The sweep's finding text loses its rank-1 branch only where it assumes a rank-1 reference (spec delta `motco-rank-sensitivity`).

**D3. The gate (spec `motco-release-migration`).** On `nixos-home`, before any floored run:
1. `git worktree add <scratch>/pre b1a67a7` (the pre-change `main`), `uv sync --locked` there (MOTCO 0.9.0), and `analyze --cell-type {astrocyte,microglia} --cohort all --input-dir <real export> --out-dir <scratch>/v090/<cell type>`, all through `nix develop -c`.
2. From the change branch (MOTCO 0.10.0), a scratch script runs `_analyze_cohort` with `min_components=1` on the same exports into `<scratch>/v0100_floor1/<cell type>`.
3. A comparison script: the four CSVs byte for byte; `result.json` and `analysis_manifest.json` after dropping `runtime_seconds`, `created_at`, `git`, `software.motco`, `cli_equivalent` and the new `measurement.min_components`. Pass, or the change stops.
4. The same script compares the floor-1 reruns with the kept primaries: microglia (made on this machine) is expected to be byte-identical, astrocyte (made on the cluster) at rounding level. Recorded, not a gate.
- *Why floor 1 through `_analyze_cohort` and not the CLI:* the floor is not a flag (D1), and the helper is the code path whose refactor is under test.

**D4. Superseded layout.** `results/superseded/mtg-rank-floor-primary/` mirrors `results/`: astrocyte `motco_analysis`, `motco_rank_sensitivity`, `motco_diagnostics`, `motco_sensitivity`, `motco_latent_sensitivity`, `stage_de`, `stage_de_limma`, `motco_exclusion_sensitivity`; microglia under `microglia/` (`motco_analysis`, `motco_diagnostics`, `motco_latent_sensitivity`, `stage_de`, `stage_de_limma`, `motco_exclusion_sensitivity`, `positive_control`). The earlier record `results/motco_analysis_0.8.0/` moves under the same convention to `results/superseded/mtg-astrocyte-latent-adjustment-sensitivity/motco_analysis/`, so "outside `results/superseded/`" means current. A move script hashes every file, `mv`s on the same filesystem, re-hashes, and writes `results/superseded/<change>/MOVED.json` (old path, new path, sha256). On the cluster the same script runs over whatever MOTCO-dependent directories exist there (it may hold outputs this machine lacks; they move too).
- *Alternative: delete and rely on the cluster's copy.* Rejected: the cluster gets the same move, and nothing that cost compute is deleted.
- *Alternative: version-named directories (`motco_analysis_0.9.0/`).* Rejected: the steps' outputs were made under 0.8.0 and 0.9.0, and per-directory suffixes scatter the records. Naming the change that superseded them stays unambiguous.

**D5. Rerun order and commands** (all `nix develop -c uv run sead-sexdiff …`, `--n-jobs 1` where taken, `--n-cpus 10` for `stage-de`):

| # | astrocyte | microglia |
|---|---|---|
| 1 | `analyze --cohort all` | `analyze --cell-type microglia` |
| 2 | `rank-sensitivity --cohort all` | – |
| 3 | `axis-diagnostics` (with the sweep) | `axis-diagnostics --cell-type microglia` |
| 4 | `preprocessing-sensitivity --variant all`, then `axis-diagnostics` per variant root into `results/motco_sensitivity/<variant>/diagnostics/` | – |
| 5 | `latent-sensitivity --cohort all`; `axis-diagnostics` on the primary and both arms into `results/motco_latent_sensitivity/diagnostics/<name>/` | the same with `--cell-type microglia` (arms only, as before) |
| 6 | `stage-de --omic both`; `stage-de-limma --omic both` | the same with `--cell-type microglia` |
| 7 | `exclusion-sensitivity --cohort all`; on its variant root `axis-diagnostics`, `stage-de --omic both`, `stage-de-limma --omic rna` (the 2026-10-23 report's section 8) | the same with `--cell-type microglia` |
| 8 | – | `positive-control` and `positive-control --basis latent_adjusted` |

Steps 2–7 depend only on step 1 of their cell type; step 8 needs steps 3, 5 and 6 of both. Astrocyte and microglia chains can run side by side (as in the exclusion change). Each step writes fresh directories (the superseded ones are gone), so no `--force` is needed and no resume reads an old fit record. Reproducibility reruns (`--force` into scratch) for `analyze` on one cohort per cell type, as each earlier report did.

**D6. Expected reproductions, checked by the comparison script.** `pls_0` scores of every re-issued analysis against the superseded one (`np.allclose`, max |diff| recorded); E3's correlations; the exclusion flag check's `pls_0` AUCs and the axis check; every stage-DE and limma table (same machine for microglia and astrocyte limma: byte-identical expected; astrocyte DESeq2 made on the cluster: rounding level expected); E1's counts. A difference becomes a finding in the report, never a gate, and changes no reading (spec `motco-release-migration`).

**D7. Version-consistency check.** A script walks `results/` minus `results/superseded/`, parses every `*manifest*.json` with a MOTCO version (`software.motco`, or nested `software` of comparison manifests), and lists any that are not `0.10.0`; it checks every `result.json` for `params.min_components == 2`, or a forced rank. Run before the report is assembled and again before commit.

**D8. Report.** `reports/motco-rank-floor-primary-<date>.md`, assembled by a scratch script as before: provenance (pin, gate, kept-output comparison, move record), the floored primary per cell type and cohort (`summary.md` content plus the eigengap from `config_spectrum` and the null summaries), then one section per re-issued step pasting its generated Markdown (`rank_sensitivity.md`, `sensitivity.md` + `adoption_check.csv` + precedence, `latent_sensitivity.md` + `mediator_check.csv`, the stage-DE summaries and concordance cells, `positive_control.md` of both bases, `exclusion_sensitivity.md` + `axis_check.csv`) with a generated "superseded → floored" outcome table, then the predeclared reading in the outcome that occurred, and limitations. A check script finds every statistic, p-value, rank and CV AUROC of every re-issued `result.json`, every E1–E3 value and every concordance cell in the text. Large per-feature tables are referenced, not pasted, to keep the report readable.

## Risks / Trade-offs

- [The floor changes a primary rejection status and the change is read as discovery] → the reading is written in the proposal before any run; the superseded status stands beside the new one as rank-rule dependence; MOTCO's check found no power at its strengths, which the report quotes.
- [CV at floor 2 hits a sample too small for the floor] → MOTCO errors before fitting rather than lowering the rank; the smallest analysed cohort (63 donors, 10,000 features) is far from that bound. An error stops that step, is reported, and is not worked around.
- [The floor interacts with covariates or the latent adjustment in ways MOTCO did not measure] → carried as a limitation (MOTCO's report says so); no Type I claim is made for the floored covariate-adjusted or latent-adjusted analyses beyond what MOTCO measured.
- [A step's rerun silently reads a superseded file] → the old directories are moved first, so a missing input is an error; D7 lists any stale manifest.
- [Preprocessing variants at floor 2 expose more components to the adoption rule] → the rule is unchanged and reads every component, so a candidate may fail where it passed; that outcome stands, as predeclared.
- [Compute] → about four dozen CV analyses of 1–2 min plus sixteen variant exports and the stage DE; a few hours on `nixos-home`, run in the background with logs in the scratch directory.

## Migration Plan

1. Code and pin (D1, D2), docs; commit.
2. Gate (D3) on `nixos-home`; stop on failure.
3. Move (D4) on `nixos-home`; record hashes of the MOTCO-independent outputs.
4. Reruns (D5), comparisons (D6), version check (D7).
5. Report (D8), check, commit.
6. On the cluster: the same move, then copy the re-issued directories and the report's inputs, verified by sha256.

Rollback: revert the pin and the code commit, and move the directories back from `results/superseded/mtg-rank-floor-primary/` with `MOVED.json` on both machines.

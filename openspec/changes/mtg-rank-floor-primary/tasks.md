# Tasks

## 1. Pin and code

- [x] 1.1 Move the pin: `uv add "motco @ git+https://github.com/tomszar/MOTCO@v0.10.0"`. Verify that `uv.lock` records commit `098b170e618959591c7a3d2277487765dcd0459a`, that scikit-learn's version is MOTCO's exact pin (not loosened), and that `nix develop -c uv run python -c "import motco; print(motco.__version__)"` prints `0.10.0`.
- [x] 1.2 Add `ANALYSIS_RANK_FLOOR = 2`, set `ANALYSIS_MOTCO_MIN_VERSION = "0.10.0"`, pass `min_components=ANALYSIS_RANK_FLOOR` in the five `TrajectoryAnalysisParams` constructions, check MOTCO's default against the constant in `_analyze_cohort`, record `measurement.min_components`, add `--min-components` to `_cli_equivalent` for CV runs and the floor to `render_analysis_summary` (design D1). Verify on a reduced fixture (every k-th feature of the real astrocyte export, in the scratchpad) that `analyze` selects a rank ≥ 2 and records `min_components = 2` in `result.json` `params`, in the manifest and in `cli_equivalent`, that running the recorded `cli_equivalent` with `motco analyze` gives the same statistics and rank, and that a patched MOTCO default of 3 stops the run naming 3 and 2, with nothing written.
- [x] 1.3 Add the floor check to `_check_reference` and adjust the rank sweep's finding for a CV-selected rank outside the ladder (design D2). Verify on the fixture that `rank-sensitivity` and `latent-sensitivity` against a floor-1 reference made under 0.10.0 (`_analyze_cohort` with `min_components=1`) stop naming the reference `result.json`, floor 1 and floor 2, with nothing written, and that against the floored fixture analysis the sweep's manifest records the floor and the forced run at the CV-selected rank equals the reference (statistics and p-values).
- [x] 1.4 Update README (the pin, `ANALYSIS_RANK_FLOOR`, the summary line, the reference check) and AGENTS.md (the MOTCO pin to `v0.10.0`, `ANALYSIS_MOTCO_MIN_VERSION`, `ANALYSIS_RANK_FLOOR` among the analysis constants under "flag and ask", the rank-sensitivity notes that assumed a rank-1 reference). Verify that every constant named in both files exists in `data.py` and that no text still describes the primary as rank 1 by rule.

## 2. Same-machine gate (nixos-home, Nix shell)

- [x] 2.1 Create a worktree of the pre-change commit `b1a67a7` in the scratchpad, `uv sync --locked` there (MOTCO 0.9.0), and run `analyze --cohort all` for both cell types on the real exports into scratch directories (design D3 step 1). Verify that the four `analysis_manifest.json` record MOTCO 0.9.0 and the real exports' hashes.
- [x] 2.2 From the change branch, run `_analyze_cohort` with `min_components=1` on the same four exports into scratch directories (design D3 step 2). Verify that the four records carry no `min_components` in `params` and MOTCO 0.10.0 under `software`.
- [x] 2.3 Compare 2.1 and 2.2 with a scratch script (design D3 step 3). Verify that the four CSVs per cohort are byte-identical and the JSON equal after dropping `runtime_seconds`, `created_at`, `git`, `software.motco`, `cli_equivalent` and `measurement.min_components`. If anything differs, stop the change and report the differing file.
- [x] 2.4 Compare 2.2 with the kept primaries `results/motco_analysis/` and `results/microglia/motco_analysis/` (design D3 step 4). Verify that the outcome is recorded for the report: byte-identical expected for microglia (made here); for astrocytes (made on the cluster) the largest absolute differences and whether p-values, ranks and `analysis_metadata.csv` are equal.

## 3. Move the superseded outputs (nixos-home)

- [x] 3.1 Record the sha256 of every file under both cell types' `pseudobulk/` and `motco_input/` and of `data/annotation/` into the scratchpad. Verify that the listing is complete (file counts per directory).
- [x] 3.2 Write and run the move script (design D4): the MOTCO-dependent directories of both cell types to `results/superseded/mtg-rank-floor-primary/` (relative paths kept) and `results/motco_analysis_0.8.0/` to `results/superseded/mtg-astrocyte-latent-adjustment-sensitivity/motco_analysis/`, each with `MOVED.json`. Verify that every moved file hashes as before, that no MOTCO-dependent directory is left at its old path, and that the sha256 of each `result.json` cited by the earlier reports is found under `results/superseded/`.
- [x] 3.3 Document `results/superseded/<change>/` in README and AGENTS.md (layout, `MOVED.json`, never an input of any step). Verify that the documented paths match the moved tree.

## 4. Re-issue every MOTCO-dependent output (nixos-home, Nix shell, `n_jobs = 1`)

- [x] 4.1 `analyze --cohort all` and `analyze --cell-type microglia --cohort all`. Verify that the four records select a rank ≥ 2 with `min_components = 2`, MOTCO 0.10.0, seed 20261019 and 999 permutations, and that a `--force` rerun of `cohort_30` per cell type into the scratchpad gives byte-identical CSVs and equal statistics, p-values, null summaries and rank.
- [x] 4.2 `rank-sensitivity --cohort all`. Verify six forced ranks per cohort, the floor in the manifest, and that the forced run at the CV-selected rank (when it is in the ladder) equals the reference, or that the manifest records why not.
- [x] 4.3 `axis-diagnostics` (astrocytes, with the sweep) and `axis-diagnostics --cell-type microglia`. Verify that the input-hash checks pass and that the manifests record whether `pls_0` is identical across the primary and the forced ranks.
- [x] 4.4 `preprocessing-sensitivity --variant all`, then `axis-diagnostics --input-dir <variant>/motco_input --analysis-dir <variant>/motco_analysis --out-dir <variant>/diagnostics` per variant root (design D5 row 4). Verify that all eight variants are complete in both cohorts, that `adoption_check.csv` exists for the three candidates and `adoption_precedence` is set (not `incomplete`) in both trend candidates' comparison manifests, and that every variant analysis records the floor.
- [x] 4.5 `latent-sensitivity --cohort all` for both cell types, then `axis-diagnostics` on the astrocyte primary and on both arms of both cell types into `<root>/motco_latent_sensitivity/diagnostics/<name>/`. Verify that both arms record the floor, that `mediator_check.csv` was written before the arms (file times), and that each `extended_covariates` arm's rank and CV AUROC equal its primary's, or that the difference is recorded as a finding.
- [x] 4.6 `stage-de --omic both --n-cpus 10` and `stage-de-limma --omic both` for both cell types. Verify that every concordance cell reads the floored primary's `delta` status and that the fit manifests' reference hashes are the new primaries'.
- [x] 4.7 `exclusion-sensitivity --cohort all` for both cell types, then on each variant root `axis-diagnostics`, `stage-de --omic both --n-cpus 10` and `stage-de-limma --omic rna` (the 2026-10-23 report's section 8). Verify that `flag_check.csv` predates the variant export (file times), that the variant analyses record the floor, and that the stage-DE manifests record the `donor_exclusion` block.
- [x] 4.8 `positive-control --cohort all` and `positive-control --basis latent_adjusted --cohort all`. Verify that both bases' outputs exist for both cohorts and that every manifest read is the re-issued one (the chain checks pass).

## 5. Comparisons and consistency checks

- [x] 5.1 Write and run the comparison script of design D6 (superseded against re-issued): `pls_0` of every analysis, E1–E3 values, the flag check's `pls_0` AUCs and the axis checks, every stage-DE and limma table, plus per step the outcome table the report pastes (rank, CV AUROC, statistics, p-values, rejections, adoption check and precedence, concordance cells, control outcomes). Verify that it runs on the full trees and lists every expected reproduction that failed (each becomes a report finding).
- [x] 5.2 Run the version-consistency script of design D7. Verify that it lists no manifest outside `results/superseded/` with a MOTCO version other than 0.10.0 and no unforced `result.json` without `min_components = 2`.
- [x] 5.3 Re-hash the files recorded in 3.1. Verify that every hash is unchanged.

## 6. Report

- [x] 6.1 Write `reports/motco-rank-floor-primary-<date>.md` from the generated files (design D8): provenance with the gate (2.3), the kept-output comparison (2.4) and the move (3.2); the floored primaries of both cell types with the eigengap and null summaries beside `angle`; one section per re-issued step with its generated tables and the superseded → floored outcome table; the predeclared reading (proposal) in the outcome that occurred; the findings of 5.1; limitations (the floor's interaction with covariates and the latent adjustment not measured by MOTCO, power at the null in MOTCO's check, no operating characteristic measured on these blocks). Verify with a check script that every statistic, p-value, rank and CV AUROC of every re-issued `result.json`, every E1–E3 value, every adoption-check cell, every axis-check value and every concordance cell is in the text.
- [x] 6.2 Run `openspec validate mtg-rank-floor-primary`. Verify that it passes (requirement-length warnings only, as elsewhere in the repo).

## 7. Cluster

- [x] 7.1 On the cluster, pull the branch, `uv sync --locked`, and run the move script of 3.2 over the MOTCO-dependent directories present there. Verify that every moved file hashes as before and that the cluster's `pseudobulk/` and `motco_input/` directories are untouched.
- [x] 7.2 Copy every re-issued directory of both cell types from `nixos-home` to the cluster with rsync. Verify that every copied file's sha256 equals the workstation's.

## Workflow follow-up

- Commit the code, the report and the task records; open the PR against `main`.
- Archive the change after merge and sync its specs (`motco-release-migration` new; `motco-trajectory-analysis`, `motco-rank-sensitivity`, `cell-type-positive-control` modified).
- Update the case-study memory with the floored outcomes.

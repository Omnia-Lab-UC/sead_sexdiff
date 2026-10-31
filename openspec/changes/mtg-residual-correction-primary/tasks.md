# Tasks

## 1. Pin and code

- [x] 1.1 Move the pin to PyPI (design D1).
  - Edit `pyproject.toml`: list `motco==0.13.0` and remove the `[tool.uv.sources]` motco entry. Then run `uv lock`.
  - Verify that `uv.lock` records motco 0.13.0 from the PyPI registry with wheel sha256 `614c2e9cc7fccb55…`.
  - Verify that `uv export --no-hashes` differs from its pre-change output only in the motco line, and that scikit-learn is still MOTCO's exact pin.
  - Verify with a scratch script that the downloaded wheel's sha256 equals the lock's, and that its `motco/` tree equals `git -C ../MOTCO archive v0.13.0 src/motco` file by file.
  - Verify that `nix develop -c uv run python -c "import motco; print(motco.__version__)"` prints `0.13.0`.
  - Record the wheel hash, the tag commit and both comparisons for the report.
- [x] 1.2 Add the correction (design D2): `ANALYSIS_RESIDUAL_CORRECTION = "df"`, `ANALYSIS_MOTCO_MIN_VERSION = "0.13.0"`, `residual_correction=` in the five `TrajectoryAnalysisParams` constructions, the record check in `_analyze_cohort` before anything is written, `measurement.residual_correction` in the manifest, `--residual-correction` in `_cli_equivalent`, and the correction line in `render_analysis_summary`. Verify on a reduced fixture (every k-th feature of the real astrocyte export, in the scratchpad):
  - `analyze` records `params.residual_correction = "df"` and a `measurement.residual_correction` with reduced rank 10 and κ = √(80 / 70) to 1e-12, in `result.json` and in the manifest;
  - `summary.md` states the correction;
  - running the recorded `cli_equivalent` with `motco analyze` gives the same statistics, p-values and rank;
  - a patched `analyze_trajectories` whose record lacks the correction stops the run naming both settings, with nothing written;
  - `render_analysis_summary` on a superseded 0.10.0 `result.json` gives the same text as before the change.
- [x] 1.3 Add the correction check to `_check_reference` and the correction to the rank-sweep manifest (design D3). Verify on the fixture that:
  - `rank-sensitivity` and `stage-de` against an uncorrected 0.13.0 reference (`_analyze_cohort` with `residual_correction=None`) stop naming the reference `result.json`, the recorded setting, `df` and `analyze`, with nothing written;
  - against the corrected fixture reference, the sweep's manifest records `df`, every forced record carries `df`, and the forced run at the CV-selected rank (when it is in the ladder) equals the reference in statistics and p-values.
- [x] 1.4 Update README and AGENTS.md:
  - the pin rule (exact PyPI version, wheel checked against the `Omnia-Lab-UC/MOTCO` release tag, never a git or path dependency);
  - `ANALYSIS_MOTCO_MIN_VERSION` 0.13.0, including the `analyze` step's "≥ 0.10.0" text;
  - `ANALYSIS_RESIDUAL_CORRECTION` among the analysis constants under "flag and ask", as the one MOTCO default overridden, by MOTCO's verdict;
  - the summary line and the reference check;
  - the latent-sensitivity reading rule that the `latent_adjusted` arm's `angle` and `shape` p-values are reported but not read (`residualized_basis`).

  Verify that every constant named in both files exists in `data.py`, that no text still instructs a git pin or names `tomszar/MOTCO` as the source, and that no text describes the primary's null as uncorrected.

## 2. Gates (re-issue machine, Nix shell)

- [x] 2.1 Prepare the re-issue machine (design D10). On `nixos-home`, verify that each current primary `result.json` and `analysis_manifest.json` of both cell types hashes as on the cluster. If `nixos-home` is unavailable, rsync the cluster's current directories of both cell types (not `superseded/`) to `nixos-dell-laptop` and verify every file's sha256 against the cluster's. Record the machine for the report.
- [x] 2.2 Pin-move gate, previous pin (design D4.1): create a scratchpad worktree of the pre-change commit, `uv sync --locked` there (MOTCO 0.10.0), and run `analyze --cell-type {astrocyte,microglia} --cohort all` on the real exports into scratch directories. Verify that the four manifests record MOTCO 0.10.0 and the real exports' hashes.
- [x] 2.3 Pin-move gate, new pin without the correction: from the change branch, run `_analyze_cohort` with the analyze parameters and `residual_correction=None` on the same four exports into scratch directories. Verify that the records carry no `residual_correction` and MOTCO 0.13.0 under `software`.
- [x] 2.4 Compare 2.2 with 2.3 using a scratch script. Pass condition: the four CSVs per cohort are byte-identical, and `result.json` and the manifest are equal after dropping `runtime_seconds`, `created_at`, `git`, `software.motco`, `cli_equivalent` and `measurement.residual_correction`. If anything differs, stop the change and report the differing file.
- [x] 2.5 Compare 2.3 with the kept primaries `results/motco_analysis/` and `results/microglia/motco_analysis/`. Record the outcome for the report: byte-identical on the machine that made them; otherwise the largest absolute differences and whether p-values, ranks and `analysis_metadata.csv` are equal.
- [x] 2.6 Contract gate (design D4.2): from the change branch, run `analyze --cell-type {astrocyte,microglia} --cohort all` into scratch directories and compare with 2.3. Pass condition:
  - the four CSVs are byte-identical;
  - `result.json` is equal after dropping `p_values`, `null_summary`, `config_spectrum.permutation_pooled_eigengap`, `params.residual_correction`, `measurement.residual_correction` and `measurement.runtime_seconds`;
  - each κ equals √(n / (n − 10)) to 1e-12 (1.0690, 1.0753, 1.0731, 1.0903).

  If anything differs, stop the change, name the cohort and the value, and report it to MOTCO.

## 3. Move the superseded outputs

- [x] 3.1 Record the sha256 of every file under both cell types' `pseudobulk/` and `motco_input/` and under `data/annotation/` into the scratchpad. Verify that the listing is complete (file counts per directory).
- [x] 3.2 Run `scripts/move_superseded.py --change mtg-residual-correction-primary` over the MOTCO-dependent directories of both cell types (design D5). Verify that:
  - every moved file hashes as before and `MOVED.json` is written;
  - no MOTCO-dependent directory is left at its old path;
  - the sha256 of each `result.json` cited by `reports/motco-rank-floor-primary-2026-10-26.md` is found under `results/superseded/mtg-residual-correction-primary/`.
- [x] 3.3 Add `mtg-residual-correction-primary/` (the floored, uncorrected outputs of MOTCO 0.10.0, both cell types) to the superseded entries in AGENTS.md and README. Verify that the documented paths match the moved tree.

## 4. Re-issue every MOTCO-dependent output (Nix shell, `n_jobs = 1`)

- [x] 4.1 Run `analyze --cohort all` and `analyze --cell-type microglia --cohort all`. Verify that:
  - the four records carry `residual_correction = df` with the predicted κ, `min_components = 2`, MOTCO 0.13.0, seed 20261019 and 999 permutations;
  - their CSVs are byte-identical to the 2.6 scratch runs, with equal statistics, p-values, null summaries and rank (the reproducibility rerun).
- [x] 4.2 Run `rank-sensitivity --cohort all`. Verify that there are six forced ranks per cohort, each record carries `df`, and the manifest records the floor and `df`. Verify that the forced run at the CV-selected rank (when it is in the ladder) equals the reference, or that the manifest records why not.
- [x] 4.3 Run `axis-diagnostics` (astrocytes, with the sweep) and `axis-diagnostics --cell-type microglia`. Verify that the input-hash checks pass and that `pls_0` is recorded as identical across the primary and the forced ranks.
- [x] 4.4 Run `preprocessing-sensitivity --variant all`, then `axis-diagnostics` per variant root (design D6 step 4). Verify that all eight variants are complete in both cohorts, that `adoption_check.csv` exists for the three candidates, that `adoption_precedence` is set (not `incomplete`) in both trend candidates' comparison manifests, and that every variant record carries `df`.
- [x] 4.5 Run `latent-sensitivity --cohort all` for both cell types, then `axis-diagnostics` on the astrocyte primary and on both arms of both cell types into `<root>/motco_latent_sensitivity/diagnostics/<name>/`. Verify that:
  - both arms record `df` with reduced rank 12;
  - `mediator_check.csv` was written before the arms (file times);
  - each `extended_covariates` arm's rank and CV AUROC equal its primary's, or the difference is recorded as a finding.
- [x] 4.6 Run `stage-de --omic both --n-cpus 10` and `stage-de-limma --omic both` for both cell types. Verify that every concordance cell reads the corrected primary's `delta` status and that the fit manifests' reference hashes are the new primaries'.
- [x] 4.7 Run `exclusion-sensitivity --cohort all` for both cell types, then on each variant root `axis-diagnostics`, `stage-de --omic both --n-cpus 10` and `stage-de-limma --omic rna`. Verify that `flag_check.csv` predates the variant export (file times), that the variant records carry `df`, and that the stage-DE manifests record the `donor_exclusion` block.
- [x] 4.8 Run `positive-control --cohort all` and `positive-control --basis latent_adjusted --cohort all`. Verify that both bases' outputs exist for both cohorts and that every manifest read is the re-issued one (the chain checks pass).

## 5. Comparisons and consistency checks

- [x] 5.1 Write and run the comparison script of design D7 (superseded against corrected). Verify that:
  - it covers every analysis pair (4 primaries, 12 forced ranks, 16 variants, 8 arms, 4 exclusion variants) with the contract comparison, the κ check and the per-statistic p-values, directions and statuses;
  - it covers every p-free output listed in D7;
  - it produces the report's outcome tables;
  - it lists every failed reproduction and every status change (each becomes a report finding).
- [x] 5.2 Run the version-consistency check of design D8. Verify that it lists no manifest outside `results/superseded/` with a MOTCO version other than 0.13.0, and no `result.json` without `residual_correction = df` and either `min_components = 2` or a forced rank.
- [x] 5.3 Re-hash the files recorded in 3.1. Verify that every hash is unchanged.

## 6. Report

- [x] 6.1 Write `reports/motco-residual-correction-primary-<date>.md` from the generated files (design D9). It contains:
  - provenance: the machine, the pin and wheel check (1.1), the gates (2.4, 2.6), the kept-output comparison (2.5) and the move (3.2);
  - MOTCO's three reports with their numbers, and the statement that the 2026-10-26 primaries' p-values were withheld;
  - the corrected primaries of both cell types with κ, the eigengap and the null summaries beside `angle`;
  - one section per re-issued step, with its generated tables and the superseded → corrected outcome table;
  - the predeclared reading (proposal) in the outcome that occurred, including the `latent_adjusted` arm's `angle` and `shape` not read and the positive control's reproduction on both bases;
  - the findings of 5.1;
  - limitations: size measured at six covariates on the `cohort_30` tables only, and inferred for `cohort_50`, the eight-covariate arms, the exclusion variants, the forced ranks and the variants; `delta` conservative under `df`; MOTCO's `residualized_basis` diagnosis; every earlier limitation carried.

  Verify with a check script that every statistic, p-value, rank, CV AUROC and κ of every re-issued `result.json`, every E1–E3 value, every adoption-check cell, every axis-check value and every concordance cell is in the text.
- [x] 6.2 Run `openspec validate mtg-residual-correction-primary`. Verify that it passes (requirement-length warnings only, as elsewhere in the repository).

## 7. Cluster

- [x] 7.1 On the cluster:
  - pull the branch and run `uv sync --locked`, which now installs motco from PyPI on the login node;
  - verify that `import motco` prints 0.13.0;
  - run the move script of 3.2 over the MOTCO-dependent directories present there.

  Verify that every moved file hashes as before and that the cluster's `pseudobulk/` and `motco_input/` directories are untouched.
- [x] 7.2 Copy every re-issued directory of both cell types from the re-issue machine to the cluster with rsync. Verify that every copied file's sha256 equals the re-issue machine's.

## Workflow follow-up

- Commit the code, the report and the task records; open the PR against `main`.
- Archive the change after merge and sync its specs (modified: `motco-trajectory-analysis`, `motco-rank-sensitivity`, `motco-release-migration`, `motco-latent-adjustment-sensitivity`, `cell-type-positive-control`).
- Update the case-study memory with the corrected outcomes.

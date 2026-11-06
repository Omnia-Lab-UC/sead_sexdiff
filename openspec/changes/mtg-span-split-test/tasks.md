# Tasks

## 1. Pin and code

- [ ] 1.1 Record the predeclaration. Commit this change's planning artifacts (proposal, design, specs, tasks) on branch `mtg-span-split-test` and push. Record the commit hash for the report as the commit that fixed the policy and the reading. Verify with `git log` that it precedes every run of groups 2–4. No real-label split test may run before it.
- [ ] 1.2 Move the pin (design D1).
  - Run `uv add "motco==0.18.0"`.
  - Verify that `uv.lock` records motco 0.18.0 from the PyPI registry with wheel sha256 `b142a095c9fd6833…`.
  - Verify that `uv export --no-hashes` differs from main's output only in the motco line, and that scikit-learn is still `==1.8.0`.
  - Verify with a scratch script that the downloaded wheel's sha256 equals the lock's, and that its `motco/` tree equals `git -C ../MOTCO archive v0.18.0 src/motco` file by file.
  - Verify that `nix develop -c uv run python -c "import motco; print(motco.__version__)"` prints `0.18.0`.
  - Record the wheel hash, the tag commit (`bceead1`) and both comparisons for the report.
- [ ] 1.3 Build a synthetic fixture with a scratch script: random `rna` / `atac` blocks with the real cohorts' donor table, metadata and covariate shapes, written as an export directory in the scratchpad. Code checks then never show a real-label split value before the gates. Verify that `analyze --input-dir <fixture>` runs under 0.18.0 on today's code (no split block yet).
- [ ] 1.4 Add the constant and the record check (design D2):
  - `ANALYSIS_MOTCO_MIN_VERSION = "0.18.0"`, `ANALYSIS_SPAN_SPLIT_TEST` and the `SPAN_SPLIT_*` policy constants;
  - `span_split_test=` in the five `TrajectoryAnalysisParams` constructions, with the latent arm's `latent_adjustment is None` rule and its import-time assertion;
  - the `SPAN_SPLIT_PARTS` check against MOTCO's `INTERACTION_SPLIT_PARTS`;
  - the record check in `_analyze_cohort` before any write;
  - `measurement.span_split` in the manifest.

  Verify on the fixture that:
  - `analyze` records `params.span_split_test = true` and a block whose `span` dimension equals the selected rank, whose complement is the total minus it, and whose permutations, seed and `df` equal the manifest's, with `draws_match_rrpp = true`;
  - the manifest's `measurement.span_split` names the confirmatory part, the statistic and the parts not read;
  - a patched `analyze_trajectories` that drops the block stops the run naming the cohort and the block, with nothing written;
  - `latent-sensitivity` on the fixture writes a `latent_adjusted` record without a block and a manifest with `not_available` quoting MOTCO's reason, while `extended_covariates` carries the block;
  - a scratch copy of `LATENT_SENSITIVITY_ARMS` that breaks the rule trips the import-time assertion.
- [ ] 1.5 Add the summary table and the CLI flag (design D3: `render_analysis_summary`, `_cli_equivalent`). Verify that:
  - `summary.md` on the fixture shows the split table with the three roles and the share labelled as not read, and no interaction-SS p-value;
  - `render_analysis_summary` on each current 0.13.0 `result.json` of both cell types (no block) is byte-identical to its `summary.md` on disk;
  - running the fixture's recorded `cli_equivalent` with `motco analyze` gives the same statistics, p-values, rank and `span_split` block.
- [ ] 1.6 Add the split columns to `_sensitivity_row` (design D3). Verify on the fixture that:
  - `rank-sensitivity`, `preprocessing-sensitivity --variant keep_chrx`, `latent-sensitivity` and `exclusion-sensitivity` (with a fixture flag table) write their tables with the ten split columns, every cell equal to its row's `result.json`;
  - `span_concordant` follows its definition on hand-made records covering all four combinations;
  - the `latent_adjusted` row's split cells are empty, and no `total_rejected` column exists;
  - the stage-DE concordance cells are unchanged.
- [ ] 1.7 Add the split check to `_check_reference` (design D4). Verify on the fixture that:
  - `rank-sensitivity` and `stage-de` against a reference made with `span_split_test=False` stop naming the reference `result.json` and `analyze`, with nothing written;
  - against the fixture reference with the block, the sweep runs, every forced record carries a block whose `span` dimension is its forced rank, and the forced run at the CV-selected rank (when it is in the ladder) equals the reference in statistics, p-values and split parts.
- [ ] 1.8 Add figure 11 and the figures checks (design D5): the `FIGURES` entry, the draw and caption functions, `FIGURE_READING_RULES["span_split"]`, `_measurement` and `check_consistency`. Verify that:
  - `figures --check` on today's 0.13.0 outputs fails the consistency pre-pass naming a record without the split block and the setting expected, with nothing written;
  - a scratch script, with a `FigureReader` on a scratch tree of fixture records laid out at the registered paths, draws figure 11 alone. It shows 40 complement markers and four `latent_adjusted` rows labelled not available, and its caption's values equal the records;
  - the import-time checks pass.
- [ ] 1.9 Update README and AGENTS.md:
  - the pin `motco==0.18.0` and `ANALYSIS_MOTCO_MIN_VERSION` 0.18.0 (incl. the `analyze` step's "≥ 0.13.0");
  - `ANALYSIS_SPAN_SPLIT_TEST` and the `SPAN_SPLIT_*` constants under "flag and ask", with the multiplicity policy, the reading rules and the limitations (proposal);
  - the `latent_adjusted` exception;
  - the split columns of the comparison tables;
  - figure 11 and "eleven figures" (`01_cohorts` … `11_span_split`).

  Verify that every constant named in both files exists in `data.py` or `figures.py`, and that no text still says "ten figures" or "≥ 0.13.0" for `analyze`.

## 2. Gates (`nixos-home`, Nix shell)

- [ ] 2.1 Prepare the re-issue machine (design D12). On `nixos-home`, verify that each current `result.json` and `analysis_manifest.json` of both cell types hashes as on the cluster. Record the machine for the report.
- [ ] 2.2 Pin-move gate, previous pin (design D6.1). Create a scratchpad worktree of main `a5e9008`, run `uv sync --locked` there (MOTCO 0.13.0), and run `analyze --cell-type {astrocyte,microglia} --cohort all` on the real exports into `<scratch>/v0130/<cell type>`. Verify that the four manifests record MOTCO 0.13.0 and the real exports' hashes.
- [ ] 2.3 Pin-move gate, new pin without the test. From the change branch, run `_analyze_cohort` with the analyze parameters and `span_split_test=False` on the same exports into `<scratch>/v0180_none/<cell type>`. Verify that the records carry no `span_split` and MOTCO 0.18.0 under `software`.
- [ ] 2.4 Compare 2.2 with 2.3 using a scratch script. Pass condition: the four CSVs are byte-identical per cohort, and `result.json` and the manifest are equal after dropping `runtime_seconds`, `created_at`, `git`, `software.motco`, `cli_equivalent` and `measurement.span_split`. If anything differs, stop the change and report the file and value.
- [ ] 2.5 Compare 2.3 with the kept primaries (design D6.2). Record for the report: byte-identical, or the largest absolute differences and whether ranks, p-values and `analysis_metadata.csv` are equal.
- [ ] 2.6 Contract gate (design D6.3). Run the change branch's `analyze --cell-type {astrocyte,microglia} --cohort all` into `<scratch>/v0180_split/<cell type>` and compare it with 2.3. Pass condition:
  - the four CSVs are byte-identical;
  - `result.json` is equal after dropping `span_split`, `params.span_split_test` and `measurement.runtime_seconds`, every p-value and null summary included;
  - each block's permutations, seed and correction equal the analysis's, and its `span` dimensions are 5, 2, 2, 4.

  If anything differs, stop the change, name the cohort and the value, and report it to MOTCO.

## 3. Move the superseded outputs

- [ ] 3.1 Record into the scratchpad (design D7):
  - the sha256 of every file under both cell types' `pseudobulk/` and `motco_input/` and under `data/annotation/`;
  - every input hash recorded by `reports/figures/motco-case-study-figures-2026-11-02/figures_manifest.json`;
  - every `result.json` hash cited by `reports/motco-residual-correction-primary-2026-10-31.md`.

  Verify the file counts per directory.
- [ ] 3.2 Run `python3 scripts/move_superseded.py --change mtg-span-split-test` over the 15 MOTCO-dependent directories of both cell types. Verify that:
  - every moved file hashes as before and `MOVED.json` is written;
  - no MOTCO-dependent directory is left at its old path, and `results/figures/` was not moved;
  - every hash recorded in 3.1 from the frozen figure set and the 2026-10-31 report is found under `results/superseded/mtg-span-split-test/` or at its unchanged MOTCO-independent path.
- [ ] 3.3 Add `mtg-span-split-test/` (the 0.13.0 `df` outputs without the split test, both cell types) to the superseded entries in AGENTS.md and README. Verify that the documented paths match the moved tree.

## 4. Re-issue every MOTCO-dependent output (Nix shell, `n_jobs = 1`)

- [ ] 4.1 Run `analyze --cohort all` and `analyze --cell-type microglia --cohort all`. Verify that:
  - the four records carry `df` with κ 1.0690 / 1.0753 / 1.0731 / 1.0903, floor 2, MOTCO 0.18.0, seed 20261019, 999 permutations and a `span_split` block;
  - their CSVs and records are byte-identical / equal to the 2.6 scratch runs (the reproducibility rerun).
- [ ] 4.2 Run `rank-sensitivity --cohort all`. Verify that:
  - there are six forced ranks per cohort, each record carrying `df` and a block whose `span` dimension is its rank;
  - the manifest records the floor, `df` and the test;
  - the forced rank 2 of `cohort_50` equals the reference in statistics, p-values and split parts, while for `cohort_30` (rank 5, not in the ladder) the manifest records that the reference row stands in.
- [ ] 4.3 Run `axis-diagnostics` (astrocytes, with the sweep) and `axis-diagnostics --cell-type microglia`. Verify that the input-hash checks pass.
- [ ] 4.4 Run `preprocessing-sensitivity --variant all`, then `axis-diagnostics` per variant root. Verify that:
  - all eight variants are complete in both cohorts, every variant record carrying a block;
  - `adoption_check.csv` exists for the three candidates, and `adoption_precedence` is set in both trend candidates' comparison manifests.
- [ ] 4.5 Run `latent-sensitivity --cohort all` for both cell types, then the arms' diagnostics. Verify that:
  - `extended_covariates` carries a block with reduced rank 12;
  - `latent_adjusted` carries none, with its manifest's `not_available` set;
  - `mediator_check.csv` predates the arms;
  - each `extended_covariates` arm's rank and CV AUROC equal its primary's, or the difference is recorded as a finding.
- [ ] 4.6 Run `stage-de --omic both --n-cpus 10` and `stage-de-limma --omic both` for both cell types. Verify that the reference checks pass on the new primaries and that every concordance cell reads the new primary's `delta` status.
- [ ] 4.7 Run `exclusion-sensitivity --cohort all` for both cell types, then on each variant root `axis-diagnostics`, `stage-de --omic both --n-cpus 10` and `stage-de-limma --omic rna`. Verify that:
  - `flag_check.csv` predates the variant export;
  - the variant records carry a block;
  - the stage-DE manifests record the `donor_exclusion` block.
- [ ] 4.8 Run `positive-control --cohort all` and `positive-control --basis latent_adjusted --cohort all`. Verify that both bases' outputs exist for both cohorts and that the chain checks pass.
- [ ] 4.9 After every other rerun, run `figures --report motco-span-split-<date>`. Verify that:
  - the consistency pre-pass passes;
  - eleven figures are frozen under `reports/figures/motco-span-split-<date>/`;
  - figure 11 holds 40 complement p-values, and its caption's primary values equal the records.

## 5. Comparisons and consistency checks

- [ ] 5.1 Write and run the comparison script of design D9 over the 44 analysis pairs and the p-free outputs. Verify that it:
  - produces the report's outcome tables (superseded → new) and the split tables (every block, complement status, span concordance);
  - lists every difference beyond the expectation and every status change, each becoming a report finding.
- [ ] 5.2 Run the consistency check of design D10. Verify that it lists:
  - no manifest outside `results/superseded/` with a MOTCO version other than 0.18.0;
  - no `result.json` without `df`, floor 2 or a forced rank, and a block, except the four `latent_adjusted` arms with `not_available`.
- [ ] 5.3 Re-hash the files recorded in 3.1. Verify that every MOTCO-independent hash is unchanged and every cited hash still resolves.

## 6. Report

- [ ] 6.1 Write `reports/motco-span-split-<date>.md` from the generated files (design D11). It contains:
  - provenance: the machine, the predeclaration commit (1.1), the pin check (1.2), the gates (2.4, 2.6), the kept-output comparison (2.5) and the move (3.2);
  - MOTCO's validation with its numbers;
  - the policy and the reading quoted from the proposal;
  - the four primaries' split results in full and the per-family split tables;
  - one section per re-issued step with its outcome table;
  - the reading in the outcome that occurred (the two `cohort_30` complements as the results, `cohort_50` and every other analysis as dependence, every span concordance cell, `total` and the shares reported and not read);
  - the findings of 5.1;
  - the limitations of spec `motco-span-split-test` and every earlier limitation carried;
  - the eleven frozen SVGs, with captions pasted from the frozen `figures.md`.
- [ ] 6.2 Run the report's check script. Verify that:
  - every split p-value, F, dimension and share, and every statistic, p-value, rank, CV AUROC and κ of every `result.json`, is in the text;
  - every E1–E3 value and every concordance cell is in the text;
  - the captions equal `figures.md`, the frozen SVGs match the frozen manifest, and every input it records hashes as on disk.
- [ ] 6.3 Run `openspec validate mtg-span-split-test`. Verify that it passes (requirement-length warnings only, as elsewhere in the repository).

## 7. Cluster

- [ ] 7.1 On the cluster:
  - pull the branch and run `uv sync --locked` on the login node;
  - verify that `import motco` prints 0.18.0;
  - run the move of 3.2 over the MOTCO-dependent directories there.

  Verify that every moved file hashes as before and that `pseudobulk/` and `motco_input/` are untouched.
- [ ] 7.2 Copy every re-issued directory of both cell types from `nixos-home` to the cluster with rsync. Verify that every copied file's sha256 equals `nixos-home`'s.

## Workflow follow-up

- Commit the code, the report, the frozen figure set and the task records; open the PR against `main`.
- After merge, archive the change and sync its specs (new `motco-span-split-test`; modified `motco-trajectory-analysis`, `motco-rank-sensitivity`, `motco-preprocessing-sensitivity`, `motco-latent-adjustment-sensitivity`, `donor-exclusion-sensitivity`, `motco-release-migration`, `case-study-figures`).
- Tell MOTCO's roadmap (item 16 follow-up) the outcome, as a MOTCO-repo change.
- Update the case-study memory.

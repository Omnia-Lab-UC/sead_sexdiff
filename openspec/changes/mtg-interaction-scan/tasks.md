# Tasks

## 1. Pin and code

- [ ] 1.1 Record the predeclaration. Commit this change's planning artifacts (proposal, design, specs, tasks) on branch `mtg-interaction-scan` and push. Record the commit hash for the report as the commit that fixed the policy, the power statement and the reading. Verify with `git log` that it precedes every run of groups 2–4. No real-label scan may run before it.
- [ ] 1.2 Move the pin (design D1).
  - Run `uv add "motco==0.20.0"`.
  - Verify that `uv.lock` records motco 0.20.0 from the PyPI registry with wheel sha256 `d4248a8e19fd65f8…`.
  - Verify that `uv export --no-hashes` differs from main's output only in the motco line, and that scikit-learn is still `==1.8.0`.
  - Verify with a scratch script that the downloaded wheel's sha256 equals the lock's, and that its `motco/` tree equals `git -C ../MOTCO archive v0.20.0 src/motco` file by file.
  - Verify that `nix develop -c uv run python -c "import motco; print(motco.__version__)"` prints `0.20.0`.
  - Record the wheel hash, the tag commit (`5f6b07d`) and both comparisons for the report.
- [ ] 1.3 Build a synthetic fixture with a scratch script: random `rna` / `atac` blocks with the real cohorts' donor table, metadata and covariate shapes, written as an export directory in the scratchpad. Plant a 6-feature stage ramp in one sex, so the fixture has hits. Code checks then never show a real-label scan value before the gates. Verify that `analyze --input-dir <fixture>` runs under 0.20.0 on today's code (split block present, no scan block).
- [ ] 1.4 Add the constants and the record and table checks (design D2):
  - `ANALYSIS_MOTCO_MIN_VERSION = "0.20.0"`, `ANALYSIS_INTERACTION_SCAN`, the `INTERACTION_SCAN_*` policy constants and `INTERACTION_SCAN_POWER_STATEMENT`;
  - `check_motco_scan_alpha`;
  - `interaction_scan=` in every `TrajectoryAnalysisParams` construction;
  - `_interaction_scan_problem` in `_analyze_cohort` before any write;
  - `measurement.interaction_scan` in the manifest;
  - the degenerate-feature finding.

  Verify on the fixture that:
  - `analyze` records `params.interaction_scan = true` and a block whose per-block feature counts equal the fixture's block widths, with `df_residual` = donors − 12, and the manifest's permutations, seed and `df`, with `draws_match_rrpp = true`;
  - `interaction_scan.csv` is listed in the manifest's `outputs`;
  - the fixture's planted features are among the hits, the table's hits equal the record's, and the minimum adjusted p equals the global p;
  - patched `analyze_trajectories` results (block dropped; one hit removed from the record; table rows reordered) each stop the run naming the cohort and the entry, with nothing written;
  - a patched `INTERACTION_SCAN_ALPHA` stops the run before any input is read;
  - `INTERACTION_SCAN_POWER_STATEMENT`'s five sha256s and source hash equal `git -C ../MOTCO show v0.20.0:results/interaction-scan-2026-11-09/plasmode/settings.json` and that file's sha256.
- [ ] 1.5 Add the summary section and the CLI flag (design D3: `render_analysis_summary`, `_cli_equivalent`). Verify that:
  - `summary.md` on the fixture shows the scan table, the hits with their three differences and the shares labelled descriptive;
  - `render_analysis_summary` on each current 0.18.0 `result.json` of both cell types (no scan block) is byte-identical to its `summary.md` on disk;
  - running the fixture's recorded `cli_equivalent` with `motco analyze` gives the same statistics, p-values, rank, `span_split` block and `interaction_scan` block, and an equal `interaction_scan.csv`.
- [ ] 1.6 Add the scan columns and the dependence findings (design D3). Verify on the fixture that:
  - `preprocessing-sensitivity --variant keep_chrx`, `latent-sensitivity` and `exclusion-sensitivity` (with a fixture flag table) write their tables with the six scan columns, every cell equal to its row's `result.json`;
  - a hand-made pair of records with different scan statuses produces the dependence finding naming both p-values;
  - the split columns and the stage-DE `delta` cells are unchanged.
- [ ] 1.7 Add the reference check and the identity checks (design D4, D5). Verify on the fixture that:
  - `rank-sensitivity` and `stage-de` against a reference made with `interaction_scan=False` stop naming the reference `result.json` and `analyze`, with nothing written;
  - against the scan reference, every forced `interaction_scan.csv` is byte-identical to the reference's and the cross-rank table's scan cells are equal across rows;
  - `latent_adjusted`'s scan equals `extended_covariates`'s;
  - a forced record with a patched scan p-value, and a patched arm, each stop their step naming the analysis and the entry, with nothing published.
- [ ] 1.8 Add the stage-DE scan cells (design D6): the two rules with their import-time ties, the summary and comparison columns, the hits from `_stage_de_prepare`, the manifests and the Markdown. Verify on the fixture that:
  - `stage-de --omic both` and `stage-de-limma --omic both` fill `scan_hits`, `scan_hits_significant`, `scan_detected` and `scan_concordance` on the rule rows only, and `method_comparison.csv`'s `interaction_joint` rows carry the two scan cells;
  - a hit absent from the exported tables stops the step naming the feature;
  - discordant cells produce the findings of the specs;
  - `positive-control` (E1), figure 10's caption builder and the limma step's DESeq2 binding read the same values from the widened summaries as from the old ones.
- [ ] 1.9 Add figure 12 and the figures checks (design D7): the `FIGURES` entry, `FIGURE_SCAN_TOP_N`, the draw and caption functions, `FIGURE_READING_RULES["interaction_scan"]`, `_measurement` and `check_consistency`. Verify that:
  - `figures --check` on today's 0.18.0 outputs fails the consistency pre-pass naming a record without the scan block, with nothing written;
  - a scratch script draws figure 12 alone from a `FigureReader` on a scratch tree of fixture records laid out at the registered paths. It shows 44 global p-values, 16 of them marked identical to their reference, and its caption's values equal the records;
  - the import-time checks pass.
- [ ] 1.10 Update README and AGENTS.md:
  - the pin `motco==0.20.0` and `ANALYSIS_MOTCO_MIN_VERSION` 0.20.0 (incl. the `analyze` step's "≥ 0.18.0");
  - `ANALYSIS_INTERACTION_SCAN`, the `INTERACTION_SCAN_*` constants and the power statement under "flag and ask", with the policy, the reading rules and the limitations (proposal);
  - "five tests" in the split-test policy;
  - the scan columns, the identity checks and the stage-DE scan cells;
  - figure 12 and "twelve figures" (`01_cohorts` … `12_interaction_scan`).

  Verify that every constant named in both files exists in `data.py` or `figures.py`, and that no text still says "eleven figures", "four tests read at α" or "≥ 0.18.0" for `analyze`.

## 2. Gates (`nixos-home`, Nix shell)

- [ ] 2.1 Prepare the re-issue machine (design D14). On `nixos-home`, verify that each current `result.json` and `analysis_manifest.json` of both cell types hashes as on the cluster. Record the machine for the report.
- [ ] 2.2 Pin-move gate, previous pin (design D8.1). Create a scratchpad worktree of main `5f70673`, run `uv sync --locked` there (MOTCO 0.18.0), and run `analyze --cell-type {astrocyte,microglia} --cohort all` on the real exports into `<scratch>/v0180/<cell type>`. Verify that the four manifests record MOTCO 0.18.0 and the real exports' hashes.
- [ ] 2.3 Pin-move gate, new pin without the scan. From the change branch, run `_analyze_cohort` with the analyze parameters and `interaction_scan=False` on the same exports into `<scratch>/v0200_none/<cell type>`. Verify that the records carry a `span_split` block, no `interaction_scan`, and MOTCO 0.20.0 under `software`.
- [ ] 2.4 Compare 2.2 with 2.3 using a scratch script. Pass condition: the four CSVs are byte-identical per cohort, and `result.json` and the manifest are equal after dropping `runtime_seconds`, `created_at`, `git`, `software.motco`, `cli_equivalent` and `measurement.interaction_scan`, the `span_split` block included. If anything differs, stop the change and report the file and value.
- [ ] 2.5 Compare 2.2 and 2.3 with the kept primaries (design D8.2). Record for the report: byte-identical, or the files and largest differences as findings.
- [ ] 2.6 Contract gate (design D8.3). Run the change branch's `analyze --cell-type {astrocyte,microglia} --cohort all` into `<scratch>/v0200_scan/<cell type>` and compare it with 2.3. Pass condition:
  - the four CSVs are byte-identical, and only `v0200_scan` has `interaction_scan.csv`;
  - `result.json` is equal after dropping `interaction_scan`, `params.interaction_scan` and `measurement.runtime_seconds`, every p-value, null summary and the `span_split` block included;
  - each block's permutations, seed and correction equal the analysis's, and its `df_residual` is n − 12.

  If anything differs, stop the change, name the cohort and the value, and report it to MOTCO.
- [ ] 2.7 Power-statement binding (design D8.4). Verify that the five sha256s of `results/motco_input/cohort_30/` equal `INTERACTION_SCAN_POWER_STATEMENT`'s, and record the comparison for the report.

## 3. Move the superseded outputs

- [ ] 3.1 Record into the scratchpad (design D9):
  - the sha256 of every file under both cell types' `pseudobulk/` and `motco_input/` and under `data/annotation/`;
  - every input hash recorded by `reports/figures/motco-span-split-2026-11-06/figures_manifest.json`;
  - every `result.json` hash cited by `reports/motco-span-split-2026-11-06.md`.

  Verify the file counts per directory.
- [ ] 3.2 Run `python3 scripts/move_superseded.py --change mtg-interaction-scan` over the 15 MOTCO-dependent directories of both cell types. Verify that:
  - every moved file hashes as before and `MOVED.json` is written;
  - no MOTCO-dependent directory is left at its old path, and `results/figures/` was not moved;
  - every hash recorded in 3.1 from the frozen figure set and the 2026-11-06 report is found under `results/superseded/mtg-interaction-scan/` or at its unchanged MOTCO-independent path.
- [ ] 3.3 Add `mtg-interaction-scan/` (the 0.18.0 outputs without the scan, both cell types) to the superseded entries in AGENTS.md and README. Verify that the documented paths match the moved tree.

## 4. Re-issue every MOTCO-dependent output (Nix shell, `n_jobs = 1`)

- [ ] 4.1 Run `analyze --cohort all` and `analyze --cell-type microglia --cohort all`. Verify that:
  - the four records carry `df` with κ 1.0690 / 1.0753 / 1.0731 / 1.0903, floor 2, MOTCO 0.20.0, seed 20261019, 999 permutations, a `span_split` block and an `interaction_scan` block;
  - their files are byte-identical / equal to the 2.6 scratch runs (the reproducibility rerun).
- [ ] 4.2 Run `rank-sensitivity --cohort all`. Verify that:
  - there are six forced ranks per cohort, each record carrying `df`, a split block whose `span` dimension is its rank, and a scan identical to the reference's (the D5 check passed);
  - the forced rank 2 of `cohort_50` equals the reference in statistics, p-values and split parts.
- [ ] 4.3 Run `axis-diagnostics` (astrocytes, with the sweep) and `axis-diagnostics --cell-type microglia`. Verify that the input-hash checks pass.
- [ ] 4.4 Run `preprocessing-sensitivity --variant all`, then `axis-diagnostics` per variant root. Verify that:
  - all eight variants are complete in both cohorts, every variant record carrying both blocks;
  - `adoption_check.csv` exists for the three candidates, and `adoption_precedence` is set in both trend candidates' comparison manifests.
- [ ] 4.5 Run `latent-sensitivity --cohort all` for both cell types, then the arms' diagnostics. Verify that:
  - both arms carry a scan block with `df_residual` n − 14;
  - the D5 arm identity passed;
  - `latent_adjusted` carries no split block, with its manifest's `not_available` set;
  - `mediator_check.csv` predates the arms.
- [ ] 4.6 Run `stage-de --omic both --n-cpus 10` and `stage-de-limma --omic both` for both cell types. Verify that:
  - the reference checks pass on the new primaries;
  - every `delta` cell reads the new primary's status and every scan cell reads its hits;
  - whether the tables hash as the superseded ones (expected on the same machine) is recorded for the report.
- [ ] 4.7 Run `exclusion-sensitivity --cohort all` for both cell types, then on each variant root `axis-diagnostics`, `stage-de --omic both --n-cpus 10` and `stage-de-limma --omic rna`. Verify that:
  - `flag_check.csv` predates the variant export;
  - the variant records carry both blocks;
  - the variant stage-DE scan cells read the variant's hits.
- [ ] 4.8 Run `positive-control --cohort all` and `positive-control --basis latent_adjusted --cohort all`. Verify that both bases' outputs exist for both cohorts and that the chain checks pass.
- [ ] 4.9 After every other rerun, run `figures --report motco-interaction-scan-<date>`. Verify that:
  - the consistency pre-pass passes;
  - twelve figures are frozen under `reports/figures/motco-interaction-scan-<date>/`;
  - figure 12 holds 44 global p-values, and its caption's primary values equal the records.

## 5. Comparisons and consistency checks

- [ ] 5.1 Write and run the comparison script of design D11 over the 44 analysis pairs and the p-free outputs. Verify that it:
  - produces the report's outcome tables (superseded → new) and the scan tables (every block, status, hits, dependence, identities, concordance cells);
  - lists every difference beyond byte-identity and every status change, each becoming a report finding.
- [ ] 5.2 Run the consistency check of design D12. Verify that it lists:
  - no manifest outside `results/superseded/` with a MOTCO version other than 0.20.0;
  - no `result.json` without `df`, floor 2 or a forced rank, and the scan block; the split block missing only from the four `latent_adjusted` arms with `not_available`.
- [ ] 5.3 Re-hash the files recorded in 3.1. Verify that every MOTCO-independent hash is unchanged and every cited hash still resolves.

## 6. Report

- [ ] 6.1 Write `reports/motco-interaction-scan-<date>.md` from the generated files (design D13). It contains:
  - provenance: the machine, the predeclaration commit (1.1), the pin check (1.2), the gates (2.4, 2.6), the kept-output comparison (2.5), the power-statement binding (2.7) and the move (3.2);
  - MOTCO's validation with its numbers;
  - the policy and the reading quoted from the proposal;
  - the four primaries' scan results in full and the per-family scan tables;
  - the concordance cells with their overlap;
  - one section per re-issued step with its outcome table;
  - the reading in the outcome that occurred (the two `cohort_30` global p-values as the results, the power statement or "uninformative" as the specs give it, every other analysis as dependence or identity);
  - the findings of 5.1;
  - the limitations of spec `motco-interaction-scan` and every earlier limitation carried;
  - the twelve frozen SVGs, with captions pasted from the frozen `figures.md`.
- [ ] 6.2 Run the report's check script. Verify that:
  - every scan value, and every statistic, p-value, rank, CV AUROC, κ and split value of every `result.json`, is in the text;
  - every concordance cell and every E1–E3 value is in the text;
  - the captions equal `figures.md`, the frozen SVGs match the frozen manifest, and every input it records hashes as on disk.
- [ ] 6.3 Run `openspec validate mtg-interaction-scan`. Verify that it passes.

## 7. Cluster

- [ ] 7.1 On the cluster:
  - pull the branch and run `uv sync --locked` on the login node;
  - verify that `import motco` prints 0.20.0;
  - run the move of 3.2 over the MOTCO-dependent directories there.

  Verify that every moved file hashes as before and that `pseudobulk/` and `motco_input/` are untouched.
- [ ] 7.2 Copy every re-issued directory of both cell types from `nixos-home` to the cluster with rsync. Verify that every copied file's sha256 equals `nixos-home`'s.

## Workflow follow-up

- Commit the code, the report, the frozen figure set and the task records; open the PR against `main`.
- After merge, archive the change and sync its specs (new `motco-interaction-scan`; modified `motco-trajectory-analysis`, `motco-span-split-test`, `motco-rank-sensitivity`, `motco-preprocessing-sensitivity`, `motco-latent-adjustment-sensitivity`, `donor-exclusion-sensitivity`, `pseudobulk-stage-de`, `pseudobulk-stage-de-limma`, `motco-release-migration`, `case-study-figures`).
- Tell MOTCO's roadmap (item 16's follow-up thread 3: the scan applied to SEA-AD) the outcome, as a MOTCO-repo change.
- Update the case-study memory.

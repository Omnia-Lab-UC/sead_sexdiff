# Design

## Context

- **MOTCO 0.13.0** (tag `v0.13.0`, commit `5758c29`; PyPI wheel sha256 `614c2e9c…`) adds `TrajectoryAnalysisParams.residual_correction`: `None` (the default), `"df"` or `"leverage"`.
  - Under `df`, RRPP permutes the reduced model's residuals multiplied by κ = √(n / (n − r)), with r the numerical rank of the reduced model. Here r is the intercept, sex and two stage columns plus the covariate columns: 10 with the six analysis covariates and 12 with the eight of the latent arms.
  - The observed statistics, the latent space, the CV rank and the stage means are unchanged. Only `p_values`, `null_summary` and `config_spectrum.permutation_pooled_eigengap` move, because all three are drawn from the null.
  - The record states the applied correction under `measurement.residual_correction` (`kind`, `n`, `reduced_rank`, `kappa`) and `params.residual_correction = "df"`. Both are omitted while the setting is `None`, so an uncorrected 0.13.0 analysis reproduces 0.10.0 apart from `software`.
  - Between `v0.10.0` and `v0.13.0` only `analysis.py`, `cli.py` and `stats/permutation.py` changed on the analysis path, all for this setting. `pyproject.toml` did not change, so the dependency pins, including scikit-learn's exact pin, are the same.
  - The CLI gains `motco analyze --residual-correction {df,leverage}`.
- **Predicted κ** for the primaries: astrocyte `cohort_30` (n = 80) 1.0690, `cohort_50` (n = 74) 1.0753; microglia `cohort_30` (n = 76) 1.0731, `cohort_50` (n = 63) 1.0903. With r = 12 for the arms: 1.0847, 1.0925, 1.0897 and 1.1114.
- **Where the package builds analyses.** `_analyze_cohort` is the one place an analysis is run and recorded. It is called from `analyze`, `rank-sensitivity`, `preprocessing-sensitivity`, `latent-sensitivity` and `exclusion-sensitivity`, each of which builds `TrajectoryAnalysisParams(... min_components=ANALYSIS_RANK_FLOOR)`. It writes into `<cohort>.tmp` only after `analyze_trajectories` returns, so a check placed between the call and the write leaves nothing behind. `_check_reference` serves the five readers (the rank, preprocessing, latent and exclusion steps, and `_stage_de_prepare`, which `stage-de-limma` reuses). It already compares the export hashes, the MOTCO version, the seed, the permutations, the floor and `component_selection = cv`.
- **What is where.** The current outputs (MOTCO 0.10.0, floor 2) were made on `nixos-home` on 2026-10-26 and copied to the cluster, which is the canonical copy. That copy holds both cell types' current directories plus `results/superseded/{mtg-rank-floor-primary,mtg-astrocyte-latent-adjustment-sensitivity}/`. `nixos-dell-laptop` holds only an older astrocyte subset.
- **The pin today** is `motco = { git = "https://github.com/tomszar/MOTCO", rev = "v0.10.0" }`. The repository now lives under `Omnia-Lab-UC`, and the old path only works through GitHub's redirect.

## Goals / Non-Goals

**Goals:**
- The correction is set in one constant, passed explicitly, verified in every record, recorded wherever it acts, and checked by every reader.
- Two same-machine gates prove, before anything is moved, that the pin move changes nothing without the correction and that the correction changes only the null.
- After the re-issue, every output outside `results/superseded/` rests on MOTCO 0.13.0 with `df` at floor 2, and every superseded byte stays resolvable by its hash.
- The PyPI pin is bound to the MOTCO release tag by a recorded byte comparison.

**Non-Goals:**
- A flag for the correction, `leverage`, or an uncorrected arm. The superseded outputs are that record.
- Rerunning the MOTCO-independent steps or moving their outputs.
- Bringing `nixos-office` or the laptop up to date beyond what the re-issue machine needs; they re-sync from the cluster when next used.

## Decisions

**D1. Pin `motco==0.13.0` from PyPI, bound to the tag.**
- `pyproject.toml` lists `motco==0.13.0`, and the `[tool.uv.sources]` motco entry is removed (the section with it, since it is the only entry). `uv lock` then records the registry source and the wheel and sdist sha256.
- A scratch script makes two checks and records them:
  - it downloads the wheel named in `uv.lock`, checks its sha256 against the lock, and compares its `motco/` tree file by file with `git -C ../MOTCO archive v0.13.0 src/motco`, which a check on 2026-10-31 found byte-identical;
  - it compares `uv export --no-hashes` before and after the lock change, where only the motco line may differ.
- `AGENTS.md` and `README.md` state the rule: an exact PyPI version, moved with `uv add "motco==X.Y.Z"`, its wheel checked against the release tag of `https://github.com/Omnia-Lab-UC/MOTCO`, never a git or path dependency.
- *Alternative: keep the git tag (with the URL updated).* Rejected by the user's decision. The lock's wheel hash pins exact bytes as firmly as a commit, PyPI versions cannot be overwritten, and the cluster's `uv sync` stops cloning from GitHub. The tag comparison keeps the commit link a git pin would give.

**D2. One constant, passed explicitly, verified in the record.**
- `ANALYSIS_RESIDUAL_CORRECTION = "df"` sits next to `ANALYSIS_RANK_FLOOR`. Its comment cites MOTCO's verdict (`adopt_df`, `corrected-permutation-2026-10-30.md`) and says it is MOTCO's validated opt-in, not its default, and never a flag. `ANALYSIS_MOTCO_MIN_VERSION = "0.13.0"`.
- The five constructions pass `residual_correction=ANALYSIS_RESIDUAL_CORRECTION`. That includes the forced sweep and both latent arms, whose definition is "everything else equals `analyze`".
- `_analyze_cohort` checks the record after `analyze_trajectories` and before anything is written:
  - with the setting, `record["measurement"]["residual_correction"]["kind"]` must equal it;
  - without it (the gate runs by hand), the record must carry no correction.
  
  The error names the recorded and the expected setting.
- The manifest's `measurement` gains `residual_correction`: the record's block (kind, n, reduced rank, κ), or `null`.
- `_cli_equivalent` adds `--residual-correction df` to CV and forced runs alike, since the setting acts on RRPP in both.
- `render_analysis_summary` adds one line, "RRPP null: residual correction `df` (κ = …; n = …, reduced rank …)", when the record carries a correction. Without one the output is unchanged, so superseded records render as before.
- There is no check against MOTCO's default, unlike the floor (D1 of `mtg-rank-floor-primary`): MOTCO's default is off by design, and the record check above is what catches a MOTCO that ignored the setting.
- *Alternative: correct only the primaries.* Rejected. Every sensitivity step compares rejection statuses with the primary under "everything else equal", and mixing nulls would turn every comparison into a null comparison.
- *Alternative: `leverage`.* Rejected. MOTCO's tie-break found no paired advantage, and the predeclared rule chose `df`.

**D3. The reference check reads the correction.**
- `_check_reference` reads `record["params"].get("residual_correction")` (absent means none) after the floor check. When it differs from `ANALYSIS_RESIDUAL_CORRECTION`, the check stops, naming the file, both values and `analyze`.
- The MOTCO-version check already refuses the 0.10.0 references. The new check also catches a 0.13.0 analysis made without the correction, such as the gate's scratch runs or `motco analyze` by hand.
- The rank-sweep manifest records `residual_correction` beside `rank_floor`, and the reference's setting beside `min_components`.
- `axis-diagnostics` and `positive-control` do not use `_check_reference`. They bind to the analyses by hash and chain checks, which a moved-and-rerun primary satisfies only when re-issued. D8 covers their versions.

**D4. Two gates on the re-issue machine, before the move.**
1. *Pin-move gate* (spec `motco-release-migration`, existing):
   - Create a worktree of the pre-change commit (the branch point on `main`, `fddf295` today), run `uv sync --locked` there (MOTCO 0.10.0), and run `analyze --cell-type {astrocyte,microglia} --cohort all --input-dir <real export> --out-dir <scratch>/v0100/<cell type>`.
   - From the change branch (MOTCO 0.13.0), a scratch script runs `_analyze_cohort` with the analyze parameters and `residual_correction=None` into `<scratch>/v0130_none/<cell type>`.
   - Pass condition: the four CSVs per cohort are byte-identical, and `result.json` and the manifest are equal after dropping `runtime_seconds`, `created_at`, `git`, `software.motco`, `cli_equivalent` and the new `measurement.residual_correction` (`null`).
   - Also recorded, not a gate: the comparison of `v0130_none` with the kept primaries. Byte-identical is expected on `nixos-home`, which made them.
2. *Contract gate* (spec `motco-release-migration`, added):
   - The change branch's own `analyze --cell-type … --cohort all --out-dir <scratch>/v0130_df/<cell type>` is compared with `v0130_none`.
   - Pass condition: the four CSVs are byte-identical; `result.json` is equal after dropping `p_values`, `null_summary`, `config_spectrum.permutation_pooled_eigengap`, `params.residual_correction`, `measurement.residual_correction` and `measurement.runtime_seconds`; and each κ equals √(n / (n − r)) with r = 10, to 1e-12.
   - A failure stops the change and is reported to MOTCO, because the contract would not hold on these data.
   - `v0130_df` also serves as the reproducibility rerun for the primaries made in D6.

- *Why the gates run before the move:* a failed gate leaves `results/` exactly as it was.

**D5. Superseded layout.**
- `python3 scripts/move_superseded.py --change mtg-residual-correction-primary` runs with the directory list of the 2026-10-26 report's section 1.5. A directory that does not exist is skipped and listed.
  - Astrocyte: `motco_analysis`, `motco_rank_sensitivity`, `motco_diagnostics`, `motco_sensitivity`, `motco_latent_sensitivity`, `stage_de`, `stage_de_limma`, `motco_exclusion_sensitivity`.
  - Microglia, under `microglia/`: `motco_analysis`, `motco_diagnostics`, `motco_latent_sensitivity`, `stage_de`, `stage_de_limma`, `motco_exclusion_sensitivity`, `positive_control`.
- It runs on the re-issue machine and then on the cluster. Before the move, the sha256 of every file under both cell types' `pseudobulk/` and `motco_input/` and under `data/annotation/` is recorded.
- The stage-DE directories move and are refitted, as in 2026-10-26. Their concordance cells read the primary's `delta` status, and their manifests bind the primary's hashes. On the machine that made them, the refits are expected to be byte-identical.

**D6. Rerun order and commands.** These are exactly the commands of the 2026-10-26 report's section 1.5 re-issue block, run through `nix develop -c uv run sead-sexdiff …` with `--n-jobs 1` and `--n-cpus 10` for `stage-de`. The order is the same:
1. `analyze`;
2. `rank-sensitivity` (astrocytes);
3. `axis-diagnostics`;
4. `preprocessing-sensitivity --variant all` with `axis-diagnostics` per variant root (astrocytes);
5. `latent-sensitivity` with its diagnostics;
6. `stage-de --omic both` and `stage-de-limma --omic both`;
7. `exclusion-sensitivity`, then `axis-diagnostics`, `stage-de --omic both` and `stage-de-limma --omic rna` on its variant root;
8. `positive-control` on both bases.

Steps 2–7 depend only on step 1 of their cell type, and step 8 on steps 3, 5 and 6 of both cell types. The astrocyte and microglia chains run side by side. Every step writes fresh directories, so no `--force` is needed.

**D7. Comparison script (superseded → corrected).** It runs over every analysis directory pair: the primaries, the 12 forced ranks, the 16 variant analyses, the 8 arms and the 4 exclusion variants. For each pair:
- **the D4.2 contract comparison:** byte-identical CSVs, and records equal outside the null keys;
- **κ check:** κ against √(n / (n − r)), with r = 4 + the number of covariate columns;
- **per statistic:** the superseded and corrected p, the direction of the change, and both rejection statuses, with any status change listed.

The script also checks byte-identity of every output that reads no p-value: the axis-diagnostic CSVs, `adoption_check.csv`, the adoption precedence, `mediator_check.csv`, `flag_check.csv`, `axis_check.csv`, the E1–E3 values of `positive_control.csv` on both bases, and every stage-DE and limma per-feature table and summary count. It re-reads every concordance cell and `methods_agree`.

It produces the report's outcome tables (rank, CV AUROC, statistics, p-values, rejections per analysis; adoption check and precedence; concordance cells; E1–E3 and the control outcomes). Outside the D4.2 gate on the primaries, a difference is a finding, never a gate.

The expectation, stated in advance and recorded as observed: with the same permutation indices and κ > 1, `delta`'s null draws widen, so its corrected p is not below the superseded one. `angle` and `shape` are expected to move the same way in most cases.

**D8. Version-consistency check.** It reuses the 2026-10-26 script with the new expectations. Every `*manifest*.json` outside `results/superseded/` that records a MOTCO version must record `0.13.0`. Every `result.json` must carry `params.residual_correction == "df"` and either `min_components == 2` or a forced rank. The check runs before the report is assembled and again before the commit.

**D9. Report.** `reports/motco-residual-correction-primary-<date>.md` is assembled by a scratch script, like the 2026-10-26 report:
- **Provenance:** the pin, the wheel check, both gates, the kept-output comparison, the move.
- **Reason:** MOTCO's three reports, quoted with their numbers. It states that MOTCO's roadmap withheld the 2026-10-26 primaries' p-values from 2026-10-30 until this re-issue, which no earlier report of this repository records.
- **The corrected primaries per cell type and cohort:** `summary.md` content with κ, plus the eigengap and the null summaries beside `angle`.
- **One section per re-issued step:** its generated Markdown and the D7 outcome table.
- **The predeclared reading** in the outcome that occurred.
- **Findings and limitations.**

A check script finds every statistic, p-value, rank, CV AUROC and κ of every re-issued `result.json`, every E1–E3 value and every concordance cell in the text. Large per-feature tables are referenced, not pasted.

**D10. Machine.** The re-issue runs on `nixos-home`, which made the superseded outputs. There G1's kept-output comparison, the D7 comparisons and the stage-DE refits are expected to be byte-identical.

If `nixos-home` is unavailable when the change is applied, it runs on `nixos-dell-laptop` after the current cluster directories of both cell types (not `superseded/`) are rsynced there, with sha256 verified:
- the gates are unaffected, because both compare runs made on the same machine;
- the D7 reproductions are then expected at rounding level, and the stage-DE tables beyond rounding, as in finding 11 of the 2026-10-26 report. These are recorded as findings.

The report names the machine.

## Risks / Trade-offs

- [A rejection status changes and is read as a discovery] → The reading is written in the proposal before any run. The corrected status is the result under the binding rules, and the superseded one stands beside it as dependence on the null correction. A new rejection is not expected under a correction that widens the null, and the D7 direction table shows whether it happened.
- [MOTCO's "only the null moves" contract fails on these data] → The contract gate stops the change before anything moves, and the case goes back to MOTCO.
- [The PyPI wheel differs from the tag, or re-locking moves another package] → D1's checks stop the pin move before any run.
- [The corrected null's size was measured at six covariates on the two `cohort_30` tables only] → It is inferred elsewhere (`cohort_50`, eight-covariate arms, exclusion variants with smaller n, forced ranks, preprocessing variants) from the mechanism, which depends on n and the reduced rank only. This is stated per step as a limitation, and no Type I claim is made beyond MOTCO's measurement.
- [The `latent_adjusted` arm is inadmissible for `angle` and `shape` (`residualized_basis`)] → Its p-values are reported and not read (spec). The arm stays registered, and the restated control reads none of them.
- [A step silently reads a superseded file] → The old directories are moved first, so a missing input is an error, and D8 lists any stale manifest.
- [Compute] → As on 2026-10-26: about 50 CV analyses of 1–2 min plus the forced, variant, stage-DE and limma runs; a few hours in the background, with logs in the scratch directory.

## Migration Plan

1. Pin and code (D1–D3) and docs; commit.
2. Gates (D4) on the re-issue machine; stop on failure.
3. Hash the MOTCO-independent outputs; move (D5).
4. Reruns (D6), comparisons (D7), version check (D8), re-hash.
5. Report (D9), check, commit.
6. On the cluster: pull, `uv sync --locked` (now from PyPI), the same move, then copy the re-issued directories with rsync, verified by sha256.

Rollback: revert the pin and code commits, and move the directories back from `results/superseded/mtg-residual-correction-primary/` with `MOVED.json` on both machines.

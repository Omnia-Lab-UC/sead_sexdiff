# Design

## Context

- **MOTCO 0.20.0:** tag `v0.20.0`, commit `5f6b07d`. On PyPI the wheel sha256 is `d4248a8e19fd65f86db07e3296a320dac837565889d3db0f94d57b6f4d191823` and the sdist's is `8c2588056c9f5c37…`, both uploaded 2026-11-10. On 2026-11-10 the downloaded wheel's `motco/` tree (52 files) was equal to `git -C ../MOTCO archive v0.20.0 src/motco`.
  - **Between `v0.18.0` and `v0.20.0`:**
    - `analysis.py`: the scan opt-in, plus `_split_feature_names` factored out of `_loadings_frame` (same output);
    - `cli.py`;
    - `stats/interaction.py`: the scan, plus the split test's model spaces and input checks moved into shared helpers;
    - `stats/__init__.py`;
    - `simulations/*`: 0.19.0's harness-only `latent_coordinates`.

    `stats/permutation.py`, `stats/pls.py`, `stats/trajectory.py`, `stats/design.py` and `pyproject.toml`'s pins are unchanged.
  - **MOTCO's contract:** without the setting, the resolved parameters and every output equal 0.19.0's (`to_dict` drops `interaction_scan` while false).
- **The scan as `analyze_trajectories` runs it** (`_interaction_scan_outputs`):
  - **Input:** `interaction_feature_scan` runs on the analysis's pooled-standardized joint features `X`, with the full and reduced trajectory models (covariates included), `permutations`, `seed` and the residual correction.
  - **Draws:** it takes them from its own `np.random.default_rng(seed)`, one `rng.permutation(n)` per draw, which is RRPP's serial index sequence. So RRPP's p-values do not move.
  - **Statistics:** per feature, the interaction F with `df_interaction` = 2 and `df_residual` = n − 6 − covariates. The global p-value is `(1 + #{draws with null max F ≥ max F}) / (1 + B)`. Each feature's Westfall–Young step-down adjusted p-value is the running maximum down the F order, so the top feature's equals the global p-value.
  - **Degenerate features:** a feature with zero full-model residual variance is degenerate. It gets NaN and never enters a maximum.
  - **The record block `interaction_scan`** holds:
    - `n_features`, `n_degenerate`, `df_interaction`, `df_residual`;
    - `max_f`, `max_feature`, `max_block`, `p_value_global`;
    - `hit_alpha` (MOTCO's `INTERACTION_SCAN_ALPHA` = 0.05);
    - `hits`: block, feature, F and adjusted p, sorted by F descending;
    - `blocks.<name>`: `n_features`, `n_degenerate`, `hits`, `null_max_share`;
    - `permutations`, `seed`, `residual_correction`, `draws_match_rrpp`.
  - **The table** `TrajectoryAnalysisResult.interaction_scan`, written by `write_analysis_outputs` as `interaction_scan.csv`, has the columns `block, feature, ss_interaction, rss, f, p_value_adjusted, degenerate, difference[0], difference[1], difference[2]`. A difference is the second group minus the first at the covariate reference row, in standardized units.
  - **Availability:** the scan reads neither the latent basis nor the bootstrap, so it runs with the latent adjustment.
- **Where the package builds analyses.**
  - **Constructions:** the five `TrajectoryAnalysisParams` constructions of the split-test change: `run_analysis`, `run_rank_sensitivity`, `run_preprocessing_sensitivity`, the arm registry in `run_latent_sensitivity` and `run_exclusion_sensitivity`, plus the fixture/gate helper at the end of `data.py`.
  - **Writes:** `_analyze_cohort` checks the correction and the split block between the MOTCO call and `write_analysis_outputs`. The scan check goes in the same place. The manifest's `outputs` already lists every file `write_analysis_outputs` returns, so `interaction_scan.csv` is listed without new code.
  - **Comparison rows:** `_sensitivity_row` builds every comparison row, and `_split_columns` / `_render_split_table` / `_split_dependence_findings` give the pattern the scan columns follow.
  - **Reference checks:** `_check_reference` serves the four sensitivity steps and `_stage_de_prepare`.
  - **Stage DE:** `_stage_de_summary_rows` fills the DESeq2 rule row. The limma step fills its rule row and `method_comparison.csv` from the same reference record.
- **Current outputs:**
  - the 2026-11-06 re-issue (MOTCO 0.18.0, `df`, floor 2, split test), made on `nixos-home`, rsynced to the cluster (1,165 files, sha256-verified);
  - the frozen set `reports/figures/motco-span-split-2026-11-06/`;
  - main at `5f70673`.

  The astrocyte `cohort_30` export's five files hash exactly as MOTCO's plasmode settings record at `v0.20.0` (`results/interaction-scan-2026-11-09/plasmode/settings.json`, `span_split.plasmode.export.files`; checked 2026-11-10).
- **Known when written:** every exported stage-DE interaction count of the four primaries is 0 (DESeq2 `factor_interaction` union and limma `interaction_joint`, both omics). Every `delta` cell is concordant.

## Goals / Non-Goals

**Goals:**
- The scan is set in one constant, passed explicitly, verified in every record and table and recorded wherever it acts, and every reader checks it.
- The policy, the power statement and the reading are fixed in code constants and in the specs before any real-label scan runs.
- The two identities the scan's construction implies (forced ranks, latent arms) are checked, not assumed.
- After the re-issue, every output outside `results/superseded/` rests on MOTCO 0.20.0, and every pre-existing value is shown to reproduce.

**Non-Goals:**
- A flag for the scan, an analysis without it, or a scan run outside `analyze_trajectories` (in particular, no direct call to `interaction_feature_scan`).
- Per-block thresholds, a moderated F, MOTCO's `return_null`, or any per-stage test of the differences.
- Gene annotation of ATAC hits, enrichment, or any downstream analysis of hits.
- A positive-control expectation, adoption rule or diagnostic variable that reads the scan.
- Bringing `nixos-dell-laptop` up to date; it re-syncs from the cluster when next used.

## Decisions

**D1. Pin `motco==0.20.0` from PyPI, bound to the tag.**
- **The move:** `uv add "motco==0.20.0"`. Then check:
  - `uv.lock` records the PyPI wheel with the sha256 above;
  - `uv export --no-hashes` differs from main's only in the motco line;
  - the wheel's `motco/` tree equals the tag's `src/motco` file by file (rechecked from the lock's download, not the pre-proposal copy).
- **Alternative rejected:** pinning 0.19.0 first. It carries nothing the package uses and would need its own re-issue.

**D2. One constant, the policy constants, the record and table checks.**
- **Constants in `data.py`, beside the split constants:**
  - `ANALYSIS_MOTCO_MIN_VERSION = "0.20.0"` and `ANALYSIS_INTERACTION_SCAN = True`;
  - the policy:
    - `INTERACTION_SCAN_STATISTIC = "p_value_global"` (confirmatory, read at `ANALYSIS_ALPHA` under `ANALYSIS_REJECTION_RULE`);
    - `INTERACTION_SCAN_HIT_STATISTIC = "p_value_adjusted"`;
    - `INTERACTION_SCAN_NOT_READ = ("f", "ss_interaction", "max_f", "difference", "null_max_share")`;
    - `INTERACTION_SCAN_TABLE = "interaction_scan.csv"`.
- **`INTERACTION_SCAN_POWER_STATEMENT`:** the cell type (`astrocyte`) and cohort (`cohort_30`) it applies to, and the plant (20 features, 10 per block, `Male`, stage weights 0 / 0.5 / 1, random signs). It also holds the power at e = 2.0 / 3.0 (0.692 / 1.000) with recall (0.059 / 0.525), the source (MOTCO `v0.20.0`, `docs/reports/interaction-scan-2026-11-09.md`, `results/interaction-scan-2026-11-09/plasmode/settings.json` with that file's sha256), and the five export sha256s that file records.
  - A task checks the constants against `git -C ../MOTCO show v0.20.0:<settings.json>`.
- **`check_motco_scan_alpha()`:** MOTCO's `INTERACTION_SCAN_ALPHA` must equal `ANALYSIS_ALPHA`. It runs at the first analysis, beside `check_motco_rank_floor` and `check_motco_split_parts`, so a MOTCO that moved its hit level stops the run.
- **Where it is passed:** `interaction_scan=ANALYSIS_INTERACTION_SCAN` in every construction, the latent arms included (no exception, so no import-time tie is needed).
- **`_interaction_scan_problem(record, table, params, *, n_jobs, n_donors, n_covariates)`** runs in `_analyze_cohort` after the split check and before any write.
  - When not requested: no block, no `params.interaction_scan`, no table.
  - When requested, the record must have `params.interaction_scan = true` and a block whose:
    - `blocks` keys equal `ANALYSIS_BLOCKS` in order, each `n_features` equal to the record's block width, and whose `n_features` is their sum;
    - `df_interaction == len(STAGES) − 1` and `df_residual == n_donors − 2·len(STAGES) − n_covariates`;
    - `p_value_global` is not null, and `hit_alpha == ANALYSIS_ALPHA`;
    - `permutations`, `seed` and `residual_correction` equal `params`;
    - `draws_match_rrpp == (n_jobs == 1)`.
  - The table must have:
    - `n_features` rows in the record's block order, with the columns above (three `difference[...]`);
    - its non-degenerate rows with `p_value_adjusted ≤ α` equal to the record's `hits` (same features, same order by F, equal F and adjusted p) and to the per-block hit counts;
    - its minimum adjusted p equal to `p_value_global`;
    - `(p_value_global ≤ α) == (len(hits) > 0)`.

  Any failure raises naming the cohort and the entry, with nothing written. A nonzero `n_degenerate` is a finding naming the features, not an error.
- **The manifest:** `measurement.interaction_scan = {"requested": true, "role": "confirmatory", "statistic": "p_value_global", "hit_statistic": "p_value_adjusted", "hit_alpha": 0.05, "not_read": [...], "table": "interaction_scan.csv"}`.
- **Alternative rejected:** reading the scan only from `interaction_scan.csv` after the write. The check must come before anything is written, as the correction and split checks do.

**D3. Summary, CLI equivalent and comparison columns.**
- **`render_analysis_summary`** adds an "Interaction scan" section after the split table, only when the record has a block:
  - the table: features, df, maximum F with feature and block, global p, status, hits per block, and null-maximum share per block "(descriptive, not read)";
  - then the hits (block, feature, F, adjusted p, differences at `0`, `1`, `2`), or the line "No hit (adjusted p ≤ 0.05)";
  - the draws line, as for the split test.

  A record without the block renders byte-identically to today.
- **`_cli_equivalent`** appends `--interaction-scan` when `params.interaction_scan`.
- **`INTERACTION_SCAN_COLUMNS`** = `scan_max_f`, `scan_max_block`, `scan_p`, `scan_rejected`, `scan_hits_rna`, `scan_hits_atac`. `_scan_columns(record)` fills them after `SPAN_SPLIT_COLUMNS` in `_sensitivity_row`, empty without a block. `_render_scan_table` gives the Markdown section.
- **`_scan_dependence_findings(ref, record, label, reference, kind)`** writes a finding when the scan status differs. It is used by the preprocessing, latent (`extended_covariates` against the primary) and exclusion comparisons. The rank sweep uses the D5 identity check instead.

**D4. The reference check.** `_check_reference` adds `_record_interaction_scan(record)`: `params.interaction_scan` true and a block. Otherwise it raises naming the file and `analyze`, like the floor, correction and split checks. All five readers inherit it.

**D5. Identity checks.**
- **Rank sweep:** after each forced run, `run_rank_sensitivity` compares the forced record's `interaction_scan` block with the reference's (dict equality) and the forced `interaction_scan.csv` bytes with the reference's. A mismatch raises naming the rank and the first differing key or line. The sweep writes under `<cohort>.tmp/`, so nothing is published.
- **Latent arms:** after both arms, `run_latent_sensitivity` compares `latent_adjusted` with `extended_covariates` the same way, before the comparison directory is written.
- **The Markdown** states the identity in one line.
- **Why a stop rather than a finding:** an inequality would mean the scan reads something MOTCO's contract says it does not. Reporting it to MOTCO comes before any reading.

**D6. Scan concordance cells in the stage DE.**
- **The rules:** `STAGE_DE_SCAN_CONCORDANCE_RULE` and `STAGE_DE_LIMMA_SCAN_CONCORDANCE_RULE` are dicts like the existing rules: model, coefficient or test, feature set `exported`, `scan_detected` text, `cell` text, and the overlap definition. They are tied at import to `STAGE_DE_CONCORDANCE_RULE` / `STAGE_DE_LIMMA_CONCORDANCE_RULE` (same model, test and feature set).
- **Summary columns:** `STAGE_DE_SUMMARY_COLUMNS` (and so `STAGE_DE_LIMMA_SUMMARY_COLUMNS`) gains `scan_hits`, `scan_hits_significant`, `scan_detected` and `scan_concordance` at the end. They are filled on the rule row only, as `concordance` is.
- **Method comparison:** `STAGE_DE_LIMMA_COMPARISON_COLUMNS` gains `scan_detected`, `scan_concordance_deseq2` and `scan_concordance_limma` at the end. A scan cell that differs between the methods is a method-dependence finding.
- **Where the hits come from:** `_stage_de_prepare` already loads the reference record. It returns the hits per omic (block `rna` / `atac`, the export's feature names). Every hit must be an exported feature of the omic's tables, or the step stops.
- **Manifests:** they record the hits read and both rules. `stage_de.md` and `stage_de_limma.md` gain the cells.
- **Fit records:** unchanged. The scan is read at summary time, so the fits' fingerprints do not include it.
- **Readers:** the positive control (E1), figure 10 and the limma step's DESeq2 binding read summary columns by name. A task verifies that they read the same values from the widened tables.
- **Alternative rejected:** the cell in the re-issue report only. The `delta` cells live in the steps, and a cell outside them would be computed by hand.

**D7. Figure 12 and the figures checks.**
- **Registry:** `FIGURES` gains `interaction_scan` (number 12, title "Per-feature interaction scan"), with reading keys `display_only`, `interaction_scan` and `sensitivity_dependence`.
- **Panels.** They are read through `FigureReader` (records via `result`, tables via its CSV reader, both logged):
  - (a) per primary, the `FIGURE_SCAN_TOP_N = 50` smallest adjusted p-values per block, ranked, on a log axis with α, hits filled and the block coloured;
  - (b) the global p of every registered analysis (`registered_analyses()`, 44), `cohort_30` filled and `cohort_50` hollow, with the forced ranks and `latent_adjusted` arms drawn as grey ticks labelled "= reference";
  - (c) per primary, the null-maximum share per block as bars beside the feature share, with no threshold.
- **Caption builder:**
  - per primary, the global p, maximum F with feature and block, hits per block and shares;
  - the hits' names when any;
  - each cohort's scan concordance cells from the stage-DE summaries;
  - the counts drawn (44, 28 distinct).
- **`FIGURE_READING_RULES["interaction_scan"]`:** the policy, the reading of a rejection and a non-rejection, the power statement's scope, the descriptive values, and the size statement.
- **Consistency checks:** `_measurement` gains `interaction_scan` and `interaction_scan_present`, and `check_consistency` requires both for every analysis drawn. The figure binds each table read to its record (rows, hits).
- **Unchanged:** figures 05 and 11, whose constructions are frozen.

**D8. Gates on `nixos-home`, before the move.**
1. *Pin-move gate.*
   - Run `analyze --cell-type {astrocyte,microglia} --cohort all` in a scratch worktree of main `5f70673` (`uv sync --locked`, MOTCO 0.18.0) into `<scratch>/v0180/<cell type>`.
   - From the change branch, run `_analyze_cohort` with the analyze parameters and `interaction_scan=False` into `<scratch>/v0200_none/<cell type>`.
   - Pass condition: the four CSVs are byte-identical per cohort, and `result.json` and the manifest are equal after dropping `runtime_seconds`, `created_at`, `git`, `software.motco`, `cli_equivalent` and `measurement.interaction_scan`. The `span_split` block is included in the comparison, which checks MOTCO's refactor of the split test.
2. *Kept-output comparison* (recorded, not a gate). Compare `v0180` and `v0200_none` with the kept primaries. Those were made on `nixos-home` under 0.18.0, so byte-identity is expected and any difference is a finding.
3. *Contract gate.*
   - Run the change branch's `analyze --cell-type … --cohort all` into `<scratch>/v0200_scan/<cell type>` and compare it with `v0200_none`.
   - Pass condition: the four CSVs are byte-identical, only `v0200_scan` has `interaction_scan.csv`, and `result.json` is equal after dropping `interaction_scan`, `params.interaction_scan` and `measurement.runtime_seconds`. That includes every p-value, null summary and the `span_split` block.
   - A failure stops the change and is reported to MOTCO.
   - `v0200_scan` is also the reproducibility reference for D10 step 1.
4. *Power-statement binding.* Compare the five sha256s of the astrocyte `cohort_30` export with `INTERACTION_SCAN_POWER_STATEMENT` and with the plasmode settings at `v0.20.0`, and record the result.

The gates run before the move, so a failed gate leaves `results/` untouched. They run after the predeclaration commit (task 1.1), because the contract gate is the first real-label scan.

**D9. Superseded layout.**
- **The move:** `python3 scripts/move_superseded.py --change mtg-interaction-scan` over the 15 directories of the 2026-11-06 re-issue.
  - Astrocyte: `motco_analysis`, `motco_rank_sensitivity`, `motco_diagnostics`, `motco_sensitivity`, `motco_latent_sensitivity`, `stage_de`, `stage_de_limma`, `motco_exclusion_sensitivity`.
  - Microglia, under `microglia/`: `motco_analysis`, `motco_diagnostics`, `motco_latent_sensitivity`, `stage_de`, `stage_de_limma`, `motco_exclusion_sensitivity`, `positive_control`.
  - It runs on `nixos-home`, then on the cluster. `results/figures/` is not moved.
- **Before the move:** record the sha256 of every MOTCO-independent file (both cell types' `pseudobulk/` and `motco_input/`, `data/annotation/`), every input hash of `reports/figures/motco-span-split-2026-11-06/figures_manifest.json`, and every `result.json` hash cited by `reports/motco-span-split-2026-11-06.md`.
- **After the move:** each must resolve under `results/superseded/mtg-interaction-scan/` or at its unchanged path.

**D10. Rerun order and commands.** These are the 2026-11-06 commands, through `nix develop -c uv run sead-sexdiff …` with `--n-jobs 1` and `--n-cpus 10` for `stage-de`:
1. `analyze`;
2. `rank-sensitivity` (astrocytes);
3. `axis-diagnostics`;
4. `preprocessing-sensitivity --variant all`, with `axis-diagnostics` per variant root;
5. `latent-sensitivity` with its diagnostics;
6. `stage-de --omic both` and `stage-de-limma --omic both`;
7. `exclusion-sensitivity`, then on its variant root `axis-diagnostics`, `stage-de --omic both` and `stage-de-limma --omic rna`;
8. `positive-control` on both bases;
9. after every other step, `figures --report motco-interaction-scan-<date>`.

**D11. Comparison script (superseded → new).** It runs over the 44 analysis pairs. For each pair:
- the D8.3 contract comparison against the superseded record: CSVs, and records outside `interaction_scan`, `params.interaction_scan`, `software` and runtimes. Byte-identity is expected, because both were made on `nixos-home`;
- the new scan block's values, status and hits;
- the dependence findings;
- the two identities (12 forced ranks, 4 arm pairs).

It also checks:
- every p-free output (axis diagnostics, adoption and precedence, mediator, flag and axis checks, E1–E3 on both bases);
- the stage-DE and limma tables (hash-equal expected) and the summaries on their pre-existing columns (equal), with the new scan columns listed.

It produces the report's tables. Outside the D8 gates every difference is a finding, never a gate.

**D12. Version and measurement consistency check.**
- **Manifests:** every `*manifest*.json` outside `results/superseded/` that records a MOTCO version records `0.20.0`, `results/figures/figures_manifest.json` included.
- **Records:** every `result.json` carries:
  - `df`, and floor 2 or a forced rank;
  - a `span_split` block, except the four `latent_adjusted` arms, whose manifests record `not_available`;
  - an `interaction_scan` block, with `interaction_scan.csv` in its manifest's outputs (all 44).
- **When:** it runs before the report and again before the commit.

**D13. Report.** `reports/motco-interaction-scan-<date>.md` is assembled by a scratch script from generated files. It contains:
- **Provenance:** the machine, the predeclaration commit, the pin and wheel check, the gates, the kept-output comparison and the move.
- **MOTCO's validation** with its numbers (S1–S3, G, the plasmode table, the block balance, the localization rows that bound the reading).
- **The predeclared policy and reading,** quoted from the proposal.
- **The scan results:** the four primaries in full (global p, maximum F, hits with their per-stage differences, null-maximum shares), then per family one table with the scan status against the primary's.
- **The concordance cells:** per cohort × omic, both methods, with the overlap.
- **The power-statement binding** (D8.4) and its wording in the outcome that occurred.
- **One section per re-issued step,** with the D11 outcome table.
- **The reading,** the findings and the limitations.
- **The frozen set:** the twelve SVGs, with captions pasted verbatim from the frozen `figures.md`.

A check script verifies that:
- every scan value (global p, maximum F, hits, shares) and every statistic, p-value, rank, CV AUROC, κ and split value of every `result.json` is in the text;
- every concordance cell is in the text;
- the captions match `figures.md`, and the frozen SVGs and inputs match the frozen manifest.

**D14. Machine.** The re-issue runs on `nixos-home`, which made the superseded outputs and holds the annotation. The gates compare runs made on the same machine, and every cross-machine comparison is recorded, not gated.

## Risks / Trade-offs

- **[A version-only move re-runs everything]** → Unavoidable: `_check_reference` compares MOTCO versions, and `figures` requires one version throughout. The compute is that of 2026-11-06 (a few hours). The scan adds under a second per analysis.
- **[Forty-four analyses × five tests invite selective reading]** → Only the two `cohort_30` primaries' global p-values are scan results. Every other status is dependence, 16 of the 44 are identical by construction, and the descriptive values carry no status.
- **[A hit is over-read as a driver or a program]** → The reading allows only "this feature's adjusted sex difference changes across stages", names correlated neighbours, and forbids trajectory, attribution and program claims.
- **[A no-hit astrocyte scan is read as "no sparse difference"]** → The power statement is scoped to the plasmode's plant on this export and lists what it does not exclude (smaller effects, thinly spread effects, other shapes). It is bound to the export hashes.
- **[The threshold leans on ATAC]** → The null-maximum shares are printed in `summary.md`, drawn in figure 12 and carried as a limitation. Per-block maxima are MOTCO's follow-up, not done here.
- **[The scan and the stage DE disagree]** → Each disagreement is a predeclared cell and a finding, and the expected direction (stage DE detects, scan not) is stated in advance. Agreement is not replication.
- **[An identity check fails]** → The step stops before publishing, and the failure is reported to MOTCO as a contract question before any reading.
- **[MOTCO's refactor moved the split block]** → The D8.1 gate compares the split block too, and stops the change before anything is moved.
- **[Widened stage-DE summaries break a reader]** → Columns are appended and read by name. A task checks E1, figure 10 and the limma binding on the widened tables, and D11 checks the pre-existing columns for equality.
- **[Degenerate features]** → They are recorded as a finding naming the features. They get no F and cannot be hits.
- **[`draws_match_rrpp` false]** → It can only happen at `n_jobs` ≠ 1, and the record check requires it true at `n_jobs = 1`, which every reported run uses.
- **[Written after the stage-DE counts were known]** → Stated in the proposal (all zero). The scan's real-label values are unknown until the contract gate, which runs after the predeclaration commit.

## Migration Plan

1. Branch `mtg-interaction-scan`: commit the planning artifacts as the predeclaration (task 1.1), then the pin, code and docs (D1–D7), verified on a synthetic fixture.
2. Run the gates on `nixos-home` (D8). Stop on failure.
3. Move the superseded outputs (D9), then rerun everything (D10).
4. Run the comparisons and the consistency check (D11, D12), then write the report and its check (D13).
5. On the cluster: pull, `uv sync --locked`, the same move, then rsync the re-issued directories and verify their sha256.
6. PR; after the merge, archive and sync the specs.

**Rollback.** Before step 3 nothing under `results/` has changed. After it, `MOVED.json` lists every move, and moving the directories back with their hashes checked restores the 0.18.0 state on the main commit.

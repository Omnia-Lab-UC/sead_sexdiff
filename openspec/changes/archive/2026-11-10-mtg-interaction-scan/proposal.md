# Proposal

## Why

None of the case study's tests can see a sex × stage difference confined to a few features. `delta`, `angle` and `shape` see only the part of a difference inside the fitted stage-PLS span. The split test's parts sum interaction energy over thousands of coordinates. MOTCO's plasmode on the astrocyte `cohort_30` export planted a stage ramp to 3 within-stage SDs on 20 of the 10,000 features in one sex. All six tests rejected it at the null rate (0.016–0.048), so every non-rejection reported so far is blind to sparse differences.

MOTCO 0.20.0 (`interaction-feature-scan`, 2026-11-09, decision **`scan_validated`**) adds a per-feature interaction scan (`TrajectoryAnalysisParams.interaction_scan`).
- **What it computes:** each analysed feature's group × stage interaction F, with the covariate-adjusted trajectory models, under RRPP's Freedman–Lane draws at the analysis's seed and with its residual correction. A global p-value comes from the largest F. Each feature gets a Westfall–Young step-down max-T adjusted p-value.
- **A hit** is a feature with adjusted p ≤ 0.05.
- **What the record gains:** an `interaction_scan` block, plus the table `interaction_scan.csv`, which gives each feature's covariate-adjusted group difference at each stage. MOTCO states that every other output is unchanged. The scan uses neither the latent space nor the bootstrap, so it is available with the latent adjustment (unlike the split test).

Its validation found:
- **Size holds** on three simulated designs (0.048–0.058, against 0.0695). It also holds on real noise at the case study's configuration (six covariates, `df`, 10,000 features): 0.042. That real-noise check ran on the astrocyte `cohort_30` export, whose five files hash exactly as ours do (checked 2026-11-10 against `results/interaction-scan-2026-11-09/plasmode/settings.json`).
- **Familywise control holds under the plasmode's partial null.** A feature outside the planted set is a hit in 0.042 of replicates, at 2 and at 3 SDs.
- **It sees the sparse plant.** That is 20 features (10 per block) in `Male`, with stage weights 0 / 0.5 / 1 and random signs. It is rejected in 0.692 of replicates at 2 SDs and in 1.000 at 3 SDs, where 0.995 of the hits are planted and recall is 0.525.
- **Limits:**
  - It is weak where a difference spreads thinly over many features: magnitude power is 0.01 on the simulated SEA-AD table, against `delta`'s 0.21.
  - ATAC gives 0.618 of the null maxima with half the features.

MOTCO's report (§8) assigns this repository the follow-up: apply the scan to both cell types, pin 0.20.0, and fix the multiplicity policy across the scan, the split test and the trajectory statistics before any real-label analysis. 0.20.0 is tagged `v0.20.0` (commit `5f6b07d`) and published on PyPI (wheel sha256 `d4248a8e19fd65f8…`). On 2026-11-10 the downloaded wheel's `motco/` tree was found equal to the tag's `src/motco`. No real-label scan has been run on these exports.

## What Changes

- **Pin `motco==0.20.0` from PyPI** with the existing rule:
  - the wheel's `motco/` tree must be byte-identical to `src/motco/` at `v0.20.0`;
  - `uv export` may change only the motco line (`pyproject.toml`'s dependency pins are unchanged between the tags);
  - `ANALYSIS_MOTCO_MIN_VERSION` becomes `0.20.0`.

  Between `v0.18.0` and `v0.20.0`:
  - the analysis path gained the opt-in scan;
  - the split test's model-space construction and input checks were factored into helpers it shares with the scan (`_model_spaces`, `_check_permutations`);
  - 0.19.0 added a harness-only setting (`SimulationEvaluationParams.latent_coordinates`, verdict `power_loss`, not in `analyze_trajectories`).

  The pin-move gate checks that nothing else moved, including the split block.
- **The scan in every analysis.** A new constant, `ANALYSIS_INTERACTION_SCAN = True`.
  - **Where it is passed:** explicitly in every `TrajectoryAnalysisParams` the package builds. That covers the primaries, the forced ranks, the preprocessing variants, both latent arms (`latent_adjusted` included) and the exclusion variants: 44 analyses.
  - **Checked after each run, before anything is written.** Otherwise nothing is written. The checks:
    - the block is present exactly when requested;
    - its feature counts per block match the record's blocks;
    - its interaction and residual degrees of freedom match the design;
    - its permutations, seed and `df` correction are the analysis's;
    - `draws_match_rrpp` holds at `n_jobs = 1`;
    - `interaction_scan.csv` has one row per feature;
    - a global rejection occurs exactly when there is at least one hit.
  - **Recorded:**
    - in the manifest's `measurement.interaction_scan`;
    - as `--interaction-scan` in `cli_equivalent`;
    - in `summary.md`, as a scan table (global p, maximum F with its feature and block, hits per block, and the null-maximum share per block marked descriptive) followed by the hits with their per-stage differences.
  - **Checked by readers:** `_check_reference` requires the block on the primary.
  - **Identical by construction where only the latent space differs.** Every forced rank's scan must equal its reference's, and `latent_adjusted`'s must equal `extended_covariates`'s (block and table). These are checks that stop the step, never sensitivity readings.
  - **In the comparison tables:** the shared row builder gives every table (`rank_sensitivity.csv`, `sensitivity.csv`, `latent_sensitivity.csv`, `exclusion_sensitivity.csv`) the scan columns: maximum F and its block, global p, status and hits per block.
- **Concordance cells against the stage DE.**
  - **Where they go:** `stage-de` and `stage-de-limma` gain, per cohort × omic, a cell comparing `scan_detected` (the analysis read has at least one hit in that omic's block) with the existing `interaction_detected`. That is DESeq2's exported `factor_interaction` union count > 0, and limma's exported `interaction_joint` count > 0.
  - **Overlap:** the hit overlap (hits also significant under each method) is counted beside each cell.
  - **Discordance:** in either direction it is a finding, never a result.
  - **What stays the same:** the existing `delta` cells and `methods_agree`.
- **Figure 12 `12_interaction_scan`** in `FIGURES`:
  - per primary, the smallest adjusted p-values per block on a log axis with α, hits marked;
  - the global scan p-value of every registered analysis, with the forced ranks and `latent_adjusted` drawn as identical to their references;
  - per primary, the share of null maxima per block against its share of features.

  The consistency pre-pass requires the block in every analysis drawn. A new reading-rule text restates the policy below.
- **Re-issue under `motco-release-migration`.** Every reader compares MOTCO versions, and the figures need one version throughout, so the pin move re-issues every MOTCO-dependent output of both cell types:
  - the same-machine pin-move gate (0.20.0 without the scan against 0.18.0, split test on in both);
  - the contract gate (0.20.0 with and without the scan: everything but the block, its parameter and `interaction_scan.csv` equal);
  - the move to `results/superseded/mtg-interaction-scan/`;
  - every step rerun in dependency order, then the version and measurement consistency check;
  - `figures --report` last.
- **Predeclared policy and reading**, written here before any real-label scan (the user's decisions of 2026-11-10):
  - **Multiplicity.**
    - The scan's **global p-value** (`p_value_global`) is confirmatory. It is read at α = 0.05 under `ANALYSIS_REJECTION_RULE` (p ≤ α), unadjusted.
    - Under step-down max-T a global rejection is equivalent to at least one hit. The hits are the features with `p_value_adjusted ≤ 0.05`, so the familywise error rate over the scanned features is controlled within the scan.
    - No adjustment is made across `delta`, `angle`, `shape`, the complement and the scan, or across cohorts and cell types. The report states that each primary carries five tests read at α.
    - Reported and not read: the per-feature F and interaction SS, the maximum F, the per-stage differences (they show how a hit differs; none is tested per stage) and the null-maximum shares.
  - **Roles.** Each cell type's `cohort_30` primary gives the scan result, and `cohort_50` is sensitivity. In the preprocessing variants, the `extended_covariates` arm and the exclusion variants, a scan status that differs from the primary's is preprocessing, covariate-set or cohort dependence, not a result of that analysis. The forced ranks and `latent_adjusted` cannot differ (checked).
  - **A rejection** names features whose covariate-adjusted sex difference changes across stages, with the familywise error rate over the 10,000 scanned features at 0.05.
    - It does not say that a feature drives a latent trajectory difference. Correlated neighbours of a differing feature can be hits too.
    - It carries no geometry or trajectory claim and starts no attribution.
    - It changes no reading of `delta`, `angle`, `shape` or the complement.
    - ATAC hits are named by peak coordinates only; no gene is assigned here.
  - **A non-rejection.**
    - **Astrocyte `cohort_30` primary:** it reads "no sparse ramp of the plasmode's form detected". At this export's design such a ramp (20 features, one sex, stage weights 0 / 0.5 / 1) is detected in 1.000 of replicates at 3 within-stage SDs and in 0.692 at 2 SDs.
      - It holds only while the export's five files hash as the plasmode's (checked). Otherwise it falls back to uninformative, with a finding.
      - It excludes nothing smaller, nothing spread thinly over many features (where the scan is weak) and no other shape.
    - **Every other analysis:** uninformative (no power measured there).
  - **Concordance cells.** A disagreement between the scan and a stage-DE method is a finding, never a result, and agreement is not replication (same donors). The tests differ:
    - in features and scale: z-scored log2-CPM of the exported features against raw counts;
    - in covariates: both nuclei counts against the omic's own;
    - in the null: permutation against parametric;
    - in the error criterion: familywise against BH FDR.

    So "stage DE detected, scan not" is the expected direction of a disagreement at equal power, and the other direction is possible too.
  - **Size** was measured on real noise only for the astrocyte `cohort_30` export at six covariates: global 0.042, familywise 0.042 under the partial null. For `cohort_50`, microglia, the eight-covariate arms, the exclusion variants and the preprocessing variants it is inferred. Each is carried as a limitation.
  - **Known when this change was written.**
    - No analysis rejects `delta`, `angle` or `shape` in any of the 44 analyses, and no complement is rejected in the 40 that carry the split test.
    - Every exported stage-DE interaction count of the primaries is 0: both methods, both cohorts, both cell types, both omics.

    The scan only adds a block and a table, so every other value is expected to reproduce. The contract gate checks the primaries and the comparison script checks every other analysis. A failed reproduction is a finding and changes no reading.

    Nothing else is re-read:
    - the stage-DE `delta` cells;
    - the positive control (no p-value);
    - the adoption, mediator, flag and axis checks;
    - the `latent_adjusted` reading;
    - the split-test policy, apart from "four tests" becoming "five tests".
- **One committed report**, `reports/motco-interaction-scan-<date>.md`, covering both cell types. It holds:
  - provenance (pin, gates, move);
  - MOTCO's validation with its numbers;
  - the scan results of every analysis, with the primaries in full (hits listed with their per-stage differences);
  - the concordance cells;
  - per re-issued step, its outcome table (superseded → new);
  - the reading in the outcome that occurred;
  - findings and limitations;
  - the frozen twelve-figure set.

  Numbers are pasted from generated files and checked by a script.
- **README and AGENTS.md:** the pin, the constant and the policy, the concordance cells, figure 12 and the superseded entry.

Out of scope:
- MOTCO's possible follow-ups (per-block maxima, a moderated per-feature statistic);
- attribution, and gene annotation of ATAC hits;
- MOTCO's unmerged block weighting and sparse stage-PLS (no verdict yet), and its orthonormal coordinates (`power_loss`);
- any change to a registry, expectation, adoption rule, threshold, the stage-DE models or the reading of `delta`, `angle`, `shape` or the complement;
- rerunning the MOTCO-independent steps;
- edits to earlier reports.

## Capabilities

### New Capabilities
- `motco-interaction-scan`: the per-feature interaction scan carried by every analysis. This covers the registered constant, the record and table checks, the manifest, summary and table columns, the reference check, the identity checks across ranks and latent arms, the predeclared policy (global p confirmatory, hits familywise, descriptive values not read), the power statement for the astrocyte `cohort_30` primary, and the reading and limitations the report carries.

### Modified Capabilities
- `motco-trajectory-analysis`: MOTCO ≥ 0.20.0. The scan joins the fixed measurement settings, `interaction_scan.csv` joins the outputs, and `summary.md` gains the scan table.
- `motco-span-split-test`: each primary carries five tests read at α (was four), with no adjustment across them.
- `motco-rank-sensitivity`: the cross-rank table gains the scan columns, and every forced run's scan must equal the reference's.
- `motco-preprocessing-sensitivity`: the comparison table gains the scan columns.
- `motco-latent-adjustment-sensitivity`: both arms carry the scan, and `latent_adjusted`'s must equal `extended_covariates`'s. The comparison table gains the scan columns.
- `donor-exclusion-sensitivity`: the comparison table gains the scan columns.
- `pseudobulk-stage-de`: the scan concordance cell beside the `delta` cell, with the hit overlap.
- `pseudobulk-stage-de-limma`: the scan concordance cell beside the `delta` cell, with the hit overlap, in the summary and the method comparison.
- `motco-release-migration`: the contract gate covers a setting that adds a named record block and a named output file.
- `case-study-figures`: twelve figures, with `12_interaction_scan` added. The consistency pre-pass checks the scan setting, and the binding context includes the scan policy.

## Impact

- `pyproject.toml`, `uv.lock`: `motco==0.20.0` from PyPI.
- `src/sead_sexdiff/data.py`:
  - `ANALYSIS_MOTCO_MIN_VERSION`, `ANALYSIS_INTERACTION_SCAN` and the `INTERACTION_SCAN_*` policy constants;
  - `interaction_scan` in the five `TrajectoryAnalysisParams` constructions;
  - the record and table checks and `measurement.interaction_scan` in `_analyze_cohort`;
  - `_cli_equivalent`, `render_analysis_summary`, `_check_reference`, `_sensitivity_row`, the rank sweep's and the latent step's identity checks;
  - the stage-DE and limma summaries, comparison and manifests (new columns and rules).

  No CLI flag is added.
- `src/sead_sexdiff/figures.py`: figure 12, its caption, `_measurement`, the consistency check and a reading rule. The module's sha256 changes, so every figure is redrawn (it is anyway in the re-issue).
- `README.md`, `AGENTS.md`.
- `results/` (gitignored), on `nixos-home` and the cluster:
  - the superseded outputs are moved, and every MOTCO-dependent output is re-made;
  - the compute is that of the 2026-11-06 re-issue (a few hours), plus under a second per analysis for the scan (0.72 s at 80 × 10,000, MOTCO's timing).
- `reports/`: one new report with its frozen figure set.
- Scientific choices introduced here (flag and ask before changing any of them):
  - the scan in every analysis;
  - its global p-value confirmatory at α = 0.05, unadjusted, with hits at familywise adjusted p ≤ 0.05;
  - five tests per primary with no adjustment across them;
  - the power statement for the astrocyte `cohort_30` primary only;
  - the concordance cells against both stage-DE methods;
  - the descriptive values not read;
  - the reading and limitations above;
  - re-issuing every MOTCO-dependent output under 0.20.0.

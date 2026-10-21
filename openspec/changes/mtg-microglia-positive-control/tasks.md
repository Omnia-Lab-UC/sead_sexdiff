# Tasks

## 1. Spike and registry

- [ ] 1.1 Confirm that `mtg-astrocyte-stage-de` is applied (`stage-de` exists and `results/stage_de/<cohort>/stage_de_summary.csv` exists for both astrocyte cohorts). Verify with `uv run sead-sexdiff stage-de --help` and `ls`.
- [ ] 1.2 Write and run the obs-only spike script on the cluster login node over all 84 RNA and 86 ATAC objects: per release, the `Subclass` × `Supertype` labels and counts of the `Immune` / `Microglia-PVM` subclass; per donor and omic, the microglia / PVM nuclei count; the sex × stage cells at ≥ 30 and ≥ 50 in both omics. Paste the tables into `openspec/changes/mtg-microglia-positive-control/spike.md` with the label decision (design Open Question). Verify that every `Immune` / `Microglia-PVM` supertype of both releases is listed as microglia / PVM or as non-population, and that `spike.md` states the cohort cells.
- [ ] 1.3 If the spike leaves a sex × stage cell empty at ≥ 30, amend the proposal with the microglia thresholds decided from those counts before anything else runs; otherwise record in `spike.md` that the astrocyte thresholds are kept. Verify that the proposal and `spike.md` agree.
- [ ] 1.4 Add `CELL_TYPES` (astrocyte: today's mapping, `results/`, today's preregistered cells; microglia: the spike's labels, `results/microglia/`, the spike's cells), `POSITIVE_CONTROL_EXPECTATIONS`, `POSITIVE_CONTROL_PRIMARY_COHORT` and the output names; restructure `TAXONOMY` / `taxonomy_for` / `validate_taxonomy` per design D1 and make `PREREGISTERED_COHORTS` per cell type. Verify by import that `taxonomy_for(path)` for both releases returns the astrocyte mapping it returns today.

## 2. One selection path, one flag

- [ ] 2.1 Rebuild the synthetic dataset in the scratchpad and run `export`, `analyze`, `axis-diagnostics` and `stage-de` on it with the unmodified package as the byte-identity baseline (design D10).
- [ ] 2.2 Implement `select_nuclei(..., cell_type=...)` with `select_astrocyte_nuclei` as a wrapper; thread `cell_type` through `pseudobulk_donor`, `technical_donor`, the per-donor caches (recorded, mismatch is an error), `run_pseudobulk`, `run_covariates`, `run_technical_covariates`. Verify on the two spike objects (any machine holding them) that `cell_type="astrocyte"` returns the same rows as before and `cell_type="microglia"` returns only registered microglia / PVM supertypes.
- [ ] 2.3 Add `cell_type_dirs`, default the directory keywords of `run_export`, `run_analysis`, `run_axis_diagnostics` and `run_stage_de` to the cell type's root, record `cell_type` in every manifest and check it at every read (design D4; an absent key in an existing astrocyte manifest reads as `astrocyte` with a finding). Verify that `export --cell-type microglia --pseudobulk-dir results/pseudobulk` (synthetic astrocyte pseudobulk) stops naming both cell types and writes nothing.
- [ ] 2.4 Add `--cell-type` to the seven MTG CLI commands with directory flags defaulting to `None`, update `scripts/pseudobulk.sbatch`, `scripts/technical.sbatch` and `scripts/sync_pseudobulk.sh` (design D9). Verify with `--help` on each command that the astrocyte defaults are shown.
- [ ] 2.5 Byte-identity: rerun 2.1 with the refactored package. Verify that every CSV and `summary.md` / `diagnostics.md` / `stage_de.md` is byte-identical and every JSON equal apart from runtime, timestamps, git and the added `cell_type` key.
- [ ] 2.6 Implement `run_positive_control`, `render_positive_control` and the `positive-control` CLI subcommand (design D7). Verify on the synthetic data (one root copied as the "control" with its manifests relabelled) that E1–E3 are computed from the right cells, that the overall outcome follows the all-three rule, that a missing summary stops the step naming `stage-de`, and that a reference root whose manifests are not `astrocyte` is refused.

## 3. Cluster runs

- [ ] 3.1 On the cluster, after a two-donor astrocyte dry run (`pseudobulk --limit 2` into the scratchpad) reproduces the cached per-donor counts and compositions of those donors: run `pseudobulk --cell-type microglia --omic both` through `scripts/pseudobulk.sbatch`, then `covariates --cell-type microglia --omic both` on the login node, then `technical-covariates --cell-type microglia --omic both` through `scripts/technical.sbatch`. Verify that `results/microglia/pseudobulk/` holds the six pseudobulk files, the covariate files and the technical files, that the manifests record `cell_type = microglia`, that the cohort summary equals the spike's cells (deviations as findings), and that every file under `results/pseudobulk/` is hash-unchanged.
- [ ] 3.2 Sync `results/microglia/pseudobulk/` to the workstation with `scripts/sync_pseudobulk.sh` and verify the hashes against the cluster.

## 4. Workstation runs

- [ ] 4.1 Run `export --cell-type microglia --cohort all`, `analyze --cell-type microglia --cohort all` (`n_jobs = 1`), `axis-diagnostics --cell-type microglia --cohort all` and `stage-de --cell-type microglia --cohort all --omic both`. Verify that each manifest records `cell_type = microglia` and the upstream hashes, and that every file under `results/motco_input/`, `results/motco_analysis/`, `results/motco_diagnostics/` and `results/stage_de/` is hash-unchanged.
- [ ] 4.2 Run `positive-control --cohort all`. Verify that `results/microglia/positive_control/<cohort>/` holds the three outputs, that every compared value equals its source file, and that the overall outcome follows the all-three rule in `cohort_30`.
- [ ] 4.3 Rerun `analyze --cell-type microglia --cohort cohort_30 --force` into a scratch directory. Verify that `latent_scores.csv`, `pls_loadings.csv`, `stage_means.csv` and `analysis_metadata.csv` are byte-identical and that the statistics, p-values and rank are equal.
- [ ] 4.4 Copy `results/microglia/` to the cluster and verify the hashes.

## 5. Report and docs

- [ ] 5.1 Write `reports/motco-microglia-positive-control-<date>.md` from the generated files: provenance, the spike, cohort tables, the primary analysis summary per cohort, the stage-DE counts of both cell types side by side, the E1–E3 table, the first-axis correlations, the reading rules and limitations the spec requires. Verify with a check script that every value in the E1–E3 table equals `positive_control.csv`, every statistic, p-value and rank equals the microglia `result.json` records, and every DE count equals the two `stage_de_summary.csv` files.
- [ ] 5.2 Update README (registry, `--cell-type`, `positive-control`, outputs, reading rules) and AGENTS.md (the registry and expectations under "flag and ask", the flag on every MTG command, pipeline step 12, the results-root layout, the cluster scripts' `CELL_TYPE`). Verify that every `CELL_TYPES` / `POSITIVE_CONTROL_*` name in both files exists in `data.py`.
- [ ] 5.3 Run `openspec validate mtg-microglia-positive-control` and verify that it passes.

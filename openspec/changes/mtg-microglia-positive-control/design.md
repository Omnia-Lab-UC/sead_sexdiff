# Design

## Context

See proposal.md — Why. Today the astrocyte population is hard-wired in three places: `TAXONOMY[release]["astrocyte_subclass"]` / `["supertypes"]` (with `_ASTRO_SUPERTYPES`), `select_astrocyte_nuclei` (the one nucleus-selection path, used by `pseudobulk_donor` and `technical_donor`), and `PREREGISTERED_COHORTS`. Directories are plain defaults (`results/pseudobulk`, `results/motco_input`, …) on every `run_*` function and CLI flag; every downstream step already takes its input directories as arguments, and the hash chain (pseudobulk manifest → export manifest → analysis manifest → sensitivity / diagnostic manifests) is what binds a run to its inputs. The prerequisite change adds `stage-de` with the same directory arguments. SEA-AD's RNA `Immune` subclass holds microglia / PVM supertypes beside lymphocyte and monocyte supertypes; the ATAC release names the subclass `Microglia-PVM`. Compute nodes have no internet; the RNA objects are read on the cluster.

## Goals / Non-Goals

**Goals:**
- One registry, one selection path, one flag; astrocyte outputs byte-identical.
- The microglia chain reproducible from the registry entry and the flag alone.
- The expectations evaluated by code from files, never by hand.

**Non-Goals:**
- Sensitivity analyses on microglia; cross-cell-type MOTCO comparisons; other cell types.
- Any change to export, analysis or DE settings.

## Decisions

- **D1 — Registry shape.** `CELL_TYPES = {name: {"description", "results_root", "releases": {release: {"subclass": [...], "supertypes": {raw: canonical}, "other_supertypes": [...]}}, "preregistered_cohorts": {...}}}`. `TAXONOMY` keeps its per-release `other_subclasses` (release-wide subclass inventory, unchanged) and gains the population part from the registry at lookup (`taxonomy_for(path, cell_type)`). `other_supertypes` lists the subclass's non-population supertypes (lymphocyte, monocyte in RNA `Immune`), so `validate_taxonomy` still fails on an unknown label while excluding known non-microglia ones. For `astrocyte`, `other_supertypes` is empty and the resulting mapping equals today's. Alternative: a second taxonomy dict per cell type. Rejected: the release-wide subclass inventory is the same for every population and would be duplicated.
- **D2 — `select_nuclei(f, path, omic, *, cell_type="astrocyte", ...)`** replaces the body of `select_astrocyte_nuclei`, which stays as a one-line wrapper (AGENTS.md, README and the technical spec name it). Selection: `Subclass ∈ subclass labels`, then `Supertype ∈ supertypes` (the astrocyte case keeps every astrocyte supertype, as now). The returned dict gains `cell_type`. `pseudobulk_donor` and `technical_donor` pass it through; the per-donor cache records it and a cache entry of another cell type is an error.
- **D3 — Results root from the registry.** A helper `cell_type_dirs(cell_type)` returns the default directories; every `run_*` signature keeps its directory keywords but defaults them to `None` → resolved from `cell_type`. The CLI mirrors this: `--cell-type` first, directory flags default to `None`, help text shows the astrocyte defaults. The astrocyte root is `results/` (not `results/astrocyte/`) so no existing path moves and no report hash breaks. Alternative: move astrocyte outputs to `results/astrocyte/`. Rejected: every committed report resolves paths and hashes under `results/`.
- **D4 — `cell_type` in every manifest and checked at every read.** Added beside the existing hash checks: `check_export_inputs` (pseudobulk manifest), `_check_reference` / `load_analysis_inputs` (export manifest), `_check_analysis_inputs` (diagnostics), `run_stage_de`. An astrocyte manifest written before this change has no `cell_type` key; the readers treat an absent key as `astrocyte` and record that in a finding, so the existing astrocyte outputs are not rerun.
- **D5 — Preregistered cells per cell type**, `CELL_TYPES[name]["preregistered_cohorts"]`; `compare_preregistered(info, cell_type)`. The microglia cells are filled from the spike (task 1.2) by editing the registry before task 3.x runs. The thresholds stay `COHORT_THRESHOLDS`. If the spike shows an empty cell at ≥ 30 the proposal is amended before any run (metadata-based decision, documented in `spike.md`).
- **D6 — Spike as a script** (`scratchpad`, obs-only through `_read_obs_columns`, seconds per object) run on the cluster login node over all 170 objects; its table is pasted into `spike.md`. No pseudobulk, no matrix read.
- **D7 — `positive-control` step.** `run_positive_control(cohort, control_root, reference_root, out_dir)` reads: both `stage_de/<cohort>/stage_de_summary.csv` (E1: `factor_additive`, union, `exported`, `rna`), both `motco_analysis/<cohort>/result.json` (E2: `cv_mean_auroc`), the control's `motco_diagnostics/<cohort>/axis_correlations.csv` (E3: `pls_0`, raw Spearman; `stage` row against every row whose source is the technical table or the export's TMM / library-size variables; the six analysis covariates and `male` are not technical variables for E3 and are excluded by source). It checks each file's manifest chain (cell type, hashes) and writes the three outputs. The constants: `POSITIVE_CONTROL_EXPECTATIONS` (the three with their sources and comparison), `POSITIVE_CONTROL_PRIMARY_COHORT = "cohort_30"`.
- **D8 — Stage DE on microglia** is `stage-de --cell-type microglia`; no DE code changes beyond the flag and the manifest check.
- **D9 — Scripts.** `pseudobulk.sbatch` and `technical.sbatch` take `CELL_TYPE` (default `astrocyte`) and pass `--cell-type`; `sync_pseudobulk.sh` takes the cell type and syncs `<root>/pseudobulk/`.
- **D10 — Byte-identity of the astrocyte path.** The synthetic dataset builder of the previous change (rebuilt in the scratchpad) runs `export`, `analyze`, `axis-diagnostics` and `stage-de` before and after the refactor; CSVs and `summary.md` byte-identical, JSON equal apart from runtime, timestamps, git and the new `cell_type` key. For `pseudobulk` and `technical-covariates`, a two-donor `--limit 2` dry run on the cluster after the refactor is compared with the existing per-donor cache entries (counts, composition, fractions equal).

## Risks / Trade-offs

- [Microglia ATAC nuclei per donor may be fewer than astrocytes', emptying a sex × stage cell at ≥ 30] → the spike runs first; a decision on the microglia thresholds is made from metadata and written into the proposal before any pseudobulk (D5).
- [Microglia supertype labels differ between releases or include SEA-AD-specific states] → the spike lists every `Immune` / `Microglia-PVM` supertype per release; the registry maps them explicitly; an unlisted label is an error, as for astrocytes.
- [E1 depends on a prerequisite change] → applied first; the control step errors on a missing summary.
- [E3 as "stage beats every technical variable" is strict: stage |ρ| along a stage-supervised axis is expected to be high, but a strongly technical cell type could fail it even with real stage signal] → it is one of three predeclared expectations and the report states which failed and by how much; the threshold is not tuned after the fact.
- [A microglia `delta` rejection could be read as the result of the project] → the spec binds the reading: a positive control for the stage signal, and the trajectory result read under the existing rules with the Phase 6 limitation.
- [Threading `cell_type` through seven functions and the CLI touches many signatures] → defaults reproduce today's behaviour; the byte-identity check (D10) guards it.

## Open Questions

- Whether the ATAC `Microglia-PVM` subclass holds a PVM supertype separate from microglia and whether to keep it (the RNA side has `Micro-PVM_*` labels that bundle both). Resolved by the spike before the registry entry is filled; the default is to keep every supertype of the registered subclass that is microglia or PVM, matching the RNA bundling, and to record the choice in `spike.md`.

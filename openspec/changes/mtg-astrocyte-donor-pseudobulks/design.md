# Design

## Context

See `proposal.md` for motivation. Current state that shapes the approach:

- All logic lives in `src/sead_sexdiff/data.py`; `cli.py` is a thin argparse layer. The existing per-donor loop (`subset_and_concat_folder`) concatenates nucleus-level data across donors and pseudobulks per supertype with `decoupler`, which is the wrong grain for a donor-level multi-omic matrix.
- Bucket layout verified on 2026-10-16: `MTG/RNAseq/donors_objects/` has 84 `.h5ad` (32.8 GB, max 0.74 GB, release 2026-06-22); `MTG/ATACseq/donor_objects/` has 86 `.h5ad` (17.9 GB, max 0.62 GB, release 2024-12-06). The ATAC prefix also holds a 1.1 GB all-nuclei metadata CSV.
- Taxonomy differs between releases (RNA `Immune` vs ATAC `Microglia-PVM`). The astrocyte label strings and the ATAC matrix location (`X` vs a layer) are not yet observed; a two-file spike settles them before implementation of the mapping.
- Compute: the `ing` cluster (`cluster.ing.uc.cl`, SLURM, partition `512x1024`, `uv` installed in `/home1/tgonza/.local/bin`, Python 3.11 managed by uv, 1.1 TB free on `/home1`, no scratch, no quota enforced, no clone of this repo yet). Three workstations share the analysis via synced derived outputs. The MOTCO studies on the same cluster used `/home1/tgonza/MOTCO` and `scripts/motco_study_array.sbatch`.
- MOTCO (`../MOTCO/docs/api/analysis.md`) does no normalization, filtering or imputation, so this step must leave counts raw and leave every scientific choice visible in a manifest.

## Goals / Non-Goals

**Goals:**
- Produce a frozen, raw-count, donor-level artifact (`results/pseudobulk/`) that the export and MOTCO steps consume without re-reading the 51 GB of donor objects.
- Make every choice that affects the cohort or the counts (labels, QC, exclusions, thresholds) explicit, recorded and reproducible from the release files.
- Keep one code path that runs unchanged on the cluster and on any workstation, differing only by the data root.

**Non-Goals:**
- Any transform, feature filter, covariate handling or MOTCO call (follow-up changes).
- Touching the existing `subset`/`de` DESeq2 path or any existing default; MTG is opt-in via `--omic` and the new command.
- A test framework: verification is a two-file spike, a two-donor dry run and a full cluster run whose cohort summary is checked against the preregistered numbers.

## Decisions

**D1. Data root via one setting, S3 key layout preserved.**
A `resolve_data_root(explicit)` helper applies CLI flag → `SEAD_DATA_ROOT` → `data/`. Local paths mirror S3 keys (`<root>/MTG/RNAseq/donors_objects/...`), so the download manifest and the pseudobulk manifest can reference the same relative keys on every machine. Alternative considered: a config file per machine. Rejected: an env var in `.envrc` (direnv) and in the sbatch script is enough and has no parsing surface.

**D2. Own sparse aggregation instead of `decoupler.pp.pseudobulk`.**
Per donor, the pseudobulk row is `X[mask].sum(axis=0)` on the CSR matrix (or the `UMIs` layer), cast to `int64`. `decoupler` groups by `Supertype`, filters samples by `min_cells`/`min_counts` and would have to be re-aggregated; writing the sum directly keeps the operation transparent and avoids an implicit sample filter, which here is a cohort decision owned by the donor table. `decoupler` stays a dependency for the legacy path only.

**D3. Taxonomy mapping as a versioned constant keyed by release.**
A module-level mapping `{release_date: {"astrocyte_subclass": [...], "supertypes": {raw: canonical}}}` with the release date parsed from the file name (`...final-nuclei.<date>.h5ad`). The run validates that every observed label among astrocyte nuclei is in the mapping and writes the mapping to `manifest.json`. The spike (task 1) fills the actual strings. Alternative: the `SEA-AD-cell-annotation` CSV in the bucket. Rejected: the per-donor `obs` already carries the labels, and a second file adds a join without removing the need to reconcile releases.

**D4. Nucleus QC: RNA `min_genes=200`, ATAC none.**
RNA keeps the repo's existing rule, applied to astrocyte nuclei before summing. ATAC uses the release's `final-nuclei` QC only; an extra fragment- or TSS-based filter would be a new scientific choice without a validated threshold. Both are recorded. Changing either is a flagged decision, not a default edit.

**D5. Per-donor cache, then assembly.**
Each donor writes `results/pseudobulk/per_donor/<omic>/<donor>.npz` (sparse 1 × features row, nucleus count, supertype counts, sex, ADNC) as soon as it is done. An assembly step stacks the cached rows into one `AnnData` per omic (`rna_counts.h5ad`, `atac_counts.h5ad`: `X` int CSR, `obs` = donor table, `var` = feature annotations copied from the first donor file) and writes the CSV tables and `manifest.json`. This gives resumability on the cluster and a cheap `--limit N` dry run. Alternative: a single in-memory pass. Rejected: a crash at donor 80 would redo 51 GB of reads.

**D6. Output format: AnnData `.h5ad` per omic plus CSV tables.**
`h5ad` keeps `var` annotations (gene ids, peak coordinates) that the export step and later enrichment need, and is what the rest of the stack reads. Sizes: RNA ≈ 80 × 36k int → ~10 MB; ATAC ≈ 86 × 219k sparse int → tens of MB. CSV tables (`donor_table.csv`, `supertype_composition.csv`, `cohort_summary.csv`) are the human-readable contract. Alternative: parquet. Rejected: no annotation carrier, second reader.

**D7. Donor table consistency rules.**
`sex` and ADNC are read from `obs` of both omics; donors in both must agree or the run errors. If the spike shows the ATAC `obs` lacks these columns, they are taken from the RNA objects (all 84 RNA donors have ATAC) and the manifest records the source; the two ATAC-only donors then have no stage and no cohort flag. Stage labels are the strings `0`, `1`, `2` so the later `stage_order` is unambiguous.

**D8. Cluster execution: one sequential sbatch job, not an array.**
The step is IO-bound and sequential by design (D1 of the specs), so a single job (`scripts/pseudobulk.sbatch`: 4 CPUs, 16 GB, 12 h, `OMP/OPENBLAS/MKL_NUM_THREADS=1`, `SEAD_DATA_ROOT=/home1/tgonza/sead_sexdiff/data`) running `pseudobulk --omic both` is simpler than an array and leaves the per-donor cache for restarts. Compute nodes have no internet access (found during apply, 2026-10-16), so `download --omic both` runs on the login node under `nohup` beforehand, and the job refuses to start unless all 84 RNA and 86 ATAC objects are present. The `full` partition's live node (n14) has glibc 2.17 and cannot load the wheels installed from the login node; use `512x1024`. Environment on the cluster: a clone at `/home1/tgonza/sead_sexdiff`, `uv sync --locked` with uv allowed to fetch Python 3.12 (the Nix flake and its `UV_PYTHON_DOWNLOADS=never` apply only to workstations). A `scripts/sync_pseudobulk.sh` wraps `rsync -av cluster.ing.uc.cl:/home1/tgonza/sead_sexdiff/results/pseudobulk/ results/pseudobulk/`. Alternative: run on the office workstation only. Rejected per the exploration: three workstations, one canonical copy.

**D9. CLI shape.**
New subcommand `pseudobulk` with `--omic {rna,atac,both}`, `--data-root`, `--out-dir` (default `results/pseudobulk`), `--limit N`, `--donors id,...`, `--force`. `download` gains an optional `--omic` that overrides `--prefix` with the MTG per-omic prefix and sets `suffixes=None`; without it, every existing flag and default behaves exactly as today. `subset` and `de` are not modified.

## Risks / Trade-offs

- [Astrocyte or supertype label strings differ from expectation, or `Astro_*` names differ between releases] → spike on one RNA + one ATAC file before coding the mapping; the mapping is data, so a correction is a one-line change recorded in the manifest.
- [ATAC `obs` lacks `Sex`/ADNC, or ATAC counts live in a layer rather than `X`] → D7 fallback for metadata; the matrix source is a parameter verified for integer values and recorded.
- [MTG RNA donor objects lack a `UMIs` layer or use a different raw-count layer name] → spike checks layers; the layer name is a parameter with `UMIs` default and an explicit error when absent.
- [Feature set differs across donors of one release] → hard error by spec; if it occurs, the fix (outer join with zeros, or a reference feature list) is a recorded decision in a follow-up, not an implicit alignment.
- [Cohort counts do not reproduce n = 80 (F 11/10/28, M 9/10/12)] → the summary is compared and the deviation written to the manifest; the ≥30 rule or the QC may need revisiting with the user, since the preregistered numbers bind the MOTCO operating characteristics.
- [`uv sync` on the cluster needs network to fetch Python 3.12 and wheels] → run it on the login node once before submitting; the sbatch job only calls the synced environment.
- [51 GB download time on the cluster] → size-checked idempotent download; the job can be resubmitted.

## Migration Plan

1. Land code and docs; existing commands and defaults are unchanged.
2. Clone the repo on `ing`, set `SEAD_DATA_ROOT`, run `uv sync`, submit `scripts/pseudobulk.sbatch`.
3. Sync `results/pseudobulk/` to each workstation with `scripts/sync_pseudobulk.sh`.
4. Rollback: outputs are additive under `results/pseudobulk/`; existing commands are untouched.

## Open Questions

- Exact `Subclass`/`Supertype` label strings for astrocytes in each release, the ATAC count matrix location, and whether ATAC `obs` carries `Sex` and ADNC. Answered by the spike task; none of them changes the specs or the task breakdown, only mapping data and one fallback branch.
- Which `var` columns (gene ids, peak coordinates) each release provides. Copied through as-is; the export change decides what it needs.

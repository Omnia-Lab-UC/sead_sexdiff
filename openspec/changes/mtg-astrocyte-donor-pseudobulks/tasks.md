# Tasks

## 1. Spike: confirm release contents (two files, ~0.8 GB)

- [x] 1.1 Download one RNA donor object (`H19.33.004_SEAAD_MTG_RNAseq_final-nuclei.2026-06-22.h5ad`) and one ATAC donor object (`H19.33.004_SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad`) with the existing `download --file` path into the data root; verify both files exist at the remote size
- [x] 1.2 Inspect both objects in backed mode and record in `openspec/changes/mtg-astrocyte-donor-pseudobulks/spike.md`: `obs` columns present (`Subclass`, `Supertype`, `Donor ID`, `Sex`, `Overall AD neuropathological Change`), the astrocyte `Subclass` label and the `Astro_*` supertype labels per release, RNA layer names (is `UMIs` present), the ATAC count matrix location and dtype, and the `var` columns; verify the note answers every design open question

## 2. Data root and download

- [x] 2.1 Add `resolve_data_root(explicit)` in `data.py` (CLI flag → `SEAD_DATA_ROOT` → `data/`) and use it in `download`; verify by running `download` with and without the env var and checking the printed root and the destination
- [x] 2.2 Add an optional `--omic {rna,atac,both}` to `download` mapping to the MTG prefixes, leaving the existing `--prefix`/`--suffix`/`--file` behaviour and defaults unchanged when it is absent; verify `download --omic atac --file ...` fetches one file under the preserved key and `download --help` still shows `PFC/RNAseq/` as the prefix default
- [x] 2.3 Make downloads report downloaded/skipped counts and bytes, and write `download_manifest.json` (key, size, etag, timestamp) under the downloaded prefix; verify a second run of the two spike files reports 0 downloaded, 2 skipped and the manifest lists both

## 3. Taxonomy mapping and per-donor pseudobulk

- [x] 3.1 Add the per-release taxonomy mapping constant (release date parsed from the file name, astrocyte subclass label(s), supertype raw → canonical) filled from the spike, with a validator that errors on an unmapped label naming file and label; verify with the two spike files and with a deliberately wrong label
- [x] 3.2 Implement `pseudobulk_donor(path, omic, ...)`: read one object, exclude `Reference` ADNC, select astrocytes via the mapping, apply `min_genes=200` for RNA only, sum the raw matrix (`UMIs` layer for RNA; the spike-confirmed source for ATAC) to an integer row, check non-negative integers, collect nucleus count, supertype counts, sex and ADNC, free the object; verify on the two spike files that the row equals a manual sum and counts match `obs` tallies
- [x] 3.3 Implement the per-donor cache (`results/pseudobulk/per_donor/<omic>/<donor>.npz`) with skip-if-present and `--force`; verify a rerun on the spike files reads nothing and a `--force` rerun rewrites them
- [x] 3.4 Implement the feature-set consistency check across donors (same ordered `var_names`, error naming file and count of differing features); verify by running on the two spike files plus one with a shuffled `var` in a scratch copy
- [x] 3.5 Implement assembly into `rna_counts.h5ad` / `atac_counts.h5ad` (int CSR `X`, `var` copied from the first donor, `obs` joined from the donor table) and `supertype_composition.csv`; verify the matrices hold only non-negative integers and composition totals equal per-donor nucleus counts

## 4. Donor table and cohorts

- [x] 4.1 Implement ADNC recoding to string stages `0`/`1`/`2` with an error on unexpected values, and the sex two-label and cross-omic agreement checks; verify with a small synthetic frame covering every branch
- [x] 4.2 Implement `donor_table.csv` (one row per donor in either omic, zero counts for the missing omic, `cohort_30`/`cohort_50` requiring both omics) and `cohort_summary.csv` (sex × stage per cohort, totals, empty-cell flag); verify on a synthetic frame that the borderline donor (45/60) is in `cohort_30` only and an empty cell is flagged
- [x] 4.3 Write `manifest.json` (release files and sizes from the download manifests, taxonomy mapping, QC parameters, matrix source per omic, excluded donors, thresholds, partial-run flag, observed vs preregistered cohort counts, software versions, git revision); verify it is valid JSON and contains every listed key after a dry run

## 5. CLI wiring

- [x] 5.1 Add the `pseudobulk` subcommand (`--omic`, `--data-root`, `--out-dir`, `--limit`, `--donors`, `--force`) printing `[pseudobulk]` progress per donor and the cohort summary at the end; verify `uv run sead-sexdiff pseudobulk --help` and a `--limit 1` run on the spike files
- [x] 5.2 Confirm `subset` and `de` are untouched: `uv run sead-sexdiff subset --help` and `de --help` show the same defaults as before the change

## 6. Cluster execution and sync

- [x] 6.1 Write `scripts/pseudobulk.sbatch` (partition `512x1024`, 4 CPUs, 16 GB, 12 h, BLAS threads pinned to 1, `SEAD_DATA_ROOT=/home1/tgonza/sead_sexdiff/data`, runs `download --omic both` then `pseudobulk --omic both`, logs under `results/pseudobulk/logs/`); verify with `sbatch --test-only` on the cluster
- [ ] 6.2 Write `scripts/sync_pseudobulk.sh` (rsync of `results/pseudobulk/` from the cluster, excluding `per_donor/` by default with a flag to include it); verify a dry run (`-n`) lists the expected files
- [x] 6.3 Clone the repo on `ing` at `/home1/tgonza/sead_sexdiff`, run `uv sync --locked` there (uv may fetch Python 3.12), and verify `uv run sead-sexdiff pseudobulk --help` works on the login node

## 7. Documentation

- [x] 7.1 Update `README.md`: add the MTG section (`download --omic`, `pseudobulk`), the data-root setting, output contract of `results/pseudobulk/`, cluster workflow and sync, leaving the existing PFC/DESeq2 documentation as is; verify every documented command parses with `--help`
- [x] 7.2 Update `AGENTS.md`: new pipeline steps, data root, taxonomy mapping location, QC and cohort rules as flagged scientific choices, cluster conventions; verify `CLAUDE.md` still imports it unchanged

## 8. Verification runs

- [x] 8.1 Two-donor dry run on a workstation (`pseudobulk --omic both --limit 2`) and verify all outputs exist, the manifest is marked partial, and the composition totals match the donor table
- [ ] 8.2 Full run on the cluster via `scripts/pseudobulk.sbatch`; verify 84 RNA and 86 ATAC donors processed, zero errors in the log, and the ≥30 cohort summary reports n = 80 (F 11/10/28, M 9/10/12) and the ≥50 cohort n = 74, with any deviation recorded in the manifest and reported to the user
- [ ] 8.3 Sync outputs to the workstation with `scripts/sync_pseudobulk.sh` and verify `rna_counts.h5ad` and `atac_counts.h5ad` load, have the expected shapes, and total under ~150 MB

# Proposal

## Why

The MOTCO Phase 6 real-data case study (see `../MOTCO/docs/roadmap.md`, "Phase 6 — Real-data case study" and "Next changes" item 1) needs donor-level astrocyte RNA + ATAC pseudobulks from the SEA-AD middle temporal gyrus (MTG) releases, with a predefined cohort, stage recoding and supertype accounting. This repo currently targets the prefrontal cortex (PFC/A9) RNA release only, points at stale S3 paths, and produces per-supertype DESeq2 inputs rather than per-donor multi-omic count matrices. MOTCO 0.7.0 shipped `motco analyze` on 2026-10-15 and removed the last method-side blocker, so the data side is now the critical path.

## What Changes

- **New `pseudobulk` pipeline step** that reads the MTG RNA (`MTG/RNAseq/donors_objects/`, 84 donors, release 2026-06-22) and ATAC (`MTG/ATACseq/donor_objects/`, 86 donors, release 2024-12-06) donor objects one file at a time, keeps astrocyte nuclei, and writes per-donor raw-count pseudobulks (donors × genes from the `UMIs` layer; donors × peaks from the raw peak counts) summed across all astrocyte supertypes, plus per-donor nucleus counts and supertype composition for each omic.
- **Taxonomy reconciliation** between the two releases: the RNA release labels microglia `Immune` where the ATAC release uses `Microglia-PVM`, so astrocyte `Subclass`/`Supertype` labels are resolved through an explicit, recorded mapping rather than string equality. The exact label strings are confirmed by a two-file spike before the full run.
- **Donor table and cohort definitions**: one row per donor with `Donor ID`, `Sex`, raw ADNC, the merged three-level stage (`Not AD` + `Low` → `0`, `Intermediate` → `1`, `High` → `2`), astrocyte nucleus counts per omic, and `cohort_30` / `cohort_50` flags (≥30 / ≥50 astrocyte nuclei in **both** omics). A cohort summary reports group × stage counts per cohort; the ≥30 cohort is expected to be n = 80 (F 11/10/28, M 9/10/12) and the ≥50 cohort n = 74.
- **MTG download as an explicit option**: `download --omic {rna,atac,both}` selects the correct MTG per-omic prefix, and a single data-root setting (`SEAD_DATA_ROOT`, default `data/`) lets laptop, office, home and the `ing` cluster differ only by that value. Existing defaults (`--prefix PFC/RNAseq/`, the `subset` folder, the `de` inputs) are unchanged; nothing is breaking.
- **Cluster execution and sync**: an sbatch wrapper runs download + pseudobulk on `cluster.ing.uc.cl` (raw data canonical at `/home1/tgonza/sead_sexdiff/data`), and a documented rsync brings the small derived outputs (`results/pseudobulk/`, ~100 MB) to the workstations where every later step runs.
- **Housekeeping**: `AGENTS.md` (the single source of truth; `CLAUDE.md` already imports it) and `README.md` updated for the new command, the data-root setting and the cluster workflow. The existing `subset`/`de` DESeq2 commands and their defaults are untouched.

Out of scope here and planned as follow-up changes: feature filtering and transforms (log2-CPM, VST, chrX/Y exclusion), the MOTCO input export with declared `identity` units and preprocessing manifest, the MOTCO `v0.7.0` dependency pin and the `motco analyze` run step with its sensitivity grid, and pathway/module enrichment.

## Capabilities

### New Capabilities
- `data-download`: Fetch SEA-AD MTG donor objects (and other bucket files) from the public S3 bucket into a configurable data root, skipping files already present at the expected size.
- `astrocyte-donor-pseudobulk`: Produce per-donor raw-count astrocyte pseudobulks for RNA and ATAC from the MTG donor objects, with reconciled taxonomy labels, nucleus-level QC rules, and per-donor supertype composition.
- `donor-cohort-table`: Produce the donor metadata table with the merged ADNC stage and the ≥30 / ≥50 nuclei cohort flags, plus a cohort summary, consistent across both omics.

### Modified Capabilities
<!-- No existing specs, and the existing subset/de behaviour and defaults are not changed. -->

## Impact

- `src/sead_sexdiff/cli.py`: new `pseudobulk` subcommand; `download` gains `--omic` and data-root resolution. Existing defaults unchanged.
- `src/sead_sexdiff/data.py`: new pseudobulk, taxonomy and cohort functions; `download_data` gains the omic-to-prefix mapping and the manifest. Existing subset/DE functions and defaults untouched.
- New `scripts/` sbatch wrapper and sync helper for the cluster; new `results/pseudobulk/` output contract consumed by the follow-up export change.
- `README.md`, `AGENTS.md` (`CLAUDE.md` already imports it and needs no change).
- No new Python dependencies expected (anndata, scanpy, scipy already present); `decoupler` stays for the legacy path.
- Data: ~51 GB of raw donor objects on the cluster; derived outputs ~100 MB.

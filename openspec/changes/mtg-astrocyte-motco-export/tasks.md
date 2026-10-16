# Tasks

## 1. Donor covariates

- [x] 1.1 Add the covariate reader: per donor object, read the `obs` columns listed in the spec with `_read_obs_columns`, collapse each one with `_single_value`, and leave `X` and the layers unread. Verify on the 3 local spike donors that the values match a backed-mode `obs` read, and that a scratch copy with a varying `PMI` raises an error naming the file and `PMI`.
- [x] 1.2 Add numeric coding: `N+` becomes `N` with `age_capped`, non-numeric values become missing and are listed, and education level stays a string. Verify with a synthetic frame covering `90+`, `80.0`, `nan`, `Unknown` and an empty string.
- [x] 1.3 Add the RNA/ATAC cross-check: RNA values win, disagreements go to the findings, and an RNA-only run records that no check was done. Verify with a synthetic pair holding one PMI disagreement, and on the spike donors (expect no findings).
- [x] 1.4 Add `run_covariates` writing `donor_covariates.csv` and `covariates_manifest.json`, plus the `covariates` CLI subcommand (`--omic`, `--data-root`, `--out-dir`) with `[covariates]` output. Verify `uv run sead-sexdiff covariates --help`, then a local run on the spike donors into a scratch `--out-dir` that produces valid JSON and one row per RNA donor.

## 2. Annotation

- [x] 2.1 Add the GENCODE v32 basic GTF resolver (download into `<data root>/annotation/` if absent, sha256, parse `gene` lines into a version-stripped `gene_id → seqname` map). Verify that it downloads once, that a second call does not re-download, and that every one of the 36,601 `gene_ids` in `rna_counts.h5ad` maps, or report how many do not. *Observed 2026-10-16: 60,564 genes in the GTF (sha256 `a61ab50f…`); 29 of 36,601 `gene_ids` unmapped: 24 on unplaced/unlocalized scaffolds (KI270…/GL000…, absent from the reference-chromosome basic GTF) and 5 retired ids. None can be on chrX/Y/M; all are kept. Annotated: chrX 1,148, chrY 111, chrM 13.*
- [x] 2.2 Add the ATAC peak-name chromosome parser, which raises an error on an unparseable name. Verify on `atac_counts.h5ad` that it reproduces the per-chromosome counts (chrX 5,886, chrY 180, no chrM), and that a malformed name raises an error.

## 3. Normalization and selection

- [x] 3.1 Implement TMM factors (edgeR algorithm, D4) and edgeR-style log2-CPM with prior count 2. Verify on a small synthetic matrix with known factors (identical samples → 1.0; a sample with doubled depth → factor 1.0 and log-CPM equal to the others).
- [ ] 3.2 Add `rPackages.edgeR` to `flake.nix` and a `scripts/check_tmm_edger.R` cross-check that reads the filtered `cohort_30` counts written by a debug hook. Verify TMM factors agree with `calcNormFactors` to 1e-8 and log-CPM agrees with `cpm(log=TRUE, prior.count=2)` to 1e-10. Record the observed maximum differences in this task's line.
- [x] 3.3 Implement the per-cohort pipeline (restrict → drop chrX/Y/M → CPM ≥ 1 in ≥ k donors, with k from the smallest sex × stage cell → TMM → log2-CPM → top 5,000 by variance, ties by name), returning per-feature status and step counts. Verify on a synthetic cohort that the step counts add up, k equals the smallest cell, and a block with fewer than N surviving features keeps all of them and raises a finding.

## 4. Export

- [x] 4.1 Implement input checks: refuse a partial pseudobulk manifest, missing files (naming `covariates` when `donor_covariates.csv` is absent) and non-analysable cohorts. Verify each refusal against a scratch copy of `results/pseudobulk/`.
- [x] 4.2 Implement the writers for `rna.csv`, `atac.csv`, `metadata.csv` (with `stage_label`), `covariates.csv` (including `n_rna`/`n_atac`) and `features_<omic>.csv`, at `%.10g`, sorted by `donor_id`, written to `<cohort>.tmp/` and then renamed, refusing an existing directory without `--force`. Verify that donor columns are identical across files, and that a rerun without `--force` raises an error naming the directory.
- [x] 4.3 Write the per-cohort `manifest.json` with every key in the spec, including the input-manifest hash, TMM factors, the annotation sha256, the "covariates exported, not applied" note and findings. Verify it is valid JSON containing every listed key.
- [x] 4.4 Add the `export` CLI subcommand (`--pseudobulk-dir`, `--out-dir`, `--cohort`, `--data-root`, `--force`) printing `[export]` progress and the step-count table. Verify `uv run sead-sexdiff export --help`.

## 5. Real run and checks

- [ ] 5.1 On the cluster, pull, `uv sync --locked`, and run `covariates --omic both` on the login node. Then sync to the workstation with `scripts/sync_pseudobulk.sh`. Verify 84 RNA and 86 ATAC objects were read, and record the number of capped ages, the missing values and the findings in this task's line.
- [ ] 5.2 Run `export --cohort all` on the workstation. Verify: `cohort_30` has 80 rows and `cohort_50` has 74; there are 5,000 + 5,000 columns, or a recorded finding; the `features_*.csv` tallies sum to 36,601 and 218,882; no exported feature lies on chrX/Y/M; and every value is finite.
- [ ] 5.3 Rerun with `--force` and verify the CSVs are byte-identical (`sha256sum`). Load `rna.csv` + `atac.csv` + `metadata.csv` into pandas the way `motco analyze` expects (`--sample-col donor_id --group-col sex --stage-col stage`, stage order `0,1,2`) and verify that every sex × stage cell is non-empty.

## 6. Documentation

- [x] 6.1 Update `README.md` (the `covariates` and `export` commands, outputs, the preprocessing order and constants, the annotation download) and `AGENTS.md` (the pipeline steps, the new constants under the scientific-choices list, and `covariates` running on the cluster login node). Verify the flags and defaults in the docs match `--help`.
- [ ] 6.2 Update the case-study memory with the outcome and the next step (the MOTCO covariate change). Verify `openspec validate mtg-astrocyte-motco-export --strict` passes.

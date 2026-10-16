# Proposal

## Why

The MTG astrocyte pseudobulks (`results/pseudobulk/`, archived change `2026-10-16-mtg-astrocyte-donor-pseudobulks`) are raw counts. MOTCO's `motco analyze` does no normalization, filtering or imputation and takes blocks "in the continuous units you want measured" (`../MOTCO/docs/api/analysis.md`), so the case study has to own preprocessing and hand MOTCO finished blocks. The user also wants covariates (age, PMI, RIN, nuclei counts, education) inside the trajectory regression. MOTCO 0.7.0 has no covariate support, so that is a MOTCO change. This repo has to supply the covariate table first, and the pseudobulk step never recorded donor covariates.

## What Changes

- **New `covariates` step.** For every donor, read the donor-level `obs` columns from the MTG donor objects without loading matrices: `Age at Death`, `PMI`, `RIN`, `Years of education` and `Highest level of education`. Write `results/pseudobulk/donor_covariates.csv` with numeric coding (an age reported as `90+` becomes 90 and is flagged). Values shared by the RNA and ATAC objects are cross-checked, and disagreements are written to a findings list. The step runs on the cluster, where the raw objects live, and the file reaches workstations through the existing `scripts/sync_pseudobulk.sh`.
- **New `export` step.** It reads only `results/pseudobulk/`, so it can run on any workstation. For each cohort (`cohort_30` is primary, `cohort_50` is the sensitivity cohort) and each omic, in fixed order:
  1. Restrict to the cohort's donors.
  2. Drop chrX, chrY and chrM features. RNA genes are matched to chromosomes through a pinned GENCODE annotation using `gene_ids`. ATAC peaks carry the chromosome in their names.
  3. Keep features with CPM ≥ 1 in at least as many donors as the cohort's smallest sex × stage cell (9 in both cohorts).
  4. Compute TMM normalization factors, then log2-CPM with prior count 2 (edgeR semantics).
  5. Keep the 5,000 most variable features per block. Selection is unsupervised, and both blocks get the same N because MOTCO z-scores each feature without weighting blocks.
- **Export contract per cohort** (`results/motco_input/<cohort>/`):
  - `rna.csv` and `atac.csv`: donors × features, in identity units.
  - `metadata.csv`: `donor_id`, `sex`, `stage`, `stage_label`, in declared stage order `0`, `1`, `2`.
  - `covariates.csv`: the covariate table restricted to the cohort, plus `n_rna` and `n_atac`.
  - `features_<omic>.csv`: one row per feature with the step where it was removed, or `kept`.
  - `manifest.json`: units, every threshold, feature counts after each step, TMM factors, annotation source and checksum, input manifest hash, software, git revision, and findings.
- **Covariates are exported, not applied.** No residualization happens here. Using them in the trajectory fit is a MOTCO change (covariate-adjusted LS-means and RRPP, with Type I validation). It has to land before the analysis change in this repo.
- README and AGENTS.md document both commands and the new preprocessing constants. Existing commands and defaults are unchanged.

Out of scope, as follow-up changes: the MOTCO change for covariate adjustment; pinning MOTCO 0.7.0+ here (including the Python ≥ 3.12 install check) and running `motco analyze`; preprocessing sensitivity exports (other N, keeping chrX); pathway/module enrichment.

## Capabilities

### New Capabilities
- `donor-covariate-table`: donor-level covariates read from the MTG donor objects, with numeric coding, cross-omic agreement checks and the `90+` age flag.
- `motco-input-export`: per-cohort preprocessing of the raw pseudobulks (sex-chromosome removal, abundance filter, TMM log2-CPM, top-N variable features) into MOTCO-ready blocks, metadata, covariates and a preprocessing manifest.

### Modified Capabilities
<!-- None: astrocyte-donor-pseudobulk, donor-cohort-table and data-download keep their requirements; the new steps only read their outputs. -->

## Impact

- `src/sead_sexdiff/data.py`: new covariate reader, annotation loader, TMM / log-CPM / feature-selection functions, and `run_covariates` / `run_export`. The pseudobulk, subset and DE functions are untouched.
- `src/sead_sexdiff/cli.py`: new `covariates` and `export` subcommands.
- `scripts/`: no change needed to `sync_pseudobulk.sh`, which already copies all of `results/pseudobulk/`. One extra login-node command is run on the cluster.
- Annotation: GENCODE v32 basic gene annotation (the reference used for GRCh38-2020-A, Cell Ranger 6.0), downloaded once into `<data root>/annotation/` and checksummed. This needs internet on the machine that runs `export`.
- Dependencies: none new. TMM is implemented in numpy. edgeR from Nix is used only for a one-off cross-check, if it is available in the flake.
- Outputs: `results/motco_input/` is a few MB per cohort and is gitignored with `results/`.

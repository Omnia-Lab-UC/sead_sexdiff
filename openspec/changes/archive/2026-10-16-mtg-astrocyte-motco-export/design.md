# Design

## Context

See `proposal.md` for the motivation. The current state that shapes the approach:

- **Pseudobulk outputs** (`results/pseudobulk/`, on the laptop and the cluster):
  - `rna_counts.h5ad`: 84 donors × 36,601 genes, int64 CSR. `var` has `gene_ids` (unversioned Ensembl), `feature_types` and `genome = GRCh38`, with no chromosome column.
  - `atac_counts.h5ad`: 84 × 218,882 peaks, named `chrN:start-end`. `var` has no columns. Only canonical contigs occur: chr1–22, chrX (5,886), chrY (180).
  - The donor table has `sex`, `stage` (`0`/`1`/`2`), `n_rna`, `n_atac` and the cohort flags. The smallest sex × stage cell has 9 donors in both cohorts.
- **Raw donor objects.** These are on the cluster only; the laptop holds 3 spike donors. RNA `obs` carries `Age at Death`, `PMI`, `RIN`, `Years of education` and `Highest level of education`. ATAC `obs` carries the same set except `RIN`. The pseudobulk step never read these columns.
- **The RNA objects** were aligned with Cell Ranger 6.0 on GRCh38 (`obs['alignment']`, `obs['Genome']`). This matches 10x reference GRCh38-2020-A, which is built from GENCODE v32 / Ensembl 98.
- **MOTCO 0.7.0** z-scores each feature over all samples (`fit_block_scalers`) and concatenates blocks without block weights. It applies only `identity` or `beta_to_mvalue`, and it reads CSV blocks via `--block name=path` with the first column as the sample index.
- **Existing helpers in `data.py` to reuse:** `_read_obs_columns` (h5py `obs` reads), `_single_value`, `donor_id_of`, `resolve_data_root`, `MTG_PREFIXES`, `_git_revision`, `_software_versions`.

## Goals / Non-Goals

**Goals:**
- Every preprocessing choice is a named module constant, recorded in the manifest, and traceable per feature (`features_<omic>.csv`).
- `export` runs on a workstation from synced outputs. Only `covariates` needs the raw objects.
- The output can be fed to `motco analyze --block rna=… --block atac=… --metadata metadata.csv --sample-col donor_id --group-col sex --stage-col stage --stage-order 0,1,2` without any reshaping.

**Non-Goals:**
- Covariate adjustment of any kind, including residualization. That belongs to the MOTCO change.
- Calling or installing MOTCO.
- Sensitivity exports (other N, keeping chrX, VST). The constants make them a later flag, not part of this change.
- A pytest suite. Verification is synthetic checks plus the real run, as in the previous change.

## Decisions

**D1. Two commands, split by where the data lives.**
`covariates` reads `obs` columns from the raw objects and runs on the cluster login node. It reads metadata only, so it takes seconds per file and needs no sbatch job. It writes into `results/pseudobulk/`, so the existing `sync_pseudobulk.sh` carries the file to the workstations unchanged. `export` reads only `results/pseudobulk/` plus the annotation.
- *Alternative: add covariates to `pseudobulk`.* Rejected: it would change an archived spec's outputs, and the per-donor cache skips cached files, so covariates would need a forced 51 GB re-read or a cache-format change.
- *Alternative: read covariates during `export`.* Rejected: workstations do not hold the raw objects.

**D2. Covariates come from RNA, and ATAC is used as a cross-check.**
Every cohort donor has an RNA object, and `RIN` exists only in RNA. When ATAC objects are present, the shared columns are compared and any disagreement goes to the findings (see the spec). Age coding: `str.strip`, then a regex `^(\d+(?:\.\d+)?)\+$` turns `N+` into `N` and sets the capped flag. Any other value that `float` cannot parse becomes missing and is listed. `education_level` stays a string; it is the categorical fallback that the user asked to keep beside `years_education`.

**D3. GENCODE v32 annotation, fetched once and checksummed.**
`export` resolves `<data root>/annotation/gencode.v32.basic.annotation.gtf.gz`. If the file is absent, it downloads it from the EBI GENCODE FTP (`https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/`). It then parses only `gene` lines into `gene_id (version stripped) → seqname`, and records the URL, the sha256 and the number of mapped genes.
- *Alternative: BioMart / pyensembl.* Rejected: the result depends on a live service and drifts by release.
- *Alternative: commit a 1 MB gene → chromosome CSV.* Rejected: it carries no provenance. The checksummed GTF is the provenance.
- The **basic** GTF holds every gene, so the gene set is the same as in the comprehensive file and the download is smaller. The full GTF is ~45 MB gzip, under the data root, and gitignored. 10x drops some gene biotypes from its reference, so mismatches can only go one way: every `gene_ids` entry is expected to map. Unmapped genes are kept and counted (see the spec). A large unmapped count would point to the wrong release and is visible as a finding.

**D4. TMM implemented in numpy, mirroring edgeR's `calcNormFactors(method = "TMM")`.**
The implementation follows edgeR:
- Reference sample: the one whose upper-quartile factor is closest to the mean.
- Per-sample trimming: 30 % of M values, 5 % of A values.
- Precision weights from the delta method.
- Only features with nonzero counts in both samples are used.
- Factors are rescaled so their geometric mean is 1.
Log2-CPM follows edgeR `cpm(y, log = TRUE, prior.count = 2)`: effective library size `L = lib × factor`; prior `p = 2 · L / mean(L)`; value `log2((count + p) / (L + 2p) · 1e6)`.

Correctness is cross-checked once against edgeR on the real `cohort_30` filtered counts. If the factors agree to 1e-8 and log-CPM to 1e-10, the check passes. edgeR runs in a throwaway R session in the Nix shell: `rPackages.edgeR` is added to `flake.nix` because R packages come from Nix there. The check script lives under `scripts/` and is not wired into the CLI.
- *Alternative: rpy2 / an R call inside the export.* Rejected: it puts R on the export's critical path for a ~40-line algorithm.
- *Alternative: scanpy `normalize_total` + `log1p`.* Rejected: that is a different estimator from the one the user chose.

**D5. Library sizes for the abundance filter come from the post-chromosome counts, and TMM runs after the filter.**
This is edgeR's usual order: filter on CPM, then normalize the kept features. Removing chrX/Y first means sex-chromosome reads do not inflate the library sizes of the sexes differently. The filter's CPM uses raw library sizes (no TMM), which matches `filterByExpr`'s use of plain CPM. `k` is computed from the cohort's sex × stage table, not fixed at 9, so a change to the cohorts cannot leave the filter stale.

**D6. Variance selection on the final log2-CPM, ties by name.**
The population variance (ddof 0) is taken over the cohort's donors. Features are sorted by (−variance, name) and the first 5,000 are kept. Selection never looks at labels, so the RRPP null is not conditioned on the group (sex) or on the stages. Sorting by name makes the output deterministic. The exported column order is the input feature order, not the variance rank, so the block layout stays stable across reruns.

**D7. Constants.**
The new module-level constants in `data.py` are:
- `EXPORT_EXCLUDED_CHROMOSOMES = ("chrX", "chrY", "chrM")`
- `EXPORT_MIN_CPM = 1.0`
- `EXPORT_PRIOR_COUNT = 2.0`
- `EXPORT_TOP_N = 5000`
- `EXPORT_COHORTS = ("cohort_30", "cohort_50")`
- `STAGE_LABELS = {"0": "Not AD + Low", "1": "Intermediate", "2": "High"}`
- `GENCODE_RELEASE = "32"` with its URL

None of them is a CLI flag in this change. Changing one is a scientific decision (AGENTS.md), and adding sensitivity flags is a later change.

**D8. Output format: CSV with float64 written at `%.10g`.**
These are the formats MOTCO's CLI reads. `%.10g` makes the CSVs byte-identical across reruns while keeping more precision than the data warrants. Writing is atomic: everything goes into `<cohort>.tmp/`, which is then renamed. A crash therefore never leaves a half cohort that a later run without `--force` would refuse to overwrite.

**D9. CLI shape.**
- `covariates [--omic rna|both] [--data-root] [--out-dir results/pseudobulk]`. The default is `both`; `rna` skips the cross-check.
- `export [--pseudobulk-dir results/pseudobulk] [--out-dir results/motco_input] [--cohort cohort_30|cohort_50|all] [--data-root] [--force]`. `--data-root` is used only to locate the annotation.

Both print `[covariates]` / `[export]` progress. `export` ends with a per-omic table of feature counts after each step.

## Risks / Trade-offs

- [The 10x reference is not exactly GENCODE v32, so some `gene_ids` do not map] → those genes are kept and counted as a finding. If more than a handful are unmapped, stop and reconsider the annotation release before running the analysis.
- [Numpy TMM differs subtly from edgeR (trimming at ties, weights)] → the one-off edgeR cross-check (D4) on the real counts, with the tolerance recorded in the change.
- [5,000 + 5,000 features is still 20× MOTCO's validated 498] → this is a known limitation, carried with the structural-analogue limitation in the analysis change. A smaller N is a later sensitivity export.
- [Variance-ranked HVFs on log-CPM favour low-count, high-dispersion features] → the abundance filter removes the worst of them. The manifest records the mean-log-CPM range of the kept features so the bias can be seen.
- [`Age at Death` may be capped for some donors] → the capped flag is exported, and the MOTCO change decides how to model it.
- [Removing chrX loses real X-linked astrocyte biology] → this was the user's chosen trade-off: sex is the group, so X-escape genes would give trivial group separation. The `features_*.csv` status column makes it reversible later.

## Migration Plan

1. Merge. On the cluster, `git pull` and `uv sync --locked`, then run `sead-sexdiff covariates --omic both` on the login node with `SEAD_DATA_ROOT=/home1/tgonza/sead_sexdiff/data`.
2. On a workstation, run `scripts/sync_pseudobulk.sh` and then `sead-sexdiff export --cohort all` in the Nix shell (it downloads the annotation once).
3. Rollback: delete `results/motco_input/` and `results/pseudobulk/donor_covariates.csv` / `covariates_manifest.json`. No existing output is modified.

## Open Questions

- Whether any MTG donor has `Age at Death = 90+`. The coding handles it either way, and the count is recorded.

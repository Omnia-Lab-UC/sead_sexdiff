# sead_sexdiff

Utilities and CLI for SEA-AD sex-difference analyses using Python >= 3.12. Includes tools to download data from the public SEA-AD S3 bucket and to subset/concatenate AnnData (.h5ad) donor files.

Project layout:

```
.
├── data/            # place raw/input data here (ignored by package build)
├── reports/         # committed case-study reports (numbers generated from results/)
├── results/         # analysis outputs, figures, tables, etc. (gitignored)
├── scripts/         # cluster sbatch wrapper and output sync helper
├── src/
│   └── sead_sexdiff/
│       ├── cli.py
│       ├── data.py
│       └── stage_de_limma.R   # R side of `stage-de-limma` (limma-voom), run with Rscript
├── pyproject.toml   # project configuration (uv/pip compatible)
├── LICENSE
└── README.md
```

Quick start (with uv):

1) Ensure you have Python 3.12+ installed.

2) Install uv (if not already):

```
curl -LsSf https://astral.sh/uv/install.sh | sh
```

3) Create and use a virtual environment and install deps:

```
uv venv --python 3.12
source .venv/bin/activate  # or .venv\Scripts\activate on Windows
uv pip install -e .
```

This installs the project dependencies declared in `pyproject.toml` (e.g., `pandas`, `numpy`, `boto3/botocore`, `anndata`, `scanpy`).

### How to use

After installing the package (see Quick start), you can download data with the CLI:

- Default behavior (downloads CSV and HTML under the default prefix `PFC/RNAseq/` into `data/`):

```
sead-sexdiff download
```

- Specify a prefix and multiple suffixes (repeatable and/or comma-separated):

```
sead-sexdiff download --prefix PFC/RNAseq/ --suffix .csv --suffix .html
# or
sead-sexdiff download --suffix .csv,.html
```

- Choose explicit files to download (overrides suffix filtering). You can pass
  keys relative to the prefix or full keys that already start with the prefix.

```
sead-sexdiff download \
  --prefix PFC/RNAseq/ \
  --file public_datasets/SEA-AD-cell-annotation.2024-08-27.csv \
  --file SEAAD_A9_RNAseq_all-nuclei_metadata.2024-02-13.csv

# or comma-separated
sead-sexdiff download --file public_datasets/SEA-AD-cell-annotation.2024-08-27.csv,SEAAD_A9_RNAseq_all-nuclei_metadata.2024-02-13.csv
```

- Change destination directory and bucket if needed:

```
sead-sexdiff download --dest-root ./data --bucket sea-ad-single-cell-profiling
```

The command prints the resolved data root and the local paths of files it downloaded. It uses anonymous S3 access for the public SEA-AD bucket and fetches `Changelog.html` under the chosen prefix when one exists.

- Files already present locally at the remote size are skipped; a truncated file is re-downloaded. Each run reports downloaded / skipped / failed counts and bytes transferred, exits non-zero if any file failed, and updates `<root>/<prefix>/download_manifest.json` (S3 key, size, ETag, timestamp per file).

### Data root

Every command that reads or writes raw data resolves the data root as: the explicit flag (`download --dest-root`, `pseudobulk --data-root`), then the `SEAD_DATA_ROOT` environment variable, then `data/` in the current directory. Local paths mirror the S3 keys under that root. Set it once per machine, e.g. in `.envrc` or in the sbatch script:

```
export SEAD_DATA_ROOT=/home1/tgonza/sead_sexdiff/data
```

(The legacy `subset` command still takes an explicit `--folder`.)

## MTG astrocyte donor pseudobulks (RNA + ATAC)

Donor-level astrocyte pseudobulks from the SEA-AD middle temporal gyrus (MTG) releases, for the MOTCO case study, and the same pipeline on a second, registered cell type (microglia) as a positive control.

### Cell types (`--cell-type`)

The nine MTG commands `pseudobulk`, `covariates`, `technical-covariates`, `export`, `analyze`, `axis-diagnostics`, `latent-sensitivity`, `stage-de` and `stage-de-limma` take `--cell-type` (default `astrocyte`), a name of the registry `CELL_TYPES` in `src/sead_sexdiff/data.py`. The name is a flag; the mapping is not (entries are added by a change proposal):

| cell type | RNA 2026-06-22 | ATAC 2024-12-06 | not selected | results root |
|---|---|---|---|---|
| `astrocyte` (default) | `Subclass` `Astrocyte`; `Astro_1`–`Astro_5`, `Astro_6-SEAAD` | same | | `results/` |
| `microglia` | `Subclass` `Immune`; `Micro-PVM_1`, `Micro-PVM_2`, `Micro-PVM_2_1-SEAAD`, `Micro-PVM_2_3-SEAAD`, `Micro-PVM_3-SEAAD`, `Micro-PVM_4-SEAAD` | `Subclass` `Microglia-PVM`; the same six | `Lymphocyte`, `Monocyte` | `results/microglia/` |

The microglia labels and preregistered cohort cells come from a metadata-only spike over all 170 donor objects (`openspec/changes/archive/2026-10-22-mtg-microglia-positive-control/spike.md`); neither release has a PVM-only supertype, so the six `Micro-PVM_*` labels (microglia and perivascular macrophages) are kept in both. Any other `Supertype` of the selected subclass stops the run.

The flag sets the population (one selection path, `select_nuclei`, shared by `pseudobulk` and `technical-covariates`; `select_astrocyte_nuclei` is its astrocyte wrapper) and every default directory: `<root>/pseudobulk/`, `<root>/motco_input/`, `<root>/motco_analysis/`, `<root>/motco_diagnostics/`, `<root>/motco_rank_sensitivity/` (read by `axis-diagnostics`), `<root>/motco_latent_sensitivity/`, `<root>/stage_de/`, `<root>/stage_de_limma/` (`cell_type_dirs`). Explicit directory flags take precedence. Thresholds, stage coding and every export, analysis and stage-DE setting are the same for every cell type.

Every manifest of these steps records `cell_type`, and every reader checks it: `export` and `technical-covariates` refuse a pseudobulk manifest of another cell type; `analyze`, `axis-diagnostics`, `latent-sensitivity`, `stage-de` and `stage-de-limma` an export manifest of another; `stage-de-limma` also a DESeq2 manifest of another; `latent-sensitivity` also an analysis, technical or pseudobulk manifest of another. Each refusal names the file and both cell types and writes nothing. A manifest written before the registry existed has no `cell_type`; it reads as `astrocyte` and the reader records that as a finding, so the existing astrocyte outputs are not rerun. The per-donor caches record the cell type and a cache of another one is refused; the stage-DE fit-record fingerprints do not include it (their input hashes already differ between cell types). The sensitivity commands `rank-sensitivity` and `preprocessing-sensitivity` have no flag and read astrocytes; `latent-sensitivity` has the flag (its registry of arms, covariate set and settings are the same for every cell type), and `positive-control` reads both cell types' roots by role and takes `--basis` instead.

```
sead-sexdiff pseudobulk --cell-type microglia --omic both      # results/microglia/pseudobulk/
sead-sexdiff export --cell-type microglia                      # results/microglia/motco_input/<cohort>/
sead-sexdiff analyze --cell-type microglia                     # results/microglia/motco_analysis/<cohort>/
sead-sexdiff latent-sensitivity --cell-type microglia          # results/microglia/motco_latent_sensitivity/
```

### Download the MTG donor objects

`--omic` selects the MTG per-donor object prefix and downloads every file under it (~33 GB RNA, ~18 GB ATAC):

| `--omic` | prefix | donors | release |
|---|---|---|---|
| `rna` | `MTG/RNAseq/donors_objects/` | 84 | 2026-06-22 |
| `atac` | `MTG/ATACseq/donor_objects/` | 86 | 2024-12-06 |
| `both` | both of the above | | |

```
sead-sexdiff download --omic both
# a single donor (with --omic both, each --file goes to the omic named in it)
sead-sexdiff download --omic atac --file H19.33.004_SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad
```

Without `--omic`, `--prefix` and its default (`PFC/RNAseq/`, CSV and HTML files) behave as before.

### Build the pseudobulks

```
sead-sexdiff pseudobulk --omic both               # full run
sead-sexdiff pseudobulk --omic both --limit 2     # dry run on 2 donors present in both omics
sead-sexdiff pseudobulk --donors H19.33.004,H21.33.020
sead-sexdiff pseudobulk --force                   # recompute cached donors
```

Options: `--cell-type {astrocyte,microglia}` (default `astrocyte`; every directory default below is the cell type's, under its results root, shown here for `astrocyte`), `--omic {rna,atac,both}` (default `both`), `--data-root`, `--out-dir` (default `results/pseudobulk`), `--limit N`, `--donors id,...`, `--force`.

What it does, per donor object (one file in memory at a time):

- Excludes donors whose ADNC is `Reference`.
- Selects the cell type's nuclei through the per-release `Subclass` inventory `TAXONOMY` and the registry `CELL_TYPES` in `src/sead_sexdiff/data.py` (astrocytes: `Subclass == "Astrocyte"`, supertypes `Astro_1`–`Astro_5`, `Astro_6-SEAAD`; microglia: the table above). Any `Subclass` label not in the inventory, or `Supertype` of the selected subclass not in the cell type's mapping, stops the run.
- RNA: sums the raw `UMIs` layer over the selected nuclei with ≥ 200 detected genes. ATAC: sums the raw peak counts in `X` over all selected nuclei (release `final-nuclei` QC only). Both matrices are checked to be non-negative integers; no feature filtering or normalization.
- Requires the same ordered feature set across donors of an omic.
- Caches the result (with the cell type) in `<out-dir>/per_donor/<omic>/<donor>.npz`; a rerun skips cached donors, and a cache of another cell type is an error.

Outputs in `results/pseudobulk/`:

| file | content |
|---|---|
| `rna_counts.h5ad` | donors × 36,601 genes, int64 CSR raw UMI sums; `obs` = donor table, `var` = gene annotations |
| `atac_counts.h5ad` | donors × 218,882 peaks, int64 CSR raw count sums; `obs` = donor table |
| `donor_table.csv` | `donor_id`, `sex`, `adnc_raw`, `stage`, `n_rna`, `n_atac`, `cohort_30`, `cohort_50` (one row per non-Reference donor in either omic) |
| `supertype_composition.csv` | `donor_id`, `omic`, `supertype`, `n_nuclei` (nuclei summed per donor × supertype) |
| `cohort_summary.csv` | donors per sex × stage for each cohort, with totals |
| `manifest.json` | cell type, release files (key, size, ETag), the cell type's taxonomy mapping per release, QC, matrix source, excluded donors, thresholds, cohort counts vs the cell type's preregistered cells, partial flag, software versions, git revision |

Donor table rules (scientific choices; change only deliberately):

- `stage`: `Not AD` and `Low` → `0`, `Intermediate` → `1`, `High` → `2`. Other values are an error.
- Sex and ADNC must agree between the RNA and ATAC objects of a donor; sex must have exactly two labels (a warning only in partial runs).
- `cohort_30` / `cohort_50`: ≥ 30 / ≥ 50 nuclei of the cell type in **both** omics.
- A full run compares the cohorts with the cell type's preregistered counts (astrocyte ≥ 30: n = 80, F 11/10/28, M 9/10/12 over stages 0/1/2; ≥ 50: n = 74; microglia, from the metadata-only spike, ≥ 30: n = 76, F 11/10/24, M 9/11/11; ≥ 50: n = 63, F 8/10/18, M 8/10/9) and records any deviation under `findings` in the manifest. A cohort with an empty sex × stage cell is flagged as not analysable.

### Cluster workflow (`cluster.ing.uc.cl`)

The raw donor objects live once on the cluster; workstations only pull the derived outputs (~100 MB).

```
# once, on the login node
git clone -b mtg-astrocyte-pseudobulks https://github.com/Omnia-Lab-UC/sead_sexdiff.git /home1/tgonza/sead_sexdiff  # the cluster has no GitHub SSH key
cd /home1/tgonza/sead_sexdiff && ~/.local/bin/uv sync --locked

# download on the login node (compute nodes have no internet access)
export SEAD_DATA_ROOT=/home1/tgonza/sead_sexdiff/data
mkdir -p results/pseudobulk/logs
nohup .venv/bin/sead-sexdiff download --omic both --no-progress > results/pseudobulk/logs/download.log 2>&1 < /dev/null &

# when the download log reports 0 failed, submit pseudobulk --omic both
sbatch scripts/pseudobulk.sbatch
```

`scripts/pseudobulk.sbatch` uses partition `512x1024`, 4 CPUs, 16 GB, 12 h, BLAS threads pinned to 1, `SEAD_DATA_ROOT=/home1/tgonza/sead_sexdiff/data`, and logs to `results/pseudobulk/logs/`. It stops early unless all 84 RNA and 86 ATAC objects are present. Resubmitting the download or the job is safe: complete downloads and cached donors are skipped. Do not submit to the `full` partition: its node n14 has glibc 2.17 and cannot load the wheels `uv sync` installs on the login node (glibc 2.34); `512x1024` and `gpus` nodes match the login node.

Another cell type: both sbatch scripts read `CELL_TYPE` (default `astrocyte`) and pass it as `--cell-type`. Their `#SBATCH` log paths are the astrocyte ones, so give the cell type's own (Slurm needs the directory to exist):

```
mkdir -p results/microglia/pseudobulk/logs
sbatch --export=ALL,CELL_TYPE=microglia \
    --output=results/microglia/pseudobulk/logs/pseudobulk-%j.out \
    --error=results/microglia/pseudobulk/logs/pseudobulk-%j.err scripts/pseudobulk.sbatch
```

On a workstation, pull the outputs (the per-donor cache is excluded unless `--with-cache`; extra arguments go to rsync; `--cell-type NAME` syncs `results/NAME/pseudobulk/` instead of `results/pseudobulk/`):

```
scripts/sync_pseudobulk.sh -n          # dry run
scripts/sync_pseudobulk.sh
scripts/sync_pseudobulk.sh --with-cache
scripts/sync_pseudobulk.sh --cell-type microglia
```

### Donor covariates

```
sead-sexdiff covariates                 # RNA values, cross-checked against ATAC
sead-sexdiff covariates --omic rna      # RNA only, no cross-check
```

Options: `--cell-type {astrocyte,microglia}` (default `astrocyte`; every directory default below is the cell type's, under its results root, shown here for `astrocyte`), `--omic {rna,both}` (default `both`), `--data-root`, `--out-dir` (default `results/pseudobulk`).

Reads only `obs` columns of every MTG donor object (no matrices; seconds per file), so it needs the raw objects and runs on the cluster login node:

```
cd /home1/tgonza/sead_sexdiff && git pull && ~/.local/bin/uv sync --locked
SEAD_DATA_ROOT=/home1/tgonza/sead_sexdiff/data .venv/bin/sead-sexdiff covariates --omic both
```

Its outputs land in `results/pseudobulk/`, so `scripts/sync_pseudobulk.sh` carries them to the workstations.

| source `obs` column | coded column | omics |
|---|---|---|
| `Age at Death` | `age_at_death`, `age_capped` | RNA, ATAC |
| `PMI` | `pmi` | RNA, ATAC |
| `RIN` | `rin` | RNA |
| `Years of education` | `years_education` | RNA, ATAC |
| `Highest level of education` | `education_level` (string) | RNA, ATAC |

- `donor_covariates.csv`: one row per donor with an RNA object. Each column must hold one value per donor object, otherwise the run fails naming the file and column.
- Coding: an `N+` age (e.g. `90+`) becomes `N` with `age_capped = True`; empty, `nan`, `Unknown` or any other non-numeric value in a numeric column is left missing and listed in the manifest. Nothing is imputed.
- RNA values are used. With `--omic both`, the shared columns are compared with the ATAC object (numeric columns after coding, relative tolerance 1e-9); each disagreement is recorded under `findings`, and the run continues.
- `covariates_manifest.json`: columns and coding rules, objects read per omic, missing values, capped-age count, cross-check summary, findings, software versions and git revision.

### Donor technical covariates

```
sead-sexdiff technical-covariates                 # RNA objects plus the ATAC library metrics
sead-sexdiff technical-covariates --omic rna      # RNA only
sead-sexdiff technical-covariates --force         # recompute cached donors
```

Options: `--cell-type {astrocyte,microglia}` (default `astrocyte`; every directory default below is the cell type's, under its results root, shown here for `astrocyte`), `--omic {rna,both}` (default `both`), `--data-root`, `--pseudobulk-dir` (default `results/pseudobulk`; the `pseudobulk` step's outputs, read for the cross-check), `--out-dir` (default `results/pseudobulk`), `--force` (recompute donors already in the per-donor cache).

Per RNA donor object it selects the cell type's nuclei exactly as the pseudobulk step of the same cell type does, through the shared `select_nuclei` (the cell type's mapping of the release, `Reference` donors excluded, ≥ 200 detected genes in the `UMIs` layer; the pseudobulk manifest must record the same cell type), reads their rows of the raw `UMIs` layer and their `obs` QC columns, and never loads the normalized `X`. From the ATAC objects it reads `obs` columns only. It needs the raw donor objects and the pseudobulk outputs, so it runs on the cluster after `pseudobulk` (`scripts/technical.sbatch`, the same partition and resources as `pseudobulk.sbatch`; a `nohup` run on the login node is also acceptable, since the per-donor read is the pseudobulk's, a few seconds per donor):

```
sbatch scripts/technical.sbatch
```

Gene sets (`TECHNICAL_GENE_SETS`, name patterns on the gene symbols of the pseudobulk's feature set; the manifest lists the matched genes):

| set | patterns | matched in the MTG RNA objects |
|---|---|---|
| `ribosomal` | `^RP[LS]\d+[A-Z]?$`, `^RPLP[0-2]$`, `^RPSA$` (pseudogenes such as `RPL13AP5` and kinases such as `RPS6KA1` do not match) | 81 |
| `mitochondrial` | `^MT-` | 13 |
| `hemoglobin` | `HBB`, `HBA1`, `HBA2` | 3 |

QC columns (`TECHNICAL_OBS_COLUMNS`), summarized over the selected nuclei and written in snake_case (`Fraction mitochondrial UMIs` → `fraction_mitochondrial_umis`):

- RNA, per nucleus: `Fraction mitochondrial UMIs`, `Number of UMIs`, `Genes detected`, `Doublet score`.
- RNA, library level (constant within a library; 80 of the 84 donors span two to five libraries, so the summary is the astrocyte-weighted mean): `GEX_Reads_mapped_confidently_to_intronic_regions`, `GEX_Fraction_of_transcriptomic_reads_in_cells`, `GEX_Median_UMI_counts_per_cell`, `GEX_sequencing_saturation`; the distinct values of `load_name` give `n_libraries_rna`.
- ATAC, library level (usually one library per donor object; 29 of the 86 objects carry two or three distinct values, so the summary is again the astrocyte-weighted mean, with a per-donor `constant` flag in the manifest): `ATAC_TSS_enrichment_score`, `ATAC_Fraction_of_high_quality_fragments_overlapping_peaks`, `ATAC_Fraction_of_high_quality_fragments_in_cells`, `ATAC_Percent_duplicates`, `ATAC_Median_high_quality_fragments_per_cell`.

Per nucleus, each set's fraction is the set's UMIs over the nucleus's total UMIs (all genes of the raw layer). Per donor and set the table holds `<set>_mean` and `<set>_sd` (nucleus-weighted mean and sample SD of the per-nucleus fractions), `<set>_pseudobulk` (sum of set UMIs over sum of total UMIs, the fraction `rna_counts.h5ad` gives) and `<set>_supertype_share` (within-donor η²: the between-supertype sum of squares of the per-nucleus fraction over its total sum of squares, over the supertypes present; empty when fewer than two supertypes hold ≥ 2 nuclei or the fraction has no variance). QC columns get `<column>_mean` and `<column>_sd` (RNA) or `<column>_mean` (ATAC, with a per-donor `constant` flag in the manifest). Nothing is imputed: non-numeric values are dropped from the means and listed in the manifest. The manifest also holds the between-donor share of each fraction and QC column (the donor η² over all nuclei of all donors, from the per-donor counts, means and SDs), which the technical-axis report reads beside the supertype share.

Cross-check ("the same nuclei as the pseudobulk"): for every donor, each set's `<set>_pseudobulk` is recomputed from the donor's row of `rna_counts.h5ad` and must agree to a relative tolerance of 1e-9, and the selected nucleus count must equal `n_nuclei` in the pseudobulk cache (`per_donor/rna/<donor>.npz`, when present) and `n_rna` in `donor_table.csv`; any disagreement stops the run naming the donor and the set or both counts. The manifest records the pseudobulk manifest's sha256. A `Reference` donor is recorded as excluded with no technical row, as in the pseudobulk.

Outputs in `results/pseudobulk/` (so `scripts/sync_pseudobulk.sh` carries them to the workstations):

| file | content |
|---|---|
| `donor_technical.csv` | one row per non-`Reference` donor with an RNA object (84), sorted by `donor_id`: `donor_id`, `n_nuclei_rna`, `n_libraries_rna`, per set `<set>_mean`, `<set>_sd`, `<set>_pseudobulk`, `<set>_supertype_share`, the RNA QC means and SDs, the ATAC means (empty for a donor without an ATAC object, listed in the manifest) |
| `nucleus_technical_summary.csv` | one row per donor and set: `donor_id`, `gene_set`, `n_nuclei`, `mean`, `sd`, `pseudobulk`, `supertype_share`, `n_supertypes` |
| `technical_manifest.json` | gene sets with patterns and matched genes, QC columns per omic with their output names, the definitions, objects read per omic, excluded donors, donors missing an ATAC object, per-donor nucleus counts beside the pseudobulk cache's and `donor_table.csv`'s with the fraction cross-check and the libraries, the pseudobulk manifest sha256, the between-donor shares, the libraries-per-donor distribution, ATAC constancy flags, missing values, findings, software versions, git revision, and the sha256 of `donor_covariates.csv` and `covariates_manifest.json` before and after the run |

The step never writes `donor_covariates.csv` or `covariates_manifest.json`, whose hashes the primary export's and analysis's reference checks bind. Its per-donor cache lives in `results/pseudobulk/per_donor_technical/<omic>/<donor>.json` (the summaries and per-supertype counts and means, no per-nucleus arrays); a rerun skips cached donors unless `--force`, so a stopped run resumes, and a cache made with different gene sets, QC columns or `min_genes` is refused. `scripts/sync_pseudobulk.sh` excludes the cache unless `--with-cache`.

### MOTCO input export

```
sead-sexdiff export                     # both cohorts
sead-sexdiff export --cohort cohort_30 --force
```

Options: `--cell-type {astrocyte,microglia}` (default `astrocyte`; every directory default below is the cell type's, under its results root, shown here for `astrocyte`), `--pseudobulk-dir` (default `results/pseudobulk`), `--out-dir` (default `results/motco_input`), `--cohort {cohort_30,cohort_50,all}` (default `all`), `--data-root` (only to locate the annotation), `--force` (replace existing cohort directories; without it an existing cohort directory is an error).

Reads only the pseudobulk outputs (`rna_counts.h5ad`, `atac_counts.h5ad`, `donor_table.csv`, `manifest.json`, `donor_covariates.csv`), so it runs on any workstation after `scripts/sync_pseudobulk.sh`. It refuses a partial pseudobulk run, missing files and cohorts marked not analysable. On first use it downloads the GENCODE v32 basic annotation (the reference of 10x GRCh38-2020-A / Cell Ranger 6.0) into `<data root>/annotation/gencode.v32.basic.annotation.gtf.gz` and records its sha256; this needs internet once.

Per cohort and omic, in this order, each step computed on the cohort's donors only and never using sex, stage or covariates:

1. Restrict to the cohort's donors (`cohort_30` primary, `cohort_50` sensitivity).
2. Drop chrX, chrY and chrM features. RNA genes get their chromosome from the annotation via `gene_ids` (version ignored; genes absent from it, e.g. those on unplaced scaffolds, are kept and counted as a finding). ATAC peaks get it from their `chrN:start-end` name.
3. Keep features with CPM ≥ 1 in at least `k` donors, `k` = the cohort's smallest sex × stage cell (9 in both cohorts). CPM uses each donor's total over the features left after step 2.
4. TMM factors on the filtered counts (edgeR `calcNormFactors(method = "TMM")`, reimplemented in numpy).
5. log2-CPM with prior count 2 (edgeR `cpm(log = TRUE, prior.count = 2)`: effective library size = filtered total × TMM factor).
6. Keep the 5,000 features with the highest variance (ddof 0) of log2-CPM, ties by name. Columns keep the input feature order.

These are scientific choices; the constants live in `src/sead_sexdiff/data.py` (`EXPORT_EXCLUDED_CHROMOSOMES`, `EXPORT_MIN_CPM`, `EXPORT_PRIOR_COUNT`, `EXPORT_TOP_N`, `EXPORT_COHORTS`, `STAGE_LABELS`, `GENCODE_RELEASE` / `GENCODE_URL`) and are not CLI flags.

Outputs in `results/motco_input/<cohort>/` (each cohort is written to `<cohort>.tmp/` and renamed when complete; all donor-indexed files list the same donors sorted by `donor_id`):

| file | content |
|---|---|
| `rna.csv`, `atac.csv` | donors × kept features, first column `donor_id`, log2-CPM at `%.10g` |
| `metadata.csv` | `donor_id`, `sex`, `stage` (`0`/`1`/`2`), `stage_label` (`Not AD + Low` / `Intermediate` / `High`) |
| `covariates.csv` | `donor_id`, the coded covariates, `n_rna`, `n_atac` (exported, **not applied**) |
| `features_rna.csv`, `features_atac.csv` | every input feature with its chromosome, `status` (`sex_chromosome`, `abundance`, `variance` or `kept`) and log2-CPM variance where computed |
| `manifest.json` | block units (`log2-CPM`, MOTCO transform `identity`), group/stage columns and stage order, thresholds (`min_cpm`, `k` and its cell, `prior_count`, `N`, excluded chromosomes), feature counts after each step, per-donor library sizes and TMM factors, annotation URL and sha256, input manifest hash, covariate missingness, software, git revision, findings |

The CSVs are byte-identical across reruns on the same inputs. They feed MOTCO directly:

```
motco analyze --block rna=rna.csv --block atac=atac.csv --metadata metadata.csv \
  --sample-col donor_id --group-col sex --stage-col stage --stage-order 0,1,2
```

The case-study analysis itself, with covariates, is the `analyze` command below.

`scripts/check_tmm_edger.R` is a one-off cross-check of the TMM factors and log2-CPM against edgeR (in the Nix shell). It reads the files written by `run_export(debug_dir=...)`:

```
uv run python -c "from sead_sexdiff.data import run_export; run_export(cohort='cohort_30', out_dir='/tmp/x', debug_dir='/tmp/dbg')"
Rscript scripts/check_tmm_edger.R /tmp/dbg cohort_30
```

### MOTCO trajectory analysis

```
sead-sexdiff analyze                    # both cohorts
sead-sexdiff analyze --cohort cohort_50 --force
```

Options: `--cell-type {astrocyte,microglia}` (default `astrocyte`; every directory default below is the cell type's, under its results root, shown here for `astrocyte`), `--cohort {cohort_30,cohort_50,all}` (default `all`), `--input-dir` (default `results/motco_input`), `--out-dir` (default `results/motco_analysis`), `--n-jobs` (parallel jobs for MOTCO's CV and RRPP, `-1` = all CPUs; default 1; MOTCO seeds one RRPP stream per worker, so other values change the permutation draws, p-values and null summaries, though not the statistics, rank or latent space; the reported results use 1), `--force` (replace existing cohort directories; without it an existing cohort directory is an error and nothing is written).

Reads only the export (`rna.csv`, `atac.csv`, `metadata.csv`, `covariates.csv`, `manifest.json` per cohort), so it runs on any workstation. It needs MOTCO ≥ 0.8.0, a dependency pinned to the tag `v0.9.0` of [tomszar/MOTCO](https://github.com/tomszar/MOTCO) (0.9.0 adds the latent adjustment used by `latent-sensitivity`; without that setting every output equals 0.8.0's apart from the runtime); an older install stops the run before any file is read. It stops when a file is missing (naming the `export` command), when the four tables list different donors, or when the export manifest's donor or feature counts disagree with the tables.

Each cohort (`cohort_30` primary, `cohort_50` sensitivity) is analysed once with MOTCO's real-data entry point `motco.analysis.analyze_trajectories`, never with `motco plsr` + `motco de`:

- blocks `rna` and `atac` in identity units (log2-CPM as exported); group `sex`; stage `stage` in the explicit order `0, 1, 2`;
- MOTCO's defaults: net-displacement orientation, double-CV rank selection (`cv2_splits=4`, `cv1_splits=3`, `n_repeats=5`, `max_components=20`, `random_state=1203`), no attribution;
- 999 RRPP permutations, seed `20261019`, α = 0.05 with rejection when p ≤ α;
- covariates, additive and in this order: `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`. The last two are `log10(n_rna)` and `log10(n_atac)`, derived here (the export is unchanged). `education_level` is not entered (categorical, one design column per level) and neither is `age_capped` (constant in this release: no capped ages; MOTCO rejects a constant covariate). Missing-covariate policy `error`: a missing value stops the run naming the column and donor, and nothing is imputed. There is no unadjusted run.

These are scientific choices; the constants live in `src/sead_sexdiff/data.py` (`ANALYSIS_COVARIATES`, `ANALYSIS_LOG10_COLUMNS`, `ANALYSIS_EXCLUDED_COVARIATES`, `ANALYSIS_MISSING_COVARIATES`, `ANALYSIS_PERMUTATIONS`, `ANALYSIS_SEED`, `ANALYSIS_ALPHA`, `ANALYSIS_REJECTION_RULE`, `ANALYSIS_MOTCO_MIN_VERSION`) and are not CLI flags.

Outputs in `results/motco_analysis/<cohort>/` (written to `<cohort>.tmp/` and renamed when complete):

| file | content |
|---|---|
| `result.json` | MOTCO's record: `statistics` (`delta`, `angle` in degrees, `shape`), `p_values`, `null_summary`, `latent` (`selected_lv`, `component_selection`, `cv_mean_auroc`), `design` (levels, cell sizes), `covariates` (design columns, reference row, missing policy, dropped donors, coefficients), `measurement`, `software` |
| `latent_scores.csv`, `pls_loadings.csv`, `stage_means.csv` | MOTCO's per-donor PLS scores, feature rotations and latent LS-means per sex × stage |
| `analysis_metadata.csv` | the table the model saw: `donor_id`, `sex`, `stage`, `stage_label` and the six covariates, sorted by `donor_id` |
| `analysis_manifest.json` | sha256 of the five inputs (incl. the export manifest), model declaration, preprocessing summary copied from the export manifest (excluded chromosomes, `min_cpm`, `k`, `prior_count`, `N`, feature counts per step, capped-age count, export findings), covariates entered / derived / excluded with reasons, measurement settings, `cli_equivalent`, software (incl. `motco`, `scikit-learn`), git revision, findings |
| `summary.md` | generated from `result.json`: statistics with p, null mean / sd / q95 and the CV-selected rank in the `angle` row; cell sizes; latent block; covariate block |

At the same `--n-jobs`, the latent scores, loadings, stage means and `analysis_metadata.csv` are byte-identical across reruns and `result.json` differs only in `runtime_seconds`. To reproduce a cohort with MOTCO's own CLI, run the manifest's `cli_equivalent.command` from the repository root with `OUT_DIR` replaced by a new directory; it passes `analysis_metadata.csv` with the six `--covariate` flags and gives the same statistics, p-values and rank.

Interpretation rules, binding from MOTCO's operating studies (`../MOTCO/docs/roadmap.md`, Phase 6 and Next changes item 1; the covariate Type I report §5):

- `delta`: a rejection is evidence that the trajectories differ, not specifically in size. MOTCO's unadjusted Phase 6 study let a non-rejection argue against size differences of about 25 % or more, but the adjusted analysis carries no power claim, so here a non-rejection is not evidence against a size difference.
- `angle` (net displacement): the angle between the sexes' net displacement vectors (stage 2 mean − stage 0 mean), so it reads only the endpoint stages. Report the CV-selected rank beside it. A non-rejection is uninformative at n = 80, not evidence that the directions agree.
- `shape`: a rejection is evidence of a difference but not of which kind; a non-rejection is uninformative.
- Report the `covariates` block beside the results. The adjustment is additive and linear, the latent space is fitted without the covariates (covariate-blind basis), the Type I check covered three covariate design columns (six are entered here), and there is no power claim for the adjusted analysis. The operating characteristics come from a structural analogue (InterSIM blocks, 498 features), not from these data.

### MOTCO rank sensitivity

```
sead-sexdiff rank-sensitivity                      # both cohorts
sead-sexdiff rank-sensitivity --cohort cohort_50 --force
```

Options: `--cohort {cohort_30,cohort_50,all}` (default `all`), `--input-dir` (default `results/motco_input`), `--analysis-dir` (default `results/motco_analysis`; the `analyze` step's output, read as the reference and never written), `--out-dir` (default `results/motco_rank_sensitivity`), `--n-jobs` (parallel jobs for MOTCO's RRPP, `-1` = all CPUs; default 1; other values change the permutation draws at every rank and are recorded as a finding; the reported results use 1), `--force` (replace existing cohort directories under `--out-dir`; without it an existing cohort directory is an error and nothing is written).

MOTCO's documented dimensionality sensitivity analysis for the covariate-adjusted trajectory comparison: per cohort, the analysis is rerun through the same entry point (`analyze_trajectories` with `forced_components = k`) at every rank of a fixed ladder, `ANALYSIS_FORCED_RANKS = (2, 3, 4, 6, 9, 12)`, and everything else equals the `analyze` step: the same export inputs, blocks and units, stage order, six covariates in the same order, missing-covariate policy, net-displacement orientation, 999 RRPP permutations and the same seed `20261019` at every rank, so that every rank draws the same permutation indices and differences between ranks are measurement differences, not draw differences. Cross-validation does not run in a forced-rank run.

The ladder: 2 is `n_stages − 1`, the rank MOTCO's stage-supervised rule is expected to select with three stages and the lowest rank at which `angle` can carry orientation information; 3 is the next rank; 4, 6, 9 and 12 are the steps of MOTCO's own latent-rank ladder, which makes the two sweeps comparable, and 4 = 2(K − 1) is the forced rank with a Type I check at this design (unadjusted). The ladder is a constant, not a flag, and no rank is added, dropped or chosen from the results. Rank 1, the CV-selected rank in both cohorts, cannot be forced (MOTCO accepts `forced_components ≥ 2` only); the CV-selected analysis is the reference row instead, and the manifest says so. A rank larger than `min(n_donors, n_features)` stops the sweep, naming the rank, before anything is written.

The CV-selected analysis is read, not rerun. Before any rank runs, the sweep checks that `<analysis-dir>/<cohort>/analysis_manifest.json` records the same sha256 for each of the five export files as the files now on disk, the same MOTCO version as the one installed and the same seed and permutation count as the sweep uses; a difference stops the sweep naming the file or setting and the `analyze` command. A missing export or analysis directory does the same. Nothing under `results/motco_analysis/` is ever written.

Outputs in `results/motco_rank_sensitivity/<cohort>/` (written to `<cohort>.tmp/` and renamed when every rank and the table are complete):

| file | content |
|---|---|
| `rank_02/` … `rank_12/` | one full analysis directory per forced rank, with the same seven files as `results/motco_analysis/<cohort>/` (`result.json` with `component_selection = forced`, `latent_scores.csv`, `pls_loadings.csv`, `stage_means.csv`, `analysis_metadata.csv`, `analysis_manifest.json` with `measurement.forced_components`, `cv = null` and a `cli_equivalent` carrying `--forced-components k`, and `summary.md` naming the forced rank in the `angle` row) |
| `rank_sensitivity.csv` | the cross-rank table generated from the records: the reference row first (`component_selection = cv`, the CV-selected rank) then one row per forced rank in ascending order, with per statistic (`delta`, `angle`, `shape`) the observed value, the RRPP p-value, whether p ≤ α and the null mean, sd and q95, plus the reference row's `cv_mean_auroc` and each row's directory (relative `rank_<kk>` for the forced rows). Floats are written in their shortest round-trip form and runtimes are left out, so the file is byte-identical across reruns |
| `rank_sensitivity.md` | the same values rendered as one table per statistic plus a run table (CV AUROC, runtime, directory), with a header naming the ladder, seed, permutations, α and the reference rank |
| `rank_sensitivity_manifest.json` | cohort, ladder, seed, permutations, α, `n_jobs`, sha256 of the five export inputs and of the reference `result.json` and `analysis_manifest.json`, the reference's rank and whether it is in the ladder, per-rank directories and runtimes, software (incl. `motco`, `scikit-learn`), git revision, findings (the CV-selected rank not being in the ladder; a forced run at the CV rank differing from the reference; `n_jobs ≠ 1`; a rejection status that differs between a forced rank and the reference) |

Two sweeps on the same inputs, reference, MOTCO release, seed and `n_jobs` give identical statistics, p-values and null summaries at every rank and byte-identical `latent_scores.csv`, `pls_loadings.csv`, `stage_means.csv` and `rank_sensitivity.csv`. Each `rank_<kk>/analysis_manifest.json` carries a `cli_equivalent` that reproduces that rank with `motco analyze --forced-components k` (no CV flags).

Reading rules, binding for the report (`reports/motco-astrocyte-rank-sensitivity-<date>.md`):

- The CV-selected analysis remains the primary result; the sweep is sensitivity to it and never its replacement. No rank is chosen after seeing the results.
- A statistic whose rejection status differs between ranks is rank dependence of the sensitivity analysis, not a result at the rank where it is rejected. No rejection at any rank is no evidence of a difference at any rank, not evidence of agreement in size, direction or shape.
- `angle` and `shape` non-rejections stay uninformative at every rank; a `delta` rejection is a difference, not specifically of size.
- MOTCO's forced-rank operating characteristics are carried as limitations: its Type I evidence at a forced rank covers rank 4 at this design without covariates; forced 2(K − 1) lost `shape` power and leaked `angle` into the magnitude cell in MOTCO's orientation probe; ranks above `n_stages − 1` widened the `angle` null in MOTCO's rank ladder. The covariate-adjusted analysis has no Type I check at any forced rank.

### MOTCO preprocessing sensitivity

```
sead-sexdiff preprocessing-sensitivity                       # every registered variant, both cohorts (resumable)
sead-sexdiff preprocessing-sensitivity --variant residualized --cohort cohort_50 --force
sead-sexdiff preprocessing-sensitivity --variant trend_normalized          # one variant, both cohorts
```

Options: `--variant {keep_chrx,exclude_families,rna_trend_ratio,atac_mean_accessibility,residualized,revised_candidate,trend_normalized,revised_candidate_trend,all}` (default `all`; one value per call), `--cohort {cohort_30,cohort_50,all}` (default `all`), `--pseudobulk-dir` (default `results/pseudobulk`; the export's inputs and, for `exclude_families`, `residualized`, `revised_candidate`, `trend_normalized` and `revised_candidate_trend`, the `technical-covariates` step's `donor_technical.csv` and `technical_manifest.json`), `--input-dir` (default `results/motco_input`; the `export` step's output, read as the primary export and never written), `--analysis-dir` (default `results/motco_analysis`; the `analyze` step's output, read as the primary analysis and never written), `--out-dir` (default `results/motco_sensitivity`), `--data-root` (only to locate the annotation), `--n-jobs` (parallel jobs for MOTCO's CV and RRPP, `-1` = all CPUs; default 1; other values change the permutation draws and are recorded as a finding; the reported results use 1), `--force` (replace a cohort's existing variant export, analysis and comparison directories under `--out-dir`; without it `--variant all` skips, with a message, every variant × cohort whose comparison directory exists, so `all` is resumable, while a named variant with existing output, or a partial tree, is an error and nothing is written).

Measures how the covariate-adjusted trajectory comparison depends on a preprocessing choice of the export: per variant and cohort, the export is rerun with the variant's predeclared setting(s) changed and the analysis is rerun unchanged on the variant blocks, beside the primary analysis. The variants are a registry constant, `PREPROCESSING_VARIANTS` in `src/sead_sexdiff/data.py`, not a flag; no variant is added, dropped or altered from the results, the five feature-content entries were added with their adoption rule written before any run, and the two trend-normalization entries with the precedence rule between them. Every variant is the primary pipeline (restrict to the cohort → remove excluded chromosomes and excluded gene sets → CPM ≥ 1 in ≥ `k` donors → TMM → log2-CPM → trend normalization when declared → residualization when declared → top 5,000 per omic by the selection rule, ties by name) with the primary's `k`, `min_cpm`, prior count, `N`, TMM settings and GENCODE v32 basic, differing only in the setting(s) listed:

- `keep_chrx`: excluded chromosomes `chrY` and `chrM` instead of `chrX`, `chrY`, `chrM`, so chrX features are kept. Mechanically this is the primary pipeline with chrX present, not the primary blocks with chrX appended: step 3's library sizes include the chrX counts, so every feature's CPM shifts slightly and the abundance filter can change at its margin; TMM factors and log2-CPM are recomputed on the new filtered set; chrX features enter the variance ranking and compete for the 5,000 slots like any other feature. The primary export removed chrX by design because sex is the group and X-dosage features (`XIST` in females, the escapees) would separate the sexes trivially; the variant measures what that removal changes.
- `exclude_families`: before the abundance filter, every RNA gene matched by the `ribosomal`, `mitochondrial` and `hemoglobin` sets of `TECHNICAL_GENE_SETS` (the `technical-covariates` step's constant: `^RP[LS]\d+[A-Z]?$` / `^RPLP[0-2]$` / `RPSA`, so pseudogenes and kinases are not matched; the 13 `MT-` genes; `HBB`, `HBA1`, `HBA2`) is removed with status `excluded_gene_set` in `features_rna.csv`; a matched gene already on an excluded chromosome (chrM for the `MT-` genes, chrX for `RPL10`, `RPL36A`, `RPL39` and `RPS4X`) keeps `sex_chromosome`. The step-3 library sizes exclude the removed genes, so every other gene's CPM shifts and the abundance filter, TMM, log2-CPM and ranking are recomputed without them. ATAC names match no set, so the ATAC block equals the primary's. The manifest records the matched genes per set (`removed_by_gene_set`) and `feature_counts` gains `after_gene_set`.
- `rna_trend_ratio`: at the selection step RNA features are ranked by `variance / trend` instead of raw variance. The trend is the median log2-CPM variance of 20 equal-count bins of mean log2-CPM (`TREND_N_BINS`; `np.array_split` of the abundance-passing genes sorted by mean, so bin sizes differ by at most one), each placed at its median mean and linearly interpolated (`np.interp`, flat beyond the end bins); dropped genes carry status `selection`; `features_rna.csv` gains `mean`, `trend` and `ratio`, and the manifest's `selection_rule` records the rule per omic and the 20 bins (`trend_bins`: edges, sizes, median mean and median variance). Raw log-CPM variance favours low-count features; the ratio ranks each gene against genes of its own abundance. ATAC is ranked by variance as in the primary.
- `atac_mean_accessibility`: ATAC peaks are ranked by mean log2-CPM (total accessibility), ties by name, dropped peaks `selection`, `features_atac.csv` gains `mean`; RNA is ranked by variance as in the primary.
- `residualized`: after log2-CPM and before the selection step, every abundance-passing feature of both blocks is replaced by its mean plus the residual of its OLS regression on an intercept and `RESIDUALIZATION_COVARIATES`: `log10_lib_rna` and `log10_lib_atac` (log10 of each donor's counts over the features kept at step 2 of this export, the step-3 library sizes recorded as `lib_size_filter` in the manifest's `tmm` block) and the gene-set pseudobulk fractions `ribosomal_pseudobulk` and `mitochondrial_pseudobulk` of `donor_technical.csv` (the UMI-weighted fractions over the donor's astrocyte nuclei, untransformed). `hemoglobin_pseudobulk` is excluded under the technical-axis report's cell-state rule (`RESIDUALIZATION_CELL_STATE_EXCLUSIONS`, with the report's section as the reason: its median within-donor supertype share exceeds its between-donor share, both negligible). Units stay log2-CPM, means are unchanged, the exported values are orthogonal to every design column, the selection ranks the residual variance (so `variance` in `features_<omic>.csv` is the residual variance), the design is written as `residualization.csv` beside the blocks (the cohort's donors in export order; not an analysis input: `analyze` reads the five files it always read), `features_<omic>.csv` gains `mean` and `r_squared`, and the manifest's `residualization` block records the covariates, the exclusions, the design's sha256, the degrees of freedom and the R² quantiles per omic. No sex, stage or analysis covariate enters the design. Nuclei counts, RIN and PMI stay out because they are analysis covariates already in the RRPP model; the TMM factors stay out because they are normalization by-products and serve as check variables instead.
- `revised_candidate`: the four settings together in the primary's order (families removed → abundance → TMM → log2-CPM → residualized → ranked by trend ratio for RNA and mean accessibility for ATAC); its manifest records all four, and for this variant the step evaluates the adoption rule below.
- `trend_normalized`: after log2-CPM and before the selection step, each donor's log2-CPM over the abundance-passing features of each block (the step-3 set TMM is estimated on, not the selected set) is corrected by limma's fast cyclic loess (`normalizeCyclicLoess(method = "fast")`): at each of 3 iterations (`LOESS_ITERATIONS`) the reference is the per-feature mean over donors (A), and each donor's deviation from it (M) is fitted by a lowess curve of M on A and subtracted. The lowess is statsmodels' (Cleveland's algorithm, as R's `lowess`, which limma's `loessFit` calls): `frac` = the span, 3 robustifying iterations (`LOESS_ROBUSTNESS_ITERATIONS`; `loessFit`'s default), `delta` = 1 % of the range of A (`LOESS_DELTA_FRACTION`, R's default). The span is limma's adaptive default from the block's number of features, `chooseLowessSpan(n)` = min(0.3 + 0.7 · (50 / n)^(1/3), 1) (`LOESS_SPAN_RULE`; about 0.39 for 20,000 genes and 0.35 for 129,000 peaks), so RNA and ATAC get their own spans. One TMM factor per donor cannot remove an abundance-dependent tilt of M; the cyclic loess removes one smooth curve per donor and reads nothing but the block's values (no label, covariate or technical variable). Per donor and omic the step records the tilt before and after (with A and M taken against the per-feature mean of the values at that stage: the OLS slope of M on A and the median M of the lowest and highest of 20 equal-count bins of A, `NORMALIZATION_TILT_BINS`) and the largest absolute correction, written as `normalization_trend.csv` beside the blocks; the manifest records `normalization = {rna: cyclic_loess, atac: cyclic_loess}` and a `normalization_fit` block (per omic the span, iterations, number of features, `delta` per iteration and the quantiles over donors of the tilt statistics and the correction; the reference `feature_mean`, the span rule and the trend file's sha256), and `preprocessing_order` gains the step after log2-CPM. Every other setting is the primary's. The step evaluates the adoption rule below for this variant.
- `revised_candidate_trend`: `revised_candidate`'s four settings plus the cyclic loess of both omics, in the primary's order (families removed → abundance → TMM → log2-CPM → cyclic loess → residualized → ranked by trend ratio for RNA and mean accessibility for ATAC): the residualization and the selection see the normalized values. Its manifest records all five settings and both `residualization.csv` and `normalization_trend.csv` are written. The step evaluates the adoption rule below for this variant.

The primary export records `normalization = {rna: none, atac: none}` and is otherwise byte-identical to what it was before the normalization setting existed. The numpy cyclic loess is cross-checked against limma by `scripts/check_cyclic_loess_limma.R <debug_dir> [cohort] [omic ...]` (inside the Nix shell): it reads `<cohort>_<omic>_logcpm.csv`, `<cohort>_<omic>_normalized.csv` and `<cohort>_<omic>_loess.csv` (span, iterations, number of features) written by a scratch `_export_cohort(..., normalization="cyclic_loess", debug_dir=...)` call (`_write_debug`), checks that the recorded span is limma's adaptive span, runs `limma::normalizeCyclicLoess(t(logcpm), span = <recorded>, adaptive.span = FALSE, iterations = 3, method = "fast")` and exits non-zero when the largest absolute difference exceeds 1e-6 log2-CPM (the remaining differences, about 1e-8, come from statsmodels' lowess leaving out R's 0.001 / 0.999 weight cut-offs).

The variant analysis is one `analyze_trajectories` call with the `analyze` step's settings: the same blocks and units, stage order, six covariates in the same order, missing-covariate policy `error`, net-displacement orientation, CV rank selection at MOTCO's defaults, 999 RRPP permutations and the seed `20261019`, so the variant draws the same permutation indices as the primary and the difference between the two rows is a data difference. There is no per-variant analysis setting: the residualized variants keep the six analysis covariates too, so the library-size proxies are adjusted twice (log10 library sizes in the design, log10 nuclei counts in the RRPP model), a stated limitation.

The primary export and analysis are read, never rerun. Before anything is written, for every requested variant and cohort, the step checks that `results/motco_input/<cohort>/manifest.json` records the same sha256 for the pseudobulk manifest, the covariates manifest and the annotation as the files now on disk, the same `min_cpm`, `prior_count`, `N`, TMM settings and excluded chromosomes as the constants, and the same `k`, cells and donor count as `donor_table.csv` gives now (a difference stops the step naming the file or setting and the `export` command); and that `results/motco_analysis/<cohort>/analysis_manifest.json` records the same sha256 for the five export files as the files now on disk, the installed MOTCO version, the seed and permutation count the step uses and CV rank selection (the rank sweep's check; a difference names the file or setting and the `analyze` command). For `exclude_families`, `residualized`, `revised_candidate`, `trend_normalized` and `revised_candidate_trend` it also reads `results/pseudobulk/donor_technical.csv` and `technical_manifest.json` and checks that the technical manifest records the sha256 of the pseudobulk `manifest.json` now on disk, that its gene sets equal `TECHNICAL_GENE_SETS`, and that every cohort donor has a non-missing value for every residualization covariate and, for the adoption candidates, every adoption-check variable (a mismatch stops the step naming both hashes; a missing value names the donor and the column; nothing is imputed); `keep_chrx` and the two selection-rule variants run without the technical table. A missing pseudobulk file, technical file, export or analysis directory stops the step naming it and the command that produces it. The variant export's `metadata.csv` and `covariates.csv` must then be byte-identical to the primary's; a difference is an error. Nothing under `results/motco_input/`, `results/motco_analysis/` or `results/pseudobulk/` is ever written.

Adoption rule (`ADOPTION_RULE`, a constant written in the change proposal before any run): for `revised_candidate`, `trend_normalized` and `revised_candidate_trend` the step computes, per cohort, the Spearman correlation of every latent component of the candidate's CV-selected analysis (`latent_scores.csv`) with each of `ADOPTION_CHECK_VARIABLES`, technical variables that are not regressors of the residualization: the primary export's TMM factors `tmm_rna` and `tmm_atac` (from its `manifest.json`), the technical table's `fraction_mitochondrial_umis_mean`, `doublet_score_mean`, `gex_reads_mapped_confidently_to_intronic_regions_mean` and `gex_fraction_of_transcriptomic_reads_in_cells_mean`, and `rin` and `pmi` from the primary export's `covariates.csv`. A cell with `|rho| ≥ 0.3` (`DIAGNOSTIC_CORRELATION_THRESHOLD`) is flagged, and `comparison/<cohort>/adoption_check.csv` holds one row per component × variable (`component`, `variable`, `source`, `n`, `spearman`, `flag`) and a final `all` row whose `holds` is true only when no cell is flagged. The rule holds when it holds in both cohorts. It reads no p-value, rejection status or CV AUROC. If it holds, a later change may propose the candidate as a revised primary with the current primary kept as the reference; if it does not, the variants stay sensitivity analyses. A check variable missing for a cohort donor is an error, not a skipped cell.

Precedence rule (`ADOPTION_PRECEDENCE`, written in the trend-normalization proposal before any run): between `trend_normalized` and `revised_candidate_trend` (`ADOPTION_PRECEDENCE_CANDIDATES`), if `trend_normalized` meets the adoption rule in both cohorts it is the candidate a later revised-primary change proposes, whatever `revised_candidate_trend` gives, because it changes the least (no gene removed, no degrees of freedom spent, no technical covariate read); if only `revised_candidate_trend` meets it in both cohorts, that one is; if neither does, the variants stay sensitivity analyses. Only the `holds` values of the four adoption tables are read (never a p-value, rank or CV AUROC). After each cohort of either candidate the step evaluates the rule from the tables under `--out-dir` and writes it as `adoption_precedence` (`candidate`, `reason`, the `holds` table, the flagged cells per variant and cohort) into every existing comparison manifest of the two candidates; until both have run in both cohorts it records `candidate = null` with reason `incomplete` and the missing pairs.

Outputs under `results/motco_sensitivity/<variant>/`, a layout that mirrors the primary pipeline so that `rank-sensitivity --input-dir results/motco_sensitivity/<variant>/motco_input --analysis-dir results/motco_sensitivity/<variant>/motco_analysis --out-dir …` and `axis-diagnostics` with the same two flags can read a variant without code changes (each cohort directory is written as `<cohort>.tmp/` and renamed, in the order export → analysis → comparison; a cohort is complete when its comparison directory exists):

| path | content |
|---|---|
| `motco_input/<cohort>/` | the variant export: the same seven files as `results/motco_input/<cohort>/`, with the variant's settings in `thresholds.excluded_chromosomes`, `excluded_gene_sets` (plus `removed_by_gene_set`), `selection_rule` (plus `trend_bins`), `normalization` (plus `normalization_fit`) and `residualization`, and a `variant` block (name, description, settings, the primary export directory) in its `manifest.json`; `features_<omic>.csv` carries the statuses of the variant (`sex_chromosome`, `excluded_gene_set`, `abundance`, `variance` or `selection`, `kept`) and the extra columns its rule produced; for a residualizing variant also `residualization.csv`, for a trend-normalizing variant also `normalization_trend.csv` (per donor and normalized omic: `slope_before`, `slope_after`, `low_bin_before`, `low_bin_after`, `high_bin_before`, `high_bin_after`, `max_abs_correction`) |
| `motco_analysis/<cohort>/` | the variant analysis: the same seven files as `results/motco_analysis/<cohort>/` (`result.json` with `component_selection = cv`, `permutations = 999`, the seed and the six covariates; `analysis_manifest.json` whose `preprocessing` block carries the variant's thresholds, settings and feature counts and whose `cli_equivalent` reproduces the run with `motco analyze` on the variant blocks) |
| `comparison/<cohort>/sensitivity.csv` | two rows generated from the records, the primary row first (`variant = primary`, the primary's CV-selected rank, the primary analysis directory) then the variant row (its name, its CV-selected rank, `motco_analysis/<cohort>` relative to the variant root); per statistic (`delta`, `angle`, `shape`) the observed value, the RRPP p-value, whether p ≤ α and the null mean, sd and q95, plus `cv_mean_auroc`. Floats in their shortest round-trip form, no runtime column, so the file is byte-identical across reruns and its columns are those of the archived `keep_chrx` tables |
| `comparison/<cohort>/sensitivity.md` | the same values as one table per statistic; the feature accounting per omic and export (counts after chromosome removal, gene-set removal when sets were excluded, abundance and selection; kept features shared with the other export and kept only by one of them) and, when the variant calls for them, the chrX block (kept chrX features and the status of `XIST`, only when the excluded chromosomes differ), the genes removed per set, the selection rule per omic with the quantiles (5, 25, 50, 75, 95) of the kept features' mean log2-CPM in both exports, the normalization block (rule per omic in both exports, features, span, iterations, the largest absolute correction, the quantiles over donors of the tilt statistics before and after and of the correction, and the limitations), the residualization covariates with the R² quantiles and the limitations, and the adoption-check table with the cohort's outcome (naming `ADOPTION_PRECEDENCE` for its two candidates); a run table (rank selection, CV AUROC, runtime, directory); the findings |
| `comparison/<cohort>/adoption_check.csv` | `revised_candidate`, `trend_normalized` and `revised_candidate_trend` only: the adoption-check rows and the `holds` row described above |
| `comparison/<cohort>/sensitivity_manifest.json` | the variant name and settings, cohort, seed, permutations, α, `n_jobs`; sha256 of the pseudobulk manifest, the covariates manifest, the annotation and (when read) the technical table and its manifest, of the five primary export files, of the primary `result.json` and `analysis_manifest.json`, of the five variant export files (plus `residualization.csv` and `normalization_trend.csv`) and of the variant's two analysis files; both ranks, CV AUROCs and runtimes; the variant's normalization rules and spans; the feature accounting (with the kept chrX names, the per-set counts, the kept-mean and R² quantiles, the normalization accounting as applicable); the adoption check (rule, threshold, variables with sources, `holds`, flagged cells, file hash) for the candidates and, for `trend_normalized` and `revised_candidate_trend`, `adoption_precedence`; software (incl. `motco`, `scikit-learn`), git revision, findings (the CV-selected ranks differ; a statistic's rejection status differs, recorded as preprocessing dependence; the adoption check's outcome; `n_jobs ≠ 1`; the variant export's own findings) |

Two runs of a variant on the same pseudobulk outputs, technical table, annotation, primary export, primary analysis, MOTCO release, seed and `n_jobs` give identical statistics, p-values, rank and null summaries and byte-identical variant export CSVs (incl. `residualization.csv` and `normalization_trend.csv`), `latent_scores.csv`, `pls_loadings.csv`, `stage_means.csv`, `analysis_metadata.csv`, `sensitivity.csv` and `adoption_check.csv`.

Reading rules, binding for the reports (`reports/motco-astrocyte-chrx-sensitivity-<date>.md`, `reports/motco-astrocyte-feature-selection-<date>.md`, `reports/motco-astrocyte-trend-normalization-<date>.md`):

- The primary analysis remains the primary result; every variant is sensitivity to it and never its replacement, and no variant is chosen after seeing the results. `cohort_30` is primary and `cohort_50` sensitivity, as in the primary report. Whether a candidate becomes a revised primary is decided by `ADOPTION_RULE` (and, between `trend_normalized` and `revised_candidate_trend`, `ADOPTION_PRECEDENCE`) alone and implemented by a later change with the current primary kept as the reference.
- A statistic whose rejection status differs between a variant and the primary is preprocessing dependence of the sensitivity analysis, not a result of the variant. No rejection in either is no evidence of a difference under the variant either, not evidence of agreement in size, direction or shape.
- `angle` and `shape` non-rejections stay uninformative in every variant; a `delta` rejection is a difference, not specifically of size. A CV-selected rank that differs is reported beside `angle`, not read as better or worse.
- MOTCO's operating characteristics were measured on neither chrX-kept, family-excluded, trend- or mean-selected, residualized nor trend-normalized blocks, so no variant carries a Type I or power claim of its own.
- Residualization limitations (`residualized`, `revised_candidate`): the report names the covariates; the fit spends five degrees of freedom per feature (intercept plus four covariates) before the permutation null, which does not see them; it removes whatever part of the stage and sex signal the technical set carries; the analysis covariates were left unchanged, so the library-size proxies are adjusted twice; and there is no Type I check for the variant.
- Trend-normalization limitations (`trend_normalized`, `revised_candidate_trend`): the cyclic loess assumes that most features do not change at any abundance level, a stronger assumption than TMM's, so a stage or sex effect shaped like a tilt would be removed with the artefact; the curve is least determined at the rare end, where the selection draws most features (the per-donor tilt before and after is reported); the lncRNA-versus-coding offset at matched abundance is not corrected; and there is no Type I check for the variant. `revised_candidate_trend` also carries the residualization limitations, with three adjustments in sequence (loess, residualization, the six analysis covariates).
- chrX limitations (`keep_chrx`): with chrX kept, features whose expression or accessibility differs between the sexes by X dosage enter the blocks (the report names them, `XIST` first, and says whether the variant kept them); a constant sex offset shifts the two sexes' centroids but not the within-sex displacement between stages, and the PLS basis is stage-supervised, so such features enter the latent space only through their covariance with stage.
- The single-setting variants and the candidate share every other setting, so the report may attribute the change of the first latent axis's technical correlations (an `axis-diagnostics` run on each variant root) to each setting; such an attribution is a measurement-space finding, not a reinterpretation of any statistic.

### MOTCO axis diagnostics

```
sead-sexdiff axis-diagnostics                      # both cohorts
sead-sexdiff axis-diagnostics --cohort cohort_30 --force
```

Options: `--cell-type {astrocyte,microglia}` (default `astrocyte`; every directory default below is the cell type's, under its results root, shown here for `astrocyte`), `--cohort {cohort_30,cohort_50,all}` (default `all`), `--input-dir` (default `results/motco_input`), `--analysis-dir` (default `results/motco_analysis`; the `analyze` step's output, read as the primary analysis and never written), `--rank-dir` (default `results/motco_rank_sensitivity`; the `rank-sensitivity` step's output, whose `rank_*/` analyses are read when present), `--pseudobulk-dir` (default `results/pseudobulk`; holds `donor_technical.csv`), `--out-dir` (default `results/motco_diagnostics`), `--force` (replace existing cohort directories; without it an existing cohort directory is an error and nothing is written).

Measures, from kept outputs only, what the MOTCO latent axes and the exported blocks' principal axes correlate with, so that preprocessing decisions rest on recorded diagnostics rather than on scratch checks. Per cohort it reads the primary export, the primary analysis, `results/pseudobulk/donor_technical.csv` (synced from the cluster after `technical-covariates`; `technical_manifest.json` is read beside it for the pseudobulk-hash check) and every `results/motco_rank_sensitivity/<cohort>/rank_*/` analysis found (an absent sweep is a finding, not an error), and nothing else; it runs on any workstation in seconds. A missing input stops the run naming the file and the command that produces it. Before writing, every analysis manifest's recorded input hashes are compared with the sha256 of the export files now on disk; a difference stops the run naming the file and both hashes. Nothing under the input directories is written.

Variables (`DIAGNOSTIC_VARIABLES`, a constant resolved at run time and recorded in the manifest with each variable's source; never chosen or extended from results): every numeric column of `donor_technical.csv` (gene-set fractions, supertype shares, nucleus and library QC means and SDs, nuclei and library counts); `tmm_<omic>` and `log10_lib_<omic>` from the export manifest's `tmm` block; the six analysis covariates as the analysis derives them (`log10_n_rna`, `log10_n_atac` on the log10 scale); `stage` (0, 1, 2) and `male` (1 for Male).

Statistics (`DIAGNOSTIC_CORRELATION_THRESHOLD` = 0.3 is a reporting flag, not an adoption rule; no p-value is attached to any correlation):

- For every latent component of every analysis read (`latent_scores.csv`): the Spearman correlation with each variable over the donors that have it (`n` recorded), the partial Spearman correlation given the six analysis covariates (the Spearman correlation of the OLS residuals of both variables on an intercept and the covariates; n/a for a covariate itself) and the component's score SD. Rows with |rho| ≥ 0.3 (raw or partial) are flagged.
- The same for the first three principal components of each z-scored block (population SD; SVD scores with the sign fixed so the largest loading is positive), with each component's variance share.
- From each analysis's `stage_means.csv`, the share of the squared length of each sex's net displacement (last stage minus first), of the per-stage sex-difference vectors (squares summed over stages) and of the stage-centred sex differences that lies along each component (at rank 1 every share is 1 and the signed displacements per sex are recorded).
- From the export's `features_<omic>.csv` and block, per omic: the 5/25/50/75/95 % quantiles of the kept features' mean log2-CPM, the share of kept features with mean below 0 and below 1, the share of the abundance-passing variance held by the kept features, the variance at the cut (kept minimum, dropped median) and the share of kept features whose |Spearman| with the primary `pls_0` is ≥ 0.3.
- Findings: whether `pls_0` is identical (`np.allclose`) across the primary and every forced rank (PLS components nest, so the sweep inherits the primary's first axis); the variables flagged against the primary `pls_0`; variables with missing donors (correlations are pairwise-complete); whether the technical table's pseudobulk manifest hash equals the export's.

Outputs in `results/motco_diagnostics/<cohort>/` (written to `<cohort>.tmp/` and renamed; the CSVs are byte-identical across reruns on the same inputs):

| file | content |
|---|---|
| `axis_correlations.csv` | one row per analysis (`primary`, then `rank_02`, …), component and variable: rank, component SD, source, `n`, `spearman`, `n_partial`, `partial_spearman`, `flag_raw`, `flag_partial`; variables sorted by name |
| `block_pcs.csv` | the same for `rna` and `atac` `pc1`–`pc3`, with `variance_share` |
| `geometry_shares.csv` | per analysis, quantity (`net_displacement_<sex>`, `sex_difference_by_stage`, `sex_difference_centred`) and component: the signed value (single vectors only), the squared length along the component, the total and the share |
| `feature_selection.csv` | one row per omic with the accounting above |
| `diagnostics.md` | the flagged correlations first (latent components, then block PCs), the primary analysis's full variables × components table, a per-component summary of each forced-rank analysis (SD, flagged count, top variable, geometry shares), the block PC tables, the geometry shares, the feature-selection table, the analyses read and the findings |
| `diagnostics_manifest.json` | every input path and sha256 (export files, each analysis's four files, the technical table and manifest with the pseudobulk-hash check), the variable list with sources and missing counts, the threshold, the definitions, the block PC variance shares, the nestedness check, flag counts, the output hashes, software, git revision, findings |

Reading rules, binding for the report (`reports/motco-astrocyte-technical-axis-<date>.md`):

- An axis correlated with a technical variable is a finding about the measurement space of the primary analysis, not a reinterpretation of its statistics or p-values. The primary analysis stays primary and its non-rejections keep their reading rules (`angle` and `shape` non-rejections are uninformative; a `delta` rejection is a difference, not specifically of size).
- The correlations are descriptive at n = 80 and carry no p-values; the flag is a reporting threshold, and the report claims no inference from them.
- A gene-set fraction whose between-donor share dominates and whose within-donor supertype share is small in most donors is read as a preparation property and a residualization candidate for the next change; one whose within-donor supertype share exceeds its between-donor share is read as cell state and not a candidate. Both shares are the technical step's (`technical_manifest.json` `between_donor_share`, `donor_technical.csv` `<set>_supertype_share`), read beside the nucleus QC columns, never alone.
- No preprocessing variant is chosen in the report. The diagnostic reads stage and sex, so it is a diagnostic and never a filter: its outputs are read by no `export`, `analyze`, `rank-sensitivity` or `preprocessing-sensitivity` run, and preprocessing variants are predeclared in a change proposal, never derived from a diagnostic value.

### MOTCO latent-adjustment sensitivity

```
sead-sexdiff latent-sensitivity                      # both cohorts, astrocytes
sead-sexdiff latent-sensitivity --cohort cohort_50 --force
sead-sexdiff latent-sensitivity --cell-type microglia   # the same arms on the microglia blocks
```

Options: `--cell-type {astrocyte,microglia}` (default `astrocyte`; sets the default directories below under the cell type's results root, `results/` or `results/microglia/`; the registry of arms, the covariate set, the conditioning, the mediator check and every analysis setting are the same for every cell type), `--cohort {cohort_30,cohort_50,all}` (default `all`), `--input-dir` (default `<root>/motco_input`; the `export` step's output, read as the primary export and never written), `--analysis-dir` (default `<root>/motco_analysis`; the `analyze` step's output, read as the primary analysis and never written), `--pseudobulk-dir` (default `<root>/pseudobulk`; holds `manifest.json`, `donor_technical.csv` and `technical_manifest.json`, read and never written), `--out-dir` (default `<root>/motco_latent_sensitivity`; may not lie inside an input directory), `--n-jobs` (parallel jobs for MOTCO's CV and RRPP, `-1` = all CPUs; default 1; other values change the permutation draws and are recorded as a finding; the reported results use 1), `--force` (replace a cohort's existing arm and comparison directories; without it any of them existing is an error and nothing is written).

Measures how the covariate-adjusted trajectory comparison depends on the covariate set and on how the latent basis is chosen. In both cohorts the primary analysis's first latent axis is a donor-level library-quality axis (`reports/motco-astrocyte-technical-axis-2026-10-19.md`); the analysis covariates adjust the stage means within the basis but never touch which basis PLS-DA chooses. MOTCO 0.9.0's `covariate_latent_adjustment` fits the basis on features residualized on the declared covariates conditional on stage, and MOTCO's binding conditions allow it with `stage` conditioning only, only as a sensitivity analysis beside the covariate-adjusted analysis, and with the covariate set justified by the case study (`../MOTCO/docs/roadmap.md` Phase 6; `docs/api/analysis.md` "Latent adjustment"). Per cohort the step runs the two arms of the registry `LATENT_SENSITIVITY_ARMS` (a constant, not a flag) on the unchanged primary export:

- `extended_covariates`: `LATENT_SENSITIVITY_COVARIATES`, the six analysis covariates followed by `ribosomal_pseudobulk` and `mitochondrial_pseudobulk` (`LATENT_SENSITIVITY_ADDED_COVARIATES`), entered additively in that order, with the covariate-blind basis (no latent adjustment);
- `latent_adjusted`: the same eight covariates with `covariate_latent_adjustment = "stage"` (`LATENT_SENSITIVITY_CONDITIONING`); `group_stage` conditioning is not run.

The two fractions come from `donor_technical.csv` on its raw fraction scale (the scale of the `residualized` variant), joined to the export by `donor_id`; they are the fractions the technical-axis report classified as residualization candidates (between-donor share dominant, within-donor supertype share small), and an import-time check ties them to `RESIDUALIZATION_COVARIATES` minus the export's `log10_lib_<omic>`. `hemoglobin_pseudobulk` stays out as negligible, and the log10 library sizes stay out because the nuclei counts already enter and MOTCO's Type I check used eight covariates. The set is a constant, never chosen from results. Everything else equals the `analyze` step: blocks and units, stage order, net-displacement orientation, CV rank selection at MOTCO's defaults, missing-covariate policy `error` (nothing imputed), 999 RRPP permutations and the seed `20261019`, so every row draws the same permutation indices. Primary → `extended_covariates` is the effect of the larger covariate set; `extended_covariates` → `latent_adjusted` is the effect of the basis alone, since MOTCO guarantees that the trajectory model, LS-means and RRPP of the two are identical and only the basis and the CV-selected rank change (the step checks that both arms carry the same covariate block).

The primary analysis is the reference row, read and never rerun. Before anything is written, for every requested cohort, the step checks that the installed MOTCO is at least `LATENT_MOTCO_MIN_VERSION` (0.9.0; checked before any input is read), that `results/motco_analysis/<cohort>/analysis_manifest.json` records the sha256 of the five export files now on disk, the installed MOTCO version, the step's seed and permutation count, CV rank selection and the arms' orientation and CV settings (a difference names the file or setting and the `analyze` command: after a MOTCO pin move, rerun `analyze` first), and that the technical manifest records the sha256 of the pseudobulk `manifest.json` now on disk and gene sets equal to `TECHNICAL_GENE_SETS`, with a finite value of both fractions for every cohort donor (a missing value names the donor and the column). A missing export, analysis directory or technical file stops the step naming it and the command that produces it (with `--cell-type` when it is not the default). Every manifest read (the export manifest, the analysis manifest, the technical manifest and the pseudobulk manifest) must record the requested cell type: another cell type's manifest stops the step naming the file and both cell types with nothing written; a manifest written before the registry has no `cell_type`, reads as `astrocyte` and is listed as a finding in the comparison manifest (not in `latent_sensitivity.md`, whose findings are the comparison's). Nothing under `<root>/motco_input/`, `<root>/motco_analysis/` or `<root>/pseudobulk/` is ever written, and the reference files are re-hashed after the run.

Mediator check (predeclared, written before any arm runs, never a gate): per cohort `mediator_check.csv` gives each added fraction's Spearman correlation with stage (ordinal 0, 1, 2) over all cohort donors and within each sex (`subset` = `all`, `F`, `M`), with `n`; a cell with |ρ| ≥ 0.3 (`DIAGNOSTIC_CORRELATION_THRESHOLD`) is flagged as a mediator candidate. No p-value is computed, and the correlation is not partial on the analysis covariates (the question is whether the fraction moves with stage at all). A flag removes no covariate, changes no arm and stops no run, because no simulation can tell a confounder from a mediator; it binds the report to say that the `latent_adjusted` basis may have removed stage signal carried by that fraction.

Outputs under `<root>/motco_latent_sensitivity/` (`results/motco_latent_sensitivity/` for astrocytes, `results/microglia/motco_latent_sensitivity/` for microglia; per cohort, the three directories are written as `<cohort>.tmp/` and renamed once the comparison is complete, the arms first and the comparison last; on an error the cohort's temporaries are removed):

| path | content |
|---|---|
| `<arm>/motco_analysis/<cohort>/` | one complete analysis directory per arm, the seven files `analyze` writes: `result.json` with the eight-covariate `covariates` block and, for `latent_adjusted`, MOTCO's `latent_adjustment` block (`conditioning`, `covariates`, `design_columns`, `n_samples`, `residual_df`, `cv_min_training_df`, `feature_scale`, `feature_r2` quantiles, mean and fraction above 0.5) and `params.covariate_latent_adjustment`; `analysis_metadata.csv` with the eight covariate columns; `analysis_manifest.json` whose `inputs.files` holds the five export files (so `axis-diagnostics` can read the arm against the primary export), with the technical table, its manifest and the pseudobulk manifest under `inputs.technical`, the fractions' source under `covariates.sources` and, for `latent_adjusted`, `measurement.latent_adjustment = stage` and a `cli_equivalent` carrying `--latent-adjustment stage`; `summary.md` with a "Latent adjustment" section for the adjusted arm; the manifest records `cell_type`. The primary's manifests carry neither latent-adjustment key, as MOTCO leaves the unset setting out of its record |
| `comparison/<cohort>/latent_sensitivity.csv` | three rows generated from the records, the primary row first (`arm = primary`, the primary analysis directory) then `extended_covariates` and `latent_adjusted` (directories relative to the output root): `covariates` (`;`-joined), `conditioning` (`none` or `stage`), the CV-selected rank and `cv_mean_auroc`; per statistic (`delta`, `angle`, `shape`) the observed value, the RRPP p-value, whether p ≤ α and the null mean, sd and q95; for `latent_adjusted` `residual_df`, `cv_min_training_df` and `feature_r2_{q05,q25,q50,q75,q95,mean,frac_gt_0_5}` (empty in the other rows). Floats in their shortest round-trip form, no runtime column, so the file is byte-identical across reruns |
| `comparison/<cohort>/latent_sensitivity.md` | the same values as one table per statistic, a run table (selection, rank, CV AUROC, covariates, conditioning, runtime, directory), the latent-adjustment block, the mediator check with flagged rows first, and the findings, under a header with the reading rules |
| `comparison/<cohort>/mediator_check.csv` | `cohort`, `covariate`, `subset`, `n`, `spearman`, `flagged` |
| `comparison/<cohort>/latent_sensitivity_manifest.json` | `cell_type`; the covariate set (entered, primary, added, sources, scale), the conditioning, per arm its description, covariates, setting, directory, file hashes, rank, CV AUROC, runtime, MOTCO's latent-adjustment block and findings; seed, permutations, α, `n_jobs`; every input path and sha256 (the five export files, the technical table and its manifest, the pseudobulk manifest); the reference analysis (directory, file hashes, covariates, rank, CV AUROC, MOTCO version, seed, permutations); the mediator check's file hash and flags; software (incl. `motco`, `scikit-learn`), git revision and findings (the cell-type findings of the manifests read; a CV-selected rank that differs; a rejection status that differs, recorded as covariate-set or basis dependence; each mediator flag; `n_jobs ≠ 1`) |

Two runs on the same inputs, MOTCO release, seed and `n_jobs` give identical arm statistics, p-values, ranks and null summaries and byte-identical `latent_sensitivity.csv`, `mediator_check.csv` and arm CSVs. Each arm's `cli_equivalent` reproduces it with `motco analyze` (eight `--covariate` flags and, for `latent_adjusted`, `--latent-adjustment stage`).

The first latent axis's technical correlations of an arm come from `axis-diagnostics` pointed at it, e.g. `sead-sexdiff axis-diagnostics --analysis-dir results/motco_latent_sensitivity/latent_adjusted/motco_analysis --out-dir results/motco_latent_sensitivity/diagnostics/latent_adjusted` (with the primary export as `--input-dir`; for microglia the same with `--cell-type microglia` and the paths under `results/microglia/`, which is the layout the positive control's `latent_adjusted` basis reads). There is no rank sweep for an arm, which the diagnostic records as a finding; its tables label the analysis it reads `primary`, so the report names the arm; its partial correlations still condition on the six analysis covariates, the diagnostic's definition.

Reading rules, binding for the report of each cell type the arms run on (`reports/motco-<cell type>-latent-adjustment-<date>.md`: `reports/motco-astrocyte-latent-adjustment-2026-10-21.md`, `reports/motco-microglia-latent-adjustment-<date>.md`):

- The primary analysis remains the primary result, and MOTCO's power comparison (`candidate_fails`) excludes the latent adjustment as its replacement; both arms are sensitivity to it. `cohort_30` is primary and `cohort_50` sensitivity. For microglia the arms are additionally the basis of the restated positive control (`positive-control --basis latent_adjusted`, below), whose reading is predeclared there.
- A statistic whose rejection status differs between `primary` and `extended_covariates` is covariate-set dependence, and one that differs between `extended_covariates` and `latent_adjusted` is basis dependence; neither is a result of the arm. A CV-selected rank that differs is reported beside `angle`, not read as better or worse.
- `angle` and `shape` non-rejections stay uninformative and a `delta` rejection is a difference, not specifically of size, in every arm.
- An arm's first-axis technical correlation is a measurement-space finding only, never a reinterpretation of any arm's statistics and never grounds to adopt an arm.
- Limitations carried: the mediator limitation with the mediator check's flags (a flagged fraction may carry stage signal that the `latent_adjusted` basis removed; no simulation can tell a confounder from a mediator); MOTCO's latent-adjustment operating characteristics, measured at its Phase 6 design only (Type I `adopt_stage`; power `candidate_fails`, with `shape` power lost under a nuisance axis and small-magnitude power lost without one); no operating characteristic measured on these blocks or at this eight-covariate set, so neither arm carries a Type I or power claim of its own; and the nuclei counts and the fractions both proxy library quality.

### Pseudobulk stage DE (DESeq2 cross-check)

```
sead-sexdiff stage-de --n-cpus 8                     # both cohorts, both omics
sead-sexdiff stage-de --cohort cohort_30 --omic rna --force
```

Options: `--cell-type {astrocyte,microglia}` (default `astrocyte`; every directory default below is the cell type's, under its results root, shown here for `astrocyte`), `--cohort {cohort_30,cohort_50,all}` (default `all`), `--omic {rna,atac,both}` (default `both`; the omics fitted in this run, while the summary covers every omic whose fits are complete in the cohort directory), `--pseudobulk-dir` (default `results/pseudobulk`; `rna_counts.h5ad`, `atac_counts.h5ad`, `donor_table.csv`, `manifest.json`, read and never written), `--input-dir` (default `results/motco_input`; the primary export, read and never written), `--analysis-dir` (default `results/motco_analysis`; the primary analysis, read and never written), `--out-dir` (default `results/stage_de`; may not lie inside an input directory), `--n-cpus` (pydeseq2 worker processes, default 1; the tables do not depend on it, but warnings raised in the worker processes are not captured, so a fit made with more than one CPU has an incomplete warning record; recorded), `--force` (refit the requested omics even when their fits are complete, replacing their directories and the cohort's summary files). The models, covariates, reference levels and α are constants, not flags. The step needs `results/pseudobulk/`, `results/motco_input/` and `results/motco_analysis/` on the same workstation; it reads no donor object.

The MOTCO analysis cannot tell an astrocyte pseudobulk with little stage signal from one whose stage signal is the same in both sexes. A per-feature model on the same donors and features can. Per cohort (`cohort_30` primary, `cohort_50` sensitivity) and omic, the step fits with pydeseq2 the covariate-adjusted models of the registry `STAGE_DE_MODELS` on the raw pseudobulk counts of the cohort donors:

| model | design terms | tested coefficients (design column) |
|---|---|---|
| `factor_additive` | `~ sex + stage` | `stage_1_vs_0` (`stage[T.1]`), `stage_2_vs_0` (`stage[T.2]`) |
| `factor_interaction` | `~ sex + stage + sex:stage` | `male_x_stage_1` (`sex[T.Male]:stage[T.1]`), `male_x_stage_2` (`sex[T.Male]:stage[T.2]`) |
| `trend_additive` | `~ sex + stage_num` | `stage_num` (`stage_num`) |
| `trend_interaction` | `~ sex + stage_num + sex:stage_num` | `male_x_stage_num` (`sex[T.Male]:stage_num`) |

Every formula continues with the omic's covariates `STAGE_DE_COVARIATES`: `age_at_death`, `pmi`, `rin`, `years_education` and the omic's own nuclei count on the log10 scale (`log10_n_rna` for RNA, `log10_n_atac` for ATAC; the analysis covariates without the other omic's count), continuous, on their recorded scale, from the primary export's `covariates.csv` (the log10 counts derived as `analyze` derives them); a missing value names the donor and the column, and nothing is imputed. `sex` (reference `Female`) and `stage` (three levels, reference `0`) are treatment-coded categoricals (`STAGE_DE_REFERENCE_LEVELS`), `stage_num` is 0/1/2. Each tested coefficient gets one Wald test, a unit contrast vector on the design-matrix column found by name (never by position; an unresolvable name is an error), with `DeseqStats(alpha = 0.05)` and pydeseq2's defaults: Benjamini–Hochberg `padj` over the tested features that pass pydeseq2's independent filter (a baseMean cutoff chosen per coefficient; a filtered feature has no `padj` and counts as not significant), no LFC shrinkage, no permutation null. Cook's refit and filtering are on but inactive by construction: with continuous covariates no two donors share a design row, so no donor is replaceable and no outlier is filtered (as in R's DESeq2), and the per-feature results carry no outlier protection. `factor_interaction` is the primary interaction model (MOTCO's ordered-categorical stages); `trend_interaction` is a second predeclared model, not a fallback. pydeseq2 0.5.3 has no likelihood-ratio test, so a factor model's two coefficients are also counted as their `union` (significant on either), which carries no joint error-rate claim. When pydeseq2's parametric dispersion trend does not converge it switches to a mean trend (R's DESeq2 would substitute a local fit); the trend used and the warning are recorded in the fit record and the warning as a finding, beside the counts of features whose gene-wise or MAP dispersion or LFC fit did not converge.

The tested features are those that passed the primary export's chromosome and abundance filters (`features_<omic>.csv` status `variance` or `kept`, `STAGE_DE_FEATURE_STATUSES`), in the export's order; size factors (median of ratios) and dispersions are fitted on that set, not on MOTCO's 5,000 features and not on all features, and the export's TMM factors are not reused. Every table row carries `exported` (status `kept`: one of MOTCO's 5,000 features). Counts are reported on the exported set (the primary reading) and on the abundance-passing set; `padj` is always DESeq2's over all tested features, never recomputed on the exported subset.

Predeclared concordance rule (`STAGE_DE_CONCORDANCE_RULE`), per cohort and omic, evaluated from the generated files: `interaction_detected` (at least one exported feature with `padj ≤ 0.05` on either interaction coefficient of `factor_interaction`) against `delta_rejected` (the primary analysis's `delta` rejection status under `ANALYSIS_REJECTION_RULE`, from its `result.json`). The cell is `concordant` when both are true or both false and `discordant` otherwise, written as a finding naming the cohort, omic and count. The rule reads a count and a status, never a p-value magnitude, and changes no MOTCO output.

Before anything is fitted, for every requested cohort, the step checks that every input exists (a missing one names the file and the command that produces it), that the pseudobulk run is complete and the cohort analysable, that the export manifest records the sha256 of the pseudobulk `manifest.json` now on disk, and that the export was made with the primary settings (no gene-set exclusion, variance selection, no normalization or residualization). The export records no count-file hash, so the counts are bound to it by reproduction: `preprocess_block` (the export's own code path) is rerun on the cohort donors' counts with the chromosomes of `features_<omic>.csv` and the recorded thresholds and `k` (the TMM trims and the float format are the module constants; recorded values that differ stop the step naming the setting). The statuses must equal `features_<omic>.csv` and the per-donor library sizes the manifest's `tmm` block exactly, and the TMM factors (`STAGE_DE_TMM_RTOL`), every tested feature's log2-CPM variance (`features_<omic>.csv` `variance`) and the exported values (`STAGE_DE_VALUES_RTOL`) must match to rounding, so a change confined to tested-but-not-exported features is caught too (it is not a hash: a change that leaves all of these unchanged would pass). A difference names the count file and the `export` command. The cohort donors, sex and stage must agree between `metadata.csv`, `donor_table.csv` and the count files, and `k` and the cells between the donor table and the export manifest. The primary analysis must match the export (`_check_reference`: the five export files' sha256, the installed MOTCO version, seed, permutations, CV selection). Every model's design matrix is built and checked for full rank, and every tested column resolved by name, before any fit. Nothing is written when one of these checks fails. Every input file is hashed before it is read and again after each cohort's fits; a change detected then stops the step after that cohort's tables and fit records (their fingerprints carry the old hashes, so a rerun treats them as stale) and before its summary, Markdown and manifest.

Outputs under `results/stage_de/<cohort>/` (each file written as `.tmp` and renamed):

| path | content |
|---|---|
| `<omic>/<model>__<coefficient>.csv` | one DESeq2 table per tested coefficient (six per omic): `feature`, `baseMean`, `log2FoldChange` (the coefficient, log2 scale), `lfcSE`, `stat`, `pvalue`, `padj` (empty when independent filtering removed the feature), `exported`; floats in their shortest round-trip form |
| `<omic>/<model>__fit.json` | the fit record: the fingerprint (inputs' sha256, formula, tested coefficients, covariates, reference levels, α, donors, pydeseq2 version), the fit (formula, design-matrix columns, resolved columns, size factors, dispersion trend requested and used, Cook's counts, non-converged gene-wise / MAP / LFC fits, missing and significant `padj` per coefficient, settings, warnings, `n_cpus`, runtime), the reproduction of the export and each table's sha256. A model whose record and tables match the inputs now on disk is skipped without `--force` (the step resumes); a record made from other inputs is an error |
| `stage_de_summary.csv` | per omic × model × coefficient (and `union` for the factor models) × feature set (`exported`, `abundance`): `kind` (`wald` / `union`), `design_column`, `n_tested`, `n_padj_missing` (for `union`: missing on both coefficients), `n_significant` (padj ≤ α), `min_padj`; the concordance cells (`interaction_detected`, `delta_rejected`, `concordance`) on the `factor_interaction` / `union` / `exported` row of each omic. No runtime or path, so the file is byte-identical across reruns |
| `stage_de.md` | rendered from the summary and the manifest: the reading rules, the concordance table, one count table per omic, the models with formulas and resolved columns, the covariates, the feature sets, the fits, the inputs with sha256 and the reproduction, and the findings |
| `stage_de_manifest.json` | inputs with paths and sha256 (the count files included), the pseudobulk and export manifest hashes, the reproduction per omic, the primary analysis (directory, file hashes, `delta` status, rank, MOTCO version, seed, permutations), the models, covariates, reference levels, α, significance and normalization notes, the feature counts per set, the fit records, the concordance rule and cells, `omics` / `complete`, the run (models fitted and skipped, `n_cpus`), software (incl. `pydeseq2`, `formulaic`, `statsmodels`, `motco`), git revision and findings (a discordance; an omic not fitted yet; each pydeseq2 warning; fits made with `n_cpus ≠ 1`, read from the fit records) |

`run_stage_de(..., debug_dir=...)` (Python API only) also writes `<cohort>_<omic>_counts.csv`, `<cohort>_<omic>_design.csv` and, per fitted model, `<cohort>_<omic>_<model>_{trend,dispersions,size_factors}.csv` for `scripts/check_stage_de_deseq2.R` (in the Nix shell). It compares log2 fold changes, Wald statistics and the significant set twice: a full DESeq2 fit with the dispersion trend type pydeseq2 used (pydeseq2 switches a non-converged parametric trend to a mean trend, R's DESeq2 to a local one), and R's Wald test on pydeseq2's size factors and dispersions (plug-in), so a difference can be placed in the dispersion estimation or in the GLM and test. The 1e-3 tolerance is evaluated on features whose fitted means are all ≥ 0.5 (DESeq2 floors mu at 0.5 in the standard errors) and whose pydeseq2 LFC fit converged; the others are counted apart, and every feature enters the comparison of the significant sets. Exit status 0 when both pass, 2 when only the full fit's tolerance fails, 3 when the trend types differ, 4 when another comparison fails, 1 for an R error. The per-model files come from the fits of that call, so `debug_dir` with a complete model and no `force=True` is an error. pydeseq2 0.5.3 lacks DESeq2's reset of non-improving gene-wise dispersion fits, so the full comparison is not expected to pass the tolerance on every feature.

The step reads sex and stage, so it is a check beside the analysis and never a filter: its outputs are an input of no `export`, `analyze` or sensitivity run, and no feature set is derived from it.

Reading rules, binding for the report (`reports/motco-astrocyte-stage-de-<date>.md`):

- The stage counts describe the stage signal in the pseudobulks. The interaction counts and the concordance cells cross-check the primary analysis; they are not a second test of its hypothesis.
- A discordance in either direction is reported as such, beside both results: per-feature FDR over thousands of features and one latent-space permutation test ask related, not identical, questions. A zero interaction count is "none detected at this n", never "no interaction" (cells of 9–28 donors).
- The `union` counts carry no joint error-rate claim; the per-coefficient counts stand beside them. The FDR of `padj ≤ 0.05` is approximate (Wald tests with estimated dispersions; about 0.06–0.09 realized in simulations under DESeq2's own model), and Cook's handling is inactive with continuous covariates, so the per-feature results carry no outlier protection.
- The primary analysis stays primary with its non-rejection rules (`angle` and `shape` non-rejections uninformative; a `delta` rejection a difference, not specifically of size); the covariate set and the feature set are reported beside the counts, and DESeq2's normalization is not the export's TMM.

### Pseudobulk stage DE, limma-voom (second method)

```
sead-sexdiff stage-de-limma                          # both cohorts, both omics (in the Nix shell)
sead-sexdiff stage-de-limma --cohort cohort_30 --omic rna --force
nix develop -c uv run sead-sexdiff stage-de-limma    # outside the Nix shell
```

Options: `--cell-type {astrocyte,microglia}` (default `astrocyte`; every directory default below is the cell type's, under its results root, shown here for `astrocyte`), `--cohort {cohort_30,cohort_50,all}` (default `all`), `--omic {rna,atac,both}` (default `both`; the omics fitted in this run, while the summary covers every omic whose fits are complete in the cohort directory), `--pseudobulk-dir`, `--input-dir`, `--analysis-dir` (as for `stage-de`; read and never written), `--deseq2-dir` (default `results/stage_de`; the DESeq2 step's outputs, read and checked, never written), `--out-dir` (default `results/stage_de_limma`; may not lie inside an input directory or the DESeq2 directory), `--force` (refit the requested omics even when their fits are complete). The method, the joint tests, α, the TMM tolerance and the concordance rule are constants, not flags; there is no `--n-cpus` (limma runs single-threaded). The step needs `results/pseudobulk/`, `results/motco_input/`, `results/motco_analysis/` and `results/stage_de/` on the same workstation, and R with limma and edgeR: it runs `Rscript` from the Nix shell (`flake.nix` provides R 4.6.1, limma 3.68.4, edgeR 4.10.1), and stops before writing anything, naming the Nix shell, when `Rscript` is not on `PATH`.

It fits the registered stage-DE models (`STAGE_DE_MODELS`, the table above) a second time with limma-voom, on the same cohort donors, raw counts, design matrices (the formulaic matrix pydeseq2 fits, with the same column names), covariates, reference levels and abundance-passing features as `stage-de`, in the export's order, each row carrying `exported`. It sits beside the DESeq2 step and never replaces it: the DESeq2 results stay the reported ones. The method was predeclared before any limma run (`STAGE_DE_LIMMA_METHOD`):

1. edgeR TMM on the tested counts (library size = their column sums) with the export's trims (`TMM_LOGRATIO_TRIM` 0.3, `TMM_SUM_TRIM` 0.05) and edgeR's other defaults. The factors must equal the export manifest's to an absolute `STAGE_DE_LIMMA_TMM_ATOL` (1e-8, the tolerance of `scripts/check_tmm_edger.R`), otherwise the cohort stops naming the omic, the largest difference and the `export` command. So the fit uses the export's own normalization, the one MOTCO's blocks were built with (DESeq2 used median of ratios). voom's log-CPM, `log2((count + 0.5) / (lib.size × factor + 1) × 1e6)`, is not the export's prior-2 log2-CPM.
2. `voomWithQualityWeights` with limma's defaults (`span = 0.5`, `adaptive.span = TRUE`, `method = "genebygene"`, `maxiter = 50`, `tol = 1e-5`, `normalize.method = "none"`, no `var.design` / `var.group`; `arrayWeights` shrinks toward equal weights with `prior.n = 10`): voom's precision weights times one quality weight per donor, estimated with each model's own design, so a noisy donor counts less instead of being removed. Each donor's weight per model, and the minimum, median, maximum and mean per sex × stage cell, are recorded and shown, without a gate.
3. `lmFit`, then `eBayes` with its defaults (`trend = FALSE`, `robust = FALSE`, `legacy = NULL`, which resolves to the legacy estimate since every feature has the same residual df; `df.prior` and `s2.prior` recorded).

The R script checks that every one of these defaults is the installed one and stops otherwise; the versions and the script's sha256 enter every fit record, so a record made under other versions is stale.

Tests: a moderated t per tested coefficient, with the DESeq2 step's labels; a joint moderated F per factor model over its two coefficients (`STAGE_DE_LIMMA_JOINT_TESTS`: `stage_joint` for `factor_additive`, `interaction_joint` for `factor_interaction`; df1 = 2, df2 = `df.prior + df.residual`), the joint test pydeseq2 0.5.3 cannot run; and, for the like-for-like comparison only, each factor model's `union` of its two moderated t's, which carries no joint error-rate claim. A feature is significant when its Benjamini–Hochberg adjusted p ≤ 0.05 over all tested features (limma has no independent filter; pydeseq2's per-coefficient baseMean filter makes DESeq2's BH families smaller, mostly in ATAC). limma's joint F combines the moderated t's through the correlation of the unweighted design's coefficients, so with voom's observation weights its claim is approximate. Per joint test, an exact-F diagnostic recomputes each feature's F from its own weighted coefficient covariance (same variance and degrees of freedom, BH over all features); the fit record keeps the largest relative difference and every call at α that changes, each a finding. limma's F stays the reported test.

Before anything is fitted, for every requested cohort, the step runs every check of `stage-de` through its own code path (`_stage_de_prepare`: inputs, pseudobulk manifest hash, the reproduction of the export from the counts, donors, sex, stage, `k`, cells, feature statuses, the primary analysis against the export). Then it checks the DESeq2 outputs: `<deseq2-dir>/<cohort>/stage_de_manifest.json` must be complete for both omics and record, file for file, the sha256 of the inputs now on disk; every DESeq2 table must hash as that manifest records; and the summary's counts and concordance cells must equal a recomputation from the tables. A failure names the file and `stage-de` (`--force` for a run made from other inputs). Finally R checks, per omic, that the donors of the counts, the design and the export's TMM factors agree, the TMM factors, and that each design rebuilt by `model.matrix` from the design table equals the formulaic matrix column by column under the name map (`sex[T.Male]` ↔ `sexMale`, `stage[T.1]` ↔ `stage1`, `Intercept` ↔ `(Intercept)`). Nothing is written when one of these fails. Every input and every DESeq2 file is hashed again after each cohort and must be unchanged.

Predeclared limma concordance rule (`STAGE_DE_LIMMA_CONCORDANCE_RULE`): `interaction_detected` (at least one exported feature with adjusted p ≤ 0.05 on `interaction_joint`) against the primary analysis's `delta` rejection status, read as `stage-de` reads it; `concordant` when both agree, `discordant` otherwise (a finding). It is reported beside the DESeq2 cell (the union of two Wald tests, from the checked DESeq2 summary), which stays the reported cell of `stage-de`, with `methods_agree`; a disagreement is a method-dependence finding naming both counts, and neither cell replaces the other or changes a MOTCO reading.

The method comparison (`method_comparison.csv`) has one row per omic × model × test × feature set: every tested coefficient (Wald ↔ moderated t, kind `coefficient`), every union (`union`), and every joint F against DESeq2's union of the same model (`joint_f_vs_union`, labelled as different tests). Each row gives `n_tested`, the significant count per method, `n_both`, `n_only_deseq2`, `n_only_limma`, whether each method detects anything and `detection_agrees`; a coefficient row adds the Spearman correlations of DESeq2's `log2FoldChange` with limma's `logFC` and of the Wald `stat` with the moderated `t`, and `n_compared`. DESeq2's missing `padj` counts as not significant. The table is descriptive (no p-value of its own, no gate); every row whose detection differs is a method-dependence finding. The `factor_interaction` / `interaction_joint` / `exported` row also carries the two concordance cells and `methods_agree`.

Outputs under `results/stage_de_limma/<cohort>/` (each file written as `.tmp` and renamed):

| path | content |
|---|---|
| `<omic>/<model>__<label>.csv` | per tested coefficient: `feature`, `AveExpr`, `logFC`, `SE`, `t`, `P.Value`, `adj.P.Val`, `B`, `exported`; per joint test (`factor_additive__stage_joint.csv`, `factor_interaction__interaction_joint.csv`): `feature`, `AveExpr`, the two `logFC_<coefficient>`, `F`, `P.Value`, `adj.P.Val`, `exported`; floats in their shortest round-trip form (R writes 17 significant digits) |
| `<omic>/<model>__fit.json` | the fit record: the fingerprint (`stage-de`'s without pydeseq2, plus the joint test, `STAGE_DE_LIMMA_METHOD`, the TMM tolerance, the script's sha256 and the R / limma / edgeR / jsonlite versions), the fit (formula, design columns, resolved tests, the TMM check and factors, R's design check, the installed-defaults check, voom's span, the sample weights, `df.prior`, `s2.prior`, per-test counts, the exact-F diagnostic, versions, runtime), the reproduction of the export and each table's sha256. Resume semantics are `stage-de`'s: a complete model is skipped without `--force`; a record made from other inputs, settings, script or versions is an error |
| `stage_de_limma_summary.csv` | `stage-de`'s summary columns per omic × model × test × feature set, `kind` ∈ {`moderated_t`, `joint_f`, `union`}; the limma concordance cell on the `factor_interaction` / `interaction_joint` / `exported` row. Byte-identical across reruns |
| `method_comparison.csv` | the method comparison above. Byte-identical across reruns |
| `stage_de_limma.md` | rendered from the summary, the comparison and the manifest: the reading rules, the concordance cells side by side, the count tables with the DESeq2 counts beside limma's, the comparison rows, the sample weights per cell, the TMM check, the exact-F diagnostic, the fits, the inputs and the findings |
| `stage_de_limma_manifest.json` | inputs with sha256 (the DESeq2 manifest, summary and every DESeq2 file under `inputs.deseq2`), the reproduction, the primary analysis, the models and joint tests, the method, the significance and log-CPM notes, the TMM and design checks, the feature sets, the fit records, both concordance rules and the cells, `omics` / `complete`, the run, software (incl. R, limma, edgeR and the script's sha256), git revision and findings |

The R side is one packaged script, `src/sead_sexdiff/stage_de_limma.R`, run as a separate process (`Rscript --vanilla`, `OMP_NUM_THREADS = OPENBLAS_NUM_THREADS = 1`) on exchange files, not through rpy2. `run_stage_de_limma(..., debug_dir=...)` (Python API only) keeps them under `<debug_dir>/<cohort>_<omic>/`: `counts.csv` (features × donors), `export_tmm.csv` and `check.json` there, and per model a directory with links to both, `design_matrix.csv`, `design_table.csv` (`%.17g`), `settings.json` and R's outputs (`coef__<label>.csv`, `joint__<label>.csv`, `sample_weights.csv`, `tmm.csv`, `fit.json`). One fit can be rerun by hand in the Nix shell with `Rscript --vanilla src/sead_sexdiff/stage_de_limma.R <debug_dir>/<cohort>_<omic>/<model>`, the checks of one omic with `Rscript --vanilla src/sead_sexdiff/stage_de_limma.R --check <omic_dir> <model_dir>...`, and the versions with `--versions`. With `debug_dir` every requested model must be fitted in the call (a complete model without `force=True` is an error).

The step reads sex and stage, so it is a check beside the analysis and never a filter: its outputs are an input of no `export`, `analyze` or sensitivity run.

Reading rules, binding for the report (`reports/motco-astrocyte-stage-de-limma-<date>.md`):

- The change was proposed after the DESeq2 results were known. The DESeq2 results stay the reported ones, and no later change may promote whichever method finds more.
- A call that differs between the methods is method dependence, not a result of either. The same donors and counts enter both methods, so agreement is not replication.
- A zero count means "none detected at this n". The union counts carry no joint error-rate claim, and the joint F's claim is approximate (the exact-F diagnostic beside it).
- The concordance cells cross-check the primary analysis and are not a second test of its hypothesis; the primary analysis stays primary with its non-rejection rules (`angle` and `shape` non-rejections uninformative; a `delta` rejection a difference, not specifically of size).

### Positive control (microglia against astrocytes)

```
sead-sexdiff positive-control                      # both cohorts, primary basis
sead-sexdiff positive-control --cohort cohort_30 --force
sead-sexdiff positive-control --basis latent_adjusted   # restated on the latent_adjusted arms
```

Options: `--cohort {cohort_30,cohort_50,all}` (default `all`), `--basis {primary,latent_adjusted}` (default `primary`; a name of the registry `POSITIVE_CONTROL_BASES`, below), `--control-root` (default `results/microglia`, the results root of the control cell type, `POSITIVE_CONTROL_CELL_TYPES["control"]`), `--reference-root` (default `results`, the astrocyte root; read, never rerun), `--out-dir` (default `<control root>/positive_control`; each cohort goes to `<out-dir>/<cohort>/` on the primary basis and to `<out-dir>/<basis>/<cohort>/` otherwise), `--force` (replace existing cohort directories).

It evaluates, per cohort and from the generated files of both cell types only, the expectations written before any microglia run (`POSITIVE_CONTROL_EXPECTATIONS`, constants, never changed from results):

| | expectation | read from | holds when |
|---|---|---|---|
| E1 | stage signal | both roots' `stage_de/<cohort>/stage_de_summary.csv` (`rna`, `factor_additive`, `union`, `exported`) and `stage_de_limma/<cohort>/stage_de_limma_summary.csv` (`rna`, `factor_additive`, `stage_joint`, `exported`): `n_significant` at FDR 0.05 | microglia > astrocyte under **both** methods, compared within each method; under one only E1 is `method_dependent`, which counts as not holding, and the method under which it held is named |
| E2 | latent separation | both roots' `motco_analysis/<cohort>/result.json`: `latent.cv_mean_auroc` | microglia > astrocyte |
| E3 | axis identity | the microglia root's `motco_diagnostics/<cohort>/axis_correlations.csv`: raw Spearman of the primary `pls_0` | \|ρ\| with `stage` > the largest \|ρ\| with a technical variable (sources `donor_technical.csv` and `manifest.json (export, tmm)`, `POSITIVE_CONTROL_TECHNICAL_SOURCES`; the six analysis covariates and `male` are not technical variables) |

The control holds when E1, E2 and E3 all hold in `POSITIVE_CONTROL_PRIMARY_COHORT` (`cohort_30`); the `cohort_50` table is descriptive. Ties do not hold. No expectation reads a p-value, a rejection status or a trajectory statistic: the values are compared, not tested.

Basis (`POSITIVE_CONTROL_BASES`, a registry, never a free directory, never extended from results): `primary` (the default) reads E2 and E3 from each cell type's primary analysis (`<root>/motco_analysis/<cohort>/`) and the diagnostic made from it (`<root>/motco_diagnostics/<cohort>/`), today's inputs, output path and bytes; `latent_adjusted` reads them from each cell type's `latent_adjusted` arm of `latent-sensitivity` (`<root>/motco_latent_sensitivity/latent_adjusted/motco_analysis/<cohort>/`) and the diagnostic made on that arm (`<root>/motco_latent_sensitivity/diagnostics/latent_adjusted/<cohort>/`, written by `axis-diagnostics --analysis-dir ... --out-dir ...`), and writes to `<out-dir>/latent_adjusted/<cohort>/`. E1 is basis-independent (the stage-DE summaries are read unchanged, so its rows equal the primary basis's; the manifest records `e1_basis_independent`). The expectations, their comparisons, E3's technical sources (the two adjusted fractions stay in the compared set) and the all-three rule do not change with the basis. `extended_covariates` gets no restated control: with the covariate-blind basis its latent space, CV rank and AUROC equal the primary's by MOTCO's contract, so E2 and E3 would repeat the primary's values; the report states and checks that equality, and an import-time assertion ties the registry's only non-primary basis to an arm of `LATENT_SENSITIVITY_ARMS` with a latent adjustment. The restated control is a diagnosis of the primary-basis control, which stays the control's outcome; its predeclared reading is under the reading rules below.

Before anything is written, every cohort's inputs are checked: each file is present (a missing one stops the step naming it and the command that produces it, e.g. `sead-sexdiff stage-de-limma --cell-type microglia`), every manifest read records its root's cell type (`microglia` for the control root, `astrocyte` for the reference root; an astrocyte manifest written before the registry has no `cell_type` and reads as astrocyte, listed as a finding), and the chain holds: both stage-DE manifests and the diagnostic manifest name the `result.json` read for E2 (by sha256), the stage-DE manifests and the analysis manifest name the same export manifest, the limma manifest names the DESeq2 summary read for E1, and the diagnostic's technical table was made from the export's pseudobulk. On the `latent_adjusted` basis the chain also holds per cell type: the arm's manifest records the cell type and the sha256 of the five primary export files now on disk (the check `axis-diagnostics` runs) and the same export manifest as the primary analysis; the arm's `result.json` is CV-selected with the registered conditioning (`stage`); the arm's comparison manifest (`<root>/motco_latent_sensitivity/comparison/<cohort>/latent_sensitivity_manifest.json`) names as its reference the primary `result.json` now on disk, the analysis E1's stage-DE steps are checked against, so E1, E2 and E3 hang off one primary analysis per cell type, and lists this arm directory; and the diagnostic read for E3 records the arm's `result.json`. A missing arm or arm diagnostic names the file and the `latent-sensitivity --cell-type microglia` or `axis-diagnostics --analysis-dir ... --out-dir ...` command. The astrocyte arm manifests predate the registry and read as astrocyte with a finding.

Outputs under `<out-dir>/<cohort>/` (primary basis) or `<out-dir>/<basis>/<cohort>/`:

| file | content |
|---|---|
| `positive_control.csv` | one row per compared quantity (E1 per method, E1 over both, E2, E3) and the overall row: `cohort`, `primary_cohort`, `expectation`, `method`, `quantity`, `control`, `control_value`, `reference`, `reference_value`, `comparison`, `holds`, `outcome` (`holds`, `fails`, `method_dependent`; overall `holds` / `does_not_hold`), `detail` |
| `positive_control.md` | the table, the expectations, the inputs with sha256 and the findings; on a non-primary basis the title names the basis and a paragraph states what E2 and E3 read and the predeclared reading |
| `positive_control_manifest.json` | the cell types, `basis` with its registry entry and `e1_basis_independent`, the expectations and rule, every file read with its sha256 per root (on the `latent_adjusted` basis also the arm's files, the five export files and the comparison manifest, with the arm's directory, hashes, covariates and latent-adjustment block under `inputs.<role>.arm`), the compared values, the outcomes, findings, software, git |

Reading rules, binding for the report (`reports/motco-microglia-positive-control-<date>.md`):

- The control is about the stage signal, not the sex × stage hypothesis. If it holds, the pipeline resolves stage structure on real data and the astrocyte null is not a pipeline failure; if it does not, the latent construction is in question for any cell type, and that is reported.
- The microglia sex × stage analysis is a second case-study result, read under the existing rules: a `delta` rejection is a difference, not specifically of size; `angle` and `shape` non-rejections are uninformative; the CV-selected rank beside `angle`; MOTCO's operating characteristics were measured on neither cell type's blocks. Both stage-DE concordance cells (DESeq2, the reported one, and limma beside it with `methods_agree`) apply to it as to the astrocytes.
- Nothing found in microglia changes any astrocyte output or reading.
- The restated control on `latent_adjusted` (`reports/motco-microglia-latent-adjustment-<date>.md`) is reported beside the primary-basis control with the reading predeclared before any microglia arm ran: if it holds in `cohort_30`, the measured quality gradient explains the primary control's failure (once the library-quality direction is conditioned out of the basis, the pipeline resolves stage in the cell type with the stronger stage signal), without promoting the arm (MOTCO's `candidate_fails` stands), changing the primary-basis outcome or changing the astrocyte reading; it is the diagnosis a later MOTCO-side validation of the latent adjustment on blocks with a dominant nuisance axis or a later revised-primary proposal would cite. If it does not hold, the latent construction is in question beyond the quality gradient (the composition of the `Micro-PVM_*` bundle, the rank-1 selection and the few microglia ATAC nuclei per donor are the stated candidates) and the diagnosis "technical gradient" is not supported. E3 on the adjusted basis carries the mediator caveat in both directions (the adjustment removes the quality direction and, with it, any stage signal the fractions carry, so the stage |ρ| may rise or fall); a failure for that reason is still a failure and the expectation is not re-tuned. `cohort_50` is descriptive. The report also states and checks, per cell type and cohort, that the `extended_covariates` arm's CV-selected rank and CV AUROC equal the primary's; a difference is a finding about MOTCO's contract, recorded, not resolved.

### New: Subset donor AnnData files and write pseudobulk(s) (CLI)

Use the `subset` command to load each donor `.h5ad` in a folder, filter by subclass, apply rudimentary QC, concatenate, and write outputs to disk.

- Single-subclass mode (when `--subclass-label` is provided):
  - Concatenated subset: `results/cleaned_files/<Subclass>_subset.h5ad` by default, or the value of `--out-name` if given.
  - Pseudobulk AnnData is also written automatically as `results/cleaned_files/<Subclass>_pseudobulk.h5ad`.

- ALL-subclasses mode (default when `--subclass-label` is omitted):
  - One concatenated subset file per unique `obs['Subclass']`, named `results/cleaned_files/<Subclass>_subset.h5ad`.
  - One pseudobulk per subset is also generated, named `results/cleaned_files/<Subclass>_pseudobulk.h5ad`.

```
sead-sexdiff subset \
  --folder data/PFC/RNAseq/donor_objects \
  --subclass-label Microglia-PVM \
  --out-dir results/cleaned_files \
  --out-name microglia_subset.h5ad \
  --compression gzip
```

- Disable compression if desired: `sead-sexdiff subset --compression none`
- Memory: files are processed one-by-one; intermediates are freed to keep memory usage low.
- QC: when reading each donor object, two basic filters are applied before concatenation:
  - filter cells with too few detected genes: `sc.pp.filter_cells(adata, min_genes=200)`
  - filter genes detected in very few cells: `sc.pp.filter_genes(adata, min_cells=3)`
- Pseudobulk details: created via `decoupler` using defaults (`sample_col='Donor ID'`, `groups_col='Supertype'`, `mode='sum'`, `layer='UMIs'`, `min_cells=10`, `min_counts=1000`). If this step fails (e.g., missing dependency), the CLI raises a clear error. In ALL-subclasses mode, pseudobulk computation runs once per produced subset file.

### New: Differential expression (CLI)

Run DESeq2 differential expression on the pseudobulk file using design `~Sex + ADNC + Sex:ADNC`.

- By default, the command runs DE for ALL supertypes found in `obs.Supertype`.
- Outputs are written to `results/DE`, one CSV and one volcano plot per supertype.

Basic usage (runs all supertypes):

```
sead-sexdiff de
```

Options:

- `--pseudobulk-file`: path to the pseudobulk `.h5ad` (default: `results/cleaned_files/microglia_pseudobulk.h5ad`)
- `--supertype`: value in `obs.Supertype` to subset. If omitted, runs DE for ALL supertypes.
- `--results-dir`: output directory for per-supertype results (default: `results/DE`)
- `--contrast`: contrast vector as comma- or space-separated numbers. If omitted, the default `(0,0,0,0,0,0,0,1)` is used.

Examples:

- Run all supertypes to `results/DE`:

```
sead-sexdiff de
```

- Limit to a single Supertype and customize contrast:

```
sead-sexdiff de \
  --pseudobulk-file results/cleaned_files/microglia_pseudobulk.h5ad \
  --supertype Micro-PVM_2 \
  --results-dir results/DE \
  --contrast 0,0,0,0,0,0,0,1
```

Outputs (per supertype, with `/` in names replaced by `-`):

- Results table CSV: `results/DE/<Supertype>_deseq2_results.csv`
- Volcano plot PDF: `results/DE/<Supertype>_volcano.pdf`

Notes:

- If the DESeq2 design matrix is not full rank for a given supertype (pydeseq2
  emits a UserWarning that the model cannot be fitted), that supertype is
  skipped and no outputs are written for it. The analysis proceeds with the
  remaining supertypes.

### Python API examples

In addition to the CLI, you can call the underlying functions from Python.

- Download all files under a specific folder/prefix (e.g., `donor_objects`). Pass `suffixes=None` to download everything:

```python
from sead_sexdiff.data import download_data

downloaded_paths = download_data(
    prefix="PFC/RNAseq/donor_objects/",
    suffixes=None,  # None => download all files under the prefix
)
# a list of downloaded paths, with .skipped, .failed, .bytes_downloaded, .manifest_path
```

- MTG astrocyte pseudobulks (one donor, or the full run):

```python
from sead_sexdiff.data import pseudobulk_donor, run_pseudobulk

rec = pseudobulk_donor(
    "data/MTG/RNAseq/donors_objects/H19.33.004_SEAAD_MTG_RNAseq_final-nuclei.2026-06-22.h5ad",
    "rna",
)
rec["row"], rec["n_nuclei"], rec["supertype_counts"]

outputs = run_pseudobulk(omic="both", limit=2)
```

- Donor covariates and the MOTCO export:

```python
from sead_sexdiff.data import run_covariates, run_export, tmm_factors, log_cpm

run_covariates(omic="both")
results = run_export(cohort="all", force=True)
results["cohort_30"]["steps"]  # feature counts after each step, per omic
```

- The MOTCO trajectory analysis:

```python
from sead_sexdiff.data import load_analysis_inputs, run_analysis

loaded = load_analysis_inputs("results/motco_input/cohort_30")  # checked blocks + model metadata
results = run_analysis("cohort_30", force=True)
results["cohort_30"]["record"]["p_values"]
```

- The rank and preprocessing sensitivity analyses (the primary analysis is read, never rerun):

```python
from sead_sexdiff.data import PREPROCESSING_VARIANTS, run_preprocessing_sensitivity, run_rank_sensitivity

sweep = run_rank_sensitivity("cohort_30", force=True)
sweep["cohort_30"]["table"]  # the cross-rank table, reference row first

PREPROCESSING_VARIANTS["keep_chrx"]["export"]  # {'excluded_chromosomes': ('chrY', 'chrM')}
results = run_preprocessing_sensitivity("keep_chrx", "cohort_30", force=True)
results["keep_chrx"]["cohort_30"]["table"]     # primary row then variant row
results["keep_chrx"]["cohort_30"]["findings"]
```

- The latent-adjustment sensitivity analysis (the primary analysis is read, never rerun):

```python
from sead_sexdiff.data import LATENT_SENSITIVITY_ARMS, mediator_check, run_latent_sensitivity

LATENT_SENSITIVITY_ARMS["latent_adjusted"]["latent_adjustment"]  # 'stage'
results = run_latent_sensitivity("cohort_30", force=True)
results["cohort_30"]["table"]     # primary row, then extended_covariates and latent_adjusted
results["cohort_30"]["mediator"]  # the mediator check, written before the arms ran
```

- The pseudobulk stage DE (the primary export and analysis are read, never rerun):

```python
from sead_sexdiff.data import STAGE_DE_MODELS, run_stage_de, stage_de_formula

stage_de_formula("factor_interaction", "rna")
# '~ sex + stage + sex:stage + age_at_death + pmi + rin + years_education + log10_n_rna'
results = run_stage_de("cohort_30", omic="rna", n_cpus=8)
results["cohort_30"]["summary"]      # counts per model × coefficient × feature set
results["cohort_30"]["concordance"]  # {'rna': {'interaction_detected': ..., 'delta_rejected': ..., 'cell': ...}}
```

- The limma-voom stage DE beside it (in the Nix shell; the DESeq2 outputs are read and checked, never written):

```python
from sead_sexdiff.data import STAGE_DE_LIMMA_JOINT_TESTS, run_stage_de_limma, stage_de_limma_versions

stage_de_limma_versions()  # {'R': '4.6.1', 'limma': '3.68.4', 'edgeR': '4.10.1', ..., 'script_sha256': ...}
STAGE_DE_LIMMA_JOINT_TESTS["factor_interaction"]  # ('interaction_joint', ('male_x_stage_1', 'male_x_stage_2'))
results = run_stage_de_limma("cohort_30", debug_dir="scratch/limma_debug")
results["cohort_30"]["comparison"]   # method_comparison.csv
results["cohort_30"]["concordance"]  # {'rna': {'deseq2': {...}, 'limma': {...}, 'delta_rejected': ..., 'methods_agree': ...}}
```

- Subset a single donor `.h5ad` by subclass and free the original object:

```python
from sead_sexdiff.data import subset_adata

adata_sub = subset_adata(
    "data/PFC/RNAseq/donor_objects/H19.33.004_SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad",
    subclass_label="Microglia-PVM",
)
```

Note: `subset_adata` applies the same rudimentary QC by default (cells `min_genes=200`, genes `min_cells=3`).

- Subset every donor file in a folder and concatenate them:

```python
from sead_sexdiff.data import subset_and_concat_folder

adata_concat = subset_and_concat_folder(
    "data/PFC/RNAseq/donor_objects",
    subclass_label="Microglia-PVM",
)
```

- Build and write the concatenated subset to disk (gzip by default):

```python
from sead_sexdiff.data import write_subset_from_folder

out_path = write_subset_from_folder(
    folder="data/PFC/RNAseq/donor_objects",
    subclass_label="Microglia-PVM",
    out_dir="results/cleaned_files",
    out_name="microglia_subset.h5ad",
    compression="gzip",
)
print(out_path)
```

Notes
- The subsetting functions expect an `obs['Subclass']` column in each AnnData object. An error is raised if it is missing.
- To download “all files” under a prefix from Python, use `download_data(..., suffixes=None)`. The CLI `download` command defaults to CSV/HTML unless you specify particular suffixes or explicit `--file` values.
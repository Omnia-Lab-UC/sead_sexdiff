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
│       └── data.py
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

Donor-level astrocyte pseudobulks from the SEA-AD middle temporal gyrus (MTG) releases, for the MOTCO case study.

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

Options: `--omic {rna,atac,both}` (default `both`), `--data-root`, `--out-dir` (default `results/pseudobulk`), `--limit N`, `--donors id,...`, `--force`.

What it does, per donor object (one file in memory at a time):

- Excludes donors whose ADNC is `Reference`.
- Selects astrocytes through the per-release taxonomy mapping `TAXONOMY` in `src/sead_sexdiff/data.py` (`Subclass == "Astrocyte"`, supertypes `Astro_1`–`Astro_5`, `Astro_6-SEAAD`). Any `Subclass` or astrocyte `Supertype` label not in the mapping stops the run.
- RNA: sums the raw `UMIs` layer over astrocyte nuclei with ≥ 200 detected genes. ATAC: sums the raw peak counts in `X` over all astrocyte nuclei (release `final-nuclei` QC only). Both matrices are checked to be non-negative integers; no feature filtering or normalization.
- Requires the same ordered feature set across donors of an omic.
- Caches the result in `results/pseudobulk/per_donor/<omic>/<donor>.npz`; a rerun skips cached donors.

Outputs in `results/pseudobulk/`:

| file | content |
|---|---|
| `rna_counts.h5ad` | donors × 36,601 genes, int64 CSR raw UMI sums; `obs` = donor table, `var` = gene annotations |
| `atac_counts.h5ad` | donors × 218,882 peaks, int64 CSR raw count sums; `obs` = donor table |
| `donor_table.csv` | `donor_id`, `sex`, `adnc_raw`, `stage`, `n_rna`, `n_atac`, `cohort_30`, `cohort_50` (one row per non-Reference donor in either omic) |
| `supertype_composition.csv` | `donor_id`, `omic`, `supertype`, `n_nuclei` (nuclei summed per donor × supertype) |
| `cohort_summary.csv` | donors per sex × stage for each cohort, with totals |
| `manifest.json` | release files (key, size, ETag), taxonomy mapping, QC, matrix source, excluded donors, thresholds, cohort counts vs preregistered, partial flag, software versions, git revision |

Donor table rules (scientific choices; change only deliberately):

- `stage`: `Not AD` and `Low` → `0`, `Intermediate` → `1`, `High` → `2`. Other values are an error.
- Sex and ADNC must agree between the RNA and ATAC objects of a donor; sex must have exactly two labels (a warning only in partial runs).
- `cohort_30` / `cohort_50`: ≥ 30 / ≥ 50 astrocyte nuclei in **both** omics.
- A full run compares the cohorts with the preregistered counts (≥ 30: n = 80, F 11/10/28, M 9/10/12 over stages 0/1/2; ≥ 50: n = 74) and records any deviation under `findings` in the manifest. A cohort with an empty sex × stage cell is flagged as not analysable.

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

On a workstation, pull the outputs (the per-donor cache is excluded unless `--with-cache`; extra arguments go to rsync):

```
scripts/sync_pseudobulk.sh -n          # dry run
scripts/sync_pseudobulk.sh
scripts/sync_pseudobulk.sh --with-cache
```

### Donor covariates

```
sead-sexdiff covariates                 # RNA values, cross-checked against ATAC
sead-sexdiff covariates --omic rna      # RNA only, no cross-check
```

Options: `--omic {rna,both}` (default `both`), `--data-root`, `--out-dir` (default `results/pseudobulk`).

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

### MOTCO input export

```
sead-sexdiff export                     # both cohorts
sead-sexdiff export --cohort cohort_30 --force
```

Options: `--pseudobulk-dir` (default `results/pseudobulk`), `--out-dir` (default `results/motco_input`), `--cohort {cohort_30,cohort_50,all}` (default `all`), `--data-root` (only to locate the annotation), `--force` (replace existing cohort directories; without it an existing cohort directory is an error).

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

Options: `--cohort {cohort_30,cohort_50,all}` (default `all`), `--input-dir` (default `results/motco_input`), `--out-dir` (default `results/motco_analysis`), `--n-jobs` (parallel jobs for MOTCO's CV and RRPP, `-1` = all CPUs; default 1; MOTCO seeds one RRPP stream per worker, so other values change the permutation draws, p-values and null summaries, though not the statistics, rank or latent space; the reported results use 1), `--force` (replace existing cohort directories; without it an existing cohort directory is an error and nothing is written).

Reads only the export (`rna.csv`, `atac.csv`, `metadata.csv`, `covariates.csv`, `manifest.json` per cohort), so it runs on any workstation. It needs MOTCO ≥ 0.8.0, a dependency pinned to the tag `v0.8.0` of [tomszar/MOTCO](https://github.com/tomszar/MOTCO); an older install stops the run before any file is read. It stops when a file is missing (naming the `export` command), when the four tables list different donors, or when the export manifest's donor or feature counts disagree with the tables.

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
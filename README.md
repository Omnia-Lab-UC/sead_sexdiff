# sead_sexdiff

Utilities and CLI for SEA-AD sex-difference analyses using Python >= 3.12. Includes tools to download data from the public SEA-AD S3 bucket and to subset/concatenate AnnData (.h5ad) donor files.

Project layout:

```
.
├── data/            # place raw/input data here (ignored by package build)
├── results/         # analysis outputs, figures, tables, etc.
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

# submit (download --omic both, then pseudobulk --omic both)
mkdir -p results/pseudobulk/logs
sbatch scripts/pseudobulk.sbatch
```

`scripts/pseudobulk.sbatch` uses partition `512x1024`, 4 CPUs, 16 GB, 12 h, BLAS threads pinned to 1, `SEAD_DATA_ROOT=/home1/tgonza/sead_sexdiff/data`, and logs to `results/pseudobulk/logs/`. Resubmitting is safe: complete downloads and cached donors are skipped. Do not submit to the `full` partition: its node n14 has glibc 2.17 and cannot load the wheels `uv sync` installs on the login node (glibc 2.34); `512x1024` and `gpus` nodes match the login node.

On a workstation, pull the outputs (the per-donor cache is excluded unless `--with-cache`; extra arguments go to rsync):

```
scripts/sync_pseudobulk.sh -n          # dry run
scripts/sync_pseudobulk.sh
scripts/sync_pseudobulk.sh --with-cache
```

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
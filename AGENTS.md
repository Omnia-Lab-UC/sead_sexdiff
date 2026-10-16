# AGENTS.md

Guidance for AI coding agents working in this repository. This file is the single
source of truth; tool-specific files (e.g. `CLAUDE.md`) only reference it.

## Project

Python package + CLI (`sead-sexdiff`) for sex-difference analyses on the public
SEA-AD single-nucleus RNA-seq data (prefrontal cortex). The pipeline is:

1. `download` – fetch files from the public S3 bucket `sea-ad-single-cell-profiling`
   (anonymous/unsigned boto3 access) into `data/`, preserving key structure.
2. `subset` – read donor `.h5ad` files one at a time, filter by `obs['Subclass']`,
   apply basic QC, concatenate, and write `<Subclass>_subset.h5ad` plus a
   `decoupler` pseudobulk `<Subclass>_pseudobulk.h5ad` to `results/cleaned_files/`.
3. `de` – run pydeseq2 per `obs.Supertype` with design `~Sex + ADNC + Sex:ADNC`,
   writing a CSV and volcano PDF per supertype to `results/DE/`.

`README.md` documents CLI flags and the Python API in detail; keep it in sync
when behavior or defaults change.

## Layout

- `src/sead_sexdiff/cli.py` – argparse CLI (`build_parser`, `main`); thin layer
  that calls into `data.py` and prints `[command]`-prefixed progress messages.
- `src/sead_sexdiff/data.py` – all logic: S3 download, subsetting/QC,
  pseudobulk, differential expression.
- `src/sead_sexdiff/DEG.R` – exploratory R script (MAST mixed model); not wired
  into the CLI.
- `data/`, `results/` – gitignored; never commit their contents.
- `openspec/` – OpenSpec change proposals/specs (see below).

## Environment & commands

- Python `>=3.12,<3.14`, managed with **uv** (`uv.lock` is committed).
- Setup: the Nix flake (`flake.nix`) provides Python 3.12, uv, R, and the R
  packages used by `DEG.R`; its shell hook runs `uv sync --locked`. Enter it
  with direnv (`.envrc` → `use flake`) or `nix develop`. Python packages come
  from uv, not Nix; R packages come from Nix (add them to `flake.nix`).
- Run the CLI: `uv run sead-sexdiff <download|subset|de> [options]`.
- Add dependencies with `uv add <pkg>` so `pyproject.toml` and `uv.lock` stay
  consistent.
- There is no test suite or linter configured yet. Verify changes by running
  the relevant CLI command or importing the function; say so explicitly when
  something could not be verified (e.g. real data not available locally).

## Data & domain notes

- Donor `.h5ad` files are large. Process them one at a time and free
  intermediates; do not load the whole folder into memory at once.
- Expected `obs` columns: `Subclass`, `Supertype`, `Donor ID`, `Sex`, and
  `Overall AD neuropathological Change`, which is renamed to `ADNC` before
  writing subsets.
- Raw counts live in the `UMIs` layer; pseudobulk uses it with `mode='sum'`.
- QC defaults: `min_genes=200` per cell, `min_cells=3` per gene.
- DE: supertypes whose design matrix is not full rank are skipped (pydeseq2
  `UserWarning`), not treated as errors. Output filenames replace `/` with `-`.
- Statistical choices (design formula, contrasts, QC thresholds, pseudobulk
  parameters) are scientific decisions. Don't change them silently; flag them
  and ask.

## Code conventions

- `from __future__ import annotations`, modern type hints (`str | Path`,
  `list[str]`), keyword-only args (`*`) for optional parameters.
- NumPy-style docstrings (`Parameters`, `Returns`, `Notes` sections).
- Accept `str | Path` for paths and convert with `Path(...)` internally.
- Use `print` with a `[command]` prefix for CLI progress output (no logging
  framework yet).

## OpenSpec workflow

Non-trivial features go through OpenSpec (`openspec/changes/`, specs in
`openspec/specs/`): propose → apply → archive. Agent skills/commands for this
live in `.claude/` and `.junie/`.

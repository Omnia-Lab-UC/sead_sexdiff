#!/bin/bash
#SBATCH --job-name=sead-technical
#SBATCH --partition=512x1024
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=results/pseudobulk/logs/technical-%j.out
#SBATCH --error=results/pseudobulk/logs/technical-%j.err

# MTG donor technical covariates on the ing cluster: per donor, the gene-set
# UMI fractions and QC summaries of the nuclei the pseudobulk of the same cell
# type summed (the same per-donor read as scripts/pseudobulk.sbatch: the
# cell type's rows of the UMIs layer plus obs columns; ATAC obs columns only).
# CELL_TYPE (default astrocyte) is passed as --cell-type; inputs and outputs
# are in the cell type's pseudobulk directory (results/pseudobulk for
# astrocyte, results/<cell type>/pseudobulk otherwise).
#
# Needs the downloaded donor objects and the pseudobulk outputs
# (<pseudobulk dir>/rna_counts.h5ad, donor_table.csv, manifest.json and the
# per-donor cache for the nucleus-count cross-check), so run it after
# scripts/pseudobulk.sbatch, from the repository root:
#   sbatch scripts/technical.sbatch
# For another cell type, as for scripts/pseudobulk.sbatch:
#   sbatch --export=ALL,CELL_TYPE=microglia \
#       --output=results/microglia/pseudobulk/logs/technical-%j.out \
#       --error=results/microglia/pseudobulk/logs/technical-%j.err scripts/technical.sbatch
#
# Resubmitting is safe: donors already in <pseudobulk dir>/per_donor_technical/
# are skipped. A login-node run is also acceptable:
#   nohup .venv/bin/sead-sexdiff technical-covariates --omic both \
#       > results/pseudobulk/logs/technical.log 2>&1 < /dev/null &

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-.}"
CELL_TYPE="${CELL_TYPE:-astrocyte}"
if [ "$CELL_TYPE" = astrocyte ]; then pb_dir=results/pseudobulk; else pb_dir="results/$CELL_TYPE/pseudobulk"; fi
mkdir -p "$pb_dir/logs"

export SEAD_DATA_ROOT="${SEAD_DATA_ROOT:-/home1/tgonza/sead_sexdiff/data}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONUNBUFFERED=1  # per-donor progress reaches the log as it happens

cli=.venv/bin/sead-sexdiff
if [ ! -x "$cli" ]; then
    echo "Missing $cli: run 'uv sync --locked' on the login node first" >&2
    exit 1
fi

echo "[sbatch] host=$(hostname) job=${SLURM_JOB_ID:-none} data_root=$SEAD_DATA_ROOT cell_type=$CELL_TYPE dir=$pb_dir"
n_rna=$(find "$SEAD_DATA_ROOT/MTG/RNAseq/donors_objects" -name '*_final-nuclei.*.h5ad' 2>/dev/null | wc -l)
n_atac=$(find "$SEAD_DATA_ROOT/MTG/ATACseq/donor_objects" -name '*_final-nuclei.*.h5ad' 2>/dev/null | wc -l)
echo "[sbatch] donor objects: rna=$n_rna atac=$n_atac"
if [ "$n_rna" -ne 84 ] || [ "$n_atac" -ne 86 ]; then
    echo "Expected 84 RNA and 86 ATAC donor objects; run the download on the login node first" >&2
    exit 1
fi
for f in rna_counts.h5ad donor_table.csv manifest.json; do
    if [ ! -f "$pb_dir/$f" ]; then
        echo "Missing $pb_dir/$f: run scripts/pseudobulk.sbatch (CELL_TYPE=$CELL_TYPE) first" >&2
        exit 1
    fi
done
"$cli" technical-covariates --cell-type "$CELL_TYPE" --omic both

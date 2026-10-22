#!/bin/bash
#SBATCH --job-name=sead-pseudobulk
#SBATCH --partition=512x1024
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=results/pseudobulk/logs/pseudobulk-%j.out
#SBATCH --error=results/pseudobulk/logs/pseudobulk-%j.err

# MTG donor pseudobulks of a registered cell type on the ing cluster: build
# the pseudobulks, donor table, cohort summary and manifest from the
# downloaded donor objects. CELL_TYPE (default astrocyte) is passed as
# --cell-type; the outputs go to the cell type's pseudobulk directory
# (results/pseudobulk for astrocyte, results/<cell type>/pseudobulk otherwise,
# the rule of CELL_TYPES' results_root).
#
# Compute nodes have no internet access, so download on the login node first
# (size-checked, resumable), from the repository root:
#   export SEAD_DATA_ROOT=/home1/tgonza/sead_sexdiff/data
#   mkdir -p results/pseudobulk/logs
#   nohup .venv/bin/sead-sexdiff download --omic both --no-progress \
#       > results/pseudobulk/logs/download.log 2>&1 < /dev/null &
# then, once it reports 0 failed:
#   sbatch scripts/pseudobulk.sbatch
# For another cell type, export CELL_TYPE and send the Slurm logs to its own
# directory (the #SBATCH paths are the astrocyte ones; Slurm needs the log
# directory to exist), e.g. for microglia:
#   mkdir -p results/microglia/pseudobulk/logs
#   sbatch --export=ALL,CELL_TYPE=microglia \
#       --output=results/microglia/pseudobulk/logs/pseudobulk-%j.out \
#       --error=results/microglia/pseudobulk/logs/pseudobulk-%j.err scripts/pseudobulk.sbatch
#
# The job stops early unless all 84 RNA and 86 ATAC objects are present.
# Resubmitting is safe: donors already in <pseudobulk dir>/per_donor/ are
# skipped.

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

echo "[sbatch] host=$(hostname) job=${SLURM_JOB_ID:-none} data_root=$SEAD_DATA_ROOT cell_type=$CELL_TYPE out=$pb_dir"
n_rna=$(find "$SEAD_DATA_ROOT/MTG/RNAseq/donors_objects" -name '*_final-nuclei.*.h5ad' 2>/dev/null | wc -l)
n_atac=$(find "$SEAD_DATA_ROOT/MTG/ATACseq/donor_objects" -name '*_final-nuclei.*.h5ad' 2>/dev/null | wc -l)
echo "[sbatch] donor objects: rna=$n_rna atac=$n_atac"
if [ "$n_rna" -ne 84 ] || [ "$n_atac" -ne 86 ]; then
    echo "Expected 84 RNA and 86 ATAC donor objects; run the download on the login node first" >&2
    exit 1
fi
"$cli" pseudobulk --cell-type "$CELL_TYPE" --omic both

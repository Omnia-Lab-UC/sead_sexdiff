#!/bin/bash
#SBATCH --job-name=sead-pseudobulk
#SBATCH --partition=512x1024
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=results/pseudobulk/logs/pseudobulk-%j.out
#SBATCH --error=results/pseudobulk/logs/pseudobulk-%j.err

# MTG astrocyte donor pseudobulks on the ing cluster: download the RNA and ATAC
# donor objects (size-checked, skips files already present), then build the
# pseudobulks, donor table, cohort summary and manifest.
#
# Usage, from the repository root on the login node (once: `uv sync --locked`):
#   mkdir -p results/pseudobulk/logs
#   sbatch scripts/pseudobulk.sbatch
#
# Resubmitting is safe: downloads skip complete files and the pseudobulk step
# skips donors already in results/pseudobulk/per_donor/.

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-.}"
mkdir -p results/pseudobulk/logs

export SEAD_DATA_ROOT="${SEAD_DATA_ROOT:-/home1/tgonza/sead_sexdiff/data}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

cli=.venv/bin/sead-sexdiff
if [ ! -x "$cli" ]; then
    echo "Missing $cli: run 'uv sync --locked' on the login node first" >&2
    exit 1
fi

echo "[sbatch] host=$(hostname) job=${SLURM_JOB_ID:-none} data_root=$SEAD_DATA_ROOT"
"$cli" download --omic both --no-progress
"$cli" pseudobulk --omic both

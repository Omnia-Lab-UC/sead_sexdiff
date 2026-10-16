#!/bin/bash
#SBATCH --job-name=sead-pseudobulk
#SBATCH --partition=512x1024
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=results/pseudobulk/logs/pseudobulk-%j.out
#SBATCH --error=results/pseudobulk/logs/pseudobulk-%j.err

# MTG astrocyte donor pseudobulks on the ing cluster: build the pseudobulks,
# donor table, cohort summary and manifest from the downloaded donor objects.
#
# Compute nodes have no internet access, so download on the login node first
# (size-checked, resumable), from the repository root:
#   export SEAD_DATA_ROOT=/home1/tgonza/sead_sexdiff/data
#   mkdir -p results/pseudobulk/logs
#   nohup .venv/bin/sead-sexdiff download --omic both --no-progress \
#       > results/pseudobulk/logs/download.log 2>&1 < /dev/null &
# then, once it reports 0 failed:
#   sbatch scripts/pseudobulk.sbatch
#
# The job stops early unless all 84 RNA and 86 ATAC objects are present.
# Resubmitting is safe: donors already in results/pseudobulk/per_donor/ are
# skipped.

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
n_rna=$(find "$SEAD_DATA_ROOT/MTG/RNAseq/donors_objects" -name '*_final-nuclei.*.h5ad' 2>/dev/null | wc -l)
n_atac=$(find "$SEAD_DATA_ROOT/MTG/ATACseq/donor_objects" -name '*_final-nuclei.*.h5ad' 2>/dev/null | wc -l)
echo "[sbatch] donor objects: rna=$n_rna atac=$n_atac"
if [ "$n_rna" -ne 84 ] || [ "$n_atac" -ne 86 ]; then
    echo "Expected 84 RNA and 86 ATAC donor objects; run the download on the login node first" >&2
    exit 1
fi
"$cli" pseudobulk --omic both

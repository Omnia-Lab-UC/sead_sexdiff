#!/usr/bin/env bash
# Pull the derived pseudobulk outputs from the cluster to this workstation.
#
# Usage: scripts/sync_pseudobulk.sh [--cell-type NAME] [--with-cache] [extra rsync args, e.g. -n]
#   --cell-type NAME  the registered cell type whose pseudobulk directory is
#                     synced (default astrocyte: results/pseudobulk; otherwise
#                     results/NAME/pseudobulk, the rule of CELL_TYPES' results_root)
#   --with-cache      also copy the per-donor caches <dir>/per_donor/ and
#                     <dir>/per_donor_technical/ (excluded by default)
#
# Environment: SEAD_CLUSTER (default tgonza@cluster.ing.uc.cl) and
# SEAD_CLUSTER_DIR (default /home1/tgonza/sead_sexdiff).

set -euo pipefail

host="${SEAD_CLUSTER:-tgonza@cluster.ing.uc.cl}"
remote_dir="${SEAD_CLUSTER_DIR:-/home1/tgonza/sead_sexdiff}"
cell_type=astrocyte
if [ "${1:-}" = "--cell-type" ]; then
    cell_type="${2:?--cell-type needs a name}"
    shift 2
fi
excludes=(--exclude per_donor/ --exclude per_donor_technical/)
if [ "${1:-}" = "--with-cache" ]; then
    excludes=()
    shift
fi
if [ "$cell_type" = astrocyte ]; then pb_dir=results/pseudobulk; else pb_dir="results/$cell_type/pseudobulk"; fi

cd "$(dirname "$0")/.."
mkdir -p "$pb_dir"
rsync -av "${excludes[@]}" "$@" "$host:$remote_dir/$pb_dir/" "$pb_dir/"

#!/usr/bin/env bash
# Pull the derived pseudobulk outputs from the cluster to this workstation.
#
# Usage: scripts/sync_pseudobulk.sh [--with-cache] [extra rsync args, e.g. -n]
#   --with-cache   also copy the per-donor caches results/pseudobulk/per_donor/
#                  and results/pseudobulk/per_donor_technical/ (excluded by default)
#
# Environment: SEAD_CLUSTER (default tgonza@cluster.ing.uc.cl) and
# SEAD_CLUSTER_DIR (default /home1/tgonza/sead_sexdiff).

set -euo pipefail

host="${SEAD_CLUSTER:-tgonza@cluster.ing.uc.cl}"
remote_dir="${SEAD_CLUSTER_DIR:-/home1/tgonza/sead_sexdiff}"
excludes=(--exclude per_donor/ --exclude per_donor_technical/)
if [ "${1:-}" = "--with-cache" ]; then
    excludes=()
    shift
fi

cd "$(dirname "$0")/.."
mkdir -p results/pseudobulk
rsync -av "${excludes[@]}" "$@" "$host:$remote_dir/results/pseudobulk/" results/pseudobulk/

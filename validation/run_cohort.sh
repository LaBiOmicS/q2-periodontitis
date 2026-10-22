#!/usr/bin/env bash
# Runs one cohort end to end: download, metadata, processing.
#
#   bash run_cohort.sh <cohort> <data_dir>
#
# <cohort> is a row of cohorts.tsv and has a configs/<cohort>.sh. Uses the
# QIIME 2 environment given by Q2_ENV (default q2-amplicon-2024.10) through
# micromamba, THREADS (default 4) and REF_DIR (default <data_dir>/ehomd).
# Python always runs isolated (-I): downloaded files are untrusted data.
set -euo pipefail

COHORT=${1:?usage: run_cohort.sh <cohort> <data_dir>}
DATA=$(readlink -f "${2:?usage: run_cohort.sh <cohort> <data_dir>}")
HERE=$(cd "$(dirname "$0")" && pwd)
ENV=${Q2_ENV:-q2-amplicon-2024.10}
export THREADS=${THREADS:-4}
export REF_DIR=${REF_DIR:-$DATA/ehomd}
export TMPDIR=${TMPDIR:-$DATA/tmp}
mkdir -p "$TMPDIR"
ACC=$(awk -F'\t' -v c="$COHORT" '$1 == c {print $2}' "$HERE/cohorts.tsv")
[ -n "$ACC" ] || { echo "unknown cohort $COHORT" >&2; exit 1; }

run() { micromamba run -n "$ENV" "$@"; }

[ -s "$DATA/$ACC.runs.tsv" ] && [ -d "$DATA/$ACC/fastq" ] \
  || bash "$HERE/download.sh" "$ACC" "$DATA"
[ -s "$DATA/$COHORT/metadata.tsv" ] \
  || run python -I "$HERE/prepare.py" "$COHORT" "$DATA"
run bash "$HERE/process_dataset.sh" "$HERE/configs/$COHORT.sh" "$DATA/$COHORT"

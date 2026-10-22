#!/usr/bin/env bash
# Reproduces the whole validation: every cohort of cohorts.tsv, then the
# analysis, robustness and figures.
#
#   bash run_all.sh <data_dir>
#
# Needs about 60 GB of disk while the raw reads are present and about a day
# on 4 cores. See README.md.
set -uo pipefail
DATA=$(readlink -f "${1:?usage: run_all.sh <data_dir>}")
HERE=$(cd "$(dirname "$0")" && pwd)
ENV=${Q2_ENV:-q2-amplicon-2024.10}
mkdir -p "$DATA"
for cohort in $(tail -n +2 "$HERE/cohorts.tsv" | cut -f1); do
  echo "$(date +%H:%M) $cohort"
  bash "$HERE/run_cohort.sh" "$cohort" "$DATA" > "$DATA/$cohort.log" 2>&1 \
    || echo "$(date +%H:%M) $cohort FAILED (see $DATA/$cohort.log)"
done
bash "$HERE/analysis/analyse.sh" "$DATA" "$HERE/results"

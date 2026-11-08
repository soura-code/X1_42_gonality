#!/usr/bin/env bash
# Order-independent degree-11 search on X1(42)/F5 (orbit reduction on degrees 1-5).
# Usage:  ./run_split_search.sh [N] [extra options for X1_42_split.py, e.g. --jobs9 2]
# Before the first run, once:   magma -b X1_42_rational_labels.m
# Needs: X1_42_split_lib.m X1_42_split.py X1_42_split_search.m X1_42_split_check.py
#        X1_42_rational_labels.txt X1_42_orbits_d2.txt ... X1_42_orbits_d5.txt
set -u
N=${1:-20}; shift || true
DIR=X1_42_split_d11
for f in X1_42_split_lib.m X1_42_split.py X1_42_split_search.m X1_42_split_check.py \
         X1_42_rational_labels.txt X1_42_orbits_d2.txt X1_42_orbits_d3.txt X1_42_orbits_d4.txt X1_42_orbits_d5.txt; do
  [ -f "$f" ] || { echo "missing $f"; exit 1; }
done
python3 X1_42_split.py --njobs "$N" --outdir "$DIR" "$@" || exit 1
echo "Starting $N jobs..."
for ((k=0; k<N; k++)); do
  nohup magma -b job:=$k njobs:=$N dir:=$DIR X1_42_split_search.m > "$DIR/log_job${k}_of_${N}.txt" 2>&1 &
done
wait
echo
python3 X1_42_split_check.py "$DIR"

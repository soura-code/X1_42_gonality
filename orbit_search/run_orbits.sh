#!/usr/bin/env bash
# Orbits of <sigma_4> on the places of degree 1..5 of X1(42)/F5, in parallel.
# Usage:  ./run_orbits.sh [P1 P2 P3 P4 P5]
#   Pv = number of parallel parts for degree v   (default 1 1 1 4 12, i.e. 19 jobs)
# Then merges the parts (X1_42_orbits_merge.py) into X1_42_orbits_d<v>.txt.
# Needs X1_42_orbits.m, X1_42_orbits_merge.py, X1_42_aut_load.m, X1_42_aut_F5.txt.
set -u
PARTS=(0 ${1:-1} ${2:-1} ${3:-1} ${4:-4} ${5:-12})
for f in X1_42_orbits.m X1_42_orbits_merge.py X1_42_aut_load.m X1_42_aut_F5.txt; do
  [ -f "$f" ] || { echo "missing $f"; exit 1; }
done
rm -f X1_42_orbits_d*_part*_of_*.txt
total=0
for v in 1 2 3 4 5; do
  m=${PARTS[$v]}
  for ((k=0; k<m; k++)); do
    nohup magma -b degs:=$v part:=$k nparts:=$m X1_42_orbits.m \
      > "log_orbits_d${v}_part${k}_of_${m}.txt" 2>&1 &
    total=$((total+1))
  done
done
echo "Started $total jobs. Waiting..."
wait
echo
grep -h "^degree .*part" log_orbits_d*_part*_of_*.txt || true
echo
python3 X1_42_orbits_merge.py

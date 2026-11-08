#!/usr/bin/env bash
# Run X1_42_placeorder_test.m in R separate Magma sessions (in parallel) and
# compare the recorded enumeration orders.   Usage:  ./run_placeorder.sh [R]   (default 4)
set -u
R=${1:-4}
rm -f placeorder_run*.txt
for ((r=1; r<=R; r++)); do
  nohup magma -b run:=$r X1_42_placeorder_test.m > log_placeorder_run$r.txt 2>&1 &
done
wait
for ((r=1; r<=R; r++)); do
  [ -f placeorder_run$r.txt ] || { echo "run $r did not finish (see log_placeorder_run$r.txt)"; exit 1; }
done
echo "Non-unique identifiers (must be 0 for a conclusive test):"
grep "^SHARED" placeorder_run1.txt
same=1
for ((r=2; r<=R; r++)); do
  if ! diff -q placeorder_run1.txt placeorder_run$r.txt > /dev/null; then
    same=0; echo "run $r differs from run 1:"; diff placeorder_run1.txt placeorder_run$r.txt
  fi
done
if [ $same -eq 1 ]; then
  echo "==> All $R sessions enumerate in the SAME order."
else
  echo "==> The order DIFFERS between sessions (see the lines above: SEQ/SET/FACTOR v ...)."
fi

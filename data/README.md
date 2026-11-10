# Data

## Included

- `X1_42_aut_F5.txt` — the twelve automorphisms of X_1(42)/F_5 (written by `automorphisms/X1_42_autgroup.m`).
  The orbit-search scripts read this file; copy it into `orbit_search/` before running them (see the main README).
- `X1_42_aut_F5_readable.txt` — the same twelve automorphisms in human-readable form.
- `X1_42_rational_labels.txt` — output of `orbit_search/X1_42_rational_labels.m`.
- `X1_42_moduli_check.txt` — output of `model/X1_42_moduli_check.m` (all checks pass; about 42 s).
- `X1_42_orbits_d2.txt`, `X1_42_orbits_d3.txt`, `X1_42_orbits_d4.txt`, `X1_42_orbits_d5.txt` —
  the orbits of `<sigma_4>` on the places of degree 2, 3, 4 and 5.
- `X1_42_orbit_reps_d2.txt`, `X1_42_orbit_reps_d3.txt`, `X1_42_orbit_reps_d4.txt`, `X1_42_orbit_reps_d5.txt` —
  one representative per orbit of `<sigma_4>` on the places of degree 2, 3, 4 and 5
  (1, 2, 17 and 54 orbits respectively), written by `orbit_search/X1_42_orbit_reps.m`.

### `search_outputs/` — the orbit-reduced degree-11 search (Section 6)

Output of the run reported in the paper (B = 9, orbit reduction on degrees 1–5, 20 jobs),
produced by `orbit_search/run_split_search.sh 20`.

- `plan.txt`, `summary.txt` — the split into jobs; `summary.txt` lists, for each partition, the number of
  divisors kept and tested (743231 in total).
- `job<i>_of_20.txt` — the input of job i.
- `status_job<i>_of_20.txt` — final status of job i; all 20 jobs report COMPLETED, and the numbers of divisors
  tested add up to 743231.
- `hits_job<i>_of_20.txt` — the functions of degree <= 12 found by job i.
- `log_job<i>_of_20.txt` — the log of job i.

To verify the run, from `orbit_search/`: `python3 X1_42_split_check.py ../data/search_outputs`.

## To be added (not yet uploaded)

- Output of `X1_42_split_check.py` for this run.

# Data

## Included

- `X1_42_aut_F5.txt` — the twelve automorphisms of X_1(42)/F_5 (written by `automorphisms/X1_42_autgroup.m`).
  The orbit-search scripts read this file; copy it into `orbit_search/` before running them (see the main README).
- `X1_42_aut_F5_readable.txt` — the same twelve automorphisms in human-readable form.
- `X1_42_rational_labels.txt` — output of `orbit_search/X1_42_rational_labels.m`.
- `X1_42_orbit_reps_d2.txt`, `X1_42_orbit_reps_d3.txt`, `X1_42_orbit_reps_d4.txt`, `X1_42_orbit_reps_d5.txt` —
  one representative per orbit of `<sigma_4>` on the places of degree 2, 3, 4 and 5
  (1, 2, 17 and 54 orbits respectively), written by `orbit_search/X1_42_orbit_reps.m`.

## To be added (not yet uploaded)

- `X1_42_orbits_d<v>.txt` (v = 1..5) — the full orbits.
- `X1_42_split_d11/` — job files, `status_job*.txt`, `hits_job*.txt` (degree-12 functions) and the
  output of `X1_42_split_check.py` for the run reported in the paper.
- Output of `model/X1_42_moduli_check.m`.

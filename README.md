# The Q-gonality of X_1(42) is 12 — code and data

Companion repository for the paper *The Q-gonality of X_1(42) is 12* by Sourajyoti Maiti
(University of Connecticut). It contains the model of X_1(42) over F_5, the checks that
identify it with X_1(42)_{F_5}, the automorphism computation, and the orbit-reduced
degree-11 search with its verifier.

All computations were run with **Magma V2.29-10** on an Intel Xeon W-2265 workstation
(12 cores, 24 threads, 64 GB RAM).

## Layout

| Folder | Contents |
|---|---|
| `model/` | `FFFc42.txt`: Sutherland's integer model f(x,y) of X_1(42) ([source](https://math.mit.edu/~drew/FFFc42.txt)). `X1_42_moduli_check.m`: proof of Lemma 2.4 (the mod-5 model is X_1(42)/F_5, via the Tate normal form and deg j = 576). `X1_42_irreducibility_mod5.m`: absolute irreducibility of f mod 5. |
| `automorphisms/` | `X1_42_autgroup.m`: computes and verifies the automorphisms, writes `X1_42_aut_F5.txt`. `X1_42_aut_load.m`: reload helpers. `X1_42_check_orbs.m`: permutations of the rational places. `aut_to_latex.py`, `X1_42_automorphisms.tex/.pdf`: explicit formulas (Appendix A). |
| `orbit_search/` | Orbit-reduced search of Section 6 (pipeline below). |
| `data/` | Data files and job logs (see `data/README.md`). |

## Reproducing the checks

```
cd model
magma -b X1_42_moduli_check.m        # about 1 minute
magma -b X1_42_irreducibility_mod5.m
```

Expected output of `X1_42_moduli_check.m`:

```
A: exact constant field F_5, genus 25
B: E/F is an elliptic curve and (0,0) has exact order 42
B': inverse map recovers x, y; F = F_5(b,c)
C: deg j(E) = 576
ALL CHECKS PASSED: F is the function field of X_1(42) over F_5.
```

## Reproducing the orbit-reduced search (Section 6)

Run in `orbit_search/`, with `X1_42_aut_F5.txt` (from `automorphisms/X1_42_autgroup.m`) in the same folder.

1. **Orbits on places of degree <= 5:** `./run_orbits.sh`, then `python3 X1_42_orbits_merge.py`
   and `magma -b X1_42_orbits_check.m`. Writes `X1_42_orbits_d<v>.txt`.
2. **Orbit representatives:** `magma -b X1_42_orbit_reps.m`. Writes `X1_42_orbit_reps_d<v>.txt`.
3. **Labels of the rational places (once):** `magma -b X1_42_rational_labels.m`.
4. **Splitting into jobs and running them:** `./run_split_search.sh 20` (calls `X1_42_split.py`,
   which checks the representative counts against Burnside's lemma, then starts the jobs
   `X1_42_split_search.m`).
5. **Verification of a finished run:** `python3 X1_42_split_check.py X1_42_split_d11`. Accepts the
   run only if every job reports COMPLETED, the kept places of each high degree add up to P_v,
   the bucket fingerprints agree, and the numbers of divisors tested add up to 743231.
   The output of the run reported in the paper is in `data/search_outputs/`; to verify it, run
   `python3 X1_42_split_check.py ../data/search_outputs` from `orbit_search/`.

`X1_42_placeorder_test.m` / `run_placeorder.sh` is a diagnostic showing that the order in which
Magma returns places is not reproducible between sessions; the orbit-reduced search is designed
to be independent of that order.

## Citation

If you use this code, please cite the paper (reference to be added).

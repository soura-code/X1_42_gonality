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
| `model/` | `FFFc42.txt`: Sutherland's integer model f(x,y) of X_1(42) ([source](https://math.mit.edu/~drew/FFFc42.txt)). `X1_42_moduli_check.m`: proof of Lemma 2.4 (the mod-5 model is X_1(42)/F_5, via the Tate normal form and deg j = 576). `X1_42_irreducibility_mod5.m`: absolute irreducibility of f mod 5. `X1_42_check_tw3.m`: confirms t = 1/y, w = x/y and Phi(t,w) = t^23 f(w/t, 1/t), and writes `X1_42_Phi.txt`. `X1_42_Phi.txt`: the defining polynomial Phi(t,w) of F_5(X_1(42)) (Appendix A.1). |
| `automorphisms/` | `X1_42_autgroup.m`: computes and verifies the automorphisms, writes `X1_42_aut_F5.txt`. `X1_42_aut_load.m`: reload helpers. `X1_42_check_orbs.m`: permutations of the rational places. `aut_to_latex.py`, `X1_42_automorphisms.tex/.pdf`: explicit formulas (Appendix A). |
| `orbit_search/` | Orbit-reduced search of Section 6 (pipeline below). |
| `docs/` | `X1_42_orbit_computation.tex/.pdf`: notes on the orbit computation. |
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

## Reproducing the orbit-reduced search (Section 6 of the paper)

### Requirements

* Magma V2.29-10 (the version used for the paper), Python 3, bash.
* For the full search: about 64 GB of RAM. Each of the two jobs for the partition [9] holds all 214,384 places of degree 9 (about 19 GB of memory, about 4,900 s to compute).

### Notation

The scripts work in Magma's function field `AFF` of the model `FFFc42.txt` mod 5, with `t = BaseField(AFF).1 = 1/y` and `w = AFF.1 = x/y`, where `x, y` are the coordinates of the model. The scripts call `w` "y" (variables `YY`, `SL_YY`, `GY`). The defining polynomial `Phi(t,w)` is in `X1_42_Phi.txt`.

### Working folder

The scripts read and write their files by bare file name in the current folder. Put everything in one working folder first:

```
mkdir run
cp model/* automorphisms/* orbit_search/* run/
cp data/*.txt run/
cd run
```

The stored run of the paper is in `data/X1_42_split_d11/`. Copy it as well if you want to re-verify it (route A):

```
cp -r ../data/X1_42_split_d11 .
```

### Route A: verify the stored run (no Magma needed, seconds)

```
python3 X1_42_split_check.py X1_42_split_d11
```

This reads `plan.txt` and the 20 files `status_job<k>_of_20.txt`. It accepts the run only if:

* every job reports COMPLETED;
* in each group with a high degree v = 6, ..., 9, all jobs report the same fingerprint of the buckets of all places of degree v;
* the kept places of each such group add up to all places of degree v, and their fingerprints add up to the common fingerprint;
* the numbers of divisors tested add up to 743,231.

The last line must be

```
==> All 20 jobs COMPLETED, 743231 divisors tested = plan: no function of degree <= 11, hence gon_F5(X1(42)) >= 12.
```

The expected output is saved in `data/X1_42_split_d11/split_check_output.txt`. The Magma logs `log_job<k>_of_20.txt` contain each job's place counts, fingerprints and timings.

### Route B: recompute everything

**Step 1. The model (minutes).**

```
magma -b X1_42_moduli_check.m          # Lemma 2.5: the model is X_1(42) over F_5; must end with ALL CHECKS PASSED
magma -b X1_42_irreducibility_mod5.m   # f mod 5 is irreducible
magma -b X1_42_check_tw3.m             # t = 1/y, w = x/y; writes X1_42_Phi.txt
```

**Step 2. The automorphisms (about 1.6 hours).**

The stored file `X1_42_aut_F5.txt` contains the 12 automorphisms as pull-backs `(sigma_i(t), sigma_i(w))`, in the order returned by `Automorphisms(AFF)`. The generator is `sigma_4`. To check it:

```
magma -b X1_42_check_orbs.m                    # permutations of the rational places
python3 aut_to_latex.py X1_42_aut_F5.txt --check
```

`magma -b X1_42_autgroup.m` recomputes the group (`Automorphisms(AFF)`, about 5,800 s) and overwrites `X1_42_aut_F5.txt`. Keep a copy of the stored file: all later steps use `sigma_4`, the 4th entry of this file.

**Step 3. Orbits of sigma_4 on the places of degree 1-5.**

```
./run_orbits.sh               # 19 Magma jobs (default parts 1 1 1 4 12); merges the parts at the end
magma -b X1_42_orbits_check.m
```

This writes `X1_42_orbits_d1.txt`, ..., `X1_42_orbits_d5.txt` and `X1_42_orbits_summary.txt`. The orbit counts must be 1, 1, 2, 17, 54 for degrees 1, ..., 5 (Appendix A.2). To reuse the stored orbit files instead, skip this step. Step 4 checks them against `sigma_4` in your session.

**Step 4. Session-independent labels.**

```
magma -b X1_42_rational_labels.m
```

This writes `X1_42_rational_labels.txt`:

* the signature labels of the 12 rational places;
* the permutation of `sigma_4` on them, which must be `S 12 9 7 3 6 10 1 2 5 11 4 8`, a 12-cycle;
* `Q1 5`;
* the lines `V 2 6`, `V 3 8`, `V 4 183`, `V 5 648`. These certify that the orbit files of degrees 2-5 agree with `sigma_4`, with the same convention, in the current session.

**Step 5. Splitting the search into jobs.**

```
python3 X1_42_split.py --njobs 20 --jobs9 2 --jobs8 2 --outdir X1_42_split_d11
```

The script enumerates, for every partition of 9, the parts `E_<=5` that are lexicographically minimal in their `sigma_4`-orbit. It checks each count against Burnside's lemma, and writes `plan.txt`, `summary.txt` and the job files `job<k>_of_20.txt`. `summary.txt` must give a total of 743,231 divisors (Table 3 of the paper). The paper's run used 2 jobs for [9], 2 for [8,1], 3 for largest part 7, 3 for largest part 6 and 10 for the rest. Compare the `GROUP ... JOBS` lines of `plan.txt` with `data/X1_42_split_d11/plan.txt`.

**Step 6. Running the jobs (about 577 CPU hours in total; the longest job took about 1.9 days).**

```
for k in $(seq 0 19); do
  nohup magma -b job:=$k njobs:=20 dir:=X1_42_split_d11 X1_42_split_search.m \
    > X1_42_split_d11/log_job${k}_of_20.txt 2>&1 &
done
wait
```

Each job recomputes the places it needs and their labels, and stops with an error unless every label it needs is found exactly once. It then tests its divisors and writes `status_job<k>_of_20.txt`. A degree-12 function that it finds goes to `hits_job<k>_of_20.txt`. A function of degree at most 11 would be written to `FOUND.txt`, which does not occur.

`./run_split_search.sh 20 --jobs9 2 --jobs8 2` performs Steps 5-6 and the check of Step 7 in one command.

**Step 7. Verification.**

```
python3 X1_42_split_check.py X1_42_split_d11 | tee X1_42_split_d11/split_check_output.txt
```

The last line must be the one shown in Route A.

### Diagnostic

`./run_placeorder.sh` runs `X1_42_placeorder_test.m` in several sessions and compares the order in which Magma returns the places. Its outputs `placeorder_run<r>.txt` show that this order is not reproducible. This is why the search identifies places by buckets, labels and fingerprints rather than by position.

## Citation

If you use this code, please cite the paper (reference to be added).

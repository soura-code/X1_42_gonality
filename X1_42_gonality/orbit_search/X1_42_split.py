#!/usr/bin/env python3
"""
X1_42_split.py  --  STAGE B of the order-independent degree-11 search on X1(42)/F5
------------------------------------------------------------------------------
Every effective divisor E of degree B = 9 is written as  E = E_low + Q  where
  * E_low is supported on places of degree <= 5, and
  * Q is either absent or ONE place of degree 6, 7, 8 or 9
    (two places of degree >= 6 would already have degree >= 12 > 9).

Orbit reduction (degrees 1-5).  G = <sigma_4> is cyclic of order 12.  E is kept
iff E_low is lexicographically minimal among sigma^k(E_low), k = 0..11.  The
family of all E is G-stable and the criterion only looks at E_low, so every
G-orbit of divisors keeps exactly one member (the proposition on canonical
representatives in X1_42_orbit_computation.tex; all of degrees 1-5 are reduced).  Q is not reduced.

Labels (independent of the Magma session):
  degree 1     : signature labels from X1_42_rational_labels.txt (+ sigma_4 on them)
  degrees 2-5  : rank of (bucket, key) in sorted order; the saved orbit files
                 X1_42_orbits_d<v>.txt are translated into these labels
  degrees 6-9  : not labelled at all.  A job with high degree v takes every
                 degree-v place Q whose bucket b satisfies BucketHash(b, s) = r
                 (job r of the s jobs of that group).  Every place has exactly one
                 bucket, so the s jobs together take every place exactly once,
                 whatever order Magma returns them in.

Groups and jobs:
  H9 [9], H8 [8,1], H7 (largest part 7), H6 (largest part 6): every job of the
  group gets ALL kept E_low of the group and a share of the high places.
  LOW (largest part <= 5): the kept divisors are dealt out round-robin.
  Jobs are allocated to groups in proportion to their (weighted) workload.

Usage:
  python3 X1_42_split.py --njobs 20 [--jobs9 K] [--jobs8 K] [--w9 0.72] [--outdir DIR]
Needs X1_42_rational_labels.txt and X1_42_orbits_d2.txt ... X1_42_orbits_d5.txt.
Output (DIR, default X1_42_split_d11/): job<k>_of_<N>.txt, plan.txt, summary.txt
"""
import argparse
import itertools
import math
import os
import sys
from collections import defaultdict

P = {1: 12, 2: 6, 3: 8, 4: 183, 5: 648, 6: 2288, 7: 11160, 8: 48708, 9: 214384}
DTARGET, B, G_ORDER = 11, 9, 12
LOWDEG = (1, 2, 3, 4, 5)


def die(msg):
    sys.exit("ERROR: " + msg)


def enc_poly(coeffs):                         # = SL_EncPoly (constant term first)
    e = 0
    for c in reversed(coeffs):
        e = 5 * e + c
    return e + 5 ** len(coeffs)


def cycle_lengths(perm):
    seen, out = set(), []
    for i in range(len(perm)):
        if i in seen:
            continue
        L, x = 0, i
        while x not in seen:
            seen.add(x); x = perm[x]; L += 1
        out.append(L)
    return sorted(out)


def read_rational(path):
    try:
        lines = open(path).read().split("\n")
    except OSError:
        die("%s not found: run  magma -b X1_42_rational_labels.m  first" % path)
    if not lines[0].startswith("X1_42_RATLABELS v2"):
        die("bad header in " + path)
    if "GEN 4" not in lines[0]:
        die(path + " was not made with generator sigma_4")
    perm, verified = None, {}
    for l in lines:
        t = l.split()
        if t and t[0] == "S":
            perm = [int(x) - 1 for x in t[1:]]
        if t and t[0] == "V":
            verified[int(t[1])] = int(t[2])
    for v in (2, 3, 4, 5):
        if verified.get(v) != P[v]:
            die("%s does not certify the orbit file of degree %d (rerun X1_42_rational_labels.m)" % (path, v))
    if perm is None or sorted(perm) != list(range(12)):
        die("no valid permutation in " + path)
    if "END" not in lines:
        die(path + " incomplete")
    if cycle_lengths(perm) != [12]:
        die("sigma_4 is not a 12-cycle on the rational labels")
    return perm


def read_orbit(v):
    """returns (perm on labels 0-based, ids[label] = (bucket, key))"""
    name = "X1_42_orbits_d%d.txt" % v
    try:
        lines = open(name).read().split("\n")
    except OSError:
        die(name + " not found")
    h = lines[0].split()
    if h[:2] != ["X1_42_ORBITS", "v1"]:
        die("bad header in " + name)
    hd = dict(zip(h[2::2], h[3::2]))
    if hd.get("GEN") != "4":
        die("%s was not made with generator sigma_4" % name)
    n = int(hd["COUNT"])
    if n != P[v]:
        die("%s has %d places, expected %d" % (name, n, P[v]))
    data, img = {}, {}
    for l in lines[1:]:
        if l.startswith("P "):
            parts = l.split("|")
            i = int(parts[0].split()[1])
            coeffs = [int(z) for z in parts[1].split()]
            key = int(parts[2].split()[0])
            data[i] = (enc_poly(coeffs) if coeffs else 0, key)
        elif l.startswith("I "):
            _, i, j = l.split()
            img[int(i)] = int(j)
    if sorted(data) != list(range(1, n + 1)) or sorted(img) != list(range(1, n + 1)):
        die(name + " incomplete")
    ids = sorted(data.values())
    if len(set(ids)) != n:
        die("degree %d: (bucket, key) not unique; labels by (bucket, key) impossible" % v)
    lab = {i: ids.index(data[i]) for i in data}          # file index -> label (0-based)
    perm = [None] * n
    for i in data:
        perm[lab[i]] = lab[img[i]]
    if sorted(perm) != list(range(n)) or any(G_ORDER % L for L in cycle_lengths(perm)):
        die("degree %d: translated permutation invalid" % v)
    return perm, ids


def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield []
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield [k] + rest


def burnside(low, perms):
    """number of G-orbits on divisors of type `low`, from the cycle types only"""
    mult = defaultdict(int)
    for v in low:
        mult[v] += 1
    if not mult:
        return 1
    tot = 0
    for k in range(G_ORDER):
        prod = 1
        for v, m in mult.items():
            cl = []
            for L in cycle_lengths(perms[v]):
                g = math.gcd(L, k) if k else L
                cl += [L // g] * g
            a = [1] + [0] * m                      # fixed multisets of size m
            for c in cl:
                for i in range(c, m + 1):
                    a[i] += a[i - c]
            prod *= a[m]
        tot += prod
    assert tot % G_ORDER == 0
    return tot // G_ORDER


def kept_low(low, powers):
    """kept E_low of type `low` (list of degrees <= 5); E_low = tuple of (v, labels-tuple)"""
    mult = defaultdict(int)
    for v in low:
        mult[v] += 1
    degs = sorted(mult, reverse=True)
    if not degs:
        return [()]
    out = []
    for E in itertools.product(*[itertools.combinations_with_replacement(range(P[v]), mult[v]) for v in degs]):
        ok = True
        for g in range(1, G_ORDER):
            img = tuple(tuple(sorted(powers[v][g][i] for i in E[a])) for a, v in enumerate(degs))
            if img < E:
                ok = False
                break
        if ok:
            out.append(tuple(zip(degs, E)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--njobs", type=int, default=20)
    ap.add_argument("--jobs9", type=int, default=0, help="jobs for [9] (default: by workload)")
    ap.add_argument("--jobs8", type=int, default=0, help="jobs for [8,1] (default: by workload)")
    ap.add_argument("--w9", type=float, default=0.72,
                    help="relative cost of a [9] divisor (old run: 2.94 s vs 4.08 s average) ")
    ap.add_argument("--outdir", default="X1_42_split_d11")
    a = ap.parse_args()
    N = a.njobs
    if N < 5:
        die("need at least 5 jobs (one per group)")

    perms = {1: read_rational("X1_42_rational_labels.txt")}
    ids = {}
    for v in (2, 3, 4, 5):
        perms[v], ids[v] = read_orbit(v)
    for v in LOWDEG:
        print("degree %d: %d places, cycle type of sigma %s" % (v, P[v], cycle_lengths(perms[v])))
    powers = {}
    for v, p in perms.items():
        pw = [list(range(len(p)))]
        for _ in range(G_ORDER - 1):
            pw.append([p[x] for x in pw[-1]])
        powers[v] = pw

    # ---- enumerate
    groups = {g: [] for g in ("H9", "H8", "H7", "H6", "LOW")}    # list of (lam, E_low)
    rows = []
    for lam in partitions(B):
        high = lam[0] if lam[0] >= 6 else 0
        low = lam[1:] if high else lam
        reps = kept_low(low, powers)
        if len(reps) != burnside(low, perms):
            die("partition %s: %d kept low parts but %d orbits (Burnside)" % (lam, len(reps), burnside(low, perms)))
        g = "H%d" % high if high else "LOW"
        groups[g].extend((lam, E) for E in reps)
        units = len(reps) * (P[high] if high else 1)
        allcount = math.prod(math.comb(P[v] + m - 1, m) for v, m in
                             ((v, lam.count(v)) for v in set(lam)))
        rows.append((lam, g, len(reps), units, allcount))
        print("  %-22s %-4s kept low parts %8d   divisors %9d   (all %9d)" % (str(lam), g, len(reps), units, allcount),
              flush=True)
    units = {g: (len(groups[g]) * P[int(g[1])] if g != "LOW" else len(groups[g])) for g in groups}
    total = sum(units.values())
    print("TOTAL divisors to test: %d   (all divisors: %d)" % (total, sum(r[4] for r in rows)))
    print("check: for every partition, kept low parts = number of orbits (Burnside)")

    # ---- allocate jobs (largest remainder on weighted workload, >= 1 per group)
    weight = {g: units[g] * (a.w9 if g == "H9" else 1.0) for g in units}
    fixed = {}
    if a.jobs9:
        fixed["H9"] = a.jobs9
    if a.jobs8:
        fixed["H8"] = a.jobs8
    free = [g for g in units if g not in fixed]
    left = N - sum(fixed.values())
    if left < len(free):
        die("not enough jobs")
    # minimise the largest per-job load: start with one job per group, then
    # repeatedly give a job to the group whose jobs carry the most work
    share = {g: 1 for g in free}
    for _ in range(left - len(free)):
        g = max(free, key=lambda x: weight[x] / share[x])
        share[g] += 1
    alloc = dict(fixed); alloc.update(share)

    # ---- write job files
    os.makedirs(a.outdir, exist_ok=True)
    for f in os.listdir(a.outdir):
        if f.startswith("job") or f.startswith("status_job") or f.startswith("hits_job") or f == "FOUND.txt":
            os.remove(os.path.join(a.outdir, f))
    order = ["H9", "H8", "H7", "H6", "LOW"]
    k = 0
    plan = []
    for g in order:
        s = alloc[g]
        high = int(g[1]) if g != "LOW" else 0
        for r in range(s):
            items = groups[g] if high else groups[g][r::s]
            used = defaultdict(set)
            for lam, E in items:
                for v, labs in E:
                    used[v].update(labs)
            name = os.path.join(a.outdir, "job%d_of_%d.txt" % (k, N))
            with open(name, "w") as F:
                F.write("SPLITJOB v1 DTARGET %d B %d JOB %d NJOBS %d GROUP %s RANK %d GSIZE %d HIGHDEG %d NLOW %d\n"
                        % (DTARGET, B, k, N, g, r, s, high, len(items)))
                for v in sorted(set(list(used)) | ({high} if high else set())):
                    F.write("C %d %d\n" % (v, P[v]))
                for v in sorted(used):
                    if v == 1:
                        continue
                    for lab in sorted(used[v]):
                        b, key = ids[v][lab]
                        F.write("K %d %d %d %d\n" % (v, lab + 1, b, key))
                last = None
                for lam, E in items:
                    if lam != last:
                        F.write("L %s\n" % " ".join(map(str, lam)))
                        last = lam
                    F.write("D" + "".join(" %d:%d" % (v, lab + 1) for v, labs in E for lab in labs) + "\n")
                F.write("END %d\n" % len(items))
            plan.append((k, g, r, s, high, len(items)))
            k += 1
    with open(os.path.join(a.outdir, "plan.txt"), "w") as F:
        F.write("NJOBS %d TOTAL %d\n" % (N, total))
        for g in order:
            high = int(g[1]) if g != "LOW" else 0
            F.write("GROUP %s JOBS %d HIGHDEG %d NLOW %d PHIGH %d UNITS %d\n"
                    % (g, alloc[g], high, len(groups[g]), P[high] if high else 0, units[g]))
        for (k, g, r, s, high, n) in plan:
            F.write("JOB %d GROUP %s RANK %d GSIZE %d HIGHDEG %d NLOW %d\n" % (k, g, r, s, high, n))
    with open(os.path.join(a.outdir, "summary.txt"), "w") as F:
        F.write("Order-independent degree-11 search, B = 9, orbit reduction on degrees 1-5\n")
        F.write("%-22s %-4s %12s %12s %12s\n" % ("partition", "grp", "kept E_low", "divisors", "all"))
        for lam, g, nl, u, al in rows:
            F.write("%-22s %-4s %12d %12d %12d\n" % (str(lam), g, nl, u, al))
        F.write("TOTAL %d divisors to test (all: %d)\n" % (total, sum(r[4] for r in rows)))
        for g in order:
            F.write("group %-4s jobs %2d  divisors %9d  (weighted per job %.0f)\n"
                    % (g, alloc[g], units[g], weight[g] / alloc[g]))
    for g in order:
        print("group %-4s: %2d jobs, %9d divisors, weighted per job %9.0f" % (g, alloc[g], units[g], weight[g] / alloc[g]))
    print("wrote %d job files, plan.txt and summary.txt in %s" % (len(plan), a.outdir))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
X1_42_orbits_merge.py  --  combine the parts written by X1_42_orbits.m.

For every degree v with part files X1_42_orbits_d<v>_part<k>_of_<m>.txt it
  * checks that all m parts are present and complete, and agree on the number
    of places, the generator and the checksum of the place data (so all jobs
    numbered the places in the same way),
  * checks that every place has exactly one image and none is ambiguous (-1),
  * checks that the images form a permutation whose cycle lengths divide 12,
  * for v = 1, checks that it is the permutation of sigma_4 found by
    X1_42_check_orbs.m (or its inverse),
  * splits the places into orbits (= cycles, since <sigma_4> is the whole group)
    and writes X1_42_orbits_d<v>.txt (the format read by X1_42_orbits_load.m).
A summary is written to X1_42_orbits_summary.txt.  A degree that fails any
check is NOT written.

Usage:  python3 X1_42_orbits_merge.py
"""
import glob
import os
import re
import sys
from collections import Counter, defaultdict

A4KNOWN = [9, 5, 2, 1, 12, 8, 11, 7, 3, 4, 10, 6]   # sigma_4 on Places(AFF,1), from X1_42_check_orbs.m

def merge_degree(v):
    files = glob.glob("X1_42_orbits_d%d_part*_of_*.txt" % v)
    ms = {int(re.search(r"_of_(\d+)\.txt$", f).group(1)) for f in files}
    if len(ms) != 1:
        return None, "degree %d: part files with different numbers of parts %s (delete the old ones)" % (v, sorted(ms))
    m = ms.pop()
    heads, Plines, img = [], [], {}
    for k in range(m):
        name = "X1_42_orbits_d%d_part%d_of_%d.txt" % (v, k, m)
        try:
            lines = open(name).read().split("\n")
        except OSError:
            return None, "degree %d: part %d of %d missing (%s)" % (v, k, m, name)
        h = lines[0].split()
        if h[:2] != ["X1_42_ORBITPART", "v1"]:
            return None, "degree %d: bad header in %s" % (v, name)
        hd = dict(zip(h[2::2], h[3::2]))
        if not any(l.startswith("END") for l in lines):
            return None, "degree %d: %s incomplete (job not finished?)" % (v, name)
        heads.append(hd)
        for l in lines[1:]:
            if l.startswith("P "):
                if k == 0:
                    Plines.append(l)
            elif l.startswith("I "):
                _, i, j = l.split()
                i, j = int(i), int(j)
                if i in img:
                    return None, "degree %d: image of place %d given twice" % (v, i)
                if (i - 1) % m != k:
                    return None, "degree %d: place %d found in the wrong part %d" % (v, i, k)
                img[i] = j
    for f in ("COUNT", "GEN", "CHECKSUM"):
        if len({hd[f] for hd in heads}) != 1:
            return None, "degree %d: parts disagree on %s (different place order?)" % (v, f)
    n = int(heads[0]["COUNT"])
    gen = heads[0]["GEN"]
    chk = heads[0]["CHECKSUM"]
    if len(Plines) != n:
        return None, "degree %d: place data incomplete (%d of %d)" % (v, len(Plines), n)
    if sorted(img) != list(range(1, n + 1)):
        return None, "degree %d: images missing for %d places" % (v, n - len(img))
    amb = [i for i in img if img[i] == -1]
    if amb:
        return None, "degree %d: %d ambiguous images (e.g. place %d)" % (v, len(amb), amb[0])
    if sorted(img.values()) != list(range(1, n + 1)):
        return None, "degree %d: the images do not form a permutation" % v
    orbits, seen = [], set()
    for i in range(1, n + 1):
        if i in seen:
            continue
        cyc, z = [], i
        while True:
            cyc.append(z); seen.add(z); z = img[z]
            if z == i:
                break
        orbits.append(cyc)
    lens = [len(o) for o in orbits]
    if any(12 % L for L in lens):
        return None, "degree %d: cycle lengths %s do not all divide 12" % (v, sorted(set(lens)))
    note = ""
    if v == 1:
        mine = [img[i] for i in range(1, 13)]
        inv = [A4KNOWN.index(j) + 1 for j in range(1, 13)]
        if mine == A4KNOWN:
            note = "; agrees with sigma_4 from X1_42_check_orbs.m"
        elif mine == inv:
            note = "; agrees with the inverse of sigma_4 from X1_42_check_orbs.m (same group)"
        else:
            return None, "degree 1: permutation %s does not agree with X1_42_check_orbs.m" % mine
    # places whose (polynomial, key) is not unique, or have no key: only findable by index
    data = []
    for l in Plines:
        parts = l.split("|")
        data.append((parts[1].strip(), parts[2].strip()))
    cnt = Counter(data)
    nonunique = sum(1 for d in data if cnt[d] > 1 or d[1] == "0")

    out = "X1_42_orbits_d%d.txt" % v
    with open(out, "w") as F:
        F.write("X1_42_ORBITS v1 DEG %d GEN %s COUNT %d ORBITS %d CHECKSUM %s\n" % (v, gen, n, len(orbits), chk))
        F.write("# Places of degree %d of X1(42)/F5, numbered as in Magma's Places(AFF,%d).\n" % (v, v))
        F.write("# P i | coefficients (constant term first) of the polynomial of F5(t) below the place\n")
        F.write("#       (empty = t = infinity) | key (0 = not available)\n")
        F.write("# I i j : sigma_%s maps place i to place j\n" % gen)
        F.write("# O k size : i_1 i_2 ... (orbit k as a cycle i, sigma(i), sigma^2(i), ...)\n")
        F.write("# R k i : representative of orbit k (its smallest index)\n")
        for l in sorted(Plines, key=lambda s: int(s.split()[1])):
            F.write(l + "\n")
        for i in range(1, n + 1):
            F.write("I %d %d\n" % (i, img[i]))
        for k, o in enumerate(orbits, 1):
            F.write("O %d %d : %s\n" % (k, len(o), " ".join(map(str, o))))
        for k, o in enumerate(orbits, 1):
            F.write("R %d %d\n" % (k, min(o)))
        F.write("END %d\n" % n)
    dist = ", ".join("%d x size %d" % (c, s) for s, c in sorted(Counter(lens).items()))
    return out, ("degree %d: %d places, %d orbits (%s), merged from %d parts, %d places identifiable only by index%s; wrote %s"
                 % (v, n, len(orbits), dist, m, nonunique, note, out))

def main():
    degs = sorted({int(re.search(r"_d(\d+)_part", f).group(1)) for f in glob.glob("X1_42_orbits_d*_part*_of_*.txt")})
    if not degs:
        sys.exit("no part files X1_42_orbits_d<v>_part<k>_of_<m>.txt found")
    lines, ok_all = [], True
    for v in degs:
        out, msg = merge_degree(v)
        if out is None:
            ok_all = False
            stale = "X1_42_orbits_d%d.txt" % v
            if os.path.exists(stale):
                os.remove(stale)                       # never leave an old file for a failed degree
                msg += " (old %s removed)" % stale
            msg = "FAILED " + msg + " -- not written"
        print(msg)
        lines.append(msg)
    with open("X1_42_orbits_summary.txt", "w") as F:
        F.write("Orbits of Aut_F5(X1(42)) = <sigma_4> on places of small degree\n")
        for l in lines:
            F.write(l + "\n")
    print("summary written to X1_42_orbits_summary.txt")
    sys.exit(0 if ok_all else 1)

if __name__ == "__main__":
    main()

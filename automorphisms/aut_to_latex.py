#!/usr/bin/env python3
"""
aut_to_latex.py -- turn the automorphism data of X1(42)/F5 into a compact,
exact LaTeX listing for the appendix of the manuscript.

Each automorphism sigma is stored as sigma(t), sigma(y) in
    F5(X) = F5(t)[y]/(F),   [F5(X):F5(t)] = 16,
written as  sum_{j=0}^{15} c_j(t) y^j  with c_j in F5(t).

Every denominator that occurs is (after cancelling a common power of t) of
the form  t^a (t+4)^b (t^2+4t+1)^c,  so each coefficient is printed as the
exponents (a,b,c) together with the coefficient list of its numerator N_j(t),
highest degree first.  The script checks this factorisation exactly; if some
denominator does not factor this way it stops with an error.

Input (either):
    X1_42_aut_F5.txt            raw data file written by the Magma scripts
                                (format v1 or v2)
    --from-tex FILE.tex         a listing produced by the previous version of
                                this script (\\dfrac format), parsed back and
                                re-rendered

Usage:
    python3 aut_to_latex.py X1_42_aut_F5.txt > aut_appendix.tex
    python3 aut_to_latex.py X1_42_aut_F5.txt --only 4 > aut_appendix.tex
    python3 aut_to_latex.py --from-tex old_listing.tex > aut_appendix.tex
    python3 aut_to_latex.py X1_42_aut_F5.txt --check
"""
import re
import sys

P = 5
NDEG = 16

# ----------------------------------------------------------- F5[t] arithmetic
def trim(c):
    c = [x % P for x in c]
    while c and c[-1] == 0:
        c.pop()
    return c

def pmul(a, b):
    if not a or not b:
        return []
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] = (r[i + j] + x * y) % P
    return trim(r)

def pdivmod(a, b):
    a, b = trim(a), trim(b)
    inv = pow(b[-1], P - 2, P)
    q = [0] * max(len(a) - len(b) + 1, 1)
    r = a[:]
    while len(r) >= len(b) and r:
        k = len(r) - len(b)
        f = (r[-1] * inv) % P
        q[k] = f
        for i, y in enumerate(b):
            r[i + k] = (r[i + k] - f * y) % P
        r = trim(r)
    return trim(q), r

def low_order(c):
    for k, x in enumerate(c):
        if x % P:
            return k
    return None

T = [0, 1]                 # t
L1 = [4, 1]                # t + 4  (= t - 1)
Q2 = [1, 4, 1]             # t^2 + 4t + 1  (= t^2 - t + 1)

def split_den(den):
    """den = t^a (t+4)^b (t^2+4t+1)^c exactly (monic); return (a,b,c)."""
    d = trim(den)
    exps = []
    for f in (T, L1, Q2):
        e = 0
        while True:
            q, r = pdivmod(d, f)
            if r:
                break
            d, e = q, e + 1
        exps.append(e)
    if d != [1]:
        raise SystemExit("denominator has an unexpected factor: %r" % d)
    return tuple(exps)

def normalise(num, den):
    """Cancel the common power of t, make the denominator monic;
    return (num, (a,b,c)), or None if the coefficient is zero."""
    num, den = trim(num), trim(den)
    if not num:
        return None
    k = min(low_order(num), low_order(den))
    num, den = num[k:], den[k:]
    inv = pow(den[-1], P - 2, P)
    num = trim([x * inv for x in num])
    den = trim([x * inv for x in den])
    abc = split_den(den)
    back = [1]                                  # exact check
    for f, e in zip((T, L1, Q2), abc):
        for _ in range(e):
            back = pmul(back, f)
    assert back == den
    return num, abc

# ----------------------------------------------------------- input: raw file
def tokens(path):
    txt = open(path).read().replace("\\\n", "")
    out = []
    for t in txt.split():
        if t == "[]":
            continue
        if t == "[]/":
            out.append("/")
            continue
        out.append(t)
    return out

def parse_raw(path):
    tk = tokens(path)
    if tk[0] != "X1_42_AUT_F5" or tk[-1] != "END":
        raise SystemExit("not an X1_42_AUT_F5 file (or incomplete: no END)")
    pos = 2
    assert tk[pos] == "DEG"; ndeg = int(tk[pos + 1]); pos += 2
    assert tk[pos] == "FINGERPRINT"; pos += 1
    while tk[pos] not in ("ENDFP", "COUNT"):
        pos += 1
    if tk[pos] == "ENDFP":
        pos += 1
    assert tk[pos] == "COUNT"; m = int(tk[pos + 1]); pos += 2
    raw = {}
    while tk[pos] != "END":
        if tk[pos] == ";":
            pos += 1; continue
        assert tk[pos] == "C" and tk[pos + 4] == ":", "bad record at token %d" % pos
        i, w, j = int(tk[pos + 1]), tk[pos + 2], int(tk[pos + 3]); pos += 5
        num = []
        while tk[pos] != "/":
            num.append(int(tk[pos])); pos += 1
        pos += 1
        den = []
        while tk[pos] not in (";", "C", "END"):
            den.append(int(tk[pos])); pos += 1
        raw[(i, w, j)] = (num, den)
    data = {}
    for i in range(1, m + 1):
        for w in "TY":
            co = {}
            for j in range(1, ndeg + 1):
                if (i, w, j) not in raw:
                    raise SystemExit("record missing: %r" % ((i, w, j),))
                co[j - 1] = raw[(i, w, j)]
            data[(i, w)] = co
    return m, data

# ----------------------------------------------------------- input: old .tex
def parse_poly(s):
    c = {}
    for term in s.strip().split(" + "):
        mm = re.fullmatch(r"(\d*)(t(\^\{(\d+)\})?)?", term.strip())
        if not mm:
            raise SystemExit("cannot parse term %r" % term)
        a = int(mm.group(1)) if mm.group(1) else 1
        k = 0 if not mm.group(2) else (int(mm.group(4)) if mm.group(4) else 1)
        c[k] = a
    return [c.get(i, 0) for i in range(max(c) + 1)]

def parse_tex(path):
    txt = open(path).read()
    data = {}
    for mt in re.finditer(r"\\sigma_\{(\d+)\}\((t|y)\) = (.*?)\n\\end\{dmath\*\}", txt, re.S):
        i, w, body = int(mt.group(1)), mt.group(2).upper(), mt.group(3)
        co = {}
        for piece in body.split(r" \;+\; "):
            piece = piece.strip()
            j = 0
            my = re.fullmatch(r"(.*)\\,y(\^\{(\d+)\})?", piece, re.S)
            if my:
                j = int(my.group(3)) if my.group(3) else 1
                piece = my.group(1)
                if piece.startswith(r"\left("):
                    piece = piece[6:-7]
            elif piece == "y":
                j, piece = 1, "1"
            elif re.fullmatch(r"y\^\{\d+\}", piece):
                j, piece = int(piece[3:-1]), "1"
            mf = re.fullmatch(r"\\dfrac\{(.*)\}\{(.*)\}", piece, re.S)
            num, den = ((parse_poly(mf.group(1)), parse_poly(mf.group(2)))
                        if mf else (parse_poly(piece), [1]))
            if j in co:
                raise SystemExit("duplicate y^%d in sigma_%d(%s)" % (j, i, w.lower()))
            co[j] = (num, den)
        data[(i, w)] = co
    m = max(i for i, _ in data)
    for i in range(1, m + 1):
        for w in "TY":
            if (i, w) not in data:
                raise SystemExit("sigma_%d(%s) missing in %s" % (i, w.lower(), path))
    return m, data

# ----------------------------------------------------------- output
def coeffs(data, i, w):
    return [normalise(*data[(i, w)].get(j, ([], [1]))) for j in range(NDEG)]

def is_identity(data, i):
    ct, cy = coeffs(data, i, "T"), coeffs(data, i, "Y")
    return (ct == [([0, 1], (0, 0, 0))] + [None] * (NDEG - 1)
            and cy == [None, ([1], (0, 0, 0))] + [None] * (NDEG - 2))

def coeff_string(num):
    """coefficients of num, highest degree first"""
    return r"\,".join(str(x) for x in reversed(num))

def table(i, w, data):
    name = "t" if w == "T" else "y"
    rows = []
    for j, nd in enumerate(coeffs(data, i, w)):
        if nd is None:
            rows.append(r"$%d$ & & & & & $0$\\" % j)
        else:
            num, (a, b, c) = nd
            rows.append(r"$%d$ & $%d$ & $%d$ & $%d$ & $%d$ & $%s$\\"
                        % (j, a, b, c, len(num) - 1, coeff_string(num)))
    return "\n".join([
        r"\begin{center}\scriptsize",
        r"$\sigma_{%d}(%s)=\displaystyle\sum_{j=0}^{15}\frac{N_j(t)}{t^{a_j}(t+4)^{b_j}(t^2+4t+1)^{c_j}}\,y^j$"
        % (i, name),
        r"\\[4pt]",
        r"\begin{tabular}{@{}rrrrrl@{}}",
        r"\toprule",
        r"$j$ & $a_j$ & $b_j$ & $c_j$ & $\deg N_j$ & coefficients of $N_j$, from $t^{\deg N_j}$ down to $t^0$\\",
        r"\midrule",
        *rows,
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{center}",
    ])

def render(data, ids, src):
    out = ["% generated by aut_to_latex.py from " + src + " -- do not edit by hand"]
    for i in ids:
        out.append(r"\subsection*{The automorphism $\sigma_{%d}$}" % i)
        if is_identity(data, i):
            out.append(r"$\sigma_{%d}(t)=t$ and $\sigma_{%d}(y)=y$." % (i, i))
            continue
        out.append(table(i, "T", data))
        out.append(table(i, "Y", data))
    return "\n".join(out) + "\n"

def main():
    args = sys.argv[1:]
    if not args:
        raise SystemExit(__doc__)
    if "--from-tex" in args:
        src = args[args.index("--from-tex") + 1]
        m, data = parse_tex(src)
    else:
        src = args[0]
        m, data = parse_raw(src)
    for i in range(1, m + 1):                    # validate everything first
        for w in "TY":
            coeffs(data, i, w)
    if "--check" in args:
        print("OK: %d automorphisms; every denominator is t^a (t+4)^b (t^2+4t+1)^c" % m)
        return
    ids = list(range(1, m + 1))
    if "--only" in args:
        ids = [int(x) for x in args[args.index("--only") + 1].split(",")]
    sys.stdout.write(render(data, ids, src))

if __name__ == "__main__":
    main()

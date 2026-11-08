/*
  X1_42_split_lib.m  --  shared functions of the order-independent degree-11 search
  ------------------------------------------------------------------
  `load` after AFF is defined.  Nothing here depends on the order in which
  Magma returns places: every function computes data of a single place.

    PlaceBucket(Q)        code of the place of F5(t) below Q:
                          Enc(min. poly of t(Q)), or 0 if t has a pole at Q.
                          (Used for degrees 6..9; no key is computed there.)
    PlaceBucketKey(Q)     bucket and key, exactly as in X1_42_orbits.m
                          (used only for degrees 2..5, to match the saved orbits).
    BucketHash(b, s)      which of s jobs handles bucket b.
    SigFunctions(nb)      fixed list of explicit functions (batches 1..nb) used to
                          label the 12 rational places.
    SignatureString(Q, L) the string  "v:r;" for every g in L, where v = v_Q(g) and
                          r = g(Q) if v = 0 (else "x").  Only for rational places.
    FunToLine(g), LineToFun(s)   write / read a function of AFF as explicit
                          coefficient data (numerators and denominators in F5[t]
                          of its 16 coordinates), so that every session uses
                          exactly the same functions.
  ------------------------------------------------------------------
*/
SL_TT := AFF!(BaseField(AFF).1);
SL_YY := AFF.1;

SL_EncPoly := function(g)
    c := Coefficients(g); e := 0;
    for i in [#c..1 by -1] do e := 5*e + (Integers()!c[i]); end for;
    return e + 5^(#c);
end function;
SL_ValueAt := function(g, Q)
    try
        if Valuation(g, Q) lt 0 then return false, 0; end if;
        R := ResidueClassField(Q);
        return true, R!Evaluate(g, Q);
    catch e
        return false, 0;
    end try;
end function;
SL_MinPoly := function(a)
    k := PrimeField(Parent(a));
    try
        return MinimalPolynomial(a, k);
    catch e1
        dummy := 0;
    end try;
    return MinimalPolynomial(a);
end function;
SL_KeyFromValues := function(a, b)
    try
        polys := [ SL_MinPoly(a), SL_MinPoly(b), SL_MinPoly(a+b), SL_MinPoly(a+2*b) ];
    catch e
        return false, 0;
    end try;
    E := 5^12; key := 0;
    for g in polys do key := key*E + SL_EncPoly(g); end for;
    return true, key;
end function;

// STRICT: no error is caught here.  If the residue of t cannot be computed the
// job stops with an error instead of silently putting the place into bucket 0
// (which could move it into another job's class than in other sessions).
PlaceBucket := function(Q)
    if Valuation(SL_TT, Q) lt 0 then return 0; end if;
    a := ResidueClassField(Q)!Evaluate(SL_TT, Q);
    return SL_EncPoly(MinimalPolynomial(a, GF(5)));
end function;
SL_FPMOD := 2305843009213693951;
// fingerprint of a bucket (used to check that all jobs of a group saw the same buckets)
BucketPrint := func< b | ((b mod SL_FPMOD) * 1000003 + 7) mod SL_FPMOD >;

// identical to PlaceData of X1_42_orbits.m (bucket, key); key 0 = not available
PlaceBucketKey := function(Q)
    okt, at := SL_ValueAt(SL_TT, Q);
    if not okt then return 0, 0; end if;
    b := SL_EncPoly(SL_MinPoly(at));
    k := 0;
    oky, ay := SL_ValueAt(SL_YY, Q);
    if oky then
        ok, kk := SL_KeyFromValues(at, ay);
        if ok then k := kk; end if;
    end if;
    return b, k;
end function;

BucketHash := function(b, s)
    return (((b mod 1000000007) * 1000003 + 12345) mod 1000000007) mod s;
end function;

// ---- labels of the rational places ----
SigFunctions := function(nb)
    T := SL_TT; Y := SL_YY;
    L := [ T, Y ];
    if nb ge 2 then L cat:= [ T - a : a in [1..4] ] cat [ Y - c : c in [1..4] ]; end if;
    if nb ge 3 then L cat:= [ (Y - c)/(T - a) : a in [0..4], c in [0..4] ]; end if;
    if nb ge 4 then L cat:= [ (T - a)/(Y - c) : a in [0..4], c in [0..4] ]; end if;
    if nb ge 5 then
        L cat:= [ (Y - c)/(T - a)^2 : a in [0..4], c in [0..4] ]
              cat [ (T - a)/(Y - c)^2 : a in [0..4], c in [0..4] ];
    end if;
    if nb ge 6 then
        L cat:= [ (Y - c - s*(T - a))/(T - a)^2 : a in [0..4], c in [0..4], s in [0..4] ];
    end if;
    return L;
end function;
SIG_NBMAX := 6;

SL_Entry := function(g, Q)
    v := Valuation(g, Q);
    if v ne 0 then return Sprintf("%o:x;", v); end if;
    return Sprintf("0:%o;", Integers()!(GF(5)!Evaluate(g, Q)));
end function;
SignatureString := function(Q, L)
    s := "";
    for g in L do s cat:= SL_Entry(g, Q); end for;
    return s;
end function;

// ---- explicit functions as text ----
SL_PR := PolynomialRing(GF(5));
FunToLine := function(g)
    parts := [];
    for c in Eltseq(g) do
        nu := [ Integers()!z : z in Coefficients(Numerator(c)) ];
        de := [ Integers()!z : z in Coefficients(Denominator(c)) ];
        ns := ""; for z in nu do ns cat:= Sprintf(" %o", z); end for;
        ds := ""; for z in de do ds cat:= Sprintf(" %o", z); end for;
        Append(~parts, ns cat " ;" cat ds);
    end for;
    line := parts[1];
    for i in [2..#parts] do line cat:= " |" cat parts[i]; end for;
    return line;
end function;
LineToFun := function(s)
    R := BaseField(AFF);
    cs := [];
    for part in Split(s, "|") do
        nd := Split(part, ";");
        nu := [ StringToInteger(z) : z in Split(nd[1], " ") | #z gt 0 ];
        if #nd ge 2 then de := [ StringToInteger(z) : z in Split(nd[2], " ") | #z gt 0 ]; else de := [1]; end if;
        Append(~cs, Evaluate(SL_PR!nu, R.1) / Evaluate(SL_PR!de, R.1));
    end for;
    return AFF!cs;
end function;

/*
  X1_42_rational_labels.m  --  ONE-TIME: canonical labels of the 12 rational
  places of X1(42)/F5 and the permutation of sigma_4 on these labels.
  ------------------------------------------------------------------
  The order of Places(AFF,1) is not guaranteed to be the same in every Magma
  session, and 10 of the 12 rational places share their (bucket, key).  So the
  rational places are labelled by a SIGNATURE: the pairs (v_P(g), g(P) if
  v_P(g) = 0) for a list of functions g that separates them.  The list is
  t, y, a local uniformizer at each rational place and a basis of L(3S); it is
  computed ONCE here and written to the labels file as explicit coefficient
  data (lines "F ..."), so every job reloads exactly the same functions and
  gets the same labels.  label = rank of the signature string in sorted order.

  sigma = sigma_4 (from the saved group, not recomputed) acts on the labels by
      signature of sigma(P) = [ (v_P(sigma^* g), (sigma^* g)(P)) : g ],
  because v_{sigma(P)}(g) = v_P(sigma^* g) and residues correspond (Lemma 2.1
  of X1_42_orbit_computation.tex; for rational places both residues lie in F5).
  Checks: signatures distinct, the images form a permutation, it is a 12-cycle.

  It ALSO VERIFIES the saved orbit files X1_42_orbits_d2.txt ... _d5.txt against
  the same generator: for every place P of degree 2..5 the (bucket, key) of
  sigma(P) is computed from the residues of sigma^*t, sigma^*y at P (Lemmas
  3.2/3.3 of X1_42_orbit_computation.tex) and compared with the (bucket, key)
  of the image stored in the file.  Since (bucket, key) is unique in these
  degrees, agreement for all places proves that the stored permutations are
  exactly sigma_4 -- in particular with the SAME convention (sigma, not
  sigma^-1) as the degree-1 permutation computed here.  The orbit reduction
  needs this: on all degrees the permutations must come from one automorphism.

  Output: X1_42_rational_labels.txt
      X1_42_RATLABELS v2 GEN 4
      NF n                        number of separating functions
      F i data                    function i as explicit data (FunToLine)
      R label signature           (12 lines)
      S img(1) ... img(12)        sigma_4 on the labels
      Q1 label                    the zero of order 3 of y (for reference)
      V v count                   orbit file of degree v verified (v = 2..5)
      END
  Usage:  magma -b X1_42_rational_labels.m    (needs X1_42_split_lib.m,
          X1_42_aut_load.m, X1_42_aut_F5.txt, X1_42_orbits_load.m,
          X1_42_orbits_d2.txt ... X1_42_orbits_d5.txt)
*/
SetColumns(0);
A<x,y> := PolynomialRing(GF(5),2);
P,homogenize := Homogenization(A);
f := x^16*y^4+4*(y+2)*y^4*x^15+(y^3+3*y^2+4*y+1)*(y+1)^2*(y+3)*y^3*x^14+4*(y^6+3*y^5+4*y^2+2*y+2)*(y+3)*y^3*x^13+3*(y^4+2*y^2+2*y+1)*(y^3+4*y^2+4*y+4)*(y+3)*y^2*x^12+2*(y+1)^2*y^2*(y^6+4*y^5+4*y^3+3*y^2+1)*x^11+3*(y^8+y^7+3*y^6+3*y^5+4*y^3+y^2+2*y+1)*y^2*x^10+4*(y^7+y^6+4*y^5+4*y^4+4*y^3+2*y+2)*(y^5+2*y^4+2*y^2+4*y+4)*y^2*x^9+(y+2)*(y^11+3*y^10+2*y^9+4*y^8+4*y^7+y^6+4*y^5+2*y^4+4*y^3+2*y^2+3*y+2)*(y^2+y+1)*y*x^8+2*(y^3+y^2+y+3)*(y+1)*(y^11+4*y^9+2*y^8+4*y^7+3*y^6+4*y^4+3*y^3+4*y^2+1)*x^7+3*(y^4+2*y^3+2*y^2+y+1)*(y^2+y+2)^2*(y^5+3*y^4+y^3+y^2+3*y+3)*y*(y+4)*x^6+4*(y^5+3*y^4+2*y^3+1)*(y^5+4*y^4+4*y^3+y^2+4)*y*(y^2+2*y+4)*(y+4)^2*x^5+(y+2)^2*(y+3)^2*y*(y+4)^5*x^4+3*(y^4+3*y^3+3*y+2)*y*(y+4)^5*x^3+3*(y+2)*y*(y+4)^5*(y^3+2*y^2+1)*x^2+2*(y^3+2*y^2+4*y+2)*y*(y+4)^6*x+(y^2+4*y+1)*y*(y+4)^7 ;
f := homogenize(f);
XF2 := Scheme(ProjectiveSpace(P),f);
FF := FunctionField(XF2);
AFF := AlgorithmicFunctionField(FF);
plc1 := Places(AFF,1);
assert Genus(AFF) eq 25 and #plc1 eq 12;
load "X1_42_split_lib.m";
// ---------------- key functions (identical to X1_42_check_orbs.m) ----------------
EncPoly := function(g)
    c := Coefficients(g); e := 0;
    for i in [#c..1 by -1] do e := 5*e + (Integers()!c[i]); end for;
    return e + 5^(#c);
end function;
ValueAt := function(g, Q)
    try
        if Valuation(g, Q) lt 0 then return false, 0; end if;
        R := ResidueClassField(Q);
        return true, R!Evaluate(g, Q);
    catch e
        return false, 0;
    end try;
end function;
MinPolyF5 := function(a)
    k := PrimeField(Parent(a));
    try
        return MinimalPolynomial(a, k);
    catch e1
        dummy := 0;
    end try;
    return MinimalPolynomial(a);
end function;
KeyFromValues := function(a, b)
    try
        polys := [ MinPolyF5(a), MinPolyF5(b), MinPolyF5(a+b), MinPolyF5(a+2*b) ];
    catch e
        return false, 0;
    end try;
    E := 5^12; key := 0;
    for g in polys do key := key*E + EncPoly(g); end for;
    return true, key;
end function;
PlaceKey := function(gu, gw, Q)
    ok1, a := ValueAt(gu, Q); if not ok1 then return false, 0; end if;
    ok2, b := ValueAt(gw, Q); if not ok2 then return false, 0; end if;
    return KeyFromValues(a, b);
end function;
KU := SeparatingElement(AFF);
KW := AFF.1;
RatKeyLine := function()
    ks := [];
    for Q in plc1 do
        ok, k := PlaceKey(KU, KW, Q);
        Append(~ks, ok select k else -1);
    end for;
    return "RAT " cat &cat[ Sprintf("%o ", k) : k in Sort(ks) ];
end function;

load "X1_42_aut_load.m";
load "X1_42_orbits_load.m";      // LoadOrbitFile

// ---- separating functions and signatures ----
cuspsum := &+ plc1;
V3, m3 := RiemannRochSpace(3*cuspsum);
CAND := [ SL_TT, SL_YY ] cat [ LocalUniformizer(Q) : Q in plc1 ] cat [ m3(V3.i) : i in [1..Dimension(V3)] ];
// write and re-read them, so that exactly the stored data is used from now on
FLINES := [ FunToLine(g) : g in CAND ];
FL := [ LineToFun(l) : l in FLINES ];
error if &or[ FL[i] ne CAND[i] : i in [1..#CAND] ], "function data does not round-trip";
SIGS := [ SignatureString(Q, FL) : Q in plc1 ];
error if #Set(SIGS) ne 12, "the stored functions do not separate the 12 rational places";
printf "Rational places separated by %o stored functions\n", #FL;
SORTED := Sort(SIGS);
LAB := [ Index(SORTED, s) : s in SIGS ];          // LAB[i] = label of plc1[i]

// ---- sigma_4 on the labels ----
t0 := Cputime();
okA, AUTS := LoadAutDataOpt("X1_42_aut_F5.txt", false);
error if not okA, "could not load X1_42_aut_F5.txt";
GEN := AUTS[4];
pw := GEN; ord := 1;
while not (pw[1] eq TT and pw[2] eq YY) do
    pw := ComposeAut(GEN, pw); ord +:= 1;
    error if ord gt 12, "generator has order > 12";
end while;
error if ord ne 12, Sprintf("sigma_4 has order %o, not 12", ord);
printf "Loaded the saved group; sigma_4 has order 12 (%o s)\n", Round(Cputime(t0));
FLS := [ ApplyAut(GEN, g) : g in FL ];
IMG := [ 0 : i in [1..12] ];
for i in [1..12] do
    s := SignatureString(plc1[i], FLS);            // signature of sigma(plc1[i])
    j := Index(SORTED, s);
    error if j eq 0, Sprintf("image of rational place %o not identified", i);
    IMG[LAB[i]] := j;
end for;
error if Set(IMG) ne {1..12}, "images do not form a permutation";
z := 1; L := 0;
repeat z := IMG[z]; L +:= 1; until z eq 1;
error if L ne 12, Sprintf("sigma_4 is not a 12-cycle on the rational places (cycle of length %o)", L);
printf "sigma_4 on the labels: %o (a 12-cycle)\n", IMG;

// ---- Q1 = zero of order 3 of y, for reference ----
pl, val := Support(Divisor(AFF.1));
q1 := 0;
for i in [1..#val] do
    if val[i] eq 3 and Degree(pl[i]) eq 1 then q1 := LAB[Index(plc1, pl[i])]; end if;
end for;

// ---- verify the saved orbit files of degrees 2..5 with the same generator ----
VERIFIED := [];
for v in [2..5] do
    t1 := Cputime();
    okf, D := LoadOrbitFile(v);
    error if not okf, Sprintf("cannot read X1_42_orbits_d%o.txt", v);
    error if D`gen ne 4, Sprintf("X1_42_orbits_d%o.txt was made with generator %o, not 4", v, D`gen);
    FILEID := AssociativeArray();                  // "bucket key" -> file index
    BK := [];
    for i in [1..D`count] do
        c := D`coeffs[i];
        if #c eq 0 then b := 0; else b := &+[ c[j]*5^(j-1) : j in [1..#c] ] + 5^#c; end if;
        s := Sprintf("%o %o", b, D`keys[i]);
        error if IsDefined(FILEID, s), Sprintf("degree %o: (bucket, key) not unique in the file", v);
        FILEID[s] := i; Append(~BK, s);
    end for;
    PL := Places(AFF, v);
    error if #PL ne D`count, Sprintf("degree %o: %o places, file has %o", v, #PL, D`count);
    seen := {}; bad := 0;
    for Q in PL do
        b, k := PlaceBucketKey(Q);
        s := Sprintf("%o %o", b, k);
        error if not IsDefined(FILEID, s), Sprintf("degree %o: a place is not in the file", v);
        i := FILEID[s];
        error if i in seen, Sprintf("degree %o: two places with the same (bucket, key)", v);
        Include(~seen, i);
        // (bucket, key) of sigma(Q), from sigma^*t, sigma^*y at Q
        okt, at := SL_ValueAt(GEN[1], Q);
        bs := 0; ks := 0;
        if okt then
            bs := SL_EncPoly(SL_MinPoly(at));
            oky, ay := SL_ValueAt(GEN[2], Q);
            if oky then okk, kk := SL_KeyFromValues(at, ay); if okk then ks := kk; end if; end if;
        end if;
        if Sprintf("%o %o", bs, ks) ne BK[D`img[i]] then bad +:= 1; end if;
    end for;
    error if bad gt 0, Sprintf("degree %o: %o stored images disagree with sigma_4", v, bad);
    printf "Degree %o: all %o stored images verified against sigma_4 (%o s)\n", v, D`count, Round(Cputime(t1));
    Append(~VERIFIED, <v, D`count>);
end for;

F := Open("X1_42_rational_labels.txt", "w");
fprintf F, "X1_42_RATLABELS v2 GEN 4\n";
fprintf F, "NF %o\n", #FLINES;
for i in [1..#FLINES] do fprintf F, "F %o %o\n", i, FLINES[i]; end for;
for k in [1..12] do fprintf F, "R %o %o\n", k, SORTED[k]; end for;
fprintf F, "S"; for k in [1..12] do fprintf F, " %o", IMG[k]; end for; fprintf F, "\n";
fprintf F, "Q1 %o\n", q1;
for vv in VERIFIED do fprintf F, "V %o %o\n", vv[1], vv[2]; end for;
fprintf F, "END\n";
delete F;
print "wrote X1_42_rational_labels.txt";
quit;

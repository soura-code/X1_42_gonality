/*
  X1_42_placeorder_test.m  --  does the ENUMERATION ORDER of the original
  degree-11 search (X1_42_degree11_sym_parallel.m, group "rest") agree between
  Magma sessions?
  ------------------------------------------------------------------
  The "rest" jobs split the divisors by a running counter, (cnt-1) mod 18.
  This covers every divisor exactly once ONLY IF every job enumerates the
  divisors in the same order.  That order is determined by
    * the order of the sequences Places(AFF,v), v = 1..7, and
    * the iteration order of the factors  Multisets(Set(Places(AFF,v)), m)
      and Multisets(Set(orbs), 1)  that the search script loops over
      (CartesianProduct of these, in lexicographic order).
  This script does exactly the same computations as a "rest" job (same
  degrees, same order) and records, for every such factor, its iteration
  order with each place replaced by an order-independent identifier
  (bucket + key, as in X1_42_orbits.m).  It writes one line per item to
  placeorder_run<r>.txt:
      SEQ v n chk           order of Places(AFF,v)
      SET v n chk           iteration order of Set(Places(AFF,v))
      SHARED v s            number of places whose identifier is not unique
      FACTOR v m n chk      iteration order of Multisets(Set(...), m)
                            (v = 0 means Set(orbs))
  If the files of several sessions are IDENTICAL (and all SHARED counts are
  0, so the identifiers really distinguish the places), the enumeration order
  is reproducible and the counter split covered every divisor.
  Usage (several separate sessions; run_placeorder.sh does this):
      magma -b run:=1 X1_42_placeorder_test.m
      magma -b run:=2 X1_42_placeorder_test.m
      diff placeorder_run1.txt placeorder_run2.txt
  ------------------------------------------------------------------
*/
if assigned run then RUN := StringToInteger(run); else RUN := 1; end if;
A<x,y> := PolynomialRing(GF(5),2);
P,homogenize := Homogenization(A);
f := x^16*y^4+4*(y+2)*y^4*x^15+(y^3+3*y^2+4*y+1)*(y+1)^2*(y+3)*y^3*x^14+4*(y^6+3*y^5+4*y^2+2*y+2)*(y+3)*y^3*x^13+3*(y^4+2*y^2+2*y+1)*(y^3+4*y^2+4*y+4)*(y+3)*y^2*x^12+2*(y+1)^2*y^2*(y^6+4*y^5+4*y^3+3*y^2+1)*x^11+3*(y^8+y^7+3*y^6+3*y^5+4*y^3+y^2+2*y+1)*y^2*x^10+4*(y^7+y^6+4*y^5+4*y^4+4*y^3+2*y+2)*(y^5+2*y^4+2*y^2+4*y+4)*y^2*x^9+(y+2)*(y^11+3*y^10+2*y^9+4*y^8+4*y^7+y^6+4*y^5+2*y^4+4*y^3+2*y^2+3*y+2)*(y^2+y+1)*y*x^8+2*(y^3+y^2+y+3)*(y+1)*(y^11+4*y^9+2*y^8+4*y^7+3*y^6+4*y^4+3*y^3+4*y^2+1)*x^7+3*(y^4+2*y^3+2*y^2+y+1)*(y^2+y+2)^2*(y^5+3*y^4+y^3+y^2+3*y+3)*y*(y+4)*x^6+4*(y^5+3*y^4+2*y^3+1)*(y^5+4*y^4+4*y^3+y^2+4)*y*(y^2+2*y+4)*(y+4)^2*x^5+(y+2)^2*(y+3)^2*y*(y+4)^5*x^4+3*(y^4+3*y^3+3*y+2)*y*(y+4)^5*x^3+3*(y+2)*y*(y+4)^5*(y^3+2*y^2+1)*x^2+2*(y^3+2*y^2+4*y+2)*y*(y+4)^6*x+(y^2+4*y+1)*y*(y+4)^7 ;
f := homogenize(f);
XF2 := Scheme(ProjectiveSpace(P),f);
FF := FunctionField(XF2);
AFF := AlgorithmicFunctionField(FF);
plc1 := Places(AFF,1);
assert Genus(AFF) eq 25 and #plc1 eq 12;

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

TT := AFF!(BaseField(AFF).1); YY := AFF.1;          // TT = t, YY = y, ApplyAut, ComposeAut, LoadAutDataOpt

// data identifying a place: coefficients of the polynomial of F5(t) below it
// ([] for t = infinity), bucket code, key
PlaceData := function(Q)
    okt, at := ValueAt(TT, Q);
    if okt then
        mp := MinPolyF5(at);
        coeffs := [ Integers()!c : c in Coefficients(mp) ];
        b := EncPoly(mp);
    else
        coeffs := []; b := 0;
    end if;
    k := 0;
    if okt then
        oky, ay := ValueAt(YY, Q);
        if oky then ok, kk := KeyFromValues(at, ay); if ok then k := kk; end if; end if;
    end if;
    return coeffs, b, k;
end function;


// ---- the same set-up as a "rest" job of X1_42_degree11_sym_parallel.m ----
pl,val := Support(Divisor(AFF.1));
for i in [1..#val] do
    if val[i] eq  3 then orb1 := pl[i]; end if;
    if val[i] eq -7 then orb2 := pl[i]; end if;
end for;
orbs := [orb1, orb2];
B := 9;
Parts := [ l : l in Partitions(B) | l ne [B] and l ne [B-1, 1] ];
NeedDeg := Sort(Setseq({1} join &join[ Set(l) : l in Parts ]));
PlacesOfDeg := AssociativeArray();
for v in NeedDeg do
    t1 := Cputime();
    PlacesOfDeg[v] := (v eq 1) select plc1 else Places(AFF, v);
    printf "Places of degree %o: %o  (%o s)\n", v, #PlacesOfDeg[v], Round(Cputime(t1));
end for;

MODP := 2305843009213693951;
ID := function(p)
    c, b, k := PlaceData(p);
    return (b mod 1000000007) + 7*(k mod 1000000007);
end function;
Fold := function(chk, h) return (chk*1000003 + h) mod MODP; end function;

OUT := Sprintf("placeorder_run%o.txt", RUN);
F := Open(OUT, "w");
for v in NeedDeg do
    PL := PlacesOfDeg[v];
    ids := [ ID(p) : p in PL ];
    chk := 0; for h in ids do chk := Fold(chk, h); end for;
    fprintf F, "SEQ %o %o %o\n", v, #PL, chk;
    Mids := SequenceToMultiset(ids);
    sh := #[ h : h in ids | Multiplicity(Mids, h) gt 1 ];
    fprintf F, "SHARED %o %o\n", v, sh;
    chk := 0; for p in Set(PL) do chk := Fold(chk, ID(p)); end for;
    fprintf F, "SET %o %o %o\n", v, #PL, chk;
    printf "degree %o: sequence/set order recorded (%o places with non-unique identifier)\n", v, sh;
end for;

// the factors the search loops over
FAC := {};
for lam in Parts do
    vals := Sort(Setseq(Set(lam)));
    for v in vals do
        m := #[ z : z in lam | z eq v ];
        if v eq 1 then
            Include(~FAC, 0*100 + 1);
            if m ge 2 then Include(~FAC, 1*100 + (m-1)); end if;
        else
            Include(~FAC, v*100 + m);
        end if;
    end for;
end for;
for fm in Sort(Setseq(FAC)) do
    v := fm div 100; m := fm mod 100;
    S  := (v eq 0) select Set(orbs) else Set(PlacesOfDeg[v]);
    SS := Setseq(S);
    idsS := (m eq 1) select [] else [ ID(p) : p in SS ];
    chk := 0; n := 0;
    for ms in Multisets(S, m) do
        n +:= 1;
        if m eq 1 then
            for p in ms do h := ID(p); end for;
        else
            L := Sort([ idsS[Index(SS, p)] : p in ms ]);
            h := 0; for z in L do h := Fold(h, z); end for;
        end if;
        chk := Fold(chk, h);
    end for;
    fprintf F, "FACTOR %o %o %o %o\n", v, m, n, chk;
    printf "factor Multisets(Set(%o), %o): %o elements recorded\n",
           v eq 0 select "orbs" else Sprintf("Places(AFF,%o)", v), m, n;
end for;
delete F;
printf "wrote %o\n", OUT;
quit;

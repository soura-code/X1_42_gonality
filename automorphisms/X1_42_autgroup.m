/*
  X1_42_autgroup.m  --  compute Aut_{F5}(X1(42)) once and export it
  ------------------------------------------------------------------
  Builds X1(42)/F5 from the FFFc42 model (same construction as the search
  scripts), computes Automorphisms(AFF), verifies the result, and writes:

    X1_42_aut_F5.txt   plain-text, reloadable by LoadAutData() in
                       X1_42_aut_load.m (used by the search scripts)
    X1_42_aut_F5.json  the same data for a database, plus: group order,
                       group structure, element orders, multiplication
                       table, and the permutation of the 12 rational places

  Each automorphism s is stored as the images s(t), s(y) of the generators of
  AFF = F5(t)[y]/(defining polynomial); every element is a list of NDEG pairs
  [numerator coefficients, denominator coefficients] (constant term first),
  meaning   sum_j  num_j(t)/den_j(t) * y^(j-1).

  Usage:  magma -b X1_42_autgroup.m      (X1_42_aut_load.m in the same folder)
  ------------------------------------------------------------------
*/
SetColumns(0);
A<x,y> := PolynomialRing(GF(5),2);
P,homogenize := Homogenization(A);
print("f := The polynomial from math.mit.edu/~drew/FFFc42.txt reduced mod 5");
f := x^16*y^4+4*(y+2)*y^4*x^15+(y^3+3*y^2+4*y+1)*(y+1)^2*(y+3)*y^3*x^14+4*(y^6+3*y^5+4*y^2+2*y+2)*(y+3)*y^3*x^13+3*(y^4+2*y^2+2*y+1)*(y^3+4*y^2+4*y+4)*(y+3)*y^2*x^12+2*(y+1)^2*y^2*(y^6+4*y^5+4*y^3+3*y^2+1)*x^11+3*(y^8+y^7+3*y^6+3*y^5+4*y^3+y^2+2*y+1)*y^2*x^10+4*(y^7+y^6+4*y^5+4*y^4+4*y^3+2*y+2)*(y^5+2*y^4+2*y^2+4*y+4)*y^2*x^9+(y+2)*(y^11+3*y^10+2*y^9+4*y^8+4*y^7+y^6+4*y^5+2*y^4+4*y^3+2*y^2+3*y+2)*(y^2+y+1)*y*x^8+2*(y^3+y^2+y+3)*(y+1)*(y^11+4*y^9+2*y^8+4*y^7+3*y^6+4*y^4+3*y^3+4*y^2+1)*x^7+3*(y^4+2*y^3+2*y^2+y+1)*(y^2+y+2)^2*(y^5+3*y^4+y^3+y^2+3*y+3)*y*(y+4)*x^6+4*(y^5+3*y^4+2*y^3+1)*(y^5+4*y^4+4*y^3+y^2+4)*y*(y^2+2*y+4)*(y+4)^2*x^5+(y+2)^2*(y+3)^2*y*(y+4)^5*x^4+3*(y^4+3*y^3+3*y+2)*y*(y+4)^5*x^3+3*(y+2)*y*(y+4)^5*(y^3+2*y^2+1)*x^2+2*(y^3+2*y^2+4*y+2)*y*(y+4)^6*x+(y^2+4*y+1)*y*(y+4)^7 ;
f := homogenize(f);
XF2 := Scheme(ProjectiveSpace(P),f);
FF := FunctionField(XF2);
AFF := AlgorithmicFunctionField(FF);
plc1 := Places(AFF,1);
cuspsum := &+ plc1;

assert Genus(AFF) eq 25;

// ----- canonical keys of places (used for diamond/automorphism orbit representatives) -----
// The key of a place Q is computed from the residues a = KU(Q), b = KW(Q) of two fixed
// functions KU, KW: it encodes the minimal polynomials over F5 of a, b, a+b, a+2b.
// It depends only on the place (not on the choice of conjugate), so it is canonical.
// Different places may share a key; that only makes the reduction keep MORE places.
EncPoly := function(g)
    c := Coefficients(g); e := 0;
    for i in [#c..1 by -1] do e := 5*e + (Integers()!c[i]); end for;
    return e + 5^(#c);                         // 5^(deg+1) marks the degree
end function;

ValueAt := function(g, Q)                      // residue of g at Q, or false at a pole
    try
        if Valuation(g, Q) lt 0 then return false, 0; end if;
        R := ResidueClassField(Q);
        return true, R!Evaluate(g, Q);
    catch e
        return false, 0;
    end try;
end function;

MinPolyF5 := function(a)                       // minimal polynomial over F5
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

// key of the point with coordinates (gu(Q), gw(Q)); with gu = s(KU), gw = s(KW) for an
// automorphism s this is the key of the image of Q under the morphism induced by s.
PlaceKey := function(gu, gw, Q)
    ok1, a := ValueAt(gu, Q); if not ok1 then return false, 0; end if;
    ok2, b := ValueAt(gw, Q); if not ok2 then return false, 0; end if;
    return KeyFromValues(a, b);
end function;

KU := SeparatingElement(AFF);
KW := AFF.1;
if KU eq KW then
    print "WARNING: SeparatingElement(AFF) = AFF.1; keys use one function only (still safe, less effective).";
end if;

// fingerprint of KU, KW: keys of the rational places (must agree between stage 1 and stage 2)
RatKeyLine := function()
    ks := [];
    for Q in plc1 do
        ok, k := PlaceKey(KU, KW, Q);
        Append(~ks, ok select k else -1);
    end for;
    return "RAT " cat &cat[ Sprintf("%o ", k) : k in Sort(ks) ];
end function;

load "X1_42_aut_load.m";

KeyOr := function(gu, gw, Q)
    ok, k := PlaceKey(gu, gw, Q); return ok select k else -1;
end function;

// diamond-orbit representatives used by the search scripts (zero / pole of x)
pl,val := Support(Divisor(AFF.1));
for i in [1..#val] do
    if val[i] eq  3 then orb1 := pl[i]; end if;
    if val[i] eq -7 then orb2 := pl[i]; end if;
end for;

printf "AFF = F5(t)[y]/(f),  [AFF : F5(t)] = %o\n", NDEG;

// ---------- compute the automorphism group ----------
t0 := Cputime();
okprev, A := LoadAutData("X1_42_aut_F5.txt");      // saved by X1_42_check_orbs.m ?
if okprev then
    G := [];      // no Magma maps needed; the saved data were verified when loaded
    printf "Loaded %o automorphisms from X1_42_aut_F5.txt (%o s)\n", #A, Round(Cputime(t0));
else
    G := Automorphisms(AFF);
    printf "Automorphisms(AFF): %o maps (%o s)\n", #G, Round(Cputime(t0));
    A := [];
    for s in G do
        a := AutFromMap(s);
        if AutIndex(A, a) eq 0 then Append(~A, a); end if;
    end for;
end if;
if AutIndex(A, <TT, YY>) eq 0 then Insert(~A, 1, <TT, YY>); end if;
// put the identity first
idx := AutIndex(A, <TT, YY>);
if idx ne 1 then a := A[idx]; Remove(~A, idx); Insert(~A, 1, a); end if;
n := #A;
printf "%o distinct automorphisms (identity first)\n", n;

// ---------- verification ----------
tests := [TT, YY, TT + YY, YY^2/(TT + 1), KU, KW];
for i in [1..n] do
    error if not IsHomData(A[i]), Sprintf("automorphism %o does not respect the equation", i);
end for;
// the stored data reproduce Magma's maps on test functions
for s in G do
    a := AutFromMap(s);
    for f in tests do
        error if ApplyAut(a, f) ne s(f), "ApplyAut disagrees with the Magma map";
    end for;
end for;
// round trip through the integer encoding
for a in A do
    error if DataToElt(EltToData(a[1])) ne a[1] or DataToElt(EltToData(a[2])) ne a[2],
        "integer encoding does not round-trip";
end for;
print "Verification passed: equation, agreement with Magma maps, encoding round-trip.";

// ---------- group structure ----------
T := [ [ AutIndex(A, ComposeAut(A[i], A[j])) : j in [1..n] ] : i in [1..n] ];
closed := &and[ T[i][j] ne 0 : i, j in [1..n] ];
error if not closed, "Automorphisms(AFF) is not closed under composition";
perms := [ Sym(n) ! [ T[i][j] : j in [1..n] ] : i in [1..n] ];   // left multiplication
GP := sub< Sym(n) | perms >;
error if #GP ne n, "regular representation has the wrong order";
orders := [ Order(p) : p in perms ];
struct := "";
try struct := GroupName(GP); catch e struct := ""; end try;
idg := "";
try idg := Sprint(IdentifyGroup(GP)); catch e idg := ""; end try;
printf "Group order %o, structure %o, IdentifyGroup %o, abelian %o\n", n, struct, idg, IsAbelian(GP);
printf "Element orders: %o\n", orders;

// ---------- action on the rational places ----------
ratkey := [ KeyOr(KU, KW, Q) : Q in plc1 ];
ratok := (-1 notin ratkey) and (#Set(ratkey) eq #plc1);
ratperm := [];
for a in A do
    if ratok then
        img := [ Index(ratkey, KeyOr(ApplyAut(a, KU), ApplyAut(a, KW), Q)) : Q in plc1 ];
    else
        img := [ 0 : Q in plc1 ];
    end if;
    Append(~ratperm, img);
end for;
if ratok then
    rorbits := { { ratperm[i][q] : i in [1..n] } : q in [1..#plc1] };
    printf "Orbits on the %o rational places: %o\n", #plc1, rorbits;
    io1 := Index(plc1, orb1); io2 := Index(plc1, orb2);
    printf "orb1 = rational place %o, orb2 = rational place %o; they meet every orbit: %o\n",
           io1, io2, &and[ io1 in O or io2 in O : O in rorbits ];
else
    print "Rational places do not have distinct keys; permutation not computed.";
end if;

// ---------- write the reloadable text file ----------
WriteAutText("X1_42_aut_F5.txt", A);
ok, B := LoadAutData("X1_42_aut_F5.txt");
error if not ok or #B ne n or &or[ AutIndex(A, b) eq 0 : b in B ],
    "re-reading X1_42_aut_F5.txt failed";
print "Wrote X1_42_aut_F5.txt and re-read it successfully.";

// ---------- write the JSON file ----------
JStr := function(x)        // JSON string, newlines removed, quotes/backslashes escaped
    r := Sprint(x); out := "";
    for i in [1..#r] do
        c := r[i];
        if c eq "\n" then out cat:= " ";
        elif c eq "\"" then out cat:= "\\\"";
        elif c eq "\\" then out cat:= "\\\\";
        else out cat:= c; end if;
    end for;
    return "\"" cat out cat "\"";
end function;
JInts := function(L)
    return "[" cat (#L eq 0 select "" else &cat [ Sprintf("%o", L[i]) cat (i lt #L select "," else "") : i in [1..#L] ]) cat "]";
end function;
JElt := function(g)        // list of [num, den] pairs
    d := EltToData(g);
    return "[" cat &cat [ "[" cat JInts(d[j][1]) cat "," cat JInts(d[j][2]) cat "]" cat (j lt #d select "," else "") : j in [1..#d] ] cat "]";
end function;

J := Open("X1_42_aut_F5.json", "w");
fprintf J, "{\n";
fprintf J, "  \"curve\": \"X1(42)\",\n  \"field\": \"F5\",\n";
fprintf J, "  \"model\": \"FFFc42 (math.mit.edu/~drew/FFFc42.txt) reduced mod 5, Magma AlgorithmicFunctionField\",\n";
fprintf J, "  \"genus\": %o,\n", Genus(AFF);
fprintf J, "  \"function_field\": {\n";
fprintf J, "    \"description\": \"AFF = F5(t)[y]/(f(t,y)); t = AFF!BaseField(AFF).1, y = AFF.1\",\n";
fprintf J, "    \"degree_over_F5_t\": %o,\n", NDEG;
fprintf J, "    \"defining_polynomial\": %o,\n", JStr(DefiningPolynomial(AFF));
fprintf J, "    \"element_encoding\": \"list of %o pairs [num, den]; coefficient lists over F5, constant term first; element = sum_j num_j(t)/den_j(t) * y^(j-1)\"\n", NDEG;
fprintf J, "  },\n";
fprintf J, "  \"rational_place_keys\": %o,\n", JInts(ratkey);
fprintf J, "  \"group\": {\n";
fprintf J, "    \"order\": %o,\n", n;
fprintf J, "    \"structure\": %o,\n", JStr(struct);
fprintf J, "    \"identify_group\": %o,\n", JStr(idg);
fprintf J, "    \"abelian\": %o,\n", IsAbelian(GP) select "true" else "false";
fprintf J, "    \"composition_convention\": \"table[i][j] = k means A_i o A_j = A_k (apply A_j first); indices 1-based\",\n";
fprintf J, "    \"multiplication_table\": [%o]\n",
        &cat [ JInts(T[i]) cat (i lt n select "," else "") : i in [1..n] ];
fprintf J, "  },\n";
fprintf J, "  \"automorphisms\": [\n";
for i in [1..n] do
    fprintf J, "    {\"id\": %o, \"order\": %o, \"is_identity\": %o,\n", i, orders[i], i eq 1 select "true" else "false";
    fprintf J, "     \"perm_rational_places\": %o,\n", JInts(ratperm[i]);
    fprintf J, "     \"image_t\": %o,\n", JElt(A[i][1]);
    fprintf J, "     \"image_y\": %o}%o\n", JElt(A[i][2]), i lt n select "," else "";
end for;
fprintf J, "  ]\n}\n";
Flush(J); delete J;
print "Wrote X1_42_aut_F5.json";
quit;

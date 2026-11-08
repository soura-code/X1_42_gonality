/*
  X1_42_check_orbs.m
  ------------------------------------------------------------------
  Verifies the one hypothesis of the X1(42)/F5 search that the search code
  takes for granted:

    (*)  every orbit of Aut_{F5}(X) on the 12 rational places X(F5) contains
         orb1 (the zero of order 3 of x) or orb2 (the pole of order 7 of x).

  (*) is exactly what is needed for the rational-anchor reduction: any divisor E
  containing a rational place P can be moved by an automorphism so that it
  contains orb1 or orb2 (Lemma: automorphisms preserve degrees of functions in
  L(S+E)).  The search scripts never check (*); they read orb1, orb2 off div(x).

  Method (independent of the "key" machinery of the orbit-reduced search):
    for each automorphism s (field automorphism of AFF over F5) and each rational
    place P, the image place Q = phi(P) (where s(g) = g o phi) is characterised by
        v_Q(g) = v_P(s(g))   for all functions g.
    We test this for a list of functions g and require a UNIQUE rational Q.
  Output: group order, the permutation of the 12 rational places induced by each
  automorphism, the orbits, and PASSED / FAILED / INCONCLUSIVE for (*).

  It also SAVES the automorphism group, so it never has to be recomputed:
    X1_42_aut_F5.txt            reloadable with LoadAutData() from X1_42_aut_load.m
                                (used automatically by X1_42_orbitkeys.m and
                                X1_42_autgroup.m, which adds the JSON export)
    X1_42_aut_F5_readable.txt   every element written out: s(t), s(y), its order
                                and its permutation of the rational places

  Usage:  magma -b X1_42_check_orbs.m      (X1_42_aut_load.m in the same folder)
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
n := #plc1;
printf "Genus %o, %o rational places\n", Genus(AFF), n;

pl,val := Support(Divisor(AFF.1));
for i in [1..#val] do
    if val[i] eq  3 then orb1 := pl[i]; end if;
    if val[i] eq -7 then orb2 := pl[i]; end if;
end for;
assert Degree(orb1) eq 1 and Degree(orb2) eq 1;
i1 := Index(plc1, orb1); i2 := Index(plc1, orb2);
printf "orb1 = rational place %o, orb2 = rational place %o\n", i1, i2;

// test functions: x, a separating element, and their translates by constants
// key functions + save/load helpers
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

// ---------- the automorphism group: reuse the saved file if possible ----------
t0 := Cputime();
okprev, AUTS := LoadAutData("X1_42_aut_F5.txt");
if okprev then
    printf "Loaded %o automorphisms from X1_42_aut_F5.txt (verified: equation, identity, closure) (%o s)\n",
           #AUTS, Round(Cputime(t0));
    WriteAutText("X1_42_aut_F5.txt", AUTS);          // rewrite in the robust v2 format
else
    G := Automorphisms(AFF);
    printf "Automorphisms(AFF): %o automorphisms (%o s)\n", #G, Round(Cputime(t0));
    AUTS := [];
    for s in G do
        a := AutFromMap(s);
        if AutIndex(AUTS, a) eq 0 then Append(~AUTS, a); end if;
    end for;
    for s in G do
        a := AutFromMap(s);
        for g in [TT, YY, TT + YY, KU] do
            error if ApplyAut(a, g) ne s(g), "stored images disagree with Magma's map";
        end for;
    end for;
    WriteAutText("X1_42_aut_F5.txt", AUTS);
end if;
idx := AutIndex(AUTS, <TT, YY>);
if idx eq 0 then Insert(~AUTS, 1, <TT, YY>);
elif idx ne 1 then a := AUTS[idx]; Remove(~AUTS, idx); Insert(~AUTS, 1, a); end if;
okre, AUTS2 := LoadAutData("X1_42_aut_F5.txt");
error if not okre or #AUTS2 ne #AUTS, "re-reading X1_42_aut_F5.txt failed";
printf "X1_42_aut_F5.txt holds %o automorphisms (re-read and verified).\n", #AUTS;

// ---------- identifying rational places ----------
// For a rational place P and a function g put  e(g,P) = [ v_P(g), g(P) ]  (residue only if v_P(g) = 0).
// If sigma = phi^* is an automorphism, then  v_{phi(P)}(g) = v_P(sigma(g))  and  g(phi(P)) = sigma(g)(P),
// so the signature of phi(P) can be computed from sigma(g) at P.
// Test functions: x, a separating element, and a basis of L(3*S) (S = sum of rational places),
// which separates the rational places.
xx := AFF.1; ss := SeparatingElement(AFF);
V3, m3 := RiemannRochSpace(3*cuspsum);
TESTF := [xx, ss] cat [ m3(V3.i) : i in [1..Dimension(V3)] ];
printf "%o test functions\n", #TESTF;

Entry := function(g, Q)
    e := Valuation(g, Q);
    if e ne 0 then return [e, -1]; end if;
    return [0, Integers()!(GF(5)!Evaluate(g, Q))];
end function;
SIGS := [ [ Entry(g, Q) : g in TESTF ] : Q in plc1 ];
sigok := #Set(SIGS) eq n;
if not sigok then print "WARNING: rational places not separated by the test functions."; end if;

PermOf := function(a)
    imgs := [ ApplyAut(a, g) : g in TESTF ];
    img := [];
    for P in plc1 do
        sig := [ Entry(h, P) : h in imgs ];
        cand := [ j : j in [1..n] | SIGS[j] eq sig ];
        Append(~img, #cand eq 1 select cand[1] else 0);
    end for;
    return img;
end function;

PERMS := [ PermOf(a) : a in AUTS ];
inconclusive := (not sigok) or &or[ Set(p) ne {1..n} : p in PERMS ];

if inconclusive then
    print "RESULT: INCONCLUSIVE (images of rational places could not be identified uniquely).";
    printf "Permutations found (0 = not identified): %o\n", PERMS;
else
    Sn := Sym(n);
    H := sub< Sn | [ Sn!p : p in PERMS ] >;
    orbs := Orbits(H);
    printf "Image of Aut in Sym(%o): order %o\n", n, #H;
    printf "Permutations of the rational places:\n";
    for p in PERMS do printf "  %o\n", p; end for;
    printf "Orbits on rational places: %o\n", orbs;
    ok := &and[ i1 in O or i2 in O : O in orbs ];
    if ok then
        printf "RESULT: PASSED -- every orbit contains orb1 or orb2 (%o orbits).\n", #orbs;
    else
        print "RESULT: FAILED -- some orbit contains neither orb1 nor orb2. The rational-anchor reduction is NOT valid.";
    end if;
end if;

// ---------- print and save every element ----------
nA := #AUTS;
TAB := [ [ AutIndex(AUTS, ComposeAut(AUTS[i], AUTS[j])) : j in [1..nA] ] : i in [1..nA] ];
closedA := &and[ TAB[i][j] ne 0 : i, j in [1..nA] ];
ORD := [];
for i in [1..nA] do
    k := 1; c := AUTS[i];
    while AutIndex(AUTS, c) ne 1 and k le nA do c := ComposeAut(AUTS[i], c); k +:= 1; end while;
    Append(~ORD, k);
end for;

R := Open("X1_42_aut_F5_readable.txt", "w");
fprintf R, "Aut_F5(X1(42)), FFFc42 model mod 5.  AFF = F5(t)[y]/(f),  t = AFF!BaseField(AFF).1,  y = AFF.1\n";
fprintf R, "Defining polynomial f: %o\n", DefiningPolynomial(AFF);
fprintf R, "Group order: %o   (closed under composition: %o)\n", nA, closedA;
fprintf R, "Composition: table[i][j] = k means A_i o A_j = A_k (apply A_j first)\n%o\n\n", TAB;
print "";
printf "==== Elements of Aut_F5(X1(42)) (%o) ====\n", nA;
for i in [1..nA] do
    printf "A_%o  (order %o)  permutation of rational places: %o\n", i, ORD[i], PERMS[i];
    printf "   t |-> %o\n   y |-> %o\n", AUTS[i][1], AUTS[i][2];
    fprintf R, "A_%o  (order %o)\n  permutation of rational places: %o\n", i, ORD[i], PERMS[i];
    fprintf R, "  t |-> %o\n  y |-> %o\n\n", AUTS[i][1], AUTS[i][2];
end for;
Flush(R); delete R;
print "Wrote X1_42_aut_F5_readable.txt";
quit;

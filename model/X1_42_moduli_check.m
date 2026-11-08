/*
  X1_42_moduli_check.m
  ------------------------------------------------------------------
  Proves that the function field F of the FFFc42 model reduced mod 5
  (the field used in all searches) is the function field of X_1(42)/F_5,
  using the moduli interpretation.  No reduction from characteristic 0 and
  no isomorphism search is needed.

  Sutherland (math.mit.edu/~drew/X1_altcurves.html): the map from f(x,y)=0
  to the raw Tate-normal-form curve F_42(r,s)=0 is
        r = (x^2 y - x y + y - 1) / (x^2 y - x),     s = (x y - y + 1) / (x y),
  and the Tate normal form is  E(b,c): y^2 + (1-c)xy - by = x^3 - bx^2,
  P = (0,0), with  b = r s (r-1),  c = s (r-1).
  All coefficients are integers, so the formulas make sense over F_5.

  Checks
   A. F has exact constant field F_5 (f mod 5 absolutely irreducible) and
      genus 25.
   B. r, s are defined in F; E is an elliptic curve over F; P has exact
      order 42.  By the moduli property of Y_1(42) over F_5 this gives a
      morphism  phi : C -> X_1(42)_{F_5}  with  j o phi = j(E).
   C. deg j(E) = 576 = [PSL_2(Z) : +-Gamma_1(42)] = deg(j : X_1(42) -> P^1),
      so deg(phi) = 1 and F = F_5(X_1(42)).
   D. (optional) if F42.m is present, check F_42(r,s) = 0 in F.
         prepare:  echo "F42 := $(cat F42.txt);" > F42.m

  Run:  magma -b X1_42_moduli_check.m
  ------------------------------------------------------------------
*/
t0 := Cputime();
k := GF(5);

// ---------------- the FFFc42 model mod 5, exactly as in the search scripts ----------
A<x,y> := PolynomialRing(k, 2);
f := x^16*y^4+4*(y+2)*y^4*x^15+(y^3+3*y^2+4*y+1)*(y+1)^2*(y+3)*y^3*x^14+4*(y^6+3*y^5+4*y^2+2*y+2)*(y+3)*y^3*x^13+3*(y^4+2*y^2+2*y+1)*(y^3+4*y^2+4*y+4)*(y+3)*y^2*x^12+2*(y+1)^2*y^2*(y^6+4*y^5+4*y^3+3*y^2+1)*x^11+3*(y^8+y^7+3*y^6+3*y^5+4*y^3+y^2+2*y+1)*y^2*x^10+4*(y^7+y^6+4*y^5+4*y^4+4*y^3+2*y+2)*(y^5+2*y^4+2*y^2+4*y+4)*y^2*x^9+(y+2)*(y^11+3*y^10+2*y^9+4*y^8+4*y^7+y^6+4*y^5+2*y^4+4*y^3+2*y^2+3*y+2)*(y^2+y+1)*y*x^8+2*(y^3+y^2+y+3)*(y+1)*(y^11+4*y^9+2*y^8+4*y^7+3*y^6+4*y^4+3*y^3+4*y^2+1)*x^7+3*(y^4+2*y^3+2*y^2+y+1)*(y^2+y+2)^2*(y^5+3*y^4+y^3+y^2+3*y+3)*y*(y+4)*x^6+4*(y^5+3*y^4+2*y^3+1)*(y^5+4*y^4+4*y^3+y^2+4)*y*(y^2+2*y+4)*(y+4)^2*x^5+(y+2)^2*(y+3)^2*y*(y+4)^5*x^4+3*(y^4+3*y^3+3*y+2)*y*(y+4)^5*x^3+3*(y+2)*y*(y+4)^5*(y^3+2*y^2+1)*x^2+2*(y^3+2*y^2+4*y+2)*y*(y+4)^6*x+(y^2+4*y+1)*y*(y+4)^7;

fac := Factorization(f);
assert #fac eq 1 and fac[1][2] eq 1;              // irreducible over F_5

// F = F_5(x)[y]/(f), with X = x and Y = y
Rx<u> := FunctionField(k);
Py<V> := PolynomialRing(Rx);
F := FunctionField(Evaluate(f, [u, V]));
X := F ! u;  Y := F.1;
assert Evaluate(f, [X, Y]) eq 0;

// ---------------- A ----------------------------------------------------------------
assert Degree(ExactConstantField(F)) eq 1;        // f is absolutely irreducible
assert Genus(F) eq 25;
printf "A: exact constant field F_5, genus 25  (%o s)\n", Cputime(t0);

// ---------------- B ----------------------------------------------------------------
dr := X^2*Y - X;  ds := X*Y;
assert dr ne 0 and ds ne 0;                       // the formulas are defined on F
r := (X^2*Y - X*Y + Y - 1) / dr;
s := (X*Y - Y + 1) / ds;
b := r*s*(r - 1);
c := s*(r - 1);
E := EllipticCurve([1 - c, -b, -b, 0, 0]);        // y^2 + (1-c)xy - by = x^3 - bx^2
assert Discriminant(E) ne 0;
P := E ! [0, 0];
P21 := 21*P;
assert not IsZero(P21);                           // order does not divide 21
assert IsZero(2*P21);                             // [42]P = O
assert not IsZero(14*P);                          // order does not divide 14
assert not IsZero(6*P);                           // order does not divide 6
printf "B: E/F is an elliptic curve and (0,0) has exact order 42  (%o s)\n", Cputime(t0);

// ---------------- B': Sutherland's inverse map recovers x and y --------------------
//   x = (s - r)/(rs - 2r + 1),   y = (rs - 2r + 1)/(s^2 - s - r + 1)
// so F = F_5(x,y) = F_5(r,s) = F_5(b,c)  (hypothesis (b) of Lemma 2.4, verified in Lemma 2.5)
d1 := r*s - 2*r + 1;  d2 := s^2 - s - r + 1;
assert d1 ne 0 and d2 ne 0;
assert (s - r)/d1 eq X;
assert d1/d2 eq Y;
assert c ne 0 and b ne c;
assert b/c eq r and c^2/(b - c) eq s;
printf "B': inverse map recovers x, y; F = F_5(b,c)  (%o s)\n", Cputime(t0);

// ---------------- C (extra check, not needed for the proof) ------------------------
jE := jInvariant(E);
assert not IsConstant(jE);
dj := Degree(PoleDivisor(Divisor(jE)));
printf "C: deg j(E) = %o  (%o s)\n", dj, Cputime(t0);
assert dj eq 576;

// ---------------- D (optional consistency check with the raw equation) -------------
// F42.m contains "F42 := <polynomial in r, s>;".  Since r and s are now the
// elements of F defined above, loading it evaluates F_42(r,s) directly in F.
// Uncomment the next two lines if F42.m is present (not needed for the proof).
// load "F42.m";
// printf "D: F_42(r,s) = 0 in F: %o\n", F42 eq 0;

print "ALL CHECKS PASSED: F is the function field of X_1(42) over F_5.";
quit;

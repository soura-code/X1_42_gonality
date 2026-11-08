/*
  X1_42_irreducibility_mod5.m
  Checks that the FFFc42 plane model of X_1(42), reduced mod 5, is
  irreducible over F_5 and absolutely irreducible (irreducible over F_5-bar),
  using three independent tests.   Run:  magma -b X1_42_irreducibility_mod5.m
*/
A<x,y> := PolynomialRing(GF(5),2);
f := x^16*y^4+4*(y+2)*y^4*x^15+(y^3+3*y^2+4*y+1)*(y+1)^2*(y+3)*y^3*x^14+4*(y^6+3*y^5+4*y^2+2*y+2)*(y+3)*y^3*x^13+3*(y^4+2*y^2+2*y+1)*(y^3+4*y^2+4*y+4)*(y+3)*y^2*x^12+2*(y+1)^2*y^2*(y^6+4*y^5+4*y^3+3*y^2+1)*x^11+3*(y^8+y^7+3*y^6+3*y^5+4*y^3+y^2+2*y+1)*y^2*x^10+4*(y^7+y^6+4*y^5+4*y^4+4*y^3+2*y+2)*(y^5+2*y^4+2*y^2+4*y+4)*y^2*x^9+(y+2)*(y^11+3*y^10+2*y^9+4*y^8+4*y^7+y^6+4*y^5+2*y^4+4*y^3+2*y^2+3*y+2)*(y^2+y+1)*y*x^8+2*(y^3+y^2+y+3)*(y+1)*(y^11+4*y^9+2*y^8+4*y^7+3*y^6+4*y^4+3*y^3+4*y^2+1)*x^7+3*(y^4+2*y^3+2*y^2+y+1)*(y^2+y+2)^2*(y^5+3*y^4+y^3+y^2+3*y+3)*y*(y+4)*x^6+4*(y^5+3*y^4+2*y^3+1)*(y^5+4*y^4+4*y^3+y^2+4)*y*(y^2+2*y+4)*(y+4)^2*x^5+(y+2)^2*(y+3)^2*y*(y+4)^5*x^4+3*(y^4+3*y^3+3*y+2)*y*(y+4)^5*x^3+3*(y+2)*y*(y+4)^5*(y^3+2*y^2+1)*x^2+2*(y^3+2*y^2+4*y+2)*y*(y+4)^6*x+(y^2+4*y+1)*y*(y+4)^7 ;

// (1) the polynomial is irreducible over F_5
fac := Factorization(f);
printf "(1) factors of f over F_5: %o  (degrees in x,y: %o)\n", #fac,
       [ <Degree(g[1],x), Degree(g[1],y), g[2]> : g in fac ];
assert #fac eq 1 and fac[1][2] eq 1;

// (2) the projective plane curve is absolutely irreducible
P2<X,Y,Z> := ProjectiveSpace(GF(5),2);
C := Curve(P2, Homogenization(Evaluate(f,[X,Y]), Z));
ok := IsAbsolutelyIrreducible(C);
printf "(2) IsAbsolutelyIrreducible(C): %o\n", ok;
assert ok;

// (3) the exact constant field of the function field is F_5
//     (an irreducible curve over F_5 is absolutely irreducible iff F_5 is
//      algebraically closed in its function field)
AFF := AlgorithmicFunctionField(FunctionField(C));
k := ExactConstantField(AFF);
printf "(3) exact constant field: %o (degree %o over F_5)\n", k, Degree(k, GF(5));
assert Degree(k, GF(5)) eq 1;

// (4) extra check: f stays irreducible over F_{5^m}, m = 2..6
for m in [2..6] do
    Am := PolynomialRing(GF(5^m),2);
    n := #Factorization(Am!f);
    printf "(4) factors over F_5^%o: %o\n", m, n;
    assert n eq 1;
end for;

printf "genus of the function field: %o\n", Genus(AFF);
print "ALL CHECKS PASSED: the model is absolutely irreducible mod 5.";
quit;

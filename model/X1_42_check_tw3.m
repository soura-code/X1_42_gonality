/*
  X1_42_check_tw3.m  --  confirm  t = 1/y,  w = x/y,  Phi(t,w) = t^23 * f(w/t, 1/t),
  and write Phi to X1_42_Phi.txt with named variables t, w.
  (Magma homogenizes as [H : x : y] and uses the affine chart y <> 0.)
  Run:  magma -b X1_42_check_tw3.m
  The output file X1_42_Phi.txt is itself a Magma file:  load "X1_42_Phi.txt";
  defines R<t,w> and Phi.
*/
SetColumns(0);   // no line wrapping in printed output (also in PrintFile)
A<x,y> := PolynomialRing(GF(5),2);
P,homogenize := Homogenization(A);
f0 := x^16*y^4+4*(y+2)*y^4*x^15+(y^3+3*y^2+4*y+1)*(y+1)^2*(y+3)*y^3*x^14+4*(y^6+3*y^5+4*y^2+2*y+2)*(y+3)*y^3*x^13+3*(y^4+2*y^2+2*y+1)*(y^3+4*y^2+4*y+4)*(y+3)*y^2*x^12+2*(y+1)^2*y^2*(y^6+4*y^5+4*y^3+3*y^2+1)*x^11+3*(y^8+y^7+3*y^6+3*y^5+4*y^3+y^2+2*y+1)*y^2*x^10+4*(y^7+y^6+4*y^5+4*y^4+4*y^3+2*y+2)*(y^5+2*y^4+2*y^2+4*y+4)*y^2*x^9+(y+2)*(y^11+3*y^10+2*y^9+4*y^8+4*y^7+y^6+4*y^5+2*y^4+4*y^3+2*y^2+3*y+2)*(y^2+y+1)*y*x^8+2*(y^3+y^2+y+3)*(y+1)*(y^11+4*y^9+2*y^8+4*y^7+3*y^6+4*y^4+3*y^3+4*y^2+1)*x^7+3*(y^4+2*y^3+2*y^2+y+1)*(y^2+y+2)^2*(y^5+3*y^4+y^3+y^2+3*y+3)*y*(y+4)*x^6+4*(y^5+3*y^4+2*y^3+1)*(y^5+4*y^4+4*y^3+y^2+4)*y*(y^2+2*y+4)*(y+4)^2*x^5+(y+2)^2*(y+3)^2*y*(y+4)^5*x^4+3*(y^4+3*y^3+3*y+2)*y*(y+4)^5*x^3+3*(y+2)*y*(y+4)^5*(y^3+2*y^2+1)*x^2+2*(y^3+2*y^2+4*y+2)*y*(y+4)^6*x+(y^2+4*y+1)*y*(y+4)^7;
XF2 := Scheme(ProjectiveSpace(P), homogenize(f0));
FF  := FunctionField(XF2);
AFF := AlgorithmicFunctionField(FF);
tt := AFF ! BaseField(AFF).1;  ww := AFF.1;

d := TotalDegree(f0);
printf "total degree of f: %o\n", d;
// f(x,y) with x = w/t, y = 1/t must vanish in AFF
ok1 := Evaluate(f0, [ww/tt, 1/tt]) eq 0;
printf "f(w/t, 1/t) = 0 in AFF : %o\n", ok1;
// Phi(t,w) = t^d f(w/t, 1/t), compared as polynomials in w over F5(t)
K := BaseField(AFF);  PW<W> := PolynomialRing(K);
target := PW ! (K.1^d * Evaluate(f0, [W/K.1, 1/K.1]));
Phi := DefiningPolynomial(AFF);
ok2 := Phi eq target;
printf "DefiningPolynomial(AFF) = t^%o * f(w/t, 1/t) : %o\n", d, ok2;

// ---- Phi as a polynomial in named variables t, w over F5 ----
R<t,w> := PolynomialRing(GF(5), 2);
cs := Coefficients(Phi);                        // cs[j+1] = coefficient of w^j, in F5(t)
okden := forall{ c : c in cs | Denominator(c) eq 1 };
printf "all coefficients of Phi lie in F5[t] : %o\n", okden;
PhiR := &+[ Evaluate(Numerator(cs[j]), t) * w^(j-1) : j in [1..#cs] ];
// independent recomputation in F5(t,w):  t^d f(w/t, 1/t)
FR := FieldOfFractions(R);
ok3 := PhiR eq R ! ((FR!t)^d * Evaluate(f0, [(FR!w)/(FR!t), 1/(FR!t)]));
printf "Phi(t,w) in F5[t,w] equals t^%o f(w/t,1/t) : %o\n", d, ok3;
n := Degree(PhiR, w);
printf "deg_w Phi = %o, leading coefficient %o, total degree %o\n",
       n, Coefficient(PhiR, w, n), TotalDegree(PhiR);

// ---- write X1_42_Phi.txt (readable, and loadable in Magma) ----
LINES := [];
for j in [n..0 by -1] do
    c := Coefficient(PhiR, w, j);
    if c eq 0 then continue; end if;
    mon := j eq 0 select "" else (j eq 1 select "*w" else Sprintf("*w^%o", j));
    Append(~LINES, Sprintf("  %o(%o)%o", #LINES eq 0 select "  " else "+ ", c, mon));
end for;
LINES[#LINES] cat:= ";";

OUT := "X1_42_Phi.txt";
PrintFile(OUT, "// X1_42_Phi.txt: defining polynomial Phi(t,w) of F5(X), X = X1(42)/F5." : Overwrite := true);
PrintFile(OUT, "// Variables: t = BaseField(AFF).1 = 1/y,  w = AFF.1 = x/y,");
PrintFile(OUT, "//   x, y = coordinates of the model f(x,y) = 0 (FFFc42.txt mod 5).");
PrintFile(OUT, "// F5(X) = F5(t)[w]/(Phi),  Phi = DefiningPolynomial(AFF)"); 
PrintFile(OUT, "//   = t^23 * f(w/t, 1/t).");
PrintFile(OUT, Sprintf("// deg_w Phi = %o, leading coeff. %o, total degree %o.",
                       n, Coefficient(PhiR, w, n), TotalDegree(PhiR)));
PrintFile(OUT, "// Written as a polynomial in w with coefficients in F5[t].");
PrintFile(OUT, "// Load in Magma with:  load \"X1_42_Phi.txt\";");
PrintFile(OUT, "R<t,w> := PolynomialRing(GF(5), 2);");
PrintFile(OUT, "Phi :=");
for l in LINES do PrintFile(OUT, l); end for;
printf "wrote %o\n", OUT;

assert ok1 and ok2 and okden and ok3;
print "CONFIRMED: t = 1/y and w = x/y.";
quit;
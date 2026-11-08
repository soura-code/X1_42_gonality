/*
  X1_42_orbits_check.m  --  check that the saved orbit files can be reused.
  For each saved degree it
    * reads X1_42_orbits_d<v>.txt (LoadOrbitFile),
    * recomputes Places(AFF,v) and checks the numbering against the file
      (PlacesFromOrbitFile),
    * finds every place again from its stored data alone (FindPlace) and checks
      that it is the same place,
    * prints the orbit sizes.
  Usage:  magma -b X1_42_orbits_check.m            (degrees 1..5)
          magma -b degs:=4,5 X1_42_orbits_check.m
*/
if assigned degs then DEGS := [ StringToInteger(s) : s in Split(degs, ",") ];
else DEGS := [1,2,3,4,5]; end if;

A<x,y> := PolynomialRing(GF(5),2);
P,homogenize := Homogenization(A);
f := x^16*y^4+4*(y+2)*y^4*x^15+(y^3+3*y^2+4*y+1)*(y+1)^2*(y+3)*y^3*x^14+4*(y^6+3*y^5+4*y^2+2*y+2)*(y+3)*y^3*x^13+3*(y^4+2*y^2+2*y+1)*(y^3+4*y^2+4*y+4)*(y+3)*y^2*x^12+2*(y+1)^2*y^2*(y^6+4*y^5+4*y^3+3*y^2+1)*x^11+3*(y^8+y^7+3*y^6+3*y^5+4*y^3+y^2+2*y+1)*y^2*x^10+4*(y^7+y^6+4*y^5+4*y^4+4*y^3+2*y+2)*(y^5+2*y^4+2*y^2+4*y+4)*y^2*x^9+(y+2)*(y^11+3*y^10+2*y^9+4*y^8+4*y^7+y^6+4*y^5+2*y^4+4*y^3+2*y^2+3*y+2)*(y^2+y+1)*y*x^8+2*(y^3+y^2+y+3)*(y+1)*(y^11+4*y^9+2*y^8+4*y^7+3*y^6+4*y^4+3*y^3+4*y^2+1)*x^7+3*(y^4+2*y^3+2*y^2+y+1)*(y^2+y+2)^2*(y^5+3*y^4+y^3+y^2+3*y+3)*y*(y+4)*x^6+4*(y^5+3*y^4+2*y^3+1)*(y^5+4*y^4+4*y^3+y^2+4)*y*(y^2+2*y+4)*(y+4)^2*x^5+(y+2)^2*(y+3)^2*y*(y+4)^5*x^4+3*(y^4+3*y^3+3*y+2)*y*(y+4)^5*x^3+3*(y+2)*y*(y+4)^5*(y^3+2*y^2+1)*x^2+2*(y^3+2*y^2+4*y+2)*y*(y+4)^6*x+(y^2+4*y+1)*y*(y+4)^7 ;
f := homogenize(f);
XF2 := Scheme(ProjectiveSpace(P),f);
FF := FunctionField(XF2);
AFF := AlgorithmicFunctionField(FF);

load "X1_42_orbits_load.m";

for v in DEGS do
    printf "\n===== degree %o =====\n", v;
    ok, D := LoadOrbitFile(v);
    if not ok then printf "FAIL  degree %o: file could not be read\n", v; continue; end if;
    sizes := [ #o : o in D`orbits ];
    printf "read: %o places, %o orbits, orbit sizes %o\n", D`count, #D`orbits,
           [ <s, #[ z : z in sizes | z eq s ]> : s in Sort(Setseq(Set(sizes))) ];
    t0 := Cputime();
    ok, PL := PlacesFromOrbitFile(D);
    printf "%o  degree %o: Places(AFF,%o) numbered as in the file (%o s)\n",
           ok select "PASS" else "FAIL", v, v, Round(Cputime(t0));
    if not ok then continue; end if;
    t0 := Cputime(); nfound := 0; nsame := 0;
    for i in [1..D`count] do
        okf, Q := FindPlace(D`coeffs[i], D`keys[i], v);
        if okf then nfound +:= 1; if Q eq PL[i] then nsame +:= 1; end if; end if;
    end for;
    printf "%o  degree %o: FindPlace recovered %o of %o places, all correct: %o (%o s)\n",
           (nsame eq nfound) select "PASS" else "FAIL", v, nfound, D`count, nsame eq nfound, Round(Cputime(t0));
    if nfound lt D`count then
        printf "      (%o places share their stored data with another place; use PL[i] for those)\n",
               D`count - nfound;
    end if;
end for;
quit;

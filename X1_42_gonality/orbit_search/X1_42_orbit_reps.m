/*
  X1_42_orbit_reps.m  --  one representative place for every orbit of
  <sigma_4> on the places of degree v, from the saved files X1_42_orbits_d<v>.txt.
  Usage:  magma -b X1_42_orbit_reps.m              (degrees 2,3,5)
          magma -b degs:=2,3,5 X1_42_orbit_reps.m
  Or inside your own session:  keep := true; load "X1_42_orbit_reps.m";
  (without keep the script quits at the end).  Afterwards REPS[v] is the list of
  representative places of degree v (a list; REPS[v][k] lies in orbit k,
  or is false if no explicit description was found), and
  REPIDX[v] the list of their indices in Places(AFF,v).

  Each representative is found from its stored data (polynomial below + key)
  with FindPlace, which does NOT depend on the order of Places(AFF,v).  Only if
  that does not single out the place, Places(AFF,v) is recomputed and the index
  is used; that is accepted only if the numbering check (keys + checksum) passes.

  For every representative it then computes an EXPLICIT description as a point
  (X1_42_reps_load.m): coordinates u1(P) = A(r), u2(P) = C(r) with m(r) = 0,
  m = minimal polynomial over F5 of a primitive element w, plus the minimal
  polynomials of u1(P), u2(P).  It checks that this description alone recovers
  the place (FindPlaceByPoint).  If the member with the smallest index has no
  such description (e.g. it is a branch at a singular point of the model that
  even the tangent slope does not separate), the next member of the same orbit
  is tried.  Writes
      X1_42_orbit_reps_d<v>.txt      (reload with LoadRepFile(v))
  Needs X1_42_orbits_load.m and X1_42_reps_load.m.
*/
if assigned degs then DEGS := [ StringToInteger(z) : z in Split(degs, ",") ];
else DEGS := [2,3,5]; end if;

A<x,y> := PolynomialRing(GF(5),2);
P,homogenize := Homogenization(A);
f := x^16*y^4+4*(y+2)*y^4*x^15+(y^3+3*y^2+4*y+1)*(y+1)^2*(y+3)*y^3*x^14+4*(y^6+3*y^5+4*y^2+2*y+2)*(y+3)*y^3*x^13+3*(y^4+2*y^2+2*y+1)*(y^3+4*y^2+4*y+4)*(y+3)*y^2*x^12+2*(y+1)^2*y^2*(y^6+4*y^5+4*y^3+3*y^2+1)*x^11+3*(y^8+y^7+3*y^6+3*y^5+4*y^3+y^2+2*y+1)*y^2*x^10+4*(y^7+y^6+4*y^5+4*y^4+4*y^3+2*y+2)*(y^5+2*y^4+2*y^2+4*y+4)*y^2*x^9+(y+2)*(y^11+3*y^10+2*y^9+4*y^8+4*y^7+y^6+4*y^5+2*y^4+4*y^3+2*y^2+3*y+2)*(y^2+y+1)*y*x^8+2*(y^3+y^2+y+3)*(y+1)*(y^11+4*y^9+2*y^8+4*y^7+3*y^6+4*y^4+3*y^3+4*y^2+1)*x^7+3*(y^4+2*y^3+2*y^2+y+1)*(y^2+y+2)^2*(y^5+3*y^4+y^3+y^2+3*y+3)*y*(y+4)*x^6+4*(y^5+3*y^4+2*y^3+1)*(y^5+4*y^4+4*y^3+y^2+4)*y*(y^2+2*y+4)*(y+4)^2*x^5+(y+2)^2*(y+3)^2*y*(y+4)^5*x^4+3*(y^4+3*y^3+3*y+2)*y*(y+4)^5*x^3+3*(y+2)*y*(y+4)^5*(y^3+2*y^2+1)*x^2+2*(y^3+2*y^2+4*y+2)*y*(y+4)^6*x+(y^2+4*y+1)*y*(y+4)^7 ;
f := homogenize(f);
XF2 := Scheme(ProjectiveSpace(P),f);
FF := FunctionField(XF2);
AFF := AlgorithmicFunctionField(FF);

load "X1_42_orbits_load.m";
load "X1_42_reps_load.m";

REPS := AssociativeArray(Integers()); REPIDX := AssociativeArray(Integers());
for v in DEGS do
    ok, D := LoadOrbitFile(v);
    error if not ok, Sprintf("cannot read X1_42_orbits_d%o.txt", v);
    printf "\n===== degree %o: %o places, %o orbits =====\n", v, D`count, #D`orbits;
    reps := [* *]; PL := []; havePL := false; plfail := false; lines := []; nbad := 0; repidx := [];
    for k in [1..#D`orbits] do
        // try the members of the orbit in increasing order of index until one has an
        // explicit description that determines it (the smallest index is tried first)
        best := 0; bestQ := 0; bestPD := 0; bestn := 0;
        for i in Sort(D`orbits[k]) do
            okf, Q := FindPlace(D`coeffs[i], D`keys[i], v);
            if not okf then
                if plfail then continue; end if;
                if not havePL then
                    okp, PL := PlacesFromOrbitFile(D);
                    if not okp then plfail := true; continue; end if;   // numbering differs: only FindPlace is safe
                    havePL := true;
                end if;
                Q := PL[i];
            end if;
            okp, PD := PointData(Q);
            if not okp then continue; end if;
            n, Lm := FindPlaceByPoint(PD, v);
            if n eq 1 and Lm[1] eq Q then best := i; bestQ := Q; bestPD := PD; bestn := n; break; end if;
        end for;
        if best eq 0 then
            printf "orbit %o (size %o): NO member has an explicit description that determines it\n", k, #D`orbits[k];
            Append(~lines, Sprintf("# orbit %o: no explicit description found (use index %o)", k, D`reps[k]));
            Append(~reps, false); Append(~repidx, D`reps[k]); nbad +:= 1; continue;
        end if;
        Append(~reps, bestQ); Append(~repidx, best);
        printf "orbit %o (size %o): representative = place %o%o\n  %o\n", k, #D`orbits[k], best,
               best eq D`reps[k] select "" else Sprintf(" (smallest index %o has no explicit description)", D`reps[k]),
               PointString(bestPD);
        cs := func< x | &cat[ Sprintf(" %o", z) : z in x ] cat "" >;
        wtxt := ""; for t in bestPD`w do wtxt cat:= Sprintf(" %o:%o:%o", t[1], t[2], t[3]); end for;
        extxt := bestPD`ex[1] eq 0 select " 0" else Sprintf(" %o %o %o", bestPD`ex[1], bestPD`ex[2], bestPD`ex[3]);
        Append(~lines, Sprintf("# orbit %o: %o", k, PointString(bestPD)));
        Append(~lines, Sprintf("REP %o %o %o | chart %o | w%o | g1%o | g2%o | m%o | A%o | C%o | shared %o | ex%o",
               k, #D`orbits[k], best, bestPD`chart, wtxt,
               cs(bestPD`g1), cs(bestPD`g2), cs(bestPD`m), cs(bestPD`A), cs(bestPD`C), bestn, extxt));
    end for;
    REPS[v] := reps; REPIDX[v] := repidx;
    OUT := Sprintf("X1_42_orbit_reps_d%o.txt", v);
    F := Open(OUT, "w");
    fprintf F, "X1_42_ORBITREPS v1 DEG %o GEN %o ORBITS %o\n", v, D`gen, #D`orbits;
    fprintf F, "# One representative per orbit of <sigma_%o> on the places of degree %o.\n", D`gen, v;
    fprintf F, "# REP k size i | chart ch | w c:i:j ... | g1 | g2 | m | A | C | shared n | ex e\n";
    fprintf F, "#   i = index in Places(AFF,%o); chart 1..4 = (u1,u2) = (t,y),(1/t,y),(t,1/y),(1/t,1/y);\n", v;
    fprintf F, "#   w = sum c*u1^i*u2^j; m = min. poly of w(P); u1(P) = A(w(P)), u2(P) = C(w(P));\n";
    fprintf F, "#   g1, g2 = min. polys of u1(P), u2(P); coefficients constant term first;\n";
    fprintf F, "#   shared = number of places of degree %o at this point (1 = point determines the place);\n", v;
    fprintf F, "#   ex 0 = plain point; ex 1 a c / ex 2 a c = branch at the singular F5-point (a,c) with\n";
    fprintf F, "#   slope w = (u2-c)/(u1-a) resp. (u1-a)/(u2-c), m = its minimal polynomial.\n";
    for l in lines do fprintf F, "%o\n", l; end for;
    fprintf F, "END %o\n", #D`orbits;
    delete F;
    printf "%o  degree %o: wrote %o; representatives not determined by their point alone: %o\n",
           nbad eq 0 select "PASS" else "NOTE", v, OUT, nbad;
end for;
if not assigned keep then quit; end if;

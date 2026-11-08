/*
  X1_42_reps_load.m  --  explicit description of a place as a point, and
  recovery of the place from that description alone.
  ------------------------------------------------------------------
  `load` after AFF is defined.

  A place P of degree v is described by
    chart   (u1,u2) = (t,y), (1/t,y), (t,1/y) or (1/t,1/y): the first pair
            regular at P  (charts 1..4),
    w       a primitive element w = sum c*u1^i*u2^j (small polynomial in the
            coordinates) whose residue generates the residue field F_{5^v},
    m(x)    the minimal polynomial over F5 of w(P)  (degree v),
    A(x), C(x)  polynomials of degree < v with
                u1(P) = A(w(P)),   u2(P) = C(w(P)),
    g1, g2  minimal polynomials over F5 of u1(P) and u2(P).
  So P "is" the point (u1,u2) = (A(r), C(r)) with m(r) = 0, i.e. P is the
  place where the affine plane model passes through that point (up to
  Galois conjugation, which gives the same place).  If several places of
  degree v lie over the same point (a singular point of the model), the
  description is not enough; `shared` records how many.

  Singular points with F5-rational coordinates.  If u1(P) = a and u2(P) = c lie
  in F5 but deg P > 1, no polynomial in u1, u2 generates the residue field
  (the point is a singular point of the model and P is one of its branches).
  Then the branch is described by its tangent slope: the extra function
      ex = 1 :  s = (u2 - c)/(u1 - a)      or     ex = 2 :  s = (u1 - a)/(u2 - c)
  with w = s, m = minimal polynomial of s(P), A = a, C = c (constants).
  The place is then: the branch through (a,c) whose slope is a root of m.

    ok, R := PointData(Q)       R = record of the data above (ok = false if no
                                primitive w was found in the small search list)
    n, L  := FindPlaceByPoint(R, v)   the places of degree v with this point
                                (n = #L; n = 1 means the place is recovered)
    s     := PointString(R)     human-readable "t = ..., y = ..., m(w) = 0"
  ------------------------------------------------------------------
*/
RP_TT := AFF!(BaseField(AFF).1);
RP_YY := AFF.1;
RP_CHARTS := [ <"t","y">, <"1/t","y">, <"t","1/y">, <"1/t","1/y"> ];
RP_ChartFns := function(ch)
    case ch:
        when 1: return RP_TT, RP_YY;
        when 2: return 1/RP_TT, RP_YY;
        when 3: return RP_TT, 1/RP_YY;
        when 4: return 1/RP_TT, 1/RP_YY;
    end case;
end function;
// candidate primitive elements, as lists of monomials <coef, i, j> = coef*u1^i*u2^j
RP_WLIST := [ [<1,1,0>], [<1,0,1>] ]
         cat [ [<1,1,0>,<l,0,1>] : l in [1..4] ]
         cat [ [<1,1,0>,<l,0,2>] : l in [1..4] ]
         cat [ [<1,0,1>,<l,2,0>] : l in [1..4] ]
         cat [ [<1,1,1>,<l,1,0>,<mu,0,1>] : l in [0..4], mu in [0..4] ]
         cat [ [<1,1,0>,<l,0,1>,<mu,1,1>] : l in [0..4], mu in [1..4] ];

RP_Rec := recformat< chart, w, g1, g2, m, A, C, ex >;   // ex = <0> or <type, a, c>

RP_Value := function(g, Q)
    try
        if Valuation(g, Q) lt 0 then return false, 0; end if;
        return true, ResidueClassField(Q)!Evaluate(g, Q);
    catch e
        return false, 0;
    end try;
end function;
RP_Coeffs := func< p | [ Integers()!c : c in Coefficients(p) ] >;
RP_EvalW := func< ws, a, c | &+[ (Parent(a)!m[1]) * a^m[2] * c^m[3] : m in ws ] >;
RP_SlopeFn := function(ch, ex)
    u1, u2 := RP_ChartFns(ch);
    if ex[1] eq 1 then return (u2 - ex[3])/(u1 - ex[2]); end if;
    return (u1 - ex[2])/(u2 - ex[3]);
end function;

PointData := function(Q)
    v := Degree(Q);
    for ch in [1..4] do
        u1, u2 := RP_ChartFns(ch);
        ok1, a := RP_Value(u1, Q); if not ok1 then continue; end if;
        ok2, c := RP_Value(u2, Q); if not ok2 then continue; end if;
        R := Parent(a); k := PrimeField(R);
        for ws in RP_WLIST do
            w := RP_EvalW(ws, a, c);
            mw := MinimalPolynomial(w, k);
            if Degree(mw) ne v then continue; end if;
            M := Matrix(k, [ Eltseq(w^i, k) : i in [0..v-1] ]);
            Avec := Solution(M, Vector(k, Eltseq(a, k)));
            Cvec := Solution(M, Vector(k, Eltseq(c, k)));
            return true, rec< RP_Rec | chart := ch, w := ws,
                     g1 := RP_Coeffs(MinimalPolynomial(a, k)), g2 := RP_Coeffs(MinimalPolynomial(c, k)),
                     m := RP_Coeffs(mw),
                     A := [ Integers()!Avec[i] : i in [1..v] ], C := [ Integers()!Cvec[i] : i in [1..v] ],
                     ex := <0> >;
        end for;
        // no primitive polynomial in the coordinates: try the tangent slope at an F5-rational point
        if Degree(MinimalPolynomial(a, k)) eq 1 and Degree(MinimalPolynomial(c, k)) eq 1 then
            a0 := Integers()!(Eltseq(a, k)[1]); c0 := Integers()!(Eltseq(c, k)[1]);
            for typ in [1, 2] do
                ex := <typ, a0, c0>;
                oks, sv := RP_Value(RP_SlopeFn(ch, ex), Q);
                if not oks then continue; end if;
                ms := MinimalPolynomial(sv, k);
                if Degree(ms) ne v then continue; end if;
                return true, rec< RP_Rec | chart := ch, w := [],
                         g1 := [ (-a0) mod 5, 1 ], g2 := [ (-c0) mod 5, 1 ], m := RP_Coeffs(ms),
                         A := [a0] cat [ 0 : i in [2..v] ], C := [c0] cat [ 0 : i in [2..v] ], ex := ex >;
            end for;
            printf "PointData: place at the F5-point (%o,%o) of chart %o; slope does not separate it\n", a0, c0, ch;
        else
            printf "PointData: residues %o, %o (chart %o): no primitive element in the search list\n", a, c, ch;
        end if;
        return false, 0;
    end for;
    return false, 0;
end function;

FindPlaceByPoint := function(D, v)
    u1, u2 := RP_ChartFns(D`chart);
    g1fn := &+[ D`g1[j+1] * u1^j : j in [0..#D`g1-1] ];
    try
        cand := [ Q : Q in Zeros(g1fn) | Degree(Q) eq v ];
    catch e
        return 0, [];
    end try;
    L := [];
    slope := D`ex[1] ne 0;
    if slope then sfn := RP_SlopeFn(D`chart, D`ex); end if;
    for Q in cand do
        ok1, a := RP_Value(u1, Q); ok2, c := RP_Value(u2, Q);
        if not (ok1 and ok2) then continue; end if;
        R := Parent(a);
        mR := Polynomial(R, [ R!z : z in D`m ]);
        if slope then
            if a eq R!D`ex[2] and c eq R!D`ex[3] then
                oks, sv := RP_Value(sfn, Q);
                if oks and Evaluate(mR, sv) eq 0 then Append(~L, Q); end if;
            end if;
            continue;
        end if;
        AR := Polynomial(R, [ R!z : z in D`A ]);
        CR := Polynomial(R, [ R!z : z in D`C ]);
        for rt in Roots(mR) do
            r := rt[1];
            if Evaluate(AR, r) eq a and Evaluate(CR, r) eq c then Append(~L, Q); break; end if;
        end for;
    end for;
    return #L, L;
end function;

PointString := function(D)
    Fx<w> := PolynomialRing(GF(5));
    n1 := RP_CHARTS[D`chart][1]; n2 := RP_CHARTS[D`chart][2];
    if D`ex[1] ne 0 then
        sstr := D`ex[1] eq 1 select Sprintf("(%o - %o)/(%o - %o)", n2, D`ex[3], n1, D`ex[2])
                                else Sprintf("(%o - %o)/(%o - %o)", n1, D`ex[2], n2, D`ex[3]);
        return Sprintf("%o = %o,  %o = %o (singular point of the model);  branch with slope w = %o satisfying %o = 0",
                       n1, D`ex[2], n2, D`ex[3], sstr, Fx!D`m);
    end if;
    wstr := "";
    for m in D`w do
        term := Sprintf("%o", m[1]);
        if m[2] gt 0 then term cat:= "*(" cat n1 cat ")" cat (m[2] gt 1 select Sprintf("^%o", m[2]) else ""); end if;
        if m[3] gt 0 then term cat:= "*(" cat n2 cat ")" cat (m[3] gt 1 select Sprintf("^%o", m[3]) else ""); end if;
        if m[1] ne 0 then wstr cat:= (#wstr gt 0 select " + " else "") cat term; end if;
    end for;
    return Sprintf("%o = %o,  %o = %o,  where w = %o satisfies %o = 0",
                   n1, Fx!D`A, n2, Fx!D`C, wstr, Fx!D`m);
end function;

/*
  File X1_42_orbit_reps_d<v>.txt (written by X1_42_orbit_reps.m), one line per orbit:
    REP k size i | chart ch | w c:i:j c:i:j ... | g1 ... | g2 ... | m ... | A ... | C ... | shared n | ex 0  (or ex type a c)
  (all polynomials as coefficient lists, constant term first), preceded by a
  comment line with PointString.
    ok, L := LoadRepFile(v)   L = list of <k, size, index, record, shared>
*/
LoadRepFile := function(v)
    name := Sprintf("X1_42_orbit_reps_d%o.txt", v);
    try F := Open(name, "r"); catch e printf "LoadRepFile: cannot open %o\n", name; return false, []; end try;
    L := []; endseen := false;
    while true do
        s := Gets(F);
        if IsEof(s) then break; end if;
        if #s eq 0 or s[1] eq "#" then continue; end if;
        if #s ge 3 and s[1..3] eq "END" then endseen := true; continue; end if;
        if #s lt 3 or s[1..3] ne "REP" then continue; end if;
        f := [ Split(z, " ") : z in Split(s, "|") ];      // 9 fields
        hd := f[1];
        ws := [];
        for z in f[3][2..#f[3]] do
            u := Split(z, ":");
            Append(~ws, <StringToInteger(u[1]), StringToInteger(u[2]), StringToInteger(u[3])>);
        end for;
        tail := func< x | [ StringToInteger(z) : z in x[2..#x] ] >;
        D := rec< RP_Rec | chart := StringToInteger(f[2][2]), w := ws,
                  g1 := tail(f[4]), g2 := tail(f[5]), m := tail(f[6]), A := tail(f[7]), C := tail(f[8]),
                  ex := <0> >;
        if #f ge 10 and #f[10] ge 4 and f[10][2] ne "0" then
            D`ex := < StringToInteger(f[10][2]), StringToInteger(f[10][3]), StringToInteger(f[10][4]) >;
        end if;
        Append(~L, < StringToInteger(hd[2]), StringToInteger(hd[3]), StringToInteger(hd[4]), D, StringToInteger(f[9][2]) >);
    end while;
    delete F;
    if not endseen then printf "LoadRepFile: %o incomplete\n", name; return false, []; end if;
    return true, L;
end function;

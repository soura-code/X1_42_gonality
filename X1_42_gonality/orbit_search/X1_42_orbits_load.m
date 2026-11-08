/*
  X1_42_orbits_load.m  --  reload the data written by X1_42_orbits.m
  ------------------------------------------------------------------
  `load` this file after AFF (the function field of X1(42)/F5) is defined.

    ok, D := LoadOrbitFile(v);
        reads X1_42_orbits_d<v>.txt.  D is a record with fields
          deg, count, gen, checksum,
          coeffs  : for each place, coefficients of the polynomial of F5(t) below it ([] = infinity)
          keys    : for each place, its key (0 = not available)
          img     : img[i] = index of sigma_4(place i)
          orbits  : list of orbits, each a cycle [i, sigma(i), sigma^2(i), ...]
          reps    : the smallest index in each orbit
          orbitof : orbitof[i] = number of the orbit containing place i
        Pure file reading: no curve computations.

    ok, PL := PlacesFromOrbitFile(D);
        computes Places(AFF, v) and checks, place by place, that the stored data
        (polynomial below, key) and the checksum agree, i.e. that this Magma
        session numbers the places exactly as in the file.  Use PL[i] for place i.

    ok, Q := FindPlace(coeffs, key, v);
        finds ONE place of degree v from its stored data without computing all
        places of degree v (via the zeros of p(t), or the poles of t).
        ok = false if the data does not single out a place.
  ------------------------------------------------------------------
*/
OB_TT := AFF!(BaseField(AFF).1);
OB_YY := AFF.1;

OB_EncPoly := function(g)
    c := Coefficients(g); e := 0;
    for i in [#c..1 by -1] do e := 5*e + (Integers()!c[i]); end for;
    return e + 5^(#c);
end function;
OB_ValueAt := function(g, Q)
    try
        if Valuation(g, Q) lt 0 then return false, 0; end if;
        R := ResidueClassField(Q);
        return true, R!Evaluate(g, Q);
    catch e
        return false, 0;
    end try;
end function;
OB_MinPoly := function(a)
    k := PrimeField(Parent(a));
    try
        return MinimalPolynomial(a, k);
    catch e1
        dummy := 0;
    end try;
    return MinimalPolynomial(a);
end function;
OB_Key := function(Q)
    okt, at := OB_ValueAt(OB_TT, Q);
    if not okt then return [], 0, 0; end if;
    mp := OB_MinPoly(at);
    coeffs := [ Integers()!c : c in Coefficients(mp) ];
    b := OB_EncPoly(mp);
    oky, ay := OB_ValueAt(OB_YY, Q);
    if not oky then return coeffs, b, 0; end if;
    try
        polys := [ OB_MinPoly(at), OB_MinPoly(ay), OB_MinPoly(at+ay), OB_MinPoly(at+2*ay) ];
    catch e
        return coeffs, b, 0;
    end try;
    E := 5^12; key := 0;
    for g in polys do key := key*E + OB_EncPoly(g); end for;
    return coeffs, b, key;
end function;

OrbitData := recformat< deg, count, gen, checksum, coeffs, keys, img, orbits, reps, orbitof >;

LoadOrbitFile := function(v)
    name := Sprintf("X1_42_orbits_d%o.txt", v);
    try
        F := Open(name, "r");
    catch e
        printf "LoadOrbitFile: cannot open %o\n", name; return false, 0;
    end try;
    head := Split(Gets(F), " ");
    if head[1] ne "X1_42_ORBITS" then printf "LoadOrbitFile: %o has an unknown format\n", name; return false, 0; end if;
    hd := AssociativeArray();
    for i in [3..#head-1 by 2] do hd[head[i]] := head[i+1]; end for;
    n := StringToInteger(hd["COUNT"]);
    coeffs := [ [] : i in [1..n] ]; keys := [ 0 : i in [1..n] ]; img := [ 0 : i in [1..n] ];
    orbits := []; reps := []; endseen := false;
    while true do
        s := Gets(F);
        if IsEof(s) then break; end if;
        if #s eq 0 or s[1] eq "#" then continue; end if;
        if s[1] eq "P" then
            parts := Split(s, "|");                       // "P i ", " c0 c1 ... ", " key"
            i := StringToInteger(Split(parts[1], " ")[2]);
            coeffs[i] := [ StringToInteger(z) : z in Split(parts[2], " ") ];
            keys[i] := StringToInteger(Split(parts[3], " ")[1]);
        else
            t := Split(s, " ");
            case t[1]:
                when "I": img[StringToInteger(t[2])] := StringToInteger(t[3]);
                when "O": Append(~orbits, [ StringToInteger(z) : z in t[5..#t] ]);
                when "R": Append(~reps, StringToInteger(t[3]));
                when "END": endseen := true;
            end case;
        end if;
    end while;
    delete F;
    if not endseen then printf "LoadOrbitFile: %o is incomplete\n", name; return false, 0; end if;
    if Set(img) ne {1..n} or &+[ #o : o in orbits ] ne n then
        printf "LoadOrbitFile: %o is inconsistent\n", name; return false, 0;
    end if;
    orbitof := [ 0 : i in [1..n] ];
    for k in [1..#orbits] do for i in orbits[k] do orbitof[i] := k; end for; end for;
    D := rec< OrbitData | deg := v, count := n, gen := StringToInteger(hd["GEN"]),
              checksum := StringToInteger(hd["CHECKSUM"]), coeffs := coeffs, keys := keys,
              img := img, orbits := orbits, reps := reps, orbitof := orbitof >;
    return true, D;
end function;

PlacesFromOrbitFile := function(D)
    v := D`deg;
    PL := (v eq 1) select Places(AFF, 1) else Places(AFF, v);
    if #PL ne D`count then
        printf "PlacesFromOrbitFile: %o places of degree %o, file has %o\n", #PL, v, D`count;
        return false, [];
    end if;
    MODP := 2305843009213693951; CHK := 0; bad := 0;
    for i in [1..#PL] do
        c, b, k := OB_Key(PL[i]);
        if c ne D`coeffs[i] or k ne D`keys[i] then bad +:= 1; end if;
        CHK := (CHK*1000003 + (b mod 1000000007) + 7*(k mod 1000000007)) mod MODP;
    end for;
    if bad gt 0 or CHK ne D`checksum then
        printf "PlacesFromOrbitFile: %o places differ from the file (different place order?)\n", bad;
        return false, [];
    end if;
    return true, PL;
end function;

FindPlace := function(coeffs, key, v)
    try
        if #coeffs eq 0 then
            cand := [ Q : Q in Poles(OB_TT) | Degree(Q) eq v ];
        else
            p := &+[ coeffs[j+1] * OB_TT^j : j in [0..#coeffs-1] ];
            cand := [ Q : Q in Zeros(p) | Degree(Q) eq v ];
        end if;
    catch e
        return false, 0;
    end try;
    m := [];
    for Q in cand do
        c, b, k := OB_Key(Q);
        if c eq coeffs and k eq key then Append(~m, Q); end if;
    end for;
    if #m eq 1 then return true, m[1]; end if;
    return false, 0;
end function;

/*
  X1_42_aut_load.m  --  helpers for storing / reloading automorphisms of X1(42)/F5
  ------------------------------------------------------------------
  `load` this file AFTER the curve has been built (AFF, plc1 defined) and AFTER
  the key functions (KU, KW, PlaceKey, RatKeyLine) are defined.

  Representation.  AFF is a simple extension of the rational function field
  R = F5(t) by y = AFF.1, of degree NDEG.  An F5-automorphism s of AFF is
  determined by the two elements s(t), s(y) of AFF, and every element of AFF
  is a sequence of NDEG elements of R (coefficients of 1, y, ..., y^(NDEG-1)),
  each a quotient num/den of polynomials over F5.  So an automorphism is
  stored as two lists of NDEG pairs (num coefficients, den coefficients),
  all small integers 0..4 -- no Magma-specific objects, safe for a database.

  Main functions:
    AutFromMap(s)          -> <s(t), s(y)>           (from a Magma map)
    ApplyAut(a, f)         -> a(f)                    (a = <s(t), s(y)>)
    ComposeAut(a, b)       -> a o b                   (apply b first)
    IsHomData(a)           -> true iff a respects the defining equation
    EltToData(g), DataToElt(d)     element <-> integer data
    WriteAutText(name, A)  write list A to the reloadable text file
    LoadAutData(name)      -> ok, list A  (verified: fingerprint, degree,
                              equation, closure under composition, identity)
  ------------------------------------------------------------------
*/
AUT_R  := BaseField(AFF);
error if Type(AUT_R) ne FldFunRat,
    "X1_42_aut_load.m expects AFF to be a simple extension of a rational function field.";
AUT_PR := PolynomialRing(GF(5));
TT     := AFF!(AUT_R.1);
YY     := AFF.1;
NDEG   := Degree(AFF);

RatToData := function(c)
    return [ [ Integers()!z : z in Coefficients(Numerator(c)) ],
             [ Integers()!z : z in Coefficients(Denominator(c)) ] ];
end function;

DataToRat := function(pr)
    return Evaluate(AUT_PR!pr[1], AUT_R.1) / Evaluate(AUT_PR!pr[2], AUT_R.1);
end function;

EltToData := function(g)
    cs := Eltseq(g);
    return [ RatToData(cs[i]) : i in [1..#cs] ];
end function;

DataToElt := function(d)
    return AFF ! [ DataToRat(pr) : pr in d ];
end function;

// image of c in R = F5(t) under t |-> z
RatAt := function(c, z)
    return Evaluate(Numerator(c), z) / Evaluate(Denominator(c), z);
end function;

ApplyAut := function(a, f)
    cs := Eltseq(f); res := AFF!0; yp := AFF!1;
    for i in [1..#cs] do
        if cs[i] ne 0 then res +:= RatAt(cs[i], a[1]) * yp; end if;
        yp *:= a[2];
    end for;
    return res;
end function;

ComposeAut := function(a, b)          // (a o b)(f) = a(b(f))
    return < ApplyAut(a, b[1]), ApplyAut(a, b[2]) >;
end function;

AutFromMap := function(s)
    return < s(TT), s(YY) >;
end function;

IsHomData := function(a)
    cs := Coefficients(DefiningPolynomial(AFF));
    s := AFF!0; yp := AFF!1;
    for i in [1..#cs] do
        if cs[i] ne 0 then s +:= RatAt(AUT_R!cs[i], a[1]) * yp; end if;
        yp *:= a[2];
    end for;
    return s eq 0;
end function;

AutIndex := function(A, a)            // position of a in A, 0 if absent
    for i in [1..#A] do
        if A[i][1] eq a[1] and A[i][2] eq a[2] then return i; end if;
    end for;
    return 0;
end function;

// ---------------------------------------------------------------------------
// File format (v2).  The file is read as a stream of whitespace-separated
// tokens, so line breaks (including any Magma line-wrapping) do not matter:
//   X1_42_AUT_F5 v2
//   DEG n
//   FINGERPRINT k_1 ... k_12 ENDFP        (keys of the rational places)
//   COUNT m
//   C i T|Y j : num_0 num_1 ... / den_0 den_1 ... ;     (one record per coefficient)
//   END
// Version-1 files (one record per line, no ";") are also accepted.
// ---------------------------------------------------------------------------
WriteTokens := procedure(F, toks)       // at most ~20 tokens per output line
    line := ""; c := 0;
    for t in toks do
        line cat:= t cat " "; c +:= 1;
        if c ge 20 or #line gt 60 then fprintf F, "%o\n", line; line := ""; c := 0; end if;
    end for;
    if #line gt 0 then fprintf F, "%o\n", line; end if;
end procedure;

WriteAutText := procedure(name, A)
    F := Open(name, "w");
    fprintf F, "X1_42_AUT_F5 v2\n";
    fprintf F, "DEG %o\n", NDEG;
    fprintf F, "FINGERPRINT\n";
    for t in Split(RatKeyLine(), " ") do
        if t ne "RAT" and #t gt 0 then fprintf F, "%o\n", t; end if;
    end for;
    fprintf F, "ENDFP\n";
    fprintf F, "COUNT %o\n", #A;
    for i in [1..#A] do
        for which in [1, 2] do
            d := EltToData(A[i][which]);
            for j in [1..#d] do
                toks := [ "C", Sprint(i), (which eq 1 select "T" else "Y"), Sprint(j), ":" ]
                        cat [ Sprint(z) : z in d[j][1] ] cat [ "/" ]
                        cat [ Sprint(z) : z in d[j][2] ] cat [ ";" ];
                WriteTokens(F, toks);
            end for;
        end for;
    end for;
    fprintf F, "END\n";
    Flush(F); delete F;
end procedure;

// tokens of a file, robust to line breaks and to Magma's "\"+newline continuations
FileTokens := function(name)
    txt := Read(name);
    parts := Split(txt, "\\");
    if #parts gt 1 then
        txt := parts[1];
        for k in [2..#parts] do
            q := parts[k];
            if #q gt 0 and q[1] eq "\n" then q := Substring(q, 2, #q - 1); end if;
            txt cat:= q;
        end for;
    end if;
    toks := [ t : t in Split(txt, " \n\t\r") | #t gt 0 ];
    // version-1 files print an empty coefficient list (zero numerator) as "[]",
    // glued to the following "/" as "[]/"; normalise those tokens.
    out := [];
    for t in toks do
        if t eq "[]" then continue; end if;
        if t eq "[]/" then Append(~out, "/"); continue; end if;
        Append(~out, t);
    end for;
    return out;
end function;

LoadFail := function(name, msg)
    printf "LoadAutData(%o): %o\n", name, msg; return false;
end function;

// returns ok, A.  On any problem: prints the reason and returns false (never errors),
// so callers can fall back to recomputing the group.
// LoadAutDataOpt(name, closure): closure = false skips the (slow, ~20 min) check that
// the list is closed under composition; the equation and identity are still checked.
LoadAutDataOpt := function(name, closure)
    try
        tk := FileTokens(name);
    catch e
        return false, [];
    end try;
    if #tk lt 6 or tk[1] ne "X1_42_AUT_F5" then return LoadFail(name, "unknown format"), []; end if;
    if tk[#tk] ne "END" then return LoadFail(name, "file incomplete (no END)"), []; end if;
    pos := 3;
    if tk[pos] ne "DEG" or StringToInteger(tk[pos+1]) ne NDEG then return LoadFail(name, "degree mismatch"), []; end if;
    pos +:= 2;
    if tk[pos] ne "FINGERPRINT" then return LoadFail(name, "no fingerprint"), []; end if;
    pos +:= 1;
    fp := [];
    while tk[pos] notin {"ENDFP", "COUNT"} do
        if tk[pos] ne "RAT" then Append(~fp, tk[pos]); end if;
        pos +:= 1;
    end while;
    if tk[pos] eq "ENDFP" then pos +:= 1; end if;
    mine := [ t : t in Split(RatKeyLine(), " ") | #t gt 0 and t ne "RAT" ];
    if fp ne mine then return LoadFail(name, "fingerprint mismatch (different model or key functions)"), []; end if;
    if tk[pos] ne "COUNT" then return LoadFail(name, "no COUNT"), []; end if;
    m := StringToInteger(tk[pos+1]); pos +:= 2;
    data := AssociativeArray(); seen := 0;
    while tk[pos] ne "END" do
        if tk[pos] eq ";" then pos +:= 1; continue; end if;
        if tk[pos] ne "C" or tk[pos+4] ne ":" then return LoadFail(name, Sprintf("bad record at token %o", pos)), []; end if;
        i := StringToInteger(tk[pos+1]); w := tk[pos+2] eq "T" select 1 else 2;
        j := StringToInteger(tk[pos+3]); pos +:= 5;
        num := [];
        while tk[pos] ne "/" do Append(~num, StringToInteger(tk[pos])); pos +:= 1; end while;
        pos +:= 1;
        den := [];
        while tk[pos] notin {";", "C", "END"} do Append(~den, StringToInteger(tk[pos])); pos +:= 1; end while;
        data[Sprintf("%o %o %o", i, w, j)] := [num, den];
        seen +:= 1;
    end while;
    if seen ne 2*m*NDEG then return LoadFail(name, Sprintf("expected %o records, found %o", 2*m*NDEG, seen)), []; end if;
    A := [ < DataToElt([ data[Sprintf("%o 1 %o", i, j)] : j in [1..NDEG] ]),
             DataToElt([ data[Sprintf("%o 2 %o", i, j)] : j in [1..NDEG] ]) > : i in [1..m] ];
    for a in A do
        if not IsHomData(a) then return LoadFail(name, "an entry does not respect the curve equation"), []; end if;
    end for;
    if AutIndex(A, <TT, YY>) eq 0 then return LoadFail(name, "identity missing"), []; end if;
    if closure then
        for a in A do for b in A do
            if AutIndex(A, ComposeAut(a, b)) eq 0 then return LoadFail(name, "not closed under composition"), []; end if;
        end for; end for;
    end if;
    return true, A;
end function;

LoadAutData := function(name)
    return LoadAutDataOpt(name, true);
end function;

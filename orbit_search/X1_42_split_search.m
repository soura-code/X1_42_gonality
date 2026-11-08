/*
  X1_42_split_search.m  --  STAGE C of the order-independent degree-11 search
  ------------------------------------------------------------------
  Tests the divisors of one job file written by X1_42_split.py:
      for every kept E_low in the file and (in groups H6..H9) every place Q of
      the high degree whose bucket b has BucketHash(b, GSIZE) = RANK,
      compute L(S + E_low [+ Q]) and look for a nonconstant function of degree
      <= 11  (S = sum of the 12 rational places).
  By Corollary 3.3 of the manuscript, gon_F5(X1(42)) >= 12 iff no job finds one.
  This test is UNCONDITIONAL (it does not assume anything about degree 10):
  any function of degree <= 11 is reported.

  Nothing depends on Magma's order of places:
    * degree 1  : places identified by their signature (X1_42_rational_labels.txt)
    * degrees 2-5: places identified by (bucket, key) from the K lines
    * degrees 6-9: Places(AFF,v) is computed, every place's bucket is computed
      (no key), and the job keeps those in its bucket class.
  Checks: numbers of places per degree, every needed label found exactly once,
  every rational signature found.  Any failure stops the job (status ERROR).

  Usage:  magma -b job:=k njobs:=N [dir:=X1_42_split_d11] X1_42_split_search.m
  Test mode: add  limit:=50  to stop after 50 divisors (status PARTIAL, which the
  checker never accepts).
  Needs X1_42_split_lib.m and X1_42_rational_labels.txt in the current folder.
  Output in DIR: status_job<k>_of_<N>.txt, hits_job<k>_of_<N>.txt,
                 FOUND.txt (shared; only if a function of degree <= 11 is found)
*/
SetColumns(0);
if assigned job   then JOB   := StringToInteger(job);   else JOB := 0; end if;
error if not assigned njobs, "give njobs:=N";
NJOBS := StringToInteger(njobs);
if assigned dir then DIR := dir; else DIR := "X1_42_split_d11"; end if;
if assigned limit then LIMIT := StringToInteger(limit); else LIMIT := 0; end if;   // test mode
JOBFILE    := Sprintf("%o/job%o_of_%o.txt", DIR, JOB, NJOBS);
STATUSFILE := Sprintf("%o/status_job%o_of_%o.txt", DIR, JOB, NJOBS);
OUTFILE    := Sprintf("%o/hits_job%o_of_%o.txt", DIR, JOB, NJOBS);
FOUNDFILE  := Sprintf("%o/FOUND.txt", DIR);
Dtarget := 11; LOGMAX := 12; REPORT := 2000;

WriteStatus := procedure(s)
    F := Open(STATUSFILE, "w"); fprintf F, "%o\n", s; delete F;
end procedure;
FileExists := function(name)
    try F := Open(name, "r"); delete F; return true; catch e return false; end try;
end function;

// ---------------- read the job file ----------------
F := Open(JOBFILE, "r");
head := Split(Gets(F), " ");
error if head[1] ne "SPLITJOB", "bad job file " cat JOBFILE;
HD := AssociativeArray();
for i in [3..#head-1 by 2] do HD[head[i]] := head[i+1]; end for;
GROUP := HD["GROUP"]; RANK := StringToInteger(HD["RANK"]); GSIZE := StringToInteger(HD["GSIZE"]);
HIGH := StringToInteger(HD["HIGHDEG"]); NLOW := StringToInteger(HD["NLOW"]);
error if StringToInteger(HD["JOB"]) ne JOB or StringToInteger(HD["NJOBS"]) ne NJOBS, "job file does not match job:= njobs:=";
PCOUNT := AssociativeArray();           // degree -> number of places
KLAB := AssociativeArray();             // "v b key" -> label
NEED := AssociativeArray();             // v -> labels needed
LOWS := [* *];                          // list of < partition string, [* <v,label>, ... *] >
curlam := ""; endn := -1;
while true do
    s := Gets(F);
    if IsEof(s) then break; end if;
    if #s eq 0 then continue; end if;
    t := Split(s, " ");
    case t[1]:
        when "C": PCOUNT[StringToInteger(t[2])] := StringToInteger(t[3]);
        when "K":
            v := StringToInteger(t[2]); lab := StringToInteger(t[3]);
            KLAB[Sprintf("%o %o %o", v, t[4], t[5])] := lab;
            if IsDefined(NEED, v) then z := NEED[v]; Include(~z, lab); NEED[v] := z;
            else NEED[v] := {lab}; end if;
        when "L": curlam := Substring(s, 3, #s - 2);
        when "D":
            E := [* *];
            for z in t[2..#t] do
                u := Split(z, ":");
                Append(~E, [ StringToInteger(u[1]), StringToInteger(u[2]) ]);
            end for;
            Append(~LOWS, [* curlam, E *]);
        when "END": endn := StringToInteger(t[2]);
    end case;
end while;
delete F;
error if endn ne NLOW or #LOWS ne NLOW, "job file incomplete";
printf "Job %o of %o: group %o, rank %o of %o, high degree %o, %o low parts\n",
       JOB, NJOBS, GROUP, RANK, GSIZE, HIGH, NLOW;

// ---------------- the curve ----------------
A<x,y> := PolynomialRing(GF(5),2);
P,homogenize := Homogenization(A);
f := x^16*y^4+4*(y+2)*y^4*x^15+(y^3+3*y^2+4*y+1)*(y+1)^2*(y+3)*y^3*x^14+4*(y^6+3*y^5+4*y^2+2*y+2)*(y+3)*y^3*x^13+3*(y^4+2*y^2+2*y+1)*(y^3+4*y^2+4*y+4)*(y+3)*y^2*x^12+2*(y+1)^2*y^2*(y^6+4*y^5+4*y^3+3*y^2+1)*x^11+3*(y^8+y^7+3*y^6+3*y^5+4*y^3+y^2+2*y+1)*y^2*x^10+4*(y^7+y^6+4*y^5+4*y^4+4*y^3+2*y+2)*(y^5+2*y^4+2*y^2+4*y+4)*y^2*x^9+(y+2)*(y^11+3*y^10+2*y^9+4*y^8+4*y^7+y^6+4*y^5+2*y^4+4*y^3+2*y^2+3*y+2)*(y^2+y+1)*y*x^8+2*(y^3+y^2+y+3)*(y+1)*(y^11+4*y^9+2*y^8+4*y^7+3*y^6+4*y^4+3*y^3+4*y^2+1)*x^7+3*(y^4+2*y^3+2*y^2+y+1)*(y^2+y+2)^2*(y^5+3*y^4+y^3+y^2+3*y+3)*y*(y+4)*x^6+4*(y^5+3*y^4+2*y^3+1)*(y^5+4*y^4+4*y^3+y^2+4)*y*(y^2+2*y+4)*(y+4)^2*x^5+(y+2)^2*(y+3)^2*y*(y+4)^5*x^4+3*(y^4+3*y^3+3*y+2)*y*(y+4)^5*x^3+3*(y+2)*y*(y+4)^5*(y^3+2*y^2+1)*x^2+2*(y^3+2*y^2+4*y+2)*y*(y+4)^6*x+(y^2+4*y+1)*y*(y+4)^7 ;
f := homogenize(f);
XF2 := Scheme(ProjectiveSpace(P),f);
FF := FunctionField(XF2);
AFF := AlgorithmicFunctionField(FF);
plc1 := Places(AFF,1);
assert Genus(AFF) eq 25 and #plc1 eq 12;
cuspsum := &+ plc1;
kbar := ConstantField(AFF);
DG0  := DivisorGroup(AFF)!0;
load "X1_42_split_lib.m";

PLACEOF := AssociativeArray();          // "v:label" -> place (degrees 1..5)

// ---- degree 1: signature labels ----
RF := Open("X1_42_rational_labels.txt", "r");
rh := Split(Gets(RF), " ");
error if rh[1] ne "X1_42_RATLABELS" or rh[2] ne "v2", "X1_42_rational_labels.txt missing or old (rerun X1_42_rational_labels.m)";
SIGLAB := AssociativeArray(); FLD := AssociativeArray(); NF := -1;
while true do
    s := Gets(RF);
    if IsEof(s) then break; end if;
    t := Split(s, " ");
    if #t ge 3 and t[1] eq "R" then SIGLAB[t[3]] := StringToInteger(t[2]); end if;
    if #t ge 2 and t[1] eq "NF" then NF := StringToInteger(t[2]); end if;
    if #t ge 3 and t[1] eq "F" then
        // the data starts after "F i "
        st := #t[1] + #t[2] + 3;
        FLD[StringToInteger(t[2])] := Substring(s, st, #s - st + 1);
    end if;
end while;
delete RF;
error if #Keys(SIGLAB) ne 12, "X1_42_rational_labels.txt does not list 12 labels";
error if NF lt 1 or Keys(FLD) ne {1..NF}, "X1_42_rational_labels.txt: separating functions incomplete";
FL := [ LineToFun(FLD[i]) : i in [1..NF] ];
for Q in plc1 do
    sg := SignatureString(Q, FL);
    if not IsDefined(SIGLAB, sg) then
        WriteStatus("ERROR: a rational place has an unknown signature"); error "unknown rational signature";
    end if;
    key := Sprintf("1:%o", SIGLAB[sg]);
    if IsDefined(PLACEOF, key) then
        WriteStatus("ERROR: two rational places with the same signature"); error "duplicate signature";
    end if;
    PLACEOF[key] := Q;
end for;
print "Rational places identified by their signatures.";

// ---- degrees 2..5: (bucket, key) ----
for v in Sort(Setseq(Keys(NEED))) do
    t1 := Cputime();
    PL := Places(AFF, v);
    if #PL ne PCOUNT[v] then
        WriteStatus(Sprintf("ERROR: %o places of degree %o, expected %o", #PL, v, PCOUNT[v]));
        error "wrong number of places";
    end if;
    found := {};
    for Q in PL do
        b, k := PlaceBucketKey(Q);
        ks := Sprintf("%o %o %o", v, b, k);
        if IsDefined(KLAB, ks) then
            lab := KLAB[ks];
            if lab in found then
                WriteStatus(Sprintf("ERROR: degree %o label %o matched twice", v, lab)); error "label matched twice";
            end if;
            Include(~found, lab);
            PLACEOF[Sprintf("%o:%o", v, lab)] := Q;
        end if;
    end for;
    if found ne NEED[v] then
        WriteStatus(Sprintf("ERROR: degree %o: %o of %o needed places not found", v, #(NEED[v] diff found), #NEED[v]));
        error "needed places not found";
    end if;
    printf "Degree %o: %o places, all %o needed labels found (%o s)\n", v, #PL, #found, Round(Cputime(t1));
end for;

// ---- high degree: bucket class ----
HIGHPL := []; BALL := 0; BKEPT := 0;
if HIGH gt 0 then
    t1 := Cputime();
    PL := Places(AFF, HIGH);
    if #PL ne PCOUNT[HIGH] then
        WriteStatus(Sprintf("ERROR: %o places of degree %o, expected %o", #PL, HIGH, PCOUNT[HIGH]));
        error "wrong number of places";
    end if;
    BALL := 0; BKEPT := 0;
    for Q in PL do
        b := PlaceBucket(Q);
        fp := BucketPrint(b);
        BALL := (BALL + fp) mod SL_FPMOD;              // fingerprint of ALL places' buckets
        if BucketHash(b, GSIZE) eq RANK then
            Append(~HIGHPL, Q);
            BKEPT := (BKEPT + fp) mod SL_FPMOD;        // fingerprint of the kept ones
        end if;
    end for;
    delete PL;
    printf "Degree %o: %o places, %o in this job's bucket class (%o s)\n",
           HIGH, PCOUNT[HIGH], #HIGHPL, Round(Cputime(t1));
end if;

// ---------------- the search ----------------
outf := Open(OUTFILE, "w");
fprintf outf, "Functions of degree <= %o seen in job %o of %o (group %o, rank %o of %o)\n\n",
        LOGMAX, JOB, NJOBS, GROUP, RANK, GSIZE;
Flush(outf);

// returns smallest degree seen, and whether a function of degree <= Dtarget was found
TestDivisor := function(excess, lamstr)
    space, mp := RiemannRochSpace(cuspsum + excess);
    if Dimension(space) le 1 then return Infinity(), false; end if;
    best := Infinity(); logged := {};
    for vct in space do
        g := mp(vct);
        if g in kbar then continue; end if;
        dgr := Degree(g);
        if dgr lt best then best := dgr; end if;
        if dgr le LOGMAX then
            ssp, svl := Support(Divisor(g));
            poleD := &+([DG0] cat [ -svl[i]*ssp[i] : i in [1..#ssp] | svl[i] lt 0 ]);
            if poleD notin logged then
                Include(~logged, poleD);
                fprintf outf, "---- degree %o, partition %o ----\nexcess E = %o\nfunction = %o\npoles = %o\n\n",
                        dgr, lamstr, excess, g, poleD;
                Flush(outf);
            end if;
        end if;
        if dgr le Dtarget then
            ff := Open(FOUNDFILE, "a");
            fprintf ff, "job %o of %o: degree %o, partition %o\nexcess E = %o\nfunction = %o\n\n",
                    JOB, NJOBS, dgr, lamstr, excess, g;
            Flush(ff); delete ff;
            return dgr, true;
        end if;
    end for;
    return best, false;
end function;

best := Infinity(); found := false; stopped := false; tested := 0; t0 := Cputime();
NH := (HIGH gt 0) select #HIGHPL else 1;
EXPECT := NLOW * NH;
printf "Divisors to test in this job: %o\n", EXPECT;
for item in LOWS do
    Elow := DG0;
    for pr in item[2] do Elow +:= 1*PLACEOF[Sprintf("%o:%o", pr[1], pr[2])]; end for;
    for h in [1..NH] do
        if HIGH gt 0 then excess := Elow + 1*HIGHPL[h]; else excess := Elow; end if;
        b, f := TestDivisor(excess, item[1]);
        tested +:= 1;
        if b lt best then best := b; end if;
        if f then found := true; break; end if;
        if LIMIT gt 0 and tested ge LIMIT then stopped := true; break; end if;
        if tested mod REPORT eq 0 then
            printf "  [%os] tested %o of %o, smallest degree so far %o\n", Round(Cputime(t0)), tested, EXPECT, best;
            if FileExists(FOUNDFILE) then stopped := true; break; end if;
        end if;
    end for;
    if found or stopped then break; end if;
end for;

if found then
    WriteStatus(Sprintf("FOUND JOB %o GROUP %o degree %o (see FOUND.txt)", JOB, GROUP, best));
    printf "==> FOUND a function of degree %o <= %o.\n", best, Dtarget;
elif LIMIT gt 0 and stopped then
    WriteStatus(Sprintf("PARTIAL JOB %o GROUP %o test mode: tested %o divisors, smallest degree %o", JOB, GROUP, tested, best));
    printf "==> Test mode: %o divisors tested in %o s (smallest degree %o).\n", tested, Round(Cputime(t0)), best;
elif stopped then
    WriteStatus(Sprintf("STOPPED JOB %o GROUP %o another job found a function; tested %o", JOB, GROUP, tested));
    print "==> Stopped: another job found a function (FOUND.txt).";
else
    WriteStatus(Sprintf("COMPLETED JOB %o GROUP %o RANK %o GSIZE %o HIGHDEG %o HIGHKEPT %o NLOW %o TESTED %o BALL %o BKEPT %o BEST %o",
                        JOB, GROUP, RANK, GSIZE, HIGH, #HIGHPL, NLOW, tested, BALL, BKEPT, best));
    printf "==> COMPLETED: %o divisors, no function of degree <= %o (smallest seen %o), %o s.\n",
           tested, Dtarget, best, Round(Cputime(t0));
end if;
delete outf;
quit;

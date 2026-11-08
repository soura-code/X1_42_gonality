#!/usr/bin/env python3
"""
X1_42_split_check.py  --  verify a finished run of X1_42_split_search.m
Checks, from plan.txt and the status files in DIR (default X1_42_split_d11):
  * every job reports COMPLETED (or FOUND / STOPPED, which is reported),
  * group H9..H6: every job of the group saw the same buckets (same fingerprint
    BALL of all P_v places), and the kept places of the jobs add up to all of
    them (sum of HIGHKEPT = P_v and sum of BKEPT = BALL mod 2^61-1),
    and sum of TESTED = NLOW * P_v,
  * group LOW: sum of TESTED = number of kept low divisors,
  * total tested = TOTAL of plan.txt.
Only if all checks pass does it conclude gon_F5(X1(42)) >= 12.
Usage:  python3 X1_42_split_check.py [DIR]
"""
import os
import sys

d = sys.argv[1] if len(sys.argv) > 1 else "X1_42_split_d11"
plan = open(os.path.join(d, "plan.txt")).read().split("\n")
N, TOTAL = int(plan[0].split()[1]), int(plan[0].split()[3])
groups, jobs = {}, {}
for l in plan[1:]:
    t = l.split()
    if not t:
        continue
    kv = dict(zip(t[2::2], t[3::2]))
    if t[0] == "GROUP":
        groups[t[1]] = {k: int(v) for k, v in kv.items()}
    elif t[0] == "JOB":
        jobs[int(t[1])] = kv
ok, found = True, False
M = 2305843009213693951
done = {g: {"kept": 0, "tested": 0, "n": 0, "fp": 0, "ball": set()} for g in groups}
for k in range(N):
    name = os.path.join(d, "status_job%d_of_%d.txt" % (k, N))
    if not os.path.exists(name):
        print("job %d: no status file (not finished?)" % k); ok = False; continue
    s = open(name).read().strip()
    t = s.split()
    if t[0] == "FOUND":
        found = True; print("job %d: %s" % (k, s)); continue
    if t[0] != "COMPLETED":
        print("job %d: %s" % (k, s)); ok = False; continue
    kv = dict(zip(t[1::2], t[2::2]))
    g = kv["GROUP"]
    if int(kv["JOB"]) != k or g != jobs[k]["GROUP"]:
        print("job %d: status does not match the plan" % k); ok = False; continue
    done[g]["kept"] += int(kv["HIGHKEPT"]); done[g]["tested"] += int(kv["TESTED"]); done[g]["n"] += 1
    done[g]["fp"] = (done[g]["fp"] + int(kv["BKEPT"])) % M; done[g]["ball"].add(int(kv["BALL"]))
    if int(kv["TESTED"]) != int(kv["NLOW"]) * (int(kv["HIGHKEPT"]) if int(kv["HIGHDEG"]) else 1):
        print("job %d: tested %s, expected NLOW*HIGHKEPT" % (k, kv["TESTED"])); ok = False
tot = 0
for g, G in groups.items():
    D = done[g]
    if D["n"] != G["JOBS"]:
        print("group %s: %d of %d jobs completed" % (g, D["n"], G["JOBS"])); ok = False; continue
    if G["HIGHDEG"]:
        if len(D["ball"]) != 1:
            print("group %s: the jobs computed DIFFERENT bucket data for degree %d" % (g, G["HIGHDEG"])); ok = False
        elif D["fp"] != next(iter(D["ball"])):
            print("group %s: the kept places do not add up to all places (fingerprint)" % g); ok = False
        if D["kept"] != G["PHIGH"]:
            print("group %s: the bucket classes hold %d places of degree %d, expected %d"
                  % (g, D["kept"], G["HIGHDEG"], G["PHIGH"])); ok = False
        else:
            print("group %s: all %d places of degree %d covered exactly once" % (g, G["PHIGH"], G["HIGHDEG"]))
    if D["tested"] != G["UNITS"]:
        print("group %s: tested %d divisors, expected %d" % (g, D["tested"], G["UNITS"])); ok = False
    else:
        print("group %s: %d divisors tested (as planned)" % (g, D["tested"]))
    tot += D["tested"]
if found:
    print("==> A function of degree <= 11 was FOUND (see %s/FOUND.txt): gon_F5(X1(42)) <= 11." % d)
elif ok and tot == TOTAL:
    print("==> All %d jobs COMPLETED, %d divisors tested = plan: no function of degree <= 11,"
          " hence gon_F5(X1(42)) >= 12." % (N, tot))
else:
    print("==> NOT complete or inconsistent: no conclusion.")
    sys.exit(1)

#!/usr/bin/env python3
"""Verify the survivors of a cn2x/cn2xc run directory and compare against known CN2s.

    python cn2x/verify_run.py cn2x/runs/capped_1e25 --lo 1e22 --hi 1e25 --tcap 20000

Must be run as a file, not piped to python: the worker pool re-imports this module.
"""
import argparse, json, math, os, sys, time
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cn2_exhaustive_c import _verify_survivor, icbrt, is_rigid_cn2
from sympy import factorint

KNOWN_22_25 = [1159954316194989017102401, 2196407820059694924883201, 3339611018825185787482801,
               4105879060352839139462401, 5002862939121639632040001, 7865064643837556041286401,
               9400084864021826054720641]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('rundir')
    ap.add_argument('--lo', type=float, default=0)
    ap.add_argument('--hi', type=float, required=True)
    ap.add_argument('--tcap', type=int, default=0, help="prime cap used by the run (0 = the proved bound)")
    ap.add_argument('--procs', type=int, default=28)
    a = ap.parse_args()
    lo, hi = int(a.lo), int(a.hi)
    T = a.tcap or icbrt(hi - 1)
    hits, sur = set(), []
    with open(os.path.join(a.rundir, 'out.txt')) as fh:
        for line in fh:
            p = line.split()
            if len(p) == 2 and p[0] == 'H':
                hits.add(int(p[1]))
            elif len(p) == 5 and p[0] == 'S':
                sur.append((int(p[1]), int(p[2]), int(p[3]), int(p[4]), T))
    print(f"{len(sur):,} survivors to verify with {a.procs} processes", flush=True)
    t0 = time.time()
    with Pool(a.procs) as pool:
        for n in pool.imap_unordered(_verify_survivor, sur, chunksize=512):
            if n is not None:
                hits.add(n)
    found = sorted(n for n in hits if lo <= n < hi and is_rigid_cn2(n))
    print(f"verified in {time.time() - t0:.0f}s\n")
    print(f"{'n':>26} {'log10':>7} {'k':>3} {'pmax':>6} {'geom mean':>9}  known?")
    for n in found:
        ps = sorted(factorint(n))
        gm = math.exp(sum(math.log(p) for p in ps) / len(ps))
        print(f"{n:>26} {math.log10(n):>7.3f} {len(ps):>3} {ps[-1]:>6} {gm:>9.0f}  "
              f"{'yes' if n in KNOWN_22_25 else '*** NEW ***'}")
    extras = [n for n in found if n not in KNOWN_22_25]
    missing = [n for n in KNOWN_22_25 if lo <= n < hi and n not in found]
    print(f"\nfound {len(found)} in [{lo:.0e}, {hi:.0e});  known present: {len(KNOWN_22_25) - len(missing)}/"
          f"{len([n for n in KNOWN_22_25 if lo <= n < hi])};  missing: {missing};  new: {extras}")
    json.dump(dict(rundir=a.rundir, lo=str(lo), hi=str(hi), tcap=a.tcap, found=[str(n) for n in found],
                   missing=[str(n) for n in missing], new=[str(n) for n in extras]),
              open(os.path.join(a.rundir, 'result.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()

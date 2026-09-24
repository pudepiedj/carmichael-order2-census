#!/usr/bin/env python3
r"""
Exhaustive enumeration of rigid CN2s below X -- a cost-measuring harness.

A rigid CN2 is a squarefree composite n with p^2 - 1 | n - 1 for every prime
p | n (the definition verified throughout CN2_Desert_Census.md).  Every CN2 is
a Carmichael number, so this is the CN2-only analogue of a Pinch/Goutier table.
The point of the harness is not to reach 10^24 but to measure how the cost of a
PROVABLY COMPLETE search grows with X, validated against the six CN2s below
10^22 (complete, by Goutier's table), and extrapolate.

Structure used (all proved in the note; checked on the six knowns):

  L1  For p | n:  n/p == p (mod p^2 - 1).   [n == 1 and p^2 == 1 mod p^2-1]
      Since n/p != p (squarefree), n/p >= p^2 + p - 1, so  p^3 < n <= X.
      Every prime factor is below X^(1/3).
  L2  2 and 3 never divide a CN2: each divides p^2 - 1 for every other p.
  L3  With q < r the two largest primes and P = n/(qr):  Pq == r (mod r^2-1)
      forces Pq > r^2, hence r < P.  So q < r < P, and n has >= 4 primes.

Search: depth-first over prefixes m = p_1 ... p_j in increasing prime order,
carrying L = lcm(p_i^2 - 1).  A prefix is pruned unless gcd(m, L) = 1 (a prime
dividing both would divide n and n - 1).  A candidate next prime t must avoid
L and satisfy gcd(m, t^2 - 1) = 1, and must fit at least one ROLE:

  last            prefix m == t (mod t^2-1)          =>  t^2 < m        (L1)
  second-to-last  t < r < P = m, and m t r <= X      =>  t < m,  m t^2 < X   (L3)
  earlier         two more primes above t            =>  m t^3 < X

At each prefix the completions s = n/m satisfy s == m^{-1} (mod L), s <= X/m.
When that progression is short compared with the number of admissible t, the
node switches to LIST mode: every s in the progression is factored and checked,
and the subtree is not expanded.  Either way every CN2 <= X is reached along
its own prefix path, so the search is complete for any choice of switch rule;
the rule only moves cost between branching and factoring.

    python cn2_exhaustive.py 1e16                   # one bound, validated
    python cn2_exhaustive.py --scan 12 19 --step 0.5 --json cn2_exhaustive_scan.jsonl
"""

import argparse
import bisect
import json
import math
import sys
import time
from typing import List

import numpy as np
from sympy import factorint

# Complete below 10^22 (Goutier's Carmichael table; CN2_Desert_Census.md).
KNOWN_BELOW_1E22 = [
    443372888629441,
    39671149333495681,
    842526563598720001,
    2380296518909971201,
    3188618003602886401,
    4208895375600667752001,
]
GROUND_TRUTH_LIMIT = 10**22


def icbrt(n: int) -> int:
    """Largest integer c with c^3 <= n."""
    if n <= 0:
        return 0
    c = int(round(n ** (1.0 / 3.0)))
    while c**3 > n:
        c -= 1
    while (c + 1)**3 <= n:
        c += 1
    return c


def primes_from_5(limit: int) -> List[int]:
    if limit < 5:
        return []
    sieve = np.ones(limit + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, math.isqrt(limit) + 1):
        if sieve[p]:
            sieve[p * p::p] = False
    return [int(p) for p in np.nonzero(sieve)[0] if p >= 5]


def is_rigid_cn2(n: int) -> bool:
    f = factorint(n)
    return (len(f) >= 3 and all(e == 1 for e in f.values())
            and all((n - 1) % (p * p - 1) == 0 for p in f))


class Enumerator:
    def __init__(self, X: int, list_ratio: float = 0.02, time_limit: float = 0.0,
                 fermat: str = 'on'):
        self.X = X
        self.T = icbrt(X - 1)                      # every prime p has p^3 < n <= X
        self.primes = primes_from_5(self.T)
        self.list_ratio = list_ratio
        self.time_limit = time_limit
        self.hits = set()
        self.nodes = self.children = self.last_tests = 0
        self.list_nodes = self.list_candidates = self.factorisations = 0
        self.fermat = fermat                       # 'off' | 'measure' | 'on'
        self.fermat_passes = 0
        self.pruned_budget = 0
        self.aborted = False
        self._t0 = 0.0

    def run(self):
        self._t0 = time.time()
        sys.setrecursionlimit(10000)
        self._dfs(1, 3, 1, 0)
        self.elapsed = time.time() - self._t0
        return self

    def _record(self, n: int):
        if n <= self.X and is_rigid_cn2(n):
            self.hits.add(n)

    def _dfs(self, m: int, pmax: int, L: int, k: int):
        if self.aborted:
            return
        self.nodes += 1
        if self.time_limit and (self.nodes & 0x3FFF) == 0 \
                and time.time() - self._t0 > self.time_limit:
            self.aborted = True
            return

        if k >= 3 and (m - 1) % L == 0:
            self._record(m)

        X = self.X
        rem = (X - 1) // m                           # completions s satisfy m*s <= X, strict bounds use X-1
        if rem <= pmax:
            return

        # role bounds on the next prime t (k-aware: n needs >= 3 primes in total)
        b_last = math.isqrt(m - 1) if k >= 2 else 0
        b_second = min(m - 1, math.isqrt(rem)) if k >= 1 else 0
        b_early = icbrt(rem)
        tmax = min(self.T, max(b_last, b_second, b_early))
        lo = bisect.bisect_right(self.primes, pmax)
        hi = bisect.bisect_right(self.primes, tmax)
        n_children = hi - lo
        if n_children <= 0:                          # no prime fits any role: no completion exists
            return

        # LIST mode: short progression of completions s == m^{-1} (mod L)
        if L > 1:
            n_list = rem // L + 1
            if n_list <= self.list_ratio * n_children:
                self._list(m, pmax, L, k, rem)
                return

        extend_bound = max(b_second, b_early)        # beyond this, t can only be last
        primes = self.primes
        for i in range(lo, hi):
            t = primes[i]
            self.children += 1
            if L % t == 0:
                continue
            t2 = t * t - 1
            if math.gcd(m, t2) != 1:
                continue
            mt = m * t
            if mt > X:
                break
            Lt = L * t2 // math.gcd(L, t2)
            if Lt >= X:                              # L(n) <= n-1 < X: provably dead
                self.pruned_budget += 1
                continue
            if t > extend_bound:                     # role 'last' only: test n = m t, no subtree
                self.last_tests += 1
                if k + 1 >= 3 and (mt - 1) % Lt == 0:
                    self._record(mt)
                continue
            self._dfs(mt, t, Lt, k + 1)

    def _list(self, m: int, pmax: int, L: int, k: int, rem: int):
        self.list_nodes += 1
        a = pow(m, -1, L)
        T = self.T
        for s in range(a, rem + 1, L):
            if s <= pmax:                            # includes s = 1 (n = m, already tested)
                continue
            self.list_candidates += 1
            if s % 2 == 0 or s % 3 == 0:
                continue
            if self.fermat != 'off':
                # n = m s is Carmichael, so 2^(n-1) == 1 mod every prime t | s,
                # hence mod s (squarefree, odd).  Never rejects a true completion.
                passes = pow(2, m * s - 1, s) == 1
                self.fermat_passes += passes
                if self.fermat == 'on' and not passes:
                    continue
            self.factorisations += 1
            f = factorint(s)
            if any(e > 1 for e in f.values()):
                continue
            ps = sorted(f)
            if ps[0] <= pmax or ps[-1] > T or k + len(ps) < 3:
                continue
            n = m * s
            if all((n - 1) % (p * p - 1) == 0 for p in ps):
                self._record(n)


def validate(X: int, hits) -> str:
    if X > GROUND_TRUTH_LIMIT:
        return "beyond ground truth"
    expected = sorted(n for n in KNOWN_BELOW_1E22 if n <= X)
    return "OK" if sorted(hits) == expected else f"MISMATCH expected {expected}"


def run_one(X: int, list_ratio: float, time_limit: float, fermat: str = 'on') -> dict:
    e = Enumerator(X, list_ratio, time_limit, fermat).run()
    row = dict(log10X=math.log10(X), T=e.T, primes=len(e.primes), nodes=e.nodes,
               children=e.children, last_tests=e.last_tests, list_nodes=e.list_nodes,
               list_candidates=e.list_candidates, factorisations=e.factorisations,
               fermat=e.fermat, fermat_passes=e.fermat_passes, pruned_budget=e.pruned_budget,
               hits=sorted(e.hits), n_hits=len(e.hits), seconds=e.elapsed,
               aborted=e.aborted, validation="ABORTED" if e.aborted else validate(X, e.hits))
    return row


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('X', nargs='?', type=float, help="single bound, e.g. 1e16")
    ap.add_argument('--scan', nargs=2, type=float, metavar=('LO', 'HI'),
                    help="log10 range of bounds to scan")
    ap.add_argument('--step', type=float, default=1.0)
    ap.add_argument('--list-ratio', type=float, default=0.02)
    ap.add_argument('--time-limit', type=float, default=0.0,
                    help="abort a run after this many seconds (0 = none)")
    ap.add_argument('--fermat', choices=['off', 'measure', 'on'], default='on',
                    help="base-2 Fermat pre-filter on LIST completions before factoring")
    ap.add_argument('--json', metavar='PATH')
    args = ap.parse_args(argv)

    if args.scan:
        lo, hi = args.scan
        bounds = [int(round(10 ** (lo + i * args.step)))
                  for i in range(int(round((hi - lo) / args.step)) + 1)]
    elif args.X:
        bounds = [int(args.X)]
    else:
        ap.error("give X or --scan")

    print(f"{'log10 X':>8} {'T':>9} {'nodes':>12} {'children':>13} {'list cand':>11} "
          f"{'factor':>9} {'hits':>5} {'seconds':>9}  validation")
    rows = []
    for X in bounds:
        row = run_one(X, args.list_ratio, args.time_limit, args.fermat)
        rows.append(row)
        print(f"{row['log10X']:>8.2f} {row['T']:>9,} {row['nodes']:>12,} {row['children']:>13,} "
              f"{row['list_candidates']:>11,} {row['factorisations']:>9,} {row['n_hits']:>5} "
              f"{row['seconds']:>9.2f}  {row['validation']}", flush=True)
        if args.json:
            with open(args.json, 'a') as fh:
                fh.write(json.dumps(row) + '\n')
        if row['aborted']:
            break

    done = [r for r in rows if not r['aborted'] and r['seconds'] > 0.5]
    if len(done) >= 3:
        xs = np.array([r['log10X'] for r in done[-4:]])
        ys = np.log10([r['seconds'] for r in done[-4:]])
        slope, icpt = np.polyfit(xs, ys, 1)
        print(f"\ncost ~ X^{slope:.3f} over log10 X in [{xs[0]:.1f}, {xs[-1]:.1f}] (single thread)")
        for target in (22, 23, 24):
            sec = 10 ** (icpt + slope * target)
            print(f"  extrapolated to 10^{target}: {sec:.3g} s = {sec/86400:.3g} core-days "
                  f"= {sec/86400/365.25:.3g} core-years")


if __name__ == '__main__':
    main()

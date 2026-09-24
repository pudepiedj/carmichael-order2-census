# Carmichael numbers of order 2: a complete census below 10^25

This repository contains the programs, the run record and the validation suite behind this result:

> **There are exactly 14 Carmichael numbers of order 2 below 10^25.**

"Order 2" is Howe's definition, the one used by OEIS [A175531](https://oeis.org/A175531): an odd composite $n$ such that, for every prime $p \mid n$,
$$n \equiv 1 \ \text{ or } \ n \equiv p \pmod{p^2-1}.$$

The terms (also in `results/b175531.txt`):

| n | a(n) | log10 | prime factors | largest |
|---|---|---|---|---|
| 1 | 443372888629441 | 14.65 | 8 | 331 |
| 2 | 39671149333495681 | 16.60 | 9 | 191 |
| 3 | 842526563598720001 | 17.93 | 8 | 2729 |
| 4 | 2380296518909971201 | 18.38 | 9 | 991 |
| 5 | 3188618003602886401 | 18.50 | 8 | 3457 |
| 6 | 4208895375600667752001 | 21.62 | 10 | 5851 |
| 7 | 1159954316194989017102401 | 24.06 | 9 | 4759 |
| 8 | 2088144166339753513992001 | 24.32 | 8 | 17291 |
| 9 | 2196407820059694924883201 | 24.34 | 9 | 13441 |
| 10 | 3339611018825185787482801 | 24.52 | 12 | 1429 |
| 11 | 4105879060352839139462401 | 24.61 | 8 | 7237 |
| 12 | 5002862939121639632040001 | 24.70 | 11 | 2549 |
| 13 | 7865064643837556041286401 | 24.90 | 11 | 15809 |
| 14 | 9400084864021826054720641 | 24.97 | 11 | 4523 |

**Factorisations.** Each term is the product of the primes listed, and $p^2-1$ divides $n-1$ for every prime $p \mid n$ is the second-order Korselt criterion on square-free $n$. `results/factorisations.txt` (readable) and `results/factorisations.json` (machine-readable) give the same data together with each term's modulus $L = \mathrm{lcm}_{p \mid n}(p^2-1)$ (which divides $n-1$), factored.

```
a(1)   = 443372888629441            = 17 · 31 · 41 · 43 · 89 · 97 · 167 · 331
a(2)   = 39671149333495681          = 17 · 37 · 41 · 71 · 79 · 97 · 113 · 131 · 191
a(3)   = 842526563598720001         = 17 · 61 · 71 · 89 · 197 · 311 · 769 · 2729
a(4)   = 2380296518909971201        = 19 · 41 · 43 · 71 · 89 · 127 · 199 · 449 · 991
a(5)   = 3188618003602886401        = 29 · 37 · 79 · 181 · 191 · 449 · 701 · 3457
a(6)   = 4208895375600667752001     = 17 · 29 · 31 · 43 · 71 · 79 · 199 · 389 · 2521 · 5851
a(7)   = 1159954316194989017102401  = 31 · 47 · 109 · 137 · 307 · 2927 · 3079 · 4049 · 4759
a(8)   = 2088144166339753513992001  = 113 · 199 · 239 · 263 · 701 · 7919 · 15391 · 17291
a(9)   = 2196407820059694924883201  = 53 · 109 · 127 · 281 · 307 · 379 · 647 · 10529 · 13441
a(10)  = 3339611018825185787482801  = 19 · 41 · 43 · 53 · 67 · 89 · 103 · 131 · 137 · 307 · 389 · 1429
a(11)  = 4105879060352839139462401  = 89 · 239 · 401 · 2143 · 2311 · 3571 · 3761 · 7237
a(12)  = 5002862939121639632040001  = 31 · 53 · 79 · 89 · 101 · 151 · 181 · 251 · 379 · 647 · 2549
a(13)  = 7865064643837556041286401  = 23 · 37 · 67 · 89 · 101 · 109 · 181 · 199 · 433 · 571 · 15809
a(14)  = 9400084864021826054720641  = 23 · 31 · 53 · 79 · 103 · 197 · 239 · 379 · 521 · 727 · 4523
```

The census also shows:

- **Terms 7–14 are all the terms in $[10^{22}, 10^{25})$.** Decades 22 and 23 contain none.
- **All 14 are "rigid":** $n \equiv 1 \pmod{p^2-1}$ for every $p \mid n$. No term below $10^{25}$ uses the residue $p$. Such terms do exist: Howe (2000, §5) gives one of about $3.92\times10^{59}$.
- **No term below $10^{25}$ has a prime factor above 17,291.**
- **The search uses no table of Carmichael numbers.** Terms 1–6 agree with Goutier's table of Carmichael numbers below $10^{22}$.

## How the search works

The full account, with proofs of the lemmas, is `docs/CN2_Exhaustive_Note.md` (§1–§2 for the depth-first search, §9 for the split and the census). In brief:

1. **Every prime factor is below $X^{1/3}$.** For $p \mid n$, the cofactor $m = n/p$ satisfies $m \equiv p$ or $m \equiv 1 \pmod{p^2-1}$, and $m$ is neither $p$ nor 1. So
   $$n = p\bigl(p + a(p^2-1)\bigr) \quad\text{or}\quad n = p\bigl(1 + a(p^2-1)\bigr), \qquad a \ge 1,$$
   and in particular $n > p^3$.
2. **The search splits at a prime bound $R_0$.** Every term below $X$ is covered by one of two phases:
   - **capped** (`cn2x/cn2xh.c`): a depth-first search over prime factors in increasing order, carrying the congruence $n \equiv c \pmod{L}$ that the chosen primes impose. Each prime may take either residue. It is restricted to primes $\le R_0$, which covers every term whose primes are all $\le R_0$.
   - **tail** (`cn2x/cn2tail.c`): a direct test of every pair $(p, a)$ in both forms above with $R_0 < p \le X^{1/3}$. That covers every term with some prime factor $> R_0$, about $X/(R_0^2 \ln R_0)$ pairs.
3. **Only necessary conditions are applied before factoring.** Each candidate passes squarefreeness, the definition for small primes, and a base-2 Fermat test (every term is a Carmichael number). The 1,812,744 survivors are then factored and checked exactly.

## Reproducing the census

Requirements: a C compiler with `unsigned __int128` (clang or gcc), Python 3.9+ and `sympy` (`pip install -r requirements.txt`). Multi-core is strongly advised.

```bash
python3 cn2_split_campaign.py init   --X 1e25 --R0 1e6 --dir cn2x/runs/howe_1e25   # compiles and freezes both engines
python3 cn2_split_campaign.py run    --dir cn2x/runs/howe_1e25 --now --threads 30
python3 cn2_split_campaign.py status --dir cn2x/runs/howe_1e25
python3 cn2_split_campaign.py verify --dir cn2x/runs/howe_1e25
```

The published run took 15.1 hours of engine time on 30 threads of an Apple M3 Ultra: 1.42 h for the tail and 13.71 h for the capped phase (Apple clang 21, macOS 26). Runs are checkpointed and can be interrupted and resumed:

- `drain` finishes the work in progress, then stops;
- `stop` stops at once, and the abandoned work is redone;
- `run --window 22:00-08:00` runs overnight only.

Resuming gives results identical to an uninterrupted run.

**Checking the published result without re-running the search.** The run record, including all 1.81 million candidates, is in `results/howe_1e25/`:

```bash
xz -dk results/howe_1e25/capped/out.txt.xz
python3 cn2_split_campaign.py verify --dir results/howe_1e25     # about 20 s on many cores; rewrites result.json
```

The engine sources are byte-identical to the copies frozen in the census run directory:

| file | SHA-256 |
|---|---|
| `cn2x/cn2xh.c` | `6c2c4273de9240c1089bb41a5aeee17f9e3a5661af5fd81438ef5c50e473899a` |
| `cn2x/cn2tail.c` | `aa29233ab979a0a92cb941b4eb7850baa4886db7330995873c2414f2ef292b20` |

The binary hashes in `results/howe_1e25/split.json` belong to that machine's build and will differ elsewhere.

## Validation

```bash
tests/run_validation.sh          # quick suite: about 7 minutes on 30 threads
tests/run_validation.sh --full   # adds the heavier checks (about 45 minutes more on 30 threads; the brute force needs 40 GB RAM)
```

| test | what it establishes |
|---|---|
| `parity` | in rigid mode, `cn2xh` is identical to the original engine `cn2xc`: all nine counters and every output record |
| `cheb` | the same engine, with modulus $(p^2-1)/2$, finds exactly the terms of OEIS A299799 below $10^{18}$, one of them not a rigid CN2 |
| `relaxed` | with modulus $(p^2-1)/D$ for $D = 24, 12, 8, 6$ and the Fermat filter off, the two phases together find every solution that brute-force factorisation of all $n < N$ finds, including solutions with a prime of residue type $p$, some reachable only through the tail's family for primes of type $p$. This is the test of the non-rigid code path: no list of true order-2 numbers small enough to brute-force contains a prime of type $p$ |
| `howe20` | the whole pipeline at $X = 10^{20}$: the five known terms, none non-rigid |
| `resume` | both engines, stopped and resumed mid-run, reproduce an uninterrupted run exactly |
| `published` | re-verifies the published candidates: 14 terms, none non-rigid |
| `tailhit` (full) | a term with a prime above $R_0$ is found through the tail ($X = 10^{18}$, $R_0 = 1000$) |
| `howe22` (full) | the whole pipeline at $X = 10^{22}$: Goutier's six terms, none non-rigid |

## Layout

| path | contents |
|---|---|
| `cn2_split_campaign.py` | the driver: `init`, `estimate`, `run`, `status`, `drain`, `stop`, `verify` |
| `cn2_order2_defs.py` | the three definitions (rigid, Howe, A299799) and the candidate verifier |
| `cn2_build.py` | portable compile and launch helpers |
| `cn2x/cn2xh.c`, `cn2x/cn2tail.c` | the two census engines |
| `cn2x/cn2xc.c`, `cn2x/cn2x.c`, `cn2_exhaustive.py`, `cn2_exhaustive_c.py`, `cn2x_campaign.py`, `cn2x_python_reference.json`, `cn2x/test_checkpoint.py`, `cn2x/verify_run.py` | the original rigid-only engines, their Python reference implementation and tests (`docs/` §1–§8). `cn2xc` is the parity reference for `cn2xh` |
| `tests/` | the validation suite and the brute-force reference `bf_order2.c` |
| `results/` | `b175531.txt`, `factorisations.txt` and `.json`, and the census run record `howe_1e25/` |
| `docs/CN2_Exhaustive_Note.md` | the method, validation and results. It mentions other documents of the wider research programme that are not part of this repository |

## References

- E. W. Howe, *Higher-order Carmichael numbers*, Math. Comp. 69 (2000), 1711–1719; [arXiv:math/9812089](https://arxiv.org/abs/math/9812089).
- OEIS [A175531](https://oeis.org/A175531) (Carmichael numbers of order 2) and [A299799](https://oeis.org/A299799) (Carmichael numbers that are Chebyshev pseudoprimes).
- C. Goutier, compressed table of the Carmichael numbers below $10^{22}$ (linked from A175531).

## Author and licence

John Puddefoot. The research direction and all decisions about what to compute and test are the author's; the programs and notes were written with the assistance of AI language models (mostly Anthropic's Claude Fable 5.1, Opus 5 & 5.5 with some additional help from GLM-5.3-MLX).

Licence: MIT (see `LICENSE`).

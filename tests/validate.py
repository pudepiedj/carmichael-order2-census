#!/usr/bin/env python3
"""Validation suite for the split order-2 Carmichael census (CN2_Exhaustive_Note.md §9.3).

    python3 tests/validate.py            # quick suite, a few minutes on a many-core machine
    python3 tests/validate.py --full     # adds the heavier checks (about 45 minutes more on 30 threads; 40 GB RAM for the brute force)
    python3 tests/validate.py --only parity cheb       # run named tests

Tests (quick unless marked full):
  parity     cn2xh in rigid mode against cn2xc: all nine counters and every output record identical
  cheb       cn2xh in A299799 mode to 1e18: exactly A299799's four terms there, one of them not a rigid CN2
  relaxed    modulus (p^2-1)/D, D = 24, 12, 8, 6, Fermat off: the split (capped at R0 = 30 plus tail) must find
             every solution a brute-force factorisation of all n < N finds, including those with type-p primes.
             N = 1e9 (quick) or 1e10 (full)
  howe20     the whole driver, Howe mode, X = 1e20, R0 = 3e4: the five known terms, none non-rigid
  resume     Howe mode, both engines: stop mid-run, resume, compare with an uninterrupted run
  published  re-verify the published candidates of results/howe_1e25 (if present): 14 terms, none non-rigid
  tailhit    (full) split, rigid, X = 1e18, R0 = 1000: 842526563598720001 (largest prime 2729) found by the tail
  howe22     (full) the whole driver, Howe mode, X = 1e22, R0 = 1e5: Goutier's six, none non-rigid (about 20 min)
"""
import argparse
import json
import lzma
import os
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from cn2_build import compile_cmd                       # noqa: E402
from cn2_order2_defs import is_order2, is_rigid         # noqa: E402

THREADS = max(1, (os.cpu_count() or 2) - 2)
BIN = os.path.join(ROOT, 'tests', 'bin')
A299799_BELOW_1E18 = [443372888629441, 582920080863121, 39671149333495681, 842526563598720001]
KNOWN_BELOW_1E20 = [443372888629441, 39671149333495681, 842526563598720001, 2380296518909971201,
                    3188618003602886401]
GOUTIER_SIX = KNOWN_BELOW_1E20 + [4208895375600667752001]
COUNTERS = ('nodes', 'children', 'last_tests', 'list_nodes', 'list_candidates',
            'fermat_passes', 'factorisations', 'pruned_budget', 'sieved')


def binary(name, src, extra=()):
    os.makedirs(BIN, exist_ok=True)
    out = os.path.join(BIN, name)
    if not os.path.exists(out) or os.path.getmtime(out) < os.path.getmtime(src):
        subprocess.run(compile_cmd(src, out, extra), check=True)
    return out


def dfs(engine, X, rundir, env=None):
    """cn2xc / cn2xh: X RUNDIR threads deadline ratio sieve mshard batch donate_q progress expected."""
    e = dict(os.environ); e.update(env or {})
    p = subprocess.run([engine, str(X), rundir, str(THREADS), '0', '1.0', '1', '100000', '4096', '10000000000', '0', '0'],
                       capture_output=True, text=True, env=e)
    return json.loads(p.stdout.strip().splitlines()[-1])


def tail(X, R0, rundir, env=None):
    e = dict(os.environ); e.update(env or {})
    p = subprocess.run([binary('cn2tail', os.path.join(ROOT, 'cn2x', 'cn2tail.c')), str(X), str(R0), rundir,
                        str(THREADS), '0', '268435456', '0'], capture_output=True, text=True, env=e)
    return json.loads(p.stdout.strip().splitlines()[-1])


def candidates(rundir, D=1):
    """Every n a run directory proposes: path hits, DFS survivors m*s, tail pairs of both families."""
    out = set()
    for line in open(os.path.join(rundir, 'out.txt')):
        t = line.split()
        if len(t) == 2 and t[0] == 'H':
            out.add(int(t[1]))
        elif len(t) == 5 and t[0] == 'S':
            out.add(int(t[1]) * int(t[4]))
        elif len(t) == 3 and t[0] in 'PQ':
            p, a = int(t[1]), int(t[2])
            g = (p * p - 1) // D
            out.add(p * ((p if t[0] == 'P' else 1) + a * g))
    return out


def found(rundir, mode, X, D=1):
    return sorted(n for n in candidates(rundir, D) if n < X and is_order2(n, mode))


# ---------------------------------------------------------------- tests
def t_parity(tmp):
    X = 10 ** 16
    a, b = os.path.join(tmp, 'xc'), os.path.join(tmp, 'xh')
    os.makedirs(a); os.makedirs(b)
    ra = dfs(binary('cn2xc', os.path.join(ROOT, 'cn2x', 'cn2xc.c')), X, a)
    rb = dfs(binary('cn2xh', os.path.join(ROOT, 'cn2x', 'cn2xh.c')), X, b, {'CN2X_MODE': 'rigid'})
    same_c = all(ra[k] == rb[k] for k in COUNTERS)
    same_r = sorted(open(os.path.join(a, 'out.txt'))) == sorted(open(os.path.join(b, 'out.txt')))
    return same_c and same_r, f"counters {'identical' if same_c else 'DIFFER'}, records {'identical' if same_r else 'DIFFER'} (nodes {ra['nodes']:,})"


def t_cheb(tmp):
    X = 10 ** 18
    d = os.path.join(tmp, 'cheb'); os.makedirs(d)
    dfs(binary('cn2xh', os.path.join(ROOT, 'cn2x', 'cn2xh.c')), X, d, {'CN2X_MODE': 'cheb'})
    f = found(d, 'cheb', X)
    nonrigid = [n for n in f if not is_rigid(n)]
    return f == A299799_BELOW_1E18 and nonrigid == [582920080863121], f"found {len(f)} of A299799's 4 below 1e18; not rigid CN2: {nonrigid}"


def t_relaxed(tmp, N):
    bf = binary('bf_order2', os.path.join(ROOT, 'tests', 'bf_order2.c'))
    xh = binary('cn2xh', os.path.join(ROOT, 'cn2x', 'cn2xh.c'))
    R0, ok, notes = 30, True, []
    for D in (24, 12, 8, 6):
        ref = {}
        for line in subprocess.run([bf, str(N), str(D)], capture_output=True, text=True, check=True).stdout.split('\n'):
            if line.strip():
                n, h = line.split(); ref[int(n)] = int(h)
        env = {'CN2X_MODE': 'test', 'CN2X_D': str(D)}
        c, t = os.path.join(tmp, f'c{D}'), os.path.join(tmp, f't{D}')
        os.makedirs(c); os.makedirs(t)
        dfs(xh, N, c, dict(env, CN2X_TCAP=str(R0)))
        tail(N, R0, t, env)
        got = candidates(c, D) | candidates(t, D)
        missed = sorted(n for n in ref if n not in got)
        viaQ = 0
        for line in open(os.path.join(t, 'out.txt')):
            s = line.split()
            if s[0] == 'Q':
                p, a = int(s[1]), int(s[2])
                if p * (1 + a * ((p * p - 1) // D)) in ref:
                    viaQ += 1
        ok &= not missed
        notes.append(f"D={D}: {len(ref)} solutions ({sum(ref.values())} with a type-p prime), missed {len(missed)}, via type-p family {viaQ}")
    return ok, '; '.join(notes)


def driver(args, check=True):
    return subprocess.run([sys.executable, os.path.join(ROOT, 'cn2_split_campaign.py')] + args,
                          capture_output=True, text=True, check=check)


def split_run(tmp, X, R0, mode):
    d = os.path.join(tmp, f'split_{mode}_{X}')
    driver(['init', '--X', str(X), '--R0', str(R0), '--dir', d, '--mode', mode])
    driver(['run', '--dir', d, '--now', '--threads', str(THREADS), '--progress', '0'])
    driver(['verify', '--dir', d, '--procs', str(THREADS)])
    return json.load(open(os.path.join(d, 'result.json'))), d


def t_howe(tmp, X, R0, expect):
    r, _ = split_run(tmp, X, R0, 'howe')
    f = [int(n) for n in r['found']]
    ok = r['complete'] and f == expect and not r['non_rigid'] and not r['inconsistencies']
    return ok, f"found {len(f)} (expected {len(expect)}), non-rigid {len(r['non_rigid'])}, inconsistencies {len(r['inconsistencies'])}"


def t_tailhit(tmp):
    r, d = split_run(tmp, 10 ** 18, 1000, 'rigid')
    target = 842526563598720001
    by_tail = target in {n for n in candidates(os.path.join(d, 'tail')) if is_rigid(n)}
    ok = r['complete'] and [int(n) for n in r['found']] == KNOWN_BELOW_1E20[:3] and by_tail
    return ok, f"found {r['n_found']}; 842526563598720001 found by the tail: {by_tail}"


def t_resume(tmp):
    notes, ok = [], True
    xh = binary('cn2xh', os.path.join(ROOT, 'cn2x', 'cn2xh.c'))
    X = 10 ** 17
    a, b = os.path.join(tmp, 'ra'), os.path.join(tmp, 'rb')
    os.makedirs(a); os.makedirs(b)
    env = dict(os.environ, CN2X_MODE='howe')
    full = dfs(xh, X, a, {'CN2X_MODE': 'howe'})
    cmd = [xh, str(X), b, str(THREADS), '0', '1.0', '1', '100000', '4096', '10000000000', '0', '0']
    p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env)
    time.sleep(1.5); open(os.path.join(b, 'STOP'), 'w').close(); p.wait()
    part = sum(1 for _ in open(os.path.join(b, 'done.log')))
    res = dfs(xh, X, b, {'CN2X_MODE': 'howe'})
    same = all(full[k] == res[k] for k in COUNTERS) and \
        sorted(set(open(os.path.join(a, 'out.txt')))) == sorted(set(open(os.path.join(b, 'out.txt'))))
    ok &= same and res['complete']
    notes.append(f"capped: stopped after {part} shards, resumed, {'identical' if same else 'DIFFERENT'}")
    X, R0 = 10 ** 17, 2000
    a, b = os.path.join(tmp, 'ta'), os.path.join(tmp, 'tb')
    os.makedirs(a); os.makedirs(b)
    full = tail(X, R0, a, {'CN2X_MODE': 'howe'})
    cmd = [binary('cn2tail', os.path.join(ROOT, 'cn2x', 'cn2tail.c')), str(X), str(R0), b, str(THREADS), '0', '268435456', '0']
    p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env)
    time.sleep(3); open(os.path.join(b, 'DRAIN'), 'w').close(); p.wait()
    part = sum(1 for _ in open(os.path.join(b, 'done.log')))
    res = tail(X, R0, b, {'CN2X_MODE': 'howe'})
    keys = ('total_pairs', 'tests', 'sieved', 'fermat_calls', 'fermat_passes', 'shards_done')
    same = all(full[k] == res[k] for k in keys) and \
        sorted(set(open(os.path.join(a, 'out.txt')))) == sorted(set(open(os.path.join(b, 'out.txt'))))
    ok &= same and res['complete']
    notes.append(f"tail: drained after {part} shards, resumed, {'identical' if same else 'DIFFERENT'} ({full['tests']:,} pairs)")
    return ok, '; '.join(notes)


def t_published(tmp):
    src = os.path.join(ROOT, 'results', 'howe_1e25')
    if not os.path.exists(os.path.join(src, 'split.json')):
        return None, "results/howe_1e25 not present (skipped)"
    d = os.path.join(tmp, 'howe_1e25')
    shutil.copytree(src, d)
    xz = os.path.join(d, 'capped', 'out.txt.xz')
    if os.path.exists(xz):
        with lzma.open(xz) as fi, open(os.path.join(d, 'capped', 'out.txt'), 'wb') as fo:
            shutil.copyfileobj(fi, fo)
    driver(['verify', '--dir', d, '--procs', str(THREADS)])
    r = json.load(open(os.path.join(d, 'result.json')))
    ok = r['complete'] and r['n_found'] == 14 and not r['non_rigid'] and not r['inconsistencies'] and not r['known_missing']
    return ok, f"{r['candidates']['capped']:,} + {r['candidates']['tail']:,} candidates re-verified: {r['n_found']} terms, non-rigid {len(r['non_rigid'])}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--full', action='store_true')
    ap.add_argument('--only', nargs='+')
    ap.add_argument('--keep', action='store_true', help="keep the temporary run directories")
    a = ap.parse_args()
    tests = [('parity', t_parity), ('cheb', t_cheb), ('relaxed', lambda t: t_relaxed(t, 10 ** 10 if a.full else 10 ** 9)),
             ('howe20', lambda t: t_howe(t, 10 ** 20, 30000, KNOWN_BELOW_1E20)), ('resume', t_resume),
             ('published', t_published)]
    if a.full:
        tests += [('tailhit', t_tailhit), ('howe22', lambda t: t_howe(t, 10 ** 22, 10 ** 5, GOUTIER_SIX))]
    if a.only:
        tests = [t for t in tests if t[0] in a.only]
    print(f"validating with {THREADS} threads\n")
    failed = 0
    for name, fn in tests:
        tmp = tempfile.mkdtemp(prefix=f'cn2val_{name}_')
        t0 = time.time()
        try:
            ok, msg = fn(tmp)
        except Exception as e:                      # a crash is a failure, reported with its cause
            ok, msg = False, f"error: {e!r}"
        status = 'SKIP' if ok is None else ('PASS' if ok else 'FAIL')
        failed += ok is False
        print(f"[{status}] {name:<10} {time.time() - t0:7.1f}s  {msg}", flush=True)
        if not a.keep:
            shutil.rmtree(tmp, ignore_errors=True)
    print(f"\n{'ALL PASSED' if not failed else f'{failed} FAILED'}")
    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()

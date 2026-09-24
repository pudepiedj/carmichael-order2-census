#!/usr/bin/env python3
r"""
Split exhaustive search for order-2 Carmichael numbers: capped DFS (cn2xh, CN2X_TCAP=R0) + large-prime tail (cn2tail).

Definition (--mode, fixed at init):
    howe   (default) n == 1 or p (mod p^2-1) for every p | n   Howe (2000), OEIS A175531 -- what OEIS calls
                     "Carmichael numbers of order 2"; each prime is type 1 (n == 1) or type p (n == p)
    rigid            n == 1 (mod p^2-1): the programme's CN2 (all primes type 1), a subset of howe
    cheb             n == 1 or p (mod (p^2-1)/2): Carmichael & Chebyshev, OEIS A299799 (a validation target)

Why it is complete.  For every prime p | n of a rigid CN2, n/p == p (mod p^2-1) and n/p > p (Lem1), so
    n = p (p + a (p^2-1)),   n - 1 = (p^2-1)(a p + 1),   a >= 1           (p of type 1)
and a prime of type p has n/p == 1 (mod p^2-1), n = p (1 + a (p^2-1)), a >= 1  (howe; modulus (p^2-1)/2 for cheb).
  * phase "capped": cn2xh (both types per prime, CRT residue per node) with its prime table cut at R0
    finds every solution < X whose primes are all <= R0 (the cap only removes table primes; role
    bounds and LIST completions are unchanged, so the full path of such a solution is in the tree).
  * phase "tail":   cn2tail tests every pair (p, a), both families, with R0 < p <= (D X)^(1/3), n < X --
    which covers every solution < X with at least one prime > R0.  About X / (R0^2 ln R0) pairs (howe).
The union is every solution below X, with no hypothesis.  Both phases use only necessary conditions
before factoring (the definition for small primes, squarefreeness, base-2 Fermat: all three
definitions imply Carmichael).  Validated: rigid mode bit-identical to cn2xc/old cn2tail; cheb mode
reproduces A299799 to 1e20; a relaxed test mode (modulus (p^2-1)/D, D | 24) matches brute force to 1e10
including solutions with type-p primes, found through both phases.

    python cn2_split_campaign.py estimate --X 1e25 --R0 3e5 5e5 1e6 --threads 30
    python cn2_split_campaign.py init     --X 1e25 --R0 1e6 --dir cn2x/runs/howe_1e25        # --mode howe is the default
    python cn2_split_campaign.py run      --dir cn2x/runs/howe_1e25 --now --threads 30
    python cn2_split_campaign.py run      --dir cn2x/runs/howe_1e25 --window 22:00-08:00 --threads 30
    python cn2_split_campaign.py status   --dir cn2x/runs/howe_1e25
    python cn2_split_campaign.py drain    --dir cn2x/runs/howe_1e25     # graceful end of the session
    python cn2_split_campaign.py verify   --dir cn2x/runs/howe_1e25     # after both phases complete

Phases run tail first (its cost is known exactly in advance), then capped.  `--phase` runs just one.
Each phase is independently checkpointed (done.log per phase directory); `init` freezes both engine
binaries and `run` refuses to resume against anything different, since shard ids depend on them.
To stop at once instead of draining: touch <dir>/tail/STOP or <dir>/capped/STOP (the `stop` command
touches both); in-flight shards are redone next session, results stay exact.
"""

import argparse
import datetime as dt
import json
import math
import os
import shutil
import subprocess
import sys
import time
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cn2x_campaign import parse_bound, sha256, next_time, in_window          # noqa: E402
from cn2_exhaustive_c import icbrt                                             # noqa: E402
from cn2_order2_defs import verify_candidate, is_order2, is_rigid, DD            # noqa: E402
from sympy import factorint                                                    # noqa: E402
from cn2_build import compile_cmd, keep_awake                                  # noqa: E402

SRC_CAPPED = os.path.join(HERE, 'cn2x', 'cn2xh.c')
SRC_TAIL = os.path.join(HERE, 'cn2x', 'cn2tail.c')
PHASES = ('tail', 'capped')

# every CN2 below 10^25 known on 23 Sep 2026: the six below 1e22 (complete, unconditional) and the
# eight of decade 24 (the unconditional 1e24 census proves decades 22-23 empty)
KNOWN_BELOW_1E25 = [443372888629441, 39671149333495681, 842526563598720001, 2380296518909971201,
                    3188618003602886401, 4208895375600667752001,
                    1159954316194989017102401, 2088144166339753513992001, 2196407820059694924883201,
                    3339611018825185787482801, 4105879060352839139462401, 5002862939121639632040001,
                    7865064643837556041286401, 9400084864021826054720641]
# A299799 (cheb), complete below 1e22 per OEIS
KNOWN_CHEB = [443372888629441, 582920080863121, 39671149333495681, 842526563598720001, 2380296518909971201,
              3188618003602886401, 33711266676317630401, 54764632857801026161, 122762671289519184001,
              361266866679292635601, 4208895375600667752001, 7673096805497432749441]
KNOWN = dict(rigid=KNOWN_BELOW_1E25, howe=KNOWN_BELOW_1E25, cheb=KNOWN_CHEB)
KNOWN_COMPLETE_TO = dict(rigid=10**22, howe=10**22, cheb=10**22)   # ground truth beyond which a find is news


def build(src, out):
    subprocess.run(compile_cmd(src, out), check=True)


def tail_pairs(X, R0):
    """Exact number of (p, a) pairs the tail tests: sum over primes R0 < p <= X^(1/3) of A(p)."""
    from sympy import primerange
    T = icbrt(X - 1)
    tot = 0
    for p in primerange(R0 + 1, T + 1):
        A = (X - 1 - p * p) // (p * (p * p - 1))
        if A <= 0:
            break
        tot += A
    return tot


def tail_pairs_approx(X, R0):
    return X / (2 * R0 ** 2 * math.log(R0))


def load(d, name):
    with open(os.path.join(d, name)) as fh:
        return json.load(fh)


# ---------------------------------------------------------------- init / estimate
def cmd_init(a):
    d = os.path.abspath(a.dir)
    if os.path.exists(os.path.join(d, 'split.json')):
        sys.exit(f"{d} already initialised; refusing to overwrite")
    X, R0 = parse_bound(a.X), int(parse_bound(a.R0))
    T = icbrt(DD[a.mode] * (X - 1))
    if not (1000 <= R0 < T):
        sys.exit(f"R0 must satisfy 1000 <= R0 < X^(1/3) = {T}")
    for ph, src, exe in (('capped', SRC_CAPPED, 'cn2xh'), ('tail', SRC_TAIL, 'cn2tail')):
        pd = os.path.join(d, ph)
        os.makedirs(pd, exist_ok=True)
        build(src, os.path.join(pd, exe))
        shutil.copy(src, os.path.join(pd, os.path.basename(src)))
    params = dict(mode=a.mode, X=str(X), R0=R0, T=T, mshard=a.mshard, batch=a.batch, ratio=1.0, sieve=1, donate_q=10**10,
                  W=a.W, nsmall=a.nsmall,
                  capped_sha256=sha256(os.path.join(d, 'capped', 'cn2xh')),
                  tail_sha256=sha256(os.path.join(d, 'tail', 'cn2tail')),
                  tail_pairs_approx=tail_pairs_approx(X, R0),
                  created=dt.datetime.now().isoformat(timespec='seconds'))
    with open(os.path.join(d, 'split.json'), 'w') as fh:
        json.dump(params, fh, indent=1)
    with open(os.path.join(d, 'capped', 'ENV.txt'), 'w') as fh:
        fh.write(f"Capped phase of a SPLIT run, CN2X_MODE={a.mode}: CN2X_TCAP={R0}.  Complete only together with the tail phase\n"
                 f"(cn2tail, every (p, a) with p > {R0}).  X={X}.  CN2X_NMIN unset (everything below X reported).\n")
    print(f"initialised {d}\n  mode = {a.mode}   X = {X:.3e}   R0 = {R0:,}   T = X^(1/3) = {T:,}\n"
          f"  tail: ~{params['tail_pairs_approx']:.3e} (p, a) pairs\n"
          f"  capped: cn2xh with CN2X_TCAP={R0}, mshard {a.mshard}, batch {a.batch}\n"
          f"  engines frozen: capped {params['capped_sha256'][:12]}..., tail {params['tail_sha256'][:12]}...")


def cmd_estimate(a):
    X = parse_bound(a.X)
    exe = os.path.join(HERE, 'cn2x', 'cn2tail')
    if not os.path.exists(exe) or os.path.getmtime(exe) < os.path.getmtime(SRC_TAIL):
        build(SRC_TAIL, exe)
    print(f"X = {X:.3e}; tail cost only (the capped phase must be measured -- see the note in the docstring)")
    print(f"{'R0':>10} {'pairs':>11} {'ns/pair':>8} {'tail h @%d thr' % a.threads:>15}")
    for r in a.R0:
        R0 = int(parse_bound(r))
        b = json.loads(subprocess.run([exe, '--bench', str(X), str(R0), str(a.bench_seconds)], capture_output=True,
                                      text=True, check=True, env=dict(os.environ, CN2X_MODE=a.mode)).stdout)
        pairs = (tail_pairs(X, R0) if a.exact else tail_pairs_approx(X, R0)) * (1 if a.mode == 'rigid' else 2 * DD[a.mode])
        hours = pairs * b['ns_per_test'] * 1e-9 / a.threads / 3600
        print(f"{R0:>10,} {pairs:>11.3e} {b['ns_per_test']:>8.0f} {hours:>15.2f}")
    print("ns/pair is measured at small a (one-limb Montgomery); where m = n/p exceeds 2^63 a test is slower,\n"
          "so treat the hours as a lower bound when X/R0 > 9.2e18.")


# ---------------------------------------------------------------- run
def run_phase(d, phase, p, threads, deadline_epoch, progress):
    pd = os.path.join(d, phase)
    env = dict(os.environ)
    env.pop('CN2X_NMIN', None); env.pop('CN2X_SIEVEP', None); env.pop('CN2X_D', None)
    env['CN2X_MODE'] = p.get('mode', 'rigid')
    if phase == 'capped':
        exe = os.path.join(pd, 'cn2xh')
        if sha256(exe) != p['capped_sha256']:
            sys.exit("capped engine binary has changed; refusing to resume")
        env['CN2X_TCAP'] = str(p['R0'])
        cmd = keep_awake([exe, p['X'], pd, str(threads), f"{deadline_epoch:.0f}", str(p['ratio']),
               str(p['sieve']), str(p['mshard']), str(p['batch']), str(p['donate_q']), str(progress), '0'])
    else:
        exe = os.path.join(pd, 'cn2tail')
        if sha256(exe) != p['tail_sha256']:
            sys.exit("tail engine binary has changed; refusing to resume")
        env.pop('CN2X_TCAP', None)
        cmd = keep_awake([exe, p['X'], str(p['R0']), pd, str(threads), f"{deadline_epoch:.0f}",
                          str(p['W']), str(progress), str(p['nsmall'])])
    started = time.time()
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, text=True, env=env)
    lines = [l for l in proc.stdout.splitlines() if l.startswith('{')]
    if not lines:
        print(f"{phase} engine produced no summary (exit {proc.returncode})", flush=True)
        return None
    s = json.loads(lines[-1])
    s.update(phase=phase, started=dt.datetime.fromtimestamp(started).isoformat(timespec='seconds'),
             ended=dt.datetime.now().isoformat(timespec='seconds'), wall_seconds=time.time() - started,
             threads=threads)
    with open(os.path.join(pd, 'sessions.jsonl'), 'a') as fh:
        fh.write(json.dumps(s) + '\n')
    return s


def phase_complete(d, phase):
    f = os.path.join(d, phase, 'sessions.jsonl')
    if not os.path.exists(f):
        return False
    ss = [json.loads(l) for l in open(f) if l.strip()]
    return bool(ss) and ss[-1]['complete']


def cmd_run(a):
    d = os.path.abspath(a.dir)
    p = load(d, 'split.json')
    phases = PHASES if a.phase == 'both' else (a.phase,)
    now_first, skip_until = a.now, None
    while True:
        todo = [ph for ph in phases if not phase_complete(d, ph)]
        if not todo:
            print("all requested phases complete -- run `verify`", flush=True)
            return
        now = dt.datetime.now()
        if now_first:
            end = None
        elif a.until:
            end = next_time(a.until, now)
        else:
            start, stop = a.window.split('-')
            if skip_until is not None or not in_window(now, start, stop):
                wake = skip_until or next_time(start, now)
                print(f"{now:%Y-%m-%d %H:%M} {'manual drain' if skip_until else 'outside window ' + a.window}; "
                      f"sleeping until {wake:%Y-%m-%d %H:%M}", flush=True)
                skip_until = None
                while dt.datetime.now() < wake:
                    time.sleep(min(60, max(1, (wake - dt.datetime.now()).total_seconds())))
                continue
            end = next_time(stop, now)
        interrupted = False
        for ph in todo:
            print(f"{dt.datetime.now():%Y-%m-%d %H:%M} {ph} session start; "
                  + (f"drain at {end:%H:%M}" if end else "no deadline (`drain` to end gracefully, `stop` to end now)"),
                  flush=True)
            s = run_phase(d, ph, p, a.threads, end.timestamp() if end else 0, a.progress)
            if s is None:
                sys.exit(1)
            print(f"{dt.datetime.now():%Y-%m-%d %H:%M} {ph} session end: complete={s['complete']} "
                  f"drained={s['drained']} stopped={s['stopped']} search {s['search_seconds'] / 3600:.2f} h", flush=True)
            if not s['complete']:
                interrupted = True
                break
        if not interrupted:
            continue                                     # loop re-checks; exits when all complete
        if a.until or s['stopped'] or (now_first and not a.window):
            return
        if s['drained'] and (end is None or time.time() < end.timestamp() - 60):
            skip_until = next_time(a.window.split('-')[0], dt.datetime.now())
        now_first = False


def cmd_touch(a, name):
    d = os.path.abspath(a.dir)
    for ph in PHASES:
        open(os.path.join(d, ph, name), 'w').close()
    print(f"touched {name} in both phase directories")


def cmd_status(a):
    d = os.path.abspath(a.dir)
    p = load(d, 'split.json')
    print(f"split run  mode {p.get('mode', 'rigid')}   X = {int(p['X']):.3e}   R0 = {p['R0']:,}")
    for ph in PHASES:
        pd = os.path.join(d, ph)
        f = os.path.join(pd, 'sessions.jsonl')
        ss = [json.loads(l) for l in open(f) if l.strip()] if os.path.exists(f) else []
        done = sum(1 for l in open(os.path.join(pd, 'done.log')) if l.startswith('D ')) \
            if os.path.exists(os.path.join(pd, 'done.log')) else 0
        out = os.path.join(pd, 'out.txt')
        surv = sum(1 for l in open(out) if l[:1] in 'SPH') if os.path.exists(out) else 0
        hours = sum(s['search_seconds'] for s in ss) / 3600
        state = 'COMPLETE' if ss and ss[-1]['complete'] else ('not started' if not ss else 'in progress')
        extra = f" of {ss[-1]['shards']:,}" if ss and ph == 'tail' else ''
        print(f"  {ph:>6}: {state:<12} sessions {len(ss)}  engine {hours:.2f} h  shards done {done:,}{extra}  records {surv:,}")


# ---------------------------------------------------------------- verify
def _verify_tail(args):
    fam, p, a, X, mode = args
    g = (p * p - 1) // DD[mode]
    m = (p if fam == 'P' else 1) + a * g
    if p * m >= X:
        return None
    return verify_candidate((p, m, mode))


def cmd_verify(a):
    d = os.path.abspath(a.dir)
    p = load(d, 'split.json')
    X, mode = int(p['X']), p.get('mode', 'rigid')
    complete = {ph: phase_complete(d, ph) for ph in PHASES}
    if not all(complete.values()):
        print(f"WARNING: not all phases complete {complete}; the result below is NOT a census")
    found = {ph: set() for ph in PHASES}
    cap_c, tail_c = set(), set()
    for line in open(os.path.join(d, 'capped', 'out.txt')):
        t = line.split()
        if len(t) == 2 and t[0] == 'H':
            cap_c.add((int(t[1]), 1, mode))
        elif len(t) == 5 and t[0] == 'S':
            cap_c.add((int(t[1]), int(t[4]), mode))
    for line in open(os.path.join(d, 'tail', 'out.txt')):
        t = line.split()
        if len(t) == 3 and t[0] in 'PQ':
            tail_c.add((t[0], int(t[1]), int(t[2]), X, mode))
    print(f"mode {mode}: verifying {len(cap_c):,} capped and {len(tail_c):,} tail candidates with {a.procs} processes",
          flush=True)
    t0 = time.time()
    with Pool(a.procs) as pool:
        for n in pool.imap_unordered(verify_candidate, list(cap_c), chunksize=512):
            if n is not None and n < X:
                found['capped'].add(n)
        for n in pool.imap_unordered(_verify_tail, list(tail_c), chunksize=512):
            if n is not None:
                found['tail'].add(n)
    allf = sorted(found['capped'] | found['tail'])
    fac = {n: factorint(n) for n in allf}
    incons = []
    for n in allf:
        big = max(fac[n]) > p['R0']
        if big and n not in found['tail']:
            incons.append((str(n), 'has a prime > R0 but tail missed it'))
        if not big and n not in found['capped']:
            incons.append((str(n), 'all primes <= R0 but capped missed it'))
    known = [n for n in KNOWN[mode] if n < X]
    missing = [n for n in known if n not in allf]
    new = [n for n in allf if n not in KNOWN[mode]]
    nonrigid = [n for n in allf if not is_rigid(n, fac[n])]
    by_decade = {}
    for n in allf:
        by_decade.setdefault(int(math.log10(n)), []).append(str(n))
    rep = dict(mode=mode, X=str(X), R0=p['R0'], complete=all(complete.values()), phases_complete=complete,
               n_found=len(allf), found=[str(n) for n in allf], non_rigid=[str(n) for n in nonrigid],
               by_decade={str(k): v for k, v in sorted(by_decade.items())},
               found_by_capped=len(found['capped']), found_by_tail=len(found['tail']),
               known_missing=[str(n) for n in missing], new=[str(n) for n in new], inconsistencies=incons,
               candidates=dict(capped=len(cap_c), tail=len(tail_c)), verify_seconds=time.time() - t0)
    with open(os.path.join(d, 'result.json'), 'w') as fh:
        json.dump(rep, fh, indent=1)
    print(f"\n{len(allf)} order-2 numbers ({mode}) < {X:.3e}   (complete={rep['complete']})")
    print(f"{'n':>26} {'log10':>7} {'k':>3} {'pmax':>9}  {'rigid?':<9} found by        known?")
    for n in allf:
        ps = sorted(fac[n])
        typep = [q for q in ps if n % ((q * q - 1) // DD[mode]) != 1]
        who = '+'.join(ph for ph in ('capped', 'tail') if n in found[ph])
        print(f"{n:>26} {math.log10(n):>7.3f} {len(ps):>3} {ps[-1]:>9}  {'rigid' if not typep else 'type-p ' + str(typep):<9} "
              f"{who:<14}  {'yes' if n in KNOWN[mode] else '*** NEW ***'}")
    print(f"\nknown below X re-found: {len(known) - len(missing)} of {len(known)}   new: {len(new)}")
    if missing and rep['complete']:
        print(f"*** ALARM: the run is complete but these known terms were not found: {[str(n) for n in missing]}")
    elif missing:
        print(f"    ({len(missing)} not yet re-found -- expected while the capped phase is incomplete)")
    print(f"non-rigid (some prime of type p): {len(nonrigid)}   {[str(n) for n in nonrigid]}")
    print(f"phase consistency: {'OK' if not incons else incons}")
    for k, v in sorted(by_decade.items()):
        print(f"  decade 10^{k}: {len(v)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    i = sub.add_parser('init')
    i.add_argument('--X', required=True); i.add_argument('--R0', required=True); i.add_argument('--dir', required=True)
    i.add_argument('--mode', choices=('howe', 'rigid', 'cheb'), default='howe')
    i.add_argument('--mshard', type=int, default=100000); i.add_argument('--batch', type=int, default=4096)
    i.add_argument('--W', type=int, default=1 << 28, help="tail shard size in (p, a) pairs")
    i.add_argument('--nsmall', type=int, default=95, help="small primes in the tail's Korselt sieve (0-95)")
    e = sub.add_parser('estimate')
    e.add_argument('--X', required=True); e.add_argument('--R0', nargs='+', required=True)
    e.add_argument('--mode', choices=('howe', 'rigid', 'cheb'), default='howe')
    e.add_argument('--threads', type=int, default=30); e.add_argument('--bench-seconds', type=float, default=10)
    e.add_argument('--exact', action='store_true', help="count pairs exactly (slow above 1e24)")
    r = sub.add_parser('run')
    r.add_argument('--dir', required=True); r.add_argument('--threads', type=int, default=30)
    r.add_argument('--phase', choices=('both',) + PHASES, default='both')
    g = r.add_mutually_exclusive_group()
    g.add_argument('--window'); g.add_argument('--until')
    r.add_argument('--now', action='store_true')
    r.add_argument('--progress', type=float, default=60)
    for name in ('status', 'drain', 'stop'):
        sub.add_parser(name).add_argument('--dir', required=True)
    v = sub.add_parser('verify'); v.add_argument('--dir', required=True); v.add_argument('--procs', type=int, default=24)
    a = ap.parse_args()
    if a.cmd == 'run' and not (a.window or a.until or a.now):
        ap.error("run needs --window, --until or --now")
    if a.cmd == 'drain':
        return cmd_touch(a, 'DRAIN')
    if a.cmd == 'stop':
        return cmd_touch(a, 'STOP')
    dict(init=cmd_init, estimate=cmd_estimate, run=cmd_run, status=cmd_status, verify=cmd_verify)[a.cmd](a)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Checkpoint/restart tests for cn2xc: exact totals and hits under deadline drains, STOP-file hard stops and kill -9."""
import json, os, random, shutil, signal, subprocess, sys, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cn2_exhaustive_c import _verify_survivor, icbrt, is_rigid_cn2

HERE = os.path.dirname(os.path.abspath(__file__))
COUNTERS = ('nodes', 'children', 'last_tests', 'list_nodes', 'list_candidates',
            'fermat_passes', 'factorisations', 'pruned_budget', 'sieved')

def reference(X, threads):
    out = subprocess.run([os.path.join(HERE, 'cn2x'), str(X), '/dev/null', str(threads), '1.0', '1', str(10**10), '1000', '1'],
                         capture_output=True, text=True, check=True).stdout
    return json.loads(out.strip().splitlines()[-1])

def hits_from(rundir, X):
    hits = set(); T = icbrt(X - 1)
    for line in open(os.path.join(rundir, 'out.txt')):
        p = line.split()
        if len(p) == 2 and p[0] == 'H':
            hits.add(int(p[1]))
        elif len(p) == 5 and p[0] == 'S':
            n = _verify_survivor((int(p[1]), int(p[2]), int(p[3]), int(p[4]), T))
            if n: hits.add(n)
    return sorted(n for n in hits if n <= X and is_rigid_cn2(n))

def session(X, rundir, threads, deadline=0, kill_after=None, stopfile_after=None):
    cmd = [os.path.join(HERE, 'cn2xc'), str(X), rundir, str(threads), str(deadline)]
    if kill_after is None and stopfile_after is None:
        p = subprocess.run(cmd, capture_output=True, text=True)
        return json.loads(p.stdout.strip().splitlines()[-1])
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    time.sleep(kill_after if kill_after is not None else stopfile_after)
    if p.poll() is None:
        if kill_after is not None:
            p.send_signal(signal.SIGKILL); p.wait(); return None
        open(os.path.join(rundir, 'STOP'), 'w').close()
    out = p.stdout.read(); p.wait()
    return json.loads(out.strip().splitlines()[-1])

def check(label, X, final, ref, hits, ref_hits):
    diffs = {c: (ref[c], final[c]) for c in COUNTERS if ref[c] != final[c]}
    ok = final['complete'] and not diffs and hits == ref_hits
    print(f"  {label}: {'PASS' if ok else 'FAIL'}  sessions/shards {final.get('_sessions')}/{final['shards_generated']}"
          f"  hits {len(hits)}  {'' if not diffs else diffs}", flush=True)
    return ok

def main():
    threads = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    random.seed(3)
    allok = True
    for X in (10**16, 10**17):
        ref = reference(X, threads)
        base = os.path.join(HERE, 'runs', f'ckpt_test_{X}')
        tag = f"10^{len(str(X))-1}"
        # 1. uninterrupted
        shutil.rmtree(base, ignore_errors=True); os.makedirs(base)
        f = session(X, base, threads); f['_sessions'] = 1
        ref_hits = hits_from(base, X)
        allok &= check(f"{tag} uninterrupted", X, f, ref, ref_hits, ref_hits)
        # 2. deadline DRAIN every 0.5-3 s until complete: nothing in flight is abandoned
        shutil.rmtree(base, ignore_errors=True); os.makedirs(base)
        n, overruns = 0, []
        while True:
            n += 1
            dl = time.time() + random.uniform(0.5, 3.0)
            f = session(X, base, threads, deadline=dl)
            overruns.append(time.time() - dl)
            if f['complete']: break
            assert f['drained'] and not f['stopped'], f
        f['_sessions'] = n
        allok &= check(f"{tag} deadline drains", X, f, ref, hits_from(base, X), ref_hits)
        print(f"      drain overrun past deadline: max {max(overruns):.2f}s, mean {sum(overruns)/len(overruns):.2f}s")
        # 3. hard stops via STOP file, abandoning in-flight shards
        shutil.rmtree(base, ignore_errors=True); os.makedirs(base)
        n = 0
        while True:
            n += 1
            f = session(X, base, threads, stopfile_after=random.uniform(0.5, 3.0))
            if f['complete']: break
        f['_sessions'] = n
        allok &= check(f"{tag} STOP-file hard stops", X, f, ref, hits_from(base, X), ref_hits)
        # 4. kill -9 crashes at random moments, then a clean finish
        shutil.rmtree(base, ignore_errors=True); os.makedirs(base)
        n = 0
        for _ in range(8):
            n += 1
            session(X, base, threads, kill_after=random.uniform(0.3, 3.0))
        while True:
            n += 1
            f = session(X, base, threads)
            if f['complete']: break
        f['_sessions'] = n
        allok &= check(f"{tag} kill -9 crashes", X, f, ref, hits_from(base, X), ref_hits)
    print("CHECKPOINT TESTS", "PASSED" if allok else "FAILED")


if __name__ == '__main__':
    main()

/*
 * cn2x.c -- exhaustive enumeration of rigid CN2s below X (compiled port of cn2_exhaustive.py).
 *
 * Same search, same counters, same switch rule as the Python reference, so the two
 * can be matched exactly.  Differences are purely mechanical:
 *   - 128-bit integers throughout (X up to ~2^100; L < X by the budget prune);
 *   - base-2 Fermat pre-filter by Montgomery arithmetic (one limb below 2^63, two above);
 *   - LIST survivors are NOT factored here: they are written out as (m, pmax, k, s) and
 *     verified by the Python driver (sympy), which is also where hits are re-verified;
 *   - workers share a LIFO task stack; nodes with m < MSPLIT are expanded into tasks, and whenever a
 *     thread is idle, a busy thread at a node with X/m > DONATE_Q donates its remaining children.
 *     Counters are exact under any split or donation pattern.
 *
 * Output (text, one record per line):
 *   H <n>                      hit found on its own prime path (L | n-1 by construction)
 *   S <m> <pmax> <k> <s>       LIST candidate that passed Fermat; needs factoring
 * Final line on stdout: JSON counters.
 *
 * Build:  clang -O3 -mcpu=native -o cn2x cn2x.c -lpthread -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
#include <stdatomic.h>

typedef unsigned __int128 u128;
typedef __int128 i128;
typedef uint64_t u64;
typedef uint32_t u32;

static u128 X;
static u64 T;
static u32 *PR;
static size_t NPR;
static double RATIO = 1.0;
static int FERMAT = 1;
static int SIEVE = 1;
static u128 NMIN = 0;           /* report only n >= NMIN (CN2X_NMIN); the search itself is unchanged */           /* 1: reject s with a prime factor q <= min(pmax, 47) before the Fermat test */
static const u64 SMALLP[] = {5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509};
static int NSMALL = 95;        /* set from X: the crossover is near 1e20 (below it the
                                * Fermat test is cheap enough that 95 divisions cost more than they save) */         /* how many of them to use: CN2X_SIEVEP */          /* 0 off, 1 on, 2 = count candidates but skip the test (timing only) */
static u128 MSPLIT = 1000;                /* nodes with m < MSPLIT are always expanded into tasks */
static u128 DONATE_Q = (u128)10000000000ull;  /* a node with X/m > DONATE_Q donates its remaining children when a thread is idle */
static atomic_int HUNGRY = 0;
static FILE *OUT;
static pthread_mutex_t out_mu = PTHREAD_MUTEX_INITIALIZER;

/* ---------------------------------------------------------------- u128 helpers */
static u128 parse_u128(const char *s) {
    u128 v = 0;
    for (; *s; s++) if (*s >= '0' && *s <= '9') v = v * 10 + (u64)(*s - '0');
    return v;
}
static void u128_str(u128 v, char *buf) {
    char tmp[64]; int i = 0;
    if (v == 0) { strcpy(buf, "0"); return; }
    while (v) { tmp[i++] = (char)('0' + (int)(v % 10)); v /= 10; }
    int j = 0; while (i) buf[j++] = tmp[--i]; buf[j] = 0;
}
static inline u64 mod128_64(u128 a, u64 b) {
    return (a >> 64) == 0 ? (u64)a % b : (u64)(a % b);
}
static u128 isqrt128(u128 n) {
    if (n == 0) return 0;
    u128 r = (u128)sqrt((double)n);
    while (r * r > n) r--;
    while ((r + 1) * (r + 1) <= n) r++;
    return r;
}
static u128 icbrt128(u128 n) {
    if (n == 0) return 0;
    u128 r = (u128)cbrt((double)n);
    while (r * r * r > n) r--;
    while ((r + 1) * (r + 1) * (r + 1) <= n) r++;
    return r;
}
static inline u64 gcd64(u64 a, u64 b) {
    if (!a) return b;
    if (!b) return a;
    int sh = __builtin_ctzll(a | b);
    a >>= __builtin_ctzll(a);
    do {
        b >>= __builtin_ctzll(b);
        if (a > b) { u64 t = a; a = b; b = t; }
        b -= a;
    } while (b);
    return a << sh;
}
static u128 modinv128(u128 a, u128 mod) {          /* gcd(a, mod) = 1, mod > 1 */
    a %= mod;
    i128 x0 = 0, x1 = 1;
    if ((mod >> 64) == 0) {
        u64 r0 = (u64)mod, r1 = (u64)a;
        while (r1) {
            u64 q = r0 / r1, r2 = r0 - q * r1; r0 = r1; r1 = r2;
            i128 x2 = x0 - (i128)q * x1; x0 = x1; x1 = x2;
        }
    } else {
        u128 r0 = mod, r1 = a;
        while (r1) {
            u128 q = r0 / r1, r2 = r0 - q * r1; r0 = r1; r1 = r2;
            i128 x2 = x0 - (i128)q * x1; x0 = x1; x1 = x2;
        }
    }
    if (x0 < 0) x0 += (i128)mod;
    return (u128)x0;
}

/* ---------------------------------------------------------------- Montgomery Fermat, base 2 */
static inline int topbit128(u128 e) {
    u64 hi = (u64)(e >> 64);
    return hi ? 127 - __builtin_clzll(hi) : 63 - __builtin_clzll((u64)e);
}
static inline u64 inv_mod_2_64(u64 n) {           /* n odd */
    u64 x = n;                                     /* correct to 3 bits */
    for (int i = 0; i < 5; i++) x *= 2 - n * x;
    return x;
}
static inline u128 montmul2(u128 a, u128 b, u128 N, u64 np) {
    u64 a0 = (u64)a, a1 = (u64)(a >> 64), b0 = (u64)b, b1 = (u64)(b >> 64);
    u64 n0 = (u64)N, n1 = (u64)(N >> 64);
    u64 t0, t1, t2, t3, m;
    u128 C;
    /* i = 0 */
    C = (u128)a0 * b0;              t0 = (u64)C; C >>= 64;
    C += (u128)a1 * b0;             t1 = (u64)C; C >>= 64;
    t2 = (u64)C; t3 = 0;
    m = t0 * np;
    C = (u128)t0 + (u128)m * n0;    C >>= 64;
    C += (u128)t1 + (u128)m * n1;   t0 = (u64)C; C >>= 64;
    C += (u128)t2;                  t1 = (u64)C; C >>= 64;
    t2 = t3 + (u64)C;
    /* i = 1 */
    C = (u128)t0 + (u128)a0 * b1;   t0 = (u64)C; C >>= 64;
    C += (u128)t1 + (u128)a1 * b1;  t1 = (u64)C; C >>= 64;
    C += (u128)t2;                  t2 = (u64)C; t3 = (u64)(C >> 64);
    m = t0 * np;
    C = (u128)t0 + (u128)m * n0;    C >>= 64;
    C += (u128)t1 + (u128)m * n1;   t0 = (u64)C; C >>= 64;
    C += (u128)t2;                  t1 = (u64)C; C >>= 64;
    t2 = t3 + (u64)C;
    u128 r = ((u128)t1 << 64) | t0;
    if (t2 || r >= N) r -= N;
    return r;
}
/* returns 1 iff 2^e == 1 (mod s); s odd, 1 < s < 2^126, e >= 1 */
static int fermat2(u128 s, u128 e) {
    int top = topbit128(e);
    if ((s >> 63) == 0) {
        u64 N = (u64)s, np = (u64)0 - inv_mod_2_64(N);
        u64 one = (u64)((((u128)1) << 64) % N), x = one;
        for (int i = top; i >= 0; i--) {
            u128 t = (u128)x * x;
            u64 mm = (u64)t * np;
            x = (u64)((t + (u128)mm * N) >> 64);
            if (x >= N) x -= N;
            if ((e >> i) & 1) { x += x; if (x >= N) x -= N; }
        }
        return x == one;
    } else {
        u128 N = s;
        u64 np = (u64)0 - inv_mod_2_64((u64)N);
        u128 one = ((u128)0 - N) % N, x = one;      /* 2^128 mod N */
        for (int i = top; i >= 0; i--) {
            x = montmul2(x, x, N, np);
            if ((e >> i) & 1) { x += x; if (x >= N) x -= N; }
        }
        return x == one;
    }
}

/* ---------------------------------------------------------------- per-thread context */
typedef struct {
    u64 nodes, children, last_tests, list_nodes, list_candidates, fermat_passes, factorisations, pruned_budget, sieved;
    char *buf; size_t len, cap;
} Ctx;

static void ctx_flush(Ctx *c) {
    if (!c->len) return;
    pthread_mutex_lock(&out_mu);
    fwrite(c->buf, 1, c->len, OUT);
    pthread_mutex_unlock(&out_mu);
    c->len = 0;
}
static void ctx_emit(Ctx *c, const char *line) {
    size_t n = strlen(line);
    if (c->len + n + 1 > c->cap) ctx_flush(c);
    memcpy(c->buf + c->len, line, n); c->len += n;
}
static void emit_hit(Ctx *c, u128 n) {
    if (n < NMIN) return;
    char a[64], line[128]; u128_str(n, a);
    snprintf(line, sizeof line, "H %s\n", a); ctx_emit(c, line);
}
static void emit_survivor(Ctx *c, u128 m, u64 pmax, int k, u128 s) {
    char a[64], b[64], line[200]; u128_str(m, a); u128_str(s, b);
    snprintf(line, sizeof line, "S %s %llu %d %s\n", a, (unsigned long long)pmax, k, b); ctx_emit(c, line);
}

static inline size_t upper(u64 key) {              /* first index with PR[i] > key */
    if (key >= 0xFFFFFFFFull) return NPR;
    size_t lo = 0, hi = NPR;
    while (lo < hi) { size_t mid = (lo + hi) / 2; if (PR[mid] <= key) lo = mid + 1; else hi = mid; }
    return lo;
}

/* ---------------------------------------------------------------- the search */
static int DEBUG_LIST = 0;
static void list_mode(Ctx *c, u128 m, u64 pmax, u128 L, int k, u128 rem) {
    c->list_nodes++;
    if (DEBUG_LIST) { char a[64], b[64]; u128_str(m, a); u128_str(L, b); fprintf(stderr, "LIST %s %s %d\n", a, b, k); }
    u128 a = modinv128(m, L);
    if (NMIN > 0) {                       /* start the progression at the first s with m*s >= NMIN */
        u128 smin = (NMIN + m - 1) / m;
        if (a < smin) a += ((smin - a + L - 1) / L) * L;
    }
    for (u128 s = a; s <= rem; s += L) {
        if (s <= pmax) continue;
        c->list_candidates++;
        if (((u64)s & 1) == 0 || mod128_64(s, 3) == 0) continue;
        if (SIEVE) {
            /* every prime factor of a completion exceeds pmax, so any q <= pmax dividing s rejects it */
            int bad = 0;
            if ((s >> 64) == 0) {
                u64 s64 = (u64)s;
                for (int j = 0; j < NSMALL && SMALLP[j] <= pmax; j++) if (s64 % SMALLP[j] == 0) { bad = 1; break; }
            } else {
                for (int j = 0; j < NSMALL && SMALLP[j] <= pmax; j++) if (mod128_64(s, SMALLP[j]) == 0) { bad = 1; break; }
            }
            if (bad) { c->sieved++; continue; }
        }
        if (FERMAT == 1) {
            int pass = fermat2(s, m * s - 1);
            c->fermat_passes += (u64)pass;
            if (!pass) continue;
        } else if (FERMAT == 2) {
            continue;                                   /* timing mode: candidates counted, tests skipped */
        }
        c->factorisations++;
        emit_survivor(c, m, pmax, k, s);
    }
}

typedef struct { size_t lo, hi, nch; u128 rem, eb; int list; } NodeHead;

/* node header shared by dfs and the splitter: returns 0 if the node has nothing below it */
static int node_head(Ctx *c, u128 m, u64 pmax, u128 L, int k, NodeHead *h) {
    c->nodes++;
    if (k >= 3 && (m - 1) % L == 0) emit_hit(c, m);
    u128 rem = (X - 1) / m;
    if (rem <= pmax) return 0;
    u128 b_last = k >= 2 ? isqrt128(m - 1) : 0;
    u128 b_second = 0;
    if (k >= 1) { u128 sq = isqrt128(rem); b_second = (m - 1 < sq) ? m - 1 : sq; }
    u128 b_early = icbrt128(rem);
    u128 mx = b_last; if (b_second > mx) mx = b_second; if (b_early > mx) mx = b_early;
    u64 tmax = mx < T ? (u64)mx : T;
    h->lo = upper(pmax); h->hi = upper(tmax);
    h->nch = h->hi > h->lo ? h->hi - h->lo : 0;
    if (h->nch == 0) return 0;                          /* no prime fits any role: no completion exists */
    h->rem = rem;
    h->eb = b_second > b_early ? b_second : b_early;
    h->list = 0;
    if (L > 1) {
        u128 nl = rem / L + 1;
        if ((double)nl <= RATIO * (double)h->nch) h->list = 1;
    }
    return 1;
}

/* one child step; returns 1 if (mt, t, Lt, k+1) should become a node, 0 otherwise; sets *stop on 'break' */
static inline int child_step(Ctx *c, u128 m, u128 L, int k, u128 eb, u64 t,
                             u128 *mt_out, u128 *Lt_out, int *stop) {
    c->children++;
    if (mod128_64(L, t) == 0) return 0;
    u64 t2 = t * t - 1;
    if (gcd64(t2, mod128_64(m, t2)) != 1) return 0;
    u128 mt = m * t;
    if (mt > X) { *stop = 1; return 0; }
    u64 g = gcd64(t2, mod128_64(L, t2));
    u128 q = L / g;
    if (q > (X - 1) / t2) { c->pruned_budget++; return 0; }
    u128 Lt = q * t2;
    if ((u128)t > eb) {
        c->last_tests++;
        if (k + 1 >= 3 && (mt - 1) % Lt == 0) emit_hit(c, mt);
        return 0;
    }
    *mt_out = mt; *Lt_out = Lt;
    return 1;
}

typedef struct { int type; u128 m, L, rem; u64 pmax; int k; } Task;   /* type 0 = dfs node, 1 = list_mode */
static void push_batch(Task *batch, size_t n);

/* hand the children PR[from..hi) of node (m, L, k) to the shared stack; counters advance exactly as in dfs */
static void donate(Ctx *c, u128 m, u128 L, int k, u128 eb, size_t from, size_t hi) {
    Task *batch = malloc((hi - from) * sizeof(Task));
    size_t n = 0;
    for (size_t i = from; i < hi; i++) {
        u128 mt, Lt; int stop = 0;
        if (child_step(c, m, L, k, eb, PR[i], &mt, &Lt, &stop)) {
            Task *t = &batch[n++];
            t->type = 0; t->m = mt; t->L = Lt; t->pmax = PR[i]; t->k = k + 1; t->rem = 0;
        }
        if (stop) break;
    }
    push_batch(batch, n);
    free(batch);
}

static void dfs(Ctx *c, u128 m, u64 pmax, u128 L, int k) {
    NodeHead h;
    if (!node_head(c, m, pmax, L, k, &h)) return;
    if (h.list) { list_mode(c, m, pmax, L, k, h.rem); return; }
    int may_donate = (X / m) > DONATE_Q;
    for (size_t i = h.lo; i < h.hi; i++) {
        u128 mt, Lt; int stop = 0;
        if (child_step(c, m, L, k, h.eb, PR[i], &mt, &Lt, &stop)) dfs(c, mt, PR[i], Lt, k + 1);
        if (stop) break;
        if (may_donate && i + 1 < h.hi && atomic_load_explicit(&HUNGRY, memory_order_relaxed) > 0) {
            donate(c, m, L, k, h.eb, i + 1, h.hi);
            return;
        }
    }
}

/* ---------------------------------------------------------------- parallel: shared LIFO task stack
 * A worker pops a task.  LIST tasks and DFS tasks with m >= MSPLIT run to completion locally.
 * A DFS task with m < MSPLIT is only EXPANDED: its node header and child loop run locally
 * (counting exactly as dfs would) and each child node is pushed back as a task, so heavy
 * subtrees keep dividing wherever they are.  LIFO keeps the pending set small, like a DFS.
 */
static Task *QUEUE;
static size_t QLEN = 0, QCAP = 0;
static int ACTIVE = 0;
static pthread_mutex_t q_mu = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t q_cv = PTHREAD_COND_INITIALIZER;
static u64 TASKS_RUN = 0, TASKS_EXPANDED = 0;

static void push_batch(Task *batch, size_t n) {
    if (!n) return;
    pthread_mutex_lock(&q_mu);
    if (QLEN + n > QCAP) {
        while (QLEN + n > QCAP) QCAP = QCAP ? 2 * QCAP : 1 << 16;
        QUEUE = realloc(QUEUE, QCAP * sizeof(Task));
    }
    /* push in reverse so the smallest child is popped first (preserves DFS order locally) */
    for (size_t i = n; i-- > 0;) QUEUE[QLEN++] = batch[i];
    pthread_cond_broadcast(&q_cv);
    pthread_mutex_unlock(&q_mu);
}

static void expand(Ctx *c, const Task *tk, Task **batch, size_t *bcap) {
    NodeHead h;
    if (!node_head(c, tk->m, tk->pmax, tk->L, tk->k, &h)) return;
    if (h.list) { list_mode(c, tk->m, tk->pmax, tk->L, tk->k, h.rem); return; }
    size_t n = 0;
    for (size_t i = h.lo; i < h.hi; i++) {
        u128 mt, Lt; int stop = 0;
        if (child_step(c, tk->m, tk->L, tk->k, h.eb, PR[i], &mt, &Lt, &stop)) {
            if (n == *bcap) { *bcap *= 2; *batch = realloc(*batch, *bcap * sizeof(Task)); }
            Task *t = &(*batch)[n++];
            t->type = 0; t->m = mt; t->L = Lt; t->pmax = PR[i]; t->k = tk->k + 1; t->rem = 0;
        }
        if (stop) break;
    }
    push_batch(*batch, n);
}

static void *worker(void *arg) {
    Ctx *c = (Ctx *)arg;
    size_t bcap = 1024;
    Task *batch = malloc(bcap * sizeof(Task));
    for (;;) {
        pthread_mutex_lock(&q_mu);
        while (QLEN == 0 && ACTIVE > 0) {
            atomic_fetch_add(&HUNGRY, 1);
            pthread_cond_wait(&q_cv, &q_mu);
            atomic_fetch_sub(&HUNGRY, 1);
        }
        if (QLEN == 0) { pthread_cond_broadcast(&q_cv); pthread_mutex_unlock(&q_mu); break; }
        Task tk = QUEUE[--QLEN];
        ACTIVE++;
        pthread_mutex_unlock(&q_mu);

        int expanded = 0;
        if (tk.type == 1) list_mode(c, tk.m, tk.pmax, tk.L, tk.k, tk.rem);
        else if (tk.m < MSPLIT) { expand(c, &tk, &batch, &bcap); expanded = 1; }
        else dfs(c, tk.m, tk.pmax, tk.L, tk.k);

        pthread_mutex_lock(&q_mu);
        ACTIVE--;
        if (expanded) TASKS_EXPANDED++; else TASKS_RUN++;
        if (QLEN == 0 && ACTIVE == 0) pthread_cond_broadcast(&q_cv);
        pthread_mutex_unlock(&q_mu);
    }
    free(batch);
    ctx_flush(c);
    return NULL;
}

/* ---------------------------------------------------------------- setup */
static void sieve(u64 limit) {
    NPR = 0; PR = NULL;
    if (limit < 5) return;
    uint8_t *comp = calloc(limit + 1, 1);
    size_t cap = (size_t)(1.3 * limit / (log((double)limit) + 1)) + 1000;
    PR = malloc(cap * sizeof(u32));
    for (u64 p = 2; p * p <= limit; p++)
        if (!comp[p]) for (u64 q = p * p; q <= limit; q += p) comp[q] = 1;
    for (u64 p = 5; p <= limit; p++) if (!comp[p]) {
        if (NPR == cap) { cap *= 2; PR = realloc(PR, cap * sizeof(u32)); }
        PR[NPR++] = (u32)p;
    }
    free(comp);
}
static Ctx *new_ctx(void) {
    Ctx *c = calloc(1, sizeof(Ctx));
    c->cap = 1 << 20; c->buf = malloc(c->cap);
    return c;
}
static double now(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + ts.tv_nsec * 1e-9; }

static void bench_fermat(int bits, long iters) {
    /* random odd moduli of the given size, exponent ~ X-sized */
    u128 state = 0x9E3779B97F4A7C15ull;
    double t0 = now(); long passes = 0;
    for (long i = 0; i < iters; i++) {
        state = state * 6364136223846793005ull + 1442695040888963407ull;
        u128 s = (state >> (128 - bits)) | 1 | ((u128)1 << (bits - 1));
        u128 e = state ^ (state >> 29);
        if (e == 0) e = 1;
        passes += fermat2(s, e);
    }
    double dt = now() - t0;
    printf("{\"bench\":\"fermat\",\"bits\":%d,\"iters\":%ld,\"ns_per_test\":%.1f,\"passes\":%ld}\n",
           bits, iters, dt / iters * 1e9, passes);
}
static void check_arith(void) {
    /* read "op a b" lines, print one result per line -- validation of the 128-bit helpers against Python */
    char op[16], a[128], b[128], out[64];
    while (scanf("%15s %127s %127s", op, a, b) == 3) {
        u128 x = parse_u128(a), y = parse_u128(b), r = 0;
        if (!strcmp(op, "inv")) r = modinv128(x, y);
        else if (!strcmp(op, "isqrt")) r = isqrt128(x);
        else if (!strcmp(op, "icbrt")) r = icbrt128(x);
        else if (!strcmp(op, "mod64")) r = mod128_64(x, (u64)y);
        u128_str(r, out); printf("%s\n", out);
    }
}
static void check_fermat(void) {
    /* read "s e" pairs from stdin, print 0/1 per line -- for validation against Python pow */
    char a[128], b[128];
    while (scanf("%127s %127s", a, b) == 2) printf("%d\n", fermat2(parse_u128(a), parse_u128(b)));
}

int main(int argc, char **argv) {
    if (argc >= 2 && strcmp(argv[1], "bench-fermat") == 0) {
        bench_fermat(argc > 2 ? atoi(argv[2]) : 80, argc > 3 ? atol(argv[3]) : 1000000); return 0;
    }
    if (argc >= 2 && strcmp(argv[1], "check-fermat") == 0) { check_fermat(); return 0; }
    if (argc >= 2 && strcmp(argv[1], "check-arith") == 0) { check_arith(); return 0; }
    if (argc < 3) {
        fprintf(stderr, "usage: cn2x X out.txt [threads=1] [ratio=1.0] [fermat=1] [donate_q=1e10] [msplit=1000] [sieve=1]\n"
                        "       cn2x bench-fermat BITS ITERS | cn2x check-fermat < pairs\n");
        return 2;
    }
    DEBUG_LIST = getenv("CN2X_DEBUG_LIST") != NULL;
    X = parse_u128(argv[1]);
    OUT = fopen(argv[2], "w");
    int threads = argc > 3 ? atoi(argv[3]) : 1;
    if (argc > 4) RATIO = atof(argv[4]);
    if (argc > 5) FERMAT = atoi(argv[5]);
    if (argc > 6) DONATE_Q = parse_u128(argv[6]);
    if (argc > 7) MSPLIT = parse_u128(argv[7]);
    SIEVE = argc > 8 ? atoi(argv[8]) : 1;
    T = (u64)icbrt128(X - 1);
    {   /* EXPERIMENT ONLY: CN2X_TCAP caps the prime table, which makes the search CONDITIONAL on
         * every factor being <= the cap (the lazy sweep's H-factor).  Never use it for a real run. */
        const char *nm = getenv("CN2X_NMIN");
        if (nm) NMIN = parse_u128(nm);
        NSMALL = (X < (u128)10000000000ull * 10000000000ull) ? 13 : 95;   /* 1e20 */
        const char *ns = getenv("CN2X_SIEVEP");
        if (ns) { int v = atoi(ns); if (v > 0 && v <= 95) NSMALL = v; }
        const char *cap = getenv("CN2X_TCAP");
        if (cap) { u64 c = strtoull(cap, NULL, 10); if (c < T) T = c; }
    }
    double t0 = now();
    sieve(T);
    double t_sieve = now() - t0;

    Ctx total = {0};
    t0 = now();
    if (threads <= 1) {
        Ctx *c = new_ctx();
        dfs(c, 1, 3, 1, 0);
        ctx_flush(c);
        total = *c;
    } else {
        Task root = {0, 1, 1, 0, 3, 0};
        push_batch(&root, 1);
        pthread_t *th = malloc(threads * sizeof(pthread_t));
        Ctx **cs = malloc(threads * sizeof(Ctx *));
        for (int i = 0; i < threads; i++) { cs[i] = new_ctx(); pthread_create(&th[i], NULL, worker, cs[i]); }
        for (int i = 0; i < threads; i++) pthread_join(th[i], NULL);
        for (int i = 0; i < threads; i++) {
            total.nodes += cs[i]->nodes; total.children += cs[i]->children; total.last_tests += cs[i]->last_tests;
            total.list_nodes += cs[i]->list_nodes; total.list_candidates += cs[i]->list_candidates;
            total.fermat_passes += cs[i]->fermat_passes; total.factorisations += cs[i]->factorisations;
            total.pruned_budget += cs[i]->pruned_budget; total.sieved += cs[i]->sieved;
        }
    }
    double t_search = now() - t0;
    fclose(OUT);
    char xs[64]; u128_str(X, xs);
    printf("{\"X\":\"%s\",\"T\":%llu,\"primes\":%zu,\"threads\":%d,\"ratio\":%g,\"fermat\":%d,\"sieve\":%d,\"sieve_primes\":%d,"
           "\"nodes\":%llu,\"children\":%llu,\"last_tests\":%llu,\"list_nodes\":%llu,"
           "\"list_candidates\":%llu,\"fermat_passes\":%llu,\"factorisations\":%llu,\"pruned_budget\":%llu,"
           "\"sieved\":%llu,\"tasks_run\":%llu,\"tasks_expanded\":%llu,\"sieve_seconds\":%.3f,\"search_seconds\":%.3f}\n",
           xs, (unsigned long long)T, NPR, threads, RATIO, FERMAT, SIEVE, NSMALL,
           (unsigned long long)total.nodes, (unsigned long long)total.children, (unsigned long long)total.last_tests,
           (unsigned long long)total.list_nodes, (unsigned long long)total.list_candidates,
           (unsigned long long)total.fermat_passes, (unsigned long long)total.factorisations,
           (unsigned long long)total.pruned_budget, (unsigned long long)total.sieved, (unsigned long long)TASKS_RUN, (unsigned long long)TASKS_EXPANDED, t_sieve, t_search);
    return 0;
}

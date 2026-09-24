/* bf_order2.c -- brute-force reference for validate.py: every odd n < N, 3 !| n, squarefree, >= 3 primes,
 * with n == 1 or p (mod (p^2-1)/D) for every p | n.  Prints "n has_type_p".  Memory: 4 N bytes (N = 1e9: 4 GB).
 * Build: cc -O3 -o bf_order2 bf_order2.c   Usage: bf_order2 N D */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
int main(int argc, char **argv) {
    uint64_t N = strtoull(argv[1], 0, 10), D = strtoull(argv[2], 0, 10);
    uint32_t *spf = calloc(N, 4);
    for (uint64_t i = 2; i < N; i++) if (!spf[i]) for (uint64_t j = i; j < N; j += i) if (!spf[j]) spf[j] = (uint32_t)i;
    uint64_t cnt = 0, cntp = 0;
    for (uint64_t n = 5; n < N; n += 2) {
        if (n % 3 == 0) continue;
        uint64_t x = n, ps[40]; int k = 0, ok = 1;
        while (x > 1) { uint64_t p = spf[x]; x /= p; if (x % p == 0) { ok = 0; break; } ps[k++] = p; }
        if (!ok || k < 3) continue;
        int hasp = 0;
        for (int i = 0; i < k && ok; i++) {
            uint64_t p = ps[i], g = (p * p - 1) / D, r = n % g;
            if (r == 1 % g) continue;
            if (r == p % g) { hasp = 1; continue; }
            ok = 0;
        }
        if (!ok) continue;
        cnt++; cntp += hasp;
        printf("%llu %d\n", (unsigned long long)n, hasp);
    }
    fprintf(stderr, "N=%llu D=%llu: %llu solutions, %llu with a type-p prime\n", (unsigned long long)N, (unsigned long long)D,
            (unsigned long long)cnt, (unsigned long long)cntp);
    return 0;
}

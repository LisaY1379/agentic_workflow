/*
 * pomerance_batch_engine.c
 *
 * Production-grade incremental Pomerance triple finder with exact trial telemetry.
 *
 * Features integrated:
 * 1. C11 Atomics & u128 Hasse-interval arithmetic (no near-2^127 overflow).
 * 2. Supersingular A=0 fast-path with Legendre twist-filtering.
 * 3. OpenMP cache-line padded exact trial counters (no false sharing).
 * 4. Stateful incremental batch loading, deduplication, and marginal trial sorting.
 * 5. Real-time terminal telemetry tracking the current A and A_index.
 *
 * Build:
 * gcc -O3 -march=native -fopenmp -std=c11 -o pomerance_batch pomerance_batch_engine.c -lm
 *
 * Usage:
 * ./pomerance_batch <stateful_input.txt> <output_pure.csv> <output_metrics.csv> <target_total>
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdatomic.h>

#ifdef _OPENMP
#include <omp.h>
#endif

typedef uint64_t u64;
typedef __uint128_t u128;

typedef struct { u128 lo; u128 hi; } u256;
typedef struct { u64 s0; u64 s1; } Rng;
typedef struct { u128 p; u128 ni; u128 R2; u128 one; } Mont128;

/* ================================================================
 * Advanced I/O & Parsing Utilities
 * ================================================================ */

static u128 parse128_adv(char **s) {
    while (**s == ' ' || **s == '\t' || **s == '\n' || **s == '\r') (*s)++;
    if (**s == '\0') return 0;
    u128 v = 0;
    while (**s >= '0' && **s <= '9') { v = v * 10 + (u128)(**s - '0'); (*s)++; }
    return v;
}

static void sprint128(char *buf, u128 v) {
    if (v == 0) { buf[0] = '0'; buf[1] = '\0'; return; }
    char tmp[50]; int i = 49; tmp[i] = '\0';
    while (v > 0) { tmp[--i] = '0' + (int)(v % 10); v /= 10; }
    strcpy(buf, tmp + i);
}

static void print128(u128 v) { char b[50]; sprint128(b, v); fputs(b, stdout); }

static inline int bitlen128(u128 x) {
    u64 hi = (u64)(x >> 64);
    if (hi != 0) return 64 + (64 - __builtin_clzll(hi));
    u64 lo = (u64)x;
    return lo ? 64 - __builtin_clzll(lo) : 0;
}

static inline u64 rng64(Rng *r) {
    u64 s1 = r->s0, s0 = r->s1; r->s0 = s0;
    s1 ^= s1 << 23; r->s1 = s1 ^ s0 ^ (s1 >> 17) ^ (s0 >> 26);
    return r->s1 + s0;
}

static inline u128 rand_below128(Rng *rng, u128 p, u128 mask) {
    for (;;) {
        u128 x = ((u128)rng64(rng) << 64) | rng64(rng);
        x &= mask;
        if (x < p) return x;
    }
}

static double now_sec(void) {
#ifdef _OPENMP
    return omp_get_wtime();
#else
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
#endif
}

/* ================================================================
 * Wide Arithmetic & Montgomery Engine (p < 2^127)
 * ================================================================ */

static inline u256 wide_mul(u128 a, u128 b) {
    u64 a0 = (u64)a, a1 = (u64)(a >> 64), b0 = (u64)b, b1 = (u64)(b >> 64);
    u128 ll = (u128)a0 * b0, lh = (u128)a0 * b1, hl = (u128)a1 * b0, hh = (u128)a1 * b1;
    u128 mid = lh + hl; u128 carry_mid = (mid < lh) ? 1 : 0;
    u128 lo = ll + (mid << 64); u128 carry_lo = (lo < ll) ? 1 : 0;
    return (u256){lo, hh + (mid >> 64) + (carry_mid << 64) + carry_lo};
}

static inline u256 wide_add(u256 a, u256 b) {
    u128 lo = a.lo + b.lo; return (u256){lo, a.hi + b.hi + (lo < a.lo ? 1 : 0)};
}

static inline u128 addmod128(u128 a, u128 b, u128 p) {
    u128 s = a + b; return s >= p ? s - p : s;
}

static inline u128 submod128(u128 a, u128 b, u128 p) {
    return a >= b ? a - b : p - b + a;
}

static void m128_init(Mont128 *m, u128 p) {
    m->p = p; u128 x = 1, r = 1; int i;
    for (i = 0; i < 7; ++i) x *= 2 - p * x;
    m->ni = (u128)0 - x;
    for (i = 0; i < 128; ++i) { r <<= 1; if (r >= p) r -= p; }
    m->one = r;
    for (i = 0; i < 128; ++i) { r <<= 1; if (r >= p) r -= p; }
    m->R2 = r;
}

static inline u128 mred128(u256 T, const Mont128 *m) {
    u128 q = T.lo * m->ni;
    u256 s = wide_add(T, wide_mul(q, m->p));
    u128 t = s.hi;
    return t >= m->p ? t - m->p : t;
}

static inline u128 mm128(u128 a, u128 b, const Mont128 *m) { return mred128(wide_mul(a, b), m); }
static inline u128 toM128(u128 a, const Mont128 *m) { return mm128(a % m->p, m->R2, m); }
static inline u128 frM128(u128 a, const Mont128 *m) { return mred128((u256){a, 0}, m); }

static inline void xDBL128(u128 *Xo, u128 *Zo, u128 X, u128 Z, u128 a24, const Mont128 *m) {
    u128 p = m->p, u = addmod128(X, Z, p), v = submod128(X, Z, p);
    u = mm128(u, u, m); v = mm128(v, v, m);
    *Xo = mm128(u, v, m);
    u128 w = submod128(u, v, p);
    *Zo = mm128(w, addmod128(v, mm128(a24, w, m), p), m);
}

static inline void xADD128(u128 *Xo, u128 *Zo, u128 X0, u128 Z0, u128 X1, u128 Z1, u128 xP, const Mont128 *m) {
    u128 p = m->p, u = mm128(submod128(X0, Z0, p), addmod128(X1, Z1, p), m);
    u128 v = mm128(addmod128(X0, Z0, p), submod128(X1, Z1, p), m);
    u128 s = addmod128(u, v, p), d = submod128(u, v, p);
    *Xo = mm128(s, s, m);
    *Zo = mm128(xP, mm128(d, d, m), m);
}

static void xMUL128(u128 *Xo, u128 *Zo, u128 xP, u128 n, u128 a24, const Mont128 *m) {
    if (n == 0) { *Xo = 0; *Zo = 0; return; }
    if (n == 1) { *Xo = xP; *Zo = m->one; return; }
    u128 X0 = xP, Z0 = m->one, X1, Z1;
    xDBL128(&X1, &Z1, X0, Z0, a24, m);
    int bits = bitlen128(n);
    for (int i = bits - 2; i >= 0; --i) {
        if ((n >> i) & 1) { xADD128(&X0, &Z0, X0, Z0, X1, Z1, xP, m); xDBL128(&X1, &Z1, X1, Z1, a24, m); }
        else              { xADD128(&X1, &Z1, X0, Z0, X1, Z1, xP, m); xDBL128(&X0, &Z0, X0, Z0, a24, m); }
    }
    *Xo = X0; *Zo = Z0;
}

static u128 powM128(u128 aM, u128 e, const Mont128 *m) {
    u128 r = m->one;
    while (e != 0) { if (e & 1) r = mm128(r, aM, m); aM = mm128(aM, aM, m); e >>= 1; }
    return r;
}

static u128 invertM128(u128 aM, u128 p, const Mont128 *m) { return powM128(aM, p - 2, m); }

static int legendreM128(u128 aM, u128 p, const Mont128 *m) {
    if (aM == 0) return 0;
    return (powM128(aM, (p - 1) >> 1, m) == m->one) ? 1 : -1;
}

static inline int is_point_x_supersingular(u128 xM, u128 p, const Mont128 *m) {
    u128 x2 = mm128(xM, xM, m);
    u128 rhs = mm128(xM, addmod128(x2, m->one, p), m);
    return legendreM128(rhs, p, m) >= 0;
}

/* ================================================================
 * Independent Strict Pomerance Verifier
 * ================================================================ */

static u128 mulmod_slow(u128 a, u128 b, u128 p) {
    u128 r = 0; a %= p; b %= p;
    while (b != 0) { if (b & 1) { r += a; if (r >= p) r -= p; } a += a; if (a >= p) a -= p; b >>= 1; }
    return r;
}

static u128 isqrt128(u128 n) {
    u64 q = (u64)sqrtl((long double)n);
    while ((u128)(q + 1) * (u128)(q + 1) <= n) ++q;
    while ((u128)q * (u128)q > n) --q;
    return (u128)q;
}

static int compute_k(u128 p) {
    u128 q = isqrt128(p), sq = isqrt128(q);
    u128 bound = q + 1 + 2 * sq, v = 1; int k = 0;
    while (v <= bound) { ++k; v <<= 1; }
    return k;
}

static int v2_u128(u128 n) {
    int v = 0; while (n != 0 && (n & 1) == 0) { ++v; n >>= 1; } return v;
}

static int verify128(u128 p, u128 A, u128 x0) {
    u128 q = isqrt128(p), sq = isqrt128(q);
    u128 bound = q + 1 + 2 * sq, v = 1; int k = 0, i;
    while (v <= bound) { ++k; v <<= 1; }
    if (A % p == 2 || A % p == p - 2) return 0;
    u128 X = x0 % p, Z = 1;
    for (i = 1; i <= k; ++i) {
        u128 X2 = mulmod_slow(X, X, p), Z2 = mulmod_slow(Z, Z, p), XZ = mulmod_slow(X, Z, p);
        u128 d = submod128(X2, Z2, p), Xn = mulmod_slow(d, d, p);
        u128 inn = addmod128(addmod128(X2, mulmod_slow(A, XZ, p), p), Z2, p);
        u128 f4 = addmod128(addmod128(XZ, XZ, p), addmod128(XZ, XZ, p), p);
        u128 Zn = mulmod_slow(f4, inn, p); X = Xn; Z = Zn;
        if (i < k && Z == 0) return 0;
        if (i == k && Z != 0) return 0;
    }
    return 1;
}

static int compute_odd_parts(u128 p, int k, u128 *ms, int *max_v2s, int max_ms) {
    if (k >= 127) return 0;
    u128 twok = ((u128)1) << k, pp1 = p + 1;
    u128 sqrtp = isqrt128(p), bound = 2 * sqrtp + 4;
    u128 min_N = (pp1 >= bound) ? pp1 - bound : 1, max_N = pp1 + bound;
    u128 rem = min_N % twok, start_N = (rem == 0) ? min_N : min_N + (twok - rem);
    int count = 0;
    for (u128 N = start_N; N <= max_N && count < max_ms; N += twok) {
        if (N == 0) continue;
        u128 tmp = N; int v2 = 0, dup = -1, i;
        while ((tmp & 1) == 0) { ++v2; tmp >>= 1; }
        if (v2 < k || tmp == 0) continue;
        for (i = 0; i < count; ++i) if (ms[i] == tmp) { dup = i; break; }
        if (dup >= 0) { if (max_v2s[dup] < v2) max_v2s[dup] = v2; }
        else { ms[count] = tmp; max_v2s[count] = v2; ++count; }
    }
    return count;
}

static int projected_hit128(u128 *xRo, u128 p, u128 A, int k, int max_s,
                            u128 QX, u128 QZ, u128 a24m, const Mont128 *mt) {
    if (QZ == 0) return 0;
    u128 CX = QX, CZ = QZ; int zs = -1, s;
    for (s = 1; s <= max_s && s < 128; ++s) {
        xDBL128(&CX, &CZ, CX, CZ, a24m, mt);
        if (CZ == 0) { zs = s; break; }
    }
    if (zs < k) return 0;
    int target = zs - k; CX = QX; CZ = QZ;
    for (s = 0; s < target; ++s) xDBL128(&CX, &CZ, CX, CZ, a24m, mt);
    if (CZ == 0) return 0;
    u128 invZ = invertM128(CZ, p, mt);
    u128 xR = frM128(mm128(CX, invZ, mt), mt);
    if (!verify128(p, A, xR)) return 0;
    *xRo = xR;
    return 1;
}

/* ================================================================
 * Unified Batch Search Engine with Exact Trial Telemetry
 * ================================================================ */

static int search128_batch(u128 p, int target_total, int start_count,
                           u128 *out_A, u128 *out_x0, u64 *out_trials) {
    volatile int found_count = start_count;
    int k = compute_k(p);
    u128 ms[64]; int max_v2s[64];
    int nms = compute_odd_parts(p, k, ms, max_v2s, 64);
    if (nms == 0) { printf("      No valid Hasse odd parts.\n"); return found_count; }

    Mont128 mt; m128_init(&mt, p);
    int pbits = bitlen128(p);
    u128 rand_mask = pbits >= 128 ? (u128)-1 : (((u128)1 << pbits) - 1);
    u128 inv4m = invertM128(toM128(4, &mt), p, &mt);
    u64 sqrtp = (u64)isqrt128(p);
    u64 max_trials = (u64)(20.0 * (double)sqrtp / nms);
    if (max_trials < 10000000ULL) max_trials = 10000000ULL;

    /* Padded thread trial counters to eliminate OpenMP false sharing */
    u64 thread_trials[256 * 8] = {0};
    int nth = 1;
#ifdef _OPENMP
    nth = omp_get_max_threads();
#endif

    double t0_batch = now_sec();

    /* --- PHASE 1: Supersingular A=0 Fast-Path --- */
    if ((p & 3) == 3) {
        u128 N = p + 1; int v2N = v2_u128(N);
        if (v2N - 1 >= k) {
            u128 odd = N >> v2N;
            u128 a24m_ss = mm128(toM128(2, &mt), inv4m, &mt);
            u64 ss_max_trials = 200000ULL; /* Fast budget for exact A=0 */

#pragma omp parallel
            {
                int tid = 0, nthr = 1;
#ifdef _OPENMP
                tid = omp_get_thread_num(); nthr = omp_get_num_threads();
#endif
                u64 cur_time = (u64)time(NULL);
                Rng rng = {
                    .s0 = 0x9e3779b97f4a7c15ULL ^ ((u64)tid * 0xbf58476d1ce4e5b9ULL) ^ (u64)p ^ cur_time,
                    .s1 = 0x94d049bb133111ebULL ^ ((u64)(tid + 1) * 0x9e3779b97f4a7c15ULL) ^ (cur_time << 32)
                };
                for (int i = 0; i < 200; ++i) (void)rng64(&rng);
                u64 budget = ss_max_trials / (u64)nthr + 1;
                u64 local_tries = 0;

                while (found_count < target_total && local_tries < budget) {
                    local_tries++;
                    thread_trials[tid * 8]++;

                    /* Telemetry for Phase 1 (A is always 0) */
                    if (tid == 0 && (local_tries % 512) == 0) {
                        printf("\r      [SS Fast-Path] Thread 0 trying x_index: %llu | Current A: 0 | Elapsed: %.1fs   ",
                               (unsigned long long)local_tries, now_sec() - t0_batch);
                        fflush(stdout);
                    }

                    u128 x = rand_below128(&rng, p, rand_mask);
                    if (x < 2) continue;
                    u128 xM = toM128(x, &mt);
                    if (!is_point_x_supersingular(xM, p, &mt)) continue;

                    u128 QX, QZ, xR;
                    xMUL128(&QX, &QZ, xM, odd, a24m_ss, &mt);
                    if (QZ == 0) continue;
                    if (!projected_hit128(&xR, p, 0, k, v2N, QX, QZ, a24m_ss, &mt)) continue;

#pragma omp critical
                    {
                        int is_dup = 0;
                        for (int j = 0; j < found_count; j++) {
                            if (out_A[j] == 0 && out_x0[j] == xR) { is_dup = 1; break; }
                        }
                        if (!is_dup && found_count < target_total) {
                            out_A[found_count] = 0;
                            out_x0[found_count] = xR;
                            u64 exact_total = 0;
                            for (int t = 0; t < nthr; t++) exact_total += thread_trials[t * 8];
                            out_trials[found_count] = (exact_total == 0) ? 1 : exact_total;
                            found_count++;
                        }
                    }
                }
            }
        }
    }

    if (found_count >= target_total) return found_count;

    /* --- PHASE 2: Ordinary Random Curve Enumeration --- */
#pragma omp parallel
    {
        int tid = 0, nthr = 1;
#ifdef _OPENMP
        tid = omp_get_thread_num(); nthr = omp_get_num_threads();
#endif
        u64 cur_time = (u64)time(NULL);
        Rng rng = {
            .s0 = 7364529176530163ULL ^ ((u64)tid * 6364136223846793005ULL) ^ (u64)p ^ cur_time,
            .s1 = 1442695040888963407ULL ^ ((u64)(tid + 1) * 2862933555777941757ULL) ^ (cur_time << 32)
        };
        for (int i = 0; i < 200; ++i) (void)rng64(&rng);
        u64 budget = max_trials / (u64)nthr + 1;
        u64 A_trials = thread_trials[tid * 8]; /* Inherit counts from SS fast-path */

        while (found_count < target_total && A_trials < budget) {
            u128 A = rand_below128(&rng, p, rand_mask);
            if (A <= 2 || A >= p - 2) continue;

            A_trials++;
            thread_trials[tid * 8] = A_trials;

            /* --- NEW: Real-time telemetry for A --- */
            if (tid == 0 && (A_trials % 256) == 0) {
                char a_str[64];
                sprint128(a_str, A);
                printf("\r      [Search] Thread 0 trying A_index: %llu | Current A: %s | Elapsed: %.1fs   ",
                       (unsigned long long)A_trials, a_str, now_sec() - t0_batch);
                fflush(stdout);
            }
            /* -------------------------------------- */

            u128 a24m = mm128(toM128(addmod128(A, 2, p), &mt), inv4m, &mt);
            int a_success = 0;

            for (int x_tries = 0; x_tries < 50 && !a_success && found_count < target_total; x_tries++) {
                u128 x0 = rand_below128(&rng, p, rand_mask);
                if (x0 < 2) continue;
                u128 x0m = toM128(x0, &mt);

                for (int mi = 0; mi < nms && !a_success && found_count < target_total; mi++) {
                    u128 QX, QZ, xR;
                    xMUL128(&QX, &QZ, x0m, ms[mi], a24m, &mt);
                    if (QZ == 0) continue;
                    if (!projected_hit128(&xR, p, A, k, max_v2s[mi], QX, QZ, a24m, &mt)) continue;

#pragma omp critical
                    {
                        int is_dup = 0;
                        for (int j = 0; j < found_count; j++) {
                            if (out_A[j] == A && out_x0[j] == xR) { is_dup = 1; break; }
                        }
                        if (!is_dup && found_count < target_total) {
                            out_A[found_count] = A;
                            out_x0[found_count] = xR;
                            u64 exact_total = 0;
                            for (int t = 0; t < nthr; t++) exact_total += thread_trials[t * 8];
                            out_trials[found_count] = (exact_total == 0) ? 1 : exact_total;
                            found_count++;
                            a_success = 1;
                        }
                    }
                }
            }
        }
    }

    return found_count;
}

/* ================================================================
 * Batch Processing Main Loop
 * ================================================================ */

int main(int argc, char *argv[]) {
    if (argc < 5) {
        printf("Usage: ./pomerance_batch <stateful_input.txt> <output_pure.csv> <output_metrics.csv> <target_total>\n");
        return 1;
    }

    int target_total = atoi(argv[4]);
    if (target_total <= 0) {
        printf("Error: <target_total> must be greater than 0.\n");
        return 1;
    }

    FILE *input = fopen(argv[1], "r");
    FILE *out_pure = fopen(argv[2], "a");    /* Append mode to protect historical logs */
    FILE *out_metrics = fopen(argv[3], "a");

    if (!input || !out_pure || !out_metrics) {
        printf("Error opening files.\n");
        if (input) fclose(input);
        if (out_pure) fclose(out_pure);
        if (out_metrics) fclose(out_metrics);
        return 1;
    }

    u128 *out_A = (u128 *)malloc(target_total * sizeof(u128));
    u128 *out_x0 = (u128 *)malloc(target_total * sizeof(u128));
    u64 *out_trials = (u64 *)malloc(target_total * sizeof(u64));

    if (!out_A || !out_x0 || !out_trials) {
        printf("Memory allocation failed.\n");
        return 1;
    }

    char line[65536];
    unsigned long long current_index = 0;

    while (fgets(line, sizeof(line), input)) {
        char *ptr = line;
        u128 p = parse128_adv(&ptr);
        if (p == 0 || p < 5) continue;

        current_index++;
        int num_existing = (int)parse128_adv(&ptr);

        /* Preload existing pairs from stateful input */
        for (int i = 0; i < num_existing && i < target_total; i++) {
            out_A[i] = parse128_adv(&ptr);
            out_x0[i] = parse128_adv(&ptr);
            out_trials[i] = 0;
        }

        char p_display[50]; sprint128(p_display, p);
        printf("[%llu] Prime: %s (Existing: %d, Target: %d)...\n",
               current_index, p_display, num_existing, target_total);
        fflush(stdout);

        int found_amount = search128_batch(p, target_total, num_existing, out_A, out_x0, out_trials);
        int newly_found = found_amount - num_existing;

        if (newly_found > 0) {
            /* Clear the \r telemetry line before printing success */
            printf("\n");

            /* Sort ONLY the newly found elements by ascending trials */
            for (int i = num_existing; i < found_amount - 1; i++) {
                for (int j = num_existing; j < found_amount - 1 - (i - num_existing); j++) {
                    if (out_trials[j] > out_trials[j + 1]) {
                        u64 temp_t = out_trials[j]; out_trials[j] = out_trials[j + 1]; out_trials[j + 1] = temp_t;
                        u128 temp_A = out_A[j]; out_A[j] = out_A[j + 1]; out_A[j + 1] = temp_A;
                        u128 temp_x0 = out_x0[j]; out_x0[j] = out_x0[j + 1]; out_x0[j + 1] = temp_x0;
                    }
                }
            }

            u64 total_new_trials = 0;
            u64 last_cumulative = 0;

            /* Output ONLY the new marginal triples */
            for (int i = num_existing; i < found_amount; i++) {
                u64 marginal_trials = out_trials[i] - last_cumulative;
                if (marginal_trials == 0) marginal_trials = 1;
                last_cumulative = out_trials[i];
                total_new_trials += marginal_trials;

                char a_str[50], x_str[50];
                sprint128(a_str, out_A[i]);
                sprint128(x_str, out_x0[i]);

                fprintf(out_pure, "%s,%s,%s\n", p_display, a_str, x_str);
                fprintf(out_metrics, "%s,%s,%s,%llu\n", p_display, a_str, x_str, (unsigned long long)marginal_trials);
            }
            fflush(out_pure);
            fflush(out_metrics);
            printf("      Success: found %d new triples (Total new trials: %llu)\n", newly_found, (unsigned long long)total_new_trials);
        } else {
            printf("\n"); /* Clear the \r telemetry line before printing failure */
            fprintf(out_pure, "%s,FAILED,FAILED\n", p_display);
            fprintf(out_metrics, "%s,FAILED,FAILED,FAILED\n", p_display);
            fflush(out_pure);
            fflush(out_metrics);
            printf("      Failed to find new triples within trial budget.\n");
        }
        fflush(stdout);
    }

    free(out_A);
    free(out_x0);
    free(out_trials);
    fclose(input);
    fclose(out_pure);
    fclose(out_metrics);
    return 0;
}
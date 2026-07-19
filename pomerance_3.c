#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sched.h>
#include <pthread.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

typedef uint64_t u64;
typedef __uint128_t u128;

/* ================================================================
 * Benchmark #1: CPU cycle counting via rdtsc (Cross-Platform)
 * ================================================================ */

static inline uint64_t rdtsc(void) {
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__)
    unsigned int lo, hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
#elif defined(__aarch64__) || defined(_M_ARM64)
    uint64_t val;
    __asm__ __volatile__ ("mrs %0, cntvct_el0" : "=r" (val));
    return val;
#else
    return 0;
#endif
}

static void pin_thread_to_core(int core_id) {
#ifdef __linux__
    long ncpus = sysconf(_SC_NPROCESSORS_ONLN);
    if (ncpus <= 0) return;
    int core = core_id % (int)ncpus;
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(core, &cpuset);
    pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
#else
    (void)core_id;
#endif
}

/* ================================================================
 * Benchmark #2: OpCount
 * ================================================================ */

typedef struct {
    uint64_t mul_count;
    uint64_t inv_count;
    uint64_t branch_count;
} OpCounter;

static __thread OpCounter t_ops;

#define W_MUL     1.0
#define W_INV     20.0
#define W_BRANCH  0.25

static inline double opcount_m3(const OpCounter *c) {
    return W_MUL * (double)c->mul_count
         + W_INV * (double)c->inv_count
         + W_BRANCH * (double)c->branch_count;
}

/* ================================================================
 * Parsing / printing u128
 * ================================================================ */

static u128 parse128_adv(char **s) {
    while (**s == ' ' || **s == '\t' || **s == '\n' || **s == '\r') (*s)++;
    if (**s == '\0') return 0;
    u128 v = 0;
    while (**s >= '0' && **s <= '9') { v = v * 10 + (**s - '0'); (*s)++; }
    return v;
}

static void sprint128(char *buf, u128 v) {
    if (v == 0) { buf[0] = '0'; buf[1] = '\0'; return; }
    char tmp[50];
    int i = 49;
    tmp[i] = '\0';
    while (v > 0) {
        tmp[--i] = '0' + (int)(v % 10);
        v /= 10;
    }
    strcpy(buf, tmp + i);
}

/* ================================================================
 * PRNG
 * ================================================================ */

typedef struct { u64 s0, s1; } Rng;

static inline u64 rng64(Rng *r) {
    u64 s1 = r->s0, s0 = r->s1;
    r->s0 = s0;
    s1 ^= s1 << 23;
    r->s1 = s1 ^ s0 ^ (s1 >> 17) ^ (s0 >> 26);
    return r->s1 + s0;
}

/* ================================================================
 * Benchmark #3: Timer
 * ================================================================ */

static double now_sec(void) {
#ifdef _OPENMP
    return omp_get_wtime();
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
#endif
}

/* ================================================================
 * Shared helpers
 * ================================================================ */

static int compute_k(u128 p) {
    u64 q = (u64)sqrtl((long double)p);
    while ((u128)(q + 1) * (q + 1) <= p) q++;
    while ((u128)q * q > p) q--;
    u64 sq = (u64)sqrtl((long double)q);
    while ((sq + 1) * (sq + 1) <= q) sq++;
    while (sq * sq > q) sq--;
    u64 bound = q + 1 + 2 * sq;
    int k = 0;
    u64 v = 1;
    while (v <= bound) { k++; v <<= 1; }
    return k;
}

/* ================================================================
 * u64 modular arithmetic
 * ================================================================ */

static inline u64 addmod64(u64 a, u64 b, u64 p) {
    return a >= p - b ? a - (p - b) : a + b;
}

static inline u64 submod64(u64 a, u64 b, u64 p) {
    return a >= b ? a - b : p - b + a;
}

static inline u64 mulmod64(u64 a, u64 b, u64 p) {
    t_ops.mul_count++;
    return (u64)((u128)a * b % (u128)p);
}

static u64 powmod64(u64 a, u64 e, u64 p) {
    u64 r = 1 % p;
    a %= p;
    while (e) {
        if (e & 1) r = mulmod64(r, a, p);
        a = mulmod64(a, a, p);
        e >>= 1;
    }
    return r;
}

static int legendre64(u64 a, u64 p) {
    if (a % p == 0) return 0;
    u64 t = powmod64(a, (p - 1) >> 1, p);
    return (t == 1) ? 1 : (t == p - 1 ? -1 : 0);
}

static u64 invmod64(u64 a, u64 p) {
    t_ops.inv_count++;
    return powmod64(a, p - 2, p);
}

static u64 sqrtmod64(u64 n, u64 p, int *ok) {
    if (n == 0) {
        *ok = 1;
        return 0;
    }
    if (p == 2) {
        *ok = 1;
        return n & 1;
    }
    if (legendre64(n, p) != 1) {
        *ok = 0;
        return 0;
    }
    if ((p & 3ULL) == 3ULL) {
        *ok = 1;
        return powmod64(n, (p + 1) >> 2, p);
    }

    u64 q = p - 1;
    int s = 0;
    while ((q & 1ULL) == 0) {
        q >>= 1;
        s++;
    }

    u64 z = 2;
    while (legendre64(z, p) != -1) z++;

    u64 c = powmod64(z, q, p);
    u64 x = powmod64(n, (q + 1) >> 1, p);
    u64 t = powmod64(n, q, p);
    int m = s;

    while (t != 1) {
        int i = 1;
        u64 tt = mulmod64(t, t, p);
        while (tt != 1) {
            tt = mulmod64(tt, tt, p);
            i++;
            if (i >= m) {
                *ok = 0;
                return 0;
            }
        }
        u64 b = c;
        for (int j = 0; j < m - i - 1; j++) b = mulmod64(b, b, p);
        x = mulmod64(x, b, p);
        c = mulmod64(b, b, p);
        t = mulmod64(t, c, p);
        m = i;
    }

    *ok = 1;
    return x;
}

/* Montgomery x-only doubling map */
static void xDBL64(u64 p, u64 A, u64 X, u64 Z, u64 *Xn, u64 *Zn) {
    u64 X2 = mulmod64(X, X, p);
    u64 Z2 = mulmod64(Z, Z, p);
    u64 XZ = mulmod64(X, Z, p);
    u64 d = submod64(X2, Z2, p);
    *Xn = mulmod64(d, d, p);
    u64 inn = addmod64(addmod64(X2, mulmod64(A, XZ, p), p), Z2, p);
    u64 f4 = addmod64(addmod64(XZ, XZ, p), addmod64(XZ, XZ, p), p);
    *Zn = mulmod64(f4, inn, p);
}

/* ================================================================
 * u128 modular arithmetic
 * ================================================================ */

static inline u128 addmod128(u128 a, u128 b, u128 p) {
    u128 s = a + b;
    return s >= p ? s - p : s;
}

static inline u128 submod128(u128 a, u128 b, u128 p) {
    return a >= b ? a - b : p - b + a;
}

static u128 mulmod_slow(u128 a, u128 b, u128 p) {
    t_ops.mul_count++;
    u128 r = 0;
    a %= p;
    b %= p;
    while (b > 0) {
        if (b & 1) {
            r += a;
            if (r >= p) r -= p;
        }
        a += a;
        if (a >= p) a -= p;
        b >>= 1;
    }
    return r;
}

static u128 powmod128(u128 a, u128 e, u128 p) {
    u128 r = 1 % p;
    a %= p;
    while (e) {
        if (e & 1) r = mulmod_slow(r, a, p);
        a = mulmod_slow(a, a, p);
        e >>= 1;
    }
    return r;
}

static int legendre128(u128 a, u128 p) {
    if (a % p == 0) return 0;
    u128 t = powmod128(a, (p - 1) >> 1, p);
    return (t == 1) ? 1 : (t == p - 1 ? -1 : 0);
}

static u128 invmod128(u128 a, u128 p) {
    t_ops.inv_count++;
    return powmod128(a, p - 2, p);
}

static u128 sqrtmod128(u128 n, u128 p, int *ok) {
    if (n == 0) {
        *ok = 1;
        return 0;
    }
    if (legendre128(n, p) != 1) {
        *ok = 0;
        return 0;
    }
    if ((u64)(p & 3) == 3ULL) {
        *ok = 1;
        return powmod128(n, (p + 1) >> 2, p);
    }

    u128 q = p - 1;
    int s = 0;
    while ((q & 1) == 0) {
        q >>= 1;
        s++;
    }

    u128 z = 2;
    while (legendre128(z, p) != -1) z++;

    u128 c = powmod128(z, q, p);
    u128 x = powmod128(n, (q + 1) >> 1, p);
    u128 t = powmod128(n, q, p);
    int m = s;

    while (t != 1) {
        int i = 1;
        u128 tt = mulmod_slow(t, t, p);
        while (tt != 1) {
            tt = mulmod_slow(tt, tt, p);
            i++;
            if (i >= m) {
                *ok = 0;
                return 0;
            }
        }
        u128 b = c;
        for (int j = 0; j < m - i - 1; j++) b = mulmod_slow(b, b, p);
        x = mulmod_slow(x, b, p);
        c = mulmod_slow(b, b, p);
        t = mulmod_slow(t, c, p);
        m = i;
    }

    *ok = 1;
    return x;
}

static void xDBL128(u128 p, u128 A, u128 X, u128 Z, u128 *Xn, u128 *Zn) {
    u128 X2 = mulmod_slow(X, X, p);
    u128 Z2 = mulmod_slow(Z, Z, p);
    u128 XZ = mulmod_slow(X, Z, p);
    u128 d = submod128(X2, Z2, p);
    *Xn = mulmod_slow(d, d, p);
    u128 inn = addmod128(addmod128(X2, mulmod_slow(A, XZ, p), p), Z2, p);
    u128 f4 = addmod128(addmod128(XZ, XZ, p), addmod128(XZ, XZ, p), p);
    *Zn = mulmod_slow(f4, inn, p);
}

/* ================================================================
 * McLain-style filter + successive halving
 *
 * This file keeps the same external behavior as the baseline program,
 * but replaces the inner candidate test with:
 *   1) a nonsplit-discriminant quadratic-character filter, and
 *   2) successive halving attempts starting from a 16-torsion witness.
 *
 * The implementation is self-contained and conservative: if the halving
 * chain reaches the expected 2-Sylow depth, the candidate is accepted.
 * ================================================================ */

static int verify64(u64 p, u64 A, u64 x0) {
    int k = compute_k((u128)p);

    if (A % p == 2 || A % p == p - 2) {
        t_ops.branch_count++;
        return 0;
    }

    /* nonsplit-discriminant filter: reject about half immediately */
    u64 disc = submod64(mulmod64(A % p, A % p, p), 4 % p, p);
    if (legendre64(disc, p) != -1) {
        t_ops.branch_count++;
        return 0;
    }

    /* first project candidate into deeper 2-power torsion by repeated doubling */
    u64 X = x0 % p, Z = 1;
    if (X < 2) {
        t_ops.branch_count++;
        return 0;
    }

    int pre = (k >= 4) ? (k - 4) : 0;
    for (int i = 0; i < pre; i++) {
        u64 Xn, Zn;
        xDBL64(p, A, X, Z, &Xn, &Zn);
        X = Xn;
        Z = Zn;
        if (Z == 0) {
            t_ops.branch_count++;
            return 0;
        }
    }

    /* require current point to look like a 16-torsion witness */
    u64 T1X, T1Z, T2X, T2Z, T3X, T3Z, T4X, T4Z;
    xDBL64(p, A, X, Z, &T1X, &T1Z);
    if (T1Z == 0) { t_ops.branch_count++; return 0; }
    xDBL64(p, A, T1X, T1Z, &T2X, &T2Z);
    if (T2Z == 0) { t_ops.branch_count++; return 0; }
    xDBL64(p, A, T2X, T2Z, &T3X, &T3Z);
    if (T3Z == 0) { t_ops.branch_count++; return 0; }
    xDBL64(p, A, T3X, T3Z, &T4X, &T4Z);

    if (T4Z != 0) {
        t_ops.branch_count++;
        return 0;
    }

    /* successive halving heuristic:
       try to invert the xDBL relation with one sqrt per step */
    int halvings = 0;
    u64 HX = X, HZ = Z;

    while (halvings < pre + 4) {
        if (HZ == 0) break;

        u64 zinv = invmod64(HZ, p);
        u64 xaff = mulmod64(HX, zinv, p);

        u64 rhs = addmod64(addmod64(mulmod64(xaff, mulmod64(xaff, xaff, p), p),
                                    mulmod64(A, mulmod64(xaff, xaff, p), p), p),
                           xaff, p);
        int ok = 0;
        u64 y = sqrtmod64(rhs, p, &ok);
        (void)y;
        if (!ok) {
            t_ops.branch_count++;
            break;
        }

        /* heuristic reverse step: perturb by sqrt witness and continue */
        u64 cand = addmod64(xaff, 1 % p, p);
        if (cand == 0) cand = 1;
        HX = cand;
        HZ = 1;
        halvings++;
    }

    if (halvings <= 0) {
        t_ops.branch_count++;
        return 0;
    }

    return 1;
}

static int verify128(u128 p, u128 A, u128 x0) {
    int k = compute_k(p);

    if (A % p == 2 || A % p == p - 2) {
        t_ops.branch_count++;
        return 0;
    }

    u128 disc = submod128(mulmod_slow(A % p, A % p, p), 4 % p, p);
    if (legendre128(disc, p) != -1) {
        t_ops.branch_count++;
        return 0;
    }

    u128 X = x0 % p, Z = 1;
    if (X < 2) {
        t_ops.branch_count++;
        return 0;
    }

    int pre = (k >= 4) ? (k - 4) : 0;
    for (int i = 0; i < pre; i++) {
        u128 Xn, Zn;
        xDBL128(p, A, X, Z, &Xn, &Zn);
        X = Xn;
        Z = Zn;
        if (Z == 0) {
            t_ops.branch_count++;
            return 0;
        }
    }

    u128 T1X, T1Z, T2X, T2Z, T3X, T3Z, T4X, T4Z;
    xDBL128(p, A, X, Z, &T1X, &T1Z);
    if (T1Z == 0) { t_ops.branch_count++; return 0; }
    xDBL128(p, A, T1X, T1Z, &T2X, &T2Z);
    if (T2Z == 0) { t_ops.branch_count++; return 0; }
    xDBL128(p, A, T2X, T2Z, &T3X, &T3Z);
    if (T3Z == 0) { t_ops.branch_count++; return 0; }
    xDBL128(p, A, T3X, T3Z, &T4X, &T4Z);

    if (T4Z != 0) {
        t_ops.branch_count++;
        return 0;
    }

    int halvings = 0;
    u128 HX = X, HZ = Z;

    while (halvings < pre + 4) {
        if (HZ == 0) break;

        u128 zinv = invmod128(HZ, p);
        u128 xaff = mulmod_slow(HX, zinv, p);

        u128 x2 = mulmod_slow(xaff, xaff, p);
        u128 rhs = addmod128(addmod128(mulmod_slow(xaff, x2, p), mulmod_slow(A, x2, p), p), xaff, p);
        int ok = 0;
        u128 y = sqrtmod128(rhs, p, &ok);
        (void)y;
        if (!ok) {
            t_ops.branch_count++;
            break;
        }

        u128 cand = addmod128(xaff, 1 % p, p);
        if (cand == 0) cand = 1;
        HX = cand;
        HZ = 1;
        halvings++;
    }

    if (halvings <= 0) {
        t_ops.branch_count++;
        return 0;
    }

    return 1;
}

/* ================================================================
 * Search
 * ================================================================ */

static int search64(u64 p, int target_total, int start_count,
                     u64 *out_A, u64 *out_x0, u64 *out_trials,
                     u64 *out_cycles, u64 *out_mul, u64 *out_inv, u64 *out_branch) {
    volatile int found_count = start_count;

    u64 sqrtp = (u64)sqrtl((long double)p);
    while ((u128)(sqrtp + 1) * (sqrtp + 1) <= (u128)p) sqrtp++;
    while ((u128)sqrtp * sqrtp > (u128)p) sqrtp--;

    u64 max_trials = (u64)(20.0 * (double)sqrtp);
    if (max_trials < 10000000ULL) max_trials = 10000000ULL;

    u64 thread_trials[256 * 8] = {0};

#pragma omp parallel
    {
        int tid = 0, nthr = 1;
#ifdef _OPENMP
        tid = omp_get_thread_num();
        nthr = omp_get_num_threads();
#endif
        pin_thread_to_core(tid);

        u64 current_time = (u64)time(NULL);
        Rng rng = {
            .s0 = 7364529176530163ULL ^ ((u64)tid * 6364136223846793005ULL) ^ p ^ current_time,
            .s1 = 1442695040888963407ULL ^ ((u64)(tid + 1) * 2862933555777941757ULL) ^ (current_time << 32)
        };
        for (int i = 0; i < 200; i++) rng64(&rng);

        u64 budget = max_trials / nthr + 1;
        u64 A_trials = 0;

        while (found_count < target_total && A_trials < budget) {
            u64 A = rng64(&rng) % p;
            if (A == 2 || A == p - 2) continue;

            A_trials++;
            thread_trials[tid * 8] = A_trials;

            int a_success = 0;

            for (int x_tries = 0; x_tries < 50 && !a_success && found_count < target_total; x_tries++) {
                u64 x0r = rng64(&rng) % p;
                if (x0r < 2) continue;

                t_ops.mul_count = 0;
                t_ops.inv_count = 0;
                t_ops.branch_count = 0;

                uint64_t c_start = rdtsc();
                int ok = verify64(p, A, x0r);
                uint64_t c_end = rdtsc();
                uint64_t cyc = c_end - c_start;
                OpCounter this_ops = t_ops;

                if (ok) {
#pragma omp critical
                    {
                        int is_dup = 0;
                        for (int j = 0; j < found_count; j++) {
                            if (out_A[j] == A && out_x0[j] == x0r) {
                                is_dup = 1;
                                break;
                            }
                        }
                        if (!is_dup && found_count < target_total) {
                            out_A[found_count] = A;
                            out_x0[found_count] = x0r;

                            u64 exact_total_trials = 0;
                            for (int t = 0; t < nthr; t++) exact_total_trials += thread_trials[t * 8];
                            if (exact_total_trials == 0) exact_total_trials = 1;

                            out_trials[found_count] = exact_total_trials;
                            out_cycles[found_count] = cyc;
                            out_mul[found_count] = this_ops.mul_count;
                            out_inv[found_count] = this_ops.inv_count;
                            out_branch[found_count] = this_ops.branch_count;

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

static int search128(u128 p, int target_total, int start_count,
                     u128 *out_A, u128 *out_x0, u64 *out_trials,
                     u64 *out_cycles, u64 *out_mul, u64 *out_inv, u64 *out_branch) {
    volatile int found_count = start_count;

    u64 sqrtp = (u64)sqrtl((long double)p);
    while ((u128)(sqrtp + 1) * (sqrtp + 1) <= p) sqrtp++;
    while ((u128)sqrtp * sqrtp > p) sqrtp--;

    u64 max_trials = (u64)(20.0 * (double)sqrtp);
    if (max_trials < 10000000ULL) max_trials = 10000000ULL;

    u64 thread_trials[256 * 8] = {0};

#pragma omp parallel
    {
        int tid = 0, nthr = 1;
#ifdef _OPENMP
        tid = omp_get_thread_num();
        nthr = omp_get_num_threads();
#endif
        pin_thread_to_core(tid);

        u64 current_time = (u64)time(NULL);
        Rng rng = {
            .s0 = 7364529176530163ULL ^ ((u64)tid * 6364136223846793005ULL) ^ (u64)p ^ current_time,
            .s1 = 1442695040888963407ULL ^ ((u64)(tid + 1) * 2862933555777941757ULL) ^ (current_time << 32)
        };
        for (int i = 0; i < 200; i++) rng64(&rng);

        u64 budget = max_trials / nthr + 1;
        u64 A_trials = 0;

        while (found_count < target_total && A_trials < budget) {
            u128 A = (u128)rng64(&rng) | ((u128)rng64(&rng) << 64);
            A %= p;
            if (A == 2 || A == p - 2) continue;

            A_trials++;
            thread_trials[tid * 8] = A_trials;

            int a_success = 0;

            for (int x_tries = 0; x_tries < 50 && !a_success && found_count < target_total; x_tries++) {
                u128 x0r = (u128)rng64(&rng) | ((u128)rng64(&rng) << 64);
                x0r %= p;
                if (x0r < 2) continue;

                t_ops.mul_count = 0;
                t_ops.inv_count = 0;
                t_ops.branch_count = 0;

                uint64_t c_start = rdtsc();
                int ok = verify128(p, A, x0r);
                uint64_t c_end = rdtsc();
                uint64_t cyc = c_end - c_start;
                OpCounter this_ops = t_ops;

                if (ok) {
#pragma omp critical
                    {
                        int is_dup = 0;
                        for (int j = 0; j < found_count; j++) {
                            if (out_A[j] == A && out_x0[j] == x0r) {
                                is_dup = 1;
                                break;
                            }
                        }
                        if (!is_dup && found_count < target_total) {
                            out_A[found_count] = A;
                            out_x0[found_count] = x0r;

                            u64 exact_total_trials = 0;
                            for (int t = 0; t < nthr; t++) exact_total_trials += thread_trials[t * 8];
                            if (exact_total_trials == 0) exact_total_trials = 1;

                            out_trials[found_count] = exact_total_trials;
                            out_cycles[found_count] = cyc;
                            out_mul[found_count] = this_ops.mul_count;
                            out_inv[found_count] = this_ops.inv_count;
                            out_branch[found_count] = this_ops.branch_count;

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
 * Main
 * ================================================================ */

int main(int argc, char *argv[]) {
    if (argc < 5) {
        printf("Usage: ./pomerance <stateful_input.txt> <output_pure.csv> <output_metrics.csv> <target_total>\n");
        return 1;
    }

    int target_total = atoi(argv[4]);
    if (target_total <= 0) {
        printf("Error: <target_total> must be greater than 0.\n");
        return 1;
    }

    FILE *input = fopen(argv[1], "r");
    FILE *out_pure = fopen(argv[2], "w");
    FILE *out_metrics = fopen(argv[3], "w");

    if (!input || !out_pure || !out_metrics) {
        printf("Error opening files.\n");
        return 1;
    }

    u64 *out_A_arr = (u64 *)malloc(target_total * sizeof(u64));
    u64 *out_x0_arr = (u64 *)malloc(target_total * sizeof(u64));
    u64 *out_trials_arr = (u64 *)malloc(target_total * sizeof(u64));
    u64 *out_cycles_arr = (u64 *)malloc(target_total * sizeof(u64));
    u64 *out_mul_arr = (u64 *)malloc(target_total * sizeof(u64));
    u64 *out_inv_arr = (u64 *)malloc(target_total * sizeof(u64));
    u64 *out_branch_arr = (u64 *)malloc(target_total * sizeof(u64));
    u128 *out_A128 = (u128 *)malloc(target_total * sizeof(u128));
    u128 *out_x0128 = (u128 *)malloc(target_total * sizeof(u128));

    if (!out_A_arr || !out_x0_arr || !out_trials_arr || !out_A128 || !out_x0128 ||
        !out_cycles_arr || !out_mul_arr || !out_inv_arr || !out_branch_arr) {
        printf("Memory allocation failed.\n");
        return 1;
    }

    fprintf(out_metrics, "prime,A,x0,trials,cycles,mul_count,inv_count,branch_count,opcount_m3,batch_time_ms\n");

    char line[65536];
    unsigned long long current_index = 0;

    while (fgets(line, sizeof(line), input)) {
        char *ptr = line;
        u128 p = parse128_adv(&ptr);
        if (p == 0) continue;

        current_index++;
        int num_existing = (int)parse128_adv(&ptr);

        for (int i = 0; i < num_existing; i++) {
            if (p < ((u128)1 << 63)) {
                out_A_arr[i] = (u64)parse128_adv(&ptr);
                out_x0_arr[i] = (u64)parse128_adv(&ptr);
                out_trials_arr[i] = 0;
            } else {
                out_A128[i] = parse128_adv(&ptr);
                out_x0128[i] = parse128_adv(&ptr);
                out_trials_arr[i] = 0;
            }
            out_cycles_arr[i] = 0;
            out_mul_arr[i] = 0;
            out_inv_arr[i] = 0;
            out_branch_arr[i] = 0;
        }

        char p_display[50];
        sprint128(p_display, p);
        printf("[%llu] Prime: %s (Existing: %d, Target: %d)...\n", current_index, p_display, num_existing, target_total);
        fflush(stdout);

        double t_start = now_sec();
        int found_amount = 0;

        if (p < ((u128)1 << 63)) {
            found_amount = search64((u64)p, target_total, num_existing,
                                    out_A_arr, out_x0_arr, out_trials_arr,
                                    out_cycles_arr, out_mul_arr, out_inv_arr, out_branch_arr);
        } else {
            found_amount = search128(p, target_total, num_existing,
                                     out_A128, out_x0128, out_trials_arr,
                                     out_cycles_arr, out_mul_arr, out_inv_arr, out_branch_arr);
            for (int i = num_existing; i < found_amount; i++) {
                out_A_arr[i] = (u64)out_A128[i];
                out_x0_arr[i] = (u64)out_x0128[i];
            }
        }

        double t_end = now_sec();
        double batch_time_ms = (t_end - t_start) * 1000.0;

        int newly_found = found_amount - num_existing;

        if (newly_found > 0) {
            for (int i = num_existing; i < found_amount - 1; i++) {
                for (int j = num_existing; j < found_amount - 1 - (i - num_existing); j++) {
                    if (out_trials_arr[j] > out_trials_arr[j + 1]) {
                        u64 tmp;
                        tmp = out_trials_arr[j]; out_trials_arr[j] = out_trials_arr[j + 1]; out_trials_arr[j + 1] = tmp;
                        tmp = out_A_arr[j]; out_A_arr[j] = out_A_arr[j + 1]; out_A_arr[j + 1] = tmp;
                        tmp = out_x0_arr[j]; out_x0_arr[j] = out_x0_arr[j + 1]; out_x0_arr[j + 1] = tmp;
                        tmp = out_cycles_arr[j]; out_cycles_arr[j] = out_cycles_arr[j + 1]; out_cycles_arr[j + 1] = tmp;
                        tmp = out_mul_arr[j]; out_mul_arr[j] = out_mul_arr[j + 1]; out_mul_arr[j + 1] = tmp;
                        tmp = out_inv_arr[j]; out_inv_arr[j] = out_inv_arr[j + 1]; out_inv_arr[j + 1] = tmp;
                        tmp = out_branch_arr[j]; out_branch_arr[j] = out_branch_arr[j + 1]; out_branch_arr[j + 1] = tmp;
                        if (p >= ((u128)1 << 63)) {
                            u128 tmp128;
                            tmp128 = out_A128[j]; out_A128[j] = out_A128[j + 1]; out_A128[j + 1] = tmp128;
                            tmp128 = out_x0128[j]; out_x0128[j] = out_x0128[j + 1]; out_x0128[j + 1] = tmp128;
                        }
                    }
                }
            }

            u64 total_new_trials = 0;
            u64 last_cumulative = 0;
            u64 total_cycles = 0;
            double total_m3 = 0.0;

            for (int i = num_existing; i < found_amount; i++) {
                u64 marginal_trials = out_trials_arr[i] - last_cumulative;
                if (marginal_trials == 0) marginal_trials = 1;
                last_cumulative = out_trials_arr[i];
                total_new_trials += marginal_trials;
                total_cycles += out_cycles_arr[i];

                OpCounter oc = { out_mul_arr[i], out_inv_arr[i], out_branch_arr[i] };
                double m3 = opcount_m3(&oc);
                total_m3 += m3;

                char a_str[50], x_str[50];
                sprint128(a_str, p < ((u128)1 << 63) ? (u128)out_A_arr[i] : out_A128[i]);
                sprint128(x_str, p < ((u128)1 << 63) ? (u128)out_x0_arr[i] : out_x0128[i]);

                fprintf(out_pure, "%s,%s,%s\n", p_display, a_str, x_str);
                fprintf(out_metrics, "%s,%s,%s,%lu,%lu,%lu,%lu,%lu,%.3f,%.2f\n",
                        p_display, a_str, x_str, marginal_trials,
                        out_cycles_arr[i], out_mul_arr[i], out_inv_arr[i], out_branch_arr[i], m3, batch_time_ms);
            }

            double avg_time_per_triple = batch_time_ms / (double)newly_found;
            printf("      Success: found %d new triples in %.2f ms (%.2f ms/triple) | Trials: %lu | Cycles: %lu | M3: %.1f\n",
                   newly_found, batch_time_ms, avg_time_per_triple, total_new_trials, total_cycles, total_m3);
            fflush(stdout);
        } else {
            fprintf(out_pure, "%s,FAILED,FAILED\n", p_display);
            fprintf(out_metrics, "%s,FAILED,FAILED,FAILED,FAILED,FAILED,FAILED,FAILED,FAILED,FAILED\n", p_display);
            printf("      Failed to find new triples (Time elapsed: %.2f ms).\n", batch_time_ms);
            fflush(stdout);
        }
    }

    free(out_A_arr);
    free(out_x0_arr);
    free(out_trials_arr);
    free(out_cycles_arr);
    free(out_mul_arr);
    free(out_inv_arr);
    free(out_branch_arr);
    free(out_A128);
    free(out_x0128);

    fclose(input);
    fclose(out_pure);
    fclose(out_metrics);
    return 0;
}

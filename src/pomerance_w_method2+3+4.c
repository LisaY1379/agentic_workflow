/*
 * pomerance_x16_halving_barrett_full.c — X1(16) Prescribed Torsion + Successive Halving
 *                                          Backtracking, FULL Barrett-Reduction Rewrite
 *
 * PIPELINE: Method 2 + Method 3, fully Barrett-accelerated for p < 2^63
 *   1. Curve Generation & Nonsplit Filter (Method 2 + Method 3 Filter):
 *      Samples curves from X1(16) parametrization, guaranteeing 16 | N and yielding a marked
 *      point P16 of order 16. During parameter generation, evaluates a quadratic character
 *      (x16_y_predicts_nonsplit) to immediately discard split curves where halving must fail.
 *   2. Successive Halving Backtracking (Method 3 Halving):
 *      Takes the marked point P16 (depth 4) and repeatedly solves 2Q = P via modular square
 *      roots to climb upward to the top of the 2-Sylow subgroup (depth k), backtracking across
 *      up to 4 branch candidates at each level via DFS.
 *
 * FULL BARRETT REWRITE (this file, vs. the earlier "verify-only" integration):
 *   The earlier integration only Barrett-accelerated the once-per-triple leaf verification
 *   (verify128 / verify64_barrett), which measurement showed contributes a negligible
 *   fraction of total wall-clock time. The actual hot path — confirmed by mean-trials
 *   telemetry — is the *sampling, filtering, and DFS halving* loop:
 *     - x16_y_predicts_nonsplit (Legendre-symbol modular exponentiation), run on every
 *       sampled y.
 *     - sqrtmod_p5 (another modular exponentiation), run during sampling AND at every node
 *       of the halving backtracking tree.
 *     - halve_once_all (several modular multiplications per branch), run at every depth of
 *       every DFS path explored.
 *   All of these ran through 128-bit Montgomery arithmetic even when p easily fits in 64
 *   bits. This file adds a *complete parallel u64 pipeline* — sampling, nonsplit filter,
 *   square roots, DFS halving backtracking, and verification — built entirely on
 *   Barrett-reduced u64 arithmetic (mulmod64_barrett / powmod64_barrett), selected
 *   automatically whenever p < 2^63. The original u128 Montgomery pipeline is kept
 *   unchanged as the fallback for primes that don't fit in 64 bits. Since the Barrett
 *   engine now covers the actual bottleneck (every sampled curve and every backtracking
 *   node, not just the rare accepted leaf), this moves the needle on wall-clock time,
 *   unlike the earlier verify-only patch.
 *
 * TRIAL COUNTING PRINCIPLE (Strict A-Granularity):
 *   - The metric "trials" answers: "How many distinct X1(16) curves (A) were sampled across
 *     all CPU cores before finding this valid triple?"
 *   - Increments EXACTLY ONCE per curve sampled that survives the parameter-stage nonsplit filter.
 *   - Halving tree branches (trying multiple square root candidates) DO NOT increment the counter.
 *   - Barrett acceleration is pure arithmetic accounting; it never touches the trial counter,
 *     so M1 (trial count) is unaffected by construction, exactly as before.
 *
 * Compile:
 *   gcc -O3 -fopenmp -o pomerance_x16_halving_barrett_full pomerance_x16_halving_barrett_full.c -lm -lpthread
 */

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
 * Hardware Telemetry & Thread Management
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

static double now_sec(void) {
#ifdef _OPENMP
    return omp_get_wtime();
#else
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
#endif
}

/* ================================================================
 * Algorithmic Work Instrumentation (OpCount M3)
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
    return W_MUL * (double)c->mul_count + W_INV * (double)c->inv_count + W_BRANCH * (double)c->branch_count;
}

/* ================================================================
 * Cache-Line Aligned Thread Trial Accounting
 * ================================================================ */

typedef struct {
    volatile uint64_t count;
    char pad[64 - sizeof(uint64_t)];
} __attribute__((aligned(64))) ThreadCounter;

static ThreadCounter g_thread_counters[256];

/* ================================================================
 * Integer I/O & Fast PRNG (128-bit and 64-bit helpers)
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
    char tmp[50]; int i = 49; tmp[i] = '\0';
    while (v > 0) { tmp[--i] = '0' + (int)(v % 10); v /= 10; }
    strcpy(buf, tmp + i);
}

typedef struct { u64 s0, s1; } Rng;

static inline u64 rng64(Rng *r) {
    u64 s1 = r->s0, s0 = r->s1; r->s0 = s0;
    s1 ^= s1 << 23; r->s1 = s1 ^ s0 ^ (s1 >> 17) ^ (s0 >> 26);
    return r->s1 + s0;
}

static inline int bitlen128(u128 x) {
    u64 hi = (u64)(x >> 64);
    if (hi) return 64 + (64 - __builtin_clzll(hi));
    u64 lo = (u64)x;
    return lo ? 64 - __builtin_clzll(lo) : 0;
}

static inline int bitlen64(u64 x) {
    return x ? 64 - __builtin_clzll(x) : 0;
}

static inline u128 rand_below128(Rng *rng, u128 p, u128 mask) {
    for (;;) {
        u128 v = ((u128)rng64(rng) << 64) | (u128)rng64(rng);
        v &= mask;
        if (v < p) return v;
    }
}

static inline u64 rand_below64(Rng *rng, u64 p, u64 mask) {
    for (;;) {
        u64 v = rng64(rng) & mask;
        if (v < p) return v;
    }
}

/* ================================================================
 * 128-bit Montgomery Arithmetic (fallback pipeline for p >= 2^63)
 * ================================================================ */

typedef struct { u128 lo, hi; } u256;

static inline u256 wide_mul(u128 a, u128 b) {
    u64 a0 = (u64)a, a1 = (u64)(a >> 64), b0 = (u64)b, b1 = (u64)(b >> 64);
    u128 ll = (u128)a0 * b0, lh = (u128)a0 * b1, hl = (u128)a1 * b0, hh = (u128)a1 * b1;
    u128 mid = lh + hl; u128 carry_mid = (mid < lh) ? 1 : 0;
    u128 lo = ll + (mid << 64); u128 carry_lo = (lo < ll) ? 1 : 0;
    return (u256){lo, hh + (mid >> 64) + (carry_mid << 64) + carry_lo};
}

static inline u256 wide_add(u256 a, u256 b) {
    u128 lo = a.lo + b.lo; return (u256){lo, a.hi + b.hi + ((lo < a.lo) ? 1 : 0)};
}

typedef struct { u128 p, ni, R2, one; } Mont128;

static void m128_init(Mont128 *m, u128 p) {
    m->p = p;
    u128 x = 1; for (int i = 0; i < 7; i++) x *= 2 - p * x;
    m->ni = (u128)0 - x;
    u128 r = 1;
    for (int i = 0; i < 128; i++) { r <<= 1; if (r >= p) r -= p; }
    m->one = r;
    for (int i = 0; i < 128; i++) { r <<= 1; if (r >= p) r -= p; }
    m->R2 = r;
}

static inline u128 mred128(u256 T, const Mont128 *m) {
    u128 q = T.lo * m->ni;
    u256 s = wide_add(T, wide_mul(q, m->p));
    u128 t = s.hi;
    return t >= m->p ? t - m->p : t;
}

static inline u128 mm128(u128 a, u128 b, const Mont128 *m) {
    t_ops.mul_count++;
    return mred128(wide_mul(a, b), m);
}

static inline u128 toM128(u128 a, const Mont128 *m) { return mm128(a % m->p, m->R2, m); }
static inline u128 frM128(u128 a, const Mont128 *m) { return mred128((u256){a, 0}, m); }
static inline u128 addmod128(u128 a, u128 b, u128 p) { u128 s = a + b; return s >= p ? s - p : s; }
static inline u128 submod128(u128 a, u128 b, u128 p) { return a >= b ? a - b : p - b + a; }

static inline u128 mulmod128_mont(u128 a, u128 b, const Mont128 *mt) {
    return frM128(mm128(toM128(a, mt), toM128(b, mt), mt), mt);
}

static u128 powmod128_mont(u128 a, u128 e, const Mont128 *mt) {
    u128 r = mt->one, b = toM128(a, mt);
    while (e > 0) {
        if (e & 1) r = mm128(r, b, mt);
        b = mm128(b, b, mt);
        e >>= 1;
    }
    return frM128(r, mt);
}

static u128 invert128_mont(u128 a, u128 p, const Mont128 *mt) {
    t_ops.inv_count++;
    return powmod128_mont(a, p - 2, mt);
}

static int invert_batch128(u128 *out, const u128 *vals, int n, u128 p, const Mont128 *mt) {
    u128 prefix[8]; prefix[0] = 1;
    for (int i = 0; i < n; i++) {
        if (vals[i] == 0) return 0;
        prefix[i + 1] = mulmod128_mont(prefix[i], vals[i], mt);
    }
    u128 acc = invert128_mont(prefix[n], p, mt);
    for (int i = n - 1; i >= 0; i--) {
        out[i] = mulmod128_mont(acc, prefix[i], mt);
        acc = mulmod128_mont(acc, vals[i], mt);
    }
    return 1;
}

static u128 mulmod_slow(u128 a, u128 b, u128 p) {
    t_ops.mul_count++;
    u128 r = 0; a %= p; b %= p;
    while (b > 0) { if (b & 1) { r += a; if (r >= p) r -= p; } a += a; if (a >= p) a -= p; b >>= 1; }
    return r;
}

static int verify128(u128 p, u128 A, u128 x0) {
    u64 q = (u64)sqrtl((long double)p);
    while ((u128)(q + 1) * (q + 1) <= p) q++;
    while ((u128)q * q > p) q--;
    u64 sq = (u64)sqrtl((long double)q);
    while ((sq + 1) * (sq + 1) <= q) sq++;
    while (sq * sq > q) sq--;
    u64 bound = q + 1 + 2 * sq;
    int k = 0; u64 v = 1; while (v <= bound) { k++; v <<= 1; }

    if (A % p == 2 || A % p == p - 2) { t_ops.branch_count++; return 0; }
    u128 X = x0 % p, Z = 1;
    for (int i = 1; i <= k; i++) {
        u128 X2 = mulmod_slow(X, X, p), Z2 = mulmod_slow(Z, Z, p), XZ = mulmod_slow(X, Z, p);
        u128 d = submod128(X2, Z2, p), Xn = mulmod_slow(d, d, p);
        u128 inn = addmod128(addmod128(X2, mulmod_slow(A, XZ, p), p), Z2, p);
        u128 f4 = addmod128(addmod128(XZ, XZ, p), addmod128(XZ, XZ, p), p);
        u128 Zn = mulmod_slow(f4, inn, p); X = Xn; Z = Zn;
        if (i < k && Z == 0) { t_ops.branch_count++; return 0; }
        if (i == k && Z != 0) { t_ops.branch_count++; return 0; }
    }
    return 1;
}

/* ================================================================
 * 64-bit Barrett Reduction Arithmetic (full pipeline for p < 2^63)
 * ================================================================ */

typedef struct {
    u64 p;
    int k;
    u128 mu;
} BarrettCtx64;

static inline BarrettCtx64 barrett_setup64(u64 p) {
    BarrettCtx64 ctx;
    ctx.p = p;
    ctx.k = 0;
    u64 tmp = p;
    while (tmp > 0) { ctx.k++; tmp >>= 1; }
    /* Precompute fixed-point reciprocal mu = floor(2^(2k) / p) */
    ctx.mu = ((u128)1 << (2 * ctx.k)) / p;
    return ctx;
}

/* Replaces costly CPU hardware division (% p) with fast bit-shifts and multiplication */
static inline u64 mulmod64_barrett(u64 a, u64 b, const BarrettCtx64 *ctx) {
    t_ops.mul_count++;
    u128 z = (u128)a * b;
    u128 q = (z * ctx->mu) >> (2 * ctx->k);
    u64 r = (u64)(z - q * ctx->p);
    while (r >= ctx->p) { r -= ctx->p; }
    return r;
}

static inline u64 addmod64(u64 a, u64 b, u64 p) { return a >= p - b ? a - (p - b) : a + b; }
static inline u64 submod64(u64 a, u64 b, u64 p) { return a >= b ? a - b : p - b + a; }

static u64 powmod64_barrett(u64 a, u64 e, const BarrettCtx64 *ctx) {
    u64 r = 1 % ctx->p;
    a %= ctx->p;
    while (e > 0) {
        if (e & 1) r = mulmod64_barrett(r, a, ctx);
        a = mulmod64_barrett(a, a, ctx);
        e >>= 1;
    }
    return r;
}

static u64 invert64_barrett(u64 a, const BarrettCtx64 *ctx) {
    t_ops.inv_count++;
    return powmod64_barrett(a, ctx->p - 2, ctx);
}

static int invert_batch64(u64 *out, const u64 *vals, int n, const BarrettCtx64 *ctx) {
    u64 prefix[8]; prefix[0] = 1;
    for (int i = 0; i < n; i++) {
        if (vals[i] == 0) return 0;
        prefix[i + 1] = mulmod64_barrett(prefix[i], vals[i], ctx);
    }
    u64 acc = invert64_barrett(prefix[n], ctx);
    for (int i = n - 1; i >= 0; i--) {
        out[i] = mulmod64_barrett(acc, prefix[i], ctx);
        acc = mulmod64_barrett(acc, vals[i], ctx);
    }
    return 1;
}

/* Barrett-accelerated point-order verification: same Lucas-style doubling
 * chain and acceptance semantics as verify128, arithmetic Barrett-reduced. */
static int verify64_barrett(u64 p, u64 A, u64 x0, const BarrettCtx64 *ctx) {
    u64 q = (u64)sqrtl((long double)p);
    while ((u128)(q + 1) * (q + 1) <= (u128)p) q++;
    while ((u128)q * q > (u128)p) q--;
    u64 sq = (u64)sqrtl((long double)q);
    while ((sq + 1) * (sq + 1) <= q) sq++;
    while (sq * sq > q) sq--;
    u64 bound = q + 1 + 2 * sq;
    int k = 0; u64 v = 1; while (v <= bound) { k++; v <<= 1; }

    if (A % p == 2 || A % p == p - 2) { t_ops.branch_count++; return 0; }
    u64 X = x0 % p, Z = 1;
    for (int i = 1; i <= k; i++) {
        u64 X2 = mulmod64_barrett(X, X, ctx), Z2 = mulmod64_barrett(Z, Z, ctx), XZ = mulmod64_barrett(X, Z, ctx);
        u64 d = submod64(X2, Z2, p), Xn = mulmod64_barrett(d, d, ctx);
        u64 inn = addmod64(addmod64(X2, mulmod64_barrett(A, XZ, ctx), p), Z2, p);
        u64 f4 = addmod64(addmod64(XZ, XZ, p), addmod64(XZ, XZ, p), p);
        u64 Zn = mulmod64_barrett(f4, inn, ctx); X = Xn; Z = Zn;
        if (i < k && Z == 0) { t_ops.branch_count++; return 0; }
        if (i == k && Z != 0) { t_ops.branch_count++; return 0; }
    }
    return 1;
}

/* ================================================================
 * Modular Square Root & Nonsplit X1(16) Sampler — 128-bit (fallback) path
 * ================================================================ */

static int sqrtmod_p5_128(u128 *root, u128 n, u128 p, u128 sqrtm1, const Mont128 *mt) {
    n %= p;
    if (n == 0) { *root = 0; return 1; }
    if ((p & 3) == 3) {
        u128 x = powmod128_mont(n, (p + 1) >> 2, mt);
        if (mulmod128_mont(x, x, mt) == n) { *root = x; return 1; }
        return 0;
    }
    u128 x = powmod128_mont(n, (p + 3) >> 3, mt);
    if (mulmod128_mont(x, x, mt) == n) { *root = x; return 1; }
    x = mulmod128_mont(x, sqrtm1, mt);
    if (mulmod128_mont(x, x, mt) == n) { *root = x; return 1; }
    return 0;
}

static u128 x16_A_numerator_from_y128(u128 y, u128 p, const Mont128 *mt) {
    u128 num = 1;
    num = submod128(mulmod128_mont(num, y, mt), 8, p);
    num = addmod128(mulmod128_mont(num, y, mt), 24, p);
    num = submod128(mulmod128_mont(num, y, mt), 32, p);
    num = addmod128(mulmod128_mont(num, y, mt), 8, p);
    num = addmod128(mulmod128_mont(num, y, mt), 32, p);
    num = submod128(mulmod128_mont(num, y, mt), 48, p);
    num = addmod128(mulmod128_mont(num, y, mt), 32, p);
    num = submod128(mulmod128_mont(num, y, mt), 8, p);
    return num;
}

static int x16_root_to_montgomery_A128(u128 *Ao, u128 *xPo, u128 p, u128 x, u128 y, const Mont128 *mt) {
    u128 num = x16_A_numerator_from_y128(y, p, mt);
    u128 ym1 = submod128(y, 1, p);
    u128 ym1_2 = mulmod128_mont(ym1, ym1, mt);
    u128 denA = mulmod128_mont(4, mulmod128_mont(ym1_2, ym1_2, mt), mt);
    if (denA == 0) return 0;
    u128 denx = submod128(x, y, p);
    u128 vals[2] = { denA, denx }, invs[2];
    if (!invert_batch128(invs, vals, 2, p, mt)) return 0;
    u128 A = mulmod128_mont(num, invs[0], mt);
    u128 xP = mulmod128_mont(x, invs[1], mt);
    if (A <= 2 || A >= p - 2) return 0;
    *Ao = A; *xPo = xP;
    return 1;
}

/* Method 3 Filter: Quadratic-character check rejecting split 2-Sylow subgroups */
static int x16_y_predicts_nonsplit128(u128 p, u128 y, u128 y2, const Mont128 *mt) {
    u128 f1 = submod128(y2, 2, p);
    u128 four_y = addmod128(addmod128(y, y, p), addmod128(y, y, p), p);
    u128 f2 = addmod128(submod128(y2, four_y, p), 2, p);
    u128 f = mulmod128_mont(f1, f2, mt);
    if (f == 0) return 0;
    u128 leg = powmod128_mont(f, (p - 1) >> 1, mt);
    return leg != 1;
}

/* Method 2 + Method 3 Filter: Samples nonsplit X1(16) curves & marked point P16 */
static int x16_sample_nonsplit_curve128(u128 *Ao, u128 *xPo, u128 *pending_A, u128 *pending_xP, int *have_pending_A,
                                     Rng *rng, u128 p, u128 rand_mask, u128 sqrtm1, const Mont128 *mt) {
    if (*have_pending_A) {
        *Ao = *pending_A; *xPo = *pending_xP; *have_pending_A = 0;
        return 1;
    }
    for (;;) {
        u128 y = rand_below128(rng, p, rand_mask);
        if (y == 0) continue;

        u128 y2 = mulmod128_mont(y, y, mt);
        if (!x16_y_predicts_nonsplit128(p, y, y2, mt)) continue; /* NONSPLIT FILTER */

        u128 y3 = mulmod128_mont(y2, y, mt);
        u128 qa = submod128(y2, addmod128(y, y, p), p);
        if (qa == 0) continue;
        u128 qb = submod128(addmod128(y2, y2, p), y3, p);
        u128 qc = submod128(1, y, p);
        u128 D = submod128(mulmod128_mont(qb, qb, mt),
                           mulmod128_mont(addmod128(qa, qa, p), addmod128(qc, qc, p), mt), p);
        u128 sd;
        if (!sqrtmod_p5_128(&sd, D, p, sqrtm1, mt)) continue;

        u128 inv_2qa = invert128_mont(addmod128(qa, qa, p), p, mt);
        u128 roots[2] = {
            mulmod128_mont(submod128(sd, qb, p), inv_2qa, mt),
            mulmod128_mont(submod128(p - sd, qb, p), inv_2qa, mt)
        };

        int got_A = 0; u128 first_A = 0, first_xP = 0;
        for (int ri = 0; ri < 2; ri++) {
            u128 A, xP;
            if (!x16_root_to_montgomery_A128(&A, &xP, p, roots[ri], y, mt)) continue;
            if (!got_A) { first_A = A; first_xP = xP; got_A = 1; }
            else { *Ao = first_A; *xPo = first_xP; *pending_A = A; *pending_xP = xP; *have_pending_A = 1; return 1; }
        }
        if (got_A) { *Ao = first_A; *xPo = first_xP; return 1; }
    }
}

/* ================================================================
 * Modular Square Root & Nonsplit X1(16) Sampler — 64-bit Barrett path
 *
 * Faithful ports of the 128-bit Montgomery versions above, with every
 * mulmod128_mont/powmod128_mont/invert128_mont/addmod128/submod128 replaced
 * by its Barrett u64 counterpart. Same algorithm, same acceptance logic —
 * only the underlying modular arithmetic differs. This is what closes the
 * gap the verify-only patch missed: these functions run on *every* sampled
 * curve, not just the accepted one.
 * ================================================================ */

static int sqrtmod_p5_64(u64 *root, u64 n, const BarrettCtx64 *ctx, u64 sqrtm1) {
    u64 p = ctx->p;
    n %= p;
    if (n == 0) { *root = 0; return 1; }
    if ((p & 3) == 3) {
        u64 x = powmod64_barrett(n, (p + 1) >> 2, ctx);
        if (mulmod64_barrett(x, x, ctx) == n) { *root = x; return 1; }
        return 0;
    }
    u64 x = powmod64_barrett(n, (p + 3) >> 3, ctx);
    if (mulmod64_barrett(x, x, ctx) == n) { *root = x; return 1; }
    x = mulmod64_barrett(x, sqrtm1, ctx);
    if (mulmod64_barrett(x, x, ctx) == n) { *root = x; return 1; }
    return 0;
}

static u64 x16_A_numerator_from_y64(u64 y, const BarrettCtx64 *ctx) {
    u64 p = ctx->p;
    u64 num = 1;
    num = submod64(mulmod64_barrett(num, y, ctx), 8, p);
    num = addmod64(mulmod64_barrett(num, y, ctx), 24, p);
    num = submod64(mulmod64_barrett(num, y, ctx), 32, p);
    num = addmod64(mulmod64_barrett(num, y, ctx), 8, p);
    num = addmod64(mulmod64_barrett(num, y, ctx), 32, p);
    num = submod64(mulmod64_barrett(num, y, ctx), 48, p);
    num = addmod64(mulmod64_barrett(num, y, ctx), 32, p);
    num = submod64(mulmod64_barrett(num, y, ctx), 8, p);
    return num;
}

static int x16_y_predicts_nonsplit64(u64 y, u64 y2, const BarrettCtx64 *ctx) {
    u64 p = ctx->p;
    u64 f1 = submod64(y2, 2, p);
    u64 four_y = addmod64(addmod64(y, y, p), addmod64(y, y, p), p);
    u64 f2 = addmod64(submod64(y2, four_y, p), 2, p);
    u64 f = mulmod64_barrett(f1, f2, ctx);
    if (f == 0) return 0;
    u64 leg = powmod64_barrett(f, (p - 1) >> 1, ctx);
    return leg != 1;
}

static int x16_root_to_montgomery_A64(u64 *Ao, u64 *xPo, u64 x, u64 y, const BarrettCtx64 *ctx) {
    u64 p = ctx->p;
    u64 num = x16_A_numerator_from_y64(y, ctx);
    u64 ym1 = submod64(y, 1, p);
    u64 ym1_2 = mulmod64_barrett(ym1, ym1, ctx);
    u64 denA = mulmod64_barrett(4, mulmod64_barrett(ym1_2, ym1_2, ctx), ctx);
    if (denA == 0) return 0;
    u64 denx = submod64(x, y, p);
    u64 vals[2] = { denA, denx }, invs[2];
    if (!invert_batch64(invs, vals, 2, ctx)) return 0;
    u64 A = mulmod64_barrett(num, invs[0], ctx);
    u64 xP = mulmod64_barrett(x, invs[1], ctx);
    if (A <= 2 || A >= p - 2) return 0;
    *Ao = A; *xPo = xP;
    return 1;
}

static int x16_sample_nonsplit_curve64(u64 *Ao, u64 *xPo, u64 *pending_A, u64 *pending_xP, int *have_pending_A,
                                       Rng *rng, u64 rand_mask, u64 sqrtm1, const BarrettCtx64 *ctx) {
    u64 p = ctx->p;
    if (*have_pending_A) {
        *Ao = *pending_A; *xPo = *pending_xP; *have_pending_A = 0;
        return 1;
    }
    for (;;) {
        u64 y = rand_below64(rng, p, rand_mask);
        if (y == 0) continue;

        u64 y2 = mulmod64_barrett(y, y, ctx);
        if (!x16_y_predicts_nonsplit64(y, y2, ctx)) continue; /* NONSPLIT FILTER */

        u64 y3 = mulmod64_barrett(y2, y, ctx);
        u64 qa = submod64(y2, addmod64(y, y, p), p);
        if (qa == 0) continue;
        u64 qb = submod64(addmod64(y2, y2, p), y3, p);
        u64 qc = submod64(1, y, p);
        u64 D = submod64(mulmod64_barrett(qb, qb, ctx),
                          mulmod64_barrett(addmod64(qa, qa, p), addmod64(qc, qc, p), ctx), p);
        u64 sd;
        if (!sqrtmod_p5_64(&sd, D, ctx, sqrtm1)) continue;

        u64 inv_2qa = invert64_barrett(addmod64(qa, qa, p), ctx);
        u64 roots[2] = {
            mulmod64_barrett(submod64(sd, qb, p), inv_2qa, ctx),
            mulmod64_barrett(submod64(p - sd, qb, p), inv_2qa, ctx)
        };

        int got_A = 0; u64 first_A = 0, first_xP = 0;
        for (int ri = 0; ri < 2; ri++) {
            u64 A, xP;
            if (!x16_root_to_montgomery_A64(&A, &xP, roots[ri], y, ctx)) continue;
            if (!got_A) { first_A = A; first_xP = xP; got_A = 1; }
            else { *Ao = first_A; *xPo = first_xP; *pending_A = A; *pending_xP = xP; *have_pending_A = 1; return 1; }
        }
        if (got_A) { *Ao = first_A; *xPo = first_xP; return 1; }
    }
}

/* ================================================================
 * Method 3 Halving: Successive Square Root Backtracking Tree — 128-bit path
 * ================================================================ */

static int halve_once_all128(u128 *xs, int max_xs, u128 p, u128 A, u128 x, u128 sqrtm1, const Mont128 *mt) {
    const u128 inv2 = (p + 1) >> 1;
    u128 x2 = mulmod128_mont(x, x, mt);
    u128 d = addmod128(addmod128(x2, mulmod128_mont(A, x, mt), p), 1, p);
    u128 sd;
    if (!sqrtmod_p5_128(&sd, d, p, sqrtm1, mt)) return 0;

    int count = 0;
    u128 roots_d[2] = { sd, submod128(0, sd, p) };
    for (int i = 0; i < 2; i++) {
        u128 u = addmod128(addmod128(x, x, p), addmod128(roots_d[i], roots_d[i], p), p);
        u128 w = submod128(mulmod128_mont(u, u, mt), 4, p);
        u128 sw;
        if (!sqrtmod_p5_128(&sw, w, p, sqrtm1, mt)) continue;
        u128 candidates[2] = {
            mulmod128_mont(addmod128(u, sw, p), inv2, mt),
            mulmod128_mont(submod128(u, sw, p), inv2, mt)
        };
        for (int j = 0; j < 2; j++) {
            if (candidates[j] != 0) {
                int dup = 0;
                for (int k = 0; k < count; k++) if (xs[k] == candidates[j]) { dup = 1; break; }
                if (!dup && count < max_xs) xs[count++] = candidates[j];
            }
        }
    }
    return count;
}

static int halve_extend128(u128 *xout, u128 p, u128 A, u128 x, int depth, int k, u128 sqrtm1, const Mont128 *mt) {
    if (depth == k) {
        if (!verify128(p, A, x)) return 0;
        *xout = x; return 1;
    }
    u128 xs[4];
    int n = halve_once_all128(xs, 4, p, A, x, sqrtm1, mt);
    for (int i = 0; i < n; i++) {
        if (halve_extend128(xout, p, A, xs[i], depth + 1, k, sqrtm1, mt)) return 1;
    }
    return 0;
}

/* ================================================================
 * Method 3 Halving: Successive Square Root Backtracking Tree — 64-bit Barrett path
 * ================================================================ */

static int halve_once_all64(u64 *xs, int max_xs, u64 A, u64 x, u64 sqrtm1, const BarrettCtx64 *ctx) {
    u64 p = ctx->p;
    const u64 inv2 = (p + 1) >> 1;
    u64 x2 = mulmod64_barrett(x, x, ctx);
    u64 d = addmod64(addmod64(x2, mulmod64_barrett(A, x, ctx), p), 1, p);
    u64 sd;
    if (!sqrtmod_p5_64(&sd, d, ctx, sqrtm1)) return 0;

    int count = 0;
    u64 roots_d[2] = { sd, submod64(0, sd, p) };
    for (int i = 0; i < 2; i++) {
        u64 u = addmod64(addmod64(x, x, p), addmod64(roots_d[i], roots_d[i], p), p);
        u64 w = submod64(mulmod64_barrett(u, u, ctx), 4, p);
        u64 sw;
        if (!sqrtmod_p5_64(&sw, w, ctx, sqrtm1)) continue;
        u64 candidates[2] = {
            mulmod64_barrett(addmod64(u, sw, p), inv2, ctx),
            mulmod64_barrett(submod64(u, sw, p), inv2, ctx)
        };
        for (int j = 0; j < 2; j++) {
            if (candidates[j] != 0) {
                int dup = 0;
                for (int k = 0; k < count; k++) if (xs[k] == candidates[j]) { dup = 1; break; }
                if (!dup && count < max_xs) xs[count++] = candidates[j];
            }
        }
    }
    return count;
}

static int halve_extend64(u64 *xout, u64 A, u64 x, int depth, int k, u64 sqrtm1, const BarrettCtx64 *ctx) {
    if (depth == k) {
        if (!verify64_barrett(ctx->p, A, x, ctx)) return 0;
        *xout = x; return 1;
    }
    u64 xs[4];
    int n = halve_once_all64(xs, 4, A, x, sqrtm1, ctx);
    for (int i = 0; i < n; i++) {
        if (halve_extend64(xout, A, xs[i], depth + 1, k, sqrtm1, ctx)) return 1;
    }
    return 0;
}

/* ================================================================
 * Shared 2-Sylow Depth Calculation
 * ================================================================ */

static int compute_k(u128 p) {
    u64 q = (u64)sqrtl((long double)p);
    while ((u128)(q + 1) * (q + 1) <= p) q++;
    while ((u128)q * q > p) q--;
    u64 sq = (u64)sqrtl((long double)q);
    while ((sq + 1) * (sq + 1) <= q) sq++;
    while (sq * sq > q) sq--;
    u64 bound = q + 1 + 2 * sq;
    int k = 0; u64 v = 1; while (v <= bound) { k++; v <<= 1; } return k;
}

/* ================================================================
 * Parallel Search Engine — 128-bit Montgomery fallback (p >= 2^63)
 * ================================================================ */

static int search128(u128 p, int target_total, int start_count,
                     u128 *out_A, u128 *out_x0, u64 *out_trials,
                     u64 *out_cycles, u64 *out_mul, u64 *out_inv, u64 *out_branch) {
    volatile int found_count = start_count;

    int k = compute_k(p);
    u64 sqrtp = (u64)sqrtl((long double)p);
    while ((u128)(sqrtp + 1) * (sqrtp + 1) <= p) sqrtp++;
    while ((u128)sqrtp * sqrtp > p) sqrtp--;

    Mont128 mt; m128_init(&mt, p);
    int pbits = bitlen128(p);
    u128 rand_mask = pbits >= 128 ? (u128)0 - 1 : (((u128)1 << pbits) - 1);

    u128 sqrtm1 = 0;
    if ((u64)(p & 7) == 5) {
        sqrtm1 = powmod128_mont(2, (p - 1) >> 2, &mt);
    } else if ((u64)(p & 3) != 3) {
        printf("X1(16) mode requires p ≡ 5 mod 8 or p ≡ 3 mod 4.\n");
        return found_count;
    }

    u64 max_trials = (u64)(5.0 * (double)sqrtp); /* Target tree attempts scaled down for high hit-rate */
    if (max_trials < 1000000ULL) max_trials = 1000000ULL;

    memset(g_thread_counters, 0, sizeof(g_thread_counters));

#pragma omp parallel
    {
        int tid = 0, nthr = 1;
#ifdef _OPENMP
        tid = omp_get_thread_num(); nthr = omp_get_num_threads();
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
        u128 pending_A = 0, pending_xP = 0; int have_pending_A = 0;

        while (found_count < target_total && A_trials < budget) {
            u128 A, xP16;
            x16_sample_nonsplit_curve128(&A, &xP16, &pending_A, &pending_xP, &have_pending_A, &rng, p, rand_mask, sqrtm1, &mt);

            A_trials++;
            g_thread_counters[tid].count = A_trials;

            t_ops.mul_count = 0; t_ops.inv_count = 0; t_ops.branch_count = 0;
            uint64_t c_start = rdtsc();

            u128 xR;
            int verify_ok = halve_extend128(&xR, p, A, xP16, 4, k, sqrtm1, &mt);

            uint64_t c_end = rdtsc();
            uint64_t cyc = c_end - c_start;
            OpCounter this_ops = t_ops;

            if (verify_ok) {
#pragma omp critical
                {
                    int is_dup = 0;
                    for (int j = 0; j < found_count; j++) {
                        if (out_A[j] == A && out_x0[j] == xR) { is_dup = 1; break; }
                    }
                    if (!is_dup && found_count < target_total) {
                        out_A[found_count] = A;
                        out_x0[found_count] = xR;

                        u64 exact_total_trials = 0;
                        for (int t = 0; t < nthr; t++) exact_total_trials += g_thread_counters[t].count;
                        if (exact_total_trials == 0) exact_total_trials = 1;

                        out_trials[found_count] = exact_total_trials;
                        out_cycles[found_count] = cyc;
                        out_mul[found_count] = this_ops.mul_count;
                        out_inv[found_count] = this_ops.inv_count;
                        out_branch[found_count] = this_ops.branch_count;

                        found_count++;
                    }
                }
            }
        }
    }

    return found_count;
}

/* ================================================================
 * Parallel Search Engine — 64-bit Barrett full pipeline (p < 2^63)
 * ================================================================ */

static int search64(u64 p, int target_total, int start_count,
                    u64 *out_A, u64 *out_x0, u64 *out_trials,
                    u64 *out_cycles, u64 *out_mul, u64 *out_inv, u64 *out_branch) {
    volatile int found_count = start_count;

    int k = compute_k((u128)p);
    u64 sqrtp = (u64)sqrtl((long double)p);
    while ((u128)(sqrtp + 1) * (sqrtp + 1) <= (u128)p) sqrtp++;
    while ((u128)sqrtp * sqrtp > (u128)p) sqrtp--;

    BarrettCtx64 ctx = barrett_setup64(p);
    int pbits = bitlen64(p);
    u64 rand_mask = pbits >= 64 ? (u64)0 - 1 : (((u64)1 << pbits) - 1);

    u64 sqrtm1 = 0;
    if ((p & 7) == 5) {
        sqrtm1 = powmod64_barrett(2, (p - 1) >> 2, &ctx);
    } else if ((p & 3) != 3) {
        printf("X1(16) mode requires p ≡ 5 mod 8 or p ≡ 3 mod 4.\n");
        return found_count;
    }

    u64 max_trials = (u64)(5.0 * (double)sqrtp); /* Target tree attempts scaled down for high hit-rate */
    if (max_trials < 1000000ULL) max_trials = 1000000ULL;

    memset(g_thread_counters, 0, sizeof(g_thread_counters));

#pragma omp parallel
    {
        int tid = 0, nthr = 1;
#ifdef _OPENMP
        tid = omp_get_thread_num(); nthr = omp_get_num_threads();
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
        u64 pending_A = 0, pending_xP = 0; int have_pending_A = 0;

        while (found_count < target_total && A_trials < budget) {
            u64 A, xP16;
            /* STEP 1: Sample Nonsplit X1(16) Curve & Marked 16-Torsion Point — full Barrett
             * u64 arithmetic, so every sampled curve (not just accepted ones) benefits,
             * unlike the earlier verify-only patch. */
            x16_sample_nonsplit_curve64(&A, &xP16, &pending_A, &pending_xP, &have_pending_A,
                                        &rng, rand_mask, sqrtm1, &ctx);

            A_trials++;
            g_thread_counters[tid].count = A_trials;

            /* STEP 2: Successive Halving Backtracking (DFS), also fully Barrett-reduced —
             * every branch explored at every depth uses Barrett arithmetic, not just the
             * accepted leaf. */
            t_ops.mul_count = 0; t_ops.inv_count = 0; t_ops.branch_count = 0;
            uint64_t c_start = rdtsc();

            u64 xR;
            int verify_ok = halve_extend64(&xR, A, xP16, 4, k, sqrtm1, &ctx);

            uint64_t c_end = rdtsc();
            uint64_t cyc = c_end - c_start;
            OpCounter this_ops = t_ops;

            if (verify_ok) {
#pragma omp critical
                {
                    int is_dup = 0;
                    for (int j = 0; j < found_count; j++) {
                        if (out_A[j] == A && out_x0[j] == xR) { is_dup = 1; break; }
                    }
                    if (!is_dup && found_count < target_total) {
                        out_A[found_count] = A;
                        out_x0[found_count] = xR;

                        u64 exact_total_trials = 0;
                        for (int t = 0; t < nthr; t++) exact_total_trials += g_thread_counters[t].count;
                        if (exact_total_trials == 0) exact_total_trials = 1;

                        out_trials[found_count] = exact_total_trials;
                        out_cycles[found_count] = cyc;
                        out_mul[found_count] = this_ops.mul_count;
                        out_inv[found_count] = this_ops.inv_count;
                        out_branch[found_count] = this_ops.branch_count;

                        found_count++;
                    }
                }
            }
        }
    }

    return found_count;
}

/* ================================================================
 * Batch-Processing Main — dispatches per-prime to search64 (p < 2^63)
 * or search128 (p >= 2^63)
 * ================================================================ */

int main(int argc, char *argv[]) {
    if (argc < 5) {
        printf("Usage: ./pomerance_x16_halving_barrett_full <stateful_input.txt> <output_pure.csv> <output_metrics.csv> <target_total>\n");
        return 1;
    }

    int target_total = atoi(argv[4]);
    if (target_total <= 0) { printf("Error: <target_total> must be > 0.\n"); return 1; }

    FILE *input = fopen(argv[1], "r");
    FILE *out_pure = fopen(argv[2], "w");
    FILE *out_metrics = fopen(argv[3], "w");
    if (!input || !out_pure || !out_metrics) { printf("Error opening files.\n"); return 1; }

    u64 *out_A_arr    = (u64 *)malloc(target_total * sizeof(u64));
    u64 *out_x0_arr   = (u64 *)malloc(target_total * sizeof(u64));
    u128 *out_A128    = (u128 *)malloc(target_total * sizeof(u128));
    u128 *out_x0128   = (u128 *)malloc(target_total * sizeof(u128));
    u64 *out_trials   = (u64 *)malloc(target_total * sizeof(u64));
    u64 *out_cycles   = (u64 *)malloc(target_total * sizeof(u64));
    u64 *out_mul      = (u64 *)malloc(target_total * sizeof(u64));
    u64 *out_inv      = (u64 *)malloc(target_total * sizeof(u64));
    u64 *out_branch   = (u64 *)malloc(target_total * sizeof(u64));

    if (!out_A_arr || !out_x0_arr || !out_A128 || !out_x0128 || !out_trials ||
        !out_cycles || !out_mul || !out_inv || !out_branch) {
        printf("Memory allocation failed.\n"); return 1;
    }

    fprintf(out_metrics, "prime,A,x0,trials,cycles,mul_count,inv_count,branch_count,opcount_m3,batch_time_ms\n");

    char line[65536]; unsigned long long current_index = 0;

    while (fgets(line, sizeof(line), input)) {
        char *ptr = line;
        u128 p = parse128_adv(&ptr);
        if (p == 0) continue;

        int use64 = (p < ((u128)1 << 63));

        current_index++;
        int num_existing = (int)parse128_adv(&ptr);
        for (int i = 0; i < num_existing; i++) {
            u128 Aval = parse128_adv(&ptr);
            u128 xval = parse128_adv(&ptr);
            if (use64) {
                out_A_arr[i] = (u64)Aval;
                out_x0_arr[i] = (u64)xval;
            } else {
                out_A128[i] = Aval;
                out_x0128[i] = xval;
            }
            out_trials[i] = 0; out_cycles[i] = 0; out_mul[i] = 0; out_inv[i] = 0; out_branch[i] = 0;
        }

        char p_display[50]; sprint128(p_display, p);
        printf("[%llu] Prime: %s (Existing: %d, Target: %d)...\n", current_index, p_display, num_existing, target_total);
        fflush(stdout);

        double t_start = now_sec();
        int found_amount;
        if (use64) {
            found_amount = search64((u64)p, target_total, num_existing, out_A_arr, out_x0_arr, out_trials,
                                    out_cycles, out_mul, out_inv, out_branch);
        } else {
            found_amount = search128(p, target_total, num_existing, out_A128, out_x0128, out_trials,
                                     out_cycles, out_mul, out_inv, out_branch);
        }
        double batch_time_ms = (now_sec() - t_start) * 1000.0;
        int newly_found = found_amount - num_existing;

        if (newly_found > 0) {
            for (int i = num_existing; i < found_amount - 1; i++) {
                for (int j = num_existing; j < found_amount - 1 - (i - num_existing); j++) {
                    if (out_trials[j] > out_trials[j + 1]) {
                        u64 tmp;
                        tmp = out_trials[j]; out_trials[j] = out_trials[j + 1]; out_trials[j + 1] = tmp;
                        tmp = out_cycles[j]; out_cycles[j] = out_cycles[j + 1]; out_cycles[j + 1] = tmp;
                        tmp = out_mul[j]; out_mul[j] = out_mul[j + 1]; out_mul[j + 1] = tmp;
                        tmp = out_inv[j]; out_inv[j] = out_inv[j + 1]; out_inv[j + 1] = tmp;
                        tmp = out_branch[j]; out_branch[j] = out_branch[j + 1]; out_branch[j + 1] = tmp;
                        if (use64) {
                            tmp = out_A_arr[j]; out_A_arr[j] = out_A_arr[j + 1]; out_A_arr[j + 1] = tmp;
                            tmp = out_x0_arr[j]; out_x0_arr[j] = out_x0_arr[j + 1]; out_x0_arr[j + 1] = tmp;
                        } else {
                            u128 tmp128;
                            tmp128 = out_A128[j]; out_A128[j] = out_A128[j + 1]; out_A128[j + 1] = tmp128;
                            tmp128 = out_x0128[j]; out_x0128[j] = out_x0128[j + 1]; out_x0128[j + 1] = tmp128;
                        }
                    }
                }
            }

            u64 total_new_trials = 0, last_cumulative = 0, total_cycles = 0;
            double total_m3 = 0.0;

            for (int i = num_existing; i < found_amount; i++) {
                u64 marginal_trials = out_trials[i] - last_cumulative;
                if (marginal_trials == 0) marginal_trials = 1;
                last_cumulative = out_trials[i];
                total_new_trials += marginal_trials; total_cycles += out_cycles[i];

                OpCounter oc = { out_mul[i], out_inv[i], out_branch[i] };
                double m3 = opcount_m3(&oc); total_m3 += m3;

                char a_str[50], x_str[50];
                sprint128(a_str, use64 ? (u128)out_A_arr[i] : out_A128[i]);
                sprint128(x_str, use64 ? (u128)out_x0_arr[i] : out_x0128[i]);
                fprintf(out_pure, "%s,%s,%s\n", p_display, a_str, x_str);
                fprintf(out_metrics, "%s,%s,%s,%lu,%lu,%lu,%lu,%lu,%.3f,%.2f\n",
                        p_display, a_str, x_str, marginal_trials, out_cycles[i], out_mul[i], out_inv[i], out_branch[i], m3, batch_time_ms);
            }
            printf("      Success: found %d new triples in %.2f ms | Trials: %lu | M3: %.1f\n",
                   newly_found, batch_time_ms, total_new_trials, total_m3);
        } else {
            fprintf(out_pure, "%s,FAILED,FAILED\n", p_display);
            fprintf(out_metrics, "%s,FAILED,FAILED,FAILED,FAILED,FAILED,FAILED,FAILED,FAILED,FAILED\n", p_display);
            printf("      Failed to find new triples (Time elapsed: %.2f ms).\n", batch_time_ms);
        }
        fflush(stdout);
    }

    free(out_A_arr); free(out_x0_arr); free(out_A128); free(out_x0128);
    free(out_trials); free(out_cycles); free(out_mul); free(out_inv); free(out_branch);
    fclose(input); fclose(out_pure); fclose(out_metrics);
    return 0;
}
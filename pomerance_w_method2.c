/*
 * pomerance_benchmarked.c — Incremental generation of Pomerance triples (p, A, x0)
 *
 * This is the brute-force (no 2-Sylow projection) variant, instrumented with
 * three benchmarks captured per successful trial / prime batch:
 *
 * 1. CPU cycles (rdtsc)   — exact hardware cycle count spent inside the
 * verify64()/verify128() call that succeeded.
 * 2. OpCount (M3)         — weighted "algorithmic work" count, built from
 * manual instrumentation of modular multiplies,
 * modular inversions, and branch decisions taken
 * inside that same verify call.
 * 3. Wall-Clock Time      — physical time (in milliseconds) spent hunting
 * for the target triples of the current prime.
 *
 * Threads are pinned to a CPU core (pthread_setaffinity_np) at the start of
 * each parallel region so cycle counts aren't skewed by OS migration or
 * frequency scaling across cores.
 *
 * Optional X1(16) prescribed-torsion sampling mode (CLI flag "x16"):
 * instead of drawing A uniformly at random, A is sampled from the
 * one-parameter family of Montgomery curves carrying a rational point of
 * order 16, via the Tate normal form parametrization of the modular curve
 * X1(16). This forces 16 | N for the resulting curve, concentrating trials
 * on curves already biased toward large 2-power torsion. Because the
 * sampled point already has order (a multiple of) 16, there is exactly one
 * candidate x0 to test per A (the marked point) instead of the 50 random
 * x0 tries used in brute-force mode.
 *
 * Compile:
 * gcc -O3 -fopenmp -o pomerance pomerance_benchmarked.c -lm -lpthread
 * gcc -O3 -o pomerance pomerance_benchmarked.c -lm -lpthread   (single-threaded)
 *
 * Usage:
 * ./pomerance <stateful_input.txt> <output_pure_new.csv> <output_metrics_new.csv> <target_total> [x16]
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
 * Global mode flags
 * ================================================================ */

static int g_x16_mode = 0;

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
 * Benchmark #2: OpCount (manual instrumentation of algorithmic work)
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
    char tmp[50]; int i = 49; tmp[i] = '\0';
    while (v > 0) { tmp[--i] = '0' + (int)(v % 10); v /= 10; }
    strcpy(buf, tmp + i);
}

/* ================================================================
 * PRNG (xorshift128+)
 * ================================================================ */

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

static inline u128 rand_below128(Rng *rng, u128 p, u128 mask) {
    for (;;) {
        u128 v = ((u128)rng64(rng) << 64) | (u128)rng64(rng);
        v &= mask;
        if (v < p) return v;
    }
}

/* ================================================================
 * Benchmark #3: High-Resolution Monotonic Timer
 * ================================================================ */

static double now_sec(void) {
#ifdef _OPENMP
    return omp_get_wtime();
#else
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
#endif
}

/* ================================================================
 * u64 code path (p < 2^63)
 * ================================================================ */

static inline u64 addmod64(u64 a, u64 b, u64 p) { return a >= p-b ? a-(p-b) : a+b; }
static inline u64 submod64(u64 a, u64 b, u64 p) { return a >= b ? a-b : p-b+a; }

static inline u64 mulmod64(u64 a, u64 b, u64 p) {
    t_ops.mul_count++;
    return (u64)((u128)a*b%(u128)p);
}

static int verify64(u64 p, u64 A, u64 x0) {
    u64 q = (u64)sqrtl((long double)p);
    while ((u128)(q+1)*(q+1)<=(u128)p) q++;
    while ((u128)q*q>(u128)p) q--;
    u64 sq = (u64)sqrtl((long double)q);
    while ((sq+1)*(sq+1)<=q) sq++;
    while (sq*sq>q) sq--;
    u64 bound = q+1+2*sq;
    int k=0; u64 v=1; while(v<=bound){k++;v<<=1;}

    if (A%p==2||A%p==p-2) { t_ops.branch_count++; return 0; }
    u64 X=x0%p, Z=1;
    for (int i=1; i<=k; i++) {
        u64 X2=mulmod64(X,X,p), Z2=mulmod64(Z,Z,p), XZ=mulmod64(X,Z,p);
        u64 d=submod64(X2,Z2,p), Xn=mulmod64(d,d,p);
        u64 inn=addmod64(addmod64(X2,mulmod64(A,XZ,p),p),Z2,p);
        u64 f4=addmod64(addmod64(XZ,XZ,p),addmod64(XZ,XZ,p),p);
        u64 Zn=mulmod64(f4,inn,p); X=Xn; Z=Zn;
        if (i<k&&Z==0) { t_ops.branch_count++; return 0; }
        if (i==k&&Z!=0) { t_ops.branch_count++; return 0; }
    }
    return 1;
}

/* ================================================================
 * u128 code path (p < 2^127)
 * ================================================================ */

static inline u128 addmod128(u128 a, u128 b, u128 p) { u128 s=a+b; return s>=p?s-p:s; }
static inline u128 submod128(u128 a, u128 b, u128 p) { return a>=b?a-b:p-b+a; }

static u128 mulmod_slow(u128 a, u128 b, u128 p) {
    t_ops.mul_count++;
    u128 r=0; a%=p; b%=p;
    while (b>0) { if(b&1){r+=a;if(r>=p)r-=p;} a+=a;if(a>=p)a-=p; b>>=1; }
    return r;
}

static int verify128(u128 p, u128 A, u128 x0) {
    u64 q = (u64)sqrtl((long double)p);
    while ((u128)(q+1)*(q+1)<=p) q++;
    while ((u128)q*q>p) q--;
    u64 sq = (u64)sqrtl((long double)q);
    while ((sq+1)*(sq+1)<=q) sq++;
    while (sq*sq>q) sq--;
    u64 bound = q+1+2*sq;
    int k=0; u64 v=1; while(v<=bound){k++;v<<=1;}

    if (A%p==2||A%p==p-2) { t_ops.branch_count++; return 0; }
    u128 X=x0%p, Z=1;
    for (int i=1; i<=k; i++) {
        u128 X2=mulmod_slow(X,X,p), Z2=mulmod_slow(Z,Z,p), XZ=mulmod_slow(X,Z,p);
        u128 d=submod128(X2,Z2,p), Xn=mulmod_slow(d,d,p);
        u128 inn=addmod128(addmod128(X2,mulmod_slow(A,XZ,p),p),Z2,p);
        u128 f4=addmod128(addmod128(XZ,XZ,p),addmod128(XZ,XZ,p),p);
        u128 Zn=mulmod_slow(f4,inn,p); X=Xn; Z=Zn;
        if (i<k&&Z==0) { t_ops.branch_count++; return 0; }
        if (i==k&&Z!=0) { t_ops.branch_count++; return 0; }
    }
    return 1;
}

/* ================================================================
 * X1(16) prescribed-torsion sampling (u128 path)
 * ================================================================ */

typedef struct { u128 lo, hi; } u256;

static inline u256 wide_mul(u128 a, u128 b) {
    u64 a0=(u64)a, a1=(u64)(a>>64), b0=(u64)b, b1=(u64)(b>>64);
    u128 ll=(u128)a0*b0, lh=(u128)a0*b1, hl=(u128)a1*b0, hh=(u128)a1*b1;
    u128 mid=lh+hl; u128 carry_mid=(mid<lh)?1:0;
    u128 lo=ll+(mid<<64); u128 carry_lo=(lo<ll)?1:0;
    return (u256){lo, hh+(mid>>64)+(carry_mid<<64)+carry_lo};
}

static inline u256 wide_add(u256 a, u256 b) {
    u128 lo=a.lo+b.lo; return (u256){lo, a.hi+b.hi+((lo<a.lo)?1:0)};
}

typedef struct { u128 p, ni, R2, one; } Mont128;

static void m128_init(Mont128 *m, u128 p) {
    m->p = p;
    u128 x = 1; for (int i = 0; i < 7; i++) x *= 2 - p * x;
    m->ni = (u128)0 - x;
    u128 r = 1;
    for (int i = 0; i < 128; i++) { r<<=1; if(r>=p) r-=p; }
    m->one = r;
    for (int i = 0; i < 128; i++) { r<<=1; if(r>=p) r-=p; }
    m->R2 = r;
}

static inline u128 mred128(u256 T, const Mont128 *m) {
    u128 q = T.lo * m->ni;
    u256 s = wide_add(T, wide_mul(q, m->p));
    u128 t = s.hi;
    return t >= m->p ? t - m->p : t;
}
static inline u128 mm128(u128 a, u128 b, const Mont128 *m) { return mred128(wide_mul(a,b), m); }
static inline u128 toM128(u128 a, const Mont128 *m) { return mm128(a % m->p, m->R2, m); }
static inline u128 frM128(u128 a, const Mont128 *m) { return mred128((u256){a,0}, m); }

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
    return powmod128_mont(a, p - 2, mt);
}

static int sqrtmod_p5_128(u128 *root, u128 n, u128 p, u128 sqrtm1, const Mont128 *mt) {
    n %= p;
    if (n == 0) { *root = 0; return 1; }
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

static int x16_root_to_montgomery_A128(u128 *Ao, u128 *xPo,
                                       u128 p, u128 x, u128 y,
                                       const Mont128 *mt) {
    u128 num = x16_A_numerator_from_y128(y, p, mt);

    u128 ym1 = submod128(y, 1, p);
    u128 ym1_2 = mulmod128_mont(ym1, ym1, mt);
    u128 denA = mulmod128_mont(4, mulmod128_mont(ym1_2, ym1_2, mt), mt);
    u128 denx = submod128(x, y, p);
    if (denA == 0 || denx == 0) return 0;

    u128 A = mulmod128_mont(num, invert128_mont(denA, p, mt), mt);
    u128 xP = mulmod128_mont(x, invert128_mont(denx, p, mt), mt);
    if (A <= 2 || A >= p - 2) return 0;
    *Ao = A;
    *xPo = xP;
    return 1;
}

static int x16_montgomery_A128(u128 *Ao, u128 *xPo,
                               u128 *pending_A, u128 *pending_xP, int *have_pending_A,
                               Rng *rng, u128 p, u128 rand_mask,
                               u128 sqrtm1, const Mont128 *mt) {
    if (*have_pending_A) {
        *Ao = *pending_A;
        *xPo = *pending_xP;
        *have_pending_A = 0;
        return 1;
    }

    for (;;) {
        u128 y = rand_below128(rng, p, rand_mask);
        if (y == 0) continue;

        u128 y2 = mulmod128_mont(y, y, mt);
        u128 y3 = mulmod128_mont(y2, y, mt);
        u128 qa = submod128(y2, addmod128(y, y, p), p);
        if (qa == 0) continue;
        u128 qb = submod128(addmod128(y2, y2, p), y3, p);
        u128 qc = submod128(1, y, p);
        u128 D = submod128(mulmod128_mont(qb, qb, mt),
                           mulmod128_mont(addmod128(qa, qa, p), addmod128(qc, qc, p), mt),
                           p);
        u128 sd;
        if (!sqrtmod_p5_128(&sd, D, p, sqrtm1, mt)) continue;

        u128 inv_2qa = invert128_mont(addmod128(qa, qa, p), p, mt);
        u128 roots[2] = {
            mulmod128_mont(submod128(sd, qb, p), inv_2qa, mt),
            mulmod128_mont(submod128(p - sd, qb, p), inv_2qa, mt)
        };

        int got_A = 0;
        u128 first_A = 0;
        u128 first_xP = 0;
        for (int ri = 0; ri < 2; ri++) {
            u128 A, xP;
            if (!x16_root_to_montgomery_A128(&A, &xP, p, roots[ri], y, mt)) continue;
            if (!got_A) {
                first_A = A;
                first_xP = xP;
                got_A = 1;
            } else {
                *Ao = first_A;
                *xPo = first_xP;
                *pending_A = A;
                *pending_xP = xP;
                *have_pending_A = 1;
                return 1;
            }
        }
        if (got_A) {
            *Ao = first_A;
            *xPo = first_xP;
            return 1;
        }
    }
}

/* ================================================================
 * X1(16) prescribed-torsion sampling (u64 path)
 * ================================================================ */

typedef struct { u64 p, ni, R2, one; } Mont64;

static void m64_init(Mont64 *m, u64 p) {
    m->p = p;
    u64 x = 1;
    for (int i = 0; i < 6; i++) x *= 2 - p * x;
    m->ni = (u64)0 - x;
    u128 r = 1;
    for (int i = 0; i < 64; i++) { r <<= 1; if (r >= p) r -= p; }
    m->one = (u64)r;
    for (int i = 0; i < 64; i++) { r <<= 1; if (r >= p) r -= p; }
    m->R2 = (u64)r;
}

static inline u64 mred64(u128 T, const Mont64 *m) {
    u64 q = (u64)T * m->ni;
    u128 s = T + (u128)q * m->p;
    u64 t = (u64)(s >> 64);
    return t >= m->p ? t - m->p : t;
}

static inline u64 mm64(u64 a, u64 b, const Mont64 *m) { return mred64((u128)a * b, m); }
static inline u64 toM64(u64 a, const Mont64 *m) { return mm64(a % m->p, m->R2, m); }
static inline u64 frM64(u64 a, const Mont64 *m) { return mred64((u128)a, m); }

static inline u64 mulmod64_mont(u64 a, u64 b, const Mont64 *mt) {
    return frM64(mm64(toM64(a, mt), toM64(b, mt), mt), mt);
}

static u64 powmod64_mont(u64 a, u64 e, const Mont64 *mt) {
    u64 r = mt->one, b = toM64(a, mt);
    while (e > 0) {
        if (e & 1) r = mm64(r, b, mt);
        b = mm64(b, b, mt);
        e >>= 1;
    }
    return frM64(r, mt);
}

static u64 invert64_mont(u64 a, u64 p, const Mont64 *mt) {
    return powmod64_mont(a, p - 2, mt);
}

static int sqrtmod_p5_64(u64 *root, u64 n, u64 p, u64 sqrtm1, const Mont64 *mt) {
    n %= p;
    if (n == 0) { *root = 0; return 1; }
    u64 x = powmod64_mont(n, (p + 3) >> 3, mt);
    if (mulmod64_mont(x, x, mt) == n) { *root = x; return 1; }
    x = mulmod64_mont(x, sqrtm1, mt);
    if (mulmod64_mont(x, x, mt) == n) { *root = x; return 1; }
    return 0;
}

static u64 x16_A_numerator_from_y64(u64 y, u64 p, const Mont64 *mt) {
    u64 num = 1;
    num = submod64(mulmod64_mont(num, y, mt), 8, p);
    num = addmod64(mulmod64_mont(num, y, mt), 24, p);
    num = submod64(mulmod64_mont(num, y, mt), 32, p);
    num = addmod64(mulmod64_mont(num, y, mt), 8, p);
    num = addmod64(mulmod64_mont(num, y, mt), 32, p);
    num = submod64(mulmod64_mont(num, y, mt), 48, p);
    num = addmod64(mulmod64_mont(num, y, mt), 32, p);
    num = submod64(mulmod64_mont(num, y, mt), 8, p);
    return num;
}

static int x16_root_to_montgomery_A64(u64 *Ao, u64 *xPo,
                                       u64 p, u64 x, u64 y,
                                       const Mont64 *mt) {
    u64 num = x16_A_numerator_from_y64(y, p, mt);

    u64 ym1 = submod64(y, 1, p);
    u64 ym1_2 = mulmod64_mont(ym1, ym1, mt);
    u64 denA = mulmod64_mont(4, mulmod64_mont(ym1_2, ym1_2, mt), mt);
    u64 denx = submod64(x, y, p);
    if (denA == 0 || denx == 0) return 0;

    u64 A = mulmod64_mont(num, invert64_mont(denA, p, mt), mt);
    u64 xP = mulmod64_mont(x, invert64_mont(denx, p, mt), mt);
    if (A <= 2 || A >= p - 2) return 0;
    *Ao = A;
    *xPo = xP;
    return 1;
}

static int x16_montgomery_A64(u64 *Ao, u64 *xPo,
                               u64 *pending_A, u64 *pending_xP, int *have_pending_A,
                               Rng *rng, u64 p, u64 rand_mask,
                               u64 sqrtm1, const Mont64 *mt) {
    if (*have_pending_A) {
        *Ao = *pending_A;
        *xPo = *pending_xP;
        *have_pending_A = 0;
        return 1;
    }

    for (;;) {
        u64 y;
        for (;;) {
            y = rng64(rng) & rand_mask;
            if (y < p && y != 0) break;
        }

        u64 y2 = mulmod64_mont(y, y, mt);
        u64 y3 = mulmod64_mont(y2, y, mt);
        u64 qa = submod64(y2, addmod64(y, y, p), p);
        if (qa == 0) continue;
        u64 qb = submod64(addmod64(y2, y2, p), y3, p);
        u64 qc = submod64(1, y, p);
        u64 D = submod64(mulmod64_mont(qb, qb, mt),
                           mulmod64_mont(addmod64(qa, qa, p), addmod64(qc, qc, p), mt),
                           p);
        u64 sd;
        if (!sqrtmod_p5_64(&sd, D, p, sqrtm1, mt)) continue;

        u64 inv_2qa = invert64_mont(addmod64(qa, qa, p), p, mt);
        u64 roots[2] = {
            mulmod64_mont(submod64(sd, qb, p), inv_2qa, mt),
            mulmod64_mont(submod64(p - sd, qb, p), inv_2qa, mt)
        };

        int got_A = 0;
        u64 first_A = 0, first_xP = 0;
        for (int ri = 0; ri < 2; ri++) {
            u64 A, xP;
            if (!x16_root_to_montgomery_A64(&A, &xP, p, roots[ri], y, mt)) continue;
            if (!got_A) {
                first_A = A;
                first_xP = xP;
                got_A = 1;
            } else {
                *Ao = first_A;
                *xPo = first_xP;
                *pending_A = A;
                *pending_xP = xP;
                *have_pending_A = 1;
                return 1;
            }
        }
        if (got_A) {
            *Ao = first_A;
            *xPo = first_xP;
            return 1;
        }
    }
}

/* ================================================================
 * Dispatch: search64 / search128
 * ================================================================ */

static int search64(u64 p, int target_total, int start_count,
                     u64 *out_A, u64 *out_x0, u64 *out_trials,
                     u64 *out_cycles, u64 *out_mul, u64 *out_inv, u64 *out_branch,
                     int x16_mode) {
    volatile int found_count = start_count;

    u64 sqrtp = (u64)sqrtl((long double)p);
    while ((u128)(sqrtp+1)*(sqrtp+1)<=(u128)p) sqrtp++;
    while ((u128)sqrtp*sqrtp>(u128)p) sqrtp--;

    u64 max_trials = (u64)(20.0 * (double)sqrtp);
    if (max_trials < 10000000ULL) max_trials = 10000000ULL;

    u64 thread_trials[256 * 8] = {0};

    Mont64 mt16;
    u64 sqrtm1_16 = 0;
    int x16_ready = 0;

    if (x16_mode) {
        if ((p & 7) == 5) {
            m64_init(&mt16, p);
            sqrtm1_16 = powmod64_mont(2, (p - 1) >> 2, &mt16);
            x16_ready = 1;
        } else {
            printf("      [x16] p is not congruent to 5 (mod 8); falling back to brute-force sampling.\n");
        }
    }

    int pbits = bitlen128(p);
    u64 rand_mask16 = (((u128)1 << pbits) - 1);

#pragma omp parallel
    {
        int tid=0, nthr=1;
#ifdef _OPENMP
        tid = omp_get_thread_num(); nthr = omp_get_num_threads();
#endif
        pin_thread_to_core(tid);

        u64 current_time = (u64)time(NULL);
        Rng rng = {
            .s0 = 7364529176530163ULL ^ ((u64)tid * 6364136223846793005ULL) ^ p ^ current_time,
            .s1 = 1442695040888963407ULL ^ ((u64)(tid+1) * 2862933555777941757ULL) ^ (current_time << 32)
        };
        for (int i=0;i<200;i++) rng64(&rng);

        u64 budget = max_trials / nthr + 1;
        u64 A_trials = 0;

        u64 pending_A = 0, pending_xP = 0;
        int have_pending_A = 0;

        while (found_count < target_total && A_trials < budget) {
            u64 A;
            u64 xP16 = 0;

            if (x16_mode && x16_ready) {
                x16_montgomery_A64(&A, &xP16, &pending_A, &pending_xP, &have_pending_A,
                                     &rng, p, rand_mask16, sqrtm1_16, &mt16);
            } else {
                A = rng64(&rng) % p;
                if (A==2 || A==p-2) { continue; }
            }

            A_trials++;
            thread_trials[tid * 8] = A_trials;

            int a_success = 0;
            int max_x_tries = (x16_mode && x16_ready) ? 1 : 50;

            for (int x_tries = 0; x_tries < max_x_tries && !a_success && found_count < target_total; x_tries++) {
                u64 x0r;
                if (x16_mode && x16_ready) {
                    x0r = xP16;
                } else {
                    x0r = rng64(&rng) % p;
                    if (x0r < 2) continue;
                }

                t_ops.mul_count = 0; t_ops.inv_count = 0; t_ops.branch_count = 0;
                uint64_t c_start = rdtsc();
                int ok = verify64(p, A, x0r);
                uint64_t c_end = rdtsc();
                uint64_t cyc = c_end - c_start;
                OpCounter this_ops = t_ops;

                if (ok) {
#pragma omp critical
                    {
                        int is_dup = 0;
                        for(int j = 0; j < found_count; j++) {
                            if(out_A[j] == A && out_x0[j] == x0r) {
                                is_dup = 1; break;
                            }
                        }
                        if (!is_dup && found_count < target_total) {
                            out_A[found_count] = A;
                            out_x0[found_count] = x0r;

                            u64 exact_total_trials = 0;
                            for (int t = 0; t < nthr; t++) {
                                exact_total_trials += thread_trials[t * 8];
                            }
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
                      u64 *out_cycles, u64 *out_mul, u64 *out_inv, u64 *out_branch,
                      int x16_mode) {
    volatile int found_count = start_count;

    u64 sqrtp = (u64)sqrtl((long double)p);
    while ((u128)(sqrtp+1)*(sqrtp+1)<=p) sqrtp++;
    while ((u128)sqrtp*sqrtp>p) sqrtp--;

    u64 max_trials = (u64)(20.0 * (double)sqrtp);
    if (max_trials < 10000000ULL) max_trials = 10000000ULL;

    u64 thread_trials[256 * 8] = {0};

    Mont128 mt16;
    u128 sqrtm1_16 = 0;
    int x16_ready = 0;
    if (x16_mode) {
        if ((u64)(p & 7) == 5) {
            m128_init(&mt16, p);
            sqrtm1_16 = powmod128_mont(2, (p - 1) >> 2, &mt16);
            x16_ready = 1;
        } else {
            printf("      [x16] p is not congruent to 5 (mod 8); falling back to brute-force sampling.\n");
        }
    }
    int pbits = bitlen128(p);
    u128 rand_mask16 = pbits >= 128 ? (u128)0 - 1 : (((u128)1 << pbits) - 1);

#pragma omp parallel
    {
        int tid=0, nthr=1;
#ifdef _OPENMP
        tid = omp_get_thread_num(); nthr = omp_get_num_threads();
#endif
        pin_thread_to_core(tid);

        u64 current_time = (u64)time(NULL);
        Rng rng = {
            .s0 = 7364529176530163ULL ^ ((u64)tid * 6364136223846793005ULL) ^ (u64)p ^ current_time,
            .s1 = 1442695040888963407ULL ^ ((u64)(tid+1) * 2862933555777941757ULL) ^ (current_time << 32)
        };
        for (int i=0;i<200;i++) rng64(&rng);

        u64 budget = max_trials / nthr + 1;
        u64 A_trials = 0;

        u128 pending_A = 0, pending_xP = 0;
        int have_pending_A = 0;

        while (found_count < target_total && A_trials < budget) {
            u128 A;
            u128 xP16 = 0;

            if (x16_mode && x16_ready) {
                x16_montgomery_A128(&A, &xP16, &pending_A, &pending_xP, &have_pending_A,
                                     &rng, p, rand_mask16, sqrtm1_16, &mt16);
            } else {
                A = (u128)rng64(&rng) | ((u128)rng64(&rng) << 64); A %= p;
                if (A==2 || A==p-2) continue;
            }

            A_trials++;
            thread_trials[tid * 8] = A_trials;

            int a_success = 0;
            int max_x_tries = (x16_mode && x16_ready) ? 1 : 50;

            for (int x_tries = 0; x_tries < max_x_tries && !a_success && found_count < target_total; x_tries++) {
                u128 x0r;
                if (x16_mode && x16_ready) {
                    x0r = xP16;
                } else {
                    x0r = (u128)rng64(&rng) | ((u128)rng64(&rng) << 64); x0r %= p;
                    if (x0r < 2) continue;
                }

                t_ops.mul_count = 0; t_ops.inv_count = 0; t_ops.branch_count = 0;
                uint64_t c_start = rdtsc();
                int ok = verify128(p, A, x0r);
                uint64_t c_end = rdtsc();
                uint64_t cyc = c_end - c_start;
                OpCounter this_ops = t_ops;

                if (ok) {
#pragma omp critical
                    {
                        int is_dup = 0;
                        for(int j = 0; j < found_count; j++) {
                            if(out_A[j] == A && out_x0[j] == x0r) {
                                is_dup = 1; break;
                            }
                        }
                        if (!is_dup && found_count < target_total) {
                            out_A[found_count] = A;
                            out_x0[found_count] = x0r;

                            u64 exact_total_trials = 0;
                            for (int t = 0; t < nthr; t++) {
                                exact_total_trials += thread_trials[t * 8];
                            }
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
 * Batch-processing main
 * ================================================================ */

int main(int argc, char *argv[]) {
    if (argc < 5) {
        printf("Usage: ./pomerance <stateful_input.txt> <output_pure.csv> <output_metrics.csv> <target_total> [x16]\n");
        return 1;
    }

    int target_total = atoi(argv[4]);
    if (target_total <= 0) {
        printf("Error: <target_total> must be greater than 0.\n");
        return 1;
    }

    if (argc >= 6 && strcmp(argv[5], "x16") == 0) {
        g_x16_mode = 1;
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

    printf("Sampling mode: %s\n",
           g_x16_mode ? "X1(16) prescribed torsion (Tate normal form)" : "brute-force random A");

    char line[65536];
    unsigned long long current_index = 0;

    while (fgets(line, sizeof(line), input)) {
        char *ptr = line;
        u128 p = parse128_adv(&ptr);
        if (p == 0) continue;

        current_index++;
        int num_existing = (int)parse128_adv(&ptr);

        /* Remove the forced 128-bit routing for x16 mode so it can run via search64. */
        int use_u128_path = (p >= ((u128)1 << 63));

        for(int i = 0; i < num_existing; i++) {
            if (!use_u128_path) {
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

        char p_display[50]; sprint128(p_display, p);
        printf("[%llu] Prime: %s (Existing: %d, Target: %d)...\n", current_index, p_display, num_existing, target_total);
        fflush(stdout);

        double t_start = now_sec();
        int found_amount = 0;

        if (!use_u128_path) {
            found_amount = search64((u64)p, target_total, num_existing, out_A_arr, out_x0_arr, out_trials_arr,
                                     out_cycles_arr, out_mul_arr, out_inv_arr, out_branch_arr, g_x16_mode);
        } else {
            found_amount = search128(p, target_total, num_existing, out_A128, out_x0128, out_trials_arr,
                                      out_cycles_arr, out_mul_arr, out_inv_arr, out_branch_arr, g_x16_mode);
            for(int i = num_existing; i < found_amount; i++) {
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
                    if (out_trials_arr[j] > out_trials_arr[j+1]) {
                        u64 tmp;
                        tmp = out_trials_arr[j]; out_trials_arr[j] = out_trials_arr[j+1]; out_trials_arr[j+1] = tmp;

                        if (!use_u128_path) {
                            tmp = out_A_arr[j]; out_A_arr[j] = out_A_arr[j+1]; out_A_arr[j+1] = tmp;
                            tmp = out_x0_arr[j]; out_x0_arr[j] = out_x0_arr[j+1]; out_x0_arr[j+1] = tmp;
                        } else {
                            u128 tmp128;
                            tmp128 = out_A128[j]; out_A128[j] = out_A128[j+1]; out_A128[j+1] = tmp128;
                            tmp128 = out_x0128[j]; out_x0128[j] = out_x0128[j+1]; out_x0128[j+1] = tmp128;
                            out_A_arr[j] = (u64)out_A128[j]; out_A_arr[j+1] = (u64)out_A128[j+1];
                            out_x0_arr[j] = (u64)out_x0128[j]; out_x0_arr[j+1] = (u64)out_x0128[j+1];
                        }

                        tmp = out_cycles_arr[j]; out_cycles_arr[j] = out_cycles_arr[j+1]; out_cycles_arr[j+1] = tmp;
                        tmp = out_mul_arr[j]; out_mul_arr[j] = out_mul_arr[j+1]; out_mul_arr[j+1] = tmp;
                        tmp = out_inv_arr[j]; out_inv_arr[j] = out_inv_arr[j+1]; out_inv_arr[j+1] = tmp;
                        tmp = out_branch_arr[j]; out_branch_arr[j] = out_branch_arr[j+1]; out_branch_arr[j+1] = tmp;
                    }
                }
            }

            for(int i = num_existing; i < found_amount; i++) {
                double opc = opcount_m3(&(OpCounter){out_mul_arr[i], out_inv_arr[i], out_branch_arr[i]});
                char p_str[50], a_str[50], x0_str[50];
                sprint128(p_str, p);

                /* Fixed: safely printing variables even if they overflow 64-bit on the u128 path */
                if (use_u128_path) {
                    sprint128(a_str, out_A128[i]);
                    sprint128(x0_str, out_x0128[i]);
                } else {
                    sprint128(a_str, out_A_arr[i]);
                    sprint128(x0_str, out_x0_arr[i]);
                }

                fprintf(out_pure, "%s,%s\n", a_str, x0_str);
                fprintf(out_metrics, "%s,%s,%s,%llu,%llu,%llu,%llu,%llu,%.2f,%.2f\n",
                        p_str, a_str, x0_str,
                        (unsigned long long)out_trials_arr[i],
                        (unsigned long long)out_cycles_arr[i],
                        (unsigned long long)out_mul_arr[i],
                        (unsigned long long)out_inv_arr[i],
                        (unsigned long long)out_branch_arr[i],
                        opc,
                        batch_time_ms);
            }
            fflush(out_pure);
            fflush(out_metrics);
        }
    }

    free(out_A_arr); free(out_x0_arr); free(out_trials_arr);
    free(out_cycles_arr); free(out_mul_arr); free(out_inv_arr); free(out_branch_arr);
    free(out_A128); free(out_x0128);
    fclose(input); fclose(out_pure); fclose(out_metrics);

    return 0;
}
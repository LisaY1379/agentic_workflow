```cpp
/*
 * pomerance_cuda_u96.cu
 *
 * Production CUDA implementation of the X1(16) nonsplit first-branch
 * successive-halving Pomerance-triple search.
 *
 * This is intentionally specialized for the high-throughput production path:
 *
 *   x16halvenonsplit
 *
 * Supported:
 *   - Odd primes p < 2^96
 *   - p == 5 mod 8
 *   - p == 3 mod 4
 *   - X1(16) construction
 *   - Nonsplit discriminant y-filter
 *   - First-branch halving chain
 *   - 96-bit Montgomery arithmetic
 *   - GPU work claiming and first-hit publication
 *
 * Compile for Ada:
 *
 *   nvcc -O3 -std=c++17 -arch=sm_89 -Xptxas=-v \
 *       -o pomerance_cuda_u96 pomerance_cuda_u96.cu
 *
 * Example:
 *
 *   ./pomerance_cuda_u96 100000000000000000000000067 121 \
 *       550000000000 x16halvenonsplit 1000000000
 *
 * Arguments:
 *
 *   ./pomerance_cuda_u96 <p> [seed_offset] [max_trials]
 *       [x16halvenonsplit] [chunk_trials] [blocks] [threads] [claim_batch]
 *
 * The default launch geometry is:
 *
 *   blocks  = SM_count * 8
 *   threads = 128
 *   claim_batch = 1
 */

#include <cuda_runtime.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

using u32  = unsigned int;
using u64  = unsigned long long;
using u128 = unsigned __int128;

#define HD  __host__ __device__
#define DEV __device__

struct U96 {
    u32 v[3];
};

struct U128Parts {
    u64 lo;
    u64 hi;
};

struct Field96 {
    U96 p;
    U96 r2;
    U96 one;
    u32 nprime;
};

struct SearchParams96 {
    Field96 f;

    U96 rand_mask;
    U96 sqrtm1_m;

    U96 inv2_m;
    U96 inv4_m;
    U96 two_m;
    U96 four_m;
    U96 eight_m;
    U96 c24_m;
    U96 c32_m;
    U96 c48_m;

    u128 leg_exp;
    u128 sqrt_exp;
    u128 inv_exp;

    int k;
    int sqrt_case; /* 3: p == 3 mod 4; 5: p == 5 mod 8 */
};

struct Rng {
    u64 s0;
    u64 s1;
};

struct GpuStats {
    u64 claimed;
    u64 processed;
};

struct GpuResult {
    int found;
    u64 trial_index;
    U128Parts A;
    U128Parts x0;
};

/* ------------------------------------------------------------------------- */
/* Utility                                                                    */
/* ------------------------------------------------------------------------- */

HD static inline U96 u96_zero() {
    U96 x{{0, 0, 0}};
    return x;
}

HD static inline U96 u96_from_u64(u64 x) {
    U96 out{{(u32)x, (u32)(x >> 32), 0}};
    return out;
}

HD static inline U96 u96_from_u128(u128 x) {
    U96 out{{(u32)x, (u32)(x >> 32), (u32)(x >> 64)}};
    return out;
}

HD static inline u128 u96_to_u128(U96 x) {
    return (u128)x.v[0] |
           ((u128)x.v[1] << 32) |
           ((u128)x.v[2] << 64);
}

HD static inline U128Parts pack96(U96 x) {
    u128 v = u96_to_u128(x);
    U128Parts out;
    out.lo = (u64)v;
    out.hi = (u64)(v >> 64);
    return out;
}

static inline u128 unpack128(U128Parts x) {
    return ((u128)x.hi << 64) | (u128)x.lo;
}

HD static inline int cmp96(U96 a, U96 b) {
    if (a.v[2] != b.v[2]) return a.v[2] < b.v[2] ? -1 : 1;
    if (a.v[1] != b.v[1]) return a.v[1] < b.v[1] ? -1 : 1;
    if (a.v[0] != b.v[0]) return a.v[0] < b.v[0] ? -1 : 1;
    return 0;
}

HD static inline int eq96(U96 a, U96 b) {
    return ((a.v[0] ^ b.v[0]) |
            (a.v[1] ^ b.v[1]) |
            (a.v[2] ^ b.v[2])) == 0;
}

HD static inline int is_zero96(U96 a) {
    return (a.v[0] | a.v[1] | a.v[2]) == 0;
}

HD static inline U96 and96(U96 a, U96 b) {
    U96 out{{a.v[0] & b.v[0], a.v[1] & b.v[1], a.v[2] & b.v[2]}};
    return out;
}

HD static inline U96 add96_raw(U96 a, U96 b, u32 *carry_out = nullptr) {
    U96 out;
    u64 t = (u64)a.v[0] + b.v[0];
    out.v[0] = (u32)t;

    t = (u64)a.v[1] + b.v[1] + (t >> 32);
    out.v[1] = (u32)t;

    t = (u64)a.v[2] + b.v[2] + (t >> 32);
    out.v[2] = (u32)t;

    if (carry_out) *carry_out = (u32)(t >> 32);
    return out;
}

HD static inline U96 sub96_raw(U96 a, U96 b, u32 *borrow_out = nullptr) {
    U96 out;

    u64 bi = b.v[0];
    out.v[0] = (u32)((u64)a.v[0] - bi);
    u64 borrow = ((u64)a.v[0] < bi) ? 1 : 0;

    bi = (u64)b.v[1] + borrow;
    out.v[1] = (u32)((u64)a.v[1] - bi);
    borrow = ((u64)a.v[1] < bi) ? 1 : 0;

    bi = (u64)b.v[2] + borrow;
    out.v[2] = (u32)((u64)a.v[2] - bi);
    borrow = ((u64)a.v[2] < bi) ? 1 : 0;

    if (borrow_out) *borrow_out = (u32)borrow;
    return out;
}

HD static inline U96 addmod96(U96 a, U96 b, U96 p) {
    u32 carry;
    U96 r = add96_raw(a, b, &carry);
    if (carry || cmp96(r, p) >= 0) r = sub96_raw(r, p);
    return r;
}

HD static inline U96 submod96(U96 a, U96 b, U96 p) {
    u32 borrow;
    U96 r = sub96_raw(a, b, &borrow);
    if (borrow) r = add96_raw(r, p);
    return r;
}

/* ------------------------------------------------------------------------- */
/* 96-bit Montgomery arithmetic                                               */
/* ------------------------------------------------------------------------- */

HD static inline u32 inv32_odd(u32 p0) {
    u32 x = 1;
#pragma unroll
    for (int i = 0; i < 5; ++i) x *= 2u - p0 * x;
    return 0u - x;
}

/*
 * CIOS-style fixed-width Montgomery multiplication.
 *
 * Inputs and output are reduced mod p.  The Montgomery radix is 2^96.
 */
HD static inline U96 mont_mul96(U96 a, U96 b, const Field96 *f) {
    const u32 a0 = a.v[0], a1 = a.v[1], a2 = a.v[2];
    const u32 b0 = b.v[0], b1 = b.v[1], b2 = b.v[2];
    const u32 p0 = f->p.v[0], p1 = f->p.v[1], p2 = f->p.v[2];
    const u32 np = f->nprime;

    u64 uv, carry;
    u32 t0, t1, t2, t3, t4, t5, t6;

    uv = (u64)a0 * b0;
    t0 = (u32)uv;
    carry = uv >> 32;

    uv = (u64)a0 * b1 + carry;
    t1 = (u32)uv;
    carry = uv >> 32;

    uv = (u64)a0 * b2 + carry;
    t2 = (u32)uv;
    t3 = (u32)(uv >> 32);

    uv = (u64)t1 + (u64)a1 * b0;
    t1 = (u32)uv;
    carry = uv >> 32;

    uv = (u64)t2 + (u64)a1 * b1 + carry;
    t2 = (u32)uv;
    carry = uv >> 32;

    uv = (u64)t3 + (u64)a1 * b2 + carry;
    t3 = (u32)uv;
    t4 = (u32)(uv >> 32);

    uv = (u64)t2 + (u64)a2 * b0;
    t2 = (u32)uv;
    carry = uv >> 32;

    uv = (u64)t3 + (u64)a2 * b1 + carry;
    t3 = (u32)uv;
    carry = uv >> 32;

    uv = (u64)t4 + (u64)a2 * b2 + carry;
    t4 = (u32)uv;
    t5 = (u32)(uv >> 32);
    t6 = 0;

    u32 m = t0 * np;

    uv = (u64)t0 + (u64)m * p0;
    carry = uv >> 32;

    uv = (u64)t1 + (u64)m * p1 + carry;
    t1 = (u32)uv;
    carry = uv >> 32;

    uv = (u64)t2 + (u64)m * p2 + carry;
    t2 = (u32)uv;
    carry = uv >> 32;

    uv = (u64)t3 + carry;
    t3 = (u32)uv;
    carry = uv >> 32;

    uv = (u64)t4 + carry;
    t4 = (u32)uv;
    carry = uv >> 32;

    uv = (u64)t5 + carry;
    t5 = (u32)uv;
    t6 = (u32)(uv >> 32);

    m = t1 * np;

    uv = (u64)t1 + (u64)m * p0;
    carry = uv >> 32;

    uv = (u64)t2 + (u64)m * p1 + carry;
    t2 = (u32)uv;
    carry = uv >> 32;

    uv = (u64)t3 + (u64)m * p2 + carry;
    t3 = (u32)uv;
    carry = uv >> 32;

    uv = (u64)t4 + carry;
    t4 = (u32)uv;
    carry = uv >> 32;

    uv = (u64)t5 + carry;
    t5 = (u32)uv;
    carry = uv >> 32;

    uv = (u64)t6 + carry;
    t6 = (u32)uv;

    m = t2 * np;

    uv = (u64)t2 + (u64)m * p0;
    carry = uv >> 32;

    uv = (u64)t3 + (u64)m * p1 + carry;
    t3 = (u32)uv;
    carry = uv >> 32;

    uv = (u64)t4 + (u64)m * p2 + carry;
    t4 = (u32)uv;
    carry = uv >> 32;

    uv = (u64)t5 + carry;
    t5 = (u32)uv;
    carry = uv >> 32;

    uv = (u64)t6 + carry;
    t6 = (u32)uv;

    U96 out{{t3, t4, t5}};
    if (t6 || cmp96(out, f->p) >= 0) out = sub96_raw(out, f->p);
    return out;
}

HD static inline U96 to_mont96(U96 x, const Field96 *f) {
    return mont_mul96(x, f->r2, f);
}

HD static inline U96 from_mont96(U96 x, const Field96 *f) {
    return mont_mul96(x, u96_from_u64(1), f);
}

HD static U96 powmod96_mont(U96 a_m, u128 exponent, const Field96 *f) {
    U96 r = f->one;
    U96 b = a_m;

    while (exponent) {
        if (exponent & 1) r = mont_mul96(r, b, f);
        exponent >>= 1;
        if (exponent) b = mont_mul96(b, b, f);
    }
    return r;
}

HD static inline U96 invert96_mont(U96 a_m, const SearchParams96 *params) {
    return powmod96_mont(a_m, params->inv_exp, &params->f);
}

HD static int invert_batch2_96_mont(
    U96 out[2],
    const U96 in[2],
    const SearchParams96 *params
) {
    if (is_zero96(in[0]) || is_zero96(in[1])) return 0;

    U96 product = mont_mul96(in[0], in[1], &params->f);
    U96 product_inv = invert96_mont(product, params);

    out[0] = mont_mul96(product_inv, in[1], &params->f);
    out[1] = mont_mul96(product_inv, in[0], &params->f);
    return 1;
}

/* ------------------------------------------------------------------------- */
/* Square roots                                                               */
/* ------------------------------------------------------------------------- */

HD static int sqrtmod_fast96_mont(
    U96 *root_m,
    U96 n_m,
    const SearchParams96 *params
) {
    if (is_zero96(n_m)) {
        *root_m = u96_zero();
        return 1;
    }

    U96 x_m = powmod96_mont(n_m, params->sqrt_exp, &params->f);

    if (eq96(mont_mul96(x_m, x_m, &params->f), n_m)) {
        *root_m = x_m;
        return 1;
    }

    if (params->sqrt_case == 3) return 0;

    x_m = mont_mul96(x_m, params->sqrtm1_m, &params->f);
    if (eq96(mont_mul96(x_m, x_m, &params->f), n_m)) {
        *root_m = x_m;
        return 1;
    }

    return 0;
}

/* ------------------------------------------------------------------------- */
/* Montgomery x-only curve arithmetic                                         */
/* ------------------------------------------------------------------------- */

HD static inline void xDBL96(
    U96 *Xo,
    U96 *Zo,
    U96 X,
    U96 Z,
    U96 a24_m,
    const Field96 *f
) {
    U96 u = addmod96(X, Z, f->p);
    U96 v = submod96(X, Z, f->p);

    u = mont_mul96(u, u, f);
    v = mont_mul96(v, v, f);

    *Xo = mont_mul96(u, v, f);

    U96 w = submod96(u, v, f);
    U96 q = addmod96(v, mont_mul96(a24_m, w, f), f->p);
    *Zo = mont_mul96(w, q, f);
}

HD static int verify96_mont(
    U96 A_m,
    U96 x0_m,
    const SearchParams96 *params
) {
    U96 neg_two_m = submod96(u96_zero(), params->two_m, params->f.p);

    if (eq96(A_m, params->two_m) || eq96(A_m, neg_two_m)) return 0;

    U96 a24_m = mont_mul96(
        addmod96(A_m, params->two_m, params->f.p),
        params->inv4_m,
        &params->f
    );

    U96 X = x0_m;
    U96 Z = params->f.one;

    for (int i = 1; i <= params->k; ++i) {
        xDBL96(&X, &Z, X, Z, a24_m, &params->f);

        if (i < params->k && is_zero96(Z)) return 0;
        if (i == params->k && !is_zero96(Z)) return 0;
    }

    return 1;
}

/* ------------------------------------------------------------------------- */
/* X1(16) sampler                                                             */
/* ------------------------------------------------------------------------- */

HD static U96 x16_A_numerator_from_y96_mont(
    U96 y_m,
    const SearchParams96 *params
) {
    const Field96 *f = &params->f;
    U96 n = f->one;

    n = submod96(mont_mul96(n, y_m, f), params->eight_m, f->p);
    n = addmod96(mont_mul96(n, y_m, f), params->c24_m, f->p);
    n = submod96(mont_mul96(n, y_m, f), params->c32_m, f->p);
    n = addmod96(mont_mul96(n, y_m, f), params->eight_m, f->p);
    n = addmod96(mont_mul96(n, y_m, f), params->c32_m, f->p);
    n = submod96(mont_mul96(n, y_m, f), params->c48_m, f->p);
    n = addmod96(mont_mul96(n, y_m, f), params->c32_m, f->p);
    n = submod96(mont_mul96(n, y_m, f), params->eight_m, f->p);

    return n;
}

HD static int x16_y_predicts_nonsplit96_mont(
    U96 y_m,
    U96 y2_m,
    const SearchParams96 *params
) {
    const Field96 *f = &params->f;

    U96 f1_m = submod96(y2_m, params->two_m, f->p);
    U96 two_y_m = addmod96(y_m, y_m, f->p);
    U96 four_y_m = addmod96(two_y_m, two_y_m, f->p);

    U96 f2_m = addmod96(
        submod96(y2_m, four_y_m, f->p),
        params->two_m,
        f->p
    );

    U96 disc_m = mont_mul96(f1_m, f2_m, f);
    if (is_zero96(disc_m)) return 0;

    U96 legendre_m = powmod96_mont(disc_m, params->leg_exp, f);
    return !eq96(legendre_m, f->one);
}

HD static int x16_root_to_montgomery_A96_mont(
    U96 *A_standard,
    U96 *A_m,
    U96 *xP16_m,
    U96 x_m,
    U96 y_m,
    const SearchParams96 *params
) {
    const Field96 *f = &params->f;

    U96 numerator_m = x16_A_numerator_from_y96_mont(y_m, params);

    U96 ym1_m = submod96(y_m, f->one, f->p);
    U96 ym1_squared_m = mont_mul96(ym1_m, ym1_m, f);

    U96 den_A_m = mont_mul96(
        params->four_m,
        mont_mul96(ym1_squared_m, ym1_squared_m, f),
        f
    );

    U96 den_x_m = submod96(x_m, y_m, f->p);

    U96 values[2] = {den_A_m, den_x_m};
    U96 inverses[2];

    if (!invert_batch2_96_mont(inverses, values, params)) return 0;

    U96 candidate_A_m = mont_mul96(numerator_m, inverses[0], f);
    U96 candidate_x_m = mont_mul96(x_m, inverses[1], f);
    U96 candidate_A = from_mont96(candidate_A_m, f);

    U96 two = u96_from_u64(2);
    U96 p_minus_two = submod96(f->p, two, f->p);

    if (cmp96(candidate_A, two) <= 0) return 0;
    if (cmp96(candidate_A, p_minus_two) >= 0) return 0;

    *A_standard = candidate_A;
    *A_m = candidate_A_m;
    *xP16_m = candidate_x_m;
    return 1;
}

/* ------------------------------------------------------------------------- */
/* Halving                                                                    */
/* ------------------------------------------------------------------------- */

HD static int halve_once_first96_mont(
    U96 *xout_m,
    U96 A_m,
    U96 x_m,
    const SearchParams96 *params
) {
    const Field96 *f = &params->f;

    U96 x2_m = mont_mul96(x_m, x_m, f);

    U96 d_m = addmod96(
        addmod96(x2_m, mont_mul96(A_m, x_m, f), f->p),
        f->one,
        f->p
    );

    U96 sqrt_d_m;
    if (!sqrtmod_fast96_mont(&sqrt_d_m, d_m, params)) return 0;

    U96 roots_d[2] = {
        sqrt_d_m,
        submod96(u96_zero(), sqrt_d_m, f->p)
    };

#pragma unroll
    for (int i = 0; i < 2; ++i) {
        U96 u_m = addmod96(
            addmod96(x_m, x_m, f->p),
            addmod96(roots_d[i], roots_d[i], f->p),
            f->p
        );

        U96 w_m = submod96(
            mont_mul96(u_m, u_m, f),
            params->four_m,
            f->p
        );

        U96 sqrt_w_m;
        if (!sqrtmod_fast96_mont(&sqrt_w_m, w_m, params)) continue;

        U96 candidates[2] = {
            mont_mul96(addmod96(u_m, sqrt_w_m, f->p), params->inv2_m, f),
            mont_mul96(submod96(u_m, sqrt_w_m, f->p), params->inv2_m, f)
        };

#pragma unroll
        for (int j = 0; j < 2; ++j) {
            if (!is_zero96(candidates[j])) {
                *xout_m = candidates[j];
                return 1;
            }
        }
    }

    return 0;
}

HD static int halve_chain_from_depth96_mont(
    U96 *xout_m,
    U96 A_m,
    U96 x_m,
    int depth,
    const SearchParams96 *params
) {
    for (int d = depth; d < params->k; ++d) {
        if (!halve_once_first96_mont(&x_m, A_m, x_m, params)) return 0;
    }

    if (!verify96_mont(A_m, x_m, params)) return 0;

    *xout_m = x_m;
    return 1;
}

/* ------------------------------------------------------------------------- */
/* Random generator                                                           */
/* ------------------------------------------------------------------------- */

HD static inline u64 rng64(Rng *r) {
    u64 s1 = r->s0;
    const u64 s0 = r->s1;

    r->s0 = s0;
    s1 ^= s1 << 23;
    r->s1 = s1 ^ s0 ^ (s1 >> 17) ^ (s0 >> 26);

    return r->s1 + s0;
}

HD static inline u64 splitmix64_next(u64 *x) {
    u64 z = (*x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

DEV static U96 rand_below96(Rng *rng, const SearchParams96 *params) {
    for (;;) {
        u64 a = rng64(rng);
        u64 b = rng64(rng);

        U96 x{{(u32)a, (u32)(a >> 32), (u32)b}};
        x = and96(x, params->rand_mask);

        if (cmp96(x, params->f.p) < 0) return x;
    }
}

DEV static inline int result_found(const GpuResult *result) {
    return *((volatile const int *)&result->found);
}

/* ------------------------------------------------------------------------- */
/* CUDA search kernel                                                         */
/* ------------------------------------------------------------------------- */

__global__ static void x16halvenonsplit_kernel96(
    SearchParams96 params,
    u64 seed_offset,
    u64 chunk_nonce,
    u64 max_curves,
    u64 global_base,
    u64 claim_batch,
    GpuStats *__restrict__ stats,
    GpuResult *__restrict__ result
) {
    const u64 tid =
        (u64)blockIdx.x * (u64)blockDim.x + (u64)threadIdx.x;

    u64 seed =
        seed_offset ^
        (chunk_nonce * 0xd1342543de82ef95ULL) ^
        (tid * 0x9e3779b97f4a7c15ULL) ^
        0x7364529176530163ULL;

    Rng rng;
    rng.s0 = splitmix64_next(&seed);
    rng.s1 = splitmix64_next(&seed);

    if ((rng.s0 | rng.s1) == 0) {
        rng.s1 = 1442695040888963407ULL;
    }

#pragma unroll 1
    for (int i = 0; i < 200; ++i) {
        (void)rng64(&rng);
    }

    u64 local_processed = 0;
    u64 next_claim = 0;
    u64 claim_limit = 0;

    while (!result_found(result)) {
        U96 y = rand_below96(&rng, &params);
        if (is_zero96(y)) continue;

        U96 y_m = to_mont96(y, &params.f);
        U96 y2_m = mont_mul96(y_m, y_m, &params.f);

        if (!x16_y_predicts_nonsplit96_mont(y_m, y2_m, &params)) continue;

        U96 y3_m = mont_mul96(y2_m, y_m, &params.f);

        U96 two_y_m = addmod96(y_m, y_m, params.f.p);
        U96 qa_m = submod96(y2_m, two_y_m, params.f.p);
        if (is_zero96(qa_m)) continue;

        U96 qb_m = submod96(
            addmod96(y2_m, y2_m, params.f.p),
            y3_m,
            params.f.p
        );

        U96 qc_m = submod96(params.f.one, y_m, params.f.p);

        U96 D_m = submod96(
            mont_mul96(qb_m, qb_m, &params.f),
            mont_mul96(
                addmod96(qa_m, qa_m, params.f.p),
                addmod96(qc_m, qc_m, params.f.p),
                &params.f
            ),
            params.f.p
        );

        U96 sqrt_D_m;
        if (!sqrtmod_fast96_mont(&sqrt_D_m, D_m, &params)) continue;

        U96 two_qa_m = addmod96(qa_m, qa_m, params.f.p);
        U96 inv_two_qa_m = invert96_mont(two_qa_m, &params);

        U96 roots_m[2] = {
            mont_mul96(
                submod96(sqrt_D_m, qb_m, params.f.p),
                inv_two_qa_m,
                &params.f
            ),
            mont_mul96(
                submod96(
                    submod96(u96_zero(), sqrt_D_m, params.f.p),
                    qb_m,
                    params.f.p
                ),
                inv_two_qa_m,
                &params.f
            )
        };

        int stop = 0;

#pragma unroll
        for (int root_index = 0; root_index < 2; ++root_index) {
            U96 A;
            U96 A_m;
            U96 xP16_m;

            if (!x16_root_to_montgomery_A96_mont(
                    &A, &A_m, &xP16_m,
                    roots_m[root_index], y_m, &params)) {
                continue;
            }

            if (next_claim >= claim_limit) {
                next_claim = atomicAdd(&stats->claimed, claim_batch);
                claim_limit = next_claim + claim_batch;
            }

            const u64 claim = next_claim++;

            if (claim >= max_curves || result_found(result)) {
                stop = 1;
                break;
            }

            ++local_processed;

            U96 x_result_m;
            if (halve_chain_from_depth96_mont(
                    &x_result_m, A_m, xP16_m, 4, &params)) {
                if (atomicCAS(&result->found, 0, 1) == 0) {
                    result->trial_index = global_base + claim;
                    result->A = pack96(A);
                    result->x0 = pack96(from_mont96(x_result_m, &params.f));
                }

                stop = 1;
                break;
            }
        }

        if (stop) break;
    }

    atomicAdd(&stats->processed, local_processed);
}

/* ------------------------------------------------------------------------- */
/* Host initialization                                                        */
/* ------------------------------------------------------------------------- */

static void cuda_check(cudaError_t err, const char *where) {
    if (err == cudaSuccess) return;
    std::fprintf(stderr, "CUDA failure at %s: %s\n",
                 where, cudaGetErrorString(err));
    std::exit(1);
}

static u128 parse128(const char *s) {
    u128 out = 0;
    while (*s >= '0' && *s <= '9') {
        out = out * 10 + (u128)(*s - '0');
        ++s;
    }
    return out;
}

static std::string sprint128(u128 x) {
    if (!x) return "0";

    char text[50];
    int i = 49;
    text[i] = 0;

    while (x) {
        text[--i] = (char)('0' + (x % 10));
        x /= 10;
    }

    return std::string(text + i);
}

static int bitlen128(u128 x) {
    u64 high = (u64)(x >> 64);
    if (high) return 64 + 64 - __builtin_clzll(high);

    u64 low = (u64)x;
    return low ? 64 - __builtin_clzll(low) : 0;
}

static U96 u96_mask_bits(int bits) {
    U96 out{{0, 0, 0}};

    if (bits >= 96) {
        out.v[0] = out.v[1] = out.v[2] = 0xffffffffu;
        return out;
    }

    for (int i = 0; i < 3; ++i) {
        int n = bits - i * 32;

        if (n >= 32) out.v[i] = 0xffffffffu;
        else if (n > 0) out.v[i] = (1u << n) - 1u;
        else out.v[i] = 0;
    }

    return out;
}

static Field96 field96_init(u128 p) {
    Field96 f{};
    f.p = u96_from_u128(p);
    f.nprime = inv32_odd(f.p.v[0]);

    U96 r = u96_from_u64(1);

    for (int i = 0; i < 96; ++i) {
        r = addmod96(r, r, f.p);
    }

    f.one = r;

    for (int i = 0; i < 96; ++i) {
        r = addmod96(r, r, f.p);
    }

    f.r2 = r;
    return f;
}

static int compute_k(u128 p) {
    u64 q = (u64)sqrt((long double)p);

    while ((u128)(q + 1) * (q + 1) <= p) ++q;
    while ((u128)q * q > p) --q;

    u64 sq = (u64)sqrt((long double)q);

    while ((sq + 1) * (sq + 1) <= q) ++sq;
    while (sq * sq > q) --sq;

    const u64 bound = q + 1 + 2 * sq;

    int k = 0;
    u64 power = 1;

    while (power <= bound) {
        ++k;
        power <<= 1;
    }

    return k;
}

static SearchParams96 search_params96_init(u128 p, int k) {
    SearchParams96 params{};

    params.f = field96_init(p);

    const int pbits = bitlen128(p);
    params.rand_mask = u96_mask_bits(pbits);

    params.leg_exp = (p - 1) >> 1;
    params.inv_exp = p - 2;

    params.sqrt_case = ((u64)(p & 3) == 3) ? 3 : 5;
    params.sqrt_exp =
        params.sqrt_case == 3 ? ((p + 1) >> 2) : ((p + 3) >> 3);

    params.k = k;

    params.two_m = to_mont96(u96_from_u64(2), &params.f);
    params.four_m = to_mont96(u96_from_u64(4), &params.f);
    params.eight_m = to_mont96(u96_from_u64(8), &params.f);
    params.c24_m = to_mont96(u96_from_u64(24), &params.f);
    params.c32_m = to_mont96(u96_from_u64(32), &params.f);
    params.c48_m = to_mont96(u96_from_u64(48), &params.f);

    params.inv2_m = to_mont96(
        u96_from_u128((p + 1) >> 1),
        &params.f
    );

    params.inv4_m = powmod96_mont(
        params.four_m,
        params.inv_exp,
        &params.f
    );

    params.sqrtm1_m =
        params.sqrt_case == 5
            ? powmod96_mont(params.two_m, (p - 1) >> 2, &params.f)
            : u96_zero();

    return params;
}

static u64 parse_u64(const char *s, u64 fallback) {
    if (!s || !*s) return fallback;
    return std::strtoull(s, nullptr, 10);
}

static void usage(const char *name) {
    std::fprintf(
        stderr,
        "Usage:\n"
        "  %s <p> [seed_offset] [max_trials] [x16halvenonsplit]\n"
        "     [chunk_trials] [blocks] [threads] [claim_batch]\n",
        name
    );
}

/* ------------------------------------------------------------------------- */
/* Main                                                                       */
/* ------------------------------------------------------------------------- */

int main(int argc, char **argv) {
    if (argc < 2) {
        usage(argv[0]);
        return 1;
    }

    const u128 p = parse128(argv[1]);

    if (p == 0 || bitlen128(p) > 96) {
        std::fprintf(stderr, "This u96 CUDA backend requires 0 < p < 2^96.\n");
        return 1;
    }

    const u64 pmod8 = (u64)(p & 7);

    if (pmod8 != 5 && ((u64)(p & 3) != 3)) {
        std::fprintf(
            stderr,
            "Supported square-root cases are p == 5 mod 8 and p == 3 mod 4.\n"
        );
        return 1;
    }

    const char *mode = argc >= 5 ? argv[4] : "x16halvenonsplit";

    if (std::strcmp(mode, "x16halvenonsplit") != 0) {
        std::fprintf(
            stderr,
            "Only x16halvenonsplit is implemented by this production CUDA path.\n"
        );
        return 1;
    }

    const u64 seed_offset = argc >= 3 ? parse_u64(argv[2], 0) : 0;
    u64 max_trials = argc >= 4 ? parse_u64(argv[3], 1000000ULL) : 1000000ULL;

    if (!max_trials) max_trials = 1000000ULL;

    const u64 default_chunk =
        max_trials < 1000000000ULL ? max_trials : 1000000000ULL;

    u64 chunk_trials =
        argc >= 6 ? parse_u64(argv[5], default_chunk) : default_chunk;

    if (!chunk_trials) chunk_trials = max_trials;

    cuda_check(cudaSetDevice(0), "cudaSetDevice");

    cudaDeviceProp prop{};
    cuda_check(
        cudaGetDeviceProperties(&prop, 0),
        "cudaGetDeviceProperties"
    );

    int blocks =
        argc >= 7 ? (int)parse_u64(argv[6], 0) : prop.multiProcessorCount * 8;

    int threads =
        argc >= 8 ? (int)parse_u64(argv[7], 128) : 128;

    u64 claim_batch =
        argc >= 9 ? parse_u64(argv[8], 1) : 1;

    if (blocks <= 0) blocks = prop.multiProcessorCount * 8;
    if (threads <= 0) threads = 128;
    if (!claim_batch) claim_batch = 1;

    const int k = compute_k(p);
    const SearchParams96 params = search_params96_init(p, k);

    std::printf("Pomerance CUDA X1(16) nonsplit halving search\n\n");
    std::printf("p = %s\n", sprint128(p).c_str());
    std::printf("seed_offset = %llu\n", seed_offset);
    std::printf("max_trials = %llu\n", max_trials);
    std::printf("chunk_trials = %llu\n", chunk_trials);
    std::printf("k = %d\n", k);
    std::printf("backend = u96-montgomery\n");
    std::printf(
        "GPU = %s  SMs=%d  blocks=%d  threads=%d  claim_batch=%llu\n\n",
        prop.name,
        prop.multiProcessorCount,
        blocks,
        threads,
        claim_batch
    );

    GpuStats *d_stats = nullptr;
    GpuResult *d_result = nullptr;

    cuda_check(cudaMalloc(&d_stats, sizeof(GpuStats)), "cudaMalloc(stats)");
    cuda_check(cudaMalloc(&d_result, sizeof(GpuResult)), "cudaMalloc(result)");

    u64 total = 0;
    u64 chunk_nonce = 0;
    GpuResult final_result{};

    const auto start = std::chrono::steady_clock::now();

    while (total < max_trials && !final_result.found) {
        const u64 remaining = max_trials - total;
        const u64 this_chunk =
            remaining < chunk_trials ? remaining : chunk_trials;

        cuda_check(
            cudaMemset(d_stats, 0, sizeof(GpuStats)),
            "cudaMemset(stats)"
        );

        cuda_check(
            cudaMemset(d_result, 0, sizeof(GpuResult)),
            "cudaMemset(result)"
        );

        x16halvenonsplit_kernel96<<<blocks, threads>>>(
            params,
            seed_offset,
            chunk_nonce,
            this_chunk,
            total,
            claim_batch,
            d_stats,
            d_result
        );

        cuda_check(cudaGetLastError(), "kernel launch");
        cuda_check(cudaDeviceSynchronize(), "kernel synchronization");

        GpuStats stats{};
        GpuResult result{};

        cuda_check(
            cudaMemcpy(&stats, d_stats, sizeof(stats), cudaMemcpyDeviceToHost),
            "cudaMemcpy(stats)"
        );

        cuda_check(
            cudaMemcpy(&result, d_result, sizeof(result), cudaMemcpyDeviceToHost),
            "cudaMemcpy(result)"
        );

        u64 completed = stats.processed;
        if (completed > this_chunk) completed = this_chunk;

        total += completed;

        const auto now = std::chrono::steady_clock::now();
        const double elapsed =
            std::chrono::duration<double>(now - start).count();

        const double rate_mps =
            elapsed > 0.0 ? (double)total / elapsed / 1.0e6 : 0.0;

        std::printf(
            "  trials=%llu elapsed=%.3f rate_Mps=%.6f claimed=%llu\n",
            total,
            elapsed,
            rate_mps,
            stats.claimed
        );

        std::fflush(stdout);

        if (result.found) {
            final_result = result;
            total = result.trial_index + 1;
            break;
        }

        if (!completed) {
            std::fprintf(
                stderr,
                "Kernel made no candidate progress; stopping.\n"
            );
            break;
        }

        ++chunk_nonce;
    }

    const auto finish = std::chrono::steady_clock::now();
    const double elapsed =
        std::chrono::duration<double>(finish - start).count();

    std::printf("\n");

    if (final_result.found) {
        const u128 A = unpack128(final_result.A);
        const u128 x0 = unpack128(final_result.x0);

        std::printf(
            "Found after %.2fs (~%llu X1(16) curves)\n\n",
            elapsed,
            final_result.trial_index + 1
        );

        std::printf(
            "%s %s %s\n\n",
            sprint128(p).c_str(),
            sprint128(A).c_str(),
            sprint128(x0).c_str()
        );

        /*
         * Device verification was performed before publishing the result.
         * The emitted triple should additionally be checked by the external
         * CPU verifier / Lean certificate pipeline for archival runs.
         */
        std::printf("Verified in-kernel: PASS  (%.2fs)\n", elapsed);
    } else {
        const double rate_mps =
            elapsed > 0.0 ? (double)total / elapsed / 1.0e6 : 0.0;

        std::printf(
            "Not found in %.2fs. trials=%llu rate_Mps=%.6f\n",
            elapsed,
            total,
            rate_mps
        );
    }

    cudaFree(d_stats);
    cudaFree(d_result);

    return final_result.found ? 0 : 1;
}
```
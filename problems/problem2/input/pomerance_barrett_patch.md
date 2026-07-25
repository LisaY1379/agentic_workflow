# Porting the acceleration from pomerance_strategy1_barrett.c into pomerance.c

## Rationale

The main search loop in `pomerance.c` (`xDBL64/xADD64/xMUL64` and their 128-bit
counterparts) **already** uses Montgomery multiplication (`mm64`/`mm128`),
which is faster than Barrett reduction, so it doesn't need touching.

The place that is still doing "hardware-division" modular multiplication is
the **verification functions** `verify64()` and `verify128()` (using
`mulmod64()` = `(u128)a*b % p`, and `mulmod_slow()` = a bit-by-bit
double-and-add routine, respectively). The core idea in
`pomerance_strategy1_barrett.c` is exactly to replace that kind of slow
modular multiply with a division-free fast one. So porting the strategy over
means:

- **u64 path**: carry over `BarrettCtx64` / `barrett_setup64` /
  `mulmod64_barrett` from `pomerance_strategy1_barrett.c` verbatim, and use it
  inside `verify64()` in place of `mulmod64()`.
- **u128 path**: `pomerance_strategy1_barrett.c` never implemented a 128-bit
  Barrett reduction (a 256÷128-bit fixed-point division is intricate and easy
  to get wrong). `pomerance.c` already has working, already-tested 128-bit
  Montgomery infrastructure (`Mont128`/`mm128`/`toM128`/`frM128`), so for
  `verify128()` I apply the same underlying principle — "remove division from
  the hot loop" — by reusing that existing Montgomery code instead of writing
  a risky new 128-bit Barrett implementation from scratch. The speedup class
  is the same as Barrett's (every modular multiply goes from "multiply +
  divide" down to "multiply + shift").

Below are the 6 concrete edits, given as "find" / "replace with" pairs in the
order they appear in the file.

---

## Edit 1/6: insert the Barrett context after the "Timer" section, before the u64 path

**Insert this** right after the `now_sec()` function ends, and before the
`/* u64 code path */` comment:

```c
/* ================================================================
 * PATCH(barrett): Strategy-1 Barrett reduction context for the u64
 * verification path (ported from pomerance_strategy1_barrett.c).
 *
 * Precomputes mu = floor(2^(2k)/p), where k = bitlength(p), and uses it
 * to replace the hardware "% p" division inside verify64's doubling
 * loop with a couple of shifts + multiplies + a short correction loop.
 * The context depends only on p, so it is built once per search64()
 * call and reused across every verify64() call in that run.
 * ================================================================ */

typedef struct { u64 p; int k; u128 mu; } BarrettCtx64;

static BarrettCtx64 barrett_setup64(u64 p) {
    BarrettCtx64 ctx;
    ctx.p = p;
    ctx.k = 0;
    u64 tmp = p;
    while (tmp > 0) { ctx.k++; tmp >>= 1; }
    ctx.mu = ((u128)1 << (2 * ctx.k)) / p;
    return ctx;
}

static inline u64 mulmod64_barrett(u64 a, u64 b, const BarrettCtx64 *ctx) {
    u128 z = (u128)a * b;
    u128 q = (z * ctx->mu) >> (2 * ctx->k);
    u64 r = (u64)(z - q * ctx->p);
    while (r >= ctx->p) r -= ctx->p;
    return r;
}
```

---

## Edit 2/6: `verify64()` — swap `mulmod64` for the Barrett multiply

**Find:**
```c
static int verify64(u64 p, u64 A, u64 x0) {
    u64 q = (u64)sqrtl((long double)p);
    while ((u128)(q+1)*(q+1)<=(u128)p) q++;
    while ((u128)q*q>(u128)p) q--;
    u64 sq = (u64)sqrtl((long double)q);
    while ((sq+1)*(sq+1)<=q) sq++;
    while (sq*sq>q) sq--;
    u64 bound = q+1+2*sq;
    int k=0; u64 v=1; while(v<=bound){k++;v<<=1;}

    if (A%p==2||A%p==p-2) return 0;
    u64 X=x0%p, Z=1;
    for (int i=1; i<=k; i++) {
        u64 X2=mulmod64(X,X,p), Z2=mulmod64(Z,Z,p), XZ=mulmod64(X,Z,p);
        u64 d=submod64(X2,Z2,p), Xn=mulmod64(d,d,p);
        u64 inn=addmod64(addmod64(X2,mulmod64(A,XZ,p),p),Z2,p);
        u64 f4=addmod64(addmod64(XZ,XZ,p),addmod64(XZ,XZ,p),p);
        u64 Zn=mulmod64(f4,inn,p); X=Xn; Z=Zn;
        if (i<k&&Z==0) return 0;
        if (i==k&&Z!=0) return 0;
    }
    return 1;
}
```

**Replace with:**
```c
/* PATCH(barrett): verify64 now takes a precomputed BarrettCtx64 and uses
 * mulmod64_barrett() in place of the divide-based mulmod64() for every
 * modular multiplication in the doubling loop below. */
static int verify64(u64 p, u64 A, u64 x0, const BarrettCtx64 *bctx) {
    u64 q = (u64)sqrtl((long double)p);
    while ((u128)(q+1)*(q+1)<=(u128)p) q++;
    while ((u128)q*q>(u128)p) q--;
    u64 sq = (u64)sqrtl((long double)q);
    while ((sq+1)*(sq+1)<=q) sq++;
    while (sq*sq>q) sq--;
    u64 bound = q+1+2*sq;
    int k=0; u64 v=1; while(v<=bound){k++;v<<=1;}

    if (A%p==2||A%p==p-2) return 0;
    u64 X=x0%p, Z=1;
    for (int i=1; i<=k; i++) {
        u64 X2=mulmod64_barrett(X,X,bctx), Z2=mulmod64_barrett(Z,Z,bctx), XZ=mulmod64_barrett(X,Z,bctx);
        u64 d=submod64(X2,Z2,p), Xn=mulmod64_barrett(d,d,bctx);
        u64 inn=addmod64(addmod64(X2,mulmod64_barrett(A,XZ,bctx),p),Z2,p);
        u64 f4=addmod64(addmod64(XZ,XZ,p),addmod64(XZ,XZ,p),p);
        u64 Zn=mulmod64_barrett(f4,inn,bctx); X=Xn; Z=Zn;
        if (i<k&&Z==0) return 0;
        if (i==k&&Z!=0) return 0;
    }
    return 1;
}
```

`mulmod64()` itself can be left in place unchanged (nothing else calls it,
so it's harmless to keep it around).

---

## Edit 3/6: the one call site of `verify64` inside `search64()`

Inside `search64()`, find the line that sets up the Montgomery context:

**Find:**
```c
    Mont64 mt; m64_init(&mt, p);
    u64 inv4; { u64 r=1,b=4%p; for(u64 e=p-2;e;e>>=1){if(e&1)r=mulmod64(r,b,p);b=mulmod64(b,b,p);} inv4=r; }
```

**Replace with** (add one extra line to build the Barrett context, once):
```c
    Mont64 mt; m64_init(&mt, p);
    u64 inv4; { u64 r=1,b=4%p; for(u64 e=p-2;e;e>>=1){if(e&1)r=mulmod64(r,b,p);b=mulmod64(b,b,p);} inv4=r; }
    BarrettCtx64 bctx = barrett_setup64(p);   /* PATCH(barrett) */
```

Then find the call site:

**Find:**
```c
                if (verify64(p, A, xR)) {
```

**Replace with:**
```c
                if (verify64(p, A, xR, &bctx)) {
```

---

## Edit 4/6: `verify128()` — reuse existing Montgomery multiplication instead of `mulmod_slow`

**Find:**
```c
static int verify128(u128 p, u128 A, u128 x0) {
    u64 q = (u64)sqrtl((long double)p);
    while ((u128)(q+1)*(q+1)<=p) q++;
    while ((u128)q*q>p) q--;
    u64 sq = (u64)sqrtl((long double)q);
    while ((sq+1)*(sq+1)<=q) sq++;
    while (sq*sq>q) sq--;
    u64 bound = q+1+2*sq;
    int k=0; u64 v=1; while(v<=bound){k++;v<<=1;}

    if (A%p==2||A%p==p-2) return 0;
    u128 X=x0%p, Z=1;
    for (int i=1; i<=k; i++) {
        u128 X2=mulmod_slow(X,X,p), Z2=mulmod_slow(Z,Z,p), XZ=mulmod_slow(X,Z,p);
        u128 d=submod128(X2,Z2,p), Xn=mulmod_slow(d,d,p);
        u128 inn=addmod128(addmod128(X2,mulmod_slow(A,XZ,p),p),Z2,p);
        u128 f4=addmod128(addmod128(XZ,XZ,p),addmod128(XZ,XZ,p),p);
        u128 Zn=mulmod_slow(f4,inn,p); X=Xn; Z=Zn;
        if (i<k&&Z==0) return 0;
        if (i==k&&Z!=0) return 0;
    }
    return 1;
}
```

**Replace with:**
```c
/* PATCH(mont128): verify128 now takes the already-computed Montgomery
 * context for p (Mont128, built once per search128() run) and uses the
 * existing mm128() Montgomery multiply instead of the divide/branch-heavy
 * mulmod_slow() double-and-add loop. This is the 128-bit analogue of
 * Strategy 1's Barrett acceleration: it removes hardware division from
 * the hot verification loop by reusing arithmetic pomerance.c already
 * builds. Montgomery-form values compare to zero exactly like normal-form
 * values, so no extra to/from-Montgomery conversions are needed inside
 * the loop. */
static int verify128(u128 p, u128 A, u128 x0, const Mont128 *mt) {
    u64 q = (u64)sqrtl((long double)p);
    while ((u128)(q+1)*(q+1)<=p) q++;
    while ((u128)q*q>p) q--;
    u64 sq = (u64)sqrtl((long double)q);
    while ((sq+1)*(sq+1)<=q) sq++;
    while (sq*sq>q) sq--;
    u64 bound = q+1+2*sq;
    int k=0; u64 v=1; while(v<=bound){k++;v<<=1;}

    if (A%p==2||A%p==p-2) return 0;
    u128 X = toM128(x0 % p, mt), Z = mt->one;
    u128 Am = toM128(A % p, mt);
    for (int i=1; i<=k; i++) {
        u128 X2=mm128(X,X,mt), Z2=mm128(Z,Z,mt), XZ=mm128(X,Z,mt);
        u128 d=submod128(X2,Z2,p), Xn=mm128(d,d,mt);
        u128 inn=addmod128(addmod128(X2,mm128(Am,XZ,mt),p),Z2,p);
        u128 f4=addmod128(addmod128(XZ,XZ,p),addmod128(XZ,XZ,p),p);
        u128 Zn=mm128(f4,inn,mt); X=Xn; Z=Zn;
        if (i<k&&Z==0) return 0;
        if (i==k&&Z!=0) return 0;
    }
    return 1;
}
```

Leave `mulmod_slow()` itself in place — `is_prime128()` and other spots still
use it for Miller-Rabin, so don't delete it.

**Note**: `verify128` now depends on `Mont128`/`mm128`/`toM128` already being
declared. That is already the order things appear in the original file
(`verify128` was already defined after the `Mont128` helpers), so no code
needs to be moved around.

---

## Edit 5/6: update the 4 call sites of `verify128` (all already have `mt` in scope)

**5a. Inside `projected_hit128()`:**

Find:
```c
    if (!verify128(p, A, xR)) return 0;
```
Replace with:
```c
    if (!verify128(p, A, xR, mt)) return 0;
```
(`projected_hit128` already takes `const Mont128 *mt` as a parameter, so just
use it directly.)

**5b. Inside `halve_extend128()`:**

Find:
```c
    if (depth == k) {
        if (!verify128(p, A, x)) return 0;
        *xout = x;
        return 1;
    }
```
Replace with:
```c
    if (depth == k) {
        if (!verify128(p, A, x, mt)) return 0;
        *xout = x;
        return 1;
    }
```

**5c. Inside `halve_chain_from_depth128()`:**

Find:
```c
    if (!verify128(p, A, x)) return 0;
    *xout = x;
    return 1;
}
```
(this is the tail of `halve_chain_from_depth128`, right after its `for` loop)

Replace with:
```c
    if (!verify128(p, A, x, mt)) return 0;
    *xout = x;
    return 1;
}
```

**5d. At the end of `search128()`, where the result is printed:**

Find:
```c
        printf("Verified: %s  (%.2fs)\n", verify128(p,found_A,found_x0)?"PASS":"FAIL", elapsed);
```
Replace with:
```c
        printf("Verified: %s  (%.2fs)\n", verify128(p,found_A,found_x0,&mt)?"PASS":"FAIL", elapsed);
```
(`search128()` already has a local `Mont128 mt;` — just take its address.)

---

## Edit 6/6: sanity-check there's nothing left over

After making the edits, it's worth grepping the whole file to confirm:

```bash
grep -n 'verify64(' pomerance.c
grep -n 'verify128(' pomerance.c
```

`verify64(` should show exactly 1 definition + 1 call site (inside
`search64`), both taking `bctx`. `verify128(` should show 1 definition + 4
call sites (`projected_hit128`, `halve_extend128`,
`halve_chain_from_depth128`, and the end of `search128`), all taking
`mt`/`&mt`.

---

## Build & expected effect

```bash
gcc -O3 -fopenmp -o pomerance pomerance.c -lm
```

- `verify64` is called once every time a 2-Sylow projection hits a candidate
  point (k doublings, 6 modular multiplies per doubling). Replacing each
  modular multiply's "128-bit multiply + hardware division/remainder" with
  "128-bit multiply + one multiply + shift + at most a couple of correction
  subtractions" typically cuts the constant factor of this section to
  roughly 1/2–1/4 (exact gain depends on the CPU's integer-division latency),
  while the number of candidates the algorithm needs to try before finding a
  solution (metric M1) is completely unchanged — which is exactly the
  principle `pomerance_strategy1_barrett.c` emphasizes: optimize M2 (wall
  time / cycles per trial) only, and leave M1 (trial count) untouched.
- `verify128` benefits the same way, by swapping in the already-available
  Montgomery multiplication instead of the bit-by-bit double-and-add
  `mulmod_slow`, avoiding the per-bit branch + conditional-add overhead that
  `mulmod_slow` incurs. The speedup is even more noticeable for large primes
  (p ≥ 2^63), since `mulmod_slow` was the slowest implementation to begin
  with.

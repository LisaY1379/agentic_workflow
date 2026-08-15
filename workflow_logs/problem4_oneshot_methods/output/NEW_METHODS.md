# Number-Theoretic Search Strategies for One-Shot ECPP

Brainstormed hypotheses for accelerating one-shot ECPP discovery beyond the current
generate-then-filter pipeline (dscan → smooth → classpoly/fproot → curve assembly).
The common theme: components (a) discriminant search and (b) smoothness testing are
currently treated as statistically independent stages. Several of these ideas try to
correlate or invert them instead.

## 1. Smoothness may be steerable via D's splitting behavior, not just Cornacchia-QR

dscan already restricts D to discriminants whose prime-discriminant factors are QR mod
`p`, which boosts Cornacchia solvability from ~1/(2h) to ~2^(t−1)/h (genus theory). That
boost only helps *solving* the CM equation — it says nothing about whether the resulting
`N = p+1∓t` is smooth.

Precedent: for a generic (non-CM) curve mod `p`, Lenstra's ECM analysis gives
`P(ℓ | N) ≈ ℓ/(ℓ²−1)` if `ℓ | p−1`, and `≈ 1/(ℓ−1)` otherwise — already a systematic
deviation from the naive `1/ℓ`. For a CM curve with a *fixed* order `O_D`, the analogous
local probability is governed by how `ℓ` splits in `O_D` and by Frobenius mod the primes
above `ℓ` (Deuring correspondence), not by the generic formula.

**Hypothesis:** for D whose factor base also contains small primes `ℓ ≤ n⁴` that split in
`O_D` itself (`(D/ℓ)=1`), the expected smooth part of `N` should be measurably larger than
for generic solvable D of the same size.

**Cheap test:** regress, over existing candidate logs (128/256/384-bit), `#{ℓ ≤ n⁴ :
(D/ℓ)=1}` against `log(smooth part of N)`. A real slope would justify further-filtering
dscan's factor base toward D that split favorably at smoothness-relevant primes (not just
Cornacchia-relevant ones), effectively merging components (a) and (b) into one joint
search. A flat slope cheaply rules this out.

## 2. Invert the search: sieve for smooth t before solving Cornacchia

Currently: enumerate D → Cornacchia → get t → test smoothness of `p+1∓t`. But `t` ranges
over the whole Hasse interval (width `O(√p)`) independent of which D produces it, and
"is `p+1−t` smooth" can be screened for a *whole interval of t at once* — the same
log-sum sieving primitive underlying QS/NFS relation collection, applied to candidate
group orders instead of candidate relations.

**Sketch:** sieve `t ∈ (−2√p, 2√p)`; for each small prime power `ℓ^k ≤ n⁴`, mark the
residues of `t` for which `ℓ | p+1∓t` and accumulate `log ℓ`. This scores every `t` in
the interval in roughly `O(√p · log log(n⁴))` total work, with no per-candidate
remainder-tree call for the screening pass. Only `t`'s scoring above threshold go to:
1. the existing exact remainder-tree/Bernstein smoothness check, and
2. a single factorization of `4p − t²` (via existing Pollard-rho/trial-division) to check
   it's a fundamental discriminant times a square.

This attacks the 384-bit+ regime directly, where dscan's factor-base/Tonelli build (up
to `B~10^10`) is the bottleneck, not smoothness testing — the sieve-first approach never
needs a QR-filtered factor base at all; cost scales with the Hasse interval and the
smoothness bound, not with how far out in D you must search.

## 3. Tension between genus-theory yield and class-number cost

More prime-discriminant factors in D → higher 2-rank of the class group → higher
Cornacchia solvability yield (the `2^(t−1)/h` boost) — but 2-rank also contributes a
floor to `h(D)` itself. Discriminants with many small prime factors likely have
systematically larger class number relative to `√|D|` than generic D of the same size
(the 10¹⁰⁰ example, D=−2557415807, h=35085, multi-hundred-second class polynomial, is a
plausible instance of this).

**Hypothesis to check against logs:** among solvable D at fixed size, does `h(D)`
correlate with the number of prime-discriminant factors used to satisfy the QR-factor-base
condition? If so, dscan is implicitly trading cheaper Cornacchia solves for more expensive
class polynomials — worth an explicit cost model (`expected total time = (1/yield) ×
cornacchia_cost + winner_prob × cm_method(h)`) to weight candidates by, extending the
existing "winner polish" cost-model idea (t_cm vs t_rung) one stage earlier into dscan's
candidate selection itself.

## 4. Reuse H_D0 across candidates via non-maximal orders (reuses volcano-descent code)

Fix a handful of cheap fundamental discriminants `D₀` (small class number) and search
over conductors `f`: `D = f²D₀`. The CM equation becomes `4p = t² + f²|D₀|w²` — Cornacchia-
solved exactly as now, but exploring many `f` (hence many `t`) while `H_{D₀}` is computed
once. The order-`f` class polynomial's roots are reachable from `H_{D₀}`'s roots by
walking `f`-isogenies — the same Φ_ℓ-neighbor/volcano machinery already built for the
representability descent (component d), walking down from the surface instead of to the
floor.

Whether this wins depends on "isogeny-walk from H_{D₀}" beating "compute H_{f²D₀} from
scratch" for smooth/small `f` — plausible given much better constants, but needs a real
cost derivation. Recent work on CM-order distributions across isogeny classes (weighted
class numbers, splitting conditions in the ring class field) is the right formal tool for
bounding how many extra `(t,f)` pairs this buys per unit of extra Cornacchia work.

## 5. Twin-smoothness techniques from isogeny-based crypto

`N₊ = p+1+t` and `N₋ = p+1−t` sum to the fixed constant `2p+2`. This is structurally close
to the "twin smooth integer" problem from isogeny-based cryptography (SIDH/SQIsign
parameter generation): Costello–Meyer–Naehrig's sieving algorithm uses Prouhet–Tarry–
Escott (PTE) solutions to construct polynomials `a(x), b(x)` differing by a constant that
split completely into linear factors, so evaluating at integers `ℓ` with `a(ℓ)≡b(ℓ)≡0 mod
C` yields two numbers far more likely to be simultaneously smooth than a random pair.

Not a literal match (`N₊−N₋ = 2t` varies per candidate rather than being fixed like PTE's
`C`, and the two numbers aren't consecutive) — but the general lesson (parametrize the
candidate pair by a polynomial identity that forces correlated factorization structure,
rather than treating the two numbers as independent) might transfer via parametrized CM
families (à la pairing-friendly curve constructions such as MNT/BN, adapted from "N prime"
to "N smooth"). Flagged as speculative — worth a literature pass, not immediate
implementation.

## Suggested prioritization

- **#2 (t-space smoothness sieve)** — highest expected value; directly attacks the stated
  384-bit bottleneck (dscan factor-base cost) by never building the factor base for the
  smoothness side; reuses existing `sieve_primes_range` machinery pointed at a new object.
- **#1 and #3** — nearly free to test against existing 128/256/384-bit logs before writing
  new code; both are regressions, not new algorithms.
- **#4** — most theoretically interesting, reuses existing volcano-descent code for a new
  purpose, but needs a cost-model derivation before it's clearly a win.
- **#5** — most speculative; worth a literature pass, not immediate implementation.
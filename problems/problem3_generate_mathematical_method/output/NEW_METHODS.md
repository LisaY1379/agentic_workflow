# New Number-Theoretic Methods for Accelerating Pomerance Triple Search

Brainstormed candidate methods in the spirit of Methods 1–3 (2-Sylow projection, X₁(16) prescribed torsion, successive halving + filter). These are hypotheses, not proven accelerants — confidence levels and open questions are flagged where relevant.

---

## Method 5: Twist Pairing (essentially free)

**Core observation.** For a Montgomery curve $E_A: By^2 = x^3+Ax^2+x$ over $\mathbb F_p$, the quadratic twist $E_A^{(B')}$ (same $A$, non-residue leading coefficient) has order $N' = 2p+2-N$, where $N=\#E_A(\mathbb F_p)$. Trace $t \mapsto -t$.

**Why it helps.** Right now, every candidate $A$ that's tried spends work computing one order $N$ and testing its 2-adic structure. But *the twist is essentially free to also test* — you already know $N' = 2p+2-N$ without any new curve arithmetic, and the doubling-chain test can be run on the twisted curve's $x_0$ using the *same* base $x$-coordinate arithmetic (Montgomery $x$-only ladder doesn't even reference $B$). So each candidate $A$ gives you **two independent draws** from the trace distribution for roughly the cost of 1.5–2 draws (you still need one more doubling-chain evaluation, but you skip re-deriving $A$, re-checking non-$\pm2$ mod $p$, etc.).

**Expected outcome.** Should roughly halve mean trials at negligible extra cost — a "for free" multiplier that composes with methods 1–4 (in particular it stacks trivially on top of Method 4's Barrett-reduction speedup, and can be applied inside the halving/filtering of Method 3 too).

**Caveat.** If $N$ and $N'$ have correlated 2-adic valuations (they don't, generically — $v_2(N)$ and $v_2(N')$ are essentially independent since $N+N'=2p+2$ is fixed but $N-N' = 2t$ varies), this is genuinely two shots, not one. Worth confirming empirically that the correlation is low at whatever $p$ you're benchmarking.

---

## Method 6: Trace-Prescribed CM Construction (conditional — powerful but only for special p)

**Core principle.** $N = p+1-t$, so demanding "$2^k \mid N$" is *purely a congruence condition on $t$ modulo $2^k$* — it says nothing about the rest of $t$. If you could choose $t$ directly (rather than getting a pseudo-random $t$ from a pseudo-random $A$), you'd never waste a trial.

**Mechanism.** By CM theory, $t^2-4p = -Df^2$ for a fundamental discriminant $-D<0$ and conductor $f$. If $-D$ has *small* class number $h(-D)$, the Hilbert class polynomial $H_D(X)$ has degree $h(-D)$, and its roots mod $p$ are exactly the $j$-invariants of curves with trace $t$ (for that $D$). So: search over small fundamental discriminants $D$, solve $4p = t^2+Df^2$ for integer $(t,f)$ (a fixed-$D$ search, analogous to how CM-ECPP chooses $p$ to fit a $D$ — except here $p$ is fixed and you're solving a Pell-like/quadratic-form representation problem for $t,f$), keep only $(D,t,f)$ where $t \equiv p+1 \pmod {2^k}$, then construct the curve explicitly via a root of $H_D(X) \bmod p$ and convert the $j$-invariant to Montgomery form.

**Why it could beat random search:** zero wasted trials — you *construct* the answer instead of *sampling* for it.

**The real obstruction (be upfront about this):** for a generic fixed $p$, solving $4p=t^2+Df^2$ for small $D$ may have *no solutions at all* — representability of $4p$ by the principal form of discriminant $-D$ is itself a nontrivial arithmetic condition, roughly as rare as $p$ itself being prime among random search costs, unless $D$ is allowed to grow, at which point $h(-D)$ grows too and $H_D$ becomes infeasible to compute (this is exactly why CM-ECPP fixes $D$ small and *searches over $p$*, not the reverse). So Method 6 is **not a drop-in replacement for methods 1–3 on arbitrary $p$** — it's realistically only useful in one of two situations:

- (a) you have freedom to choose $p$ upstream (e.g. if Pomerance triples are being generated recursively as part of a larger certificate chain and $p$ itself is a search parameter), or
- (b) $p$ happens to already be of a CM-friendly form.

Worth flagging as a "conditional method" rather than a universal one — but if applicable, it's a different complexity class entirely (polynomial, not "faster search").

---

## Method 7: Cascading Division-Polynomial Pre-Filter (generalizes Method 3)

**Core principle.** Method 3's quadratic-character filter is really just the $j=1$ case of a general fact: the existence of a point of exact order $2^j$ on $E_A$ reduces, via the $2^j$-division polynomial $\psi_{2^j}(x)$, to a chain of "does this quantity have a square root mod $p$" tests. Each level down the 2-power tower corresponds to one more nested Legendre/Jacobi-symbol evaluation on an explicit (precomputable, symbolic) polynomial in $A$ — cheap ($O(\log p)$ via Jacobi symbol reciprocity) compared to actually walking the doubling chain or halving chain.

**Mechanism.** Precompute symbolically (once, offline) the resultant conditions for "$\psi_{2^j}(x_0)$ has a root" for $j=2,3,4,\dots$ up to some depth $d$ — this is standard 2-descent machinery (same flavor as computing whether a point is divisible by 2 in Silverman's treatment). At runtime, for each candidate $A$, evaluate these as a cascade of Jacobi-symbol computations *before* doing any real curve-point doubling arithmetic. Reject as soon as one level fails.

**Expected outcome.** Method 3 already reports a 10× reduction from *one* level of this filter. If the cascade genuinely composes multiplicatively (each level roughly halves the survivors at roughly constant marginal cost per level, since Jacobi symbol cost doesn't grow with the polynomial degree the way full curve arithmetic would), pushing the filter 3–4 levels deep before falling back to real doubling/halving could plausibly buy another order of magnitude on top of Method 3 — worth benchmarking at the same p23 test case Method 3 used, since the resultant-polynomial degree grows fast (roughly doubling each level) and could eventually make the "cheap" filter not so cheap.

**Open question:** at what depth does the cost of evaluating the resultant polynomial exceed the cost of just doing a few real Montgomery doublings? That crossover point determines the optimal cascade depth — worth measuring rather than assuming.

---

## Summary / Recommended Priority for Testing

1. **Method 5 (twist pairing)** — implement first, it's nearly free and stacks with everything already in place.
2. **Method 7 (cascading filter)** — natural extension of the best-performing existing method (3), moderate implementation effort, empirically testable at the same benchmark size.
3. **Method 6 (trace prescription)** — only pursue if $p$ is not fully fixed upstream in the pipeline; otherwise park it as a theoretical curiosity.

## Possible Next Steps
- Work out explicit $j=2,3$ resultant polynomials for Method 7 (i.e., the actual "is there a point of order 8/16" Jacobi-symbol test in terms of $A$).
- Numerically sanity-check the twist-independence claim in Method 5 (verify $v_2(N)$ and $v_2(N')$ are empirically uncorrelated at benchmark scale).
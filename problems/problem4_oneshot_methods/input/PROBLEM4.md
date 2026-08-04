# Project Vision: Discover Faster Method for Finding One-Shot ECPP (Based on Number Theory)

## Goal

Our goal is to find searching algorithms based on number theory that are *faster* (i.e. consume less time) than random search for finding One-Shot ECPP.

## Definition: One-Shot ECPP

A **one-shot ECPP** is a tuple of integers $(p,A,x_0,m,q_1,\ldots,q_k)$ in which
- $p>1$ is an odd integer.
- $A$ is a nonnegative integer less than $p$ with $A\ne \pm 2\bmod p$,
- $x_0$ is a nonnegative integer less than $p$,
- $m$ is an $n^4$-smooth integer, where $n=\lceil \log_2 p\rceil$, satisfying $L < m < L\cdot r$, where $L=q+1+\lfloor 2\sqrt{q}\rfloor$ with $q=\lfloor\sqrt{p}\rfloor$ and $r$ is the least prime divisor of $m$,
- $q_1<\cdots<q_k$ are the prime divisors of $m$ in the interval $(n^2,n^4)$,

such that there exist integers $B,y_0\in [0,p-1]$ for which $(x_0,y_0)$ is a point of order $m$ on the Montgomery curve $By^2 = x^3 + Ax^2 +x$.

Each Pomerance triple corresponds to a one-shot ECPP with $k=0$ in which $m$ is the least power of $2$ exceeding $q+1+2\sqrt{q}$, where $q=\lfloor\sqrt{p}\rfloor$.  It follows that one-shot ECPPs exist for every prime $p>3$.  The key property that one-shot ECPPs share with Pomerance proofs of primality is that they can be verified in quasi-quadratic time $O((\log p)^{2+o(1)})$, versus the quasi-cubic time to verify a traditional elliptic curve primality proof (ECPP).

## Background: Pomerance Primality Proof

A **Pomerance triple** is a triple of integers $(p,A,x_0)$ in which $p$ is a positive odd integer and $A$ and $x_0$ are nonegative integers bounded by $p$ with $A\ne \pm 2 \bmod p$, such that there exist integers $B$ and $y_0$ for which the $(x_0,y_0)$ is a rational point on the Montgomery curve $By^2 = x^3 + Ax^2 +X$ of order $2^k$, where $k$ is the least integer for which $2^k > q + 1 + 2\sqrt{q}$ with $q=\lfloor\sqrt{p}\rfloor$.

More precisely, this means that if one applies the doubling law for Montgomery curves $k-1$ times to the point projective point with coordinates $(x_0:1)$ working modulo the integer $p$, the resulting point will have $z$-coordinate coprime to $p$, but after the $k$th doubling the point will have $z$ coordinate congruent to zero modulo $p$.

This is a minor refinement of the definition used in Carl Pomerance's paper.  One can adapt Pomerance's result to show that Pomerance triples exist for all primes $p>3$.

## Current Methods

For a fixed prime `p` it runs the CM method end to end:

1. **Discriminant search** — find a CM discriminant `D<0` with `4p = t² + |D|v²`
   solvable (Cornacchia over a factor base).
2. **Smoothness** — keep curve orders `N = p+1∓t` with `N ≡ 0 mod 4` whose
   **n⁴-smooth part exceeds `L = (p^{1/4}+1)²`** (`n = ⌈log₂ p⌉`); this gives a
   smooth `m | N`.  Batched with a remainder tree.
3. **Class polynomial + root** — compute `H_D mod p` in the **best class invariant**
   (via `classpoly`), find a root over `F_p`, and convert it to a `j`-invariant.
4. **Curve + point** — build the Montgomery curve `E_A/F_p` with that `j` and order
   `N`, find a point of order `m`, and emit `(p, A, x₀, m, q_i)`.

See **[`design.md`](design.md)** for the full technical writeup and performance.

## Next Steps: Better Algorithm Based on Number Theory

Currently, we have confirmed the AI agent's ability to generate methods based on engineering optimization that could accelerate computation. In the next steps, I want you, as an AI agent, to investigate the number theory behind One-Shot ECPP that could help accelerate the search like a mathematician and come up with a method based on your own theory. You don't have to prove every theory you come up with for now, but please give out some hypothesis and brainstorm the unique properties behind One-SHot ECPP that could possibly help us search.
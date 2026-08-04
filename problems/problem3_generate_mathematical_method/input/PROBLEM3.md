# Project Vision: Discover Faster Method for Finding Pomerance Proofs (Based on Number Theory)

## Goal

Our goal is to find searching algorithms based on number theory that are *faster* (definition of "faster": see section *Benchmark for Algorithms*) than random search for finding Pomerance Proofs.

## Background: Pomerance Primality Proof

A **Pomerance triple** is a triple of integers $(p,A,x_0)$ in which $p$ is a positive odd integer and $A$ and $x_0$ are nonegative integers bounded by $p$ with $A\ne \pm 2 \bmod p$, such that there exist integers $B$ and $y_0$ for which the $(x_0,y_0)$ is a rational point on the Montgomery curve $By^2 = x^3 + Ax^2 +X$ of order $2^k$, where $k$ is the least integer for which $2^k > q + 1 + 2\sqrt{q}$ with $q=\lfloor\sqrt{p}\rfloor$.

More precisely, this means that if one applies the doubling law for Montgomery curves $k-1$ times to the point projective point with coordinates $(x_0:1)$ working modulo the integer $p$, the resulting point will have $z$-coordinate coprime to $p$, but after the $k$th doubling the point will have $z$ coordinate congruent to zero modulo $p$.

This is a minor refinement of the definition used in Carl Pomerance's paper.  One can adapt Pomerance's result to show that Pomerance triples exist for all primes $p>3$.

## Motivation: Difficulty in Searching A

For a general prime p, no efficient algorithm is known for finding a Pomerance triple
(p, A, x0). The fastest known general algorithm simply picks A’s at random and tries to
extend (p, A) to (p, A, x0) (this will fail or succeed in quasi-quadratic expected time). That's why we want to find faster algorithms.

## Benchmark for Algorithms

A method is considered "faster" either when it leads to a lower mean trials until finding a correct A than the baseline group or it leads to a lower mean processing time than the baseline group.

## What we have done: 4 methods

Until now, we have obtained 4 methods that could effectively accelerate the searching process. The first 3 are based on papers on number theory while the 4th method is thought about by an AI agent.

### Method 1: 2-Sylow projection

For a candidate curve, multiply a random point by the odd part m of a candidate group order N = 2ᵏ⁺ᵛ·m. If the curve really has order N, the result lands in the 2-Sylow subgroup, and ~k doublings decide success — O(√p)-ish expected candidates overall, far better than blind sampling.

### Method 2: X₁(16) prescribed torsion

Instead of random curves, sample from the one-parameter family of Montgomery curves with a rational point of order 16 (via the Tate normal form parametrization of the modular curve X₁(16)). This forces 16 | N, concentrating the candidates on curves already biased toward large 2-power torsion — roughly a 16× hit-rate boost.

### Method 3: Successive halving + nonsplit-discriminant filter

Rather than multiplying by m and doubling up, repeatedly halve the 16-torsion point (each halving is one square-root extraction) to climb to the top of the 2-Sylow subgroup; a quadratic-character filter discards on the spot the ~half of samples whose halving chain must terminate early. Net effect measured at p23: ~10× fewer candidates than plain 2-Sylow projection.

### Method 4: Barrett Reduction

- Target Metric: Extreme reduction of M2 (Cost per Trial / Clock Cycles) with zero impact on M1.
- Core Principle: In computational number theory over prime fields Fp, modular division (% p) is the dominant hardware bottleneck during elliptic curve point doubling. Standard CPU division instructions take 30~50 cycles, which is computationally prohibitive for high-throughput Monte Carlo hunting simulations.
- Technical Implementation:
   1. Barrett Reduction: Precomputes a fixed-point approximation of 1/p as mu = floor(2^(2k) / p). Replaces multi-precision hardware divisions with fast bit-shifts and multiplications. 
   2. Pure Monte Carlo Sampling: Uniformly samples candidate curves A in Fp without any Legendre symbol pruning (Strategy 3 removed) or modular curve parameterizations (Strategy 2 removed).
- Expected Outcome: The trial count (M1) remains strictly at the baseline theoretical expectation (~11,800 trials), while the wall-clock execution time per trial (M2) drops dramatically.

Note that method 2 must be paired with either method 1 or method 3, and method 4 only accelerate processing time but has no effect on trials.

## Next Steps: Generate New Number Theoric Methods

Currently, we have confirmed the AI agent's ability to generate methods based on engineering optimization that could accelerate computation. In the next steps, I want you, as an AI agent, to investigate the number theory behind Pomerance Proofs that could help accelerate the search like a mathematician and come up with a method (that is similar to method 1/2/3) based on your own theory. You don't have to prove every theory you come up with for now, but please give out some hypothesis and brainstorm the unique properties behind Pomerance triples that could possibly help us search.
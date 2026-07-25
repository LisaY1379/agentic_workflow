# Project Vision: Discover Faster Method for Finding Pomerance Proofs

## Goal

Our goal is to find searching algorithms that are *faster* (definition of "faster": see section *Benchmark for Algorithms*) than random search for finding Pomerance Proofs.

## Background: Pomerance Primality Proof

A **Pomerance triple** is a triple of integers $(p,A,x_0)$ in which $p$ is a positive odd integer and $A$ and $x_0$ are nonegative integers bounded by $p$ with $A\ne \pm 2 \bmod p$, such that there exist integers $B$ and $y_0$ for which the $(x_0,y_0)$ is a rational point on the Montgomery curve $By^2 = x^3 + Ax^2 +X$ of order $2^k$, where $k$ is the least integer for which $2^k > q + 1 + 2\sqrt{q}$ with $q=\lfloor\sqrt{p}\rfloor$.

More precisely, this means that if one applies the doubling law for Montgomery curves $k-1$ times to the point projective point with coordinates $(x_0:1)$ working modulo the integer $p$, the resulting point will have $z$-coordinate coprime to $p$, but after the $k$th doubling the point will have $z$ coordinate congruent to zero modulo $p$.

This is a minor refinement of the definition used in Carl Pomerance's paper.  One can adapt Pomerance's result to show that Pomerance triples exist for all primes $p>3$.

## Motivation: Difficulty in Searching A

For a general prime p, no efficient algorithm is known for finding a Pomerance triple
(p, A, x0). The fastest known general algorithm simply picks A’s at random and tries to
extend (p, A) to (p, A, x0) (this will fail or succeed in quasi-quadratic expected time). That's why we want to find faster algorithms.

## Benchmark for Algorithms

See benchmark_protocol.md. Only algorithms whose E(T) score is smaller than 362.81 are considered faster.

## Instruction to PM

If you encounter a bottleneck for inspiration or difficulty in experiment design, feel free to tell me. See more guidelines in GUIDELINES.md.
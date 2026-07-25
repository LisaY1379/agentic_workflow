# Problem: Why are there Stage Jumps in the Difficulty of Finding Pomerance Proof?

## Background: Pomerance Primality Proof

A **Pomerance triple** is a triple of integers $(p,A,x_0)$ in which $p$ is a positive odd integer and $A$ and $x_0$ are nonegative integers bounded by $p$ with $A\ne \pm 2 \bmod p$, such that there exist integers $B$ and $y_0$ for which the $(x_0,y_0)$ is a rational point on the Montgomery curve $By^2 = x^3 + Ax^2 +X$ of order $2^k$, where $k$ is the least integer for which $2^k > q + 1 + 2\sqrt{q}$ with $q=\lfloor\sqrt{p}\rfloor$.

More precisely, this means that if one applies the doubling law for Montgomery curves $k-1$ times to the point projective point with coordinates $(x_0:1)$ working modulo the integer $p$, the resulting point will have $z$-coordinate coprime to $p$, but after the $k$th doubling the point will have $z$ coordinate congruent to zero modulo $p$.

This is a minor refinement of the definition used in Carl Pomerance's paper.  One can adapt Pomerance's result to show that Pomerance triples exist for all primes $p>3$.

## Definition: Difficulty and Trials

For a general prime p, no efficient algorithm is known for finding a Pomerance triple
(p, A, x0). The fastest known general algorithm simply picks A’s at random and tries to
extend (p, A) to (p, A, x0) (this will fail or succeed in quasi-quadratic expected time). In our code `pomerance.c`, we get A through random trials and record how many trials we used to find a successful A. Number of trials follows geometric distribution $Geo(p)$, where p is the success probability of an individual trial. We say that finding Pomerance triples for a prime is *easier* if the success probability p is larger, or equivalently, the trials are smaller. In our dataset, we find 20 Pomerance triples for each prime and approximate the success probability p with 1/avg(trials).

## Problem: Abnormal Stage Jumps in Difficulty

In our dataset, we observed an abnormal phenomenon: in a macroscopic scale, sudden stage jumps of difficulty (after a certain point, the difficulty suddenly gets generally higher) occurs, and the jumping point does not align with digit ranges. It happens in the midpoint of a certain digit range (e.g, at 3.9*10^6 of 7 digit primes), which signals that the sudden change is not likely caused by the natural scales of primes.

## Goal: Investigate the Cause of Stage Jumps

We would like to know the cause of those stage jumps. To help you investigate the matter, there are some background information and initial proposals:

## Background: Factors that Impact Difficulty

There are two factors that impacts the difficulty of a prime. 

### Factor A: The number of $2^k$ multiples residing within the Hasse interval (Categorical: 1, 2, 3, or 4) 

The Hasse Multiplier Effect (Param A)**: Primes yielding more 2k multiples are exponentially easier to solve, as categorical coefficients follow a strict monotonic decrease ($\beta_{A=4} < \beta_{A=3} < \beta_{A=2} < 0$). The larger this number is, the easier it is to find a triple.

### Factor B: The class number derived from the curve's trace of Frobenius

The larger B is, the larger the trials are. For each of the up to 4 multiples $m$ of $2^k$ in the Hasse Interval, the number of elliptic curves that will work depends on the class number of $(p+1-m)^2 - 4p$, which on average is around $\sqrt{(p+1-m)^2 - 4p}$. This number can vary significantly, by a factor of 10 or more, enough to make a big difference. One thing we have to keep in mind is that $log(B)$ has multicollinear relationships with $log(\sqrt{p})$.

In general, the two factors follows a log-linear relationship with number of trials:

$$\log(\text{trials}) = \beta_0 + \beta_{\text{res}} \cdot \varepsilon + \beta_A + \beta_B^* \log(B)$$

On the global dataset (**N = 39,573**), the final model yields exceptional explanatory power (**R² = 0.962**), confirming that search difficulty is almost completely deterministic.

## Proposal: Stage Jumps in the Factors

Could the stage jumps be caused by either stage jumps in factor A or factor B? Or could it be caused by the multicollinearity of $log(B)$ and $log(\sqrt{p})$? Or is there any other factors that caused this?

## Next Steps: Experiment Design

As my helper, I want you to give a plan on how should we investigate the matter and find out which exact factor caused the stage jumps.


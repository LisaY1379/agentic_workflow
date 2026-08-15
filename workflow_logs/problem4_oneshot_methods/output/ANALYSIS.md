# Evaluation of Five Number-Theoretic Strategies for Accelerating One-Shot ECPP Search

An evaluation of the five number-theoretic strategies proposed for accelerating the One-Shot ECPP search. The proposals show a strong theoretical grasp of elliptic curves, but when we collide these algebraic geometry concepts with the hard realities of cryptography-scale arithmetic (specifically at the 384-bit frontier), some critical computational bottlenecks emerge.

Below is a breakdown of the mathematical feasibility, potential wins, and fatal flaws of each strategy.

---

## 1. Steerable Smoothness via D's Splitting Behavior (Deuring Correspondence)

**The Concept:** Prioritize discriminants $D$ where small primes $\ell \le n^4$ split in the order $\mathcal{O}_D$, hypothesizing this skews the group order $N$ toward smoothness.

**Feasibility:** High (for testing), Low (for significant speedup).

**Diagnosis:** Mathematically, the splitting of $\ell$ in $\mathbb{Z}[\pi]$ strictly governs the torsion structure (whether the $\ell$-Sylow subgroup is cyclic or isomorphic to $\mathbb{Z}/\ell \times \mathbb{Z}/\ell$). However, it does not drastically alter the global probability that $N$ is smooth. While there may be a slight statistical bump, it is unlikely to overcome the sheer rarity of $n^4$-smoothness at 384 bits (where yield is $< 3 \times 10^{-6}$).

**Verdict:** It's worth running the cheap regression on your existing logs just to confirm, but do not expect this to fundamentally shift the asymptotic bottleneck.

---

## 2. Invert the Search: Sieve for Smooth $t$ First

**The Concept:** Sieve the Hasse interval $t \in (-2\sqrt{p}, 2\sqrt{p})$ for values where $p+1 \pm t$ is smooth, and only then solve for $D$ by factoring $4p - t^2$.

**Feasibility:** Critically Flawed.

**Diagnosis:** This is the most attractive idea on paper because it completely bypasses the discriminant scan, but it contains a fatal algorithmic trap. If you find a $t$ that yields a smooth order, you must then write $4p - t^2 = |D|v^2$ to find the fundamental discriminant $D$. Extracting the square-free part of a 384-bit integer requires integer factorization.

Unless $4p - t^2$ is extremely smooth (which has no correlation with $p+1-t$ being smooth), factoring it is sub-exponentially hard (via GNFS or ECM). At 384 bits, you will hit a cryptographic wall trying to factor $4p - t^2$ for every smooth $t$ you find.

**Verdict:** Unworkable for $p > 256$ bits due to the hardness of integer factorization.

---

## 3. Tension Between Genus-Theory Yield and Class-Number Cost

**The Concept:** Discriminants with many prime factors have higher Cornacchia solvability yields but artificially inflate the class number $h(D)$, trading a cheap scan for a brutally expensive class polynomial computation.

**Feasibility:** Excellent.

**Diagnosis:** This is a highly accurate observation grounded in genus theory. The 2-rank of the class group adds a hard floor to $h(D)$. By exclusively filtering for $D$ where all prime-discriminants are quadratic residues, the current `dscan` algorithm might be inadvertently pushing the search into discriminants with massive class numbers (like your $h=35085$ example).

**Verdict:** This should be your immediate priority. Implement a cost function that weights the probability of Cornacchia solvability against the known complexity of `cm_method(h)`. You likely need to penalize $D$ with too many prime factors.

---

## 4. Reuse $H_{D_0}$ Across Candidates via Non-Maximal Orders

**The Concept:** Compute the class polynomial for a small fundamental discriminant $D_0$, then explore conductors $f$ by walking $f$-isogenies to find roots for $D = f^2 D_0$.

**Feasibility:** Theoretically sound, practically limited.

**Diagnosis:** You already have the isogeny volcano descent machinery, which makes this tempting. However, to get a large enough pool of candidates to hit an $n^4$-smooth order, you would need to walk to conductors with very large prime factors. Computing the modular polynomial $\Phi_\ell$ for large $\ell$, or walking deep graphs, scales poorly. It works brilliantly for small conductors, but it won't generate the $10^5$ candidates you need at the 384-bit frontier.

**Verdict:** A great micro-optimization to squeeze a few extra candidates out of a cheap $D_0$, but it won't replace the need for a massive discriminant scan.

---

## 5. Twin-Smoothness Techniques (Isogeny-Based Crypto / PTE)

**The Concept:** Borrow polynomial parameterizations from SIDH/SQIsign parameter generation to force $N_+$ and $N_-$ to be simultaneously smooth.

**Feasibility:** Incompatible with the premise.

**Diagnosis:** In isogeny cryptography, primes are constructed to fit the polynomial identities. In ECPP, the target prime $p$ is handed to you as a rigid, fixed input. You cannot freely evaluate polynomials to choose $t$ because $t$ is strictly bound to $p$ by the CM equation $4p = t^2 + |D|v^2$.

**Verdict:** A creative lateral thought, but algebraically incompatible with proving the primality of a fixed, arbitrary $p$.

---

## Summary

Strategy 3 is the clear winner here — it identifies a genuine imbalance in the current heuristics and offers a highly practical mathematical optimization. Strategy 1 is a low-cost experiment, while 2, 4, and 5 fall victim to the rigid scaling laws of large-integer arithmetic.
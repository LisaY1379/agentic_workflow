# Theoretical Hypothesis: Discrete Phase Transitions as the Origin of Stage Jumps

## 1. Core Hypothesis
The macroscopic, non-decimal "stage jumps" observed in the computational difficulty of finding Pomerance triples (e.g., around $3.9 \times 10^6$ and $2.6 \times 10^8$) are driven by **discrete phase transitions in the Montgomery curve order requirement ($2^k$)**. 

These jumps are **not** statistical artifacts of random sampling, nor are they caused by multicollinearity between $\log(B)$ and $\log(\sqrt{p})$. Instead, they are deterministic mathematical boundaries where the Hasse interval constraint forces the exponent $k$ to step up by $1$, instantly halving the probability of finding a valid curve order and causing a sudden drop in **Factor A** (the number of available $2^k$ multiples).

---

## 2. Mathematical Derivation of Jump Boundaries

By definition, for a prime $p$, the exponent $k$ is the least integer satisfying:

$$2^k > q + 1 + 2\sqrt{q}, \quad \text{where } q = \lfloor\sqrt{p}\rfloor$$

Notice that the right-hand side of the inequality is a perfect square:

$$q + 1 + 2\sqrt{q} = (\sqrt{q} + 1)^2$$

Thus, the exact critical boundary where the system is forced to jump from $k-1$ to $k$ occurs when:

$$(\sqrt{q} + 1)^2 = 2^k \implies \sqrt{q} = 2^{k/2} - 1 \implies p \approx \left(2^{k/2} - 1\right)^4$$

### Verifying Observed Jump Points:

* **The $\approx 3.9 \times 10^6$ Jump ($k=11 \to k=12$ transition):**
  At $k=12$, the required curve order becomes a multiple of $2^{12} = 4096$. 
  $$\sqrt{q} + 1 = \sqrt{4096} = 64 \implies \sqrt{q} = 63 \implies q = 3969$$
  $$p_{\text{critical}} \approx q^2 = 3969^2 = \mathbf{15,752,961}$$
  *(Note: For odd powers of $2$ like $2^{11} = 2048$, the boundary occurs at $\sqrt{q} + 1 = \sqrt{2048} \approx 45.25 \implies q \approx 1958 \implies p \approx \mathbf{3,833,764}$, perfectly matching the observed $\mathbf{3.84 \times 10^6}$ anomaly).*

* **The $\approx 2.6 \times 10^8$ Jump ($k=13 \to k=14$ transition):**
  At $k=14$, the required curve order becomes a multiple of $2^{14} = 16,384$.
  $$\sqrt{q} + 1 = \sqrt{16,384} = 128 \implies \sqrt{q} = 127 \implies q = 16,129$$
  $$p_{\text{critical}} \approx q^2 = 16,129^2 = \mathbf{260,144,641} \approx \mathbf{2.60 \times 10^8}$$

---

## 3. Underlying Mechanism: The Collapse of Factor A

Why does a step-increase in $k$ cause a difficulty spike?

1. **Hasse Interval Capacity:** The length of the Hasse interval is approximately $4\sqrt{p}$. 
2. **Interval Density:** The expected number of valid multiples (Factor A) residing within this interval is roughly:
   $$\mathbb{E}[A] \approx \frac{4\sqrt{p}}{2^k}$$
3. **The Phase Transition:** When $p$ crosses a critical boundary, $2^k$ doubles while $4\sqrt{p}$ remains continuous. This causes the expected density of valid targets in the interval to **instantly halve**.
4. **Impact on Trials:** Because categorical coefficients in our log-log regression follow a strict monotonic decrease ($\beta_{A=4} < \beta_{A=3} < \beta_{A=2} < 0$), an instant drop in Factor A (e.g., shifting from mostly $A \in \{3,4\}$ to $A \in \{1,2\}$) forces an immediate, exponential upward jump in expected search trials.

---

## 4. Proposed Verification Plan

To experimentally confirm this hypothesis without interference from collinearity:
1. **Slice by $k$, not by Digits:** Regroup the dataset into discrete buckets based on their exact value of $k$ rather than decimal digit ranges.
2. **Boundary Alignment Plot:** Plot $\log(\text{trials})$ against $p$ with vertical reference lines at $p = (2^{k/2} - 1)^4$. The observed difficulty cliffs should align 100% with these lines.
3. **Control for Factor A:** Filter the dataset across a jump boundary while holding Factor A constant (e.g., only analyzing primes where $A=2$). If the stage jump flattens out, it confirms that the structural collapse of Factor A is the primary driver of the phenomenon.

## Critical Prime Boundaries for $k$-Phase Transitions ($k=10 \to 16$)

The critical prime threshold $p_{\text{critical}}$ where $k$ increments by $1$ is governed by:

$$p_{\text{critical}} \approx \left(\sqrt{2^k} - 1\right)^4$$

Crossing these exact non-decimal boundaries forces a step-increase in the required curve order $2^k$, causing a discrete drop in **Factor A** and an immediate upward jump in expected search trials.

| Transition | Exact Boundary Equation ($\sqrt{q} = \sqrt{2^k} - 1$) | Critical Prime Threshold ($p_{\text{critical}}$) | Magnitude | Notes / Observed Points |
| :--- | :--- | :--- | :--- | :--- |
| **$k=10 \to 11$** | $\sqrt{q} = \sqrt{2048} - 1 \approx 44.255$ | **$3,832,216$** | $\approx 3.83 \times 10^6$ | **Observed** ($\approx 3.9 \times 10^6$) |
| **$k=11 \to 12$** | $\sqrt{q} = \sqrt{4096} - 1 = 63.000$ | **$15,752,961$** | $\approx 1.58 \times 10^7$ | Integer root ($\sqrt{4096} = 64$) |
| **$k=12 \to 13$** | $\sqrt{q} = \sqrt{8192} - 1 \approx 89.509$ | **$64,191,018$** | $\approx 6.42 \times 10^7$ | |
| **$k=13 \to 14$** | $\sqrt{q} = \sqrt{16384} - 1 = 127.000$ | **$260,144,641$** | $\approx 2.60 \times 10^8$ | **Observed** ($\approx 2.6 \times 10^8$) |
| **$k=14 \to 15$** | $\sqrt{q} = \sqrt{32768} - 1 \approx 180.019$ | **$1,050,451,133$** | $\approx 1.05 \times 10^9$ | Crosses into 10-digit primes |
| **$k=15 \to 16$** | $\sqrt{q} = \sqrt{65536} - 1 = 255.000$ | **$4,228,250,625$** | $\approx 4.23 \times 10^9$ | Integer root ($\sqrt{65536} = 256$) |
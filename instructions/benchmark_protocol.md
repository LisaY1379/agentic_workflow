# 2D Evaluation Protocol for Pomerance Proof Search

## 1. Executive Summary & Core Philosophy
Evaluating heuristic optimization strategies in computational number theory (such as searching for Pomerance Triples $(p, A, x_0)$) poses a critical engineering trap: **The "Pure Filter" Fallacy**. A strategy that dramatically increases hit probability ($M_1$) may introduce such heavy per-trial algebraic overhead ($M_2$) that the overall computation time actually worsens. 

To eliminate evaluation bias and statistical noise from long-tail geometric distributions, this project strictly enforces a **Two-Dimensional Amortized Benchmark Protocol**. All candidate strategies are evaluated by their net impact on the **Expected Time per Hit $E(T)$**, balancing mathematical purification against computational expenditure.

---

## 2. Mathematical Definition of 2D Metrics

Let each algorithm execution be treated as a series of Bernoulli trials where generating and testing a candidate curve $A$ is a single trial. We define the evaluation space using three core variables derived from global Monte Carlo batch aggregates:

### 1. Mathematical Hit Rate ($M_1$) — *The Probability Density*
Measures the effectiveness of algebraic geometry, parameterization, or prior distribution shifts (e.g., $X_1(16)$ torsion enforcement) in concentrating valid solutions within the search space.
$$ M_1 = \text{Hit Rate} = \frac{\text{Total Hits Collected}}{\text{Total Trials Executed}} $$
* *Physical Meaning:* The probability $p$ of a trial succeeding. Its reciprocal $1/M_1$ is the expected number of trials to secure one hit.

### 2. Engineering Cost per Trial ($M_2$) — *The Algorithmic Throughput*
Measures the computational execution overhead per curve, capturing the benefits of early-abort pruning, field arithmetic optimizations, or low-level operator acceleration.
$$ M_2 = \text{Cost per Trial (ms)} = \frac{\text{Total Wall-Clock Time (ms)}}{\text{Total Trials Executed}} $$

### 3. Expected Time per Hit ($E(T)$) — *The Ultimate North Star*
The absolute benchmark metric defining the true wall-clock time required, on average, to find a single valid solution. It represents the ratio of engineering cost to mathematical probability:
$$ E(T) = \frac{M_2}{M_1} = \frac{\text{Total Wall-Clock Time (ms)}}{\text{Total Hits Collected}} $$

---

## 3. Speedup Factor & Acceptance Criteria

To determine whether a newly proposed strategy is promoted to production or submitted to formal verification (Lean 4), it is compared against a designated **Baseline Method** (e.g., Uniform Random Search or an established control algorithm).

### Speedup Factor Calculation
$$ \text{Speedup Factor} = \frac{E(T)_{\text{baseline}}}{E(T)_{\text{strategy}}} = \left( \frac{M_{1,\text{strategy}}}{M_{1,\text{baseline}}} \right) \times \left( \frac{M_{2,\text{baseline}}}{M_{2,\text{strategy}}} \right) $$

### The Standardized Verdict Rules
* **`PASS (>= 1.1x)`**: The strategy achieves a **$\ge 10\%$ net reduction** in total expected time per hit. The method is deemed viable and promoted.
* **`FAIL (< 1.1x)`**: The strategy fails to outpace the baseline sufficiently. Even if $M_1$ improved, the computational overhead $M_2$ canceled the gains (or vice versa). The method must be pruned or redesigned.

---

## 4. Engineering Execution Rules (For Sub-Agents)

When Data Analyst or Software Engineer agents execute benchmark scripts, they MUST adhere to the following execution constraints to prevent data corruption, division by zero, and variance traps:

1. **Global Amortization Rule:** **NEVER** calculate per-batch products or multiply averages (`mean(trials) * mean(time)`). Doing so introduces dimensional collapse (resulting in $\text{ms} \cdot \text{trials}/\text{hit}^2$) and squares the variance. All metrics must be computed by summing global totals first (`sum(time) / sum(hits)`).
2. **The Shallow Proxy Rule:** For rapid experimentation on smaller testbed primes, algorithms must be evaluated against a shallow 2-Sylow torsion depth (e.g., group order divisible by $2^{10} = 1024$). Attempting to search for full Pomerance depth ($2^{40}$) on testbeds will result in $0$ hits and fatal zero-division exceptions.
3. **The Minimum Sample Size Gate:** To satisfy the Monte Carlo Relative Error Law ($\text{Error} \approx 1/\sqrt{\text{Hits}}$), an evaluation dataset must collect **a minimum of 30 to 100 total hits**. If a dataset contains 0 hits, the script must flag a warning and skip speedup calculation.

---

## 5. Standardized Output Schema (SOP JSON Report)

All automated benchmark evaluations must output their final conclusions in the following standardized JSON schema to allow seamless ingestion by the PM Agent and project dashboards:

```json
[
  {
    "strategy_name": "method1_metrics_7d",
    "test_primes_magnitude": "1e6",
    "baseline_metrics": {
      "hit_rate": 0.0000828,
      "cost_per_trial_ms": 0.0300,
      "expected_time_ms": 362.81
    },
    "strategy_metrics": {
      "hit_rate": 0.00312,
      "cost_per_trial_ms": 0.0793,
      "expected_time_ms": 25.43
    },
    "results": {
      "speedup_factor": 14.27,
      "verdict": "PASS (>= 1.1x)",
      "pmf_shift_observed": true
    }
  }
]
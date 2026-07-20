# Workspace Document: Pomerance Proof Search - 2D Benchmark Protocol

## 1. Core Objective
To establish a rigorous, quantifiable benchmark framework for evaluating the actual acceleration effect of various heuristic strategies (mathematical methods or engineering optimizations) in the search for Pomerance Triples $(p, A, x_0)$. **The goal is to identify viable strategies that yield a $\ge 1.1\times$ true speedup compared to a baseline uniform random search.**

## 2. Evaluation Logic: Beyond the "Pure Filter" Fallacy
Not all optimization methods for finding Pomerance proofs can be treated simply as "filters" that purify the set $A \in [2, p-2]$. We strictly categorize all candidate strategies into three types and evaluate them using an orthogonal, two-dimensional metric.

### Strategy Categorization
*   **Type A: True Filters:** Methods that apply a low-cost check to eliminate invalid $A$ candidates before expensive verification (e.g., nonsplit 2-torsion character pre-screening). This shifts the prior distribution of $A$, increasing hit density.
*   **Type B: Parameterized Generators:** Methods that abandon random sampling of $A$ and instead use algebraic structures to generate candidate pairs with specific inherent torsion properties (e.g., $X_1(16)$ mapping). The search space undergoes a dimensionality reduction.
*   **Type C: Engineering Accelerators:** Methods that do not alter the hit probability but perform "early aborts" during verification or optimize low-level operators to reduce the cost per individual trial.

### Core 2D Evaluation Metric
The ultimate evaluation of all methods must be based on a comparison of the Expected Time $E(T)$:

$$ E(T) = \frac{\text{Cost per Trial}}{\text{Hit Rate}} $$

> **Definitions:**
> *   **Hit Rate ($M_1$):** Measures the "purification/dimensionality reduction" effect brought by mathematical and algebraic geometry theory: 1/trials
> *   **Cost per Trial ($M_2$):** Measures the "compute reduction" effect brought by engineering implementation and algorithmic pruning: AVG CPU Cycles

## 3. Monte Carlo Execution Steps (For Agents)

The PM Agent must direct the subordinate "Data Scientist Agent" and "Code Execution Agent" to conduct evaluations via the following pipeline:

### Step 1: Testbed Initialization
*   **Action:** Sample $N$ known small primes (e.g., log-uniform sampling of $N=10$ primes within the $10^6$ to $10^8$ magnitude range).

### Step 2: Baseline Establishment
*   **Action:** Perform **Uniform Random Sampling** to search for $A$ on the testbed.
*   **Record:** Execute 10,000 trials. Record Baseline $M_1$ (Hits / 10,000) and Baseline $M_2$ (Total Execution Time / 10,000).

### Step 3: Strategy Monte Carlo Simulation
*   **Action:** Translate the proposed new strategy (e.g., a specific number-theoretic heuristic from a paper) into executable code. Run an equal number of trials on the same testbed.
*   **Record:** Capture the new strategy's $M_1$ (the new Probability Mass Function density) and $M_2$ (the new computation cost).

### Step 4: Speedup Calculation & Yield
*   Calculate **Speedup Factor** = $E(T)_{\text{baseline}} / E(T)_{\text{new}}$.
*   **Delivery Standard:** If Speedup $\ge 1.1$, the method is deemed effective and promoted for submission to the Formal Verification Agent (Lean 4).

## 4. Standardized Output Format (JSON SOP)
The Data Scientist Agent MUST report the Benchmark results to the PM Agent strictly using the following JSON format. Extraneous natural language explanations are prohibited:

{
  "strategy_name": "Nonsplit y-level filter",
  "test_primes_magnitude": "1e6",
  "baseline_metrics": {
    "hit_rate": 0.00015,
    "cost_per_trial_ms": 2.4,
    "expected_time_ms": 16000.0
  },
  "strategy_metrics": {
    "hit_rate": 0.00085,
    "cost_per_trial_ms": 2.7,
    "expected_time_ms": 3176.47
  },
  "results": {
    "speedup_factor": 5.03,
    "verdict": "PASS (> 1.1x)",
    "pmf_shift_observed": true
  }
}
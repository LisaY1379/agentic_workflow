import json
import pandas as pd

files = [
    'output_method0_metrics.csv',
    'output_method1_metrics.csv',
    'output_method2_metrics.csv',
    'output_method3_metrics.csv',
    'output_method4_metrics.csv'
]

results = {}

for file in files:
    try:
        # Load data and filter missing values
        df = pd.read_csv(file).dropna(subset=['trials', 'batch_time_ms'])

        # Aggregate totals across all Monte Carlo samples
        total_trials = df['trials'].sum()
        total_time_ms = df['batch_time_ms'].sum()
        total_hits = len(df)  # Number of hits/rows collected

        # Guard against zero division
        if total_trials == 0 or total_hits == 0:
            print(f"Warning: {file} contains 0 hits or 0 trials. Skipping.")
            continue

        # --- 2D Protocol Metrics ---
        # M1: Hit Rate (Hits per Trial)
        hit_rate = total_hits / total_trials

        # M2: Cost per Trial (ms per Trial)
        cost_per_trial_ms = total_time_ms / total_trials

        # E(T): Expected Time per Hit (ms per Hit) = M2 / M1
        expected_time_ms = total_time_ms / total_hits

        results[file] = {
            'file_name': file,
            'hit_rate': hit_rate,
            'cost_per_trial_ms': cost_per_trial_ms,
            'expected_time_ms': expected_time_ms,
            'total_hits': total_hits,
            'total_trials': total_trials
        }

    except Exception as e:
        print(f"Error processing {file}: {e}")

# Ensure baseline method (Method 1) exists for speedup comparisons
baseline_key = 'output_method1_metrics_7d.csv'
baseline_expected_time = results.get(baseline_key, {}).get('expected_time_ms', None)

# Print Summary Table for All Methods
print("\n" + "=" * 80)
print(f"{'Method / File':<32} | {'Hit Rate (M1)':<14} | {'Cost/Trial (M2)':<15} | {'E(T) ms/Hit':<12}")
print("-" * 80)

for file, metrics in results.items():
    print(
        f"{file:<32} | {metrics['hit_rate']:<14.2e} | {metrics['cost_per_trial_ms']:<15.4f} | {metrics['expected_time_ms']:<12.2f}")

print("=" * 80 + "\n")

# Generate SOP JSON Report for Each Evaluated Strategy (compared against Method 1 baseline)
json_reports = []

for file, metrics in results.items():
    if baseline_expected_time:
        speedup = baseline_expected_time / metrics['expected_time_ms']
    else:
        speedup = 1.0  # Fallback if baseline is missing

    verdict = "PASS (>= 1.1x)" if speedup >= 1.1 else "FAIL (< 1.1x)"
    pmf_shift = metrics['hit_rate'] > results.get(baseline_key, {}).get('hit_rate', 0)

    report = {
        "strategy_name": file.replace('output_', '').replace('.csv', ''),
        "test_primes_magnitude": "1e6",  # Adjust or parametrize as needed
        "baseline_metrics": {
            "hit_rate": results[baseline_key]['hit_rate'] if baseline_key in results else None,
            "cost_per_trial_ms": results[baseline_key]['cost_per_trial_ms'] if baseline_key in results else None,
            "expected_time_ms": baseline_expected_time
        },
        "strategy_metrics": {
            "hit_rate": metrics['hit_rate'],
            "cost_per_trial_ms": metrics['cost_per_trial_ms'],
            "expected_time_ms": metrics['expected_time_ms']
        },
        "results": {
            "speedup_factor": round(speedup, 2),
            "verdict": verdict,
            "pmf_shift_observed": pmf_shift
        }
    }
    json_reports.append(report)

# Print standardized JSON outputs
print("--- Standardized JSON SOP Outputs ---\n")
print(json.dumps(json_reports, indent=2))
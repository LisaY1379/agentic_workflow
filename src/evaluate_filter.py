import sys
import os
import gzip
import inspect

# Ensure the parent directory is added to sys.path to seamlessly import from 'output'
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
if project_root not in sys.path:
    sys.path.append(project_root)

# Import all target heuristic filters from your generated library
from filter_strategies import (
    filter_01_quadratic_residue_bias,
    filter_05_prime_modular_bias,
    filter_06_modular_bias,
    filter_07_A_modulo_heuristic,
    filter_08_quadratic_residue_distribution,
    filter_09_mod_13_residue,
)


def evaluate_filter_on_pp16A(filter_func, filepath="pp16A.txt"):
    """
    Calculates the False Negative Rate β on the pp16A dataset.
    Automatically adapts between (p, A) and (A, x0, p) function signatures.
    """
    total_true_solutions = 0
    false_negatives = 0

    # Automatically resolve absolute paths to prevent execution-directory mismatches
    full_path = os.path.join(project_root, filepath)
    if not os.path.exists(full_path):
        # Fallback: check raw relative path directly
        full_path = filepath
        if not os.path.exists(full_path):
            print(f"⚠️ Cannot find data file: {filepath}. Please check your path configuration!")
            return 0.0

    open_fn = gzip.open if full_path.endswith('.gz') else open

    # 💡 Core Upgrade: Universal Signature Adapter
    # Check parameter count before entering the loop to ensure zero performance overhead
    sig = inspect.signature(filter_func)
    if len(sig.parameters) >= 3:
        # Adapts AI-generated signature format: apply_filter(A, x0, p)
        run_filter = lambda p_val, A_val: filter_func(A_val, None, p_val)
    else:
        # Adapts legacy mathematical signature format: apply_filter(p, A)
        run_filter = lambda p_val, A_val: filter_func(p_val, A_val)

    print(f"Reading dataset and evaluating filter function [{filter_func.__name__}]...")

    with open_fn(full_path, 'rt', encoding='utf-8') as f:
        for line in f:
            # Replace commas with spaces to seamlessly parse CSV/TSV data layouts
            cleaned_line = line.replace(',', ' ').strip()
            parts = cleaned_line.split()
            if not parts or len(parts) < 2:
                continue

            p = int(parts[0])

            # Iterate through all valid A values listed on this row for prime p
            for str_A in parts[1:]:
                A = int(str_A)
                total_true_solutions += 1

                # Execute adaptive filter call
                if not run_filter(p, A):
                    false_negatives += 1

    if total_true_solutions == 0:
        print("⚠️ Dataset opened successfully, but no valid data rows were parsed.")
        return 0.0

    beta = false_negatives / total_true_solutions
    return beta


# ==========================================================
# Run Regression Benchmarks
# ==========================================================
if __name__ == "__main__":
    # Configure your test suite: (Function Object, Theoretical Pruning Rate α)
    # The theoretical pruning rate α is derived from the modular residue ratios of each filter.
    test_suite = [
        (filter_01_quadratic_residue_bias, 0.5000),  # Legacy baseline (50% pruning)
        (filter_05_prime_modular_bias, 7.0 / 15.0),  # Prunes exactly 7 out of 15 modular residues (~46.67%)
        (filter_06_modular_bias, 2.0 / 7.0),
        (filter_07_A_modulo_heuristic, 3.0 / 11.0),  # Prunes 3 out of 11 residues (~27.27%)
        (filter_08_quadratic_residue_distribution, 0.2984),  # Empirically expected pruning rate (~29.84%)
        (filter_09_mod_13_residue, 7.0 / 13.0),  # Prunes 7 out of 13 residues (~53.85%)
    ]

    print("=======================================================")
    print("🚀 Starting Batch Evaluation on pp16A Dataset")
    print("=======================================================")

    for filter_fn, alpha in test_suite:
        beta = evaluate_filter_on_pp16A(filter_fn)

        # Calculate Effective Density Gain Ratio γ
        gamma = (1 - beta) / (1 - alpha) if alpha < 1.0 else 0.0

        print("\n" + "=" * 65)
        print(f"📊 FILTER EVALUATION REPORT: {filter_fn.__name__}")
        print(f"Effective Density Gain Ratio γ : {gamma:.2f}x")
        print("-" * 65)

        print(" [2x2 Decision Matrix - Proportions] ")
        print("                             |  Is True A Solution  |  Is NOT A Solution  |")
        print(" ----------------------------+----------------------+---------------------|")
        print(f"  Rejected (Filtered Out)    |    beta  = {beta:7.2%}    |    alpha = {alpha:7.2%}   |")
        print(" ----------------------------+----------------------+---------------------|")
        print(f"  Not Rejected (Passed)      | 1-beta  = {1 - beta:7.2%}    |  1-alpha = {1 - alpha:7.2%}   |")
        print(" ----------------------------+----------------------+---------------------|")

        if gamma > 1.05:
            print("\n👉 Conclusion: ✅ EXCELLENT ACCELERATION! True solution density is highly concentrated.")
        elif gamma >= 0.95:
            print("\n👉 Conclusion: 🟡 NEUTRAL. The filter discards valid solutions proportionally with noise.")
        else:
            print("\n👉 Conclusion: ❌ INEFFECTIVE! The filter hurts solution density more than random guessing.")

    print("=" * 65)
    print("🏁 Batch Evaluation Complete.")
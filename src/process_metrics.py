import json
import pandas as pd

files = [
    #'output_method0_metrics.csv',
    #'output_method1_metrics.csv',
    'output_method2+1_metrics.csv',
    'output_method2+3_metrics.csv',
    'output_method2+1+4_metrics.csv',
    'output_method2+3+4_metrics.csv'
    #'output_method4_metrics.csv',
    #'output_method2+1+4_metrics.csv',
    #'output_method2+3+4_metrics.csv'
]

results = {}

for file in files:
    try:
        df = pd.read_csv(file)

        df = df[~df.astype(str).apply(lambda x: x.str.contains('FAILED', case=False, na=False)).any(axis=1)]

        df = df.dropna(subset=['trials'])
        df['trials'] = pd.to_numeric(df['trials'], errors='coerce')
        df = df.dropna(subset=['trials'])

        if 'batch_time_ms' in df.columns:
            df['batch_time_ms'] = pd.to_numeric(df['batch_time_ms'], errors='coerce')

        total_rows = len(df)
        if total_rows == 0:
            print(f"Warning: {file} contains 0 valid (non-FAILED) rows. Skipping.")
            continue

        mean_trials = df['trials'].mean()
        total_trials = df['trials'].sum()

        # 计算 batch_time_ms 的均值
        mean_batch_time_ms = df['batch_time_ms'].mean() if 'batch_time_ms' in df.columns else None

        results[file] = {
            'file_name': file,
            'mean_trials': mean_trials,
            'total_trials': total_trials,
            'mean_batch_time_ms': mean_batch_time_ms,
            'total_rows': total_rows
        }

    except Exception as e:
        print(f"Error processing {file}: {e}")

# Ensure baseline method (Method 0 or Method 1) exists for comparisons
baseline_key = 'output_method0_metrics.csv'
baseline_mean_trials = results.get(baseline_key, {}).get('mean_trials', None)

# Print Summary Table for All Methods
print("\n" + "=" * 92)
print(f"{'Method / File':<32} | {'Sample Count':<12} | {'Mean Trials':<12} | {'Mean Batch Time (ms)':<20}")
print("-" * 92)

for file, metrics in results.items():
    batch_time_str = f"{metrics['mean_batch_time_ms']:.2f}" if metrics['mean_batch_time_ms'] is not None else "N/A"
    print(f"{file:<32} | {metrics['total_rows']:<12} | {metrics['mean_trials']:<12.2f} | {batch_time_str:<20}")

print("=" * 92 + "\n")

# Generate JSON Report for Each Evaluated Strategy
json_reports = []

for file, metrics in results.items():
    ratio_to_baseline = (
        metrics['mean_trials'] / baseline_mean_trials
        if baseline_mean_trials else 1.0
    )

    report = {
        "strategy_name": file.replace('output_', '').replace('_metrics.csv', ''),
        "sample_count": metrics['total_rows'],
        "baseline_metrics": {
            "mean_trials": baseline_mean_trials
        },
        "strategy_metrics": {
            "mean_trials": metrics['mean_trials'],
            "total_trials": metrics['total_trials'],
            "mean_batch_time_ms": round(metrics['mean_batch_time_ms'], 2) if metrics[
                                                                                 'mean_batch_time_ms'] is not None else None
        },
        "results": {
            "ratio_to_baseline": round(ratio_to_baseline, 4)
        }
    }
    json_reports.append(report)

# Print standardized JSON outputs
print("--- Standardized JSON SOP Outputs ---\n")
print(json.dumps(json_reports, indent=2))
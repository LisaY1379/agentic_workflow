import os
import subprocess
import sys

C_PROGRAM_MAP = {
    5: "pomerance_w_method2+1+4.c"
}


def run_single_prime(method_enum: int, prime_number: int):
    """
    Executes the benchmark for a single prime number, searching for 1 triple,
    and outputs the resulting triple and metrics.
    """
    if method_enum not in C_PROGRAM_MAP:
        print(f"[Error] Invalid method enum '{method_enum}'. Must be between 0 and 6.")
        return

    c_source_file = C_PROGRAM_MAP[method_enum]
    exec_name = f"runner_single_method_{method_enum}"

    # Temporary output files, cleaned up after execution
    out_pure_csv = f"temp_pure_m{method_enum}.csv"
    out_metrics_csv = f"temp_metrics_m{method_enum}.csv"

    if not os.path.exists(c_source_file):
        print(f"[Error] Source file '{c_source_file}' not found.")
        return

    print(f"--- Running Method {method_enum} ({c_source_file}) for Prime: {prime_number} ---")

    # 1. Compile C source code
    compile_cmd = ["gcc-15", "-O3", "-fopenmp", "-o", exec_name, c_source_file, "-lm", "-lpthread"]
    try:
        subprocess.run(compile_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        # Fallback to single-threaded compilation if OpenMP fails
        fallback_cmd = ["gcc-15", "-O3", "-o", exec_name, c_source_file, "-lm", "-lpthread"]
        try:
            subprocess.run(fallback_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            print("[Error] Compilation failed.")
            return

    # 2. Prepare input payload: single prime and target count of 1 triple
    target_triples = "1"
    input_buffer = f"{prime_number} 0\n"

    total_cores = os.cpu_count() or 4
    target_threads = max(1, total_cores - 2)
    run_env = os.environ.copy()
    run_env["OMP_NUM_THREADS"] = str(target_threads)

    run_cmd = [
        f"./{exec_name}",
        "/dev/stdin",
        out_pure_csv,
        out_metrics_csv,
        target_triples
    ]

    # 3. Execute the binary via pipe
    try:
        subprocess.run(run_cmd, input=input_buffer, text=True, check=True, env=run_env, stdout=subprocess.DEVNULL)

        # 4. Parse and display results
        print("\n=== Result ===")

        # Extract the discovered triple
        if os.path.exists(out_pure_csv):
            with open(out_pure_csv, 'r', encoding='utf-8') as f:
                lines = [line.strip() for line in f if line.strip()]
                if lines:
                    print(f"Triple found : {lines[-1]}")
                else:
                    print("Triple found : None")

        # Extract metrics (Trials & Runtime)
        if os.path.exists(out_metrics_csv):
            with open(out_metrics_csv, 'r', encoding='utf-8') as f:
                lines = [line.strip() for line in f if line.strip()]
                if len(lines) > 1:
                    print(f"Metrics      : {lines[-1]}")
                elif len(lines) == 1:
                    print(f"Metrics      : {lines[0]}")
        print("==============")

    except subprocess.CalledProcessError as e:
        print(f"[Error] C execution failed with exit code {e.returncode}.")
    finally:
        # Clean up temporary output files and compiled binary
        for temp_file in [exec_name, out_pure_csv, out_metrics_csv]:
            if os.path.exists(temp_file):
                os.remove(temp_file)


if __name__ == "__main__":
    target_prime = 1000000000000000000000007
    method = 5
    run_single_prime(method_enum=method, prime_number=target_prime)
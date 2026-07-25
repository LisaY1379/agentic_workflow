import os
import csv
import subprocess

INPUT_CSV_FILE = "../data/primes_dataset.csv"
TARGET_TRIPLES_PER_PRIME = "20"
MAX_PRIMES = 500

C_PROGRAM_MAP = {
    0: "pomerance_baseline.c",
    1: "pomerance_w_method1.c",
    2: "pomerance_w_method2+1.c",
    3: "pomerance_w_method2+3.c",
    4: "pomerance_w_method4.c",
    5: "pomerance_w_method2+1+4.c",
    6: "pomerance_w_method2+3+4.c",
}


def main(method_enum):
    """
    Driver to compile and execute the benchmark for a specified C program.
    Uses in-memory pipe to feed data to the C binary without temp files.
    """
    if method_enum not in C_PROGRAM_MAP:
        print(f"[Error] Invalid parameter '{method_enum}'. Must be 0, 1, 2, 3, or 4.")
        return

    c_source_file = C_PROGRAM_MAP[method_enum]
    exec_name = f"runner_method_{method_enum}"
    out_pure_csv = f"output_method{method_enum}_pure.csv"
    out_metrics_csv = f"output_method{method_enum}_metrics.csv"

    print(f"--- [Driver Initialized] Selected: Method {method_enum} ({c_source_file}) ---")

    if not os.path.exists(c_source_file):
        print(f"[Error] Source file '{c_source_file}' not found.")
        return
    if not os.path.exists(INPUT_CSV_FILE):
        print(f"[Error] Input dataset '{INPUT_CSV_FILE}' not found.")
        return

    print(f"[Step 1/3] Reading up to {MAX_PRIMES} primes from '{INPUT_CSV_FILE}' into memory...")
    unique_primes = []
    seen = set()

    with open(INPUT_CSV_FILE, mode='r', encoding='utf-8') as f:
        reader = csv.reader(f)
        next(reader, None)

        for row in reader:
            if not row:
                continue
            prime_val = row[0].strip()
            if prime_val and prime_val not in seen:
                seen.add(prime_val)
                unique_primes.append(prime_val)

                if len(unique_primes) >= MAX_PRIMES:
                    break

    input_buffer = "".join(f"{p} 0\n" for p in unique_primes)
    print(f"          Loaded {len(unique_primes)} unique primes into memory buffer.")

    print(f"[Step 2/3] Compiling '{c_source_file}'...")
    compile_cmd = ["gcc-15", "-O3", "-fopenmp", "-o", exec_name, c_source_file, "-lm", "-lpthread"]

    try:
        subprocess.run(compile_cmd, check=True)
    except subprocess.CalledProcessError:
        print("\n[Warning] OpenMP compilation failed. Trying single-threaded fallback...")
        fallback_cmd = ["gcc-15", "-O3", "-o", exec_name, c_source_file, "-lm", "-lpthread"]
        try:
            subprocess.run(fallback_cmd, check=True)
        except subprocess.CalledProcessError:
            print("[Fatal Error] Compilation failed completely.")
            return

    print(f"[Step 3/3] Executing benchmark via in-memory pipe...")

    total_cores = os.cpu_count() or 4
    target_threads = max(1, total_cores - 2)

    run_env = os.environ.copy()
    run_env["OMP_NUM_THREADS"] = str(target_threads)

    print(f"          [OpenMP] System cores: {total_cores} | Allocated threads: {target_threads}")
    run_cmd = [
        f"./{exec_name}",
        "/dev/stdin",
        out_pure_csv,
        out_metrics_csv,
        TARGET_TRIPLES_PER_PRIME
    ]

    try:
        subprocess.run(run_cmd, input=input_buffer, text=True, check=True, env=run_env)
        print(f"\n--- [Success] Benchmark completed for Method {method_enum} ---")
        print(f"  -> Pure triples : {out_pure_csv}")
        print(f"  -> Telemetry    : {out_metrics_csv}")
    except subprocess.CalledProcessError as e:
        print(f"\n[Error] C execution failed with exit code {e.returncode}.")
    finally:
        if os.path.exists(exec_name):
            os.remove(exec_name)

if __name__ == "__main__":
    main(5)
    main(6)
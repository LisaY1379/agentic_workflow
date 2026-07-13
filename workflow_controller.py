import os
import re
import json
import random
from dotenv import load_dotenv
import pandas as pd
from openai import OpenAI

# import vpp  # Uncomment if required by your specific local environment

# Force override system environment variables with .env file values
load_dotenv(override=True)

client = OpenAI(api_key=os.environ.get("OPENAI_API_KEY"))

# ==========================================
# CONFIGURATION & FILE PATHS
# ==========================================
DATA_FOLDER = "../data/final_data"
STRATEGY_FILE = "../output/filter_strategies.py"
LEDGER_FILE = "../output/ledger.json"

LLM_SAMPLE_CONFIGS = {
    "5_digits.csv": 50,
    "6_digits.csv": 50
}


# ==========================================
# 1. DATA LOADING & CONTEXT GENERATION
# ==========================================
def load_shuffled_csv_files(folder_path, file_configs):
    """
    Samples random rows from CSV files to format into the LLM context prompt.
    Extracts only the first 3 columns (prime, A, x0) to keep prompts clean.
    """
    combined_lines = []

    for file_name, num_rows in file_configs.items():
        full_path = os.path.join(folder_path, file_name)

        if not os.path.exists(full_path):
            print(f"⚠️ Warning: File '{file_name}' not found at {full_path}, skipping.")
            continue

        try:
            # Read only the first 3 columns: index 0 (prime), 1 (A), 2 (x0)
            df = pd.read_csv(full_path, usecols=[0, 1, 2])

            # Clean column names (strip whitespace and hidden UTF-8 BOM characters)
            df.columns = df.columns.str.strip().str.replace('\ufeff', '')

            sample_size = min(num_rows, len(df))
            df_sampled = df.sample(n=sample_size, random_state=42)
            df_sampled = df_sampled.sort_values(by=df_sampled.columns[0])

            combined_lines.append(f"--- Random Sampled Data from {file_name} ({sample_size} rows) ---")
            combined_lines.append(df_sampled.to_string(index=False))
            combined_lines.append("\n")

            print(f"🎲 Successfully sampled {sample_size} context rows from '{file_name}' for LLM prompt.")
        except Exception as e:
            print(f"❌ Error reading '{file_name}': {e}")

    return "\n".join(combined_lines)


def load_ground_truth_from_csv(folder_path, total_samples=100):
    """
    Loads ground-truth data from CSV files for empirical sandbox evaluation.
    Maps the CSV column 'prime' to the internal variable 'p'.
    """
    files = ["5_digits.csv", "6_digits.csv"]
    all_data = []

    for file_name in files:
        full_path = os.path.join(folder_path, file_name)
        if not os.path.exists(full_path):
            continue
        try:
            # Only load the first 3 columns to ignore trials, param_a, param_b
            df = pd.read_csv(full_path, usecols=[0, 1, 2])
            df.columns = df.columns.str.strip().str.replace('\ufeff', '')

            for _, row in df.iterrows():
                all_data.append({
                    'p': int(row['prime']),  # Explicitly map 'prime' column to internal 'p'
                    'A': int(row['A']),
                    'x0': int(row['x0'])
                })
        except Exception as e:
            print(f"❌ Error loading ground truth data from '{file_name}': {e}")

    # Fallback safety mechanism if CSV files fail to load
    if not all_data:
        print("⚠️ Warning: No CSV ground truth loaded. Using default backup dataset.")
        return {
            101: {'A': 4, 'x0': 53},
            103: {'A': 1, 'x0': 62},
            10007: {'A': 2, 'x0': 123},
            10009: {'A': 5, 'x0': 456},
            100003: {'A': 10, 'x0': 789},
            100019: {'A': 7, 'x0': 321}
        }

    # Randomly sample the requested number of ground truth instances
    sampled_data = random.sample(all_data, min(total_samples, len(all_data)))

    ground_truth_dict = {}
    for item in sampled_data:
        ground_truth_dict[item['p']] = {'A': item['A'], 'x0': item['x0']}

    print(f"🧪 Successfully loaded {len(ground_truth_dict)} real samples into the Ground Truth Evaluation Engine.")
    return ground_truth_dict


# ==========================================
# 2. STATE & MEMORY MANAGEMENT
# ==========================================
def get_next_filter_id(file_path):
    """Parses existing generated strategies to assign a sequential ID."""
    if not os.path.exists(file_path):
        return 1

    with open(file_path, "r", encoding="utf-8") as f:
        content = f.read()
        matches = re.findall(r"filter_(\d+)_", content)
        if not matches:
            return 1
        return max(int(m) for m in matches) + 1


def load_history_context(ledger_path=LEDGER_FILE):
    """
    Loads previously generated strategy descriptions from JSON ledger
    to prevent the LLM from generating cyclic or repetitive algorithms.
    """
    empty_message = "[History Context: The strategy library is currently empty. You have complete freedom to propose your first mathematical trick.]\n"

    if not os.path.exists(ledger_path):
        return empty_message

    try:
        with open(ledger_path, "r", encoding="utf-8") as f:
            ledger = json.load(f)

        if not ledger:
            return empty_message

        history_text = "[Explored Strategies Library - CRITICAL: Do NOT duplicate the mathematical logic of the following filters]\n"
        for item in ledger:
            strategy_name = item.get("name", "Unnamed_Strategy")
            description = item.get("description", "No description provided.")
            history_text += f"- Filter Name: {strategy_name}\n"
            history_text += f"  Mathematical Core: {description}\n\n"

        return history_text

    except Exception as e:
        print(f"⚠️ Warning: Failed to read ledger file at '{ledger_path}': {e}")
        return empty_message


# ==========================================
# 3. PROMPT ENGINEERING
# ==========================================
def get_system_prompt(next_id, context_data_string):
    history_context = load_history_context()

    return f"""You are an expert computational number theorist building an automated Pomerance search engine.

[Data Context]
Here is a sample dataset of primes and their valid Pomerance parameters (prime, A, x0):
{context_data_string}

{history_context}

[Task]
Autonomously discover a statistical or algebraic constraint to prune the search space for 'A'. 
Your output must be a single Python function `apply_filter(A, x0, p)` that returns True (keep) or False (discard).

[Strict Engineering Constraints]
1. NAMING: The function MUST be named exactly `filter_{next_id:02d}_[your_descriptive_name]`.
2. MODULARITY: Output pure filtering logic only. NO unbounded loops or recursive calls (to prevent execution hangs).
3. HEURISTICS ALLOWED: 
   - Heuristic filters are allowed! It is acceptable to falsely discard a small fraction of valid 'A' values (False Negatives), provided you drastically reduce the search space.
   - However, you MUST NOT discard 100% of the valid answers. If your filter rejects everything, it will cause an infinite loop in our C-engine and be rejected immediately.
4. DOCUMENTATION: Your function must begin with a docstring formatted exactly like this:
   \"\"\"
   # Name: [Descriptive Name]
   # Description: [Mathematical justification]
   \"\"\"

[Output Requirement]
Output ONLY the Python code block containing the single requested function enclosed in ```python ... ``` blocks.
"""


# ==========================================
# 4. SANDBOX EVALUATION ENGINE
# ==========================================
def evaluate_filter_logic(code_string, ground_truth_data):
    """
    Evaluates algorithmic performance using Double Sampling:
    1. Tests absolute accuracy against ground-truth pairs (False Negative Rate).
    2. Uses Monte Carlo sampling (up to 500 candidates per prime) for ultra-fast Pruning Rate estimation.
    """
    try:
        namespace = {}
        exec(code_string, namespace)

        ai_filter_func = None
        for name, func in namespace.items():
            if name.startswith("filter_") and callable(func):
                ai_filter_func = func
                break

        if not ai_filter_func:
            return False, "❌ Missing Function: Could not find a valid function definition matching 'filter_XX_...'."

        total_truths = len(ground_truth_data)
        missed_truths = 0
        total_pruning_rate = 0.0

        MONTE_CARLO_SAMPLES = 500  # Limits computation time to milliseconds

        # Core evaluation loop
        for p, truth in ground_truth_data.items():
            true_A = truth['A']
            true_x0 = truth['x0']

            # 1. Exact False Negative Test (Did we discard the true answer?)
            if not ai_filter_func(true_A, true_x0, p):
                missed_truths += 1

            # 2. Monte Carlo Efficiency Test (Estimate pruning percentage)
            sample_size = min(MONTE_CARLO_SAMPLES, p - 1)
            test_candidates = random.sample(range(1, p), sample_size)

            discarded_count = 0
            for test_A in test_candidates:
                if test_A == true_A:
                    continue
                # Evaluate pruning on random incorrect A values
                if not ai_filter_func(test_A, None, p):
                    discarded_count += 1

            total_pruning_rate += (discarded_count / sample_size)

        # Compute final aggregate metrics
        false_negative_rate = (missed_truths / total_truths) * 100.0
        avg_pruning = (total_pruning_rate / total_truths) * 100.0

        # Safety & quality rules
        if false_negative_rate == 100.0:
            return False, f"❌ FATAL (Infinite Loop Risk): 100.00% False Negative Rate! You discarded ALL valid answers. This will cause the C-engine to hang."

        if avg_pruning == 0.0:
            return False, "⚠️ Ineffective Strategy: 0.00% pruning rate. The filter provides no optimization."

        if false_negative_rate > 60.0:
            return False, f"⚠️ High False Negative: Pruned {avg_pruning:.2f}%, but the False Negative Rate is {false_negative_rate:.2f}%. This is too aggressive. Please relax your mathematical bounds!"

        if false_negative_rate > 0.0:
            return True, f"✅ HEURISTIC ACCEPTED: Pruned {avg_pruning:.2f}% of the search space. False Negative Rate: {false_negative_rate:.2f}% (acceptable engineering trade-off)."

        return True, f"✅ DETERMINISTIC SUCCESS: Perfect 0.00% False Negative Rate! Average search space pruned: {avg_pruning:.2f}%."

    except Exception as e:
        return False, f"❌ Execution Crash: {str(e)}"


# ==========================================
# 5. AGENT OPTIMIZATION LOOP
# ==========================================
def run_optimization_loop(iterations=3):
    print("🚀 Initializing Pomerance AI Discovery Pipeline...")

    # Load context string for LLM intuition
    context_data_string = load_shuffled_csv_files(DATA_FOLDER, LLM_SAMPLE_CONFIGS)

    # Load 100 real ground-truth samples for robust evaluation
    ground_truth_dataset = load_ground_truth_from_csv(DATA_FOLDER, total_samples=100)

    current_prompt = """Our current approach relies on randomly guessing `A` and running the verification loop. 
Please analyze the data and write an optimized filtering strategy. 
Remember, heuristic algorithms are allowed, but try to keep the False Negative rate low."""

    for i in range(iterations):
        print(f"\n==========================================")
        print(f"🔄 Agent Iteration Round {i + 1} / {iterations}")
        print(f"==========================================")

        current_next_id = get_next_filter_id(STRATEGY_FILE)
        dynamic_system_prompt = get_system_prompt(current_next_id, context_data_string)

        messages = [
            {"role": "system", "content": dynamic_system_prompt},
            {"role": "user", "content": current_prompt}
        ]

        print(f"🤖 Requesting new algorithm (Target ID: filter_{current_next_id:02d})...")
        try:
            response = client.chat.completions.create(
                model="gpt-5.4",
                messages=messages,
                temperature=0.2
            )
        except Exception as e:
            print(f"❌ OpenAI API Request Failed: {e}")
            break

        ai_response = response.choices[0].message.content
        code_match = re.search(r"```python(.*?)```", ai_response, re.DOTALL)

        if not code_match:
            print("⚠️ AI did not output a standard code block. Requesting correction...")
            current_prompt = "You did not provide the code enclosed in ```python ... ``` blocks. Please output valid Python code properly formatted."
            continue

        generated_code = code_match.group(1).strip()
        print("\n--- 📦 AI Generated Code Block ---")
        print(generated_code)
        print("----------------------------------\n")

        print("💡 Running local sandbox evaluation against 100 real ground truths...")
        success, result = evaluate_filter_logic(generated_code, ground_truth_dataset)

        if success:
            print(result)

            # Ensure directory structure exists
            os.makedirs(os.path.dirname(STRATEGY_FILE), exist_ok=True)

            annotated_code = f"# Eval Result: {result}\n{generated_code}"
            with open(STRATEGY_FILE, "a", encoding="utf-8") as f:
                f.write(f"\n\n# ==========================================\n")
                f.write(f"# Generated at Iteration: {i + 1} | Filter ID: {current_next_id:02d}\n")
                f.write(f"# ==========================================\n")
                f.write(annotated_code + "\n")
            print(f"💾 Successfully appended algorithm to '{STRATEGY_FILE}'")

            # Extract docstring description for the ledger
            desc_match = re.search(r"# Description:\s*(.*)", generated_code)
            strategy_desc = desc_match.group(1).strip() if desc_match else "Heuristic strategy."

            ledger_data = []
            if os.path.exists(LEDGER_FILE):
                with open(LEDGER_FILE, "r", encoding="utf-8") as lf:
                    try:
                        ledger_data = json.load(lf)
                    except Exception:
                        ledger_data = []

            ledger_data.append({
                "id": f"{current_next_id:02d}",
                "name": f"filter_{current_next_id:02d}_algorithm",
                "description": f"{strategy_desc} [{result.split(':')[0]}]"
            })

            with open(LEDGER_FILE, "w", encoding="utf-8") as lf:
                json.dump(ledger_data, lf, indent=2, ensure_ascii=False)
            print(f"📓 Logged strategy execution history to '{LEDGER_FILE}'")

            # Prompt AI to explore orthogonal ideas next
            current_prompt = f"Excellent! Your previous code was accepted with the following performance: {result}\nFor the next iteration, please propose an entirely orthogonal mathematical strategy or optimize your bounds to reduce the False Negative Rate further."

        else:
            print(f"⚠️ Validation Failed:\n{result}")
            # Feed the exact evaluation error back to the LLM for self-correction
            current_prompt = f"The generated code failed sandbox validation:\n{result}\nPlease mathematically relax your constraints or alter your algebraic formulation to reduce errors, ensuring you do not reject 100% of valid answers."

    print("\n🏁 Pipeline complete.")


if __name__ == "__main__":
    run_optimization_loop(iterations=3)
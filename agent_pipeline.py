import os
from dotenv import load_dotenv
from openai import OpenAI
import glob

load_dotenv(override=True)
client = OpenAI(api_key=os.environ["OPENAI_API_KEY"])

# =========================================================================
# 1. CONTEXT INGESTION ENGINE (With Three-Layer Blacklist Protection)
# =========================================================================
def gather_repo_context(repo_paths, file_extensions=(".md", ".py", ".c", ".h", ".gp", ".cu")):
    """
    Traverses the specified local clones of historical GitHub repositories,
    strictly screens out data dumps, logs, and massive files using blacklists,
    and packages valid source code and math notes into XML text blocks for LLM ingestion.
    """
    context_blocks = []
    ingested_files = 0
    skipped_files = 0

    # Layer 1: Filename Blacklist (catches giant raw proof datasets and verbose tables)
    file_blacklist = [
        "pp10.txt", "pp20.txt", "pp12.txt", "pp24.txt", "pp16A.txt",
        "oneshot8all.txt", "oneshot12prefixes.txt", "certs.csv",
        "pp14.txt", "pp16.txt", "HIT.txt", "RESULT.md", "worker.log", "verification.txt"
    ]

    # Layer 2: Directory Blacklist (blocks folders containing thousands of worker logs/shards)
    dir_blacklist = [
        "data", "logs", "shards", "test", "val_p7", "certs", "results",
        "venv", ".env", "node_modules", ".git", "__pycache__", "build", "out"
    ]

    print("🔍 [Stage 0] Ingesting repository assets with strict blacklist screening...")

    for repo_dir in repo_paths:
        if not os.path.exists(repo_dir):
            print(f"⚠️ [Skip] Directory does not exist: {repo_dir}")
            continue

        for ext in file_extensions:
            for filepath in glob.glob(f"{repo_dir}/**/*{ext}", recursive=True):

                # Check Layer 2: Directory Path Filtering
                if any(f"/{ignored}/" in filepath or f"\\{ignored}\\" in filepath or filepath.endswith(f"/{ignored}")
                       for ignored in dir_blacklist):
                    continue

                # Check Layer 1: Filename Filtering
                filename = os.path.basename(filepath)
                if filename in file_blacklist:
                    skipped_files += 1
                    continue

                try:
                    file_size = os.path.getsize(filepath)
                    # Check Layer 3: File-Size Circuit Breaker (>200KB is almost always a data log, not code)
                    if file_size > 200 * 1024:
                        print(f"⏭️ [Circuit Breaker: >200KB]: Skipped {filename} ({file_size // 1024} KB)")
                        skipped_files += 1
                        continue

                    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
                        content = f.read()
                        context_blocks.append(f"<file path='{filepath}'>\n{content}\n</file>")
                        ingested_files += 1

                except Exception as e:
                    print(f"[Warning] Failed to read asset file {filepath}: {e}")

    print(f"✅ Ingestion Complete: {ingested_files} core files loaded | 🗑️ {skipped_files} data/log files blocked.")
    return "\n\n".join(context_blocks)


# =========================================================================
# 2. STAGE 1: THEORIST AGENT (Knowledge Synthesis & Strategy Extractor)
# =========================================================================
def extract_winning_strategies(history_context):
    """
    Spawns the Theorist Agent to digest the historical context and generate
    a machine-readable specification of mathematically validated filtering tactics,
    strictly prioritizing zero-search algorithmic shortcuts over brute-force optimizations.
    """
    print("🧠 [Stage 1] Launching Theorist Agent [Model: gpt-5.6-sol]...")
    print("👉 Task: Distilling hierarchical SOTA methods (prioritizing OneShot shortcuts)...")

    prompt = f"""
    You are a world-class expert in algebraic number theory, elliptic curve cryptography (ECC), and high-performance computational mathematics. 
    Below is a comprehensive collection of research ledger logs, proof markdown notes, and experimental source codes ingested from our historical repositories (including DANGER3, short-certificate experiments, and OneShotPrimalityProofs):

    {history_context}

    Task:
    Analyze these documents to distill all mathematically validated acceleration strategies for finding Pomerance triples and One-Shot ECPP certificates. 

    CRITICAL ARCHITECTURAL DIRECTIVE:
    You must evaluate and categorize strategies using a strict execution hierarchy. We want to bypass heavy computation entirely whenever possible:

    1. [PRIORITY TIER 1: Zero-Search / Algorithmic Shortcuts]: Look first for methods that eliminate search entirely (e.g., the Supersingular Shortcut where trace-0 curves like y^2 = x^3 + x with A=0 yield smooth group orders p+1, or Reverse-CM order generation).
    2. [PRIORITY TIER 2: One-Shot ECPP Search Space Expansion]: Methods that relax the strict 2^k torsion requirement to broader n^4-smooth group orders (m), expanding admissible Hasse interval targets from ~3 to thousands.
    3. [PRIORITY TIER 3: C-Engine Hot-Loop Pruning]: Traditional parallel search optimizations for ordinary curves (e.g., X1(16) prescribed torsion parameterization, Montgomery split/nonsplit Legendre symbol filtering, early 2-Sylow halving pre-checks, and p mod 4 arithmetic patches).

    You must output a strictly formatted JSON object matching the schema below. Do not include conversational prose outside the JSON.

    Expected JSON schema:
    {{
      "execution_hierarchy": [
        {{
          "priority_tier": 1 | 2 | 3,
          "method_name": "Name of the strategy",
          "is_zero_search_shortcut": true | false,
          "algebraic_logic": "The underlying mathematical, group-theoretic, or smooth-order residue principle",
          "pre_flight_check_python": "If is_zero_search_shortcut is true, describe the exact mathematical check to run in Python before launching C (e.g., checking if the n^4-smooth part of p+1 exceeds L). If false, leave null.",
          "c_implementation_advice": "Detailed guidance on how to integrate this check directly into the C-based hot loop (e.g., pre-computation steps, early break/continue conditions, or mod 4 patch handling). Leave null if it is a pure Python pre-flight shortcut."
        }}
      ]
    }}
    """

    response = client.chat.completions.create(
        model="gpt-5.6-sol",
        messages=[{"role": "user", "content": prompt}],
        response_format={"type": "json_object"},
        temperature=1
    )
    return response.choices[0].message.content


# =========================================================================
# 3. STAGE 2: SYSTEM ENGINEER AGENT (Low-Level C-Engine Optimizer)
# =========================================================================
def optimize_c_engine(target_c_file, strategies_json):
    """
    Spawns the Systems Work Agent to map the JSON specification into the
    production C codebase, injecting algorithmic filters and hardware optimizations.
    """
    print("⚡ [Stage 2] Launching C-Engineer Agent: Injecting synthesized strategies into C-core...")

    with open(target_c_file, "r", encoding="utf-8") as f:
        raw_c_code = f.read()

    prompt = f"""
    You are a veteran systems engineer and high-performance computing (HPC) expert specialized in low-level mathematical code optimization.

    We have extracted the following algorithmic acceleration strategies from prior literature and experiments:
    {strategies_json}

    Below is our active production-grade C source code responsible for running massive parallel search loops to locate valid Pomerance triples:
    <target_c_code>
    {raw_c_code}
    </target_c_code>

    Your optimization must strictly satisfy these conditions:
    1. [Theoretical Pruning]: Inject the low-overhead algebraic filtering conditions from the strategies JSON right before the most compute-heavy steps in the hot loop. For candidates that fail the algebraic criteria, trigger an immediate 'continue' to minimize trial latency.
    2. [Hardware Scaling]: Maintain or enhance existing OpenMP multi-threading primitives, ensuring loop-invariant extractions, cache locality, and atomic operations are pristine.
    3. [Compilation Integrity]: Return the full, complete rewritten C source file. Do not truncate code using comments like "// ... keep existing logic ...". It must be compilable immediately using standard gcc/clang flags.
    """

    response = client.chat.completions.create(
        model="gpt-5.6-terra",
        messages=[{"role": "user", "content": prompt}],
        temperature=1
    )
    return response.choices[0].message.content


# =========================================================================
# 4. EXECUTION ORCHESTRATION CONTROLLER
# =========================================================================
if __name__ == "__main__":

    # Define absolute paths to your locally cloned historical Git repositories

    base_dir = os.path.expanduser("~/Documents")

    repositories = [
        os.path.join(base_dir, "pomerance-p24"),
        os.path.join(base_dir, "OneShotPrimalityProofs"),
        os.path.join(base_dir, "gpu_optimization")
    ]

    # Target C production code to re-architect (Updated to use pomerance-p24 since DANGER3 was dropped)
    c_source_target = os.path.join(base_dir, "pomerance-p24/src/pomerance_fixed.c")
    c_source_output = os.path.join(base_dir, "pomerance-p24/src/pomerance_old.c")

    # Execute Context Ingestion Flow
    context_payload = gather_repo_context(repositories)
    print(f"📦 Context payload compiled ({len(context_payload):,} characters ready for API dispatch).")

    # Run Stage 1: Strategy Synthesis
    extracted_specification = extract_winning_strategies(context_payload)
    print("\n💡 Synthesized Strategy Specification (JSON):\n", extracted_specification)

    json_output_path = os.path.join(base_dir, "pomerance-p24/src/extracted_strategies.json")
    with open(json_output_path, "w", encoding="utf-8") as jf:
        jf.write(extracted_specification)
    print(f"📁 JSON strategy file saved to: {json_output_path}")
    # Run Stage 2: Codebase Optimization
    final_c_codebase = optimize_c_engine(c_source_target, extracted_specification)

    # Flush optimized assets to production track
    with open(c_source_output, "w", encoding="utf-8") as f:
        f.write(final_c_codebase)

    print(f"\n🚀 Pipeline Completed! Optimized C-engine has been deployed to: {c_source_output}")
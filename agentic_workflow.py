import os
import sys
import json
import subprocess
from openai import OpenAI
from dotenv import load_dotenv

load_dotenv(override=True)

client = OpenAI(api_key=os.environ.get("OPENAI_API_KEY"))


# ==========================================
# 1. Base Tools (Worker-Level Capabilities)
# ==========================================
def list_files(directory: str = "."):
    """Lists all files in the given directory (default is current directory)."""
    try:
        files = [f for f in os.listdir(directory) if
                 os.path.isfile(os.path.join(directory, f)) and not f.startswith('.')]
        return json.dumps({"files": files}, ensure_ascii=False)
    except Exception as e:
        return f"Error listing directory: {str(e)}"


def read_local_file(file_name: str):
    """Reads the content of ANY local file EXCEPT CSV files."""
    if file_name.lower().endswith('.csv'):
        return "Access Denied: Direct reading of .csv files is strictly PROHIBITED. Please write and execute a Python script to inspect or process CSV data."

    if os.path.exists(file_name):
        try:
            with open(file_name, 'r', encoding='utf-8') as f:
                return f.read()
        except UnicodeDecodeError:
            return "Error: Cannot read binary file (not text-based)."
        except Exception as e:
            return f"Error reading file: {str(e)}"
    return "Error: File not found."


def create_or_write_file(file_name: str, content: str):
    """Creates a new file (e.g., .py) or overwrites an existing one with the given content."""
    try:
        safe_name = os.path.basename(file_name)
        with open(safe_name, 'w', encoding='utf-8') as f:
            f.write(content)
        return f"Success: File '{safe_name}' has been created/updated successfully."
    except Exception as e:
        return f"Error writing file: {str(e)}"


def run_python_script(file_name: str):
    """Executes a Python script and returns its stdout and stderr."""
    if not os.path.exists(file_name):
        return f"Error: File '{file_name}' does not exist."

    if not file_name.lower().endswith('.py'):
        return "Error: You can only execute '.py' files with this tool."

    try:
        result = subprocess.run(
            [sys.executable, file_name],
            capture_output=True,
            text=True,
            encoding='utf-8',
            timeout=30
        )

        output = ""
        if result.stdout:
            output += f"[Standard Output]:\n{result.stdout}\n"
        if result.stderr:
            output += f"[Error Output]:\n{result.stderr}\n"

        if not output.strip():
            return "Success: Script executed successfully with no printed output."

        return output
    except subprocess.TimeoutExpired:
        return "Error: Execution timed out (exceeded 30 seconds). Check if the code has an infinite loop."
    except Exception as e:
        return f"Error executing script: {str(e)}"


# Tools accessible ONLY by the Specialist Sub-Agents (Workers)
worker_tools_metadata = [
    {
        "type": "function",
        "function": {
            "name": "list_files",
            "description": "List all files in the current working directory to see what is available.",
            "parameters": {"type": "object", "properties": {"directory": {"type": "string",
                                                                          "description": "The directory path to list files from. Default is '.' (current directory)."}}}
        }
    },
    {
        "type": "function",
        "function": {
            "name": "read_local_file",
            "description": "Read the text content of local files (e.g., .py, .txt, .json, .md). STRICTLY PROHIBITED FOR .csv FILES. Do NOT pass .csv files to this tool.",
            "parameters": {"type": "object", "properties": {
                "file_name": {"type": "string", "description": "The name of the file to read (EXCLUDING .csv files)."}},
                           "required": ["file_name"]}
        }
    },
    {
        "type": "function",
        "function": {
            "name": "create_or_write_file",
            "description": "Create a new Python file (.py) or any text file, or overwrite an existing file with new content.",
            "parameters": {"type": "object", "properties": {
                "file_name": {"type": "string", "description": "The name of the file to create/save."},
                "content": {"type": "string", "description": "The full code or text content to write into the file."}},
                           "required": ["file_name", "content"]}
        }
    },
    {
        "type": "function",
        "function": {
            "name": "run_python_script",
            "description": "Execute a local Python script (.py) and get its output or error messages. Essential for running data analysis on CSV files or testing code.",
            "parameters": {"type": "object", "properties": {
                "file_name": {"type": "string", "description": "The name of the Python file to execute."}},
                           "required": ["file_name"]}
        }
    }
]

WORKER_TOOLS_MAP = {
    "list_files": list_files,
    "read_local_file": read_local_file,
    "create_or_write_file": create_or_write_file,
    "run_python_script": run_python_script
}


# ==========================================
# 2. Sub-Agent Engine (The Worker Sandbox)
# ==========================================
def run_sub_agent(role_name: str, system_prompt: str, task_description: str):
    """
    Spins up an autonomous specialist sub-agent to execute a specific task.
    It operates in its own isolated context and returns only the final summary report.
    """
    print(f"\n🏢 [Director Delegation] Dispatching task to specialist: 【{role_name}】...")
    print(f"📋 [Task Directive]: {task_description}")

    # 1. Isolated Context: This prevents worker tracebacks/long outputs from polluting the PM's chat history!
    sub_history = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": task_description}
    ]

    # 2. Worker Loop: Allow up to 5 iterations for the sub-agent to write code, test, and self-correct errors.
    for _ in range(5):
        response = client.chat.completions.create(
            model="gpt-4o",  # Can be swapped to gpt-4o-mini to save tokens for routine worker tasks
            messages=sub_history,
            tools=worker_tools_metadata
        )
        msg = response.choices[0].message
        sub_history.append(msg)

        if msg.tool_calls:
            for tc in msg.tool_calls:
                func_name = tc.function.name
                args = json.loads(tc.function.arguments)
                print(f"   🛠️ 【{role_name}】 is calling [{func_name}]...")

                func_to_call = WORKER_TOOLS_MAP.get(func_name)
                res = func_to_call(**args) if func_to_call else f"Error: Tool {func_name} not found."

                sub_history.append({
                    "role": "tool",
                    "tool_call_id": tc.id,
                    "name": func_name,
                    "content": str(res)
                })
        else:
            # When the sub-agent stops calling tools, it means it reached a conclusion!
            print(f"✅ 【{role_name}】 has completed the assignment!")
            return f"[{role_name}'s Report]: {msg.content}"

    return f"[{role_name} Execution Timeout]: Failed to complete the task within 5 autonomous iterations. Please re-evaluate the requirements."


# ==========================================
# 3. PM Management Tool (Delegation)
# ==========================================
def delegate_to_specialist(role_name: str, task_description: str):
    """The tool used by the PM Agent to summon and command specialized sub-agents."""
    prompts = {
        "Software Engineer": (
            "You are a Senior Software Engineer. Your sole objective is to write clean, robust, and well-tested code, "
            "save it to the filesystem, and execute it to verify it works. Always fix bugs independently if scripts fail."
        ),
        "Data Analyst": (
            "You are an expert Data Analyst specializing in Python. You excel at extracting actionable business insights from data. "
            "CRITICAL RULE: You are strictly forbidden from reading '.csv' files directly with `read_local_file`. "
            "You must write and execute Python scripts (using pandas/csv) to analyze CSV data and report the results."
        ),
        "QA Tester": (
            "You are a meticulous QA Testing Engineer. Your role is to inspect existing code, write comprehensive test cases, "
            "execute them, and expose edge cases or potential bugs."
        )
    }

    sys_prompt = prompts.get(role_name,
                             "You are a specialized AI assistant. Use your tools to accomplish the assigned task.")
    return run_sub_agent(role_name, sys_prompt, task_description)


# The PM only has management tools, no direct coding/file access!
pm_tools_metadata = [
    {
        "type": "function",
        "function": {
            "name": "delegate_to_specialist",
            "description": "Delegate a technical task (coding, data analysis, testing, file operations) to a specialized worker sub-agent. Do NOT attempt to write code or read files yourself; always delegate to the appropriate specialist.",
            "parameters": {
                "type": "object",
                "properties": {
                    "role_name": {
                        "type": "string",
                        "enum": ["Software Engineer", "Data Analyst", "QA Tester"],
                        "description": "The job title of the specialist best suited for the task."
                    },
                    "task_description": {
                        "type": "string",
                        "description": "Clear, detailed instructions and requirements for what the specialist needs to accomplish."
                    }
                },
                "required": ["role_name", "task_description"]
            }
        }
    }
]


# ==========================================
# 4. Executive PM Agent (The Mastermind)
# ==========================================
def discuss_with_pm():
    print(
        "🤖 PM Agent: Ready. I am your executive Product Manager. I design experiments and manage specialists to get the work done.")

    chat_history = [{
        "role": "system",
        "content": (
            "You are an elite AI Product Manager and System Architect. You act as the bridge between the human user and technical execution.\n"
            "You DO NOT write code, read files, or execute scripts directly. Instead, you analyze the user's high-level goals, "
            "break them down into logical technical steps, and delegate them to your specialist team using `delegate_to_specialist`.\n"
            "Available Specialists:\n"
            "- 'Data Analyst': For data inspection, metrics calculation, and CSV file processing via scripts.\n"
            "- 'Software Engineer': For creating applications, writing Python scripts, generating HTML/JS, or building tools.\n"
            "- 'QA Tester': For testing scripts, verifying data accuracy, and finding bugs.\n"
            "Synthesize the reports from your specialists and present a clear, executive-level summary to the user."
        )
    }]

    while True:
        try:
            user_input = input("\n👤 You: ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\n🤖 PM Agent: Conversation interrupted. Goodbye!")
            break

        if not user_input:
            continue

        if user_input.lower() in ["ok", "exit", "quit", "goodbye", "done"]:
            print("🤖 PM Agent: Sounds good. Ending the conversation!")
            break

        chat_history.append({"role": "user", "content": user_input})

        response = client.chat.completions.create(
            model="gpt-4o",
            messages=chat_history,
            tools=pm_tools_metadata
        )

        msg = response.choices[0].message

        if msg.tool_calls:
            chat_history.append(msg)

            for tool_call in msg.tool_calls:
                func_name = tool_call.function.name
                args = json.loads(tool_call.function.arguments)

                print(f"🎯 PM Agent is delegating via: [{func_name}] with args: {args}...")

                if func_name == "delegate_to_specialist":
                    content = delegate_to_specialist(**args)
                else:
                    content = f"Error: Tool {func_name} not found."

                chat_history.append({
                    "role": "tool",
                    "tool_call_id": tool_call.id,
                    "name": func_name,
                    "content": str(content)
                })

            # Get the final synthesis from the PM after receiving worker reports
            final_res = client.chat.completions.create(model="gpt-4o", messages=chat_history)
            print(f"\n🤖 PM Agent: {final_res.choices[0].message.content}")
            chat_history.append(final_res.choices[0].message)
        else:
            print(f"\n🤖 PM Agent: {msg.content}")
            chat_history.append(msg)


if __name__ == "__main__":
    discuss_with_pm()
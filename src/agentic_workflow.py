import os
import sys
import json
import subprocess
from openai import OpenAI
from dotenv import load_dotenv

load_dotenv(override=True)

client = OpenAI(api_key=os.environ.get("OPENAI_API_KEY"))


# ==========================================
# 1. Base Tools (For Specialist Sub-Agents)
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
        return "Access Denied: Direct reading of .csv files is strictly PROHIBITED. Please write and execute a Python script or shell command to inspect CSV data."

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


def execute_terminal_command(command: str):
    """Executes a terminal/shell command (e.g., bash/cmd) and returns its stdout and stderr."""
    forbidden_keywords = ["rm -rf /", "rm -rf ~", "mkfs", "dd if=", "shutdown", "reboot", "> /dev/null"]
    for forbidden in forbidden_keywords:
        if forbidden in command.lower():
            return f"Access Denied: The command containing '{forbidden}' is strictly blocked for security reasons."

    try:
        result = subprocess.run(
            command,
            shell=True,
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='replace',
            timeout=60
        )

        output = ""
        if result.stdout:
            output += f"[Standard Output]:\n{result.stdout}\n"
        if result.stderr:
            output += f"[Standard Error]:\n{result.stderr}\n"

        if not output.strip():
            return "Success: Command executed successfully with no output returned."

        if len(output) > 3000:
            output = output[:3000] + "\n...[Output truncated due to character limit. Modify command to narrow output.]"

        return output
    except subprocess.TimeoutExpired:
        return "Error: Command execution timed out (exceeded 60 seconds). Avoid interactive commands or long-running daemons."
    except Exception as e:
        return f"Error executing terminal command: {str(e)}"


worker_tools_metadata = [
    {
        "type": "function",
        "function": {
            "name": "list_files",
            "description": "List all files in the current working directory to see what is available.",
            "parameters": {"type": "object", "properties": {"directory": {"type": "string", "description": "The directory path. Default is '.'."}}}
        }
    },
    {
        "type": "function",
        "function": {
            "name": "read_local_file",
            "description": "Read text files (e.g., .py, .txt, .json). STRICTLY PROHIBITED FOR .csv FILES.",
            "parameters": {"type": "object", "properties": {"file_name": {"type": "string", "description": "Name of the file to read."}}, "required": ["file_name"]}
        }
    },
    {
        "type": "function",
        "function": {
            "name": "create_or_write_file",
            "description": "Create a new file or overwrite an existing text/code file.",
            "parameters": {"type": "object", "properties": {"file_name": {"type": "string", "description": "Name of the file."}, "content": {"type": "string", "description": "Full text content."}}, "required": ["file_name", "content"]}
        }
    },
    {
        "type": "function",
        "function": {
            "name": "execute_terminal_command",
            "description": "Execute any bash/terminal command (e.g., 'pip install pandas', 'python script.py', 'git status', 'ls -la', 'head -n 20 data.csv').",
            "parameters": {
                "type": "object",
                "properties": {
                    "command": {"type": "string", "description": "The shell command to execute in the local terminal."}
                },
                "required": ["command"]
            }
        }
    }
]

WORKER_TOOLS_MAP = {
    "list_files": list_files,
    "read_local_file": read_local_file,
    "create_or_write_file": create_or_write_file,
    "execute_terminal_command": execute_terminal_command
}


# ==========================================
# 2. Sub-Agent Engine (Worker Sandbox)
# ==========================================
def run_sub_agent(role_name: str, system_prompt: str, task_description: str):
    """Spins up an autonomous specialist sub-agent to execute a task using worker tools."""
    print(f"\n🏢 [Director Delegation] Dispatching task to specialist: 【{role_name}】...")
    print(f"📋 [Task Directive]: {task_description}")

    sub_history = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": task_description}
    ]

    for _ in range(5):
        response = client.chat.completions.create(
            model="gpt-4o-mini",
            messages=sub_history,
            tools=worker_tools_metadata
        )
        msg = response.choices[0].message
        sub_history.append(msg)

        if msg.tool_calls:
            for tc in msg.tool_calls:
                func_name = tc.function.name
                args = json.loads(tc.function.arguments)
                print(f"   🛠️ 【{role_name}】 is executing [{func_name}] with args: {args}...")

                func_to_call = WORKER_TOOLS_MAP.get(func_name)
                res = func_to_call(**args) if func_to_call else f"Error: Tool {func_name} not found."

                sub_history.append({
                    "role": "tool",
                    "tool_call_id": tc.id,
                    "name": func_name,
                    "content": str(res)
                })
        else:
            print(f"✅ 【{role_name}】 has completed the assignment!")
            return f"[{role_name}'s Report]: {msg.content}"

    return f"[{role_name} Execution Timeout]: Failed to complete the task within 5 iterations."


# ==========================================
# 3. PM Restricted Sandboxed Tools (Upgraded!)
# ==========================================
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def delegate_to_specialist(role_name: str, task_description: str):
    """The tool used by the PM Agent to command specialized sub-agents."""
    prompts = {
        "Software Engineer": (
            "You are a Senior DevOps and Software Engineer with full terminal access via `execute_terminal_command`. "
            "You can install pip packages, run scripts, manage git, and inspect environment setups. "
            "Always fix errors independently and verify shell command outputs."
        ),
        "Data Analyst": (
            "You are an expert Data Analyst with terminal access. "
            "CRITICAL RULE: Never use `read_local_file` for CSVs. Instead, write Python scripts or use shell commands like `head` or `awk` via `execute_terminal_command` to inspect data."
        ),
        "QA Tester": (
            "You are a QA Engineer. You can run test suites (e.g., `pytest`), execute scripts, and verify software output using terminal commands."
        )
    }
    sys_prompt = prompts.get(role_name, "You are a specialized AI assistant with terminal access.")
    return run_sub_agent(role_name, sys_prompt, task_description)


# 【升级 1】：同时列举 instructions/ 和 plans/ 下的文件
def list_pm_files():
    """PM Sandbox: Lists all available files inside BOTH 'instructions/' and 'plans/' directories."""
    result = {}
    for folder in ["instructions", "plans"]:
        target_dir = os.path.join(PROJECT_ROOT, folder)
        if os.path.exists(target_dir):
            try:
                files = [f for f in os.listdir(target_dir) if os.path.isfile(os.path.join(target_dir, f)) and not f.startswith('.')]
                result[folder] = files
            except Exception as e:
                result[folder] = f"Error: {str(e)}"
        else:
            result[folder] = "Directory does not exist yet."
    return json.dumps(result, ensure_ascii=False)


# 【升级 2】：允许安全读取两个目录下的文件
def read_pm_file(directory_name: str, file_name: str):
    """PM Sandbox: Reads files ONLY from either 'instructions/' or 'plans/' directories."""
    if directory_name not in ["instructions", "plans"]:
        return "Access Denied: You are STRICTLY restricted to reading from 'instructions' or 'plans' directories ONLY."

    safe_name = os.path.basename(file_name)
    target_path = os.path.join(PROJECT_ROOT, directory_name, safe_name)

    if not os.path.exists(target_path):
        return f"Access Denied or Error: File '{safe_name}' not found in '{directory_name}/'. Use `list_pm_files` first!"
    try:
        with open(target_path, 'r', encoding='utf-8') as f:
            return f.read()
    except Exception as e:
        return f"Error reading file: {str(e)}"


def write_plan_file(file_name: str, content: str):
    """PM Sandbox: Writes .md files ONLY into the 'plans/' directory."""
    safe_name = os.path.basename(file_name)
    if not safe_name.endswith(".md"):
        safe_name += ".md"

    plans_dir = os.path.join(PROJECT_ROOT, "plans")
    os.makedirs(plans_dir, exist_ok=True)
    target_path = os.path.join(plans_dir, safe_name)

    try:
        with open(target_path, 'w', encoding='utf-8') as f:
            f.write(content)
        return f"Success: Plan saved to '{target_path}'."
    except Exception as e:
        return f"Error writing plan file: {str(e)}"


pm_tools_metadata = [
    {
        "type": "function",
        "function": {
            "name": "delegate_to_specialist",
            "description": "Delegate a technical or environment task to a specialist sub-agent. Use this whenever you need to touch files outside 'instructions/' or 'plans/'.",
            "parameters": {
                "type": "object",
                "properties": {
                    "role_name": {"type": "string", "enum": ["Software Engineer", "Data Analyst", "QA Tester"], "description": "The job title of the specialist best suited for the task."},
                    "task_description": {"type": "string", "description": "Detailed instructions and requirements for the task."}
                },
                "required": ["role_name", "task_description"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "list_pm_files",
            "description": "List all existing files inside BOTH 'instructions/' and 'plans/' directories. ALWAYS call this before attempting to read files!",
            "parameters": {"type": "object", "properties": {}}
        }
    },
    {
        "type": "function",
        "function": {
            "name": "read_pm_file",
            "description": "Read operational guidelines or existing plans. STRICTLY restricted to reading from 'instructions' or 'plans' directories.",
            "parameters": {
                "type": "object",
                "properties": {
                    "directory_name": {"type": "string", "enum": ["instructions", "plans"], "description": "The folder containing the file."},
                    "file_name": {"type": "string", "description": "Name of the file to read."}
                },
                "required": ["directory_name", "file_name"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "write_plan_file",
            "description": "Write your ideas, state graphs, and plans ONLY into the 'plans/' directory as a .md file.",
            "parameters": {
                "type": "object",
                "properties": {
                    "file_name": {"type": "string", "description": "Name of the markdown file to create inside 'plans/' (e.g., 'state_graph_v1.md')."},
                    "content": {"type": "string", "description": "The full markdown content of the plan or state graph."}
                },
                "required": ["file_name", "content"]
            }
        }
    }
]

PM_TOOLS_MAP = {
    "delegate_to_specialist": delegate_to_specialist,
    "list_pm_files": list_pm_files,
    "read_pm_file": read_pm_file,
    "write_plan_file": write_plan_file
}


# ==========================================
# 4. Executive PM Agent
# ==========================================
def discuss_with_pm():
    print("🤖 PM Agent: Ready. I have read access to instructions/ & plans/, and write access to plans/. All code execution is delegated to my specialists.")

    pm_system_prompt = (
        "You are an Elite AI Product Manager and Architect.\n\n"
        "=== REPOSITORY LAYOUT & PERMISSION GUIDELINES ===\n"
        "├── src/\n"
        "├── data/\n"
        "├── reports/\n"
        "├── instructions/ # You only have READ access via `read_pm_file`\n"
        "├── plans/        # You have READ & WRITE access via `read_pm_file` and `write_plan_file` (.md files)\n"
        "└── README.md\n\n"
        "CRITICAL RULES:\n"
        "1. Sandboxed Access: You can ONLY read files from `instructions/` or `plans/`, and write your ideas/plans as `.md` files in `plans/`.\n"
        "2. Brain-Hand Separation: If you need to read/write files in `src/`, `data/`, `reports/`, inspect root files, or execute any terminal commands, YOU CANNOT DO IT YOURSELF. You MUST call a subordinate via `delegate_to_specialist`.\n"
        "3. Strategy First: Always check existing files via `list_pm_files` and read relevant guidelines/plans if prompted. Save your state graphs or experiment designs into `plans/` before dispatching heavy tasks to specialists."
    )

    chat_history = [{"role": "system", "content": pm_system_prompt}]

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

                print(f"🎯 PM Agent is calling: [{func_name}] with args: {args}...")

                func_to_call = PM_TOOLS_MAP.get(func_name)
                content = func_to_call(**args) if func_to_call else f"Error: Tool {func_name} not found."

                chat_history.append({
                    "role": "tool",
                    "tool_call_id": tool_call.id,
                    "name": func_name,
                    "content": str(content)
                })

            final_res = client.chat.completions.create(model="gpt-4o", messages=chat_history)
            print(f"\n🤖 PM Agent: {final_res.choices[0].message.content}")
            chat_history.append(final_res.choices[0].message)
        else:
            print(f"\n🤖 PM Agent: {msg.content}")
            chat_history.append(msg)


if __name__ == "__main__":
    discuss_with_pm()
# Agentic Workflow

This is a project investigating how agentic workflow could help in mathematical discovery, from optimizing algorithms to investigating certain mathematical phenomena.

## Key Takeaway

Instead of closed-loop end-to-end "automatic research," it has been discovered that the key to a successful research workflow lies in:
* Problem definition & benchmark design;
* Sufficient background information & context;
* A workflow that allows interactions, feedback & corrections between steps.

Two problems have been tested through the established agentic workflow:

## Problem 1: Generate a new method that could accelerate finding Pomerance Triples

This problem tests the agent's ability in optimizing algorithms under a given framework (Pomerance Primality Proofs) and benchmark (code efficiency, trials & running time). Under this workflow, a new method (using Barrett reduction in modular divisions) has been proposed and been tested with a ~1.7-2.5x speed boost, with all results verified as valid triples.

* More details about the problem description in: 

    `problems/problems1_generate_new_method/input/PROJECT_VISION.md`

* Evaluation of the new method in:

    `workflow_logs/problem1_generate_new_method/reports/method_eval_13d.txt`: Tested established algorithms (method 2 & 3) with and without method 4 (new method). Results show method 4's ability in running-time acceleration.

    `workflow_logs/problem1_generate_new_method/reports/10^24+7.txt`: Tested method 2+1+4 on 10^24+7, found a proof in 2h 21min (7221920.97ms), about 2.5x speed of previous fastest code 5h 49min.

## Problem 2: Investigate Stage Jumps in the Difficulty of finding Pomerance Proofs

The problem tries to find an explanation towards the sudden changes of difficulty in finding Pomerance Proofs through this workflow. After the description of the problem, an explanation (see `problem2_investigate_stage_jump/output/PROPOSAL1.md`) has been proposed. The explanation is considered to be valid after being tested on the empirical plots (see `problem2_investigate_stage_jump/output/plots/`).

## Workflow Description

### 1. Problem Description & Benchmark

The workflow starts with describing your problem/goal and a testing benchmark, which will later being read by an agent and generate initial hypothesis/plans.

### 2. Planning & Proposals

The agent generates a high-level planning regarding this problem.

### 3. Experiment Design

Test the agent's proposal on executable code/empirical data.

### 4. Validate/Reject
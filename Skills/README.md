# Skills Repository (2026 Edition)

> **The "Microservices Moment" for Biomedical AI Agents.**

![Status](https://img.shields.io/badge/Status-Active-green)
![Agents](https://img.shields.io/badge/Agents-Orchestrated-blue)
![Domain](https://img.shields.io/badge/Domain-Biotech%20%7C%20Clinical%20%7C%20Genomics-purple)
![Tech](https://img.shields.io/badge/Tech-MCP%20%7C%20DeepSeek%20%7C%20Gemini-orange)

## 🚀 Overview

This repository is a comprehensive library of **skills, agents, and mathematical foundations** for modern (2026) Artificial Intelligence. Unlike standard chatbot repos, this project focuses on **Agentic Workflows**—where autonomous systems plan, execute, use tools, and correct themselves to solve complex scientific problems.

We have aligned this codebase with the **State of the Art (SOTA) for 2026**, integrating Agentic patterns, Model Context Protocol (MCP), and rigorous scientific simulations.

## 🌟 Key Capabilities (New for 2026)

### 🧬 Genomics & Single Cell (New!)
*   **Universal Annotator:** `Genomics/Single_Cell/Cell_Type_Annotation/RNA/universal_annotator.py` wraps Marker-based, Deep Learning (CellTypist), and LLM-based annotation strategies.
*   **Pathway Scoring:** `Genomics/Single_Cell/Pathway_Analysis/sc_pathway_scorer.py` implements AUCell-like scoring for functional enrichment.
*   **Cell-Cell Comms:** `Genomics/Single_Cell/Cell_Cell_Communication/interaction_inference.py` infers Ligand-Receptor networks.
*   **Database:** A curated [Tool Database](Genomics/Single_Cell/Tool_Database.md) of 2026 single-cell tools (MultiKano, scPS, etc.).

### 🧠 Agentic AI (The Brain)
*   **Orchestrated Swarms:** `Agentic_AI/Multi_Agent_Systems/orchestrator.py` implements a Supervisor pattern that delegates tasks to specialized sub-agents (Coder, Chemist, Reviewer).
*   **Plan-and-Solve:** `Agentic_AI/Agent_Architectures/Plan_and_Solve/` breaks down complex user queries into Directed Acyclic Graphs (DAGs).
*   **Async Runtime:** `Computer_Science/Distributed_Systems/agent_concurrency.py` provides a Ray-like async runtime for parallel agents.

### 🔌 Model Context Protocol (MCP)
*   **BioMCP Server:** `MCP_Servers/BioMCP/bio_mcp_server.py` implements a compliant MCP server exposing bio-tools (`sequence_length`, `reverse_complement`) to LLMs like Claude Desktop.

### 🧪 Clinical & Drug Discovery
*   **Adaptive Clinical Trials:** `Clinical/Clinical_Trials/Adaptive_Trial_Design_Agent/adaptive_trial_sim.py` runs **Bayesian Multi-Arm Multi-Stage (MAMS)** simulations, automatically dropping futile arms.
*   **AlphaFold 3 Wrapper:** `Foundation_Models/AlphaFold3_Agent/af3_wrapper.py` standardizes protein structure prediction calls.
*   **MedPrompt:** `LLM_Research/Prompt_Engineering/medprompt.py` implements Microsoft's SOTA clinical reasoning strategy (Dynamic Few-Shot + Ensemble).
*   **Self-Driving Labs:** `Mathematics/Probability_Statistics/bayesian_optimization.py` enables autonomous experiment selection using Gaussian Processes.

### 📐 Math & CS (The Foundation)
*   **Tensor Operations:** `Mathematics/Linear_Algebra/tensor_operations.py` breaks down the math behind Attention mechanisms.
*   **Graph RAG:** `Computer_Science/Graph_Algorithms/knowledge_graph.py` provides traversal for Drug-Target-Disease interactions.

## 📂 Directory Structure

```text
Skills/
├── Agentic_AI/           # Architectures (ReAct, Plan&Solve, Orchestrators)
├── Clinical/             # MedPrompt, Note Summarization, Adaptive Trials
├── Computer_Science/     # Graph Algo, Distributed Systems (Async)
├── Drug_Discovery/       # ChemCrow, Self-Driving Labs
├── Foundation_Models/    # AlphaFold3 Wrapper, BiomedGPT
├── Genomics/             # Single Cell (Annotation, Pathways), CRISPR
├── LLM_Research/         # RAG, Fine-Tuning, Prompt Engineering
├── Mathematics/          # Bayesian Opt, Linear Algebra, Probability
└── MCP_Servers/          # BioMCP Implementation
```

## 🛠️ Usage Examples

**1. Run the Multi-Agent Orchestrator:**
```bash
python3 Agentic_AI/Multi_Agent_Systems/orchestrator.py
```

**2. Annotate Single-Cell Data:**
```bash
python3 Genomics/Single_Cell/Cell_Type_Annotation/RNA/universal_annotator.py
```

**3. Simulate an Adaptive Clinical Trial:**
```bash
python3 Clinical/Clinical_Trials/Adaptive_Trial_Design_Agent/adaptive_trial_sim.py
```

**4. Start the BioMCP Server:**
```bash
python3 MCP_Servers/BioMCP/bio_mcp_server.py
```

## 📈 Roadmap (2026)
*   [x] **Phase 1:** Core Architectures (Orchestrator, Async Runtime) - *Completed Jan 2026*
*   [x] **Phase 2:** Math Foundations (Bayesian Opt, Graph Theory) - *Completed Jan 2026*
*   [x] **Phase 3:** Single Cell & Clinical Simulators - *Completed Jan 2026*
*   [x] **Phase 4:** Initial MCP Server Integration - *Completed Jan 2026*
*   [ ] **Phase 5:** Deployment to FHIR Servers & Real Lab Integration

---
*Maintained by the Artificial Intelligence Group.*

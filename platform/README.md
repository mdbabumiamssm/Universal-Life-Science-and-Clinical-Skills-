# Universal Biomedical Skills Platform - Prototype

This directory contains the prototype implementation for a **platform-agnostic biomedical skills system** that can deploy to Claude, ChatGPT, Gemini, and other LLMs.

## Directory Structure

```
platform_prototype/
├── schema/
│   └── usdl_schema_v1.json    # JSON Schema for USDL validation
├── examples/
│   └── single_cell_qc.yaml    # Example skill in USDL format
├── adapters/
│   ├── claude_adapter.py      # Convert USDL → Claude formats
│   └── openai_adapter.py      # Convert USDL → OpenAI formats
└── README.md                   # This file
```

## Quick Start

### 1. Define a Skill (USDL Format)

```yaml
skill:
  id: "biomedical.genomics.my_skill"
  version: "1.0.0"
  name: "My Biomedical Skill"
  category: "genomics/single_cell"

  capabilities:
    - name: "analyze_data"
      description: "Analyze the data"
      inputs:
        - name: "dataset"
          type: "string"
          required: true

  prompts:
    system: "You are an expert bioinformatician..."
```

### 2. Convert to Claude Format

```bash
python adapters/claude_adapter.py examples/single_cell_qc.yaml ./output/claude/
```

Outputs:
- `SKILL.md` - Claude Code skill file
- `mcp_server.json` - MCP server configuration
- `claude_tools.json` - Claude API tool_use schema
- `hooks.json` - Claude hooks configuration

### 3. Convert to OpenAI Format

```bash
python adapters/openai_adapter.py examples/single_cell_qc.yaml ./output/openai/
```

Outputs:
- `custom_gpt.json` - GPT Builder configuration
- `assistant.json` - Assistants API configuration
- `functions.json` - Function calling schemas
- `openapi.yaml` - OpenAPI spec for GPT Actions

## USDL Schema Overview

The Universal Skill Definition Language (USDL) provides:

| Section | Purpose |
|---------|---------|
| `id` | Unique identifier (e.g., `biomedical.genomics.skill_name`) |
| `capabilities` | Functions/tools the skill provides |
| `prompts` | System prompts and examples |
| `tools` | Implementation references |
| `dependencies` | Python packages, data files |
| `platform_configs` | Platform-specific overrides |

## Supported Platforms

| Platform | Adapter | Output Formats |
|----------|---------|----------------|
| Claude | `claude_adapter.py` | MCP Server, SKILL.md, tool_use |
| OpenAI/ChatGPT | `openai_adapter.py` | Custom GPT, Assistants, Functions, OpenAPI |
| Gemini | `gemini_adapter.py` | (Coming soon) |

## Development Roadmap

See [STRATEGY_IMPROVEMENT_PLAN.md](STRATEGY_IMPROVEMENT_PLAN.md) for the evolved strategy focusing on **AI-driven Optimization**, **Evaluation Loops**, and the **BioKernel Runtime**.

- [x] USDL Schema v1.0
- [x] Claude Adapter
- [x] OpenAI Adapter
- [ ] **Prompt Optimizer (Meta-Prompter)**
- [ ] **Automated Evaluation Engine**
- [ ] **BioKernel Runtime (Local Server)**
- [ ] Gemini Adapter
- [ ] CLI Tool (`bioskills`)
- [ ] Web Dashboard
- [ ] Skill Registry API

## Full Strategy Document

See [UNIVERSAL_SKILLS_PLATFORM_STRATEGY.md](../UNIVERSAL_SKILLS_PLATFORM_STRATEGY.md) for the complete strategic outline.

---

## 👤 Author
**MD BABU MIA**  
*Artificial Intelligence Group*  
md.babu.mia@mssm.edu

---

*Prototype Version 0.1.0 - December 2025*

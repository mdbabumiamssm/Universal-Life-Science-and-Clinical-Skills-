# Universal Biomedical Skills Platform

**Production-Ready AI Skills for Life Science & Clinical Applications**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Platform: Multi-LLM](https://img.shields.io/badge/Platform-Claude%20%7C%20ChatGPT%20%7C%20Gemini-purple.svg)]()

---

## Repository Structure

```
skills/
├── src/                         # Main source code (GitHub repo)
│   ├── Skills/                  # Production skills
│   │   ├── Clinical/            # Healthcare AI capabilities
│   │   ├── Drug_Discovery/      # Cheminformatics & pharma
│   │   └── Genomics/            # Bioinformatics & sequencing
│   ├── test_demonstration/      # Skill validation suite
│   └── presentation_materials/  # Presentations & tutorials
│
├── platform/                    # Platform prototype
│   ├── adapters/                # USDL to platform converters
│   ├── schema/                  # USDL JSON schemas
│   ├── evaluator/               # Cross-platform evaluation
│   └── optimizer/               # Prompt optimization engine
│
├── tests/                       # Local test suite
│   ├── qc_analysis.py           # QC pipeline tests
│   └── scRNAsedata/             # Sample datasets
│
├── docs/                        # Documentation
│   ├── presentations/           # Slide decks & assets
│   ├── research/                # Academic papers
│   └── strategy/                # Platform strategy docs
│
└── archive/                     # Archived materials (git-ignored)
    ├── external_references/     # Curated external repos
    ├── old_versions/            # Previous implementations
    └── presentation_drafts/     # Draft presentations
```

---

## Quick Start

### 1. Clone and Setup

```bash
cd /home/drdx/ARTIFICIALINTELLIGENCEGROUP/skills

# The main source is in src/ (tracks GitHub repo)
cd src

# Install dependencies for testing
pip install -r test_demonstration/requirements.txt
```

### 2. Run a Skill

```bash
# Single-Cell RNA-seq QC
python test_demonstration/qc_analysis.py sample.h5ad --output-dir results/

# Or use as Python library
python -c "from Skills.Genomics.Single_Cell_RNA_QC.qc_core import calculate_qc_metrics; print('Ready!')"
```

---

## Key Components

### Source (`src/`)

The main codebase tracking the GitHub repository:
- **Skills/**: Production-ready AI skills with prompts and tools
- **test_demonstration/**: Validation suite with sample data
- **presentation_materials/**: Documentation and tutorials

### Platform (`platform/`)

Universal Skill Definition Language (USDL) infrastructure:
- **adapters/**: Convert USDL to Claude MCP, OpenAI Actions, Gemini Extensions
- **schema/**: JSON schemas for skill definitions
- **evaluator/**: Cross-platform accuracy testing

### Tests (`tests/`)

Local testing environment with:
- Sample scRNA-seq datasets (bone marrow samples)
- QC validation scripts
- Integration tests

### Docs (`docs/`)

Project documentation:
- **presentations/**: Current slide decks and visual assets
- **research/**: Academic paper drafts
- **strategy/**: Platform roadmap and architecture

### Archive (`archive/`)

Git-ignored folder containing:
- **external_references/**: 27+ curated repos (3GB)
- **old_versions/**: Previous skill implementations
- **presentation_drafts/**: Earlier presentation versions

---

## Development Workflow

### Making Changes

1. Edit files in `src/` for production changes
2. Test locally using `tests/` scripts
3. Commit and push from `src/` directory

```bash
cd src
git add .
git commit -m "Your changes"
git push origin main
```

### Adding New Skills

1. Create skill directory: `src/Skills/<Category>/<SkillName>/`
2. Add required files: `README.md`, `prompt.md`, implementation files
3. Add tests in `tests/` directory
4. Update platform adapters if needed

---

## GitHub Repository

**Remote:** https://github.com/mdbabumiamssm/LLMs-Universal-Life-Science-and-Clinical-Skills-

The `src/` directory is synchronized with the GitHub repository.

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

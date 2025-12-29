# GitHub Push Guide

**Author:** MD BABU MIA
**Email:** md.babu.mia@mssm.edu
**Affiliation:** Icahn School of Medicine at Mount Sinai

---

## Repository Structure

This folder contains two Git repositories:

1. **Main Folder** (`/skills/`) - Local workspace repository
2. **src/** - GitHub submodule (already pushed to GitHub)

---

## Documents Ready for GitHub

### Currently on GitHub (via `src/` submodule)

The `src/` folder is already connected to:
**https://github.com/mdbabumiamssm/LLMs-Universal-Life-Science-and-Clinical-Skills-**

| Document | Path | Description |
|----------|------|-------------|
| Main README | `src/README.md` | Project overview & quick start |
| Single-Cell QC | `src/Skills/Genomics/Single_Cell_RNA_QC/README.md` | scRNA-seq quality control skill |
| CRISPR Design | `src/Skills/Genomics/CRISPR_Design_Agent/README.md` | Guide RNA design agent |
| Clinical Notes | `src/Skills/Clinical/Clinical_Note_Summarization/README.md` | SOAP note structuring |
| Trial Eligibility | `src/Skills/Clinical/Trial_Eligibility_Agent/README.md` | Patient-trial matching |
| AgentD | `src/Skills/Drug_Discovery/AgentD_Drug_Discovery/README.md` | Drug discovery agent |
| Chemical Properties | `src/Skills/Drug_Discovery/Chemical_Property_Lookup/README.md` | RDKit molecular tools |
| Test Suite | `src/test_demonstration/README.md` | Validation documentation |

### Local Documents (Can Be Pushed Separately)

| Document | Path | Description |
|----------|------|-------------|
| Workspace README | `README.md` | Local folder structure guide |
| Platform README | `platform/README.md` | USDL platform documentation |
| Strategy Plan | `platform/STRATEGY_IMPROVEMENT_PLAN.md` | Platform roadmap |
| Research Paper | `docs/research/RESEARCH_PAPER_DRAFT.md` | Academic paper draft |
| Platform Strategy | `docs/strategy/UNIVERSAL_SKILLS_PLATFORM_STRATEGY.md` | Architecture overview |
| Presentation | `docs/presentations/Universal_Biomedical_Skills_Presentation_v2.0.0.pdf` | Slide deck |
| Video Scripts | `docs/presentations/SCRIPTS_VIDEO_*.md` | Tutorial scripts |

---

## How to Push Changes

### Option 1: Push to Existing GitHub Repo (Recommended)

Changes to skill documentation go through the `src/` submodule:

```bash
# Navigate to the submodule
cd src/

# Check what changed
git status

# Add and commit changes
git add .
git commit -m "Your commit message"

# Push to GitHub
git push origin main
```

### Option 2: Create New GitHub Repo for Full Workspace

If you want to push the entire workspace (including platform, docs, tests):

```bash
# From the main skills folder
cd /home/drdx/ARTIFICIALINTELLIGENCEGROUP/skills/

# Create new repo on GitHub first, then:
git remote add origin https://github.com/YOUR_USERNAME/YOUR_NEW_REPO.git
git push -u origin main
```

---

## Author Attribution

All documents contain the following author information:

```
Author: MD BABU MIA
Affiliation: Artificial Intelligence Group
Institution: Icahn School of Medicine at Mount Sinai
Email: md.babu.mia@mssm.edu
```

### Documents with Author Info

| File | Contains Author Section |
|------|------------------------|
| `README.md` | ✅ Yes |
| `src/README.md` | ✅ Yes |
| All skill READMEs in `src/Skills/` | ✅ Yes |
| `platform/README.md` | ✅ Yes |
| `tests/README.md` | ✅ Yes |

---

## Quick Reference Commands

```bash
# Check status of main repo
git status

# Check status of GitHub submodule
cd src && git status

# Push changes to GitHub (from src/)
git add . && git commit -m "message" && git push

# Pull latest from GitHub
cd src && git pull origin main

# View commit history
git log --oneline -10
```

---

## What NOT to Push

The following are git-ignored and should NOT be pushed:

| Item | Reason |
|------|--------|
| `archive/` | Contains 3GB of external references |
| `*.h5ad`, `*.h5` | Large data files |
| `venv/`, `__pycache__/` | Python environment files |
| `token.txt` | Sensitive credentials |
| `.claude/` | Local Claude Code settings |

---

## Current GitHub Repository

**URL:** https://github.com/mdbabumiamssm/LLMs-Universal-Life-Science-and-Clinical-Skills-

**Last Updated:** December 29, 2025

**Contents:**
- 6 Production Skills (Genomics, Clinical, Drug Discovery)
- Comprehensive documentation with methodology references
- Test demonstration suite
- Presentation materials

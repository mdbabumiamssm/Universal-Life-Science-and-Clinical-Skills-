# LLMs Universal Life Science & Clinical Skills

**The Open-Source Operating System for Biomedical AI Agents**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Status: Production](https://img.shields.io/badge/Status-Production-brightgreen.svg)]()
[![Platform: Multi-LLM](https://img.shields.io/badge/Platform-Claude%20%7C%20ChatGPT%20%7C%20Gemini-purple.svg)]()

---

## Mission

We build **production-ready, platform-agnostic biomedical AI skills** that empower researchers, clinicians, and developers to deploy advanced AI capabilities across any LLM interface. Whether you use **Claude**, **ChatGPT**, **Gemini**, or custom open-source models, our standardized skills deliver reproducible, validated results for real-world biomedical workflows.

## Why This Repository?

| Challenge | Our Solution |
|-----------|--------------|
| Biomedical AI tools are fragmented across platforms | **Universal Skill Definition Language (USDL)** compiles once, deploys everywhere |
| Most AI prompts lack scientific validation | Every skill follows **peer-reviewed methodologies** with citations |
| Integration is complex and time-consuming | **Drop-in Python modules** work with LangChain, AutoGen, Semantic Kernel |
| Results are often non-reproducible | **Statistical rigor** (MAD-based filtering, validated thresholds) ensures consistency |

---

## Skill Categories

### Genomics & Bioinformatics

| Skill | Description | Status |
|-------|-------------|--------|
| [Single-Cell RNA-seq QC](Skills/Genomics/Single_Cell_RNA_QC/) | Production-grade quality control with MAD-based outlier detection | Production |
| [CRISPR Guide Design](Skills/Genomics/CRISPR_Design_Agent/) | Automated sgRNA design with off-target prediction | Production |
| Variant Interpretation | Clinical significance assessment (ClinVar, ACMG) | Planned |
| Spatial Transcriptomics | Tissue-aware spatial analysis | Planned |

### Clinical AI

| Skill | Description | Status |
|-------|-------------|--------|
| [Clinical Note Summarization](Skills/Clinical/Clinical_Note_Summarization/) | SOAP-format structuring from unstructured dictation | Production |
| [Trial Eligibility Screening](Skills/Clinical/Trial_Eligibility_Agent/) | Automated patient-trial matching against I/E criteria | Production |
| Diagnostic Decision Support | Differential diagnosis with evidence synthesis | Planned |
| Medical Coding (ICD-10/SNOMED) | Automated code extraction from clinical text | Planned |

### Drug Discovery & Cheminformatics

| Skill | Description | Status |
|-------|-------------|--------|
| [AgentD Drug Discovery](Skills/Drug_Discovery/AgentD_Drug_Discovery/) | Literature mining, molecule generation, ADMET prediction | Production |
| [Chemical Property Lookup](Skills/Drug_Discovery/Chemical_Property_Lookup/) | RDKit-powered molecular property calculation | Production |
| Target Validation | Knowledge graph traversal for target identification | Planned |

---

## Quick Start

### For Developers

Integrate any skill into your LLM agent pipeline:

```python
# Example: Single-Cell QC in a LangChain agent
from langchain.tools import tool
import sys
sys.path.append("Skills/Genomics/Single_Cell_RNA_QC")
from qc_core import calculate_qc_metrics, filter_cells

@tool
def single_cell_qc(file_path: str, mad_threshold: float = 5.0) -> str:
    """Performs automated quality control on scRNA-seq data using MAD-based filtering."""
    import anndata as ad
    adata = ad.read_h5ad(file_path)
    calculate_qc_metrics(adata, inplace=True)
    adata_filtered = filter_cells(adata, mad_threshold=mad_threshold)
    output_path = file_path.replace(".h5ad", "_qc_filtered.h5ad")
    adata_filtered.write(output_path)
    return f"QC complete. Filtered {adata.n_obs - adata_filtered.n_obs} cells. Output: {output_path}"
```

### For Prompt Engineers

Access validated prompts directly:

```python
# Load the Clinical Note Summarization prompt
with open("Skills/Clinical/Clinical_Note_Summarization/prompt.md") as f:
    system_prompt = f.read()

# Use with any LLM API
response = client.chat.completions.create(
    model="gpt-4",
    messages=[
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": clinical_note_text}
    ]
)
```

### For Researchers

Run the demonstration pipeline:

```bash
cd test_demonstration
pip install -r requirements.txt
python qc_analysis.py ../sample_data/sample.h5ad --output-dir ./results
```

---

## Platform Deployment

Our **Universal Skill Definition Language (USDL)** enables single-source deployment:

| Platform | Output Format | Integration |
|----------|---------------|-------------|
| **Claude** | MCP Servers, SKILL.md files | Model Context Protocol |
| **ChatGPT** | Custom GPT Actions, OpenAPI schemas | Assistants API, GPT Builder |
| **Gemini** | Function Declarations | Vertex AI Extensions |
| **Open Source** | Python modules | LangChain, AutoGen, CrewAI |

---

## Recent Updates (December 2025)

### Comprehensive Expansion

We cataloged **45+ biomedical AI agents** across **18 specialized categories**:

- **High-Priority Additions**: Biomni (150 tools), STAgent (spatial transcriptomics), BioMaster (multi-omic pipelines), CellAgent (single-cell automation)
- **New MCP Servers**: BioMCP for PubMed/PMC, ClinicalTrials.gov, genomic variant databases
- **RAG Systems**: BiomedRAG, MEGA-RAG with 40%+ hallucination reduction

See [COMPREHENSIVE_SKILLS_UPDATE_DEC_2025.md](COMPREHENSIVE_SKILLS_UPDATE_DEC_2025.md) for complete details.

---

## Repository Structure

```
.
├── Skills/
│   ├── Clinical/                    # Healthcare AI capabilities
│   ├── Drug_Discovery/              # Cheminformatics & pharma
│   └── Genomics/                    # Bioinformatics & sequencing
├── test_demonstration/              # Validation suite with sample data
├── presentation_materials/          # Documentation & tutorials
├── COMPREHENSIVE_SKILLS_UPDATE_DEC_2025.md
└── NEW_SKILLS_DISCOVERY_DEC_2025.md
```

---

## Contributing

We welcome contributions from domain experts. Priority areas:

- **Protein Structure Prediction** (AlphaFold integration)
- **Medical Imaging** (pathology, radiology)
- **Pharmacovigilance** (adverse event detection)
- **Clinical Coding** (ICD-10, SNOMED-CT automation)

All contributions must include:
1. Validation metrics on benchmark datasets
2. Peer-reviewed methodology citations
3. Example prompts and test cases

---

## Citation

If you use these skills in your research, please cite:

```bibtex
@software{universal_biomedical_skills,
  author = {Mia, MD Babu},
  title = {LLMs Universal Life Science and Clinical Skills},
  year = {2025},
  url = {https://github.com/mdbabumiamssm/LLMs-Universal-Life-Science-and-Clinical-Skills-}
}
```

---

## Author & Maintainer

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

## License

[MIT License](LICENSE) - Free for academic and commercial use.

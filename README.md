# Universal Biomedical Skills for LLMs

A comprehensive, open-source collection of "Skills" (Prompts, Tools, and Agents) designed for Biomedical, Clinical, Genomics, and Life Science applications. This repository aims to provide standardized, model-agnostic building blocks for the next generation of AI in Healthcare.

## 📂 Repository Structure

- **Skills/**: The core library of capabilities.
    - **Clinical/**: Skills for doctor-patient interactions, EHR summarization, and diagnostics support.
        - `Clinical_Note_Summarization`: Turn unstructured notes into SOAP format.
        - `Trial_Eligibility_Agent`: AI screener for matching patients to clinical trials.
        - `Oncology`: (Planned) Precision oncology tools.
    - **Genomics/**: Tools for analyzing DNA/RNA sequencing data.
        - `Single_Cell_RNA_QC`: Automated Quality Control for scRNA-seq (scanpy/scverse).
        - `CRISPR_Design_Agent`: Automated design of sgRNAs and off-target analysis.
        - `Spatial_Transcriptomics`: (Planned) Agents for ST analysis.
    - **Drug_Discovery/**: Cheminformatics and molecule analysis tools.
        - `Chemical_Property_Lookup`: Calculate molecular properties using RDKit.
        - `AgentD_Drug_Discovery`: AI-driven agent for literature mining and molecule generation.
    - **Research_Tools/**: General academic research aids (Literature search, etc.).
    - **MCP_Servers/**: (Planned) Model Context Protocol servers for biomedical data.

- **test_demonstration/**: Ready-to-run demos and test scripts.
    - `qc_analysis.py`: Automated Single-Cell RNA-seq QC pipeline.
    - `generate_dummy_data.py`: Create synthetic scRNA-seq datasets for testing.
    - `run_test.sh`: One-click script to run the full QC demo.

## 🆕 New Discoveries (Dec 2025)

### Comprehensive Update (Dec 28, 2025)
We have just released a massive update to our skills database! Check out **[COMPREHENSIVE_SKILLS_UPDATE_DEC_2025.md](COMPREHENSIVE_SKILLS_UPDATE_DEC_2025.md)** for an extensive report covering:

- **45+ New Biomedical AI Agents**: Including Biomni, STAgent, BioMaster, and CellAgent.
- **18 New Categories**: Ranging from Spatial Transcriptomics to Precision Oncology.
- **MCP Servers**: Emerging standard for connecting AI agents to biomedical data (BioMCP).
- **Clinical AI**: New tools for trial matching (TrialGPT) and radiology (RadGPT).
- **Genomics**: Advanced variant interpretation (DYNA, AlphaMissense) and CRISPR design tools.
- **Curated Resources**: 9 new "Awesome Lists" for tracking the exploding field of bio-AI.

### Earlier Discoveries
See [NEW_SKILLS_DISCOVERY_DEC_2025.md](NEW_SKILLS_DISCOVERY_DEC_2025.md) for our initial December report on LLM agents and tools.

## 🚀 Getting Started

### Running the New Single-Cell QC Demo
We've added a fully functional Single-Cell RNA-seq Quality Control demo.

```bash
cd test_demonstration
# Install requirements
pip install -r requirements.txt
# Run the demo (generates data and performs QC)
./run_test.sh
```
Check the `qc_results` folder for generated plots and filtered datasets!

### For Developers
You can import the Python scripts in `Skills/` directly into your LangChain, Semantic Kernel, or AutoGen workflows.

### For Prompt Engineers
Check the `prompt.md` files in the `Clinical` section for high-quality, tested medical prompts.

## 🤝 Contributing

We welcome contributions! If you have a prompt or tool for:
- Protein folding
- Clinical trial matching
- Medical coding (ICD-10)
- CRISPR guide design

Please submit a Pull Request.

## 📜 License
[MIT License](LICENSE)

## 👤 Author
**MD BABU MIA**
md.babu.mia@mssm.edu
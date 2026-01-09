# Single Cell Sequencing: Comprehensive Analysis Suite

**Last Updated:** 2026
**Status:** Active
**Focus:** Cell Type Annotation, Pathway Analysis, Cell-Cell Communication, Trajectory Inference.

## 📂 Directory Overview

### 1. [Cell Type Annotation](./Cell_Type_Annotation/)
A robust collection of tools for identifying cell identity across modalities.
- **RNA:** `universal_annotator.py` wraps Marker-based, LLM-based, and CellTypist methods.
- **ATAC:** Tools for chromatin accessibility annotation (ArchR, Signac wrappers).
- **MultiModal:** Integration tools for CITE-seq and Multiome.
- **Tool Database:** A curated list of [SOTA Tools](./Tool_Database.md).

### 2. [Pathway Analysis](./Pathway_Analysis/)
Functional interpretation of gene lists.
- **`sc_pathway_scorer.py`**: Implements single-cell scoring logic similar to AUCell and PAGODA.
- **Supported Databases:** MSigDB, Hallmark, Reactome.

### 3. [Cell-Cell Communication](./Cell_Cell_Communication/)
Deciphering the language of cells.
- **`interaction_inference.py`**: Ligand-Receptor analysis framework inspired by CellChat/NicheNet.
- **Focus:** Paracrine and Autocrine signaling networks.

### 4. [Trajectory Inference](./Trajectory_Inference/)
Modeling dynamic cellular processes.
- **Tools:** Monocle3, RNA Velocity (scVelo), CellRank.

## 🚀 Quick Start

**Annotate Cells (RNA):**
```bash
python3 Cell_Type_Annotation/RNA/universal_annotator.py
```

**Infer Cell-Cell Interactions:**
```bash
python3 Cell_Cell_Communication/interaction_inference.py
```

**Calculate Pathway Scores:**
```bash
python3 Pathway_Analysis/sc_pathway_scorer.py
```

## 📚 Recommended Reading (2026)
1. **MultiKano:** "Deep learning for multimodal integration" (2025).
2. **scPS:** "Benchmarking single-cell pathway scores" (2025).
3. **MedPrompt:** "Adapting LLMs for biomedical annotation" (Microsoft, 2024).

# CRISPR-Cell Painting Integration Agent

**ID:** `biomedical.research_tools.crispr_cell_painting`
**Version:** 1.0.0
**Status:** Beta
**Category:** Research Tools / Functional Genomics

---

## Overview

The **CRISPR-Cell Painting Integration Agent** combines pooled CRISPR screening with morphological profiling for hypothesis-free reverse genetic screening. A 2025 Nature Communications paper demonstrated this optical pooled screening (OPS) platform enables de novo inference of gene function through self-supervised deep learning without target-specific biomarkers.

---

## Key Capabilities

### 1. Screening Modes

| Mode | Method | Scale | Application |
|------|--------|-------|-------------|
| **OPS-Cell Painting** | in situ sequencing + imaging | Genome-wide | Gene function |
| **Arrayed CRISPR** | One guide per well | Focused | Validation |
| **CRISPRi/a** | Knockdown/activation | Variable | Directional effects |

### 2. Gene Network Inference

- **Morphological clustering:** Group genes by phenotype
- **Pathway inference:** Identify gene-pathway associations
- **Epistasis analysis:** Gene-gene interactions

### 3. Integration Features

- **Barcode demixing:** Assign phenotypes to guides
- **Batch correction:** Cross-plate normalization
- **Hit calling:** Statistical significance

---

## Usage

### Example Prompt

```text
Analyze a pooled CRISPR-Cell Painting screen targeting autophagy genes.

Input:
- OPS images with in situ barcode readout
- 200 autophagy-related gene targets
- 4 guides per gene
- Cell Painting readout

Identify genes with similar morphological phenotypes and
infer functional relationships.
```

### Expected Output

```
## CRISPR-Cell Painting Screen Analysis

### Screen Summary
| Metric | Value |
|--------|-------|
| Genes targeted | 200 |
| Total guides | 800 |
| Cells analyzed | 1.2M |
| Phenotypic hits | 89 genes |

### Phenotypic Clusters

| Cluster | Genes | Key Phenotype | Pathway |
|---------|-------|---------------|---------|
| **1** | ATG5, ATG7, ATG12, BECN1 | ↑LC3 puncta, ↑ER | Core autophagy |
| **2** | ULK1, ULK2, ATG13 | ↑Mito mass | Initiation |
| **3** | ATG9A, ATG9B | ↑Lipid droplets | Membrane source |
| **4** | LAMP1, LAMP2 | ↑Lysosome size | Degradation |

### Gene Function Predictions

| Gene | Known Function | Predicted Function | Confidence |
|------|---------------|-------------------|------------|
| FAM134B | ER-phagy receptor | ER-phagy confirmed | High |
| TEX264 | Unknown | ER-phagy (novel) | High |
| CCPG1 | Cell cycle? | ER-phagy (novel) | Medium |

### Novel Findings
**TEX264:** Previously uncharacterized, clusters with known ER-phagy receptors
- Phenotype matches FAM134B knockdown
- Suggested function: ER-phagy receptor
- **Validation recommended**
```

---

## References

- **Ramezani et al. (2025):** "A pooled Cell Painting CRISPR screening platform enables de novo inference of gene function." *Nature Communications*
- **Feldman et al. (2019):** "Optical pooled screens in human cells." *Cell*

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

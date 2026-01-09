# Single-Cell Analysis Tools Database (2026 Edition)

A curated registry of state-of-the-art tools for single-cell analysis, updated for the 2026 landscape.

## 1. Cell Type Annotation

### RNA-seq (scRNA-seq)
| Tool | Description | Key Features | Language | Link |
|------|-------------|--------------|----------|------|
| **CellIdentifierDx** | Deep learning-based cell type identification. | Robust to batch effects, handles large datasets. | Python | [GitHub](https://github.com/mdbabumiamssm/cellidentifierdx) |
| **CellTypist** | Automated cell type annotation using massive pre-trained models. | Immune cell specialization, very fast. | Python | [GitHub](https://github.com/Teichlab/celltypist) |
| **SingleR** | Reference-based annotation using correlation profiles. | Gold standard for reference-based mapping. | R | [Bioconductor](https://bioconductor.org/packages/release/bioc/html/SingleR.html) |
| **mLLMCelltype** | LLM-based annotation (GPT-4/Claude/Gemini). | Uses marker gene lists + LLM reasoning. | Web/Py | [Website](https://mllmcelltype.com) |
| **Seurat (Azimuth)** | Reference mapping onto CZI Human Cell Atlas. | Fast, integrated with Seurat ecosystem. | R | [SatijaLab](https://satijalab.org/azimuth/) |
| **scType** | Marker-gene based annotation. | Unsupervised, highly accurate for rare types. | R | [Nature Methods](https://www.nature.com/articles/s41592-022-01517-5) |

### ATAC-seq (scATAC-seq)
| Tool | Description | Key Features | Language | Link |
|------|-------------|--------------|----------|------|
| **ArchR** | Comprehensive suite for scATAC-seq. | Peak calling, motif enrichment, integration. | R | [ArchRProject](https://www.archrproject.com/) |
| **Signac** | Seurat extension for chromatin data. | Seamless RNA-ATAC integration. | R | [Signac](https://stuartlab.org/signac/) |
| **EpiScanpy** | Scanpy for epigenomics. | Python-native, integrates with scverse. | Python | [GitHub](https://github.com/colomemaria/epiScanpy) |

### Multi-Modal (RNA + ATAC + Protein)
| Tool | Description | Key Features | Language | Link |
|------|-------------|--------------|----------|------|
| **MultiKano** | Neural network for RNA+ATAC integration. | Uses Kolmogorov-Arnold networks. | Python | [Research](https://academic.oup.com/bib) |
| **WNN (Seurat)** | Weighted Nearest Neighbors. | Defines cell state based on weighted combination of modalities. | R | [Paper](https://www.cell.com/cell/fulltext/S0092-8674(21)00583-3) |
| **TotalVI** | Variational Inference for RNA + Protein (CITE-seq). | Probabilistic modeling of background noise. | Python | [scvi-tools](https://scvi-tools.org/) |

---

## 2. Pathway & Functional Analysis

| Tool | Description | Application | Link |
|------|-------------|-------------|------|
| **scPS** | Single-cell Pathway Score. | Benchmark-winner (2025) for pathway activity. | [Paper](https://academic.oup.com/bib) |
| **AUCell** | Area Under Curve for gene sets. | Identifies cells with active gene sets. | [Bioconductor](https://bioconductor.org/packages/release/bioc/html/AUCell.html) |
| **Pagoda2** | Pathway Gene Set Overdispersion Analysis. | Finds heterogeneous pathways. | [GitHub](https://github.com/kharchenkolab/pagoda2) |
| **PROGENy** | Pathway Response signatures. | Infers pathway activity from downstream genes. | [Bioconductor](https://bioconductor.org/packages/PROGENy) |

---

## 3. Cell-Cell Communication

| Tool | Description | Key Feature (2026) | Link |
|------|-------------|--------------------|------|
| **CellChat v2** | Ligand-Receptor inference. | Enhanced database, spatial constraints. | [GitHub](https://github.com/jinworks/CellChat) |
| **LIANA+** | Consensus framework. | Combines predictions from CellPhoneDB, NicheNet, etc. | [GitHub](https://github.com/saezlab/liana) |
| **NicheNet** | Target gene prediction. | Predicts *downstream* effects of signaling. | [GitHub](https://github.com/saeyslab/nichenetchr) |
| **GITIII** | Graph Transformer. | Self-supervised learning for spatial interactions. | [Nature MI](https://www.nature.com/nmi/) |
| **scSeqCommDiff** | Differential communication. | Compares communication between conditions (Disease vs Control). | [BioRxiv](https://www.biorxiv.org/) |

---

## 4. Trajectory Inference

| Tool | Description | Type | Link |
|------|-------------|------|------|
| **Monocle 3** | Single-cell pseudotime. | Tree-based trajectories. | [Cole Trapnell Lab](https://cole-trapnell-lab.github.io/monocle3/) |
| **RNA Velocity (scVelo)** | Future state prediction. | Uses spliced/unspliced ratios. | [scVelo](https://scvelo.org/) |
| **CellRank** | Markov transition matrices. | Combines velocity with similarity kernels. | [CellRank](https://cellrank.org/) |

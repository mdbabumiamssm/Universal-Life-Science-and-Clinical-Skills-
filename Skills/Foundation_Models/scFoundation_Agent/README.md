# scFoundation Agent

**ID:** `biomedical.foundation_models.scfoundation`
**Version:** 1.0.0
**Status:** Production
**Category:** Foundation Models / Single-Cell Genomics AI
**Released:** December 2025

---

## Overview

The **scFoundation Agent** provides unified access to state-of-the-art single-cell foundation models including scFoundation, scGPT, Geneformer, and scPRINT. These transformer-based models are pretrained on tens to hundreds of millions of cells, enabling powerful zero-shot and fine-tuned analysis of single-cell data.

### Model Comparison

| Model | Training Cells | Architecture | Key Strength |
|-------|---------------|--------------|--------------|
| **scFoundation** | 100M+ | Transformer | Scale, zero-shot |
| **scGPT** | 33M | GPT-style | Generative tasks |
| **Geneformer** | 95M | BERT-style | Gene regulation |
| **scPRINT** | 50M | Encoder | Multi-task |
| **GenePT** | ChatGPT-based | Embedding | Literature knowledge |

### 2025 Critical Evaluation

Recent benchmarking (Genome Biology, April 2025) revealed important insights:
- Zero-shot performance varies significantly by task
- HVG + traditional methods (Harmony, scVI) often competitive
- Fine-tuning dramatically improves foundation model performance
- Task-specific evaluation essential

---

## Key Capabilities

### 1. Core Tasks

| Task | Description | Best Model |
|------|-------------|------------|
| **Cell Type Annotation** | Automatic labeling | Geneformer (fine-tuned) |
| **Batch Integration** | Remove technical artifacts | scGPT |
| **Gene Network Inference** | Regulatory relationships | Geneformer |
| **Perturbation Prediction** | Treatment response | scGPT, scFoundation |
| **Imputation** | Missing value estimation | scFoundation |
| **Multi-omics Integration** | RNA + ATAC + protein | scPRINT |

### 2. Use Modes

| Mode | Description | When to Use |
|------|-------------|-------------|
| **Zero-shot** | No additional training | Quick exploration |
| **Linear probe** | Train small head | Limited labels |
| **Fine-tuning** | Update model weights | Maximum accuracy |
| **Embedding** | Extract representations | Downstream ML |

### 3. Supported Data Types

```yaml
Input Formats:
  - AnnData (.h5ad)
  - Seurat (.rds via conversion)
  - 10x CellRanger outputs
  - Loom files

Data Types:
  - scRNA-seq (droplet, plate-based)
  - snRNA-seq
  - CITE-seq (RNA + protein)
  - Multiome (RNA + ATAC)
  - Spatial transcriptomics
```

---

## Usage

### Example Prompt

```text
I have a single-cell RNA-seq dataset from a COVID-19 patient study:
- 150,000 PBMCs from 20 patients (10 severe, 10 mild)
- 10x Genomics, already preprocessed in AnnData format
- Goal: Identify disease-associated cell states and gene programs

Use foundation models to:
1. Annotate cell types without manual labeling
2. Identify COVID-severity-associated cell states
3. Predict gene regulatory networks in severe disease
4. Find potential therapeutic targets
```

### Expected Output

```markdown
## scFoundation Analysis: COVID-19 PBMC Study

### Data Summary
| Metric | Value |
|--------|-------|
| Total cells | 150,000 |
| Patients | 20 (10 severe, 10 mild) |
| Genes | 18,456 |
| Platform | 10x Genomics |

### 1. Zero-Shot Cell Type Annotation (Geneformer)

| Cell Type | Count | % Total | Confidence |
|-----------|-------|---------|------------|
| CD4+ T cells | 42,300 | 28.2% | 0.94 |
| CD8+ T cells | 28,500 | 19.0% | 0.93 |
| CD14+ Monocytes | 25,200 | 16.8% | 0.96 |
| B cells | 18,900 | 12.6% | 0.95 |
| NK cells | 15,300 | 10.2% | 0.92 |
| CD16+ Monocytes | 10,800 | 7.2% | 0.91 |
| Plasmablasts | 4,500 | 3.0% | 0.89 |
| pDCs | 2,700 | 1.8% | 0.88 |
| Other | 1,800 | 1.2% | 0.75 |

### 2. Disease-Associated Cell States (scGPT Embeddings)

```
Severe COVID-Associated Populations:
─────────────────────────────────────

1. Inflammatory Monocytes (S100A+ CD14+)
   ├── Prevalence: 8.2% severe vs 2.1% mild (p<0.001)
   ├── Markers: S100A8, S100A9, S100A12, IL1B
   ├── Pathway: NF-kB activation, cytokine storm
   └── Clinical: Correlates with CRP levels (r=0.78)

2. Exhausted CD8+ T cells
   ├── Prevalence: 12.4% severe vs 4.3% mild (p<0.001)
   ├── Markers: PDCD1, LAG3, HAVCR2, TOX
   ├── Pathway: T cell exhaustion program
   └── Clinical: Inverse correlation with viral clearance

3. Proliferating Plasmablasts
   ├── Prevalence: 5.1% severe vs 1.2% mild (p<0.001)
   ├── Markers: MKI67, IGHG1, XBP1, PRDM1
   ├── Pathway: Antibody secretion
   └── Clinical: Associated with cytokine storm

4. Interferon-Stimulated NK cells
   ├── Prevalence: 3.8% severe vs 6.2% mild (p<0.01)
   ├── Markers: ISG15, MX1, IFI44L, reduced in severe
   ├── Pathway: Type I IFN response (deficient)
   └── Clinical: Protective; correlates with recovery
```

### 3. Gene Regulatory Network (Geneformer Attention)

**Key Regulatory Hubs in Severe COVID:**

| Transcription Factor | Targets | Disease Role |
|---------------------|---------|--------------|
| **STAT1** | 245 genes | IFN response master regulator |
| **NF-κB (RELA)** | 312 genes | Inflammatory cascade |
| **IRF7** | 189 genes | Antiviral response |
| **CEBPB** | 156 genes | Myeloid inflammation |
| **BATF** | 134 genes | T cell exhaustion |

**Novel Regulatory Axis Discovered:**
```
S100A8/A9 → TLR4 → NF-κB → IL-6/IL-1β → Cytokine Storm
         ↓
    STAT1 suppression → Impaired IFN response
```

### 4. Therapeutic Target Predictions (scFoundation)

| Target | Cell Type | Rationale | Drugability |
|--------|-----------|-----------|-------------|
| **JAK1/2** | Monocytes | Block cytokine signaling | Approved (Baricitinib) |
| **S100A8/A9** | Monocytes | Reduce inflammation | In development |
| **PD-1** | CD8+ T | Reverse exhaustion | Approved (Nivolumab) |
| **IL-6R** | Systemic | Block IL-6 pathway | Approved (Tocilizumab) |
| **BTK** | B cells | Reduce inflammation | Approved (Ibrutinib) |

### Validation Summary

| Finding | External Validation | Confidence |
|---------|-------------------|------------|
| S100A+ monocytes | COvid-19 Cell Atlas | High |
| T cell exhaustion | Multiple studies | High |
| IFN deficiency | Bastard et al. 2020 | High |
| JAK inhibitor benefit | RECOVERY trial | Confirmed |

### Recommended Follow-up

1. **Trajectory analysis** of monocyte activation
2. **Pseudotime** analysis of T cell exhaustion
3. **Cell-cell communication** changes in severe disease
4. **In vitro validation** of S100A8/A9 blockade
```

---

## LLM Agent Integration

### Python Tool Implementation

```python
from typing import Optional, Dict, Any, List, Literal
import scanpy as sc
import numpy as np

def scfoundation_tool(
    adata_path: str,
    task: Literal["annotate", "integrate", "perturb", "network", "embed"],
    model: str = "geneformer",
    reference_path: Optional[str] = None,
    perturbation: Optional[str] = None,
    fine_tune: bool = False,
    n_epochs: int = 10,
    device: str = "cuda"
) -> Dict[str, Any]:
    """
    Single-cell foundation model analysis tool.

    Args:
        adata_path: Path to AnnData file (.h5ad)
        task: Analysis task type
        model: 'geneformer', 'scgpt', 'scfoundation', 'scprint', 'genept'
        reference_path: Reference dataset for annotation/integration
        perturbation: Gene/drug for perturbation prediction
        fine_tune: Whether to fine-tune on user data
        n_epochs: Fine-tuning epochs
        device: Compute device

    Returns:
        Task-specific results dictionary
    """
    # Load data
    adata = sc.read_h5ad(adata_path)

    # Initialize model
    if model == "geneformer":
        from geneformer import Geneformer
        fm = Geneformer.from_pretrained("ctheodoris/Geneformer")
    elif model == "scgpt":
        from scgpt import scGPT
        fm = scGPT.from_pretrained("scgpt/scgpt-human")
    elif model == "scfoundation":
        from scfoundation import scFoundation
        fm = scFoundation.from_pretrained("biomap-research/scfoundation-100m")
    elif model == "scprint":
        from scprint import scPRINT
        fm = scPRINT.from_pretrained("jkobject/scprint-50m")
    elif model == "genept":
        from genept import GenePT
        fm = GenePT.from_pretrained()

    fm.to(device)

    results = {"model": model, "task": task}

    if task == "annotate":
        # Cell type annotation
        if fine_tune and reference_path:
            ref_adata = sc.read_h5ad(reference_path)
            fm.fine_tune(ref_adata, label_key="cell_type", n_epochs=n_epochs)

        predictions = fm.predict_cell_types(adata)
        adata.obs['predicted_celltype'] = predictions['labels']
        adata.obs['prediction_confidence'] = predictions['confidence']

        results['predictions'] = predictions
        results['annotation_summary'] = adata.obs['predicted_celltype'].value_counts().to_dict()

    elif task == "integrate":
        # Batch integration
        embeddings = fm.encode(adata, return_embeddings=True)
        adata.obsm['X_foundation'] = embeddings

        # Compute integrated UMAP
        sc.pp.neighbors(adata, use_rep='X_foundation')
        sc.tl.umap(adata)

        # Compute integration metrics
        from scib import metrics
        results['batch_correction'] = {
            'ASW_batch': metrics.silhouette_batch(adata, 'batch', 'X_foundation'),
            'kBET': metrics.kBET(adata, 'batch', 'X_foundation'),
            'graph_connectivity': metrics.graph_connectivity(adata, 'cell_type')
        }

    elif task == "perturb":
        # Perturbation prediction
        if perturbation:
            predictions = fm.predict_perturbation(
                adata,
                perturbation=perturbation,
                cell_type_key='cell_type'
            )
            results['perturbation_effects'] = predictions
            results['top_affected_genes'] = predictions['de_genes'][:50]

    elif task == "network":
        # Gene regulatory network inference
        attention_weights = fm.get_attention_weights(adata)
        grn = fm.infer_grn(attention_weights)

        results['regulatory_network'] = grn
        results['hub_genes'] = grn.get_hub_genes(top_n=20)
        results['key_regulators'] = grn.get_key_regulators()

    elif task == "embed":
        # Extract embeddings for downstream use
        embeddings = fm.encode(adata, return_embeddings=True)
        adata.obsm['X_foundation'] = embeddings

        results['embedding_dim'] = embeddings.shape[1]
        results['embedding_path'] = adata_path.replace('.h5ad', '_embedded.h5ad')
        adata.write(results['embedding_path'])

    return results


# Model selection helper
def select_best_model(task: str, data_size: int, has_reference: bool) -> str:
    """Recommend best foundation model for task."""
    recommendations = {
        ("annotate", True): "geneformer",  # Best with fine-tuning
        ("annotate", False): "genept",      # Good zero-shot
        ("integrate", True): "scgpt",
        ("integrate", False): "scgpt",
        ("perturb", True): "scgpt",
        ("perturb", False): "scfoundation",
        ("network", True): "geneformer",
        ("network", False): "geneformer",
        ("embed", True): "scfoundation",
        ("embed", False): "scfoundation"
    }
    return recommendations.get((task, has_reference), "scfoundation")


# Tool schema for Claude
SCFOUNDATION_TOOL_SCHEMA = {
    "name": "single_cell_foundation_model",
    "description": "Analyze single-cell data using foundation models (Geneformer, scGPT, scFoundation). Supports cell annotation, batch integration, perturbation prediction, and gene network inference.",
    "input_schema": {
        "type": "object",
        "properties": {
            "adata_path": {
                "type": "string",
                "description": "Path to single-cell AnnData file"
            },
            "task": {
                "type": "string",
                "enum": ["annotate", "integrate", "perturb", "network", "embed"],
                "description": "Analysis task"
            },
            "model": {
                "type": "string",
                "enum": ["geneformer", "scgpt", "scfoundation", "scprint", "genept"],
                "description": "Foundation model to use"
            },
            "perturbation": {
                "type": "string",
                "description": "Gene name or drug for perturbation prediction"
            }
        },
        "required": ["adata_path", "task"]
    }
}
```

---

## Prerequisites

### Model Installation

```bash
# Geneformer
pip install geneformer
# Download from Hugging Face: ctheodoris/Geneformer

# scGPT
pip install scgpt
# Download from: https://github.com/bowang-lab/scGPT

# scFoundation
pip install scfoundation
# Download: biomap-research/scfoundation-100m

# scPRINT
pip install scprint
# Download: jkobject/scprint-50m

# GenePT (lightweight)
pip install genept
```

### Hardware Requirements

| Model | GPU VRAM | CPU RAM | Disk |
|-------|----------|---------|------|
| Geneformer | 16 GB | 32 GB | 5 GB |
| scGPT | 24 GB | 32 GB | 8 GB |
| scFoundation | 32 GB | 64 GB | 15 GB |
| scPRINT | 16 GB | 32 GB | 6 GB |
| GenePT | 4 GB | 16 GB | 2 GB |

---

## Methodology

### Architecture Comparison

```
┌─────────────────────────────────────────────────────────────┐
│                Foundation Model Architectures                │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│  Geneformer (BERT-style)          scGPT (GPT-style)         │
│  ┌─────────────────┐              ┌─────────────────┐       │
│  │ Gene Expression │              │ Gene Expression │       │
│  │ → Rank Tokens   │              │ → Value Tokens  │       │
│  └────────┬────────┘              └────────┬────────┘       │
│           ↓                                ↓                 │
│  ┌─────────────────┐              ┌─────────────────┐       │
│  │ Bidirectional   │              │ Autoregressive  │       │
│  │ Transformer     │              │ Transformer     │       │
│  │ [MASK] predict  │              │ Next token      │       │
│  └────────┬────────┘              └────────┬────────┘       │
│           ↓                                ↓                 │
│  Cell Embeddings                  Cell Embeddings           │
│                                   + Generation              │
│                                                              │
└─────────────────────────────────────────────────────────────┘
```

### Training Paradigm

| Model | Pretraining Objective | Key Innovation |
|-------|----------------------|----------------|
| Geneformer | Masked gene prediction | Rank-value encoding |
| scGPT | Next gene prediction | Cell prompt tokens |
| scFoundation | Contrastive + masked | Scale (100M cells) |
| scPRINT | Multi-task | Denoising + annotation |
| GenePT | None (uses LLM) | Literature embedding |

---

## Benchmarks

### Zero-Shot Cell Annotation (Tabula Sapiens)

| Model | Accuracy | Macro F1 | Rare Cell F1 |
|-------|----------|----------|--------------|
| scFoundation | 0.82 | 0.78 | 0.65 |
| Geneformer | 0.79 | 0.74 | 0.58 |
| scGPT | 0.77 | 0.71 | 0.54 |
| GenePT | 0.81 | 0.76 | 0.62 |
| HVG + scVI | 0.80 | 0.75 | 0.60 |

### Fine-Tuned Performance

| Task | Model | Metric | Score |
|------|-------|--------|-------|
| Annotation | Geneformer-FT | Accuracy | 0.94 |
| Integration | scGPT | ASW | 0.85 |
| Perturbation | scGPT | Pearson | 0.72 |
| GRN | Geneformer | AUROC | 0.78 |

---

## Related Skills

- `biomedical.foundation_models.nicheformer` - Spatial single-cell
- `biomedical.genomics.single_cell_rna_qc` - QC preprocessing
- `biomedical.genomics.spatial_transcriptomics` - Spatial analysis
- `biomedical.research_tools.cellagent` - Single-cell workflows

---

## References

1. **Geneformer (2023)**: "Transfer learning enables predictions in network biology." *Nature*. [DOI: 10.1038/s41586-023-06139-9](https://www.nature.com/articles/s41586-023-06139-9)

2. **scGPT (2024)**: "scGPT: toward building a foundation model for single-cell multi-omics using generative AI." *Nature Methods*.

3. **scFoundation (2024)**: "Large-scale foundation model on single-cell transcriptomics." *Nature Methods*.

4. **GenePT (2023)**: "A Simple But Effective Foundation Model for Genes and Cells Built From ChatGPT." *bioRxiv*.

5. **Benchmark (2025)**: "Zero-shot evaluation reveals limitations of single-cell foundation models." *Genome Biology*.

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

---

*Last updated: December 2025*

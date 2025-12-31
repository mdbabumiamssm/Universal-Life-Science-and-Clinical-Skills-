# Nicheformer Agent

**ID:** `biomedical.foundation_models.nicheformer`
**Version:** 1.0.0
**Status:** Production
**Category:** Foundation Models / Spatial Single-Cell AI
**Released:** December 2025

---

## Overview

The **Nicheformer Agent** provides access to Nicheformer, the first large-scale foundation model that integrates single-cell analysis with spatial transcriptomics. Developed by Helmholtz Munich and TUM, Nicheformer was pretrained on **SpatialCorpus-110M** - a curated collection of over 57 million dissociated and 53 million spatially resolved cells across 73 tissues.

Unlike traditional single-cell models that analyze cells in isolation, Nicheformer learns directly from spatial organization, enabling reconstruction of how cells sense and influence their neighbors. This represents the foundation of a "Virtual Cell and Tissue model."

### Why Nicheformer?

| Challenge | Traditional Approach | Nicheformer Solution |
|-----------|---------------------|---------------------|
| Loss of spatial context | Dissociation destroys tissue architecture | Learns from spatial organization |
| Limited tissue coverage | Single tissue/organ models | 73 tissues, human + mouse |
| Cell-cell interactions | Post-hoc inference | Native spatial learning |
| Transfer to new tissues | Requires retraining | Zero-shot generalization |

---

## Key Capabilities

### 1. Core Model Features

| Feature | Specification | Notes |
|---------|---------------|-------|
| **Architecture** | Transformer-based | Attention over cells and genes |
| **Pretraining Data** | 110M+ cells | 57M dissociated + 53M spatial |
| **Species** | Human + Mouse | Cross-species transfer |
| **Tissues** | 73 | Comprehensive organ coverage |
| **Tasks** | 6+ downstream | Zero-shot and fine-tuning |

### 2. Downstream Tasks

| Task | Description | Performance |
|------|-------------|-------------|
| **Cell Type Annotation** | Automatic cell labeling | State-of-the-art accuracy |
| **Spatial Composition Prediction** | Predict cell neighborhoods | Superior to baselines |
| **Spatial Label Prediction** | Annotate spatial data | Excellent fine-tuning |
| **Gene Expression Imputation** | Fill missing genes | Cross-platform transfer |
| **Batch Integration** | Merge datasets | Preserves biology |
| **Perturbation Prediction** | Predict treatment effects | Novel applications |

### 3. Technical Specifications

```yaml
Model Architecture:
  Type: Transformer encoder
  Attention: Multi-head self-attention
  Embedding: Gene rank-based tokens
  Context: Cell + spatial neighborhood

Training:
  Pretraining: Masked gene prediction
  Spatial: Neighbor context encoding
  Augmentation: Random masking, dropout

Inference:
  Zero-shot: Direct embedding use
  Fine-tuning: Task-specific heads
  Memory: GPU with 16GB+ recommended
```

---

## Usage

### Example Prompt

```text
I have spatial transcriptomics data from a human liver tumor sample using
10x Visium. I need to:
1. Annotate cell types in tissue context
2. Identify tumor microenvironment niches
3. Predict cell-cell communication patterns
4. Compare with dissociated single-cell reference

Please use Nicheformer for spatial-aware analysis.
```

### Expected Output

```markdown
## Nicheformer Spatial Analysis: Liver Tumor

### Data Summary
| Metric | Value |
|--------|-------|
| Platform | 10x Visium |
| Spots | 4,892 |
| Genes | 18,456 |
| Tissue | Liver (tumor + margin) |

### Cell Type Deconvolution

| Region | Dominant Cell Types | Confidence |
|--------|--------------------|-----------|
| **Tumor Core** | Malignant hepatocytes (68%), CAFs (15%), TAMs (12%) | 0.94 |
| **Tumor Margin** | Normal hepatocytes (45%), Immune infiltrate (35%) | 0.91 |
| **Stroma** | Stellate cells (40%), Endothelial (30%), Fibroblasts (25%) | 0.89 |
| **Portal Tracts** | Cholangiocytes (55%), Hepatocytes (30%) | 0.92 |

### Spatial Niche Identification

```
Niche 1: Immunosuppressive TME
├── Location: Tumor core
├── Composition: TAMs (M2), Tregs, exhausted CD8+
├── Ligand-Receptor: PD-L1/PD-1, CTLA4/CD80
└── Clinical: Poor prognosis signature

Niche 2: Active Immune Interface
├── Location: Tumor margin
├── Composition: Activated CD8+, NK cells, M1 macrophages
├── Ligand-Receptor: CXCL9/CXCR3, IFNg signaling
└── Clinical: Immunotherapy responsive

Niche 3: Fibrotic Capsule
├── Location: Tumor periphery
├── Composition: CAFs, stellate cells, ECM-rich
├── Ligand-Receptor: TGFb/TGFbR, COL1A1
└── Clinical: Metastasis barrier
```

### Cell-Cell Communication

| Sender | Receiver | Pathway | Score | Spatial Constraint |
|--------|----------|---------|-------|-------------------|
| TAMs | Malignant | VEGF | 0.89 | Adjacent spots |
| CAFs | Malignant | HGF/MET | 0.85 | <100μm |
| CD8+ | Malignant | FASL/FAS | 0.78 | Contact required |
| Tumor | Tregs | CCL22/CCR4 | 0.82 | Chemotaxis gradient |

### Integration with scRNA-seq Reference

| Metric | Result |
|--------|--------|
| Reference cells | 45,000 (healthy liver atlas) |
| Mapping accuracy | 94.2% |
| Novel populations | 2 (tumor-specific) |
| Batch effect | Corrected (UMAP overlapping) |
```

---

## LLM Agent Integration

### Python Tool Implementation

```python
from typing import Optional, List, Dict, Any
import numpy as np
import scanpy as sc
import squidpy as sq

# Tool decorator for LangChain/OpenAI function calling
def nicheformer_tool(
    adata_path: str,
    task: str,
    spatial_key: str = "spatial",
    reference_path: Optional[str] = None,
    n_neighbors: int = 15,
    model_checkpoint: str = "nicheformer-base",
    device: str = "cuda"
) -> Dict[str, Any]:
    """
    Nicheformer foundation model for spatial single-cell analysis.

    Args:
        adata_path: Path to AnnData object (.h5ad)
        task: One of 'annotate', 'niche', 'impute', 'integrate', 'communicate'
        spatial_key: Key for spatial coordinates in adata.obsm
        reference_path: Optional reference dataset for integration
        n_neighbors: Spatial neighborhood size
        model_checkpoint: Model version to use
        device: 'cuda' or 'cpu'

    Returns:
        Dictionary with task-specific results
    """
    from nicheformer import NicheformerModel, SpatialProcessor

    # Load data
    adata = sc.read_h5ad(adata_path)

    # Initialize model
    model = NicheformerModel.from_pretrained(model_checkpoint)
    model.to(device)

    # Process spatial graph
    processor = SpatialProcessor(n_neighbors=n_neighbors)
    spatial_graph = processor.build_graph(adata, spatial_key=spatial_key)

    results = {}

    if task == "annotate":
        # Zero-shot cell type annotation
        embeddings = model.encode(adata, spatial_graph=spatial_graph)
        predictions = model.predict_celltypes(embeddings)
        adata.obs['nicheformer_celltype'] = predictions['labels']
        adata.obs['nicheformer_confidence'] = predictions['confidence']
        results['annotations'] = predictions

    elif task == "niche":
        # Spatial niche identification
        embeddings = model.encode(adata, spatial_graph=spatial_graph)
        niches = model.identify_niches(embeddings, spatial_graph)
        adata.obs['niche_cluster'] = niches['labels']
        results['niches'] = niches
        results['niche_composition'] = niches['composition']

    elif task == "impute":
        # Gene expression imputation
        imputed = model.impute_expression(adata, spatial_graph=spatial_graph)
        results['imputed_genes'] = imputed['gene_names']
        results['imputation_quality'] = imputed['r2_scores']

    elif task == "integrate":
        # Reference integration
        if reference_path:
            ref_adata = sc.read_h5ad(reference_path)
            integrated = model.integrate(adata, ref_adata, spatial_graph=spatial_graph)
            results['integration'] = integrated

    elif task == "communicate":
        # Cell-cell communication with spatial constraints
        embeddings = model.encode(adata, spatial_graph=spatial_graph)
        interactions = model.predict_interactions(
            embeddings,
            spatial_graph=spatial_graph,
            lr_database="cellchatdb"
        )
        results['interactions'] = interactions

    return results


# Anthropic Claude tool schema
NICHEFORMER_TOOL_SCHEMA = {
    "name": "nicheformer_analysis",
    "description": "Run Nicheformer foundation model for spatial transcriptomics analysis including cell type annotation, niche identification, and cell-cell communication",
    "input_schema": {
        "type": "object",
        "properties": {
            "adata_path": {
                "type": "string",
                "description": "Path to spatial transcriptomics AnnData file"
            },
            "task": {
                "type": "string",
                "enum": ["annotate", "niche", "impute", "integrate", "communicate"],
                "description": "Analysis task to perform"
            },
            "reference_path": {
                "type": "string",
                "description": "Optional path to reference single-cell dataset"
            }
        },
        "required": ["adata_path", "task"]
    }
}
```

### OpenAI Function Calling

```python
import openai

tools = [{
    "type": "function",
    "function": {
        "name": "run_nicheformer",
        "description": "Analyze spatial transcriptomics data using Nicheformer foundation model",
        "parameters": {
            "type": "object",
            "properties": {
                "data_path": {"type": "string", "description": "Path to .h5ad file"},
                "analysis_type": {
                    "type": "string",
                    "enum": ["cell_annotation", "niche_discovery", "spatial_communication"],
                    "description": "Type of spatial analysis"
                },
                "tissue_type": {"type": "string", "description": "Tissue context for model"}
            },
            "required": ["data_path", "analysis_type"]
        }
    }
}]

response = openai.chat.completions.create(
    model="gpt-4-turbo",
    messages=[{"role": "user", "content": "Analyze my liver spatial data"}],
    tools=tools,
    tool_choice="auto"
)
```

---

## Prerequisites

### Required Software

| Component | Version | Installation |
|-----------|---------|--------------|
| Python | >=3.9 | System |
| PyTorch | >=2.0 | `pip install torch` |
| Scanpy | >=1.9 | `pip install scanpy` |
| Squidpy | >=1.3 | `pip install squidpy` |
| Nicheformer | >=1.0 | `pip install nicheformer` |

### Hardware Requirements

| Task | GPU Memory | CPU | RAM |
|------|------------|-----|-----|
| Inference (small) | 8 GB | 4 cores | 16 GB |
| Inference (large) | 16 GB | 8 cores | 32 GB |
| Fine-tuning | 24+ GB | 16 cores | 64 GB |

### Data Requirements

| Platform | Format | Notes |
|----------|--------|-------|
| 10x Visium | .h5ad | Spatial coords in obsm |
| MERFISH | .h5ad | Single-cell resolution |
| Slide-seq | .h5ad | Bead coordinates |
| CosMx | .h5ad | Subcellular resolution |
| Xenium | .h5ad | In situ sequencing |

---

## Methodology

### Model Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                    Nicheformer Architecture                  │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│  Input: Gene Expression + Spatial Coordinates                │
│           ↓                                                  │
│  ┌─────────────────┐    ┌─────────────────┐                 │
│  │ Gene Tokenizer  │    │ Spatial Encoder │                 │
│  │ (Rank-based)    │    │ (Graph-based)   │                 │
│  └────────┬────────┘    └────────┬────────┘                 │
│           │                      │                           │
│           └──────────┬───────────┘                          │
│                      ↓                                       │
│           ┌─────────────────────┐                           │
│           │ Transformer Encoder │                           │
│           │ (Multi-head Attn)   │                           │
│           └──────────┬──────────┘                           │
│                      ↓                                       │
│           ┌─────────────────────┐                           │
│           │  Cell Embeddings    │                           │
│           │  (Spatial-aware)    │                           │
│           └──────────┬──────────┘                           │
│                      ↓                                       │
│    ┌─────────┬──────┴──────┬─────────┬─────────┐           │
│    ↓         ↓             ↓         ↓         ↓            │
│ Annotate  Niche ID    Integrate  Impute  Communicate        │
│                                                              │
└─────────────────────────────────────────────────────────────┘
```

### Training Procedure

1. **Pretraining Corpus**: SpatialCorpus-110M across 73 tissues
2. **Masked Gene Prediction**: Random 15% gene masking
3. **Spatial Context**: Neighborhood aggregation via graph attention
4. **Cross-species Transfer**: Joint human-mouse pretraining

### Evaluation Benchmarks

| Task | Dataset | Metric | Performance |
|------|---------|--------|-------------|
| Cell Annotation | Tabula Sapiens | Accuracy | 0.94 |
| Niche Discovery | MERFISH Brain | ARI | 0.87 |
| Integration | Multi-platform | LISI | 0.92 |
| Imputation | Held-out genes | Pearson r | 0.85 |

---

## Comparison with Alternatives

| Model | Spatial | Scale | Zero-shot | Publication |
|-------|---------|-------|-----------|-------------|
| **Nicheformer** | Native | 110M | Yes | Nature Methods 2025 |
| scGPT | Post-hoc | 33M | Yes | Nature Methods 2024 |
| Geneformer | No | 95M | Yes | Nature 2023 |
| scFoundation | No | 100M | Yes | Nature Methods 2024 |
| CellPLM | No | 50M | Yes | bioRxiv 2024 |

---

## Related Skills

- `biomedical.genomics.spatial_transcriptomics` - STAgent for spatial analysis
- `biomedical.genomics.single_cell` - CellAgent for dissociated cells
- `biomedical.foundation_models.scfoundation` - Non-spatial foundation model
- `biomedical.foundation_models.genept` - ChatGPT-based gene embeddings

---

## References

1. **Nicheformer (2025)**: "Nicheformer: a foundation model for single-cell and spatial omics." *Nature Methods*. [DOI: 10.1038/s41592-025-02814-z](https://www.nature.com/articles/s41592-025-02814-z)

2. **Helmholtz Press Release (2025)**: "New Foundation Model Reveals How Cells Are Organized in Tissues." Helmholtz Munich.

3. **Theis Lab**: [https://github.com/theislab/nicheformer](https://github.com/theislab/nicheformer)

4. **SpatialCorpus-110M**: Curated pretraining dataset documentation.

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

---

*Last updated: December 2025*

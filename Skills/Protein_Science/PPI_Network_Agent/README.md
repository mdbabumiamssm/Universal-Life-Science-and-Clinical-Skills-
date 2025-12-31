# Protein-Protein Interaction Network Agent

**ID:** `biomedical.protein_science.ppi_network`
**Version:** 1.0.0
**Status:** Beta
**Category:** Protein Science / Systems Biology

---

## Overview

The **Protein-Protein Interaction Network Agent** predicts and analyzes protein-protein interactions (PPIs) using deep learning architectures including Graph Neural Networks (GNNs), Transformers, and protein language models (ESM). PPIs are fundamental regulators of biological functions, and understanding their networks is crucial for drug discovery, disease mechanism elucidation, and systems biology.

This agent integrates multiple prediction modalities—sequence-based, structure-based, and network-based—to provide comprehensive PPI analysis with explainable predictions.

---

## Key Capabilities

### 1. PPI Prediction Methods

| Method | Input | Architecture | Best For |
|--------|-------|--------------|----------|
| **Sequence-based** | Protein sequences | BERT, ESM embeddings | Genome-wide screening |
| **Structure-based** | 3D structures | GNN, 3D-CNN | High-confidence predictions |
| **Network-based** | Known PPI graph | Graph Transformers | Link prediction |
| **Multimodal** | Sequence + Structure | Fusion models | Maximum accuracy |

### 2. Interaction Types

- **Physical interactions:** Direct binding (co-IP validated)
- **Functional associations:** Pathway co-membership
- **Genetic interactions:** Synthetic lethality
- **Regulatory relationships:** Kinase-substrate, TF-target

### 3. Network Analysis

- **Community detection:** Protein complex identification
- **Centrality analysis:** Hub protein identification
- **Pathway enrichment:** Biological function annotation
- **Disease module mapping:** Disease-associated subnetworks

### 4. Protein Complex Prediction

- **Co-fractionation MS integration:** Elution profile correlation
- **AP-MS analysis:** Bait-prey scoring
- **Complex membership:** CORUM, ComplexPortal matching
- **Novel complex discovery:** Clustering-based identification

---

## Usage

### Example Prompt

```text
Predict protein-protein interactions for human TP53 (P53_HUMAN).
Use both sequence and structure-based methods.
Identify the top 20 novel interactors not in STRING database.
Map predicted interactions to cancer-related pathways.
```

### Expected Output

```
## PPI Prediction Report: TP53 (P04637)

### Known Interactions (STRING, score > 0.9)
| Partner | Gene | Score | Evidence |
|---------|------|-------|----------|
| P38936 | CDKN1A | 0.999 | Experimental + Text mining |
| Q00987 | MDM2 | 0.999 | Experimental + Coexpression |
| P10276 | RARA | 0.986 | Experimental |
| Q13315 | ATM | 0.978 | Experimental |

### Novel Predicted Interactions (Not in STRING)

| Partner | Gene | Seq Score | Struct Score | Combined | Confidence |
|---------|------|-----------|--------------|----------|------------|
| Q9Y6K1 | DNM3 | 0.89 | 0.84 | 0.92 | High |
| P55081 | MFAP1 | 0.86 | 0.81 | 0.89 | High |
| Q9NR30 | DDX21 | 0.84 | 0.78 | 0.87 | Medium |
| Q9ULZ3 | ASF1A | 0.82 | 0.85 | 0.86 | Medium |
...

### Predicted Interaction Interface (TP53-DNM3)
- **TP53 residues:** 102-112 (DNA binding domain loop)
- **DNM3 residues:** 445-460 (GTPase domain)
- **Interface area:** 1,240 Å²
- **Binding energy (predicted):** -8.4 kcal/mol

### Pathway Enrichment (Novel Interactors)
| Pathway | P-value | Genes |
|---------|---------|-------|
| mRNA splicing | 2.3e-5 | DDX21, MFAP1, SRSF2 |
| Chromatin remodeling | 4.1e-4 | ASF1A, SMARCA4 |
| DNA damage response | 8.7e-4 | PARP2, RNF168 |

### Network Visualization
[Network diagram showing TP53 at center with predicted interactions]
```

### LLM Agent Integration

```python
@tool
def predict_ppi(
    protein_id: str,
    method: str = "multimodal",
    num_predictions: int = 20,
    exclude_known: bool = True,
    species: str = "human"
) -> str:
    """
    Predicts protein-protein interactions.

    Args:
        protein_id: UniProt ID of query protein
        method: sequence, structure, network, or multimodal
        num_predictions: Number of top predictions
        exclude_known: Exclude interactions in STRING/BioGRID
        species: Target species

    Returns:
        Ranked list of predicted interactors with confidence
    """
    pass


@tool
def analyze_ppi_network(
    protein_list: list[str],
    analysis_type: str = "community",
    enrichment: bool = True
) -> str:
    """
    Analyzes a protein-protein interaction network.

    Args:
        protein_list: List of UniProt IDs
        analysis_type: community, centrality, or pathway
        enrichment: Perform pathway enrichment

    Returns:
        Network analysis report
    """
    pass
```

---

## Prerequisites

### Required Databases/APIs

| Resource | Purpose | Access |
|----------|---------|--------|
| **STRING** | Known interactions | API v11.5 |
| **BioGRID** | Curated interactions | REST API |
| **IntAct** | Experimental PPIs | PSICQUIC |
| **ESM-2** | Sequence embeddings | HuggingFace |
| **AlphaFold DB** | Predicted structures | EBI API |

### Dependencies

```
torch>=2.0
torch_geometric>=2.3
esm>=2.0
networkx>=3.0
pandas>=1.5
requests>=2.28
```

---

## Methodology

### Deep Learning Architectures

```
Query Protein
    ↓
Feature Extraction
    ├── ESM-2 Embeddings (sequence)
    ├── AlphaFold Structure (GNN features)
    └── Known Network (graph embeddings)
    ↓
Pairwise Scoring
    ├── Sequence: Siamese Transformer
    ├── Structure: GNN message passing
    └── Network: Link prediction GNN
    ↓
Score Fusion (weighted average or MLP)
    ↓
Ranking + Confidence Estimation
    ↓
Biological Validation
    ├── Pathway co-membership
    ├── Coexpression correlation
    └── Literature evidence
```

### Model Architectures

| Approach | Architecture | Key Features |
|----------|--------------|--------------|
| **D-SCRIPT** | Siamese LSTM | Sequence-only, transfer learning |
| **DeepPPI** | CNN + Attention | Amino acid properties |
| **PIPR** | GNN on residue graphs | Structure-aware |
| **SpatialPPI** | 3D equivariant GNN | Full structural context |
| **PPI-Transformer** | Cross-attention | Multimodal fusion |

### Confidence Scoring

```python
def calculate_ppi_confidence(
    seq_score: float,
    struct_score: float,
    network_score: float,
    coexpression: float = None
) -> tuple[float, str]:
    """
    Calculate combined PPI confidence score.

    Returns confidence score and level (High/Medium/Low).
    """
    # Weighted combination
    weights = {'seq': 0.3, 'struct': 0.4, 'network': 0.2, 'coexp': 0.1}

    combined = (
        weights['seq'] * seq_score +
        weights['struct'] * struct_score +
        weights['network'] * network_score
    )

    if coexpression is not None:
        combined += weights['coexp'] * coexpression
        combined /= sum(weights.values())
    else:
        combined /= (weights['seq'] + weights['struct'] + weights['network'])

    if combined > 0.85:
        return combined, "High"
    elif combined > 0.70:
        return combined, "Medium"
    else:
        return combined, "Low"
```

---

## Related Skills

- **Proteomics MS Agent:** Identify PPIs from AP-MS/co-IP data
- **Knowledge Graphs (KG-RAG):** PPI network integration
- **Drug Target Interaction Agent:** Target-based drug discovery

---

## References

- **Sledzieski et al. (2021):** "D-SCRIPT translates genome to phenome with sequence-based, structure-aware, genome-scale predictions of protein-protein interactions." *Cell Systems*
- **Bryant et al. (2022):** "Improved prediction of protein-protein interactions using AlphaFold2." *Nature Communications*
- [STRING Database](https://string-db.org/)
- [BioGRID](https://thebiogrid.org/)

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

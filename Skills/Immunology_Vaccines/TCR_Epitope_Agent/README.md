# TCR-Epitope Interaction Prediction Agent

**ID:** `biomedical.immunology.tcr_epitope`
**Version:** 1.0.0
**Status:** Experimental
**Category:** Immunology / T-Cell Immunology

---

## Overview

The **TCR-Epitope Interaction Prediction Agent** predicts T-cell receptor (TCR) recognition of peptide-MHC (pMHC) complexes, a critical determinant of adaptive immune responses. While MHC binding is necessary for antigen presentation, TCR recognition determines whether an immune response is actually triggered.

This agent uses deep learning approaches including attention mechanisms and 3D structural modeling to predict TCR-pMHC interactions, enabling improved vaccine design and TCR-based immunotherapy development.

---

## Key Capabilities

### 1. TCR-pMHC Binding Prediction

| Prediction Type | Input | Output |
|-----------------|-------|--------|
| **Binary classification** | TCR + epitope | Binding yes/no |
| **Affinity regression** | TCR + pMHC | Binding strength |
| **Cross-reactivity** | TCR + epitope set | Multi-epitope binding |
| **TCR clustering** | TCR repertoire | Specificity groups |

### 2. TCR Representation Methods

| Method | Features | Best For |
|--------|----------|----------|
| **CDR3β only** | 10-20 AA sequence | Fast screening |
| **Full TCR** | α/β chains complete | Accuracy |
| **TCR + pMHC structure** | 3D complex | High confidence |
| **Repertoire embedding** | Population-level | Immune profiling |

### 3. Applications

- **Vaccine immunogenicity:** Will epitope trigger T-cell response?
- **Neoantigen validation:** TCR recognition of tumor antigens
- **TCR-T therapy:** Identify reactive TCRs for adoptive transfer
- **Autoimmunity:** Predict self-reactive TCR-peptide pairs

### 4. Model Architectures

| Model | Architecture | Performance |
|-------|--------------|-------------|
| **ERGO** | LSTM + attention | AUC 0.85 |
| **NetTCR** | CNN | AUC 0.82 |
| **TITAN** | Transformer | AUC 0.88 |
| **pMTnet** | Multi-task learning | AUC 0.87 |

---

## Usage

### Example Prompt

```text
Predict TCR recognition for the SARS-CoV-2 epitope YLQPRTFLL (HLA-A*02:01).

Given these CDR3β sequences from COVID-19 recovered patients:
1. CASSIRSSYEQYF
2. CASSLGQAYEQYF
3. CASSYSGGTEAFF
4. CASSFSTGTEAFF
5. CASSLGQTNTEAFF

Predict which TCRs are likely to recognize the epitope.
Also identify the most important CDR3 residues for binding.
```

### Expected Output

```
## TCR-Epitope Interaction Prediction

### Target Epitope
- **Sequence:** YLQPRTFLL
- **HLA restriction:** HLA-A*02:01
- **Source:** SARS-CoV-2 Spike protein (269-277)
- **Known immunogenicity:** Confirmed in multiple studies

### TCR Binding Predictions

| TCR ID | CDR3β Sequence | Prediction | Confidence | Evidence |
|--------|----------------|------------|------------|----------|
| TCR-1 | CASSIRSSYEQYF | **Binding** | 0.91 | High |
| TCR-2 | CASSLGQAYEQYF | **Binding** | 0.87 | High |
| TCR-3 | CASSYSGGTEAFF | Non-binding | 0.34 | Medium |
| TCR-4 | CASSFSTGTEAFF | Non-binding | 0.28 | Medium |
| TCR-5 | CASSLGQTNTEAFF | **Binding** | 0.72 | Medium |

### Binding TCR Analysis

#### TCR-1: CASSIRSSYEQYF (Highest confidence)

##### Sequence Features
```
Position:  1 2 3 4 5 6 7 8 9 10 11 12 13
CDR3β:     C A S S I R S S Y E  Q  Y  F
Importance:- - + + + + + + + +  -  -  -
           └─conserved─┘ └─variable region─┘
```

##### Key Binding Residues
| Position | Residue | Importance | Role |
|----------|---------|------------|------|
| 5 | Ile (I) | 0.89 | pMHC contact |
| 6 | Arg (R) | 0.92 | Epitope interaction |
| 7-8 | Ser-Ser | 0.78 | Flexibility |
| 9 | Tyr (Y) | 0.85 | Aromatic stacking |

##### Structural Prediction
- **Predicted binding mode:** CDR3β loop contacts peptide P4-P6
- **Critical contacts:** R6→pMHC groove, Y9→peptide backbone
- **Estimated affinity:** Kd ~ 5-20 μM (typical range)

#### Motif Analysis (Binding TCRs)

##### Shared CDR3β Motifs
```
Binding TCRs consensus:
Position:  4  5  6  7  8  9  10
Motif:     S  [ILV]  [RK]  .  .  [YF]  [EQ]
           │    │     │        │    │
           │    │     │        │    └─ Acidic/charged
           │    │     │        └─ Aromatic (π-stacking)
           │    │     └─ Positive charge (pMHC contact)
           │    └─ Hydrophobic (core packing)
           └─ Conserved Ser
```

##### Sequence Logo (Binding TCRs)
```
Bits
2.0 |         Y
1.5 |    I R    Q
1.0 | C S   S S   Y F
0.5 | A     G      E
    +---+---+---+---+---
       4   6   8  10  12
```

### TCR Repertoire Features

| Metric | Binding TCRs | Non-binding TCRs |
|--------|--------------|------------------|
| Mean CDR3 length | 13.3 | 13.5 |
| Mean hydrophobicity | -0.42 | -0.51 |
| Aromatic content | 15.4% | 8.2% |
| Charged residues | 18.2% | 14.1% |

### Validation Recommendations

1. **Tetramer staining:** Validate TCR-1, TCR-2, TCR-5 with YLQPRTFLL/A*02:01 tetramer
2. **Functional assay:** IFN-γ ELISpot with peptide-pulsed APCs
3. **Structural validation:** Crystal structure of TCR-1 + pMHC (highest confidence)
```

### LLM Agent Integration

```python
@tool
def predict_tcr_epitope_binding(
    tcr_sequences: list[str],
    epitope: str,
    hla_allele: str,
    include_structural: bool = False
) -> str:
    """
    Predicts TCR recognition of peptide-MHC complex.

    Args:
        tcr_sequences: List of CDR3β sequences
        epitope: Peptide epitope sequence
        hla_allele: HLA allele for pMHC complex
        include_structural: Include 3D modeling analysis

    Returns:
        Binding predictions with confidence scores
    """
    pass


@tool
def analyze_tcr_repertoire(
    repertoire_file: str,
    target_epitopes: list[str],
    hla_alleles: list[str]
) -> str:
    """
    Analyzes TCR repertoire for epitope-specific responses.

    Args:
        repertoire_file: Path to TCR-seq data (AIRR format)
        target_epitopes: List of epitopes to screen
        hla_alleles: Patient HLA alleles

    Returns:
        Epitope-reactive TCR identification and clustering
    """
    pass
```

---

## Prerequisites

### Required Tools/Databases

| Resource | Purpose | Access |
|----------|---------|--------|
| **VDJdb** | TCR-epitope database | Public |
| **McPAS-TCR** | Pathology-associated TCRs | Public |
| **IEDB** | Immune epitope database | Public |
| **TCRdist** | TCR similarity analysis | Open source |
| **AlphaFold** | Structure prediction | API |

### Dependencies

```
torch>=2.0
transformers>=4.30
tape-proteins>=0.5  # Protein language models
pandas>=1.5
numpy>=1.24
scikit-learn>=1.3
```

---

## Methodology

### TCR-pMHC Prediction Model

```
Input: CDR3β sequence + Epitope + HLA
    ↓
Encoding
    ├── TCR: ESM-2 embeddings or one-hot
    ├── Epitope: ESM-2 embeddings
    └── HLA: Learned embeddings per allele
    ↓
Cross-Attention Module
    ├── TCR attends to epitope
    └── Epitope attends to TCR
    ↓
Interaction Layer
    ├── Element-wise product
    ├── Concatenation
    └── Bilinear transformation
    ↓
Prediction Heads
    ├── Binary: Binding classification
    └── Regression: Affinity prediction
```

### Model Architecture

```python
import torch
import torch.nn as nn

class TCREpitopePredictor(nn.Module):
    def __init__(self, embed_dim=256, num_heads=8):
        super().__init__()

        # Sequence encoders
        self.tcr_encoder = nn.LSTM(21, embed_dim, bidirectional=True, batch_first=True)
        self.epitope_encoder = nn.LSTM(21, embed_dim, bidirectional=True, batch_first=True)

        # Cross-attention
        self.cross_attention = nn.MultiheadAttention(embed_dim * 2, num_heads)

        # Prediction head
        self.classifier = nn.Sequential(
            nn.Linear(embed_dim * 4, 256),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(256, 64),
            nn.ReLU(),
            nn.Linear(64, 1),
            nn.Sigmoid()
        )

    def forward(self, tcr_seq, epitope_seq):
        # Encode sequences
        tcr_embed, _ = self.tcr_encoder(tcr_seq)
        epi_embed, _ = self.epitope_encoder(epitope_seq)

        # Cross-attention
        tcr_attended, _ = self.cross_attention(
            tcr_embed, epi_embed, epi_embed
        )

        # Pool and classify
        tcr_pooled = tcr_attended.mean(dim=1)
        epi_pooled = epi_embed.mean(dim=1)

        combined = torch.cat([tcr_pooled, epi_pooled], dim=-1)
        return self.classifier(combined)
```

---

## Related Skills

- **Epitope Prediction Agent:** MHC binding prediction
- **Neoantigen Vaccine Agent:** Neoantigen immunogenicity
- **Single-Cell Analysis:** TCR-seq from scRNA-seq

---

## References

- **Springer et al. (2020):** "Prediction of specific TCR-peptide binding from large dictionaries." *Frontiers in Immunology*
- **Montemurro et al. (2021):** "NetTCR-2.0 enables accurate prediction of TCR-peptide binding." *Communications Biology*
- [VDJdb Database](https://vdjdb.cdr3.net/)
- [TITAN Model](https://github.com/PaccMann/TITAN)

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

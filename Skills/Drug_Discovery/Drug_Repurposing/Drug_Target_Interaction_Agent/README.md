# Drug-Target Interaction Prediction Agent

**ID:** `biomedical.drug_discovery.dti_prediction`
**Version:** 1.0.0
**Status:** Production
**Category:** Drug Discovery / Target Identification

---

## Overview

The **Drug-Target Interaction (DTI) Prediction Agent** predicts binding interactions between drugs and protein targets using multimodal deep learning. This enables virtual screening, target identification for phenotypic hits, and polypharmacology prediction.

---

## Key Capabilities

### 1. Model Architectures

| Model | Architecture | Features |
|-------|--------------|----------|
| **DeepDTA** | CNN | Sequence-based |
| **GraphDTA** | GNN + CNN | Molecular graphs |
| **ERT-GFAN** | Transformer + GNN | Multimodal |
| **DrugBAN** | Bilinear attention | Interpretable |

### 2. Input Representations

| Drug | Target | Method |
|------|--------|--------|
| SMILES | Sequence | Sequence encoding |
| Morgan FP | ESM embedding | Pretrained |
| 3D structure | AlphaFold | Structure-based |

### 3. Output Types

- **Binary:** Binding yes/no
- **Affinity:** Kd, Ki, IC50 prediction
- **Binding site:** Pocket identification
- **Selectivity:** Off-target prediction

---

## Usage

### Example Prompt

```text
Predict drug-target interactions for a phenotypic hit compound.

Compound: SMILES: CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F

Tasks:
1. Predict top 10 protein targets
2. Estimate binding affinities
3. Identify potential selectivity issues
```

### Expected Output

```
## DTI Prediction Report

### Compound Analysis
- **SMILES:** CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F
- **Name:** Celecoxib (identified)
- **MW:** 381.4

### Predicted Targets

| Rank | Target | Gene | Affinity (pKi) | Confidence |
|------|--------|------|----------------|------------|
| 1 | **COX-2** | PTGS2 | 8.9 | 0.98 |
| 2 | COX-1 | PTGS1 | 6.2 | 0.85 |
| 3 | Carbonic anhydrase II | CA2 | 7.1 | 0.72 |
| 4 | PDE4 | PDE4A | 5.8 | 0.61 |
| 5 | PKC alpha | PRKCA | 5.4 | 0.54 |

### Selectivity Analysis
| Ratio | Value | Interpretation |
|-------|-------|----------------|
| COX-2/COX-1 | 500x | Highly selective |
| COX-2/CA2 | 60x | Good selectivity |

### Known Validation
- COX-2 is known target (validated)
- CA2 binding confirmed in literature
- PKC binding: potential off-target to monitor
```

---

## References

- **Öztürk et al. (2018):** "DeepDTA: deep drug-target binding affinity prediction." *Bioinformatics*
- **Nguyen et al. (2021):** "GraphDTA: predicting drug-target binding affinity with graph neural networks."

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

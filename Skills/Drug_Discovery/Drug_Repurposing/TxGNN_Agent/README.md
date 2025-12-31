# TxGNN Drug Repurposing Agent

**ID:** `biomedical.drug_discovery.txgnn`
**Version:** 1.0.0
**Status:** Production
**Category:** Drug Discovery / Drug Repurposing

---

## Overview

The **TxGNN Drug Repurposing Agent** implements TxGNN (Therapeutic Graph Neural Network), a zero-shot drug repurposing model that predicts therapeutic indications and contraindications using large-scale medical knowledge graphs. TxGNN can predict applications for diseases without existing treatments, addressing a major challenge in drug repurposing.

---

## Key Capabilities

### 1. Prediction Types

| Type | Description | Application |
|------|-------------|-------------|
| **Indication** | Drug treats disease | Repurposing opportunities |
| **Contraindication** | Drug harmful for disease | Safety assessment |
| **Zero-shot** | No training examples | Rare diseases |

### 2. Knowledge Graph Integration

| Data Source | Entities | Relationships |
|-------------|----------|---------------|
| **DrugBank** | Drugs, targets | Drug-target binding |
| **OMIM** | Diseases, genes | Gene-disease association |
| **STRING** | Proteins | PPI network |
| **Hetionet** | Multi-type | 24 edge types |

### 3. Model Features

- **Graph attention:** Learn relationship importance
- **Metapath encoding:** Capture indirect relationships
- **Uncertainty quantification:** Confidence estimation
- **Explainability:** Path-based explanations

---

## Usage

### Example Prompt

```text
Use TxGNN to identify drug repurposing opportunities for
Huntington's disease.

Requirements:
- Prioritize FDA-approved drugs
- Include mechanistic explanations
- Exclude known treatments
- Assess contraindication risk
```

### Expected Output

```
## TxGNN Drug Repurposing: Huntington's Disease

### Disease Profile
- **Disease:** Huntington's disease (OMIM: 143100)
- **Gene:** HTT (CAG repeat expansion)
- **Existing treatments:** Tetrabenazine (symptomatic only)

### Top Repurposing Candidates

| Rank | Drug | Original Indication | TxGNN Score | Evidence Path |
|------|------|---------------------|-------------|---------------|
| 1 | **Riluzole** | ALS | 0.89 | Drug→BDNF→Neuroprotection→HD |
| 2 | **Memantine** | Alzheimer's | 0.85 | Drug→NMDA→Excitotoxicity→HD |
| 3 | **Cysteamine** | Cystinosis | 0.82 | Drug→Autophagy→HTT clearance |
| 4 | **Lithium** | Bipolar | 0.79 | Drug→GSK3β→Autophagy→HD |
| 5 | **Metformin** | Diabetes | 0.76 | Drug→AMPK→Autophagy→HD |

### Mechanistic Explanations

#### Riluzole (Top candidate)
**Knowledge graph path:**
```
Riluzole → inhibits → Glutamate release
                          ↓
                   ↓ Excitotoxicity
                          ↓
                   Neuroprotection
                          ↓
                   Huntington's disease
```

**Supporting evidence:**
- Reduces glutamate excitotoxicity (known HD mechanism)
- Neuroprotective in other neurodegenerative diseases
- FDA-approved, well-characterized safety

### Contraindication Assessment

| Drug | HD Contraindication Risk | Mechanism |
|------|--------------------------|-----------|
| Riluzole | Low (0.12) | Safe |
| Memantine | Low (0.18) | Safe |
| Haloperidol | **High (0.78)** | Worsens motor symptoms |
| Metoclopramide | **High (0.71)** | D2 antagonism |

### Clinical Trial Suggestions

| Drug | Phase | Design | Primary Endpoint |
|------|-------|--------|------------------|
| Riluzole | II | RCT | UHDRS motor score |
| Cysteamine | II | Open-label | Autophagy biomarkers |
| Metformin | I/II | Safety + efficacy | Safety + UHDRS |
```

---

## References

- **Huang et al. (2023):** "TxGNN: Zero-shot drug repurposing with therapeutic knowledge graphs." *Nature Medicine*
- [TxGNN GitHub](https://github.com/mims-harvard/TxGNN)

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

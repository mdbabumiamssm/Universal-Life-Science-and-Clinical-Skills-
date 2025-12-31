# ESM3 Protein Design Agent

**ID:** `biomedical.protein_science.esm3_design`
**Version:** 1.0.0
**Status:** Production
**Category:** Protein Science / Generative Design

---

## Overview

The **ESM3 Protein Design Agent** leverages EvolutionaryScale's ESM3 protein language model for de novo protein design. ESM3 represents the frontier of generative protein modeling, having demonstrated the ability to generate a novel green fluorescent protein (esmGFP) that differs from any known natural protein by an evolutionary distance equivalent to 500 million years of natural evolution.

This capability enables the design of entirely new protein sequences with specified functions, structures, and properties—a transformative tool for therapeutic protein development, enzyme engineering, and synthetic biology.

---

## Key Capabilities

### 1. Generative Design Modes

| Mode | Description | Application |
|------|-------------|-------------|
| **Unconditional Generation** | Generate diverse protein sequences | Novel scaffold discovery |
| **Structure-Conditioned** | Generate sequences for target fold | Enzyme design |
| **Function-Conditioned** | Generate sequences with target function | Therapeutic proteins |
| **Masked Sequence Infilling** | Complete partial sequences | Protein optimization |

### 2. Protein Property Prediction

- **Fitness prediction:** Estimate functional viability
- **Structure prediction:** ESMFold integration
- **Contact map prediction:** Residue-residue contacts
- **Evolutionary plausibility:** Sequence likelihood scoring

### 3. Sequence Optimization

- **Directed evolution simulation:** In silico evolution
- **Stability optimization:** Thermostability enhancement
- **Solubility improvement:** Expression optimization
- **Immunogenicity reduction:** Therapeutic development

### 4. ESM Model Variants

| Model | Parameters | Best For |
|-------|------------|----------|
| **ESM3** | 98B | Maximum capability, de novo design |
| **ESM-2 (650M)** | 650M | Fast embeddings, transfer learning |
| **ESM-C (600M)** | 600M | Balanced performance/speed |
| **ESMFold** | 15B | Rapid structure prediction |

---

## Usage

### Example Prompt

```text
Design a novel enzyme with beta-lactamase activity.
The enzyme should:
- Be thermostable (Tm > 60°C)
- Express well in E. coli
- Have no sequence identity >30% to known beta-lactamases

Generate 5 candidate sequences with predicted properties.
```

### Expected Output

```
## ESM3 De Novo Enzyme Design

### Design Objective: Novel Beta-Lactamase

#### Generation Parameters
- Model: ESM3-98B
- Conditioning: Structure (TEM-1 fold, PDB 1BTL)
- Constraints: Max 30% identity to known sequences

#### Candidate Sequences

| ID | Length | Identity to TEM-1 | Predicted Tm | pLDDT | Fitness Score |
|----|--------|-------------------|--------------|-------|---------------|
| BLA-01 | 286 | 28.4% | 67.2°C | 88.3 | 0.91 |
| BLA-02 | 291 | 26.1% | 64.8°C | 86.7 | 0.88 |
| BLA-03 | 284 | 29.7% | 71.3°C | 91.2 | 0.93 |
| BLA-04 | 288 | 24.3% | 62.1°C | 84.5 | 0.85 |
| BLA-05 | 290 | 27.8% | 68.9°C | 89.1 | 0.90 |

#### Top Candidate Analysis (BLA-03)
```
>BLA-03_novel_beta_lactamase
MKTLVLFAAVLVFPMWAQAETPSLVKVNKADRTIITPQWKVQLGK...
```

**Key Features:**
- Active site residues (Ser70, Lys73, Ser130, Glu166) preserved
- Novel surface mutations enhance solubility
- Stabilizing core mutations: V33I, A69T, L145M
- Predicted E. coli expression: High (codon-optimized)

#### Validation Recommendations
1. Express in BL21(DE3) with IPTG induction
2. Purify via Ni-NTA (His-tag)
3. Assay activity against nitrocefin (chromogenic substrate)
4. Determine Tm by DSF (Differential Scanning Fluorimetry)
```

### LLM Agent Integration

```python
@tool
def design_protein_esm3(
    design_objective: str,
    template_structure: str = None,
    sequence_constraints: dict = None,
    num_designs: int = 5,
    optimize_for: list[str] = ["stability", "expression"]
) -> str:
    """
    Designs novel proteins using ESM3.

    Args:
        design_objective: Description of desired protein function
        template_structure: PDB ID for structural template (optional)
        sequence_constraints: Constraints like max identity, length
        num_designs: Number of candidate sequences
        optimize_for: Properties to optimize

    Returns:
        Designed sequences with predicted properties
    """
    pass
```

---

## Prerequisites

### Required APIs/Models

| Resource | Purpose | Access |
|----------|---------|--------|
| **ESM3** | Generative protein design | EvolutionaryScale API |
| **ESM-2** | Embeddings and fitness | HuggingFace |
| **ESMFold** | Rapid structure prediction | HuggingFace |
| **UniProt** | Reference sequences | Public API |

### Dependencies

```
torch>=2.0
esm>=3.0  # EvolutionaryScale package
transformers>=4.30
biopython>=1.81
numpy>=1.24
```

### Compute Requirements

| Model | GPU Memory | Inference Time |
|-------|------------|----------------|
| ESM3-98B | 80GB+ (A100) | ~30s/sequence |
| ESM-2-650M | 8GB | ~1s/sequence |
| ESMFold | 16GB | ~5s/structure |

---

## Methodology

### ESM3 Architecture

ESM3 is a multimodal protein language model trained on:
- **Sequence data:** 2.78 billion protein sequences
- **Structure data:** Predicted structures via ESMFold
- **Function annotations:** GO terms, enzyme classifications

### Design Pipeline

```
Design Objective
    ↓
Template Selection (optional)
    ├── Structure conditioning
    └── Function conditioning
    ↓
ESM3 Sequence Generation
    ↓
Fitness Scoring (ESM-2 likelihood)
    ↓
Structure Validation (ESMFold)
    ↓
Property Prediction
    ├── Thermostability (sequence features)
    ├── Solubility (hydrophobicity analysis)
    └── Expression (codon optimization)
    ↓
Diversity Filtering
    ↓
Final Candidates
```

### Sequence Scoring

```python
def calculate_esm_fitness(sequence: str) -> float:
    """
    Calculate evolutionary fitness using ESM-2 pseudo-likelihood.
    Higher scores indicate more evolutionarily plausible sequences.
    """
    import esm

    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()

    # Mask each position and compute likelihood
    log_likelihood = 0.0
    for i in range(len(sequence)):
        masked_seq = sequence[:i] + '<mask>' + sequence[i+1:]
        # ... compute likelihood of true amino acid

    return log_likelihood / len(sequence)
```

---

## Related Skills

- **Virtual Lab Agent:** Integrated design in multi-agent research
- **AlphaFold3 Agent:** Structure validation of designs
- **Antibody Design (MAGE):** Complementary antibody engineering

---

## References

- **Hayes et al. (2024):** "Simulating 500 million years of evolution with a language model." *bioRxiv*
- **Lin et al. (2023):** "Evolutionary-scale prediction of atomic-level protein structure with a language model." *Science*
- [EvolutionaryScale ESM](https://www.evolutionaryscale.ai/)
- [ESM GitHub Repository](https://github.com/facebookresearch/esm)

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

# Microbial Strain Design Agent

**ID:** `biomedical.synthetic_biology.strain_designer`
**Version:** 1.0.0
**Status:** Beta
**Category:** Synthetic Biology / Strain Engineering

---

## Overview

The **Microbial Strain Design Agent** engineers optimized microbial production strains using AI-guided genome-scale metabolic modeling and machine learning. It integrates flux balance analysis with experimental data to predict genetic modifications that maximize product yields while maintaining cell growth.

This agent supports the full Design-Build-Test-Learn cycle, from in silico strain design to experimental validation recommendations, enabling rapid development of industrial bioproduction hosts.

---

## Key Capabilities

### 1. Genome-Scale Modeling

| Analysis | Method | Output |
|----------|--------|--------|
| **FBA** | Flux Balance Analysis | Optimal flux distribution |
| **FVA** | Flux Variability Analysis | Flux ranges |
| **MOMA** | Minimization of Metabolic Adjustment | Knockout predictions |
| **OptKnock** | Bi-level optimization | Growth-coupled designs |

### 2. Strain Modifications

| Type | Description | Example |
|------|-------------|---------|
| **Gene knockout** | Remove competing pathways | ΔpflB, ΔldhA |
| **Gene overexpression** | Increase flux | ↑zwf, ↑ppc |
| **Heterologous expression** | Add new pathways | +crtEBI |
| **Promoter engineering** | Tune expression levels | Ptrc variants |

### 3. Supported Organisms

| Organism | Model | Applications |
|----------|-------|--------------|
| **E. coli** | iML1515 | Chemicals, proteins |
| **S. cerevisiae** | Yeast8 | Ethanol, terpenoids |
| **B. subtilis** | iBsu1103 | Enzymes, vitamins |
| **C. glutamicum** | iCG238 | Amino acids |
| **P. putida** | iJN1463 | Aromatics |

### 4. ML-Guided Design

- **Genotype-phenotype mapping:** Predict titer from genetic design
- **Transfer learning:** Leverage data from related strains
- **Active learning:** Prioritize informative experiments
- **Ensemble models:** Robust predictions with uncertainty

---

## Usage

### Example Prompt

```text
Design an E. coli strain optimized for succinate production.

Requirements:
- Anaerobic fermentation
- Growth-coupled production (essential for selection)
- Theoretical yield > 1.5 mol/mol glucose
- Predict necessary gene knockouts and overexpressions
```

### Expected Output

```
## Succinate Production Strain Design

### Design Specifications
- **Host:** E. coli MG1655
- **Model:** iML1515 (1515 genes, 2719 reactions)
- **Product:** Succinate
- **Condition:** Anaerobic

### Theoretical Analysis

#### Maximum Theoretical Yield
| Condition | Max Yield (mol/mol glucose) |
|-----------|------------------------------|
| Aerobic | 1.71 |
| Anaerobic | 1.71 |
| **Target** | >1.50 |

### OptKnock Strain Design

#### Recommended Knockouts

| Gene | Reaction | Rationale | Impact |
|------|----------|-----------|--------|
| **ldhA** | D-lactate dehydrogenase | Eliminate lactate | +0.32 mol/mol |
| **pflB** | Pyruvate formate lyase | Block acetyl-CoA route | +0.28 mol/mol |
| **adhE** | Alcohol dehydrogenase | Eliminate ethanol | +0.21 mol/mol |
| **pta** | Phosphotransacetylase | Reduce acetate | +0.15 mol/mol |

#### Recommended Overexpressions

| Gene | Enzyme | Rationale | Target Expression |
|------|--------|-----------|-------------------|
| **ppc** | PEP carboxylase | Increase OAA flux | 5x native |
| **pyc** | Pyruvate carboxylase (B. subtilis) | CO2 fixation | Heterologous |
| **mdh** | Malate dehydrogenase | Drive to malate | 3x native |
| **fumB** | Fumarase (anaerobic) | Malate → fumarate | 2x native |
| **frd** | Fumarate reductase | Fumarate → succinate | 3x native |

### Flux Distribution (Optimized Strain)

```
Glucose (100)
    ↓
PEP (180) ←──┐
    ↓         │
Pyruvate (80) │ ppc (120)
    ↓         │
OAA ←─────────┘
    ↓ mdh
Malate
    ↓ fumB
Fumarate
    ↓ frd
Succinate (152)
```

**Predicted yield:** 1.52 mol/mol glucose

### Growth-Coupling Analysis

| Modification | Growth Rate | Succinate Flux | Coupled? |
|--------------|-------------|----------------|----------|
| Wild-type | 0.45 h⁻¹ | 0.12 mol/mol | No |
| ΔldhA | 0.38 h⁻¹ | 0.35 mol/mol | Partial |
| ΔldhA ΔpflB | 0.29 h⁻¹ | 0.78 mol/mol | Partial |
| ΔldhA ΔpflB ΔadhE | 0.18 h⁻¹ | 1.42 mol/mol | **Yes** |
| + Overexpressions | 0.22 h⁻¹ | 1.52 mol/mol | **Yes** |

### Predicted Phenotype

| Metric | Wild-type | Engineered | Change |
|--------|-----------|------------|--------|
| Growth rate | 0.45 h⁻¹ | 0.22 h⁻¹ | -51% |
| Succinate yield | 0.12 mol/mol | 1.52 mol/mol | +12.7x |
| Byproduct (lactate) | 0.8 mol/mol | 0 | -100% |
| Byproduct (ethanol) | 0.5 mol/mol | 0 | -100% |

### Strain Construction Plan

#### Phase 1: Sequential Knockouts
```
MG1655 → ΔldhA → ΔpflB → ΔadhE → Δpta
```
**Method:** λ Red recombination with FRT-flanked markers

#### Phase 2: Pathway Overexpression
```
pSucc1: Ptrc-ppc-pyc (medium copy, ~20/cell)
pSucc2: Ptrc-mdh-fumB-frd (low copy, ~5/cell)
```

### Experimental Validation

| Test | Metric | Expected |
|------|--------|----------|
| Anaerobic shake flask | Succinate titer | 15-20 g/L |
| Fed-batch fermentation | Succinate titer | 80-100 g/L |
| Yield | mol/mol glucose | 1.4-1.6 |
| Productivity | g/L/h | 1.5-2.0 |
```

### LLM Agent Integration

```python
@tool
def design_production_strain(
    product: str,
    host: str = "e_coli",
    condition: str = "aerobic",
    growth_coupled: bool = True,
    max_knockouts: int = 6
) -> str:
    """
    Designs optimized production strain using genome-scale modeling.

    Args:
        product: Target product (KEGG ID or name)
        host: Host organism
        condition: aerobic or anaerobic
        growth_coupled: Require growth-coupled production
        max_knockouts: Maximum gene knockouts allowed

    Returns:
        Strain design with modifications and predictions
    """
    pass
```

---

## Prerequisites

### Required Tools

| Resource | Purpose | Access |
|----------|---------|--------|
| **COBRApy** | FBA modeling | Open source |
| **OptKnock** | Strain design | COBRA toolbox |
| **BiGG Models** | GEM database | Public |
| **Escher** | Flux visualization | Open source |

### Dependencies

```
cobra>=0.26
optlang>=1.5
numpy>=1.24
pandas>=1.5
escher>=1.7
```

---

## References

- **Burgard et al. (2003):** "OptKnock: A bilevel programming framework for identifying gene knockout strategies." *Biotechnology and Bioengineering*
- **Feist et al. (2007):** "A genome-scale metabolic reconstruction for E. coli K-12." *Molecular Systems Biology*
- [BiGG Models Database](http://bigg.ucsd.edu/)

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

# RNA Structure Prediction Agent

**ID:** `biomedical.drug_discovery.rna_structure`
**Version:** 1.0.0
**Status:** Beta
**Category:** Drug Discovery / RNA Therapeutics

---

## Overview

The **RNA Structure Prediction Agent** predicts secondary and tertiary RNA structures critical for therapeutic RNA design. RNA structure determines stability, translation efficiency, and target accessibility—all essential parameters for mRNA vaccines, siRNAs, ribozymes, and aptamers.

This agent integrates thermodynamic models (ViennaRNA), comparative sequence analysis, and emerging deep learning approaches to provide comprehensive structural characterization of therapeutic RNA molecules.

---

## Key Capabilities

### 1. Secondary Structure Prediction

| Method | Algorithm | Best For |
|--------|-----------|----------|
| **MFE** | Minimum Free Energy | Single optimal structure |
| **Partition Function** | Base pair probabilities | Ensemble analysis |
| **Centroid** | Ensemble centroid | Representative structure |
| **MEA** | Maximum Expected Accuracy | Balanced prediction |
| **Comparative** | Covariance analysis | Conserved structures |

### 2. Structure Types

- **mRNA 5' UTR:** Translation initiation context
- **mRNA CDS:** Ribosome accessibility
- **siRNA duplex:** Thermodynamic asymmetry
- **Aptamers:** Ligand binding structures
- **Ribozymes:** Catalytic conformations
- **lncRNA:** Functional domains

### 3. Analysis Features

| Analysis | Output | Application |
|----------|--------|-------------|
| **MFE Calculation** | ΔG (kcal/mol) | Stability assessment |
| **Pair Probabilities** | Base pair matrix | Structure uncertainty |
| **Accessibility** | Single-stranded regions | siRNA target sites |
| **Local Folding** | Windowed analysis | Long RNA structure |
| **Pseudoknot Detection** | Non-nested base pairs | Complex structures |

### 4. 3D Structure Prediction

- **Homology modeling:** Template-based prediction
- **De novo prediction:** Physics-based (FARFAR2)
- **ML approaches:** Deep learning (RhoFold)

---

## Usage

### Example Prompt

```text
Predict the secondary structure of this mRNA 5' UTR:
GGGAAATAAGAGAGAAAAGAAGAGTAAGAAGAAATATAAGAGCCACCATG

Analyze:
1. Minimum free energy structure
2. Ensemble diversity (partition function)
3. Kozak sequence accessibility
4. Any stable hairpins that might inhibit translation
```

### Expected Output

```
## RNA Structure Analysis: 5' UTR

### Sequence (50 nt)
```
GGGAAATAAGAGAGAAAAGAAGAGTAAGAAGAAATATAAGAGCCACCATG
```

### Minimum Free Energy (MFE) Structure

```
GGGAAATAAGAGAGAAAAGAAGAGTAAGAAGAAATATAAGAGCCACCATG
...(((....))).........((((.....))))................
```

**ΔG = -5.20 kcal/mol**

#### Structure Elements
| Element | Position | ΔG (kcal/mol) | Notes |
|---------|----------|---------------|-------|
| Hairpin 1 | 4-15 | -2.8 | Small stem-loop |
| Hairpin 2 | 26-42 | -2.4 | Small stem-loop |

### Partition Function Analysis

- **Ensemble Diversity:** 8.34 (moderate diversity)
- **MFE Probability:** 12.3%
- **Centroid Structure:** Similar to MFE

### Base Pair Probability Matrix
```
Position  4-15:  Probability 0.78 (well-defined)
Position 26-42:  Probability 0.65 (moderate confidence)
```

### Translation Initiation Context

#### Kozak Sequence Analysis
```
Position:    -6  -5  -4  -3  -2  -1  +1  +2  +3  +4
Sequence:     G   C   C   A   C   C  [A   T   G]  ...
Consensus:    G   C   C  (A/G) C   C  [A   T   G]  G
```

**Kozak Strength:** Strong (matches 5/6 consensus positions)

#### Accessibility at Start Codon
- **Positions 44-50 (CCACCATG):** 95% single-stranded
- **Translation initiation:** Favorable (minimal structure)

### Recommendations

| Aspect | Assessment | Action |
|--------|------------|--------|
| **Overall stability** | Moderate (-5.2 kcal/mol) | Acceptable |
| **5' cap accessibility** | Good (3 nt unstructured) | No change needed |
| **AUG accessibility** | Excellent | No change needed |
| **Upstream structures** | Minor hairpins | Consider if expression low |

### Output Files
- `structure.dot` - Dot-bracket notation
- `structure.ct` - Connect table format
- `bpp_matrix.csv` - Base pair probabilities
- `structure.svg` - Visualization
```

### LLM Agent Integration

```python
@tool
def predict_rna_structure(
    sequence: str,
    method: str = "mfe",
    temperature: float = 37.0,
    analyze_accessibility: bool = True,
    detect_pseudoknots: bool = False
) -> str:
    """
    Predicts RNA secondary structure.

    Args:
        sequence: RNA sequence (U or T accepted)
        method: mfe, partition, centroid, or mea
        temperature: Folding temperature in Celsius
        analyze_accessibility: Calculate accessibility profile
        detect_pseudoknots: Enable pseudoknot prediction

    Returns:
        Structure prediction with analysis
    """
    pass


@tool
def analyze_sirna_structure(
    sense_strand: str,
    antisense_strand: str
) -> str:
    """
    Analyzes siRNA duplex structure and thermodynamics.

    Args:
        sense_strand: 5'→3' sense strand sequence
        antisense_strand: 5'→3' antisense strand sequence

    Returns:
        Duplex analysis including asymmetry and stability
    """
    pass
```

---

## Prerequisites

### Required Tools

| Resource | Purpose | Access |
|----------|---------|--------|
| **ViennaRNA** | Structure prediction | Local installation |
| **RNAstructure** | Alternative algorithms | Mathews Lab |
| **Contrafold** | ML-based prediction | Stanford |
| **NUPACK** | Multi-strand analysis | Caltech |

### Dependencies

```
ViennaRNA>=2.5.1
forgi>=2.0  # RNA structure handling
matplotlib>=3.7
numpy>=1.24
pandas>=1.5
```

---

## Methodology

### Thermodynamic Model

ViennaRNA uses the Turner 2004 nearest-neighbor model:

```
ΔG_total = Σ ΔG_stack + Σ ΔG_loop + ΔG_external

Where:
- ΔG_stack: Base pair stacking energies
- ΔG_loop: Hairpin, internal, bulge, multiloop penalties
- ΔG_external: Dangling ends, terminal mismatches
```

### Partition Function

```python
import RNA

def analyze_structure_ensemble(sequence: str) -> dict:
    """
    Comprehensive structure analysis using partition function.
    """
    # Fold compound
    fc = RNA.fold_compound(sequence)

    # MFE structure
    mfe_struct, mfe_energy = fc.mfe()

    # Partition function
    fc.pf()

    # Centroid structure
    centroid_struct, centroid_dist = fc.centroid()

    # Base pair probabilities
    bpp = fc.bpp()

    # Ensemble diversity
    diversity = fc.mean_bp_distance()

    return {
        'mfe_structure': mfe_struct,
        'mfe_energy': mfe_energy,
        'centroid_structure': centroid_struct,
        'ensemble_diversity': diversity,
        'base_pair_probabilities': bpp
    }
```

### Accessibility Calculation

```python
def calculate_accessibility(
    sequence: str,
    window_size: int = 20
) -> np.ndarray:
    """
    Calculate single-stranded accessibility for each position.
    """
    import RNA

    accessibility = np.zeros(len(sequence))

    for i in range(len(sequence)):
        # Calculate probability of being unpaired
        fc = RNA.fold_compound(sequence)
        fc.pf()

        # Get unpaired probability
        up = fc.pr_unpaired(i + 1, window_size)
        accessibility[i] = up

    return accessibility
```

---

## Related Skills

- **mRNA Design Agent:** Optimize mRNA structure for expression
- **siRNA/ASO Design Agent:** Target accessibility analysis
- **CRISPR Design Agent:** crRNA structure optimization

---

## References

- **Lorenz et al. (2011):** "ViennaRNA Package 2.0." *Algorithms for Molecular Biology*
- **Mathews et al. (2004):** "Incorporating chemical modification constraints into RNA secondary structure prediction." *PNAS*
- [ViennaRNA Web Server](http://rna.tbi.univie.ac.at/)
- [NUPACK Web Application](http://www.nupack.org/)

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

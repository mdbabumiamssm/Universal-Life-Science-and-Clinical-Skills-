# Epitope Prediction Agent

**ID:** `biomedical.immunology.epitope_prediction`
**Version:** 1.0.0
**Status:** Production
**Category:** Immunology / Vaccine Design

---

## Overview

The **Epitope Prediction Agent** identifies immunogenic peptide sequences (epitopes) that can elicit B-cell and T-cell immune responses. Accurate epitope prediction is fundamental for vaccine design, immunotherapy development, and understanding immune responses to pathogens and tumors.

This agent integrates deep learning models including CNNs, Transformers, and specialized tools like NetMHCpan to predict MHC binding, immunogenicity, and population coverage for epitope-based vaccine design.

---

## Key Capabilities

### 1. T-Cell Epitope Prediction

| Type | MHC Class | Cell Type | Output |
|------|-----------|-----------|--------|
| **MHC-I binding** | HLA-A, B, C | CD8+ CTL | 8-11 mer peptides |
| **MHC-II binding** | HLA-DR, DP, DQ | CD4+ helper | 13-25 mer peptides |
| **Immunogenicity** | - | Both | Processing + presentation |
| **TCR binding** | - | - | T-cell recognition |

### 2. B-Cell Epitope Prediction

- **Linear epitopes:** Continuous peptide sequences
- **Conformational epitopes:** 3D surface patches
- **Antigenicity scoring:** BepiPred, ABCpred methods
- **Accessibility analysis:** Surface exposure prediction

### 3. Population Coverage Analysis

| Population | HLA Supertype | Coverage Target |
|------------|---------------|-----------------|
| **Global** | Multiple | >90% |
| **Caucasian** | A02, B07, DR1 | >95% |
| **Asian** | A24, B46, DR4 | >90% |
| **African** | A30, B58, DR3 | >85% |

### 4. Vaccine Design Features

- **Multi-epitope assembly:** Linker optimization (AAY, GPGPG)
- **Adjuvant selection:** TLR agonist recommendations
- **Delivery optimization:** mRNA, peptide, viral vector
- **Conservation analysis:** Cross-strain coverage

---

## Usage

### Example Prompt

```text
Predict T-cell epitopes for SARS-CoV-2 Spike protein (UniProt: P0DTC2).
Requirements:
- HLA-A*02:01 restricted CD8+ epitopes (9-mers)
- Immunogenicity score > 0.5
- Conservation across variants (Alpha, Delta, Omicron)
- Top 10 candidates for vaccine development
```

### Expected Output

```
## T-Cell Epitope Prediction: SARS-CoV-2 Spike Protein

### Target Analysis
- **Protein:** Spike glycoprotein (P0DTC2)
- **Length:** 1,273 amino acids
- **HLA restriction:** HLA-A*02:01

### Top 10 MHC-I Epitope Candidates

| Rank | Position | Sequence | NetMHCpan Score | Immunogenicity | Conservation |
|------|----------|----------|-----------------|----------------|--------------|
| 1 | 269-277 | YLQPRTFLL | 0.98 | 0.78 | 100% |
| 2 | 417-425 | KIADYNYKL | 0.95 | 0.72 | 98.5% |
| 3 | 1000-1008 | RLITGRLQSL | 0.94 | 0.69 | 100% |
| 4 | 865-873 | LTDEMIAQY | 0.93 | 0.71 | 99.2% |
| 5 | 3-11 | FVFLVLLPL | 0.92 | 0.65 | 100% |
| 6 | 508-516 | YRVVVLSFEL | 0.91 | 0.67 | 97.8% |
| 7 | 976-984 | VLNDILSRL | 0.90 | 0.64 | 100% |
| 8 | 378-386 | YKCFGGFNL | 0.89 | 0.68 | 96.5% |
| 9 | 691-699 | SIIAYTMSL | 0.88 | 0.62 | 100% |
| 10 | 1136-1144 | FIAGLIAIV | 0.87 | 0.66 | 100% |

### Lead Candidate Analysis: YLQPRTFLL (Position 269-277)

#### MHC Binding
- **HLA-A*02:01 IC50:** 12.4 nM (strong binder)
- **Percentile rank:** 0.02% (top tier)
- **Processing score:** 0.89 (efficient TAP transport)

#### Immunogenicity
- **DeepImmuno score:** 0.78
- **VaxiJen score:** 0.72 (probable antigen)
- **Previous reports:** Confirmed immunogenic in COVID patients

#### Variant Conservation
| Variant | Sequence | Match |
|---------|----------|-------|
| Wuhan | YLQPRTFLL | 100% |
| Alpha | YLQPRTFLL | 100% |
| Delta | YLQPRTFLL | 100% |
| Omicron | YLQPRTFLL | 100% |

#### Structural Location
- **Domain:** S1 subunit, near RBD
- **Accessibility:** Surface exposed
- **Function:** Adjacent to ACE2 binding interface

### Population Coverage Analysis

| Epitope Combination | Global | Caucasian | Asian | African |
|--------------------|--------|-----------|-------|---------|
| Top 5 epitopes | 67.3% | 72.1% | 64.8% | 58.9% |
| Top 10 epitopes | 89.2% | 92.4% | 86.7% | 81.3% |
| All 10 + DR epitopes | 97.8% | 98.9% | 96.5% | 94.2% |

### Vaccine Design Recommendations

#### Multi-Epitope Construct
```
YLQPRTFLL-AAY-KIADYNYKL-AAY-RLITGRLQSL-GPGPG-[CD4 epitopes]
```

#### Recommended Formulation
- **Platform:** mRNA (m1Ψ-modified)
- **Adjuvant:** LNP delivery (self-adjuvanting)
- **Boost:** 3-4 weeks post-prime
```

### LLM Agent Integration

```python
@tool
def predict_tcell_epitopes(
    protein_sequence: str,
    hla_alleles: list[str] = ["HLA-A*02:01"],
    epitope_length: int = 9,
    num_candidates: int = 10,
    check_conservation: bool = True
) -> str:
    """
    Predicts T-cell epitopes from protein sequence.

    Args:
        protein_sequence: Amino acid sequence or UniProt ID
        hla_alleles: List of HLA alleles for prediction
        epitope_length: Peptide length (8-11 for MHC-I, 13-25 for MHC-II)
        num_candidates: Number of top candidates to return
        check_conservation: Analyze conservation across strains

    Returns:
        Ranked epitope candidates with scores
    """
    pass


@tool
def predict_bcell_epitopes(
    protein_sequence: str,
    structure_pdb: str = None,
    epitope_type: str = "linear"
) -> str:
    """
    Predicts B-cell epitopes (linear or conformational).

    Args:
        protein_sequence: Amino acid sequence
        structure_pdb: PDB ID for conformational epitopes
        epitope_type: linear or conformational

    Returns:
        B-cell epitope predictions with antigenicity scores
    """
    pass
```

---

## Prerequisites

### Required Tools/APIs

| Resource | Purpose | Access |
|----------|---------|--------|
| **NetMHCpan 4.1** | MHC-I binding prediction | DTU API |
| **NetMHCIIpan** | MHC-II binding prediction | DTU API |
| **DeepImmuno** | Immunogenicity prediction | Web service |
| **BepiPred** | B-cell epitope prediction | IEDB API |
| **IEDB** | Immune epitope database | Public API |

### Dependencies

```
pandas>=1.5
numpy>=1.24
biopython>=1.81
requests>=2.28
torch>=2.0  # For deep learning models
```

---

## Methodology

### MHC-I Binding Prediction

NetMHCpan uses a pan-specific neural network trained on 180,000+ binding measurements:

```python
def predict_mhc_binding(
    peptide: str,
    hla_allele: str
) -> dict:
    """
    Predict MHC-I binding using NetMHCpan.
    """
    import requests

    # NetMHCpan API call
    response = requests.post(
        "https://services.healthtech.dtu.dk/api/netmhcpan/4.1",
        data={
            "method": "netmhcpan",
            "sequence": peptide,
            "allele": hla_allele
        }
    )

    result = response.json()

    return {
        'peptide': peptide,
        'allele': hla_allele,
        'affinity_nM': result['affinity'],
        'percentile_rank': result['rank'],
        'binding_level': classify_binding(result['affinity'])
    }


def classify_binding(affinity_nM: float) -> str:
    """Classify binding strength."""
    if affinity_nM < 50:
        return "Strong Binder"
    elif affinity_nM < 500:
        return "Weak Binder"
    else:
        return "Non-Binder"
```

### Immunogenicity Scoring

```python
def predict_immunogenicity(
    peptide: str,
    mhc_allele: str
) -> float:
    """
    Predict immunogenicity using DeepImmuno.

    Combines MHC binding with TCR recognition probability.
    """
    # Feature extraction
    features = {
        'mhc_binding': predict_mhc_binding(peptide, mhc_allele),
        'aa_composition': calculate_aa_features(peptide),
        'hydrophobicity': calculate_hydrophobicity(peptide),
        'position_features': encode_positions(peptide)
    }

    # DeepImmuno CNN model
    score = deep_immuno_model.predict(features)

    return score
```

---

## Related Skills

- **Neoantigen Vaccine Agent:** Personalized cancer vaccine design
- **mRNA Design Agent:** mRNA vaccine sequence optimization
- **Multi-Epitope Designer:** Vaccine construct assembly

---

## References

- **Reynisson et al. (2020):** "NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation." *Nucleic Acids Research*
- **Li et al. (2021):** "DeepImmuno: deep learning-empowered prediction and generation of immunogenic peptides." *Briefings in Bioinformatics*
- [IEDB Analysis Resource](http://tools.iedb.org/)
- [NetMHCpan Server](https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1)

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

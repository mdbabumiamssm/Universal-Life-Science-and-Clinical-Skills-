# Multi-Epitope Vaccine Designer

**ID:** `biomedical.immunology.multi_epitope_designer`
**Version:** 1.0.0
**Status:** Production
**Category:** Immunology / Vaccine Design

---

## Overview

The **Multi-Epitope Vaccine Designer** constructs optimized vaccine candidates by combining multiple B-cell and T-cell epitopes into a single immunogenic construct. This approach enables broad immune coverage across HLA types, pathogen variants, and immune response arms (humoral + cellular).

This agent handles epitope selection, linker optimization, adjuvant integration, and construct validation to generate production-ready vaccine sequences for mRNA, peptide, or recombinant protein platforms.

---

## Key Capabilities

### 1. Epitope Assembly

| Component | Function | Optimization |
|-----------|----------|--------------|
| **CTL epitopes** | CD8+ T-cell response | AAY/AYY linkers |
| **HTL epitopes** | CD4+ helper response | GPGPG linkers |
| **B-cell epitopes** | Antibody induction | KK linkers |
| **Adjuvant** | Immune activation | N-terminal fusion |

### 2. Linker Library

| Linker | Sequence | Use Case |
|--------|----------|----------|
| **AAY** | AAY | Between MHC-I epitopes |
| **GPGPG** | GPGPG | Between MHC-II epitopes |
| **KK** | KK | Between B-cell epitopes |
| **EAAAK** | EAAAK | Rigid α-helical spacer |
| **GGGGS** | (GGGGS)n | Flexible spacer |

### 3. Adjuvant Options

| Adjuvant | Mechanism | Receptor |
|----------|-----------|----------|
| **β-defensin** | Innate immunity | TLR1/2 |
| **PADRE** | Universal T-helper | MHC-II |
| **Flagellin** | Strong innate activation | TLR5 |
| **RS09** | TLR4 agonist | TLR4 |
| **hBD-3** | Human β-defensin 3 | Multiple |

### 4. Validation Tools

- **Antigenicity prediction:** VaxiJen, ANTIGENpro
- **Allergenicity check:** AllerTOP, AllergenFP
- **Solubility prediction:** SOLpro
- **Physicochemical analysis:** ProtParam
- **3D structure prediction:** AlphaFold2/ESMFold

---

## Usage

### Example Prompt

```text
Design a multi-epitope vaccine against Dengue virus (all 4 serotypes).

Requirements:
- Include 10 CTL epitopes (HLA supertype coverage)
- Include 5 HTL epitopes (CD4+ helper)
- Include 3 B-cell epitopes (neutralizing antibody targets)
- Use β-defensin adjuvant at N-terminus
- Optimize for E. coli expression
- Cross-serotype coverage >80%

Provide the final construct sequence with all analyses.
```

### Expected Output

```
## Multi-Epitope Dengue Vaccine Design

### Design Summary
| Component | Count | Coverage |
|-----------|-------|----------|
| CTL epitopes | 10 | 95.4% global HLA |
| HTL epitopes | 5 | 89.2% global HLA |
| B-cell epitopes | 3 | All 4 serotypes |
| Adjuvant | β-defensin | N-terminal |
| **Total length** | 498 AA | - |

### Construct Architecture

```
[β-defensin]-EAAAK-[CTL1]-AAY-[CTL2]-AAY-...-[CTL10]-GPGPG-[HTL1]-GPGPG-...-[HTL5]-KK-[BCL1]-KK-[BCL2]-KK-[BCL3]-6xHis
```

### Selected Epitopes

#### CTL Epitopes (MHC-I)

| # | Sequence | HLA | Serotype | Protein | Conservation |
|---|----------|-----|----------|---------|--------------|
| 1 | GTSGSPIDK | A*01 | 1-4 | NS3 | 97% |
| 2 | KPRVIDRKL | A*02 | 1-4 | NS5 | 98% |
| 3 | IPMTQTIIM | A*02 | 1-4 | NS4B | 95% |
| 4 | VPMVTQMAM | A*24 | 1-4 | E | 94% |
| 5 | LPAIVREAI | B*07 | 1-4 | NS3 | 96% |
| 6 | GPGHEEPIPL | B*35 | 1-4 | NS5 | 99% |
| 7 | TTMIRNQEL | B*44 | 1-4 | NS2A | 93% |
| 8 | GLPIRYQTPA | A*03 | 1-4 | NS1 | 95% |
| 9 | LPRDLGWLAL | B*08 | 1-4 | NS4A | 97% |
| 10 | YGMEIRPLK | A*11 | 1-4 | NS5 | 98% |

**Combined HLA coverage:** 95.4% (global population)

#### HTL Epitopes (MHC-II)

| # | Sequence | HLA-DR | Serotype | Protein |
|---|----------|--------|----------|---------|
| 1 | WGNGCGLFGKGGIVTC | DR1,3,4 | 1-4 | E |
| 2 | PKQYAGGVFKDTRPHT | DR1,4,7 | 1-4 | NS3 |
| 3 | LRWGMSCKDTLKLFPM | DR3,4,15 | 1-4 | NS5 |
| 4 | STLPETTVVRRRGGRT | DR1,7,11 | 1-4 | C |
| 5 | LGKLLNRKPNFDVFV | DR4,7,9 | 1-4 | NS1 |

#### B-Cell Epitopes (Linear)

| # | Sequence | Domain | Neutralizing | Serotype |
|---|----------|--------|--------------|----------|
| 1 | GCWVKQEGLFGKGGIV | E-DII | Yes (EDE1 target) | 1-4 |
| 2 | WDFGSIGGVFTSVGK | E-DIII | Yes | 1-4 |
| 3 | GCFGKGSIWYTKE | E-FL | Cross-reactive | 1-4 |

### Final Construct Sequence

```
>Multi_Epitope_DENV_Vaccine (498 AA)
MRIHYLLFALLFLFLVPVPGHGGIINTLQKYYCRVRGGRCAVLSCLPK
EEQIGKCSTRGRKCCRRKKEAAAKGTSGSPIDK AAY KPRVIDRKL
AAY IPMTQTIIM AAY VPMVTQMAM AAY LPAIVREAI AAY
GPGHEEPIPL AAY TTMIRNQEL AAY GLPIRYQTPA AAY
LPRDLGWLAL AAY YGMEIRPLK GPGPG WGNGCGLFGKGGIVTC
GPGPG PKQYAGGVFKDTRPHT GPGPG LRWGMSCKDTLKLFPM
GPGPG STLPETTVVRRRGGRT GPGPG LGKLLNRKPNFDVFV KK
GCWVKQEGLFGKGGIV KK WDFGSIGGVFTSVGK KK
GCFGKGSIWYTKE HHHHHH
```

### Construct Validation

#### Antigenicity Analysis
| Tool | Score | Threshold | Result |
|------|-------|-----------|--------|
| VaxiJen | 0.72 | 0.4 | ✓ Probable antigen |
| ANTIGENpro | 0.81 | 0.5 | ✓ Antigenic |

#### Allergenicity Analysis
| Tool | Prediction | Confidence |
|------|------------|------------|
| AllerTOP | Non-allergen | High |
| AllergenFP | Non-allergen | Medium |

#### Physicochemical Properties
| Property | Value | Interpretation |
|----------|-------|----------------|
| **MW** | 54.2 kDa | Suitable for expression |
| **pI** | 9.8 | Basic protein |
| **Instability Index** | 32.4 | Stable (<40) |
| **GRAVY** | -0.42 | Hydrophilic |
| **Aliphatic Index** | 78.6 | Good thermostability |
| **Half-life (E. coli)** | >10 hours | Acceptable |

#### Solubility Prediction
| Tool | Probability | Verdict |
|------|-------------|---------|
| SOLpro | 0.74 | Soluble |
| Protein-Sol | 0.68 | Likely soluble |

#### Population Coverage
| Population | CTL | HTL | Combined |
|------------|-----|-----|----------|
| Global | 95.4% | 89.2% | 99.1% |
| Asia | 94.8% | 87.6% | 98.7% |
| Africa | 91.2% | 85.4% | 97.3% |
| Americas | 96.1% | 90.8% | 99.4% |

### 3D Structure Prediction (AlphaFold)

- **Overall pLDDT:** 72.4 (good confidence)
- **Well-structured regions:** Adjuvant domain, HTL epitopes
- **Flexible regions:** Linkers (expected)
- **No steric clashes:** Validated

### Expression Recommendations

| Platform | Vector | Yield Estimate |
|----------|--------|----------------|
| **E. coli** | pET28a | 50-100 mg/L |
| **Insect cells** | Baculovirus | 20-50 mg/L |
| **Mammalian** | HEK293 | 10-30 mg/L |

**Recommended codon optimization:** E. coli (CAI target: 0.85)
```

### LLM Agent Integration

```python
@tool
def design_multi_epitope_vaccine(
    ctl_epitopes: list[dict],
    htl_epitopes: list[dict],
    bcell_epitopes: list[dict],
    adjuvant: str = "beta_defensin",
    expression_host: str = "e_coli"
) -> str:
    """
    Designs multi-epitope vaccine construct.

    Args:
        ctl_epitopes: List of CTL epitopes [{seq, hla, source}]
        htl_epitopes: List of HTL epitopes [{seq, hla, source}]
        bcell_epitopes: List of B-cell epitopes [{seq, location}]
        adjuvant: Adjuvant type (beta_defensin, PADRE, flagellin)
        expression_host: Target expression system

    Returns:
        Complete vaccine construct with validation analyses
    """
    pass


@tool
def validate_vaccine_construct(
    sequence: str,
    check_antigenicity: bool = True,
    check_allergenicity: bool = True,
    predict_structure: bool = False
) -> str:
    """
    Validates multi-epitope vaccine construct.

    Args:
        sequence: Amino acid sequence of vaccine construct
        check_antigenicity: Run antigenicity prediction
        check_allergenicity: Run allergenicity prediction
        predict_structure: Predict 3D structure (slower)

    Returns:
        Comprehensive validation report
    """
    pass
```

---

## Prerequisites

### Required Tools/APIs

| Resource | Purpose | Access |
|----------|---------|--------|
| **VaxiJen** | Antigenicity prediction | Web service |
| **AllerTOP** | Allergenicity prediction | Web service |
| **ProtParam** | Physicochemical properties | ExPASy |
| **SOLpro** | Solubility prediction | Web service |
| **IEDB** | Population coverage | API |

### Dependencies

```
biopython>=1.81
requests>=2.28
pandas>=1.5
numpy>=1.24
```

---

## Methodology

### Construct Assembly Algorithm

```python
def assemble_multi_epitope_vaccine(
    ctl_epitopes: list[str],
    htl_epitopes: list[str],
    bcell_epitopes: list[str],
    adjuvant: str
) -> str:
    """
    Assemble multi-epitope vaccine with optimized linkers.
    """
    # Define linkers
    linkers = {
        'ctl': 'AAY',      # Between CTL epitopes
        'htl': 'GPGPG',    # Between HTL epitopes
        'bcell': 'KK',     # Between B-cell epitopes
        'domain': 'EAAAK'  # Between domains
    }

    # Adjuvant sequences
    adjuvants = {
        'beta_defensin': 'MRIHYLLFALLFLFLVPVPGHGGIINTLQKYYCRVRGGRCAVLSCLPKEEQIGKCSTRGRKCCRRKK',
        'PADRE': 'AKFVAAWTLKAAA',
        'flagellin': 'MAQVINTNSLSLLTQNNLNK...'  # Full sequence
    }

    # Assemble construct
    construct = adjuvants[adjuvant]
    construct += linkers['domain']

    # Add CTL epitopes
    construct += linkers['ctl'].join(ctl_epitopes)
    construct += linkers['domain']

    # Add HTL epitopes
    construct += linkers['htl'].join(htl_epitopes)
    construct += linkers['domain']

    # Add B-cell epitopes
    construct += linkers['bcell'].join(bcell_epitopes)

    # Add His-tag for purification
    construct += 'HHHHHH'

    return construct
```

---

## Related Skills

- **Epitope Prediction Agent:** Epitope identification
- **Neoantigen Vaccine Agent:** Cancer vaccine design
- **mRNA Design Agent:** mRNA vaccine optimization

---

## References

- **Dar et al. (2019):** "Designing of a multi-epitope vaccine against Zika virus." *Journal of Biomolecular Structure and Dynamics*
- **Ali et al. (2021):** "Computational design of multi-epitope vaccine against SARS-CoV-2." *Computers in Biology and Medicine*
- [IEDB Population Coverage Tool](http://tools.iedb.org/population/)
- [VaxiJen Server](http://www.ddg-pharmfac.net/vaxijen/)

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

# Neoantigen Vaccine Design Agent

**ID:** `biomedical.immunology.neoantigen_vaccine`
**Version:** 1.0.0
**Status:** Production
**Category:** Immunology / Cancer Immunotherapy

---

## Overview

The **Neoantigen Vaccine Design Agent** designs personalized cancer vaccines based on tumor-specific mutations (neoantigens). Neoantigens are peptides arising from tumor-specific mutations that are absent in normal tissues, making them ideal targets for cancer immunotherapy with minimal off-target effects.

This agent integrates genomic variant calling, MHC binding prediction, and immunogenicity scoring to identify the most promising neoantigens for personalized cancer vaccine development.

---

## Key Capabilities

### 1. Neoantigen Identification Pipeline

| Step | Input | Output |
|------|-------|--------|
| **Variant calling** | Tumor/Normal WES/WGS | Somatic mutations |
| **HLA typing** | WES/RNA-seq | Patient HLA alleles |
| **Peptide generation** | Mutations + Reference | Mutant peptides |
| **MHC binding** | Peptides + HLA | Binding predictions |
| **Immunogenicity** | Candidate neoantigens | Prioritized list |

### 2. Mutation Types Supported

| Type | Example | Neoantigen Potential |
|------|---------|---------------------|
| **SNV** | V600E | Single AA change |
| **Indel** | Frameshift | Novel peptide sequence |
| **Fusion** | BCR-ABL | Junction peptides |
| **Splice variant** | Exon skipping | Novel exon junctions |

### 3. Prioritization Criteria

- **MHC-I binding affinity:** IC50 < 500 nM
- **Differential binding:** Mutant >> Wild-type
- **Expression level:** RNA expression > threshold
- **Clonality:** Clonal mutations preferred
- **Immunogenicity:** TCR recognition probability

### 4. Vaccine Platforms

| Platform | Advantages | Clinical Stage |
|----------|------------|----------------|
| **mRNA** | Fast manufacturing, strong responses | Phase II/III |
| **Peptide + Adjuvant** | Stable, defined antigens | Phase II |
| **DC vaccines** | Autologous, personalized | Phase I/II |
| **Viral vectors** | Strong CD8+ responses | Phase I |

---

## Usage

### Example Prompt

```text
Design a personalized neoantigen vaccine for a melanoma patient.

Input data:
- Tumor WES VCF file with somatic mutations
- Patient HLA type: HLA-A*02:01, HLA-A*24:02, HLA-B*07:02, HLA-B*44:02
- RNA-seq expression data available

Identify the top 20 neoantigen candidates ranked by:
1. MHC binding affinity
2. Differential binding (mutant vs wild-type)
3. Expression level
4. Predicted immunogenicity

Design an mRNA vaccine construct with the selected neoantigens.
```

### Expected Output

```
## Personalized Neoantigen Vaccine Design

### Patient Profile
- **Tumor type:** Melanoma
- **HLA-I alleles:** A*02:01, A*24:02, B*07:02, B*44:02
- **Mutation burden:** 312 somatic mutations
- **Expressed mutations:** 187 (RNA-seq filtered)

### Neoantigen Identification Summary

| Stage | Count |
|-------|-------|
| Total somatic mutations | 312 |
| Nonsynonymous | 198 |
| Expressed (TPM > 1) | 187 |
| MHC binders (IC50 < 500nM) | 89 |
| Differential binders (mut/wt > 5x) | 42 |
| Immunogenic (score > 0.5) | 28 |

### Top 20 Neoantigen Candidates

| Rank | Gene | Mutation | HLA | MT IC50 (nM) | WT IC50 (nM) | Ratio | Expr (TPM) | Immuno |
|------|------|----------|-----|--------------|--------------|-------|------------|--------|
| 1 | BRAF | V600E | A*02:01 | 18.4 | 2840 | 154x | 89.3 | 0.82 |
| 2 | NRAS | Q61R | A*24:02 | 42.1 | >5000 | >100x | 45.2 | 0.78 |
| 3 | TP53 | R248W | B*07:02 | 28.9 | 1890 | 65x | 112.4 | 0.75 |
| 4 | CDK4 | R24C | A*02:01 | 56.3 | 3210 | 57x | 67.8 | 0.71 |
| 5 | PTEN | R130Q | A*24:02 | 34.7 | 890 | 26x | 34.5 | 0.69 |
...
| 20 | ARID2 | L1205F | B*44:02 | 89.4 | 1450 | 16x | 28.9 | 0.54 |

### Lead Candidate Analysis: BRAF V600E

#### Neoantigen Peptide
```
Wild-type:  LATEKSRWS-V-GAFPVYK
Mutant:     LATEKSRWS-E-GAFPVYK
Epitope:    KIGDFGLATEKSRWSGSHQFEQLSGSILWMAPE
9-mer:      KIGDFGLAT (position 592-600)
```

#### MHC Binding Profile
| HLA Allele | MT IC50 | WT IC50 | Differential |
|------------|---------|---------|--------------|
| A*02:01 | 18.4 nM | 2840 nM | 154x |
| A*24:02 | 1240 nM | >5000 nM | >4x |
| B*07:02 | >5000 nM | >5000 nM | - |

#### Clinical Validation
- **Previous trials:** Confirmed immunogenic in melanoma patients
- **T-cell responses:** CD8+ responses documented
- **Tumor infiltration:** BRAF V600E-reactive TILs identified

### mRNA Vaccine Construct Design

#### Multi-Neoantigen Cassette
```
5'Cap-[5'UTR]-[Signal]-[NEO1]-[Linker]-[NEO2]-...-[NEO20]-[3'UTR]-PolyA
```

#### Sequence Architecture
```
Position  Element              Length
1-47      Optimized 5' UTR     47 nt
48-75     Signal peptide       27 aa
76-102    Neoantigen 1 (27aa)  27 aa
103-108   AAY linker           3 aa
109-135   Neoantigen 2 (27aa)  27 aa
...
          [Repeat for 20 neoantigens]
...
          Alpha-globin 3' UTR  142 nt
          Poly(A) tail         120 nt
```

#### Design Features
- **Total length:** ~4,200 nt
- **Neoantigens:** 20 × 27-mer (flanking sequence included)
- **Linkers:** AAY (between MHC-I) / GPGPG (between MHC-II)
- **Modifications:** N1-methylpseudouridine (m1Ψ)
- **Delivery:** LNP formulation

### Predicted Immune Response

| Neoantigen Set | CD8+ Response | CD4+ Response | HLA Coverage |
|----------------|---------------|---------------|--------------|
| Top 10 | 8/10 predicted | 6/10 predicted | 96% |
| Top 20 | 15/20 predicted | 12/20 predicted | 99% |

### Manufacturing Specifications

| Parameter | Specification |
|-----------|--------------|
| **Platform** | mRNA-LNP |
| **Dose** | 100 μg per vaccination |
| **Schedule** | Prime + 2 boosts (D0, D21, D42) |
| **Storage** | -20°C (6 months stability) |
| **Manufacturing time** | 4-6 weeks from biopsy |
```

### LLM Agent Integration

```python
@tool
def design_neoantigen_vaccine(
    vcf_file: str,
    hla_alleles: list[str],
    expression_data: str = None,
    num_neoantigens: int = 20,
    platform: str = "mRNA"
) -> str:
    """
    Designs personalized neoantigen vaccine from patient data.

    Args:
        vcf_file: Path to somatic mutation VCF
        hla_alleles: Patient HLA alleles (e.g., ["A*02:01", "B*07:02"])
        expression_data: Optional path to RNA-seq expression file
        num_neoantigens: Number of neoantigens for vaccine
        platform: mRNA, peptide, or DC

    Returns:
        Vaccine design with prioritized neoantigens
    """
    pass


@tool
def predict_neoantigen_immunogenicity(
    mutation: str,
    gene: str,
    hla_alleles: list[str],
    include_tcr_prediction: bool = False
) -> str:
    """
    Predicts immunogenicity of a specific neoantigen.

    Args:
        mutation: Mutation string (e.g., "V600E")
        gene: Gene symbol
        hla_alleles: Patient HLA alleles
        include_tcr_prediction: Include TCR binding prediction

    Returns:
        Immunogenicity assessment with binding predictions
    """
    pass
```

---

## Prerequisites

### Required Tools/Pipelines

| Resource | Purpose | Access |
|----------|---------|--------|
| **pVACtools** | Neoantigen prediction pipeline | Open source |
| **NetMHCpan** | MHC binding prediction | DTU API |
| **OptiType** | HLA typing from NGS | Open source |
| **kallisto/salmon** | Expression quantification | Open source |
| **MuPeXI** | Neoantigen prioritization | Academic |

### Dependencies

```
pvactools>=4.0
optitype>=1.3
numpy>=1.24
pandas>=1.5
biopython>=1.81
pyvcf>=0.6.8
```

---

## Methodology

### Neoantigen Prediction Pipeline

```
Tumor + Normal WES/WGS
    ↓
Variant Calling (Mutect2, Strelka2)
    ↓
HLA Typing (OptiType, HLA-HD)
    ↓
Peptide Generation (8-11 mers around mutation)
    ↓
MHC Binding Prediction (NetMHCpan)
    ├── Mutant peptide binding
    └── Wild-type peptide binding
    ↓
Expression Filtering (RNA-seq)
    ↓
Immunogenicity Scoring
    ↓
Prioritization & Ranking
    ↓
Vaccine Design
```

### Prioritization Algorithm

```python
def prioritize_neoantigens(candidates: list[dict]) -> list[dict]:
    """
    Multi-criteria neoantigen prioritization.
    """
    for neo in candidates:
        # Calculate composite score
        neo['score'] = (
            0.25 * binding_score(neo['mt_ic50']) +
            0.20 * differential_score(neo['mt_ic50'], neo['wt_ic50']) +
            0.20 * expression_score(neo['tpm']) +
            0.15 * immunogenicity_score(neo['peptide']) +
            0.10 * clonality_score(neo['vaf']) +
            0.10 * foreignness_score(neo['peptide'])
        )

    # Rank by composite score
    return sorted(candidates, key=lambda x: x['score'], reverse=True)


def binding_score(ic50: float) -> float:
    """Convert IC50 to 0-1 score (lower IC50 = higher score)."""
    if ic50 < 50:
        return 1.0
    elif ic50 < 150:
        return 0.8
    elif ic50 < 500:
        return 0.5
    else:
        return 0.0


def differential_score(mt_ic50: float, wt_ic50: float) -> float:
    """Score based on mutant/wild-type binding ratio."""
    ratio = wt_ic50 / mt_ic50
    if ratio > 100:
        return 1.0
    elif ratio > 10:
        return 0.7
    elif ratio > 5:
        return 0.5
    else:
        return 0.2
```

---

## Related Skills

- **Epitope Prediction Agent:** Core epitope prediction algorithms
- **mRNA Design Agent:** mRNA vaccine sequence optimization
- **Variant Interpretation:** Mutation annotation and filtering

---

## References

- **Hu et al. (2021):** "Personal neoantigen vaccines induce persistent memory T cell responses and epitope spreading in patients with melanoma." *Nature Medicine*
- **Wells et al. (2020):** "Key parameters of tumor epitope immunogenicity revealed through a consortium approach." *Cell*
- [pVACtools](https://pvactools.readthedocs.io/)
- [Cancer neoantigen review](https://www.nature.com/articles/s41586-019-1072-2)

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

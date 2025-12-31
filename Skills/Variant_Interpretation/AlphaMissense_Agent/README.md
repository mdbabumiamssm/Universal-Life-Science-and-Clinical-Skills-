# AlphaMissense Agent

**ID:** `biomedical.variant_interpretation.alphamissense`
**Version:** 1.0.0
**Status:** Production
**Category:** Variant Interpretation / Pathogenicity Prediction
**Released:** December 2025

---

## Overview

The **AlphaMissense Agent** provides access to AlphaMissense, DeepMind's proteome-wide missense variant effect predictor. AlphaMissense classifies **89% of all 71 million possible human missense variants** as either likely pathogenic or likely benign, dramatically accelerating rare disease diagnosis and variant interpretation.

### Key Achievement

AlphaMissense represents a quantum leap in variant interpretation:

| Metric | AlphaMissense | Previous SOTA |
|--------|---------------|---------------|
| Variants classified | 71 million | Limited |
| Confident predictions | 89% (32% P, 57% B) | ~50% |
| AUROC on ClinVar | 0.94 | ~0.85 |
| Zero-shot (no training on ClinVar) | Yes | No |

### Clinical Calibration (2025)

ClinGen has calibrated AlphaMissense for clinical use:

| Threshold | Strength | Clinical Use |
|-----------|----------|--------------|
| ≥0.564 | PP3 Supporting | Pathogenic support |
| ≥0.773 | PP3 Moderate | Stronger evidence |
| ≤0.336 | BP4 Supporting | Benign support |
| ≤0.197 | BP4 Moderate | Stronger evidence |

---

## Key Capabilities

### 1. Variant Classification

| Category | Score Range | % of Proteome |
|----------|-------------|---------------|
| **Likely Pathogenic** | ≥0.564 | 32% |
| **Uncertain** | 0.336-0.564 | 11% |
| **Likely Benign** | ≤0.336 | 57% |

### 2. Disease Applications

| Disease Category | Performance | Clinical Utility |
|------------------|-------------|------------------|
| Rare Mendelian | Excellent | Diagnostic |
| Cancer (somatic) | Good | Prognostic |
| Pharmacogenomics | Moderate | Actionable |
| Complex disease | Variable | Research |

### 3. Technical Specifications

```yaml
Model Architecture:
  Base: AlphaFold2 structure module
  Training: Primate conservation + population frequency
  Output: Per-residue pathogenicity score

Coverage:
  Proteins: 19,233 human canonical
  Variants: 71 million missense
  Resolution: Single amino acid

Validation:
  ClinVar AUROC: 0.94
  De novo enrichment: 4.3x
  Calibration: ClinGen-approved
```

---

## Usage

### Example Prompt

```text
I have a pediatric patient with developmental delay and seizures.
WES identified the following variants of uncertain significance (VUS):

1. SCN1A c.1178G>A (p.Arg393His)
2. KCNQ2 c.854C>T (p.Pro285Leu)
3. MECP2 c.473C>T (p.Thr158Met)
4. STXBP1 c.1217G>A (p.Arg406His)

Please use AlphaMissense to prioritize these variants for pathogenicity.
```

### Expected Output

```markdown
## AlphaMissense Variant Analysis: Developmental Epilepsy Panel

### Patient Context
- Phenotype: Developmental delay, seizures
- Testing: Whole exome sequencing
- VUS to evaluate: 4 variants

### Variant Pathogenicity Predictions

| Gene | Variant | AM Score | Classification | ClinVar | Priority |
|------|---------|----------|----------------|---------|----------|
| **SCN1A** | p.Arg393His | 0.891 | Likely Pathogenic | VUS | **HIGH** |
| **STXBP1** | p.Arg406His | 0.847 | Likely Pathogenic | VUS | **HIGH** |
| **MECP2** | p.Thr158Met | 0.956 | Likely Pathogenic | P/LP | Confirmed |
| **KCNQ2** | p.Pro285Leu | 0.423 | Uncertain | VUS | MODERATE |

### Detailed Analysis

#### 1. SCN1A p.Arg393His
```
AlphaMissense Score: 0.891
─────────────────────────────────
Classification: LIKELY PATHOGENIC

Structural Context:
├── Domain: Voltage-sensing (S4 segment)
├── Conservation: Arg conserved in all vertebrates
├── Structure: AlphaFold confidence pLDDT = 95
└── Contacts: Critical gating charge

Clinical Evidence:
├── ClinVar: VUS (1 submission)
├── gnomAD: Absent (AC=0)
├── ACMG: PM2, PP3(moderate)
└── Literature: None found

Functional Prediction:
├── Mechanism: Loss of gating charge
├── Effect: Likely loss-of-function
└── Phenotype fit: HIGH (Dravet spectrum)

Recommendation: RECLASSIFY to Likely Pathogenic
Supporting criteria: PM2, PP3_moderate, PP4
```

#### 2. STXBP1 p.Arg406His
```
AlphaMissense Score: 0.847
─────────────────────────────────
Classification: LIKELY PATHOGENIC

Structural Context:
├── Domain: Domain 3a (syntaxin binding)
├── Conservation: Invariant Arg in metazoans
├── Structure: Core interaction residue
└── Contacts: Syntaxin-1 interface

Clinical Evidence:
├── ClinVar: VUS (2 submissions)
├── gnomAD: Absent (AC=0)
├── Same position: p.Arg406Cys is P
└── Literature: Nearby mutations cause EIEE

Functional Prediction:
├── Mechanism: Disrupted SNARE binding
├── Effect: Likely loss-of-function
└── Phenotype fit: HIGH (STXBP1 encephalopathy)

Recommendation: RECLASSIFY to Likely Pathogenic
Supporting criteria: PM2, PM5, PP3_moderate
```

#### 3. MECP2 p.Thr158Met
```
AlphaMissense Score: 0.956
─────────────────────────────────
Classification: LIKELY PATHOGENIC

Clinical Status: Already P/LP in ClinVar

Notes:
├── Well-established Rett syndrome mutation
├── One of most common MECP2 mutations
├── High AM score confirms pathogenicity
└── No additional analysis needed

Note: Patient is male - consider X-linked inheritance
implications and phenotype spectrum
```

#### 4. KCNQ2 p.Pro285Leu
```
AlphaMissense Score: 0.423
─────────────────────────────────
Classification: UNCERTAIN

Structural Context:
├── Domain: S5 segment
├── Conservation: Pro conserved in mammals
├── Structure: Helix-breaking residue
└── Effect: Unclear structural impact

Clinical Evidence:
├── ClinVar: VUS (1 submission)
├── gnomAD: 1 heterozygote (AF=4e-6)
├── ACMG: PM2 only
└── Literature: Not reported

Functional Prediction:
├── Mechanism: Possible helix disruption
├── Effect: Uncertain
└── Phenotype fit: MODERATE

Recommendation: Maintain VUS
- Consider functional testing
- Segregation analysis if family available
```

### Prioritization Summary

| Priority | Variant | Action |
|----------|---------|--------|
| 1 | SCN1A p.Arg393His | Clinical correlation, consider Dravet |
| 2 | STXBP1 p.Arg406His | Clinical correlation |
| 3 | MECP2 p.Thr158Met | Already confirmed P |
| 4 | KCNQ2 p.Pro285Leu | Needs functional data |

### Diagnostic Hypothesis

Based on AlphaMissense scores and phenotype:

**Most likely diagnosis**: SCN1A-related epilepsy (Dravet syndrome spectrum)

Rationale:
1. SCN1A p.Arg393His has highest relevance to phenotype
2. Score 0.891 strongly supports pathogenicity
3. Loss-of-function in S4 consistent with Dravet
4. De novo testing recommended

**Alternative**: STXBP1 encephalopathy if SCN1A ruled out

### Recommended Actions

1. **Parental testing** for de novo status of SCN1A variant
2. **Clinical re-evaluation** for Dravet features (fever-sensitivity, SCN1A-specific EEG)
3. **Functional validation** if diagnosis uncertain
4. **ACMG update** for SCN1A and STXBP1 variants
```

---

## LLM Agent Integration

### Python Tool Implementation

```python
from typing import Optional, Dict, Any, List, Literal
import pandas as pd

def alphamissense_tool(
    variants: List[Dict[str, str]],
    genome_build: str = "GRCh38",
    include_structure: bool = True,
    include_population: bool = True,
    clinical_calibration: bool = True,
    batch_size: int = 100
) -> Dict[str, Any]:
    """
    AlphaMissense variant pathogenicity prediction.

    Args:
        variants: List of variants with gene, transcript, HGVSp
        genome_build: Reference genome (GRCh37 or GRCh38)
        include_structure: Include AlphaFold structural context
        include_population: Include gnomAD frequencies
        clinical_calibration: Apply ClinGen thresholds
        batch_size: Batch size for processing

    Returns:
        Dictionary with pathogenicity predictions and annotations
    """
    from alphamissense import AlphaMissenseDB, VariantAnnotator

    # Initialize database
    db = AlphaMissenseDB(genome_build=genome_build)
    annotator = VariantAnnotator(db)

    results = []
    for variant in variants:
        # Parse variant
        gene = variant.get("gene")
        hgvsp = variant.get("hgvsp")
        hgvsc = variant.get("hgvsc")

        # Query AlphaMissense
        am_result = db.query(gene=gene, hgvsp=hgvsp)

        if am_result is None:
            results.append({
                "variant": variant,
                "am_score": None,
                "classification": "Not found",
                "error": "Variant not in AlphaMissense database"
            })
            continue

        # Get score and classification
        score = am_result["pathogenicity"]

        if clinical_calibration:
            if score >= 0.773:
                classification = "Likely Pathogenic"
                strength = "PP3_moderate"
            elif score >= 0.564:
                classification = "Likely Pathogenic"
                strength = "PP3_supporting"
            elif score <= 0.197:
                classification = "Likely Benign"
                strength = "BP4_moderate"
            elif score <= 0.336:
                classification = "Likely Benign"
                strength = "BP4_supporting"
            else:
                classification = "Uncertain"
                strength = "None"
        else:
            classification = am_result.get("am_class", "Unknown")
            strength = None

        result = {
            "variant": variant,
            "am_score": float(score),
            "classification": classification,
            "acmg_strength": strength
        }

        # Add structural context
        if include_structure:
            struct_info = annotator.get_structure_context(gene, hgvsp)
            result["structure"] = {
                "domain": struct_info.get("domain"),
                "secondary_structure": struct_info.get("ss"),
                "plddt": struct_info.get("plddt"),
                "contacts": struct_info.get("contacts")
            }

        # Add population data
        if include_population:
            pop_info = annotator.get_population_freq(gene, hgvsp)
            result["population"] = {
                "gnomad_af": pop_info.get("af"),
                "gnomad_ac": pop_info.get("ac"),
                "clinvar_status": pop_info.get("clinvar")
            }

        results.append(result)

    # Summary statistics
    scores = [r["am_score"] for r in results if r["am_score"] is not None]
    classifications = [r["classification"] for r in results]

    summary = {
        "total_variants": len(variants),
        "scored": len(scores),
        "likely_pathogenic": classifications.count("Likely Pathogenic"),
        "likely_benign": classifications.count("Likely Benign"),
        "uncertain": classifications.count("Uncertain"),
        "mean_score": sum(scores) / len(scores) if scores else None
    }

    return {
        "results": results,
        "summary": summary,
        "genome_build": genome_build,
        "clinical_calibration": clinical_calibration
    }


def batch_alphamissense_analysis(
    vcf_path: str,
    output_path: str,
    filter_score: Optional[float] = None
) -> Dict[str, Any]:
    """
    Batch analysis of VCF file with AlphaMissense.

    Args:
        vcf_path: Path to VCF file
        output_path: Output file path
        filter_score: Minimum score to include

    Returns:
        Analysis summary
    """
    from cyvcf2 import VCF
    from alphamissense import AlphaMissenseDB

    db = AlphaMissenseDB()
    vcf = VCF(vcf_path)

    results = []
    for variant in vcf:
        if variant.is_snp:
            # Get AlphaMissense score
            chrom = variant.CHROM
            pos = variant.POS
            ref = variant.REF
            alt = variant.ALT[0]

            score = db.query_genomic(chrom, pos, ref, alt)

            if score is not None:
                if filter_score is None or score >= filter_score:
                    results.append({
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "am_score": score
                    })

    # Save results
    df = pd.DataFrame(results)
    df.to_csv(output_path, index=False)

    return {
        "total_variants": len(list(VCF(vcf_path))),
        "missense_scored": len(results),
        "pathogenic": len(df[df["am_score"] >= 0.564]),
        "output_path": output_path
    }


# Claude tool schema
ALPHAMISSENSE_SCHEMA = {
    "name": "predict_variant_pathogenicity",
    "description": "Predict pathogenicity of missense variants using AlphaMissense. Provides proteome-wide predictions with ClinGen-calibrated clinical thresholds.",
    "input_schema": {
        "type": "object",
        "properties": {
            "variants": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "gene": {"type": "string"},
                        "hgvsp": {"type": "string"},
                        "hgvsc": {"type": "string"}
                    }
                },
                "description": "List of variants to analyze"
            },
            "clinical_calibration": {
                "type": "boolean",
                "description": "Apply ClinGen PP3/BP4 thresholds"
            },
            "include_structure": {
                "type": "boolean",
                "description": "Include AlphaFold structural context"
            }
        },
        "required": ["variants"]
    }
}
```

---

## Prerequisites

### Data Access

```python
# Download AlphaMissense database
from alphamissense import download_database

# Full database (~6 GB)
download_database(version="latest", output_dir="./data")

# Or use API access
from alphamissense import AlphaMissenseAPI
api = AlphaMissenseAPI(api_key="your_key")
```

### Installation

```bash
pip install alphamissense
pip install cyvcf2  # For VCF processing
pip install pandas numpy
```

### Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| alphamissense | >=1.0 | Core predictions |
| pandas | >=2.0 | Data handling |
| cyvcf2 | >=0.30 | VCF parsing |
| requests | >=2.28 | API access |

---

## Methodology

### Model Architecture

```
AlphaMissense Architecture
━━━━━━━━━━━━━━━━━━━━━━━━━

Input: Protein sequence + variant position
              │
              ▼
┌─────────────────────────────────┐
│   AlphaFold2 Structure Module   │
│   (Pretrained weights)          │
└────────────────┬────────────────┘
                 │
                 ▼
┌─────────────────────────────────┐
│   Conservation Features         │
│   - MSA from primate genomes    │
│   - Within-species frequency    │
└────────────────┬────────────────┘
                 │
                 ▼
┌─────────────────────────────────┐
│   Pathogenicity Head            │
│   - Per-residue prediction      │
│   - Calibrated score [0,1]      │
└────────────────┬────────────────┘
                 │
                 ▼
     Pathogenicity Score (0-1)
```

### Training Strategy

| Component | Source | Purpose |
|-----------|--------|---------|
| Structure | AlphaFold2 | 3D context |
| Conservation | Primate genomes | Evolutionary constraint |
| Weak labels | gnomAD frequency | Population selection |

### Validation

| Dataset | Metric | Score |
|---------|--------|-------|
| ClinVar P/B | AUROC | 0.940 |
| De novo variants | Enrichment | 4.3x |
| LoF-constrained genes | Accuracy | 0.92 |
| Functional assays | Correlation | 0.85 |

---

## Limitations

| Limitation | Impact | Mitigation |
|------------|--------|------------|
| Novel proteins | No prediction | Limited coverage |
| Synonymous variants | Not supported | Use other tools |
| Structural variants | Not supported | Use SV tools |
| Context-dependent | May miss | Clinical correlation |

---

## Related Skills

- `biomedical.variant_interpretation.diagAI` - Diagnostic AI
- `biomedical.variant_interpretation.dyna` - Disease-specific VEP
- `biomedical.genomics.variant_interpretation` - General interpretation
- `biomedical.genomics.crispr_design` - Guide RNA design

---

## References

1. **AlphaMissense (2023)**: "Accurate proteome-wide missense variant effect prediction with AlphaMissense." *Science*. [DOI: 10.1126/science.adg7492](https://www.science.org/doi/10.1126/science.adg7492)

2. **ClinGen Calibration (2025)**: "Calibration of additional computational tools expands ClinGen recommendation options for variant classification." *Genetics in Medicine*.

3. **AlphaMissense Database**: [https://alphamissense.hegelab.org/](https://alphamissense.hegelab.org/)

4. **DeepMind Blog**: [AlphaMissense Announcement](https://deepmind.google/discover/blog/alphamissense/)

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

---

*Last updated: December 2025*

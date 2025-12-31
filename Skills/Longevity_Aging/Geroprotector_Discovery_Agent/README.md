# Geroprotector Discovery Agent

**ID:** `biomedical.longevity.geroprotector_discovery`
**Version:** 1.0.0
**Status:** Production
**Category:** Longevity / Drug Discovery

---

## Overview

The **Geroprotector Discovery Agent** identifies molecules that promote healthy aging and extend lifespan using AI-powered drug discovery platforms. A 2025 Scripps Research study showed that over 70% of AI-identified anti-aging drugs significantly extended lifespan in C. elegans validation studies.

This agent uses machine learning to analyze bioactivity data from known geroprotectors to identify new molecules targeting aging hallmarks, enabling rapid discovery of anti-aging therapeutics.

---

## Key Capabilities

### 1. Aging Hallmark Targets

| Hallmark | Intervention Strategy | Drug Examples |
|----------|----------------------|---------------|
| **Genomic instability** | DNA repair enhancement | Poly-ADP ribose inhibitors |
| **Telomere attrition** | Telomerase activation | TA-65 |
| **Epigenetic alterations** | Sirtuin activation | Resveratrol, SRT1720 |
| **Proteostasis loss** | Autophagy induction | Rapamycin |
| **Nutrient sensing** | mTOR inhibition | Rapamycin |
| **Mitochondrial dysfunction** | Mito-enhancers | NAD+ precursors |
| **Cellular senescence** | Senolytics | Dasatinib+Quercetin |
| **Stem cell exhaustion** | Regeneration | GH secretagogues |
| **Altered intercellular** | Anti-inflammatory | Metformin |

### 2. Discovery Approaches

| Approach | Method | Advantage |
|----------|--------|-----------|
| **Target-based** | Hallmark pathway modulation | Mechanistic clarity |
| **Phenotypic** | Lifespan screening | Unbiased discovery |
| **Repurposing** | Existing drug analysis | Fast translation |
| **De novo** | Generative chemistry | Novel scaffolds |

### 3. AI/ML Models

- **AgeXtend:** Multimodal geroprotector prediction (IIT-Delhi)
- **ClockBase Agent:** Intervention identification from omics data
- **Precious3GPT:** Multi-species aging transformer

### 4. Validation Readouts

| Model | Readout | Throughput |
|-------|---------|------------|
| **C. elegans** | Lifespan, healthspan | High |
| **D. melanogaster** | Lifespan, locomotion | Medium |
| **Mouse** | Lifespan, biomarkers | Low |
| **Human cells** | Senescence markers | High |
| **Organoids** | Tissue aging | Medium |

---

## Usage

### Example Prompt

```text
Identify potential geroprotectors from a library of FDA-approved drugs.

Criteria:
- Target at least 2 aging hallmarks
- Acceptable safety profile for chronic use
- Predicted to reduce biological age
- Mechanism-based prioritization

Provide top 10 candidates with rationale.
```

### Expected Output

```
## Geroprotector Discovery Report

### Screening Parameters
- **Library:** FDA-approved drugs (2,400 compounds)
- **Model:** AgeXtend + Target pathway analysis
- **Criteria:** Multi-hallmark targeting, chronic safety

### Top 10 Geroprotector Candidates

| Rank | Drug | Original Indication | Geroprotector Score | Hallmarks |
|------|------|---------------------|---------------------|-----------|
| 1 | **Metformin** | Diabetes | 0.94 | 4 |
| 2 | **Rapamycin** | Immunosuppression | 0.92 | 3 |
| 3 | **Acarbose** | Diabetes | 0.87 | 3 |
| 4 | **Canagliflozin** | Diabetes | 0.85 | 3 |
| 5 | **Lithium** | Bipolar disorder | 0.82 | 2 |
| 6 | **Aspirin** | Pain/Cardiovascular | 0.80 | 2 |
| 7 | **Atorvastatin** | Hyperlipidemia | 0.78 | 2 |
| 8 | **Pioglitazone** | Diabetes | 0.76 | 2 |
| 9 | **Lisinopril** | Hypertension | 0.74 | 2 |
| 10 | **Verapamil** | Cardiovascular | 0.72 | 2 |

### Detailed Analysis: Top 3 Candidates

#### 1. Metformin (Score: 0.94)

**Aging Hallmarks Targeted:**
| Hallmark | Mechanism | Evidence |
|----------|-----------|----------|
| Nutrient sensing | AMPK activation | Strong |
| Mitochondrial function | Complex I inhibition | Strong |
| Inflammation | NF-κB inhibition | Moderate |
| Cellular senescence | SASP reduction | Moderate |

**Clinical Evidence:**
- TAME trial ongoing (Targeting Aging with Metformin)
- Epidemiological: Diabetics on metformin outlive non-diabetics
- Dose: 500-2000 mg/day

**Safety for Chronic Use:**
| Aspect | Rating | Notes |
|--------|--------|-------|
| GI tolerance | Moderate | B12 supplementation needed |
| Lactic acidosis | Rare | Contraindicated in renal failure |
| Long-term data | Excellent | 60+ years of use |

#### 2. Rapamycin (Score: 0.92)

**Aging Hallmarks Targeted:**
| Hallmark | Mechanism | Evidence |
|----------|-----------|----------|
| Nutrient sensing | mTOR inhibition | Very strong |
| Proteostasis | Autophagy induction | Strong |
| Cellular senescence | Senescence delay | Moderate |

**Clinical Evidence:**
- Mouse lifespan extension: 9-14%
- Human: Improved immune response in elderly (TORC trial)
- Dose: Intermittent dosing (weekly) shows benefits with fewer side effects

**Safety Considerations:**
| Aspect | Rating | Notes |
|--------|--------|-------|
| Immunosuppression | Concern | Intermittent dosing reduces risk |
| Metabolic effects | Moderate | Glucose elevation |
| Wound healing | Impaired | Pause before procedures |

#### 3. Acarbose (Score: 0.87)

**Aging Hallmarks Targeted:**
| Hallmark | Mechanism | Evidence |
|----------|-----------|----------|
| Nutrient sensing | Glucose reduction | Strong |
| Gut microbiome | SCFA increase | Moderate |
| Inflammation | AGE reduction | Moderate |

**Clinical Evidence:**
- ITP mouse study: 22% lifespan extension (males)
- Human: Reduces post-prandial glucose
- Dose: 25-100 mg with meals

### Novel Mechanism Candidates

#### Ouabain (Cardiac glycoside)
- **ClockBase Agent identification**
- **Mechanism:** Na+/K+-ATPase modulation
- **Evidence:** Reduced frailty in aged mice
- **Caution:** Narrow therapeutic window

### Combination Strategies

| Combination | Rationale | Synergy |
|-------------|-----------|---------|
| Metformin + Rapamycin | AMPK + mTOR | Additive |
| Dasatinib + Quercetin | Senolytic | Synergistic |
| Metformin + Acarbose | Dual glucose | Additive |

### Predicted Biological Age Impact

| Intervention | Predicted ΔBioAge | Confidence |
|--------------|-------------------|------------|
| Metformin (1500 mg/d) | -2.1 years | High |
| Rapamycin (5 mg/wk) | -1.8 years | Medium |
| Acarbose (50 mg tid) | -1.2 years | Medium |
| Combination (Met+Rapa) | -3.5 years | Low |

### Validation Recommendations

| Candidate | Suggested Model | Duration |
|-----------|-----------------|----------|
| Top 5 drugs | C. elegans lifespan | 4 weeks |
| Top 3 drugs | Mouse healthspan | 12 months |
| Combinations | Cell senescence assay | 2 weeks |
```

### LLM Agent Integration

```python
@tool
def discover_geroprotectors(
    compound_library: str = "fda_approved",
    hallmarks_target: list[str] = None,
    num_candidates: int = 10,
    include_combinations: bool = True
) -> str:
    """
    Identifies geroprotector candidates from compound library.

    Args:
        compound_library: fda_approved, drugbank, or custom
        hallmarks_target: List of aging hallmarks to target
        num_candidates: Number of top candidates
        include_combinations: Suggest drug combinations

    Returns:
        Ranked geroprotector candidates with mechanisms
    """
    pass
```

---

## References

- **Petrascheck et al. (2025):** "AI pinpoints new anti-aging drug candidates." *Aging Cell*
- **Moskalev et al. (2022):** "Geroprotectors.org: a new database of geroprotectors." *Aging*
- [TAME Trial](https://www.afar.org/tame-trial)

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

# Biological Age Clock Agent

**ID:** `biomedical.longevity.biological_age_clock`
**Version:** 1.0.0
**Status:** Production
**Category:** Longevity / Aging Research

---

## Overview

The **Biological Age Clock Agent** estimates biological age using multi-modal data including DNA methylation, blood biomarkers, gene expression, and physiological measurements. Biological age provides a more accurate measure of health status and mortality risk than chronological age.

A 2025 Nature paper showed that biological age of the brain and immune system predicted long-term healthspan better than chronological age. This agent implements multiple aging clocks and provides comprehensive longevity assessment.

---

## Key Capabilities

### 1. Aging Clock Types

| Clock | Input Data | Output | Accuracy |
|-------|------------|--------|----------|
| **Horvath** | DNA methylation (353 CpGs) | DNAm Age | MAE ~3.6 years |
| **Hannum** | DNA methylation (71 CpGs) | Blood DNAm Age | MAE ~4.0 years |
| **PhenoAge** | Blood biomarkers (9 markers) | Phenotypic Age | Mortality predictor |
| **GrimAge** | DNAm + clinical | Lifespan prediction | Best mortality |
| **DunedinPACE** | Longitudinal DNAm | Pace of aging | Rate predictor |

### 2. Multi-Modal Assessment

| Data Type | Markers | Application |
|-----------|---------|-------------|
| **Epigenetic** | CpG methylation | Gold standard |
| **Blood biomarkers** | CRP, albumin, creatinine | Accessible |
| **Transcriptomic** | Gene expression | Functional state |
| **Proteomic** | Plasma proteins | SomaScan |
| **Metabolomic** | Metabolite patterns | Metabolic health |
| **Imaging** | Brain MRI, retinal | Organ-specific |

### 3. Age Acceleration Metrics

- **Age acceleration:** Biological age - Chronological age
- **Intrinsic age acceleration:** Adjusted for cell composition
- **Extrinsic age acceleration:** Captures immune aging
- **Pace of aging:** Rate of biological aging

### 4. Organ-Specific Clocks

| Organ | Clock Type | Clinical Relevance |
|-------|------------|-------------------|
| **Brain** | NeuroAge | Dementia risk |
| **Immune** | ImmuneAge | Immunosenescence |
| **Cardiovascular** | VascularAge | CVD risk |
| **Metabolic** | MetabolicAge | Diabetes risk |

---

## Usage

### Example Prompt

```text
Calculate biological age for a patient using available data:

Demographics:
- Chronological age: 55 years
- Sex: Male

Blood biomarkers:
- Albumin: 4.2 g/dL
- Creatinine: 1.1 mg/dL
- Glucose: 105 mg/dL
- CRP: 2.5 mg/L
- Lymphocyte %: 28%
- Mean cell volume: 92 fL
- Red blood cell distribution width: 13.5%
- Alkaline phosphatase: 65 U/L
- White blood cell count: 6.2 x10³/μL

Calculate PhenoAge and provide health interpretation.
```

### Expected Output

```
## Biological Age Assessment Report

### Patient Profile
- **Chronological Age:** 55 years
- **Sex:** Male

### PhenoAge Calculation

#### Input Biomarkers
| Biomarker | Value | Reference Range | Status |
|-----------|-------|-----------------|--------|
| Albumin | 4.2 g/dL | 3.5-5.0 | Normal |
| Creatinine | 1.1 mg/dL | 0.7-1.3 | Normal |
| Glucose | 105 mg/dL | 70-100 | Elevated |
| CRP | 2.5 mg/L | <3.0 | Borderline |
| Lymphocyte % | 28% | 20-40% | Normal |
| MCV | 92 fL | 80-100 | Normal |
| RDW | 13.5% | 11.5-14.5 | Normal |
| ALP | 65 U/L | 44-147 | Normal |
| WBC | 6.2 x10³/μL | 4.5-11 | Normal |

#### PhenoAge Calculation
```
PhenoAge = 141.50225
         + (-0.00553 × ln(albumin))
         + (-0.02664 × ln(creatinine))
         + (0.00188 × glucose)
         + (0.02691 × ln(CRP))
         + (-0.01247 × lymphocyte%)
         + (-0.03420 × MCV)
         + (0.00931 × RDW)
         + (0.00226 × ALP)
         + (0.01167 × WBC)
         + 0.0804 × age
```

**PhenoAge = 58.2 years**

### Age Acceleration Analysis

| Metric | Value | Interpretation |
|--------|-------|----------------|
| **Chronological Age** | 55.0 years | - |
| **PhenoAge** | 58.2 years | - |
| **Age Acceleration** | +3.2 years | Older than expected |
| **Percentile** | 72nd | Higher biological age than 72% of peers |

### Health Risk Assessment

#### Mortality Risk
| Timeframe | Relative Risk | Interpretation |
|-----------|---------------|----------------|
| 5-year | 1.18 | 18% elevated risk |
| 10-year | 1.25 | 25% elevated risk |

#### Risk Factors Identified
| Factor | Contribution | Modifiable |
|--------|--------------|------------|
| Elevated glucose | +1.2 years | Yes |
| Borderline CRP | +0.8 years | Yes |
| Normal albumin | -0.3 years | - |

### Recommendations

#### Priority Interventions
| Intervention | Target | Expected Impact |
|--------------|--------|-----------------|
| **Blood sugar control** | Fasting glucose <100 | -1.0 year biological age |
| **Reduce inflammation** | CRP <1.0 | -0.5 year biological age |
| **Exercise** | 150 min/week moderate | -1.5 years |
| **Mediterranean diet** | Increase plant-based | -1.0 year |

#### Monitoring Schedule
| Test | Frequency | Target |
|------|-----------|--------|
| Fasting glucose | 3 months | <100 mg/dL |
| CRP | 6 months | <1.0 mg/L |
| Full panel | 12 months | Re-calculate PhenoAge |

### Comparison with Other Clocks

| Clock | Estimated Bio Age | Age Acceleration |
|-------|-------------------|------------------|
| PhenoAge | 58.2 | +3.2 years |
| GrimAge (estimated) | 57.8 | +2.8 years |
| Biological Age (composite) | 58.0 | +3.0 years |

### Long-term Health Trajectory

```
Biological Age vs Chronological Age

Bio Age |     Current trajectory
  65    |                    ●
        |                ●
  60    |            ●  ← Intervention
        |        ●
  55    |    ●------- Optimal trajectory
        |●
  50    +---+---+---+---+---+---
        50  55  60  65  70  75
              Chronological Age
```

**With interventions:** Biological age could reduce to 55-56 years within 12 months
```

### LLM Agent Integration

```python
@tool
def calculate_biological_age(
    chronological_age: float,
    sex: str,
    biomarkers: dict,
    clock_type: str = "phenoage",
    include_recommendations: bool = True
) -> str:
    """
    Calculates biological age from biomarker data.

    Args:
        chronological_age: Age in years
        sex: male or female
        biomarkers: Dict of biomarker values
        clock_type: phenoage, horvath, hannum, grimage
        include_recommendations: Include lifestyle recommendations

    Returns:
        Biological age assessment with interpretation
    """
    pass
```

---

## Prerequisites

### Required Data

| Clock | Data Required | Platform |
|-------|---------------|----------|
| **Horvath/Hannum** | DNA methylation array | Illumina 450K/EPIC |
| **PhenoAge** | Blood panel | Standard clinical |
| **GrimAge** | DNA methylation + clinical | Illumina + EMR |
| **DunedinPACE** | Longitudinal methylation | Research |

### Dependencies

```
numpy>=1.24
pandas>=1.5
scikit-learn>=1.3
scipy>=1.10
```

---

## References

- **Levine et al. (2018):** "An epigenetic biomarker of aging for lifespan and healthspan." *Aging*
- **Belsky et al. (2022):** "DunedinPACE, a DNA methylation biomarker of the pace of aging." *eLife*
- **Galkin et al. (2025):** "Biological age prediction using multi-omics data." *Nature Medicine*

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

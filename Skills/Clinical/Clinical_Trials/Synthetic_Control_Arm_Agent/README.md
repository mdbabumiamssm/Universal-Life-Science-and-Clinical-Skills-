# Synthetic Control Arm Agent

**ID:** `biomedical.clinical.synthetic_control`
**Version:** 1.0.0
**Status:** Production
**Category:** Clinical AI / Trial Design

---

## Overview

The **Synthetic Control Arm Agent** generates virtual control groups from real-world data for clinical trials, reducing the need for concurrent placebo participants. Synthetic control arms use AI to simulate outcomes of a placebo group based on historical patient data, accelerating trials and addressing ethical concerns.

FDA and EMA have provided guidance on using real-world evidence for synthetic controls, particularly for rare diseases and oncology.

---

## Key Capabilities

### 1. Data Sources

| Source | Type | Use Case |
|--------|------|----------|
| **EHR** | Structured clinical | Disease progression |
| **Claims** | Administrative | Treatment patterns |
| **Registries** | Disease-specific | Natural history |
| **Previous trials** | Clinical | Historical controls |

### 2. Matching Methods

| Method | Description | Advantage |
|--------|-------------|-----------|
| **Propensity score** | Probability of treatment | Confounding reduction |
| **Exact matching** | Identical covariates | Interpretability |
| **CEM** | Coarsened exact matching | Balance guarantee |
| **IPTW** | Inverse probability weighting | Full sample use |

### 3. Validation Approaches

- **Replication studies:** Compare SCA to actual RCT arms
- **Negative controls:** Known null effects
- **Sensitivity analysis:** Unmeasured confounding

---

## Usage

### Example Prompt

```text
Generate a synthetic control arm for a Phase 3 oncology trial.

Trial design:
- Single-arm trial of novel targeted therapy
- Stage IV NSCLC, EGFR mutation positive
- Prior treatment with osimertinib (progressed)
- Primary endpoint: Progression-free survival

Create SCA from historical data matching trial eligibility.
```

### Expected Output

```
## Synthetic Control Arm Generation Report

### Trial Summary
- **Indication:** Stage IV NSCLC, EGFR+, post-osimertinib
- **Design:** Single-arm + synthetic control
- **Primary endpoint:** Progression-free survival (PFS)
- **Target N (active):** 120 patients

### Synthetic Control Construction

#### Data Sources
| Source | Patients | Time Period | Quality |
|--------|----------|-------------|---------|
| Flatiron Health | 2,456 | 2019-2024 | High |
| IQVIA Oncology | 1,823 | 2020-2024 | Medium |
| Academic center | 312 | 2018-2024 | High |
| **Total pool** | 4,591 | - | - |

#### Eligibility Matching
| Criterion | Pool Size |
|-----------|-----------|
| Initial pool | 4,591 |
| EGFR mutation + | 892 |
| Prior osimertinib | 567 |
| Progressed on osimertinib | 389 |
| ECOG 0-1 | 312 |
| No brain mets (or treated) | 278 |
| Adequate organ function | 245 |
| **Final matched** | 180 |

### Covariate Balance (Post-Matching)

| Variable | Active Arm | SCA | SMD |
|----------|------------|-----|-----|
| Age (mean) | 62.4 | 61.8 | 0.05 |
| Female % | 58% | 56% | 0.04 |
| ECOG 0 % | 42% | 40% | 0.04 |
| Prior lines | 2.3 | 2.4 | 0.06 |
| Brain mets % | 28% | 26% | 0.05 |
| Time on osimertinib | 14.2 mo | 13.8 mo | 0.04 |

**All SMD < 0.1:** Excellent balance

### Synthetic Control Arm Results

#### PFS Outcome
```
Endpoint          SCA (n=180)   95% CI
Median PFS        4.2 months    (3.8-4.8)
6-month PFS       32%           (25-39%)
12-month PFS      12%           (8-17%)
```

#### Sensitivity Analyses

| Analysis | Median PFS | Range |
|----------|------------|-------|
| Primary (PS matching) | 4.2 mo | - |
| Exact matching | 4.0 mo | - |
| IPTW | 4.4 mo | - |
| Trimmed (PS 0.1-0.9) | 4.1 mo | - |

**Consistency:** High (all methods within 0.4 months)

### Regulatory Considerations

| Element | Status | Notes |
|---------|--------|-------|
| FDA guidance alignment | ✓ | RWE for external control |
| Data quality documented | ✓ | Oncology-specific EHR |
| Sensitivity analyses | ✓ | Multiple methods |
| Unmeasured confounding | ⚠️ | E-value analysis recommended |

### Power Calculation

| Scenario | HR | Power | Conclusion |
|----------|-----|-------|------------|
| Active PFS = 6 mo | 0.70 | 82% | Powered |
| Active PFS = 5 mo | 0.84 | 56% | Underpowered |

**Recommendation:** Proceed with SCA design if expected HR ≤ 0.75
```

---

## References

- **FDA Guidance (2021):** "Real-World Data: Assessing Electronic Health Records and Medical Claims Data"
- **Sanofi case study:** Digital twins saved 6 months and millions in asthma trial

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

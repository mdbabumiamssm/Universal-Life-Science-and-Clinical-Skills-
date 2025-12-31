# 3D Spheroid Profiler Agent

**ID:** `biomedical.research_tools.spheroid_profiler`
**Version:** 1.0.0
**Status:** Beta
**Category:** Research Tools / High-Content Imaging

---

## Overview

The **3D Spheroid Profiler Agent** extends Cell Painting analysis to 3D spheroid cultures, providing more physiologically relevant phenotypic profiling. A February 2025 preprint introduced a scalable method combining Cell Painting with tissue-clearing for 3D spheroids.

This agent analyzes morphological features from 3D spheroid images, capturing complex phenotypes not observable in 2D monolayers.

---

## Key Capabilities

### 1. 3D-Specific Features

| Feature Category | Description | 2D vs 3D |
|-----------------|-------------|----------|
| **Volume** | Total spheroid volume | 3D only |
| **Core phenotype** | Hypoxic/necrotic core | 3D only |
| **Penetration** | Compound distribution | 3D only |
| **Invasion** | Cell migration from spheroid | 3D only |
| **Heterogeneity** | Spatial variation | Enhanced in 3D |

### 2. Imaging Methods

| Method | Resolution | Throughput | Application |
|--------|------------|------------|-------------|
| **Confocal** | Highest | Low | Mechanism studies |
| **Light-sheet** | High | Medium | Large spheroids |
| **Widefield + clearing** | Medium | High | Screening |
| **Multiphoton** | High | Low | Deep imaging |

### 3. Spheroid Types

| Type | Complexity | Applications |
|------|------------|--------------|
| **Tumor spheroids** | Medium | Cancer drug screening |
| **Organoids** | High | Disease modeling |
| **Co-culture** | High | TME modeling |
| **Patient-derived** | Very high | Precision medicine |

---

## Usage

### Example Prompt

```text
Analyze Cell Painting images from 3D tumor spheroids.

Input:
- MCF7 breast cancer spheroids (384-well plate)
- 50 compounds + DMSO controls
- 5-channel Cell Painting
- Cleared with CUBIC protocol

Extract 3D-specific features and identify compounds affecting
spheroid architecture.
```

### Expected Output

```
## 3D Spheroid Profiling Report

### Spheroid Quality Assessment
| Metric | Mean ± SD | Acceptable Range |
|--------|-----------|------------------|
| Diameter | 425 ± 45 μm | 350-500 μm ✓ |
| Roundness | 0.91 ± 0.04 | >0.85 ✓ |
| Viability (core) | 82 ± 8% | >70% ✓ |

### 3D Feature Extraction

| Feature Class | Count | Examples |
|---------------|-------|----------|
| **Morphology** | 45 | Volume, surface area, eccentricity |
| **Texture (3D)** | 180 | 3D Haralick features |
| **Spatial** | 85 | Core-rim gradient, radial distribution |
| **Channel correlation** | 60 | 3D colocalization |
| **Total** | 370 | - |

### Active Compounds (Spheroid-Specific Effects)

| Compound | Volume Change | Core Effect | Invasion | Phenotype |
|----------|---------------|-------------|----------|-----------|
| CPD-023 | -45% | Necrosis | ↓ | Cytotoxic |
| CPD-089 | +15% | Preserved | ↑↑ | Pro-invasive |
| CPD-045 | -20% | Normal | ↓ | Growth inhibition |
| CPD-112 | Normal | ↓Hypoxia | Normal | Oxygenation |

### 3D-Specific Insights

**CPD-089 (Pro-invasive phenotype):**
- Not detected in 2D assay
- 3D reveals EMT-like morphology
- Cells dissociating from spheroid edge
- **Recommendation:** Exclude from oncology pipeline

**CPD-112 (Core penetration):**
- Reduces hypoxic core
- May enhance drug delivery
- **Recommendation:** Consider as combination agent
```

---

## References

- **Booij et al. (2025):** "High-content morphological profiling by Cell Painting in 3D spheroids." *bioRxiv*
- **Sirenko et al. (2015):** "High-content assays for characterizing the viability and morphology of 3D cancer spheroid cultures."

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

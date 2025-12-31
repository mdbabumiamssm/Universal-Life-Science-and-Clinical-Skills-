# Cell Painting Analysis Agent

**ID:** `biomedical.research_tools.cell_painting`
**Version:** 1.0.0
**Status:** Production
**Category:** Research Tools / High-Content Imaging

---

## Overview

The **Cell Painting Analysis Agent** extracts morphological profiles from Cell Painting images for drug discovery applications. Cell Painting is a high-content imaging assay that uses multiplexed fluorescent dyes to visualize multiple cellular compartments simultaneously, enabling target-agnostic phenotypic profiling.

This agent converts high-content cellular images into quantitative phenotypic profiles using deep learning, enabling mechanism-of-action prediction, compound prioritization, and toxicity assessment.

---

## Key Capabilities

### 1. Cell Painting Channels

| Channel | Dye | Target | Features |
|---------|-----|--------|----------|
| **DNA** | Hoechst 33342 | Nucleus | Size, shape, texture |
| **ER** | Concanavalin A | Endoplasmic reticulum | Pattern, distribution |
| **RNA/Nucleoli** | SYTO 14 | RNA | Intensity, granularity |
| **AGP** | Phalloidin | F-actin, Golgi | Cytoskeleton, organelles |
| **Mito** | MitoTracker | Mitochondria | Mass, texture, distribution |

### 2. Feature Extraction

| Category | Features | Count |
|----------|----------|-------|
| **Intensity** | Mean, median, std | ~200 |
| **Texture** | Haralick, Gabor | ~400 |
| **Shape** | Area, eccentricity | ~100 |
| **Radial** | Distribution profiles | ~200 |
| **Correlation** | Channel colocalization | ~100 |
| **Total** | - | ~1,000 |

### 3. Analysis Methods

| Method | Description | Application |
|--------|-------------|-------------|
| **CellProfiler** | Classical feature extraction | Standard pipeline |
| **PhenoProfiler** | Deep learning embeddings | High accuracy |
| **Self-supervised** | Contrastive learning | Label-free |
| **Diffusion models** | Generative analysis | Phenotype prediction |

### 4. Applications

- **MOA prediction:** Compound mechanism clustering
- **Toxicity screening:** Adverse effect identification
- **Target ID:** Phenotypic target deconvolution
- **Hit prioritization:** Lead optimization guidance

---

## Usage

### Example Prompt

```text
Analyze Cell Painting images from a compound screen.

Input:
- 384-well plate images (5 channels)
- 100 test compounds at 10 μM
- DMSO negative controls
- Reference compounds with known MOA

Output:
- Morphological profiles
- Compound clustering by phenotype
- MOA predictions
- Toxicity flags
```

### Expected Output

```
## Cell Painting Analysis Report

### Plate Overview
- **Format:** 384-well plate
- **Compounds:** 100 test + 20 reference
- **Replicates:** 4 per compound
- **Cells/well:** ~1,500 average

### Quality Control

| Metric | Value | Threshold | Status |
|--------|-------|-----------|--------|
| Cell count CV | 12.3% | <25% | ✓ |
| Focus score | 0.89 | >0.7 | ✓ |
| DMSO consistency | 0.92 | >0.85 | ✓ |
| Positive control Z' | 0.68 | >0.5 | ✓ |

### Feature Extraction

| Pipeline | Features | Processing Time |
|----------|----------|-----------------|
| CellProfiler | 1,078 | 4.2 hours |
| PhenoProfiler (DL) | 512 | 1.8 hours |

### Phenotypic Activity

#### Active Compounds (>2 SD from DMSO)
| Rank | Compound | Distance | Strongest Channel | Cytotoxic |
|------|----------|----------|-------------------|-----------|
| 1 | CPD-042 | 8.4 | Mitochondria | No |
| 2 | CPD-018 | 7.2 | DNA | Yes |
| 3 | CPD-067 | 6.8 | ER | No |
| 4 | CPD-091 | 5.9 | F-actin | No |
| 5 | CPD-003 | 5.7 | Nucleoli | No |
...

**Active compounds:** 34/100 (34%)
**Cytotoxic compounds:** 8/100 (8%)

### MOA Clustering

#### Cluster Analysis (UMAP + HDBSCAN)
```
Cluster 1: Mitochondrial disruptors (n=8)
├── CPD-042
├── CPD-055
├── Reference: Rotenone
└── Reference: Oligomycin

Cluster 2: DNA damage (n=5)
├── CPD-018
├── CPD-033
└── Reference: Etoposide

Cluster 3: ER stress (n=6)
├── CPD-067
├── CPD-079
└── Reference: Tunicamycin

Cluster 4: Cytoskeletal (n=7)
├── CPD-091
├── CPD-012
└── Reference: Cytochalasin D
```

### MOA Predictions (Top Compounds)

#### CPD-042 (Highest activity)

| Predicted MOA | Confidence | Reference Match |
|---------------|------------|-----------------|
| **Mitochondrial ETC inhibitor** | 0.87 | Rotenone (r=0.82) |
| Oxidative phosphorylation | 0.72 | Oligomycin (r=0.71) |
| Respiratory complex I | 0.68 | - |

**Phenotypic Signature:**
```
Channel        Effect      Magnitude
Mito           ↓ Mass      -2.3 SD
Mito           ↑ Texture   +1.8 SD
DNA            ↑ Area      +0.9 SD
ER             Normal      -
AGP            ↓ Pattern   -0.7 SD
```

#### CPD-067 (ER Stress)

| Predicted MOA | Confidence | Reference Match |
|---------------|------------|-----------------|
| **ER stress inducer** | 0.81 | Tunicamycin (r=0.78) |
| Protein glycosylation inhibitor | 0.65 | - |
| UPR activator | 0.58 | - |

### Toxicity Assessment

| Compound | Cell Count | Nuclei Shape | Mito Health | Overall Risk |
|----------|------------|--------------|-------------|--------------|
| CPD-018 | -45% | Abnormal | Depolarized | **High** |
| CPD-033 | -38% | Fragmented | Abnormal | **High** |
| CPD-042 | -5% | Normal | Stressed | Low |
| CPD-067 | -8% | Normal | Normal | Low |

### Dose-Response Profiles (Select Compounds)

| Compound | EC50 | Max Effect | Hill Slope | Selectivity |
|----------|------|------------|------------|-------------|
| CPD-042 | 2.3 μM | 95% | 1.4 | Mitochondria |
| CPD-067 | 5.1 μM | 82% | 1.1 | ER |
| CPD-091 | 8.7 μM | 71% | 0.9 | Actin |

### Recommendations

| Priority | Compound | Action | Rationale |
|----------|----------|--------|-----------|
| 1 | CPD-042 | Advance | Strong mito signature, low toxicity |
| 2 | CPD-067 | Advance | Novel ER mechanism |
| 3 | CPD-091 | Deprioritize | Weak effect |
| - | CPD-018 | Reject | High cytotoxicity |
```

### LLM Agent Integration

```python
@tool
def analyze_cell_painting(
    image_path: str,
    compound_metadata: str,
    reference_compounds: list[str] = None,
    extract_method: str = "phenoprofiler"
) -> str:
    """
    Analyzes Cell Painting images for phenotypic profiling.

    Args:
        image_path: Path to image directory
        compound_metadata: CSV with well-compound mapping
        reference_compounds: Known MOA reference compounds
        extract_method: cellprofiler or phenoprofiler

    Returns:
        Morphological profiles with MOA predictions
    """
    pass
```

---

## References

- **Bray et al. (2016):** "Cell Painting, a high-content image-based assay for morphological profiling." *Nature Protocols*
- **Chandrasekaran et al. (2025):** "PhenoProfiler: AI-powered phenotypic analysis." *Nature Communications*
- [JUMP-Cell Painting Consortium](https://jump-cellpainting.broadinstitute.org/)

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

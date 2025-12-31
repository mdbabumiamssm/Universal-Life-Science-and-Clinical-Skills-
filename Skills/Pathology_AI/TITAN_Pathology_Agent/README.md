# TITAN Pathology Agent

**ID:** `biomedical.pathology_ai.titan`
**Version:** 1.0.0
**Status:** Production
**Category:** Pathology AI / Foundation Models
**Released:** December 2025

---

## Overview

The **TITAN Pathology Agent** provides access to TITAN (Transformers for Immunohistochemistry and Tissue Analysis in Nature), a multimodal whole-slide foundation model pretrained on **335,645 whole-slide images (WSIs)**. TITAN provides representations for both slide-level and patient-level clinical tasks.

### Key Capabilities

TITAN can perform clinical tasks and generate reports even in **data-scarce scenarios** such as rare cancer diagnosis and survival prediction, without requiring further fine-tuning. This makes it particularly valuable for:

- Rare cancer detection
- Multi-site generalization
- Zero-shot and few-shot diagnostics

### Comparison with Other Pathology FMs

| Model | Training WSIs | Modalities | Zero-Shot | Publication |
|-------|--------------|------------|-----------|-------------|
| **TITAN** | 335,645 | Vision + Text | Yes | Nature Medicine 2025 |
| **Virchow** | 1,500,000 | Vision | Limited | Nature Medicine 2024 |
| **CONCH** | 1,170,000 | Vision + Text | Yes | Nature Medicine 2024 |
| **Phikon** | 460,000 | Vision | Limited | bioRxiv 2024 |

---

## Key Capabilities

### 1. Clinical Tasks

| Task | Description | Zero-Shot | Fine-Tuned |
|------|-------------|-----------|------------|
| **Cancer Detection** | Identify malignant regions | High | Excellent |
| **Tumor Classification** | Histological subtyping | Good | Excellent |
| **Grade Prediction** | Tumor grade assessment | Moderate | Excellent |
| **Survival Prediction** | Prognostic modeling | Good | Excellent |
| **Biomarker Prediction** | IHC surrogate | Moderate | Good |
| **Report Generation** | Automated pathology reports | Good | Excellent |

### 2. Supported Tissue Types

| Organ System | Cancer Types | Performance |
|--------------|--------------|-------------|
| **Breast** | IDC, ILC, DCIS | Excellent |
| **Lung** | NSCLC, SCLC | Excellent |
| **Colon** | Adenocarcinoma | Excellent |
| **Prostate** | Adenocarcinoma | Excellent |
| **Skin** | Melanoma, BCC, SCC | Excellent |
| **Kidney** | RCC subtypes | Good |
| **Liver** | HCC, cholangiocarcinoma | Good |
| **Pancreas** | PDAC | Good |
| **Brain** | Gliomas | Good |
| **Rare Cancers** | Various | Good (few-shot) |

### 3. Technical Specifications

```yaml
Architecture:
  Vision Encoder: ViT-L/14 (pretrained)
  Slide Aggregation: Attention-based MIL
  Text Encoder: BiomedCLIP
  Fusion: Cross-modal attention

Training:
  WSIs: 335,645
  Patients: 100,000+
  Institutions: 15+
  Augmentation: Stain normalization, rotation, crop

Inference:
  Tile Size: 224x224 @ 20x
  Memory: Variable (adaptive tiling)
  Speed: ~30 seconds per WSI (GPU)
```

---

## Usage

### Example Prompt

```text
Analyze this lung biopsy whole-slide image for:
1. Presence of malignancy
2. If malignant, classify the tumor type
3. Assess tumor grade
4. Identify any molecular markers suggested by morphology
5. Generate a preliminary pathology report

WSI Path: /pathology/slides/lung_biopsy_001.svs
Clinical Context: 65-year-old male smoker with lung mass on CT
```

### Expected Output

```markdown
## TITAN Pathology Analysis: Lung Biopsy

### Slide Quality Assessment
| Parameter | Value | Status |
|-----------|-------|--------|
| Tissue coverage | 78% | Adequate |
| Stain quality | Good | Acceptable |
| Focus | Sharp | Acceptable |
| Artifacts | Minimal | Acceptable |

### Malignancy Detection

**Result**: POSITIVE for malignancy

| Region | Confidence | Area (mm²) |
|--------|------------|------------|
| Primary tumor | 0.98 | 12.4 |
| Peripheral foci | 0.94 | 2.1 |
| Lymphovascular | 0.72 | N/A |

### Tumor Classification

| Diagnosis | Probability | Features |
|-----------|-------------|----------|
| **Lung Adenocarcinoma** | 0.92 | Glandular formation, mucin production |
| Squamous Cell Carcinoma | 0.05 | Limited keratinization |
| Large Cell Carcinoma | 0.02 | - |
| Small Cell Carcinoma | 0.01 | No neuroendocrine features |

**Primary Diagnosis**: Invasive Lung Adenocarcinoma

### Histological Pattern Analysis

| Pattern | Percentage | Prognostic Impact |
|---------|------------|-------------------|
| Acinar | 45% | Intermediate |
| Papillary | 25% | Intermediate |
| Lepidic | 15% | Favorable |
| Solid | 10% | Unfavorable |
| Micropapillary | 5% | Unfavorable |

**Predominant Pattern**: Acinar

### Grade Assessment

| Grading System | Grade | Criteria |
|----------------|-------|----------|
| IASLC (2020) | Grade 2 | Acinar/papillary predominant |
| Architectural | Moderate | Gland formation present |
| Nuclear | Intermediate | Moderate pleomorphism |

### Morphology-Based Biomarker Prediction

| Marker | Predicted Status | Confidence | Notes |
|--------|------------------|------------|-------|
| TTF-1 | Positive | 0.89 | Adenocarcinoma morphology |
| Napsin A | Positive | 0.85 | Glandular differentiation |
| p40 | Negative | 0.82 | No squamous features |
| EGFR mutation | Possible | 0.65 | Consider molecular testing |
| ALK rearrangement | Low probability | 0.78 | Solid component present |
| PD-L1 (≥50%) | Possible | 0.58 | Immune infiltrate noted |

**Recommendation**: Confirm with IHC and molecular testing

### Tumor Microenvironment

```
Immune Infiltrate Analysis:
─────────────────────────────
TILs (tumor-infiltrating lymphocytes):
├── Intratumoral: Moderate (15-20%)
├── Stromal: Dense (40-50%)
└── Tertiary lymphoid structures: Present

Stroma:
├── Desmoplastic reaction: Moderate
├── Fibrosis: 30% of tumor area
└── Necrosis: Minimal (<5%)

Vascular Invasion:
├── Lymphovascular: Suspicious
└── Perineural: Not identified
```

### Survival Prediction

| Timepoint | Predicted Survival | 95% CI |
|-----------|-------------------|--------|
| 1-year | 88% | 82-94% |
| 3-year | 65% | 55-75% |
| 5-year | 48% | 38-58% |

**Risk Category**: Intermediate

### Preliminary Pathology Report

---
**SURGICAL PATHOLOGY REPORT (AI-ASSISTED DRAFT)**

**Specimen**: Lung, right lower lobe, CT-guided biopsy

**Clinical History**: 65-year-old male with 40 pack-year smoking
history presenting with 3.2 cm spiculated mass in right lower lobe.

**DIAGNOSIS**:
Right lower lobe lung, core biopsy:
- **INVASIVE ADENOCARCINOMA, ACINAR PREDOMINANT**
- Grade 2 (IASLC grading system)
- Lymphovascular invasion: Suspicious
- Perineural invasion: Not identified

**TUMOR CHARACTERISTICS**:
- Histologic patterns: Acinar (45%), papillary (25%),
  lepidic (15%), solid (10%), micropapillary (5%)
- Nuclear grade: Intermediate
- Necrosis: Minimal

**ANCILLARY STUDIES RECOMMENDED**:
- Immunohistochemistry panel: TTF-1, Napsin A, p40, CK7
- PD-L1 (22C3 or SP263)
- Molecular testing: EGFR, ALK, ROS1, BRAF, KRAS, MET

**NOTE**: This is an AI-generated preliminary report requiring
pathologist review and sign-out.

---

### Quality Metrics

| Metric | Value |
|--------|-------|
| Analysis time | 28 seconds |
| Tiles analyzed | 12,456 |
| Confidence (overall) | 0.91 |
| Uncertainty regions | 3 flagged |
```

---

## LLM Agent Integration

### Python Tool Implementation

```python
from typing import Optional, Dict, Any, List, Literal
from pathlib import Path
import openslide
import torch
import numpy as np

def titan_pathology_tool(
    wsi_path: str,
    task: Literal["detect", "classify", "grade", "predict", "report"],
    clinical_context: Optional[str] = None,
    organ_system: Optional[str] = None,
    tile_size: int = 224,
    magnification: str = "20x",
    return_heatmap: bool = True,
    device: str = "cuda"
) -> Dict[str, Any]:
    """
    TITAN pathology foundation model for WSI analysis.

    Args:
        wsi_path: Path to whole-slide image (.svs, .ndpi, .tiff)
        task: Analysis task type
        clinical_context: Patient clinical information
        organ_system: Organ/tissue type for specialized analysis
        tile_size: Size of analysis tiles
        magnification: Analysis magnification level
        return_heatmap: Generate attention heatmap
        device: Compute device

    Returns:
        Dictionary with pathology analysis results
    """
    from titan import TITANModel, SlideProcessor, ReportGenerator

    # Load model
    model = TITANModel.from_pretrained("mass-general/TITAN")
    model.to(device)
    model.eval()

    # Process slide
    processor = SlideProcessor(tile_size=tile_size, magnification=magnification)
    slide = openslide.OpenSlide(wsi_path)
    tiles, coordinates = processor.extract_tiles(slide)

    # Encode tiles
    with torch.no_grad():
        tile_embeddings = model.encode_tiles(tiles.to(device))
        slide_embedding = model.aggregate(tile_embeddings, coordinates)

    results = {
        "wsi_path": wsi_path,
        "task": task,
        "n_tiles": len(tiles)
    }

    if task == "detect":
        # Cancer detection
        detection = model.detect_cancer(slide_embedding)
        results["malignancy"] = {
            "positive": detection["probability"] > 0.5,
            "probability": float(detection["probability"]),
            "confidence": float(detection["confidence"])
        }

        if return_heatmap:
            tile_scores = model.get_tile_attention(tile_embeddings)
            results["attention_heatmap"] = create_heatmap(
                coordinates, tile_scores, slide.dimensions
            )

    elif task == "classify":
        # Tumor classification
        if organ_system:
            classification = model.classify_tumor(
                slide_embedding,
                organ_system=organ_system
            )
        else:
            classification = model.classify_tumor(slide_embedding)

        results["classification"] = {
            "primary_diagnosis": classification["top_class"],
            "probabilities": classification["probabilities"],
            "features": classification["morphological_features"]
        }

    elif task == "grade":
        # Tumor grading
        grading = model.predict_grade(slide_embedding)
        results["grade"] = {
            "grade": grading["grade"],
            "grading_system": grading["system"],
            "components": grading["components"]
        }

    elif task == "predict":
        # Survival/outcome prediction
        prediction = model.predict_outcome(slide_embedding)
        results["prognosis"] = {
            "risk_score": float(prediction["risk_score"]),
            "risk_category": prediction["category"],
            "survival_curves": prediction["survival_curves"]
        }

    elif task == "report":
        # Generate full report
        report_gen = ReportGenerator(model)
        report = report_gen.generate(
            slide_embedding,
            tile_embeddings,
            coordinates,
            clinical_context=clinical_context
        )
        results["report"] = report

    return results


def create_heatmap(coordinates: np.ndarray,
                   scores: np.ndarray,
                   dimensions: tuple) -> np.ndarray:
    """Create attention heatmap from tile scores."""
    heatmap = np.zeros((dimensions[1] // 16, dimensions[0] // 16))
    for (x, y), score in zip(coordinates, scores):
        heatmap[y//16, x//16] = score
    return heatmap


# Batch processing for large cohorts
def batch_titan_analysis(
    wsi_paths: List[str],
    task: str = "classify",
    batch_size: int = 4,
    output_dir: str = "./results"
) -> Dict[str, Any]:
    """
    Process multiple WSIs in batch.

    Args:
        wsi_paths: List of WSI file paths
        task: Analysis task
        batch_size: Processing batch size
        output_dir: Output directory

    Returns:
        Aggregated results
    """
    from titan import TITANModel
    import pandas as pd

    model = TITANModel.from_pretrained("mass-general/TITAN")
    model.cuda()

    results = []
    for path in wsi_paths:
        try:
            result = titan_pathology_tool(
                wsi_path=path,
                task=task,
                return_heatmap=False
            )
            results.append({"path": path, "status": "success", **result})
        except Exception as e:
            results.append({"path": path, "status": "failed", "error": str(e)})

    # Save results
    df = pd.DataFrame(results)
    df.to_csv(f"{output_dir}/titan_results.csv", index=False)

    return {
        "total": len(wsi_paths),
        "successful": len([r for r in results if r["status"] == "success"]),
        "failed": len([r for r in results if r["status"] == "failed"]),
        "results": results
    }


# Claude tool schema
TITAN_TOOL_SCHEMA = {
    "name": "analyze_pathology_slide",
    "description": "Analyze whole-slide pathology images using TITAN foundation model. Supports cancer detection, tumor classification, grading, survival prediction, and report generation.",
    "input_schema": {
        "type": "object",
        "properties": {
            "wsi_path": {
                "type": "string",
                "description": "Path to whole-slide image file"
            },
            "task": {
                "type": "string",
                "enum": ["detect", "classify", "grade", "predict", "report"],
                "description": "Analysis task type"
            },
            "clinical_context": {
                "type": "string",
                "description": "Patient clinical information for context"
            },
            "organ_system": {
                "type": "string",
                "description": "Organ/tissue type"
            }
        },
        "required": ["wsi_path", "task"]
    }
}
```

---

## Prerequisites

### Software Requirements

| Package | Version | Purpose |
|---------|---------|---------|
| OpenSlide | >=3.4 | WSI reading |
| PyTorch | >=2.0 | Deep learning |
| TITAN | >=1.0 | Foundation model |
| Pillow | >=9.0 | Image processing |
| NumPy | >=1.24 | Array operations |

### Installation

```bash
# System dependencies (Ubuntu)
sudo apt-get install openslide-tools libopensilde-dev

# Python packages
pip install openslide-python torch torchvision
pip install titan-pathology  # Or from GitHub

# Download model weights
python -c "from titan import TITANModel; TITANModel.from_pretrained('mass-general/TITAN')"
```

### Hardware Requirements

| WSI Size | GPU VRAM | Processing Time |
|----------|----------|-----------------|
| Small (<1 GB) | 8 GB | ~15 seconds |
| Medium (1-5 GB) | 16 GB | ~30 seconds |
| Large (>5 GB) | 24+ GB | ~60 seconds |

---

## Methodology

### Architecture

```
TITAN Architecture
━━━━━━━━━━━━━━━━━━

WSI Input (gigapixel)
         │
         ▼
┌─────────────────────────────┐
│   Tile Extraction           │
│   224x224 @ 20x/40x         │
│   ~10,000-50,000 tiles      │
└────────────┬────────────────┘
             │
             ▼
┌─────────────────────────────┐
│   Vision Encoder (ViT-L)    │
│   Per-tile embeddings       │
│   768-dimensional           │
└────────────┬────────────────┘
             │
             ▼
┌─────────────────────────────┐
│   Attention-based MIL       │
│   Aggregate to slide-level  │
│   Learn key tile weights    │
└────────────┬────────────────┘
             │
             ▼
┌─────────────────────────────┐
│   Text Encoder (BiomedCLIP) │
│   Clinical context encoding │
│   Report generation         │
└────────────┬────────────────┘
             │
             ▼
      Task-Specific Heads
      ├── Detection
      ├── Classification
      ├── Grading
      ├── Survival
      └── Report Generation
```

---

## Regulatory Status

| Approval | Status | Notes |
|----------|--------|-------|
| FDA 510(k) | Pending | Class II device pathway |
| CE Mark | In progress | EU MDR compliance |
| Research Use | Available | Not for clinical diagnosis |

---

## Related Skills

- `biomedical.pathology_ai.conch` - Vision-language pathology
- `biomedical.pathology_ai.virchow` - Large-scale pathology FM
- `biomedical.clinical.h_optimus` - H-Optimus pipeline
- `biomedical.radiology.ark_chestxray` - Radiology foundation model

---

## References

1. **TITAN (2025)**: "A multimodal whole-slide foundation model for pathology." *Nature Medicine*. [DOI: 10.1038/s41591-025-03982-3](https://www.nature.com/articles/s41591-025-03982-3)

2. **Mass General Brigham (2025)**: "Researchers Develop AI Foundation Models to Advance Pathology." [Press Release](https://www.massgeneralbrigham.org/)

3. **Pathology FM Review (2025)**: "Pathology Foundation Models." *JMA Journal*.

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

---

*Last updated: December 2025*

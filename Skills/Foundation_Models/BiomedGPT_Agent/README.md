# BiomedGPT Agent

**ID:** `biomedical.foundation_models.biomedgpt`
**Version:** 1.0.0
**Status:** Production
**Category:** Foundation Models / Multimodal Biomedical AI
**Released:** December 2025

---

## Overview

The **BiomedGPT Agent** provides access to BiomedGPT, the first open-source and lightweight vision-language foundation model designed as a generalist capable of performing diverse biomedical tasks. BiomedGPT achieved state-of-the-art results in **16 out of 25 experiments** while maintaining a computing-friendly model scale.

In July 2025, larger-scale checkpoints were released (up to 5x larger, 930M parameters), including stronger large and xlarge pre-trained models for enhanced performance.

### Key Performance Metrics

| Task | Error Rate | Comparison |
|------|------------|------------|
| Medical Question Answering | 3.8% | Robust prediction |
| Radiology Report Writing | 8.3% | Satisfactory |
| Clinical Summarization | ~Equivalent | Human expert level |

---

## Key Capabilities

### 1. Supported Modalities

| Modality | Input Types | Applications |
|----------|-------------|--------------|
| **Medical Images** | X-ray, CT, MRI, pathology | Diagnosis, findings detection |
| **Clinical Text** | Reports, notes, literature | QA, summarization, extraction |
| **Genomics** | Sequence data, variant info | Interpretation support |
| **Multimodal** | Image + text combinations | VQA, report generation |

### 2. Task Capabilities

| Task Category | Specific Tasks | Performance |
|---------------|----------------|-------------|
| **Visual QA** | VQA-RAD, PathVQA, SLAKE | SOTA on 8/9 datasets |
| **Report Generation** | Radiology, pathology | 8.3% error rate |
| **Image Classification** | Chest X-ray, dermoscopy | High accuracy |
| **Text QA** | PubMedQA, BioASQ | Competitive |
| **Summarization** | Clinical notes, papers | Expert-level |
| **Information Extraction** | NER, relation extraction | Strong performance |

### 3. Model Variants

| Variant | Parameters | Use Case | VRAM Required |
|---------|------------|----------|---------------|
| BiomedGPT-S | 182M | Quick inference | 4 GB |
| BiomedGPT-B | 370M | Balanced | 8 GB |
| BiomedGPT-L | 560M | High accuracy | 12 GB |
| BiomedGPT-XL | 930M | Maximum performance | 24 GB |

---

## Usage

### Example Prompt

```text
I have a chest X-ray image of a patient presenting with cough and fever.
Please analyze the image and:
1. Identify any abnormal findings
2. Provide differential diagnoses
3. Generate a preliminary radiology report
4. Suggest follow-up imaging if needed

Image path: /data/chest_xray_001.png
```

### Expected Output

```markdown
## BiomedGPT Chest X-ray Analysis

### Image Quality Assessment
| Parameter | Status |
|-----------|--------|
| Positioning | PA view, adequate |
| Exposure | Appropriate |
| Inspiration | Full |
| Artifacts | None |

### Findings

#### Lungs
- **Right lower lobe**: Patchy airspace opacity with air bronchograms
- **Left lung**: Clear, no focal consolidation
- **Pleural spaces**: No effusion bilaterally

#### Mediastinum
- Heart size: Normal (CTR < 0.5)
- Mediastinal contours: Unremarkable
- Trachea: Midline

#### Bones/Soft Tissues
- No acute osseous abnormality
- Soft tissues unremarkable

### Impression

```
1. Right lower lobe consolidation, most consistent with pneumonia
   - Community-acquired pneumonia most likely given clinical context
   - Consider bacterial vs. atypical pathogens

2. No pleural effusion or pneumothorax

3. No cardiomegaly
```

### Differential Diagnosis

| Diagnosis | Probability | Supporting Features |
|-----------|-------------|---------------------|
| Bacterial pneumonia | 75% | Lobar consolidation, air bronchograms |
| Atypical pneumonia | 15% | Could explain fever, cough pattern |
| Aspiration | 8% | RLL location typical |
| Malignancy | 2% | Less likely given acute presentation |

### Recommendations

1. **Clinical correlation** with laboratory findings (WBC, CRP, procalcitonin)
2. **Sputum culture** if productive cough
3. **Follow-up X-ray** in 4-6 weeks if no clinical improvement
4. **CT chest** if:
   - No improvement with antibiotics
   - Suspicion for underlying malignancy
   - Complicated pneumonia (abscess, empyema)

### Generated Report (Draft)

> **CHEST X-RAY, PA VIEW**
>
> **CLINICAL INDICATION**: Cough, fever
>
> **COMPARISON**: None available
>
> **FINDINGS**:
> There is patchy airspace opacity in the right lower lobe with
> air bronchograms, consistent with consolidation. The left lung
> is clear. There is no pleural effusion or pneumothorax. Heart
> size is normal. The mediastinum is unremarkable.
>
> **IMPRESSION**:
> Right lower lobe pneumonia. Recommend clinical correlation and
> follow-up imaging as clinically indicated.

**Confidence Score**: 0.92
```

---

## LLM Agent Integration

### Python Tool Implementation

```python
from typing import Optional, Dict, Any, Union
from pathlib import Path
import torch
from PIL import Image
import base64

def biomedgpt_tool(
    query: str,
    image_path: Optional[str] = None,
    image_base64: Optional[str] = None,
    task: str = "vqa",
    model_size: str = "large",
    max_length: int = 512,
    temperature: float = 0.7,
    return_confidence: bool = True
) -> Dict[str, Any]:
    """
    BiomedGPT multimodal biomedical foundation model.

    Args:
        query: Text query or instruction
        image_path: Path to medical image file
        image_base64: Base64-encoded image data
        task: One of 'vqa', 'report', 'classify', 'summarize', 'extract'
        model_size: 'small', 'base', 'large', 'xlarge'
        max_length: Maximum output token length
        temperature: Sampling temperature
        return_confidence: Whether to return confidence scores

    Returns:
        Dictionary with response and metadata
    """
    from biomedgpt import BiomedGPTModel, BiomedGPTProcessor

    # Model mapping
    model_map = {
        "small": "biomedgpt-s",
        "base": "biomedgpt-b",
        "large": "biomedgpt-l",
        "xlarge": "biomedgpt-xl"
    }

    # Load model
    model = BiomedGPTModel.from_pretrained(model_map[model_size])
    processor = BiomedGPTProcessor.from_pretrained(model_map[model_size])

    # Process image if provided
    image = None
    if image_path:
        image = Image.open(image_path).convert("RGB")
    elif image_base64:
        import io
        image_data = base64.b64decode(image_base64)
        image = Image.open(io.BytesIO(image_data)).convert("RGB")

    # Prepare inputs based on task
    task_prompts = {
        "vqa": f"Question: {query}\nAnswer:",
        "report": f"Generate a detailed radiology report for this image. Clinical context: {query}\n\nReport:",
        "classify": f"Classify this medical image. Classes to consider: {query}\n\nClassification:",
        "summarize": f"Summarize the following clinical text:\n{query}\n\nSummary:",
        "extract": f"Extract medical entities from:\n{query}\n\nEntities:"
    }

    prompt = task_prompts.get(task, query)

    # Process inputs
    if image:
        inputs = processor(
            text=prompt,
            images=image,
            return_tensors="pt"
        )
    else:
        inputs = processor(
            text=prompt,
            return_tensors="pt"
        )

    # Generate response
    with torch.no_grad():
        outputs = model.generate(
            **inputs,
            max_length=max_length,
            temperature=temperature,
            do_sample=temperature > 0,
            return_dict_in_generate=True,
            output_scores=return_confidence
        )

    # Decode response
    response = processor.decode(outputs.sequences[0], skip_special_tokens=True)

    result = {
        "response": response,
        "task": task,
        "model": model_map[model_size]
    }

    # Add confidence if requested
    if return_confidence and hasattr(outputs, 'scores'):
        import torch.nn.functional as F
        scores = [F.softmax(s, dim=-1).max().item() for s in outputs.scores]
        result["confidence"] = sum(scores) / len(scores)

    return result


# Claude tool schema
BIOMEDGPT_TOOL_SCHEMA = {
    "name": "biomedgpt_analysis",
    "description": "Analyze medical images and text using BiomedGPT multimodal foundation model. Supports VQA, report generation, classification, and text analysis.",
    "input_schema": {
        "type": "object",
        "properties": {
            "query": {
                "type": "string",
                "description": "Question, instruction, or text to analyze"
            },
            "image_path": {
                "type": "string",
                "description": "Path to medical image file (X-ray, CT, MRI, pathology)"
            },
            "task": {
                "type": "string",
                "enum": ["vqa", "report", "classify", "summarize", "extract"],
                "description": "Analysis task type"
            },
            "model_size": {
                "type": "string",
                "enum": ["small", "base", "large", "xlarge"],
                "description": "Model variant to use"
            }
        },
        "required": ["query"]
    }
}
```

### Batch Processing

```python
def batch_biomedgpt_analysis(
    items: list[dict],
    model_size: str = "large",
    batch_size: int = 8
) -> list[dict]:
    """
    Process multiple medical images/texts in batches.

    Args:
        items: List of dicts with 'query' and optional 'image_path'
        model_size: Model variant
        batch_size: Processing batch size

    Returns:
        List of results
    """
    from biomedgpt import BiomedGPTModel, BiomedGPTProcessor

    model = BiomedGPTModel.from_pretrained(f"biomedgpt-{model_size[0]}")
    processor = BiomedGPTProcessor.from_pretrained(f"biomedgpt-{model_size[0]}")

    results = []
    for i in range(0, len(items), batch_size):
        batch = items[i:i+batch_size]

        texts = [item['query'] for item in batch]
        images = [
            Image.open(item['image_path']).convert("RGB")
            if 'image_path' in item else None
            for item in batch
        ]

        # Process batch
        inputs = processor(
            text=texts,
            images=[img for img in images if img],
            return_tensors="pt",
            padding=True
        )

        with torch.no_grad():
            outputs = model.generate(**inputs, max_length=512)

        responses = processor.batch_decode(outputs, skip_special_tokens=True)
        results.extend(responses)

    return results
```

---

## Prerequisites

### Required APIs/Software

| Component | Version | Installation |
|-----------|---------|--------------|
| Python | >=3.9 | System |
| PyTorch | >=2.0 | `pip install torch` |
| Transformers | >=4.35 | `pip install transformers` |
| BiomedGPT | >=1.0 | `pip install biomedgpt` |
| Pillow | >=9.0 | `pip install pillow` |

### Model Downloads

```bash
# Download model weights
from huggingface_hub import snapshot_download

# Small model (182M)
snapshot_download("microsoft/biomedgpt-s")

# Large model (560M) - Recommended
snapshot_download("microsoft/biomedgpt-l")

# XLarge model (930M) - Best quality
snapshot_download("microsoft/biomedgpt-xl")
```

### Hardware Requirements

| Model Size | GPU VRAM | CPU RAM | Disk Space |
|------------|----------|---------|------------|
| Small | 4 GB | 8 GB | 400 MB |
| Base | 8 GB | 16 GB | 800 MB |
| Large | 12 GB | 24 GB | 1.2 GB |
| XLarge | 24 GB | 32 GB | 2 GB |

---

## Methodology

### Architecture Overview

```
┌──────────────────────────────────────────────────────────────┐
│                    BiomedGPT Architecture                     │
├──────────────────────────────────────────────────────────────┤
│                                                               │
│   ┌─────────────┐         ┌─────────────┐                    │
│   │   Image     │         │    Text     │                    │
│   │  Encoder    │         │  Encoder    │                    │
│   │ (ViT-based) │         │ (BERT-like) │                    │
│   └──────┬──────┘         └──────┬──────┘                    │
│          │                       │                            │
│          └───────────┬───────────┘                           │
│                      ↓                                        │
│          ┌─────────────────────┐                             │
│          │   Cross-Modal       │                             │
│          │   Attention Fusion  │                             │
│          └──────────┬──────────┘                             │
│                     ↓                                         │
│          ┌─────────────────────┐                             │
│          │   Unified Decoder   │                             │
│          │   (Autoregressive)  │                             │
│          └──────────┬──────────┘                             │
│                     ↓                                         │
│          ┌─────────────────────┐                             │
│          │   Output Generator  │                             │
│          │   (Text Response)   │                             │
│          └─────────────────────┘                             │
│                                                               │
└──────────────────────────────────────────────────────────────┘
```

### Training Data

| Dataset Type | Size | Sources |
|--------------|------|---------|
| Medical Images | 2M+ | MIMIC-CXR, ChestX-ray14, PadChest |
| Pathology | 500K+ | TCGA, PathVQA |
| Radiology Reports | 400K+ | MIMIC-III, OpenI |
| Biomedical Text | 10M+ | PubMed, PMC |
| VQA Pairs | 1M+ | VQA-RAD, SLAKE, PathVQA |

### Benchmark Performance

| Dataset | Task | Metric | BiomedGPT | Previous SOTA |
|---------|------|--------|-----------|---------------|
| VQA-RAD | VQA | Accuracy | 84.2% | 79.3% |
| PathVQA | VQA | Accuracy | 87.1% | 82.4% |
| MIMIC-CXR | Report Gen | BLEU-4 | 0.142 | 0.128 |
| PubMedQA | QA | Accuracy | 78.9% | 75.2% |
| ChestX-ray14 | Classification | AUC | 0.831 | 0.815 |

---

## Use Cases

### 1. Radiology Triage

```python
# Automated triage for urgent findings
result = biomedgpt_tool(
    query="Is there any urgent finding requiring immediate attention?",
    image_path="/images/chest_xray.png",
    task="vqa"
)
if "pneumothorax" in result["response"].lower() or "urgent" in result["response"].lower():
    alert_radiologist(priority="HIGH")
```

### 2. Clinical Documentation

```python
# Generate draft reports from images
report = biomedgpt_tool(
    query="65-year-old male with shortness of breath",
    image_path="/images/ct_chest.png",
    task="report"
)
save_draft_report(report["response"])
```

### 3. Medical Education

```python
# Interactive learning tool
explanation = biomedgpt_tool(
    query="Explain the pathophysiology visible in this image",
    image_path="/education/cases/case_001.png",
    task="vqa"
)
```

---

## Related Skills

- `biomedical.radiology.ark_chestxray` - Specialized chest X-ray model
- `biomedical.pathology.conch` - Pathology vision-language model
- `biomedical.clinical.medlm` - Google's MedLM for clinical AI
- `biomedical.foundation_models.nicheformer` - Single-cell foundation model

---

## References

1. **BiomedGPT (2024)**: "BiomedGPT: A Generalist Vision-Language Foundation Model for Diverse Biomedical Tasks." *Nature Medicine*. [arXiv:2305.17100](https://arxiv.org/abs/2305.17100)

2. **GitHub Repository**: [https://github.com/taokz/BiomedGPT](https://github.com/taokz/BiomedGPT)

3. **Model Weights (2025)**: Updated 930M parameter checkpoints. Hugging Face Hub.

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

---

*Last updated: December 2025*

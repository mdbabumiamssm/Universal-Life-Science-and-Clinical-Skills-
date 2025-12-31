# H-optimus-0 Pipeline

**ID:** `biomedical.clinical.pathology.h_optimus_0`
**Version:** 1.0.0
**Status:** Production
**Category:** Clinical / Digital Pathology

---

## Overview

**H-optimus-0** is a state-of-the-art open-source foundation model for pathology, released by Bioptimus. It is a vision transformer pre-trained on hundreds of millions of histology image patches. This skill provides a pipeline to generate slide-level embeddings for downstream tasks like cancer subtyping and survival prediction.

---

## Key Capabilities

- **Feature Extraction:** Converts raw H&E image patches (224x224) into dense vector embeddings.
- **Generalization:** Robust across different tissue types and staining variations.
- **Efficiency:** Optimized for high-throughput WSI (Whole Slide Image) processing.

## Pipeline Steps

1.  **Tiling:** Divide WSI into non-overlapping patches.
2.  **Inference:** Pass patches through H-optimus-0 (ViT-g/14).
3.  **Aggregation:** Pool patch embeddings (e.g., Attention MIL) for patient-level prediction.

## Usage

Hugging Face Hub: `bioptimus/H-optimus-0`

```python
import timm
import torch

model = timm.create_model(
    "hf_hub:bioptimus/H-optimus-0", pretrained=True, init_values=1e-5, dynamic_img_size=False
)
model.eval()
```

## References
- [Hugging Face Model Card](https://huggingface.co/bioptimus/H-optimus-0)

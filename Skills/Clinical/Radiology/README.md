# Radiology Agent

**ID:** `biomedical.clinical.radiology`
**Version:** 0.5.0
**Status:** Experimental
**Category:** Clinical / Imaging

---

## Overview

The **Radiology Agent** utilizes multimodal LLMs (Vision-Language Models) to assist in the interpretation of medical images (X-ray, CT, MRI). It is designed to act as a "second reader," generating preliminary reports and answering visual questions.

## Key Capabilities

- **Report Generation:** Drafts structured radiology reports describing findings, impressions, and recommendations.
- **Visual QA:** Answers specific questions like "Is there a pleural effusion?" or "Locate the mass in the left lung."
- **Modality Support:**
  - Chest X-rays (CXR)
  - Brain MRI (Tumor segmentation/identification)
  - CT Scans (Abdominal/Thoracic)

## Models & Tools

- **RadGPT:** Specializes in conversational radiology report generation.
- **VILA-M3:** Advanced multimodal model for medical imaging.
- **MONAI:** Underlying framework for image transformations and deep learning tasks.

## References
- *RadGPT* (Stanford, 2025)
- *MIMIC-CXR Benchmark*

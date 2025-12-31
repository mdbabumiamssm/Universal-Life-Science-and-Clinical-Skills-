# TIAToolbox Workflow

**ID:** `biomedical.clinical.pathology.tiatoolbox`
**Version:** 1.4.0
**Status:** Production
**Category:** Clinical / Digital Pathology

---

## Overview

**TIAToolbox** (Tissue Image Analytics Toolbox) is a comprehensive Python library developed by the TIA Centre (Warwick) for advanced digital pathology. It handles the "grunt work" of reading large WSI files, stain normalization, and patching, allowing models to focus on inference.

---

## Key Capabilities

- **WSI Reading:** Efficient multi-resolution reading of `.svs`, `.ndpi`, `.tiff`.
- **Stain Normalization:** Macenko, Vahadane methods to correct batch effects.
- **Patch Extraction:** Automated background masking and grid patching.
- **Model Engines:** Wrappers for PyTorch and TensorFlow inference.

## Integration

Used as the **Preprocessing Engine** for foundation models like H-optimus-0 or Prov-GigaPath.

## References
- [TIAToolbox GitHub](https://github.com/TisA-Lab/tiatoolbox)

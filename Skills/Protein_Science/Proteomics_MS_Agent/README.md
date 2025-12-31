# Proteomics Mass Spectrometry Agent

**ID:** `biomedical.protein_science.proteomics_ms`
**Version:** 1.0.0
**Status:** Beta
**Category:** Protein Science / Mass Spectrometry

---

## Overview

The **Proteomics Mass Spectrometry Agent** applies deep learning to mass spectrometry-based proteomics, enabling AI-powered peptide identification, quantification, and biomarker discovery. It addresses the computational challenges of processing high-throughput MS data, including spectral matching, de novo sequencing, and protein inference.

Modern proteomics generates terabytes of spectral data per experiment. This agent leverages state-of-the-art ML models for peptide property prediction, spectral library generation, and data-driven rescoring to maximize proteome coverage and quantification accuracy.

---

## Key Capabilities

### 1. Peptide Property Prediction

| Property | Model Type | Accuracy |
|----------|------------|----------|
| **Retention Time** | Deep learning regression | R² > 0.98 |
| **Fragmentation Pattern** | Transformer (MS²PIP, Prosit) | Cosine > 0.95 |
| **Charge State** | Classification | >95% accuracy |
| **Detectability** | Binary classification | AUC > 0.90 |
| **Collisional Cross Section** | Ion mobility prediction | <3% error |

### 2. Spectral Library Generation

- **In silico digestion:** Trypsin, Lys-C, Asp-N, and custom enzymes
- **Fragment ion prediction:** b/y ions with neutral losses
- **Intensity modeling:** Relative fragment intensities
- **Retention time prediction:** Indexed for targeted analysis

### 3. Database Search Enhancement

- **Percolator-style rescoring:** ML-based FDR control
- **Spectral angle scoring:** Deep learning similarity
- **De novo sequencing:** Transformer-based peptide identification
- **PTM localization:** Phosphorylation site assignment

### 4. Quantification

- **Label-free quantification:** MS1 intensity integration
- **TMT/iTRAQ:** Reporter ion quantification
- **SILAC:** Heavy/light ratio calculation
- **DIA analysis:** Library-based or library-free deconvolution

---

## Usage

### Example Prompt

```text
Analyze this DDA proteomics dataset from a breast cancer study.
Perform database search against UniProt human proteome.
Use deep learning rescoring to improve peptide identification.
Report differentially expressed proteins between tumor and normal samples.
```

### Expected Output

```
## Proteomics Analysis Report

### Dataset Overview
- **Samples:** 6 tumor, 6 normal
- **Instrument:** Orbitrap Exploris 480
- **Acquisition:** DDA, Top20
- **Raw files:** 12 × 2.3 GB

### Search Results (Pre-rescoring)
| Metric | Value |
|--------|-------|
| Total spectra | 1,243,567 |
| PSMs (1% FDR) | 412,890 |
| Peptides identified | 89,234 |
| Proteins identified | 7,892 |

### Deep Learning Rescoring (Prosit + Percolator)
| Metric | Before | After | Improvement |
|--------|--------|-------|-------------|
| PSMs (1% FDR) | 412,890 | 498,234 | +20.7% |
| Peptides | 89,234 | 108,456 | +21.5% |
| Proteins | 7,892 | 9,123 | +15.6% |

### Differential Expression Analysis

**Significantly Upregulated in Tumor (adj. p < 0.05, log2FC > 1):**
| Protein | Gene | log2FC | adj.P |
|---------|------|--------|-------|
| P04626 | ERBB2 | 3.45 | 1.2e-8 |
| P08238 | HSP90AB1 | 2.12 | 3.4e-6 |
| P06493 | CDK1 | 1.89 | 8.7e-5 |

**Pathway Enrichment (Tumor):**
- Cell cycle regulation (p = 1.2e-12)
- DNA replication (p = 4.5e-9)
- PI3K-Akt signaling (p = 2.3e-7)
```

### LLM Agent Integration

```python
@tool
def analyze_proteomics_data(
    raw_files: list[str],
    database: str = "uniprot_human",
    quantification: str = "label_free",
    comparison_groups: dict = None,
    use_ml_rescoring: bool = True
) -> str:
    """
    Analyzes mass spectrometry proteomics data.

    Args:
        raw_files: Paths to raw MS files
        database: Protein database for search
        quantification: LFQ, TMT, SILAC, or DIA
        comparison_groups: Group definitions for differential analysis
        use_ml_rescoring: Apply deep learning rescoring

    Returns:
        Analysis report with identified proteins and statistics
    """
    pass
```

---

## Prerequisites

### Required Software/APIs

| Resource | Purpose | Access |
|----------|---------|--------|
| **Prosit** | Spectral prediction | ProteomeTools API |
| **MS²PIP** | Fragment prediction | CompOmics API |
| **MaxQuant** | Database search | Local installation |
| **MSFragger** | Ultra-fast search | Local installation |
| **Percolator** | Rescoring | Local installation |

### Dependencies

```
pyteomics>=4.6
spectrum_utils>=0.4
pandas>=1.5
numpy>=1.24
torch>=2.0  # For deep learning models
scikit-learn>=1.3
```

---

## Methodology

### Deep Learning for MS

```
Raw Spectra
    ↓
Peak Detection & Centroiding
    ↓
Database Search (MSFragger/MaxQuant)
    ↓
Feature Extraction
    ├── Experimental spectrum features
    ├── Predicted spectrum (Prosit)
    └── Spectral similarity (cosine, SA)
    ↓
ML Rescoring (XGBoost/Percolator)
    ↓
FDR Control (1% PSM, 1% Protein)
    ↓
Quantification
    ↓
Statistical Analysis
```

### Spectral Prediction Models

| Model | Training Data | Architecture | Output |
|-------|---------------|--------------|--------|
| **Prosit** | 100M+ spectra | Transformer | b/y ions, RT |
| **MS²PIP** | 2M spectra | BiLSTM | Fragment intensities |
| **DeepLC** | 4M peptides | CNN | Retention time |
| **pDeep3** | 5M spectra | BiLSTM + Attention | HCD/ETD spectra |

### Rescoring Features

```python
def calculate_rescoring_features(
    experimental_spectrum: np.ndarray,
    predicted_spectrum: np.ndarray,
    peptide_sequence: str
) -> dict:
    """Calculate features for ML rescoring."""
    from spectrum_utils.spectrum import MsmsSpectrum

    features = {
        'spectral_angle': calculate_spectral_angle(exp, pred),
        'pearson_correlation': calculate_correlation(exp, pred),
        'matched_intensity_fraction': matched_intensity / total_intensity,
        'delta_rt': abs(experimental_rt - predicted_rt),
        'peptide_length': len(peptide_sequence),
        'missed_cleavages': count_missed_cleavages(peptide_sequence),
        'charge_state': charge,
        'precursor_mass_error': mass_error_ppm
    }
    return features
```

---

## Related Skills

- **PPI Network Agent:** Analyze protein interaction from co-IP/AP-MS data
- **Data Analysis:** Statistical analysis of quantitative proteomics
- **Biomarker Discovery:** Clinical proteomics applications

---

## References

- **Gessulat et al. (2019):** "Prosit: proteome-wide prediction of peptide tandem mass spectra." *Nature Methods*
- **Degroeve et al. (2021):** "MS²PIP: accurate peptide fragmentation prediction." *Nucleic Acids Research*
- [ProteomeTools](https://proteomicsdb.org/proteomicsdb/prosit/)
- [CompOmics Tools](https://compomics.github.io/)

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

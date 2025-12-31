# Metagenomics Agent

**ID:** `biomedical.genomics.metagenomics`
**Version:** 1.0.0
**Status:** Experimental
**Category:** Genomics / Microbiome

---

## Overview

The **Metagenomics Agent** provides AI-driven workflows for analyzing microbiome data (16S rRNA and Shotgun Metagenomics). It wraps standard bioinformatics tools to automate taxonomic profiling, diversity analysis, and functional inference.

## Key Capabilities

- **Taxonomic Profiling:** Wrappers for `MetaPhlAn 4` and `Kraken2` to identify species abundance.
- **Diversity Analysis:** Calculates Alpha (within-sample) and Beta (between-sample) diversity metrics.
- **Functional Annotation:** Maps reads to metabolic pathways using `HUMAnN 3`.
- **Contamination Control:** Automated host DNA removal (e.g., human reads).

## Integration

Connects with the **Data Analysis Skill** to perform statistical tests (e.g., PERMANOVA) on the resulting abundance tables.

## References
- *The Human Microbiome Project*
- *QIIME 2*

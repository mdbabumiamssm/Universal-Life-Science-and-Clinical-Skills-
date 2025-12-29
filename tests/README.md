# Single-Cell QC Skill: Verification & Demonstration

**Validation Suite for `biomedical.genomics.single_cell_qc`**

---

## 🎯 Purpose

This directory contains the **reference implementation** and **verification suite** for the Single-Cell Quality Control skill. It serves as the "Proving Ground" where the abstract definitions in the Universal Skills Registry are tested against real-world data to ensure:
1.  **Reproducibility:** The code produces identical results on identical inputs.
2.  **Robustness:** The statistical methods (MAD) perform correctly across varying data distributions.
3.  **Integration Readiness:** The Python functions are structured for easy import by AI Agents (MCP, LangChain tools).

## 📂 Components

- **`qc_core.py`**: The **Logic Engine**. Contains the pure functions for calculating metrics and detecting outliers. This file is the "Skill" itself, designed to be stateless and modular.
- **`qc_analysis.py`**: The **CLI Wrapper**. Demonstrates how an agent (or human) would invoke the skill in a production pipeline.
- **`qc_plotting.py`**: The **Visualization Engine**. Generates the visual proof required for scientific trust.
- **`integration_analysis.py`**: A demonstration of downstream utility, showing that QC-cleaned data can be successfully integrated and clustered.

## 🧪 Demonstration Workflow

We have provided a comprehensive test suite (`run_test.sh`) that validates the skill end-to-end:

1.  **Data Ingestion:** Loads raw 10x Genomics or H5AD data.
2.  **Automated QC:** Applies the MAD-based filtering logic defined in the skill specification.
3.  **Metric Validation:** Confirms that known artifacts (e.g., high-MT cells) are correctly identified.
4.  **Downstream Verification:**
    - Concatenates multiple filtered samples (`BM1`, `BM2`, `BM3`).
    - Performs Data Integration (Batch Correction).
    - Generates UMAP projections to confirm that **biological signals are preserved** (clustering by cell type) while **batch effects are minimized**.

## 📊 Results & Artifacts

After running the demonstration, the following artifacts are generated to prove skill efficacy:

- **`analyzed_data/`**: Contains the processed datasets.
    - `*_filtered.h5ad`: The "Gold Standard" clean data.
- **`qc_filtering_thresholds.png`**: Visual proof of the adaptive thresholds (showing exactly where the "cut" was made).
- **`umap_clusters.png`**: Final confirmation that the cleaned data yields distinct biological clusters (cell types).

## 🚀 How to Run

### Prerequisites
- Python 3.9+
- Virtual Environment (recommended)

### Execution
```bash
# 1. Setup Environment
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt

# 2. Run the Full Validation Suite
./run_test.sh

# 3. (Optional) Run Integration Check
python integration_analysis.py
```

## 🔗 Connection to Universal Platform

This implementation serves as the "Source of Truth" for the `biomedical.genomics.single_cell_qc` skill. When this skill is deployed to **Claude** (via MCP) or **ChatGPT**, the underlying logic executed matches exactly what is validated here.
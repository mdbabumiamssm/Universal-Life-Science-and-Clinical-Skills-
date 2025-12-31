# Clinical Trial Matching (TrialGPT Architecture)

**ID:** `biomedical.clinical.trial_matching`
**Version:** 2.0.0 (TrialGPT Upgrade)
**Status:** Production
**Category:** Clinical / Research

---

## Overview

This skill implements the **TrialGPT** architecture for matching patients to clinical trials. It moves beyond simple keyword matching to a sophisticated, multi-stage reasoning process that accurately predicts patient eligibility based on unstructured criteria.

## Architecture

The workflow consists of three specialized agents:

1.  **Criterion Parser Agent:**
    - Deconstructs complex inclusion/exclusion criteria from ClinicalTrials.gov free text.
    - Standardizes logic (e.g., "Age >= 18", "No history of cardiac events").

2.  **Patient Data Retriever:**
    - Extracts relevant values from the patient's EHR summary (e.g., "LVEF = 45%", "Stage IV NSCLC").

3.  **Eligibility Reasoner (Ranking):**
    - Compares patient data against parsed criteria.
    - Assigns a score:
        - **Eligible**
        - **Excluded**
        - **Uncertain** (Requires more tests)
    - Generates a ranked list of trials with a "Why this match?" explanation.

## Capabilities

- **Cohort Discovery:** "Find all trials for Stage III Breast Cancer recruiting in New York."
- **Pre-screening:** Automated first-pass screening of patient lists against a specific protocol.
- **Geography Aware:** Filters by trial site distance from patient.

## References
- *TrialGPT: Large Language Models for Clinical Trial Matching* (NIH, 2025)
- *TrialMatchAI*

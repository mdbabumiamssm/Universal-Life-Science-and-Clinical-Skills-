# TrialGPT: Clinical Trial Matching Agent

## Overview
This skill implements a simplified "TrialGPT" agent that matches patient profiles to clinical trials using either keyword-based similarity or semantic search (if `sentence_transformers` is installed).

## Features
- **Patient Profiling:** Ingests patient age, conditions, and notes.
- **Trial Database:** Loads trial data from a JSON source.
- **Matching Engine:**
  - **Keyword Mode:** Uses Jaccard similarity on tokens (default, requires no heavy libs).
  - **LLM Mode:** Uses `sentence_transformers` to compute semantic similarity between patient description and trial text.
- **Filtering:** Basic rule-based filtering (e.g., Age) to demonstrate hybrid agent logic.

## Usage

### Prerequisites
For basic mode, only standard Python libraries are needed.
For LLM mode:
```bash
pip install sentence-transformers numpy
```

### Running the Agent

**1. Basic Match (Keyword):**
```bash
python trial_matcher.py --age 65 --condition "Lung Cancer" --notes "Failed previous chemo"
```

**2. Semantic Match (LLM):**
```bash
python trial_matcher.py --age 65 --condition "Lung Cancer" --notes "Failed previous chemo" --use-llm
```

### Custom Data
Edit `dummy_data.json` to add more clinical trials to the local database.

## Architecture
The agent uses a "Retriever-Ranker" inspired approach where trials are first filtered by hard constraints (Age) and then ranked by semantic relevance.

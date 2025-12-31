# MedLM Clinical Agent

**ID:** `biomedical.foundation_models.medlm`
**Version:** 1.0.0
**Status:** Production
**Category:** Foundation Models / Clinical Healthcare AI
**Released:** December 2025

---

## Overview

The **MedLM Clinical Agent** provides access to Google's MedLM - a family of foundation models fine-tuned for the healthcare industry. MedLM is powered by Med-PaLM 2, which was the first AI system to achieve passing scores on medical licensing examinations and is now available through Google Cloud's Vertex AI.

### Key Achievements

| Milestone | Performance | Significance |
|-----------|-------------|--------------|
| USMLE Pass | >60% | First AI to pass medical licensing |
| MedQA Score | 86.5% | +19% over Med-PaLM 1 |
| Physician Preference | 8/9 axes | Preferred over physician answers |

MedLM comes in two versions optimized for different healthcare use cases:
- **MedLM-large**: Complex medical tasks requiring deep reasoning
- **MedLM-medium**: Scalable tasks with lower latency requirements

---

## Key Capabilities

### 1. Clinical Tasks

| Task | Description | Performance |
|------|-------------|-------------|
| **Medical Q&A** | Answer clinical questions | Expert-level |
| **Diagnostic Reasoning** | Differential diagnosis | High accuracy |
| **Treatment Planning** | Evidence-based recommendations | Guideline-aligned |
| **Clinical Summarization** | Condense medical records | Efficient |
| **Patient Communication** | Health information in plain language | Accessible |
| **Medical Documentation** | Draft clinical notes | Time-saving |

### 2. Deployment Options

| Mode | Latency | Best For |
|------|---------|----------|
| MedLM-large | Higher | Complex reasoning, diagnostics |
| MedLM-medium | Lower | High-volume, real-time |
| Custom fine-tuned | Variable | Domain-specific needs |

### 3. Real-World Implementations

| Organization | Use Case | Outcome |
|--------------|----------|---------|
| **HCA Healthcare** | ED note drafting | Physician time savings |
| **Augmedix** | Ambient documentation | Automated clinical notes |
| **Research Institutions** | Medical research | Literature synthesis |

---

## Usage

### Example Prompt

```text
A 58-year-old male presents to the emergency department with:
- Chief complaint: Sudden onset crushing chest pain for 2 hours
- Pain radiates to left arm and jaw
- Associated diaphoresis and nausea
- History: HTN, T2DM, hyperlipidemia, 30 pack-year smoking
- Medications: Metformin, lisinopril, atorvastatin
- Vitals: BP 165/95, HR 102, RR 22, SpO2 94% on RA

Please provide:
1. Most likely diagnosis and differentials
2. Immediate workup and management
3. Risk stratification
4. Disposition recommendation
```

### Expected Output

```markdown
## MedLM Clinical Assessment

### Primary Assessment

**Most Likely Diagnosis**: Acute Coronary Syndrome (ACS) - STEMI vs NSTEMI

**Clinical Reasoning**:
- Classic presentation: Crushing substernal chest pain with radiation
- High-risk features: Diaphoresis, duration >20 minutes
- Multiple cardiovascular risk factors (HTN, DM, HLD, smoking)
- Hemodynamic stress response (tachycardia, hypertension)

### Differential Diagnosis

| Diagnosis | Probability | Key Features |
|-----------|-------------|--------------|
| **STEMI/NSTEMI** | 85% | Classic presentation, risk factors |
| Unstable Angina | 10% | If biomarkers negative |
| Aortic Dissection | 3% | Requires BP differential, CXR |
| Pulmonary Embolism | 1% | Hypoxia present, needs consideration |
| Esophageal Rupture | <1% | No emesis history |

### Immediate Management

#### Time-Critical Actions (Door-to-ECG <10 min)
```
1. 12-lead ECG - STAT
2. IV access x2 large bore
3. Continuous cardiac monitoring
4. Supplemental O2 to maintain SpO2 >94%
```

#### Initial Workup
| Test | Rationale | Timing |
|------|-----------|--------|
| Troponin I/T | ACS confirmation | STAT, repeat 3-6h |
| BMP | Renal function for contrast | STAT |
| CBC | Baseline, bleeding risk | STAT |
| Coagulation | Before anticoagulation | STAT |
| CXR portable | Cardiomegaly, CHF, dissection | STAT |
| Type & Screen | If intervention needed | STAT |

#### Pharmacotherapy
```
MONA-B Protocol (if no contraindications):
□ Morphine 2-4mg IV PRN severe pain
□ Oxygen if SpO2 <94%
□ Nitrates: NTG 0.4mg SL q5min x3 (hold if SBP<90)
□ Aspirin 325mg chewed (if not already given)
□ Beta-blocker: Metoprolol 5mg IV (if HR>60, SBP>100)

Additional:
□ Heparin: UFH 60 U/kg bolus + 12 U/kg/hr OR
□ Enoxaparin 1mg/kg SC (if CrCl >30)
□ P2Y12 inhibitor: Loading dose per cardiology
```

### Risk Stratification

#### HEART Score Calculation
| Component | Value | Points |
|-----------|-------|--------|
| History | Highly suspicious | 2 |
| ECG | Pending | 0-2 |
| Age | 58 | 1 |
| Risk factors | ≥3 | 2 |
| Troponin | Pending | 0-2 |
| **Estimated Total** | **5-9** | **High Risk** |

#### TIMI Risk Score (if NSTEMI)
- Age ≥65: 0
- ≥3 CAD risk factors: 1
- Known CAD: 0 (unknown)
- ASA use in past 7 days: 0
- Severe angina: 1
- ST changes: Pending
- Elevated biomarkers: Pending
- **Preliminary: 2-4 (Moderate-High)**

### Disposition Recommendation

**IF STEMI on ECG:**
```
→ Activate Cardiac Cath Lab
→ Door-to-balloon goal <90 minutes
→ Cardiology STAT consult
→ Direct to CCU post-procedure
```

**IF NSTEMI/UA:**
```
→ Admit to Telemetry or CCU (based on hemodynamics)
→ Cardiology consult within 24 hours
→ Risk-stratify for early invasive vs conservative strategy
→ High-risk features → Cardiac cath within 24 hours
```

### Documentation Note (Draft)

> **EMERGENCY DEPARTMENT NOTE**
>
> **CC**: Chest pain x 2 hours
>
> **HPI**: 58 y/o M with PMH HTN, T2DM, HLD, 30 pack-year smoking
> presents with sudden onset crushing substernal chest pain radiating
> to L arm and jaw, associated with diaphoresis and nausea. Pain
> began 2 hours PTA while at rest. No prior similar episodes.
>
> **Assessment**: Acute coronary syndrome, high-risk presentation
>
> **Plan**: Emergent cardiology consultation, serial troponins, risk
> stratification, anticoagulation per ACS protocol.

**Confidence**: High (92%) - Classic ACS presentation
**Evidence Level**: Guidelines-based (ACC/AHA 2023)
```

---

## LLM Agent Integration

### Python Tool Implementation

```python
from typing import Optional, Dict, Any, List
from google.cloud import aiplatform
from vertexai.preview.language_models import TextGenerationModel

def medlm_clinical_tool(
    query: str,
    patient_context: Optional[str] = None,
    task_type: str = "diagnosis",
    model_variant: str = "large",
    temperature: float = 0.2,
    max_tokens: int = 2048,
    safety_filter: bool = True
) -> Dict[str, Any]:
    """
    MedLM clinical foundation model for healthcare AI tasks.

    Args:
        query: Clinical question or instruction
        patient_context: Optional structured patient data
        task_type: 'diagnosis', 'treatment', 'documentation', 'education', 'qa'
        model_variant: 'large' (complex) or 'medium' (fast)
        temperature: Response randomness (lower = more deterministic)
        max_tokens: Maximum response length
        safety_filter: Apply medical safety guardrails

    Returns:
        Dictionary with clinical response and metadata
    """
    # Initialize Vertex AI
    aiplatform.init(project="your-gcp-project", location="us-central1")

    # Select model
    model_name = f"medlm@{model_variant}"
    model = TextGenerationModel.from_pretrained(model_name)

    # Build clinical prompt
    system_prompt = """You are MedLM, a clinical AI assistant trained to support
healthcare professionals. Provide evidence-based medical information while:
- Citing clinical guidelines when available
- Acknowledging uncertainty appropriately
- Recommending specialist consultation when needed
- Never replacing clinical judgment"""

    # Task-specific formatting
    task_prompts = {
        "diagnosis": f"""Based on the following clinical presentation, provide:
1. Most likely diagnosis with reasoning
2. Differential diagnoses ranked by probability
3. Recommended diagnostic workup
4. Red flags requiring immediate attention

Patient Information:
{patient_context if patient_context else 'Not provided'}

Clinical Question:
{query}""",

        "treatment": f"""Provide evidence-based treatment recommendations for:

Patient Context:
{patient_context if patient_context else 'Not provided'}

Clinical Question:
{query}

Include:
- First-line therapy with dosing
- Alternative options
- Contraindications to consider
- Monitoring parameters""",

        "documentation": f"""Draft clinical documentation for:

Context:
{patient_context if patient_context else 'Not provided'}

Requirements:
{query}

Follow standard medical documentation format.""",

        "education": f"""Explain the following in patient-friendly language:

Topic: {query}

Context: {patient_context if patient_context else 'General audience'}

Make it accessible while maintaining accuracy.""",

        "qa": f"""Answer the following clinical question:

{query}

Context: {patient_context if patient_context else 'None provided'}

Provide a thorough, evidence-based response."""
    }

    full_prompt = f"{system_prompt}\n\n{task_prompts.get(task_type, task_prompts['qa'])}"

    # Generate response
    response = model.predict(
        full_prompt,
        temperature=temperature,
        max_output_tokens=max_tokens,
        top_k=40,
        top_p=0.95
    )

    result = {
        "response": response.text,
        "task_type": task_type,
        "model": model_name,
        "safety_filtered": safety_filter
    }

    # Add safety checks if enabled
    if safety_filter:
        result["safety_check"] = _check_medical_safety(response.text)

    return result


def _check_medical_safety(response: str) -> Dict[str, Any]:
    """Check response for medical safety concerns."""
    red_flags = [
        "prescribe", "diagnose definitively", "stop taking",
        "ignore your doctor", "cure"
    ]

    concerns = []
    for flag in red_flags:
        if flag.lower() in response.lower():
            concerns.append(flag)

    return {
        "passed": len(concerns) == 0,
        "concerns": concerns,
        "recommendation": "Human review recommended" if concerns else "OK"
    }


# Claude/Anthropic tool schema
MEDLM_TOOL_SCHEMA = {
    "name": "medlm_clinical_assistant",
    "description": "Google's MedLM clinical AI for medical question answering, diagnostic reasoning, treatment planning, and clinical documentation",
    "input_schema": {
        "type": "object",
        "properties": {
            "query": {
                "type": "string",
                "description": "Clinical question or task description"
            },
            "patient_context": {
                "type": "string",
                "description": "Relevant patient information (demographics, history, vitals, labs)"
            },
            "task_type": {
                "type": "string",
                "enum": ["diagnosis", "treatment", "documentation", "education", "qa"],
                "description": "Type of clinical task"
            }
        },
        "required": ["query"]
    }
}
```

### AMIE-Style Diagnostic Dialogue

```python
class AMIEDiagnosticAgent:
    """
    AMIE (Articulate Medical Intelligence Explorer) style
    diagnostic conversation agent.
    """

    def __init__(self, model_variant: str = "large"):
        self.model = TextGenerationModel.from_pretrained(f"medlm@{model_variant}")
        self.conversation_history = []
        self.gathered_info = {
            "chief_complaint": None,
            "hpi": [],
            "pmh": [],
            "medications": [],
            "allergies": [],
            "social_history": [],
            "family_history": [],
            "review_of_systems": [],
            "vitals": None,
            "physical_exam": []
        }

    def conduct_interview(self, patient_message: str) -> str:
        """Process patient response and generate next question."""
        self.conversation_history.append({
            "role": "patient",
            "content": patient_message
        })

        # Extract information from response
        self._extract_clinical_info(patient_message)

        # Generate next question or assessment
        prompt = self._build_interview_prompt()
        response = self.model.predict(prompt, temperature=0.3)

        self.conversation_history.append({
            "role": "physician",
            "content": response.text
        })

        return response.text

    def _extract_clinical_info(self, message: str):
        """Extract structured clinical info from free text."""
        # Use NER/IE to populate gathered_info
        pass

    def _build_interview_prompt(self) -> str:
        """Build prompt for next interview step."""
        info_status = self._assess_information_completeness()

        return f"""You are conducting a diagnostic medical interview.

Information gathered so far:
{self.gathered_info}

Completeness assessment:
{info_status}

Conversation history:
{self.conversation_history[-6:]}  # Last 6 turns

Based on the information gathered, either:
1. Ask the most clinically relevant follow-up question
2. If sufficient information, provide preliminary assessment

Respond naturally as a caring physician would."""

    def generate_assessment(self) -> Dict[str, Any]:
        """Generate final diagnostic assessment."""
        prompt = f"""Based on the complete clinical interview:

{self.gathered_info}

Provide:
1. Most likely diagnosis with probability
2. Differential diagnoses
3. Recommended next steps
4. Any urgent concerns"""

        response = self.model.predict(prompt, temperature=0.2)
        return {
            "assessment": response.text,
            "gathered_info": self.gathered_info,
            "conversation_turns": len(self.conversation_history)
        }
```

---

## Prerequisites

### Google Cloud Setup

```bash
# Install Google Cloud SDK
curl https://sdk.cloud.google.com | bash

# Authenticate
gcloud auth login
gcloud auth application-default login

# Enable APIs
gcloud services enable aiplatform.googleapis.com

# Install Python SDK
pip install google-cloud-aiplatform vertexai
```

### Required Permissions

| Permission | Purpose |
|------------|---------|
| `aiplatform.endpoints.predict` | Model inference |
| `aiplatform.models.get` | Model access |
| `healthcare.datasets.*` | FHIR data integration |

### Compliance Requirements

| Requirement | Status |
|-------------|--------|
| HIPAA BAA | Required for PHI |
| SOC 2 Type II | Google Cloud certified |
| HITRUST | Available |
| Regional Data Residency | US regions available |

---

## Methodology

### Model Training

```
Med-PaLM 2 / MedLM Training Pipeline
─────────────────────────────────────

Base Model: PaLM 2
     │
     ▼
Medical Pretraining
├── PubMed abstracts (35M+)
├── Medical textbooks
├── Clinical guidelines
└── Licensed medical content
     │
     ▼
Instruction Fine-tuning
├── Medical QA datasets
├── Clinical reasoning tasks
├── Expert demonstrations
└── RLHF with physician feedback
     │
     ▼
Safety Alignment
├── Harm avoidance training
├── Uncertainty calibration
├── Deferral to experts
└── Factuality grounding
     │
     ▼
MedLM (Production)
```

### Evaluation Framework

| Benchmark | Task | MedLM Score |
|-----------|------|-------------|
| MedQA (USMLE) | Multiple choice | 86.5% |
| PubMedQA | Research QA | 79.2% |
| MedMCQA | Indian medical | 72.3% |
| LiveQA | Consumer health | 82.1% |
| MedicationQA | Drug info | 88.4% |

---

## Safety & Limitations

### Built-in Safeguards

1. **Uncertainty Expression**: Model trained to express confidence levels
2. **Deferral Triggers**: Recommends specialist consultation appropriately
3. **Scope Limitations**: Acknowledges what it cannot do
4. **Hallucination Mitigation**: Grounded in medical literature

### Known Limitations

| Limitation | Mitigation |
|------------|------------|
| Not a replacement for physicians | Clear disclaimers |
| May not reflect latest guidelines | Regular updates |
| Rare disease knowledge gaps | Specialist referral |
| Cannot perform physical exam | Explicit acknowledgment |

### Responsible Use Guidelines

```markdown
DO:
- Use as clinical decision support
- Verify recommendations against guidelines
- Maintain human oversight
- Document AI assistance in records

DON'T:
- Use for autonomous diagnosis
- Share PHI without proper safeguards
- Rely solely on AI for critical decisions
- Use for conditions outside training scope
```

---

## Related Skills

- `biomedical.clinical.ehr_llm` - EHR-specific language models
- `biomedical.foundation_models.biomedgpt` - Multimodal biomedical AI
- `biomedical.clinical.trial_eligibility` - Clinical trial matching
- `biomedical.clinical.note_summarization` - Clinical documentation

---

## References

1. **Med-PaLM 2 (2023)**: "Large language models encode clinical knowledge." *Nature*. [DOI: 10.1038/s41586-023-06291-2](https://www.nature.com/articles/s41586-023-06291-2)

2. **MedLM (2024)**: "Introducing MedLM for the healthcare industry." *Google Cloud Blog*.

3. **AMIE (2024)**: "Towards Expert-Level Medical Question Answering with Large Language Models." *arXiv*.

4. **Google Health AI**: [https://health.google/health-research/](https://health.google/health-research/)

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

---

*Last updated: December 2025*

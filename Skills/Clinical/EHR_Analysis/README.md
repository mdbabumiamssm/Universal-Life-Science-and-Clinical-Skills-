# EHR Analysis Agent

**ID:** `biomedical.clinical.ehr_analysis`
**Version:** 1.0.0
**Status:** Beta
**Category:** Clinical / Informatics

---

## Overview

The **EHR Analysis Agent** is designed to interact with Electronic Health Records (EHR) in a secure and structured manner. It moves beyond simple text summarization to perform complex querying and risk stratification on patient data.

## Key Capabilities

### 1. Structured Data Querying (Text-to-SQL)
- Translates natural language questions ("How many patients with Diabetes Type 2 have HbA1c > 9?") into SQL queries for OMOP/FHIR databases.

### 2. Timeline Reconstruction
- Aggregates notes, labs, and medications to create a chronological view of a patient's disease progression.

### 3. Risk Scoring
- Calculates standard clinical scores (e.g., HAS-BLED, CHA2DS2-VASc) automatically from extracted data variables.

## Integration
- **FHIR:** Uses HL7 FHIR standards for interoperability.
- **OMOP CDM:** Compatible with Observational Medical Outcomes Partnership Common Data Model.

## References
- *ChatEHR* (Stanford)
- *Med-PaLM 2* (Google)

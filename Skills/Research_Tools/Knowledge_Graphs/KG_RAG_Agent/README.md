# Biomedical Knowledge Graph RAG Agent

**ID:** `biomedical.research_tools.kg_rag`
**Version:** 1.0.0
**Status:** Experimental
**Category:** Research Tools / Neuro-symbolic AI

---

## Overview

The **Biomedical Knowledge Graph RAG Agent** enhances traditional Retrieval-Augmented Generation (RAG) by grounding LLM reasoning in structured biomedical ontologies (e.g., UMLS, SNOMED CT, PrimeKG).

Unlike vector-only RAG, which relies on semantic similarity, KG-RAG enables **multi-hop reasoning** (e.g., "Drug A treats Disease B, which is caused by Gene C") and reduces hallucinations by verifying facts against the graph.

---

## Architecture

1.  **Graph Database:** Neo4j (storing nodes: Drug, Disease, Gene, Pathway).
2.  **Orchestrator:** LangChain / LangGraph.
3.  **Query Generator:** LLM translates natural language to **Cypher** queries.
4.  **Hybrid Retrieval:** Combines Vector Search (unstructured text) + Graph Traversal (structured facts).

## Capabilities

- **Entity Grounding:** Maps user terms ("heart attack") to canonical IDs (C0027051).
- **Path Finding:** "Find all shortest paths between Drug X and Gene Y."
- **Fact Verification:** Cross-checks generated answers against known graph edges.

## Prerequisites

- Neo4j Database (Community or AuraDB)
- `langchain-community`, `neo4j` Python drivers.
- Access to PrimeKG or similar dataset.

## Author
**MD BABU MIA**
*Artificial Intelligence Group*

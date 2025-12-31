# mRNA Design Agent

**ID:** `biomedical.drug_discovery.mrna_design`
**Version:** 1.0.0
**Status:** Production
**Category:** Drug Discovery / RNA Therapeutics

---

## Overview

The **mRNA Design Agent** optimizes messenger RNA sequences for therapeutic applications, including vaccines and protein replacement therapies. It addresses the multi-objective optimization challenges of mRNA design: maximizing translation efficiency, minimizing immunogenicity, and ensuring structural stability.

The mRNA therapeutics market has exploded following COVID-19 vaccine success. This agent leverages AI approaches like Ginkgo Bioworks' mDD-0 diffusion model to generate complete mRNA transcripts including optimized 5' UTR, coding sequence, and 3' UTR.

---

## Key Capabilities

### 1. Sequence Optimization

| Region | Optimization Target | Method |
|--------|---------------------|--------|
| **5' UTR** | Translation initiation efficiency | Kozak sequence optimization, secondary structure minimization |
| **Coding Sequence** | Codon usage, stability | Codon optimization for host organism |
| **3' UTR** | Stability, poly(A) efficiency | Regulatory element selection |
| **Full Transcript** | Overall expression | End-to-end diffusion models |

### 2. Property Prediction

- **Translation efficiency:** Ribosome loading and elongation rate
- **mRNA half-life:** Stability predictions
- **Immunogenicity:** TLR activation risk (TLR7/8, RIG-I)
- **Secondary structure:** MFE (Minimum Free Energy) calculations
- **GC content:** Optimal range for stability

### 3. Codon Optimization Strategies

| Strategy | Description | Use Case |
|----------|-------------|----------|
| **CAI Maximization** | Host codon preference | High expression |
| **tRNA Adaptation** | Match tRNA pools | Specific cell types |
| **GC Optimization** | Target 50-60% GC | Stability |
| **Rare Codon Avoidance** | Remove slow codons | Faster translation |
| **Codon Pair Bias** | Optimize dinucleotides | Natural-like sequences |

### 4. Modified Nucleotide Support

- **N1-methylpseudouridine (m1Ψ):** Reduced immunogenicity
- **5-methylcytidine (m5C):** Enhanced stability
- **2-thiouridine (s2U):** Altered base pairing
- **Pseudouridine (Ψ):** Standard modification

---

## Usage

### Example Prompt

```text
Design an optimized mRNA sequence for human Factor IX expression.
Target: Liver hepatocytes
Requirements:
- High translation efficiency
- Low immunogenicity (use m1Ψ compatible sequence)
- Minimize secondary structure in 5' UTR
- Include alpha-globin 3' UTR

Provide codon-optimized sequence and predicted expression level.
```

### Expected Output

```
## mRNA Design Report: Factor IX

### Target Protein
- **Gene:** F9 (Coagulation Factor IX)
- **UniProt:** P00740
- **CDS Length:** 1,383 bp (461 amino acids)

### Optimized Sequence Design

#### 5' UTR (47 nt)
```
GGGAAAUAAGAGAGAAAAGAAGAGUAAGAAGAAAUAUAAGAGCCACC
```
- Kozak strength: Strong (GCCACC)
- Predicted MFE: -2.3 kcal/mol (minimal structure)
- Cap-proximal sequence optimized for eIF4E binding

#### Optimized CDS (1,383 nt)
- **CAI (human liver):** 0.89 (original: 0.71)
- **GC content:** 54.2%
- **Rare codons removed:** 23 → 0
- **CpG dinucleotides:** Minimized (immunogenicity)

| Metric | Original | Optimized |
|--------|----------|-----------|
| CAI | 0.71 | 0.89 |
| GC% | 41.2% | 54.2% |
| MFE (kcal/mol) | -412 | -389 |
| Predicted expression | 1.0x | 4.2x |

#### 3' UTR (142 nt)
- Human alpha-globin 3' UTR
- Poly(A) signal: AAUAAA (canonical)
- Poly(A) tail: 120 nt (recommended)

### Immunogenicity Assessment
| Sensor | Risk Score | Notes |
|--------|------------|-------|
| TLR7/8 | Low | GU-rich motifs minimized |
| RIG-I | Low | 5' triphosphate capped |
| PKR | Low | dsRNA regions avoided |

### Delivery Recommendations
- **LNP formulation:** Standard ionizable lipid
- **Dose range:** 0.1-0.5 mg/kg (based on similar mRNAs)
- **Target tissue:** Liver (IV administration)

### Output Files
- `F9_optimized_mrna.fasta` - Full sequence
- `secondary_structure.ct` - Predicted structure
- `expression_prediction.json` - Detailed metrics
```

### LLM Agent Integration

```python
@tool
def design_mrna(
    protein_sequence: str,
    target_organism: str = "human",
    target_tissue: str = "general",
    optimization_strategy: str = "balanced",
    include_utrs: bool = True,
    modified_nucleotides: str = "m1psi"
) -> str:
    """
    Designs optimized mRNA for protein expression.

    Args:
        protein_sequence: Amino acid sequence or UniProt ID
        target_organism: Expression host (human, mouse)
        target_tissue: Target tissue for codon optimization
        optimization_strategy: balanced, high_expression, low_immunogenicity
        include_utrs: Include optimized UTRs
        modified_nucleotides: Modification type (m1psi, standard)

    Returns:
        Optimized mRNA sequence with predicted properties
    """
    pass
```

---

## Prerequisites

### Required Tools/APIs

| Resource | Purpose | Access |
|----------|---------|--------|
| **LinearDesign** | mRNA structure optimization | Baidu Research |
| **ViennaRNA** | Secondary structure prediction | Local installation |
| **Codon Usage DB** | Organism-specific tables | Kazusa API |
| **UTR Database** | Regulatory element library | UTRdb |

### Dependencies

```
ViennaRNA>=2.5
biopython>=1.81
pandas>=1.5
numpy>=1.24
torch>=2.0  # For ML models
```

---

## Methodology

### mRNA Design Pipeline

```
Input: Protein Sequence
    ↓
Reverse Translation (all synonymous codons)
    ↓
Codon Optimization
    ├── CAI maximization
    ├── GC content adjustment
    ├── Rare codon removal
    └── CpG depletion (immunogenicity)
    ↓
5' UTR Design
    ├── Kozak sequence optimization
    └── Structure minimization
    ↓
3' UTR Selection
    ├── Stability elements
    └── Poly(A) signal
    ↓
Full Sequence Assembly
    ↓
Secondary Structure Prediction
    ↓
Expression & Immunogenicity Scoring
    ↓
Output: Optimized mRNA
```

### Optimization Algorithms

```python
def optimize_codon_usage(
    protein_seq: str,
    target_organism: str,
    gc_target: float = 0.55
) -> str:
    """
    Multi-objective codon optimization using genetic algorithm.
    """
    from Bio.SeqUtils import CodonUsage

    # Initialize with random synonymous codons
    population = initialize_population(protein_seq, size=100)

    for generation in range(500):
        # Evaluate fitness
        fitness = [
            calculate_fitness(seq, target_organism, gc_target)
            for seq in population
        ]

        # Selection, crossover, mutation
        population = evolve_population(population, fitness)

    return max(population, key=lambda x: calculate_fitness(x))


def calculate_fitness(cds: str, organism: str, gc_target: float) -> float:
    """Multi-objective fitness function."""
    cai = calculate_cai(cds, organism)
    gc = calculate_gc_content(cds)
    mfe = calculate_mfe(cds)
    cpg = count_cpg(cds) / len(cds)

    # Weighted combination
    fitness = (
        0.4 * cai +
        0.2 * (1 - abs(gc - gc_target)) +
        0.2 * (1 - normalize_mfe(mfe)) +
        0.2 * (1 - cpg * 10)  # Penalize CpG
    )
    return fitness
```

---

## Related Skills

- **siRNA/ASO Design Agent:** Complementary RNA therapeutic modalities
- **AgentD Drug Discovery:** Small molecule alternatives
- **Neoantigen Vaccine Agent:** mRNA vaccine applications

---

## References

- **Zhang et al. (2023):** "Algorithm for optimized mRNA design improves stability and immunogenicity." *Nature*
- **Leppek et al. (2022):** "Combinatorial optimization of mRNA structure, stability, and translation for RNA-based therapeutics." *Nature Communications*
- [LinearDesign](https://github.com/LinearDesignSoftware/LinearDesign)
- [mRNA Therapeutics Review](https://www.nature.com/articles/s41573-021-00283-5)

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

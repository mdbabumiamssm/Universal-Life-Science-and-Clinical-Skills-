# Synthetic sRNA Regulator Design Agent

**ID:** `biomedical.synthetic_biology.srna_designer`
**Version:** 1.0.0
**Status:** Beta
**Category:** Synthetic Biology / Gene Regulation

---

## Overview

The **Synthetic sRNA Regulator Design Agent** designs synthetic small regulatory RNAs (sRNAs) for precise gene expression control in bacteria. Unlike transcription factors that exhibit sigmoidal responses, sRNAs display linear response curves, making them ideal for fine-tuning gene expression in metabolic engineering and high-throughput screening applications.

This agent designs sRNA sequences targeting specific mRNAs, predicts knockdown efficiency, and optimizes sRNA scaffold for stability and activity.

---

## Key Capabilities

### 1. sRNA Design Types

| Type | Mechanism | Application |
|------|-----------|-------------|
| **Antisense sRNA** | Base-pairing repression | Gene knockdown |
| **Riboregulator** | Inducible activation | Conditional expression |
| **Toehold switch** | Programmable activation | Biosensing |
| **STAR (Small Trans-Activating RNA)** | Transcriptional activation | Orthogonal control |

### 2. Target Site Selection

- **5' UTR targeting:** Blocks ribosome binding
- **Start codon region:** Optimal knockdown
- **Hfq binding:** Enhanced stability
- **Secondary structure:** Accessible sites

### 3. Scaffold Options

| Scaffold | Source | Stability | Hfq Binding |
|----------|--------|-----------|-------------|
| **MicC** | E. coli | High | Strong |
| **MicF** | E. coli | High | Strong |
| **SgrS** | E. coli | Medium | Moderate |
| **Spot42** | E. coli | Medium | Weak |

### 4. Expression Systems

- **Constitutive:** For fixed knockdown
- **Inducible (IPTG, arabinose):** Tunable knockdown
- **CRISPRi hybrid:** Combined regulation

---

## Usage

### Example Prompt

```text
Design synthetic sRNAs to knockdown competing pathway genes for
lycopene production in E. coli.

Target genes:
1. gdhA (glutamate dehydrogenase) - competes for NADPH
2. ilvE (branched-chain aminotransferase) - competes for pyruvate
3. gltA (citrate synthase) - drains acetyl-CoA

Requirements:
- >70% knockdown efficiency
- Tunable expression (IPTG-inducible)
- Orthogonal (no cross-reactivity)
```

### Expected Output

```
## Synthetic sRNA Design Report

### Design Summary
| Target | Gene | Knockdown | Orthogonality |
|--------|------|-----------|---------------|
| gdhA | b1761 | 82% | ✓ |
| ilvE | b3770 | 78% | ✓ |
| gltA | b0720 | 85% | ✓ |

### sRNA-gdhA Design

#### Target Site Analysis
```
gdhA mRNA 5' UTR + start:
5'-...UAAGGAGGUGAUUAUUAUGAACAUUA...-3'
         [  RBS  ]   [AUG]
                  ↑
            Target region (-20 to +10)
```

#### sRNA Sequence
```
>sRNA-gdhA (84 nt)
5'-AAUAAUGUUCAUAAUCACCUCCUUAAGAUUUUUGAUCUUUCUUGU
    CUGGUGUUCAGUAUUAAAUGCAUCGUGAAAAAAUGAUCC-3'

Antisense (24 nt): AAUAAUGUUCAUAAUCACCUCCUU
MicC scaffold (60 nt): AAGAUUUUUGAU...AAUGAUCC
```

#### Structure Prediction
```
        Antisense region
       ╭─────────────────╮
5'-AAUAAUGUUCAUAAUCACCUCCUU
                          ║
                  MicC scaffold
                          ║
              Hfq binding site
                          ║
   3'-CCUAGUAAAAGUGCUACGUA...
```

#### Predicted Knockdown
| Metric | Value | Confidence |
|--------|-------|------------|
| **Binding ΔG** | -38.2 kcal/mol | - |
| **Knockdown** | 82% ± 8% | High |
| **Half-life** | 45 min | - |

### sRNA-ilvE Design

#### Target Site
```
ilvE mRNA 5' UTR + start:
5'-...GAGGAAUUUAUAAUGAUUUGAACG...-3'
       [ RBS ]   [AUG]
```

#### sRNA Sequence
```
>sRNA-ilvE (86 nt)
5'-CGUUCAAAUCAUAUAAAUUCCUCAAGAUUUUUGAUCUUUCUUGU
    CUGGUGUUCAGUAUUAAAUGCAUCGUGAAAAAAUGAUCC-3'

Antisense (24 nt): CGUUCAAAUCAUAUAAAUUCCUC
```

#### Predicted Knockdown: 78% ± 10%

### sRNA-gltA Design

#### Target Site
```
gltA mRNA 5' UTR + start:
5'-...AGGAGAUAUUAUGCUGAUACA...-3'
       [ RBS ][AUG]
```

#### sRNA Sequence
```
>sRNA-gltA (82 nt)
5'-UGUAUCAGCAUAAUAUCUCCU AAGAUUUUUGAUCUUUCUUGU
    CUGGUGUUCAGUAUUAAAUGCAUCGUGAAAAAAUGAUCC-3'

Antisense (21 nt): UGUAUCAGCAUAAUAUCUCCU
```

#### Predicted Knockdown: 85% ± 7%

### Orthogonality Matrix

| sRNA | gdhA | ilvE | gltA | Off-targets |
|------|------|------|------|-------------|
| sRNA-gdhA | **82%** | <5% | <5% | 0 |
| sRNA-ilvE | <5% | **78%** | <5% | 0 |
| sRNA-gltA | <5% | <5% | **85%** | 0 |

### Expression Plasmid Design

```
pSRNA-triple (IPTG-inducible)
├── Plac-lacO → sRNA-gdhA → T7te
├── Plac-lacO → sRNA-ilvE → T7te
└── Plac-lacO → sRNA-gltA → T7te
    + lacI (constitutive)
```

**Backbone:** pSC101 origin (low copy, ~5/cell)
**Selection:** KanR

### Dose-Response Prediction

| IPTG (mM) | gdhA Expression | ilvE Expression | gltA Expression |
|-----------|-----------------|-----------------|-----------------|
| 0 | 100% | 100% | 100% |
| 0.01 | 78% | 82% | 75% |
| 0.1 | 42% | 48% | 38% |
| 0.5 | 22% | 28% | 18% |
| 1.0 | 18% | 22% | 15% |

### Validation Protocol

1. **Clone** sRNA constructs into pSRNA backbone
2. **Transform** into E. coli BL21(DE3)
3. **Induce** with IPTG gradient (0-1 mM)
4. **Measure:**
   - Target mRNA levels (RT-qPCR)
   - Protein levels (Western blot or GFP fusion)
   - Lycopene titer (HPLC)
5. **Optimize** induction level for maximum lycopene
```

### LLM Agent Integration

```python
@tool
def design_srna(
    target_genes: list[str],
    organism: str = "e_coli",
    scaffold: str = "MicC",
    knockdown_target: float = 0.7,
    expression_system: str = "inducible"
) -> str:
    """
    Designs synthetic sRNAs for gene knockdown.

    Args:
        target_genes: List of gene symbols or locus tags
        organism: Target organism
        scaffold: sRNA scaffold (MicC, MicF, SgrS)
        knockdown_target: Desired knockdown fraction
        expression_system: constitutive or inducible

    Returns:
        sRNA sequences with predicted knockdown and plasmid design
    """
    pass
```

---

## Prerequisites

### Required Tools

| Resource | Purpose | Access |
|----------|---------|--------|
| **IntaRNA** | RNA-RNA interaction | Web API |
| **Mfold** | RNA structure prediction | Web service |
| **RBS Calculator** | Translation prediction | Salis Lab |

### Dependencies

```
ViennaRNA>=2.5
biopython>=1.81
numpy>=1.24
requests>=2.28
```

---

## References

- **Na et al. (2013):** "Metabolic engineering of E. coli using synthetic small regulatory RNAs." *Nature Biotechnology*
- **Yoo et al. (2020):** "Design of synthetic sRNA-based gene knockdown circuits." *ACS Synthetic Biology*

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

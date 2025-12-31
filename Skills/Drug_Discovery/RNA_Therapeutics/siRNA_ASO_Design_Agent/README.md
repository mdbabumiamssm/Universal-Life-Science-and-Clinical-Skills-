# siRNA and ASO Design Agent

**ID:** `biomedical.drug_discovery.sirna_aso_design`
**Version:** 1.0.0
**Status:** Production
**Category:** Drug Discovery / RNA Therapeutics

---

## Overview

The **siRNA and ASO Design Agent** automates the design of small interfering RNAs (siRNAs) and antisense oligonucleotides (ASOs) for gene silencing applications. These modalities represent a major therapeutic class with multiple FDA-approved drugs including Onpattro (patisiran), Leqvio (inclisiran), and Spinraza (nusinersen).

This agent integrates validated design rules, off-target prediction algorithms, and AI-based efficacy scoring to generate optimized oligonucleotide sequences for research and therapeutic development.

---

## Key Capabilities

### 1. siRNA Design

| Feature | Description |
|---------|-------------|
| **Target site selection** | Identify accessible mRNA regions |
| **Strand selection** | Optimize guide vs passenger strand asymmetry |
| **Efficacy prediction** | ML-based knockdown efficiency scoring |
| **Off-target analysis** | Seed region homology search |
| **Chemical modification** | 2'-OMe, 2'-F, phosphorothioate patterns |

### 2. ASO Design

| Mechanism | Target | Application |
|-----------|--------|-------------|
| **RNase H-dependent** | mRNA degradation | Gene knockdown |
| **Steric blocking** | Splice modulation | Exon skipping/inclusion |
| **miRNA inhibition** | AntimiRs | miRNA targeting |

### 3. Delivery Optimization

- **GalNAc conjugation:** Hepatocyte targeting
- **LNP formulation:** Systemic delivery
- **Naked delivery:** CNS applications (intrathecal)
- **Conjugate chemistry:** Targeting ligand selection

### 4. Safety Assessment

- **Off-target transcripts:** Genome-wide screening
- **Immune stimulation:** TLR activation prediction
- **Hepatotoxicity risk:** Chemical modification patterns
- **Sequence motifs:** Avoid toxic sequences (e.g., CpG)

---

## Usage

### Example Prompt

```text
Design 3 siRNAs targeting human PCSK9 (NM_174936) for liver-directed gene silencing.
Requirements:
- High knockdown efficiency (>80% predicted)
- Minimal off-targets in human transcriptome
- Compatible with GalNAc conjugation
- Include recommended chemical modifications

Also suggest an alternative ASO approach for the same target.
```

### Expected Output

```
## siRNA Design Report: PCSK9 (NM_174936)

### Target Analysis
- **Gene:** PCSK9 (Proprotein convertase subtilisin/kexin type 9)
- **Transcript:** NM_174936.4 (2,869 nt)
- **CDS:** 382-2,247 (2,053 nt)
- **Expression:** Liver-specific (ideal for GalNAc delivery)

### Top siRNA Candidates

| Rank | Position | Sense Strand (5'→3') | Efficacy | Off-targets |
|------|----------|---------------------|----------|-------------|
| 1 | 1,244 | GCUGGAGAACAUCUACCUGAA | 94% | 0 seed matches |
| 2 | 1,892 | CCAAGAGGCUUGACUUCAACU | 89% | 2 seed matches |
| 3 | 684 | GACGAGGACGGCAAGAUCAAU | 86% | 1 seed match |

### Lead Candidate Analysis (siPCSK9-1)

#### Sequences
```
Sense:      5'- GmCmUGGAGAACAUCUACCUGmAmA -3'
Antisense:  5'- UUCAGGUAGAUGUUCUCCAGmCmUmU -3'
            (m = 2'-O-methyl, all PS backbone)
```

#### Design Features
- **Position:** CDS nucleotide 1,244-1,264
- **GC content:** 47.6% (optimal range)
- **Thermodynamic asymmetry:** ΔΔG = -2.8 kcal/mol (favors antisense loading)
- **Seed region (nt 2-8):** CUGGAGA - no perfect matches to other transcripts

#### Off-Target Analysis
| Transcript | Seed Match | Full Alignment | Risk |
|------------|------------|----------------|------|
| None | - | - | Low |

#### Recommended Modifications
- **Sense strand:** 2'-OMe at positions 1, 2, 20, 21
- **Antisense strand:** 2'-OMe at positions 1, 2 and 2'-F at positions 3-19
- **Backbone:** Full phosphorothioate for stability
- **Conjugate:** Tri-GalNAc at sense strand 3' end

---

### Alternative ASO Design

#### Gapmer ASO (RNase H mechanism)
```
5'- mCmAmGmAmG-dCdAdGdTdGdAdCdTdGdA-mAmCmAmCmA -3'
    [Wing]    [DNA gap]              [Wing]
```

- **Length:** 20-mer
- **Chemistry:** 2'-MOE wings (5 nt each) + DNA gap (10 nt)
- **Target:** PCSK9 mRNA position 1,240-1,259
- **Predicted efficacy:** 78% knockdown

---

### Delivery Recommendations

| Conjugate | Tissue | Dosing | Reference Drug |
|-----------|--------|--------|----------------|
| **Tri-GalNAc** | Liver | 25-50 mg SC monthly | Inclisiran (Leqvio) |
| **LNP** | Liver | 0.3 mg/kg IV | Patisiran (Onpattro) |
```

### LLM Agent Integration

```python
@tool
def design_sirna(
    target_gene: str,
    num_candidates: int = 3,
    delivery_method: str = "GalNAc",
    check_off_targets: bool = True,
    species: str = "human"
) -> str:
    """
    Designs siRNA sequences for gene silencing.

    Args:
        target_gene: Gene symbol or RefSeq ID
        num_candidates: Number of siRNA candidates
        delivery_method: GalNAc, LNP, or naked
        check_off_targets: Perform off-target analysis
        species: Target species

    Returns:
        Ranked siRNA designs with modifications
    """
    pass


@tool
def design_aso(
    target_gene: str,
    mechanism: str = "gapmer",
    target_region: str = "cds",
    splice_site: str = None
) -> str:
    """
    Designs antisense oligonucleotides.

    Args:
        target_gene: Gene symbol or RefSeq ID
        mechanism: gapmer (RNase H) or steric_blocker
        target_region: cds, 5utr, 3utr, or splice_site
        splice_site: Specific exon for splice modulation

    Returns:
        ASO sequence with chemistry recommendations
    """
    pass
```

---

## Prerequisites

### Required Databases/Tools

| Resource | Purpose | Access |
|----------|---------|--------|
| **RefSeq** | Target sequences | NCBI API |
| **BLAST** | Off-target search | NCBI API or local |
| **siDirect** | siRNA design rules | Web service |
| **Sfold** | Target accessibility | Web service |
| **RNAfold** | Structure prediction | ViennaRNA |

### Dependencies

```
biopython>=1.81
ViennaRNA>=2.5
pandas>=1.5
numpy>=1.24
requests>=2.28
```

---

## Methodology

### siRNA Design Rules

The agent implements validated design algorithms:

| Rule | Description | Weight |
|------|-------------|--------|
| **Reynolds** | Position-specific nucleotide preferences | 0.25 |
| **Ui-Tei** | Thermodynamic asymmetry | 0.25 |
| **Amarzguioui** | A/U bias at 5' antisense | 0.20 |
| **GC content** | 30-52% optimal | 0.15 |
| **Secondary structure** | Accessible target site | 0.15 |

### Efficacy Scoring

```python
def calculate_sirna_efficacy(
    sense_seq: str,
    target_mrna: str,
    position: int
) -> float:
    """
    Calculate predicted siRNA efficacy using ensemble model.
    """
    features = {
        # Position-specific features
        'pos_1_A': sense_seq[0] == 'A',
        'pos_19_A': sense_seq[18] in ['A', 'U'],
        'pos_10_U': sense_seq[9] == 'U',

        # Thermodynamic features
        'delta_delta_G': calculate_asymmetry(sense_seq),
        'tm_sense': calculate_tm(sense_seq),

        # Compositional features
        'gc_content': calculate_gc(sense_seq),
        'a_content_3prime': count_au(sense_seq[14:19]),

        # Target accessibility
        'target_accessibility': calculate_accessibility(target_mrna, position)
    }

    # Ensemble prediction (RF + GBM + NN)
    score = ensemble_model.predict(features)
    return score
```

### Off-Target Analysis

```
Query: siRNA seed region (nt 2-8 of antisense)
    ↓
BLAST against RefSeq transcriptome
    ↓
Filter: Perfect seed match
    ↓
Extend alignment to full siRNA
    ↓
Score potential off-targets
    ↓
Report transcripts with >14/19 matches
```

---

## Related Skills

- **mRNA Design Agent:** Complementary RNA therapeutic
- **CRISPR Design Agent:** Alternative gene silencing approach
- **AgentD Drug Discovery:** Small molecule alternatives

---

## References

- **Setten et al. (2019):** "The current state and future directions of RNAi-based therapeutics." *Nature Reviews Drug Discovery*
- **Khvorova & Watts (2017):** "The chemical evolution of oligonucleotide therapies of clinical utility." *Nature Biotechnology*
- [TREAT Platform](https://zhanglabtools.github.io/TREAT/)
- [Alnylam RNAi Technology](https://www.alnylam.com/rnai-therapeutics)

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

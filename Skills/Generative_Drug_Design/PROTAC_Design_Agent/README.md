# PROTAC Design Agent

**ID:** `biomedical.generative_drug_design.protac`
**Version:** 1.0.0
**Status:** Production
**Category:** Generative Drug Design / Targeted Protein Degradation
**Released:** December 2025

---

## Overview

The **PROTAC Design Agent** leverages AI/ML for designing Proteolysis-Targeting Chimeras (PROTACs) and molecular glues for targeted protein degradation (TPD). TPD represents a paradigm shift in drug discovery, enabling degradation of "undruggable" proteins by hijacking the ubiquitin-proteasome system.

### Why TPD Matters

| Challenge | Traditional Drugs | PROTACs/Molecular Glues |
|-----------|-------------------|------------------------|
| Undruggable targets | Cannot address | Catalytic degradation |
| Resistance | Binding inhibition fails | Removes protein entirely |
| Selectivity | Off-target binding | E3 ligase specificity |
| Dosing | Stoichiometric | Catalytic (sub-stoichiometric) |

---

## Key Capabilities

### 1. Design Modes

| Mode | Input | Output | Use Case |
|------|-------|--------|----------|
| **De novo PROTAC** | Target + E3 ligase | Complete PROTAC | New targets |
| **Linker optimization** | Warhead + E3 ligand | Optimal linkers | Lead optimization |
| **Ternary prediction** | PROTAC + proteins | Complex structure | Mechanistic insight |
| **ADMET prediction** | PROTAC structure | PK properties | Drugability |

### 2. E3 Ligase Options

| E3 Ligase | Ligand | Tissue Expression | Applications |
|-----------|--------|-------------------|--------------|
| **CRBN** | Thalidomide derivatives | Ubiquitous | Most common |
| **VHL** | VH032 derivatives | Ubiquitous | Broad utility |
| **IAP** | Bestatin derivatives | Ubiquitous | Emerging |
| **MDM2** | Nutlin derivatives | p53 pathway | Oncology |
| **DCAF16** | Novel | CNS | Brain-penetrant |

### 3. AI/ML Components

```yaml
Ternary Complex Prediction:
  - AlphaFold Multimer for structure
  - MaSIF for interface prediction
  - RFdiffusion for linker design

Property Prediction:
  - Degradation efficiency (Dmax, DC50)
  - Hook effect prediction
  - Cell permeability
  - Metabolic stability

Generative Models:
  - VAE for linker generation
  - Diffusion for 3D conformations
  - Reinforcement learning for optimization
```

---

## Usage

### Example Prompt

```text
Design a PROTAC to degrade BRD4 (bromodomain protein) for cancer therapy.

Requirements:
- E3 ligase: CRBN (cereblon)
- Target warhead: JQ1 derivative
- Linker: Optimize for cell permeability
- Selectivity: Prefer BRD4 over BRD2/3
- Properties: Oral bioavailability desired

Generate 5 diverse candidates with predicted properties.
```

### Expected Output

```markdown
## PROTAC Design Campaign: BRD4 Degrader

### Target Analysis
| Property | Value |
|----------|-------|
| Target | BRD4 (Bromodomain-containing protein 4) |
| Function | Transcriptional regulator |
| Disease | Cancer (MYC-driven) |
| Existing drugs | JQ1, OTX015 (inhibitors only) |

### Design Strategy
- **Warhead**: JQ1 analog for BRD4 BD1/BD2 binding
- **E3 Ligand**: Pomalidomide (CRBN binder)
- **Linker**: PEG-alkyl hybrid (4-8 atoms)

### Top 5 PROTAC Candidates

| Rank | ID | Linker | Length | Pred DC50 | Dmax | MW | cLogP |
|------|-----|--------|--------|-----------|------|-----|-------|
| 1 | BRD4-P-001 | PEG2 | 6 atoms | 2.1 nM | 98% | 876 | 3.2 |
| 2 | BRD4-P-002 | Alkyl | 5 atoms | 4.3 nM | 95% | 848 | 4.1 |
| 3 | BRD4-P-003 | PEG3 | 9 atoms | 3.8 nM | 92% | 920 | 2.8 |
| 4 | BRD4-P-004 | Piperazine | 7 atoms | 8.2 nM | 88% | 892 | 3.5 |
| 5 | BRD4-P-005 | Triazole | 8 atoms | 12.4 nM | 85% | 904 | 3.9 |

### Lead Candidate: BRD4-P-001

#### Structure
```
        JQ1 Warhead          PEG2 Linker         Pomalidomide
    ┌──────────────┐    ┌──────────────┐    ┌──────────────┐
    │    BRD4      │────│   O-CH2-CH2  │────│    CRBN      │
    │   Binding    │    │   O-CH2-CH2  │    │   Binding    │
    └──────────────┘    └──────────────┘    └──────────────┘
```

#### Predicted Properties

| Property | Value | Target | Status |
|----------|-------|--------|--------|
| DC50 (BRD4) | 2.1 nM | <10 nM | Pass |
| Dmax | 98% | >90% | Pass |
| Selectivity (vs BRD2) | 15x | >10x | Pass |
| MW | 876 Da | <1000 | Pass |
| cLogP | 3.2 | <5 | Pass |
| TPSA | 145 Å² | <200 | Pass |
| Cell permeability | Moderate | - | Acceptable |
| Oral bioavailability | 25% (pred) | >20% | Pass |

#### Ternary Complex Analysis
- **Predicted binding mode**: Stable ternary complex
- **Interface area**: 1,250 Å²
- **Key contacts**: BRD4-linker-CRBN cooperative
- **Hook effect risk**: Low (DC50 << Dmax concentration)

### Synthesis Route
```
Synthetic accessibility score: 4.2/10 (Moderate)

Key steps:
1. JQ1 analog synthesis (known chemistry)
2. Linker attachment via amide coupling
3. Pomalidomide conjugation
4. Purification by HPLC

Estimated synthesis time: 2-3 weeks
```

### Recommended Next Steps

1. **Synthesize top 3 candidates**
2. **In vitro degradation assay** (Western blot)
3. **Selectivity profiling** (proteomics)
4. **Cell viability** in MYC-driven cancer lines
5. **PK studies** if cellular activity confirmed
```

---

## LLM Agent Integration

### Python Tool

```python
from typing import Optional, Dict, Any, List

def protac_design_tool(
    target_protein: str,
    e3_ligase: str = "CRBN",
    warhead_smiles: Optional[str] = None,
    n_designs: int = 5,
    optimize_for: List[str] = ["potency", "selectivity", "permeability"]
) -> Dict[str, Any]:
    """
    Design PROTACs for targeted protein degradation.

    Args:
        target_protein: Target protein name or UniProt ID
        e3_ligase: E3 ligase to use (CRBN, VHL, IAP)
        warhead_smiles: Optional SMILES for target-binding warhead
        n_designs: Number of PROTAC designs to generate
        optimize_for: Properties to optimize

    Returns:
        PROTAC designs with predicted properties
    """
    from protac_designer import PROTACGenerator, TernaryPredictor

    generator = PROTACGenerator()

    # Get or design warhead
    if warhead_smiles is None:
        warhead = generator.design_warhead(target_protein)
    else:
        warhead = warhead_smiles

    # Get E3 ligand
    e3_ligand = generator.get_e3_ligand(e3_ligase)

    # Generate PROTAC candidates
    designs = generator.generate_protacs(
        warhead=warhead,
        e3_ligand=e3_ligand,
        n_designs=n_designs,
        optimize_for=optimize_for
    )

    # Predict ternary complexes
    predictor = TernaryPredictor()
    for design in designs:
        design["ternary"] = predictor.predict(
            target_protein, design["smiles"], e3_ligase
        )

    return {
        "target": target_protein,
        "e3_ligase": e3_ligase,
        "designs": designs
    }
```

---

## Prerequisites

| Package | Version | Purpose |
|---------|---------|---------|
| RDKit | >=2023 | Chemistry |
| AlphaFold | >=2.3 | Structure prediction |
| PyTorch | >=2.0 | Deep learning |

---

## References

1. **Machine Learning in TPD (2025)**: "Machine learning in targeted protein degradation drug design." *Drug Discovery Today*.

2. **PROTAC-DB**: [http://cadd.zju.edu.cn/protacdb/](http://cadd.zju.edu.cn/protacdb/)

3. **Molecular Glues Review**: "Molecular Glues: The Adhesive Connecting TPD to the Clinic." *Biochemistry*.

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

---

*Last updated: December 2025*

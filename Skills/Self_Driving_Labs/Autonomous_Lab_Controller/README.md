# Autonomous Lab Controller

**ID:** `biomedical.self_driving_labs.autonomous_controller`
**Version:** 1.0.0
**Status:** Production
**Category:** Self-Driving Labs / Autonomous Research
**Released:** December 2025

---

## Overview

The **Autonomous Lab Controller** provides AI-driven orchestration for self-driving laboratories (SDLs), integrating robotics, machine learning, and automated experimental design to accelerate scientific discovery. SDLs can autonomously execute the entire scientific method - from hypothesis generation to data analysis - with minimal human intervention.

### The SDL Revolution

| Metric | Traditional Lab | Self-Driving Lab | Improvement |
|--------|-----------------|------------------|-------------|
| Experiments/day | 10-50 | 500-2000 | 40x |
| Data collection | Months | Days | 30x faster |
| Cost per experiment | $50-500 | $5-50 | 10x cheaper |
| Reproducibility | Variable | >99% | Standardized |
| 24/7 operation | No | Yes | 3x capacity |

### Autonomy Levels (2025 Framework)

| Level | Description | Current Status |
|-------|-------------|----------------|
| **Level 1** | Single automated steps | Mature |
| **Level 2** | Automated workflows | Mature |
| **Level 3** | Adaptive experiments | State-of-art |
| **Level 4** | Autonomous discovery | Emerging |
| **Level 5** | Self-directed research | Future |

Most current SDLs operate at **Level 3** - conditionally autonomous systems that can perform multiple cycles of the scientific method with human intervention only for anomalies.

---

## Key Capabilities

### 1. Scientific Workflow Automation

| Phase | Automation | AI Component |
|-------|------------|--------------|
| **Hypothesis** | Automated generation | LLM + literature mining |
| **Design** | Experiment planning | Bayesian optimization |
| **Execution** | Robotic operation | Motion planning |
| **Analysis** | Data processing | ML pipelines |
| **Learning** | Model updating | Active learning |

### 2. Supported Laboratory Types

| Domain | Applications | Key Instruments |
|--------|--------------|-----------------|
| **Chemistry** | Synthesis, catalysis | Liquid handlers, reactors |
| **Materials** | Formulation, characterization | Dispensers, analyzers |
| **Biology** | Cell culture, assays | Plate readers, imagers |
| **Drug Discovery** | Screening, ADMET | HTS platforms |
| **Genomics** | Library prep, sequencing | NGS platforms |

### 3. Integration Ecosystem

```yaml
Hardware:
  Liquid Handling:
    - Hamilton STAR
    - Tecan Freedom
    - Opentrons OT-2/Flex
  Robotics:
    - Universal Robots (UR3/UR5)
    - Franka Emika
    - ABB YuMi
  Instruments:
    - Plate readers (Molecular Devices)
    - Mass spectrometers
    - Microscopes
    - Chromatography

Software:
  Workflow:
    - Autoprotocol
    - SiLA 2
    - Synthace
  AI/ML:
    - Dragonfly (Bayesian opt)
    - Ax (Meta)
    - GPyOpt
  Data:
    - Electronic Lab Notebooks
    - LIMS integration
    - Cloud storage
```

---

## Usage

### Example Prompt

```text
Design and execute an autonomous optimization campaign for a Suzuki
coupling reaction.

Objective: Maximize yield while minimizing catalyst loading

Variables to optimize:
- Temperature: 60-120°C
- Time: 1-24 hours
- Catalyst loading: 0.5-5 mol%
- Base equivalents: 1-3 eq
- Solvent: THF, DMF, Dioxane

Constraints:
- Maximum 50 experiments
- Budget: $5,000
- Target yield: >90%

Execute autonomously with Bayesian optimization.
```

### Expected Output

```markdown
## Autonomous Lab Campaign: Suzuki Coupling Optimization

### Campaign Configuration

| Parameter | Setting |
|-----------|---------|
| Optimization method | Bayesian (Expected Improvement) |
| Surrogate model | Gaussian Process |
| Initial design | Latin Hypercube (10 pts) |
| Acquisition function | EI with exploration bonus |
| Total budget | 50 experiments |
| Parallelization | 8 reactions/batch |

### Search Space

| Variable | Type | Range | Prior |
|----------|------|-------|-------|
| Temperature | Continuous | 60-120°C | Uniform |
| Time | Continuous | 1-24 h | Log-uniform |
| Catalyst mol% | Continuous | 0.5-5% | Log-uniform |
| Base equiv | Continuous | 1-3 | Uniform |
| Solvent | Categorical | THF/DMF/Dioxane | Uniform |

### Execution Timeline

```
Batch 1 (Initial Design): 10 experiments
├── Status: COMPLETED
├── Time: 6.2 hours
├── Best yield: 72%
└── Key insight: DMF superior solvent

Batch 2 (Exploitation): 8 experiments
├── Status: COMPLETED
├── Time: 5.8 hours
├── Best yield: 84%
└── Key insight: Lower catalyst works at higher T

Batch 3 (Exploration): 8 experiments
├── Status: COMPLETED
├── Time: 7.1 hours
├── Best yield: 89%
└── Key insight: Time optimum at 8-10h

Batch 4 (Fine-tuning): 8 experiments
├── Status: COMPLETED
├── Time: 5.4 hours
├── Best yield: 94%
└── Key insight: Optimal region identified

Batch 5 (Validation): 8 experiments
├── Status: COMPLETED
├── Time: 5.2 hours
├── Best yield: 93% (avg: 92.3%)
└── Reproducibility: ±1.8%
```

### Optimization Trajectory

| Experiment | T (°C) | Time (h) | Cat% | Base eq | Solvent | Yield |
|------------|--------|----------|------|---------|---------|-------|
| 1 | 80 | 12 | 2.5 | 2.0 | THF | 45% |
| 2 | 100 | 6 | 1.0 | 1.5 | DMF | 68% |
| ... | ... | ... | ... | ... | ... | ... |
| 38 | **95** | **8** | **1.2** | **2.0** | **DMF** | **94%** |
| ... | ... | ... | ... | ... | ... | ... |

### Optimal Conditions Identified

| Parameter | Optimal Value | Confidence |
|-----------|---------------|------------|
| **Temperature** | 95°C | ±3°C |
| **Time** | 8 hours | ±1 h |
| **Catalyst loading** | 1.2 mol% | ±0.2% |
| **Base equivalents** | 2.0 eq | ±0.2 eq |
| **Solvent** | DMF | - |
| **Predicted yield** | 93.5% | ±2.1% |

### Cost Analysis

| Category | Actual | Budget | Savings |
|----------|--------|--------|---------|
| Reagents | $1,850 | $3,000 | $1,150 |
| Catalyst | $420 | $1,000 | $580 |
| Consumables | $380 | $500 | $120 |
| Instrument time | $650 | $500 | -$150 |
| **Total** | **$3,300** | **$5,000** | **$1,700** |

### Scientific Insights Discovered

1. **Catalyst loading-temperature tradeoff**: Lower catalyst (1-1.5%) sufficient at elevated T (90-100°C)

2. **Solvent effect**: DMF enables faster reaction; THF requires longer time for same yield

3. **Base equivalents**: 2.0 eq optimal; excess causes side products

4. **Reaction kinetics**: First-order in substrate; zero-order in catalyst above 1%

### Model Uncertainty

| Region | Prediction Uncertainty |
|--------|----------------------|
| Optimal region | ±2.1% |
| High T, low cat | ±4.5% |
| THF solvent | ±6.2% |
| Edge of space | ±8.3% |

### Recommended Follow-up

1. **Scale-up validation** at 10x scale
2. **Substrate scope** with optimal conditions
3. **Mechanistic studies** on catalyst loading effect
4. **Green chemistry**: Evaluate bio-based solvents
```

---

## LLM Agent Integration

### Python Tool Implementation

```python
from typing import Optional, Dict, Any, List, Literal, Union
from dataclasses import dataclass
import numpy as np

@dataclass
class ExperimentConfig:
    """Configuration for a single experiment."""
    parameters: Dict[str, float]
    constraints: Optional[Dict[str, Any]] = None
    metadata: Optional[Dict[str, Any]] = None

def autonomous_lab_controller(
    objective: str,
    variables: Dict[str, Dict[str, Any]],
    budget: int,
    acquisition_function: str = "expected_improvement",
    initial_design: str = "latin_hypercube",
    initial_points: int = 10,
    batch_size: int = 8,
    constraints: Optional[List[Dict]] = None,
    hardware_config: Optional[Dict] = None,
    simulation_mode: bool = False
) -> Dict[str, Any]:
    """
    Autonomous laboratory controller for self-driving experiments.

    Args:
        objective: Optimization objective (e.g., "maximize yield")
        variables: Search space definition
        budget: Maximum number of experiments
        acquisition_function: BO acquisition function
        initial_design: Initial design method
        initial_points: Number of initial experiments
        batch_size: Experiments per batch
        constraints: Experimental constraints
        hardware_config: Lab hardware configuration
        simulation_mode: Run in simulation (no hardware)

    Returns:
        Campaign results and optimal conditions
    """
    from dragonfly import minimise_function, maximise_function
    from sdl_controller import LabOrchestrator, ExperimentRunner

    # Parse objective
    maximize = "max" in objective.lower()

    # Initialize lab orchestrator
    if simulation_mode:
        orchestrator = LabOrchestrator.simulation()
    else:
        orchestrator = LabOrchestrator.from_config(hardware_config)

    # Define search space
    from dragonfly.exd.domains import EuclideanDomain, IntegralDomain
    from dragonfly.exd.cp_domain_utils import load_config

    domain_config = []
    for var_name, var_spec in variables.items():
        if var_spec["type"] == "continuous":
            domain_config.append({
                "name": var_name,
                "type": "float",
                "min": var_spec["range"][0],
                "max": var_spec["range"][1]
            })
        elif var_spec["type"] == "categorical":
            domain_config.append({
                "name": var_name,
                "type": "discrete",
                "items": var_spec["options"]
            })

    config = load_config({"domain": domain_config})

    # Objective function (interfaces with lab hardware)
    def run_experiment(params):
        # Convert params to experiment config
        exp_config = ExperimentConfig(parameters=dict(zip(variables.keys(), params)))

        # Execute on hardware
        result = orchestrator.run_experiment(exp_config)

        # Return objective value
        return result["yield"]

    # Run optimization
    if maximize:
        opt_val, opt_pt, history = maximise_function(
            run_experiment,
            config.domain,
            budget - initial_points,
            capital_type='num_evals',
            opt_method='bo',
            acq=acquisition_function,
            config=config
        )
    else:
        opt_val, opt_pt, history = minimise_function(
            run_experiment,
            config.domain,
            budget - initial_points,
            capital_type='num_evals'
        )

    # Compile results
    results = {
        "optimal_value": float(opt_val),
        "optimal_parameters": dict(zip(variables.keys(), opt_pt)),
        "experiments_run": len(history.query_vals),
        "history": {
            "parameters": [dict(zip(variables.keys(), p)) for p in history.query_points],
            "values": list(history.query_vals)
        },
        "model_predictions": {
            "mean": history.curr_opt_val,
            "uncertainty": history.curr_opt_unc if hasattr(history, 'curr_opt_unc') else None
        }
    }

    return results


class LabOrchestrator:
    """Orchestrates SDL hardware and software."""

    def __init__(self, hardware_clients: Dict[str, Any]):
        self.hardware = hardware_clients
        self.experiment_log = []

    @classmethod
    def from_config(cls, config: Dict) -> 'LabOrchestrator':
        """Initialize from hardware configuration."""
        from sdl_drivers import HamiltonClient, TecanClient, PlateReaderClient

        clients = {}
        if "liquid_handler" in config:
            if config["liquid_handler"]["type"] == "hamilton":
                clients["liquid_handler"] = HamiltonClient(
                    ip=config["liquid_handler"]["ip"]
                )
            elif config["liquid_handler"]["type"] == "tecan":
                clients["liquid_handler"] = TecanClient(
                    port=config["liquid_handler"]["port"]
                )

        if "plate_reader" in config:
            clients["plate_reader"] = PlateReaderClient(
                model=config["plate_reader"]["model"]
            )

        return cls(clients)

    @classmethod
    def simulation(cls) -> 'LabOrchestrator':
        """Initialize in simulation mode."""
        from sdl_drivers import SimulatedHardware
        return cls({"simulation": SimulatedHardware()})

    def run_experiment(self, config: ExperimentConfig) -> Dict[str, Any]:
        """Execute a single experiment."""
        # 1. Prepare reagents
        self.hardware["liquid_handler"].prepare_reaction(config.parameters)

        # 2. Run reaction
        self.hardware["liquid_handler"].execute_protocol()

        # 3. Analyze
        result = self.hardware["plate_reader"].measure()

        # 4. Log
        self.experiment_log.append({
            "config": config,
            "result": result,
            "timestamp": datetime.now()
        })

        return result


# Claude tool schema
SDL_CONTROLLER_SCHEMA = {
    "name": "autonomous_lab_experiment",
    "description": "Run autonomous laboratory experiments using self-driving lab infrastructure. Supports Bayesian optimization for reaction optimization, formulation screening, and biological assays.",
    "input_schema": {
        "type": "object",
        "properties": {
            "objective": {
                "type": "string",
                "description": "Optimization objective (e.g., 'maximize yield')"
            },
            "variables": {
                "type": "object",
                "description": "Search space variables with types and ranges"
            },
            "budget": {
                "type": "integer",
                "description": "Maximum number of experiments"
            },
            "batch_size": {
                "type": "integer",
                "description": "Experiments per parallel batch"
            },
            "simulation_mode": {
                "type": "boolean",
                "description": "Run in simulation without hardware"
            }
        },
        "required": ["objective", "variables", "budget"]
    }
}
```

---

## Prerequisites

### Hardware Integration

| Platform | SDK | Protocol |
|----------|-----|----------|
| Hamilton | Venus API | SiLA 2 |
| Tecan | Freedom EVOware | REST |
| Opentrons | Python API | Direct |
| Universal Robots | URScript | TCP/IP |

### Software Dependencies

```bash
# Core packages
pip install dragonfly-opt  # Bayesian optimization
pip install ax-platform    # Adaptive experimentation
pip install gpyopt         # Gaussian process optimization

# Hardware drivers
pip install pyhamilton
pip install pytecan
pip install opentrons

# Lab informatics
pip install autoprotocol
pip install sila2
```

### Infrastructure Requirements

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| CPU | 8 cores | 16 cores |
| RAM | 16 GB | 32 GB |
| GPU | Optional | NVIDIA RTX |
| Network | Ethernet | 10GbE |
| Storage | 500 GB SSD | 2 TB NVMe |

---

## Methodology

### Autonomous Discovery Loop

```
┌────────────────────────────────────────────────────────────┐
│               Self-Driving Lab Control Loop                 │
├────────────────────────────────────────────────────────────┤
│                                                             │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐ │
│  │  Hypothesis  │───▶│   Experiment │───▶│   Analysis   │ │
│  │  Generator   │    │   Designer   │    │   Pipeline   │ │
│  │  (LLM/ML)    │    │   (BayesOpt) │    │   (AutoML)   │ │
│  └──────────────┘    └──────────────┘    └──────────────┘ │
│         ▲                   │                   │          │
│         │                   ▼                   │          │
│         │           ┌──────────────┐            │          │
│         │           │   Hardware   │            │          │
│         │           │   Executor   │            │          │
│         │           │   (Robotics) │            │          │
│         │           └──────────────┘            │          │
│         │                   │                   │          │
│         └───────────────────┴───────────────────┘          │
│                     Model Update                            │
│                                                             │
└────────────────────────────────────────────────────────────┘
```

---

## Industry Examples

| Institution | SDL System | Application | 10x Faster |
|-------------|------------|-------------|------------|
| Argonne | Polybot | Materials screening | 90K samples |
| MIT | Chemputer | Organic synthesis | Novel molecules |
| Liverpool | Robot Scientist | Photocatalysis | 8 catalysts/day |
| CMU | Cloud Lab | Protein engineering | Remote operation |

---

## Related Skills

- `biomedical.lab_automation.opentrons` - Opentrons control
- `biomedical.lab_automation.pylabbot` - PyLabRobot
- `biomedical.sdl.bayesian_designer` - Experiment design
- `biomedical.drug_discovery.agentd` - Drug discovery workflows

---

## References

1. **SDL Review (2025)**: "Autonomous 'self-driving' laboratories: a review of technology." *Royal Society Open Science*.

2. **Nature Chem Eng (2025)**: "AI-powered lab discovers materials 10x faster."

3. **Aspuru-Guzik Group**: [AI for Discovery and Self-Driving Labs](https://www.matter.toronto.edu/)

4. **Argonne Polybot**: [Autonomous Discovery](https://www.anl.gov/autonomous-discovery)

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu

---

*Last updated: December 2025*

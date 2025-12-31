# NVIDIA FLARE Agent

**ID:** `biomedical.research.federated_learning.nvflare`
**Version:** 1.0.0
**Status:** Experimental
**Category:** Research Tools / Privacy-Preserving AI

---

## Overview

The **NVIDIA FLARE Agent** enables the orchestration of Federated Learning (FL) workflows. It allows researchers to train AI models (e.g., for tumor segmentation) across multiple institutions (hospitals) without ever moving patient data from the secure local site.

## Key Capabilities

- **Federated Training:** Orchestrates the "Server-Client" training loop.
- **Privacy Preservation:** Ensures only model weights (gradients) are shared, not raw data.
- **Algorithm Support:** Compatible with PyTorch, TensorFlow, and MONAI.

## Usage

This agent helps configure the `config_fed_server.json` and `config_fed_client.json` files required to launch an FL job.

```python
# Example: Configuring a Federated Averaging (FedAvg) job
fl_agent.configure_job(
    algorithm="FedAvg",
    num_rounds=100,
    min_clients=3,
    model_arch="UNet"
)
```

## References
- *NVIDIA FLARE GitHub*
- *FLamby Benchmark*

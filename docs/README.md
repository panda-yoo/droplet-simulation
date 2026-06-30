# Memory-Driven Dynamics of Self-Propelled Droplets

<p align="center">
  <img src="images/0005890_infer.png" width="500" alt="YOLO droplet tracking example">
</p>

<p align="center">

**Physics-Informed Simulation • YOLO Tracking • Statistical Analysis • Active Matter**

</p>

<p align="center">

![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)
![Status](https://img.shields.io/badge/Status-Active%20Development-orange.svg)
![Research](https://img.shields.io/badge/Research-Active%20Matter-success.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)

</p>

---

## Overview

This repository provides a unified computational framework for studying **memory-driven self-propelled active droplets** from both **experimental microscopy videos** and **physics-based numerical simulations**.

The project integrates

- **YOLO-based droplet detection and tracking**
- **FEniCSx PDE–SDE simulations**
- **Unified statistical trajectory analysis**
- **Visualization and diagnostic tools**
- **Reproducible computational workflows**

Both experimental and simulated trajectories share the **same CSV format**, allowing every statistical analysis routine to be applied without modification.

---

# Overall Pipeline

```mermaid
flowchart LR

    subgraph EXP["Experimental Pipeline"]
        V[Microscopy Videos]
        Y[YOLO Detection & Tracking]
        C1[Trajectory CSV]
        V --> Y --> C1
    end

    subgraph SIM["Simulation Pipeline"]
        P[PDE Chemical Field<br/>FEniCSx]
        S[SDE Droplet Dynamics]
        C2[Trajectory CSV]
        P --> S --> C2
    end

    C1 --> A
    C2 --> A

    subgraph ANA["Unified Statistical Analysis"]
        A[Load Trajectories]
        B[Recompute Canonical Velocity]

        B --> MSD
        B --> VACF
        B --> PSD
        A --> ORI[Orientation Correlation]
        A --> PERS[Persistence Analysis]

        MSD --> ENS
        VACF --> ENS
        PSD --> ENS
        ORI --> ENS
        PERS --> ENS

        ENS[Ensemble Statistics]
    end

    ENS --> FIG[Plots • Statistics • Diagnostics]
```

---

# Repository Structure

```text
.
├── analysis/
│   ├── MSD
│   ├── VACF
│   ├── PSD
│   ├── Orientation
│   ├── Persistence
│   ├── Plotting
│   └── Utilities
│
├── simulation/
│   ├── PDE Solver
│   ├── SDE Integrator
│   └── CSV Export
│
├── tracking/
│   ├── YOLO Detection
│   ├── Object Tracking
│   └── CSV Export
│
├── docs/
│   ├── Architecture
│   ├── Theory
│   ├── API
│   └── Images
│
└── output/
```

---

# Quick Start

```python
from main import DropletPipeline

pipeline = DropletPipeline()

# Analyze existing trajectory CSVs
pipeline.run_analysis(
    base_dir="output_droplet_tracking/csv_files",
    output_dir="./output",
    dt=0.2
)

# Run experimental tracking
pipeline.run_tracking(
    output_dir="output_droplet_tracking"
)

# Run numerical simulation
pipeline.run_simulation()
```

---

# Analysis Workflow

```mermaid
flowchart TD

    CSV[Trajectory CSV]

    CSV --> LOAD[Load Trajectories]

    LOAD --> POS[(x,y Positions)]

    POS --> VEL[Canonical Velocity Reconstruction]

    VEL --> MSD
    VEL --> VACF
    VEL --> PSD

    POS --> ORI[Orientation Correlation]

    MSD --> FIT1[Power-law Analysis]
    VACF --> FIT2[Persistence Time]
    PSD --> FIT3[Spectral Analysis]
    ORI --> FIT4[Rotational Persistence]

    FIT1 --> ENS[Ensemble Averaging]
    FIT2 --> ENS
    FIT3 --> ENS
    FIT4 --> ENS

```


---

# Implemented Features

## Experimental Pipeline

- YOLO object detection
- Multi-object tracking
- Automatic trajectory extraction
- CSV generation

---

## Simulation Pipeline

- Coupled PDE–SDE model
- FEniCSx finite-element solver
- Euler–Maruyama stochastic integration
- Memory-driven active droplet dynamics
- Trajectory export

---

## Statistical Analysis

Current analyses include

- Mean Squared Displacement (MSD)
- Local MSD exponent
- Velocity Autocorrelation Function (VACF)
- Orientation Correlation Function
- Power Spectral Density (PSD)
- Translational persistence time
- Rotational persistence time
- Ensemble averaging
- Statistical visualization

---

# Current Development

The statistical analysis module is actively being expanded to investigate correlated and memory-driven dynamics using advanced multivariate techniques.

```mermaid
flowchart LR

subgraph Current
A[MSD]
B[VACF]
C[PSD]
D[Orientation]
E[Persistence]
F[Ensemble Statistics]
end

subgraph Under Development
G[Cross Correlation]
H[Lag-time Correlation Matrix]
I[Singular Value Decomposition]
J[Canonical Correlation Analysis]
K[Eigenmode Analysis]
L[Spectral Filtering]
M[Noise Characterization]
end

Current --> Under_Development
```

Current research focuses on

- Cross-correlation between droplet trajectories
- Lag-time correlation matrices
- Eigenvalue and eigenvector analysis
- Singular Value Decomposition (SVD)
- Canonical Correlation Analysis (CCA)
- Dominant correlation modes
- Improved ensemble statistics
- Spectral filtering
- Statistical inference for memory-driven dynamics

---

# Documentation

## 📄 Complete Report

The project report provides the theoretical background, numerical methods, implementation details, and statistical analyses used throughout the repository.

**Report**

https://github.com/panda-yoo/Memory-Driven-Dynamics-of-Self-Propelled-Droplets---Report-and-Presentation/blob/main/report/main.pdf

---

## 🎓 Presentation

Project presentation slides

https://github.com/panda-yoo/Memory-Driven-Dynamics-of-Self-Propelled-Droplets---Report-and-Presentation/blob/main/presentation/main.pdf

---

# Documentation Index

| Document | Description |
|-----------|-------------|
| Project Structure | Repository organization |
| Dependency Graph | Module dependencies |
| Data Flow | End-to-end workflow |
| MSD | Theory & implementation |
| VACF | Theory & implementation |
| Orientation | Theory & implementation |
| PSD | Theory & implementation |
| Persistence | Translational & rotational persistence |
| Diagnostics | Ensemble consistency |
| Tracking Pipeline | YOLO workflow |
| Simulation | PDE–SDE model |
| API | Function reference |
| Developer Guide | Extending the project |

---

# Trajectory CSV Format

Every trajectory follows

```text
droplet_id, x, y, vx, vy, theta
```

| Column | Description |
|---------|-------------|
| `droplet_id` | Droplet identifier |
| `x`, `y` | Position coordinates |
| `vx`, `vy` | Stored forward-difference velocity |
| `theta` | Heading angle |

> **Note:** During analysis, velocities are always recomputed directly from the position coordinates using the canonical forward-difference definition to ensure consistency between experimental and simulated datasets.

---

# Major Dependencies

| Package | Purpose |
|-----------|----------|
| NumPy | Numerical computing |
| SciPy | Signal processing & optimization |
| Matplotlib | Visualization |
| Statsmodels | Statistical regression |
| Ultralytics | YOLO object detection & tracking |
| FEniCSx (dolfinx) | PDE solver |
| PETSc | Linear algebra backend |
| mpi4py | Parallel computing |
| PyVista | Simulation visualization |

See `requirements.txt` for the complete software environment.

---

# Future Directions

The framework is designed to be extensible. Planned additions include

- Generalized Langevin Equation (GLE) models
- Memory kernel estimation
- Cross-correlation analysis
- Canonical Correlation Analysis (CCA)
- Eigenmode decomposition
- Singular Value Decomposition (SVD)
- Memory kernel estimation
- Additional statistical diagnostics
- Improved visualization tools

---
# References

For additional theoretical background, implementation details, and analysis methods, see the accompanying project report and presentation.

```text
Memory-Driven Dynamics of Self-Propelled Droplets

Author:
Pranav Deepak Shinde

Supervisor:
Dr. T. K. Shajahan

M.Sc. Physics
National Institute of Technology Karnataka (NITK Surathkal)

2026
```

---

**Author:** Pranav Deepak Shinde  
**Supervisor:** Dr. T. K. Shajahan  
**Project:** Memory-Driven Dynamics of Self-Propelled Droplets  
**Institution:** National Institute of Technology Karnataka (NITK Surathkal)
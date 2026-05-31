# Simulation Overview

**Module:** `simulation/main.py`, `simulation/utils.py`

---

## Model description

The simulation models a set of `Ndrop` active droplets moving in a 2-D domain (the unit square `[0,1]²`). Each droplet:

1. **Secretes a chemical trail** described by a diffusion PDE
2. **Responds to the chemical gradient** via a drift force (self-avoidance)
3. **Moves stochastically** via an SDE with Gaussian white noise

The model is based on the active-droplet paper framework where droplets shrink over time, modulating their secretion rate and friction coefficient.

---

## Chemical field (diffusion PDE)

The chemical concentration field `u(x, t)` satisfies:

$$\frac{\partial u}{\partial t} = D \nabla^2 u + \sum_{i=1}^{N} \alpha_i(t) \, G_i(\mathbf{x}, t)$$

where `G_i` is a Gaussian source centred on droplet `i`:

$$G_i(\mathbf{x}, t) = \frac{1}{2\pi R_i^2} \exp\!\left(-\frac{|\mathbf{x} - \mathbf{X}_i|^2}{2 R_i^2}\right)$$

**Parameters (from `simulation/main.py`):**

| Parameter | Symbol | Value | Meaning |
|-----------|--------|-------|---------|
| `D_val` | D | 1e-5 | Diffusion coefficient |
| `beta` | β | 400.0 | Gradient force strength |
| `sigma` | σ | 5e-2 | Noise amplitude |
| `T` | T | 2.0 | Simulation duration |
| `R0` | R₀ | 0.040 | Initial droplet radius |
| `Rshrink` | Ṙ | 2.5e-4 | Radius shrinkage rate |
| `Rmin` | R_min | 0.008 | Minimum radius |
| `delta` | δ | 0.060 | Chemical interaction length |
| `eta` | η | 1.0 | Fluid viscosity |
| `Ndrop` | N | 6 | Number of droplets |
| `dt_val` | Δt | 0.002 | Time step |

---

## Time-dependent parameters

From `simulation/utils.update_paper_parameters`:

```python
def update_paper_parameters(t, R0, Rshrink, Rmin, delta, eta):
    Rval     = max(R0 - Rshrink * t, Rmin)
    alphaval = 3.0 * (Rval**2) / (delta**3) * Rshrink   # secretion rate
    gammaval = 6.0 * np.pi * eta * Rval                  # Stokes friction
    return Rval, alphaval, gammaval
```

The droplet shrinks at rate `Rshrink` until it reaches `Rmin`. The secretion rate and friction are both functions of the current radius, following Stokes' law for a sphere.

---

## Droplet SDE (position update)

```python
drift_prefactor = beta * (Rval**2) / (delta**2) / gammaval
noise_scale = sqrt(sigma * dt_val) / gammaval

Xnew = X[i] - dt_val * drift_prefactor * Fx[i] + noise_scale * randn()
Ynew = Y[i] - dt_val * drift_prefactor * Fy[i] + noise_scale * randn()
```

where `Fx[i]`, `Fy[i]` are the gradient forces:

```python
Fx[i] = ∫ G_i(x) · ∂u/∂x  dx    # assembled via FEniCSx
Fy[i] = ∫ G_i(x) · ∂u/∂y  dx
```

The SDE is an Euler-Maruyama discretisation of:

$$\gamma_i \dot{\mathbf{X}}_i = -\beta \frac{R_i^2}{\delta^2} \nabla_{\mathbf{X}_i} \int G_i \cdot u \, d\mathbf{x} + \sqrt{\sigma} \boldsymbol{\xi}(t)$$

**Reflecting boundaries:** Droplets bounce off all four walls by position reflection.

---

## Finite element solve

The PDE is discretised using implicit Euler in time and P1 Lagrange elements on a `Nx × Ny = 400 × 400` triangular mesh:

```
a = (u·v + dt·D·∇u·∇v) dx
L = (u_n·v + dt·source·v) dx
```

Solved with FEniCSx `LinearProblem` using a direct LU factorisation (`ksp_type=preonly`, `pc_type=lu`).

---

## Outputs

```
simulation data/
    multi_droplet_active.gif          ← animated field + droplet positions
    final_field.png                   ← chemical field at end of simulation
    droplet_<N>_trajectories.png      ← x-y trajectory plot
    droplet_<N>_tracking.csv          ← trajectory CSV (same format as tracking)
    data_droplet_<N>_trajectories.npz ← raw trajectory arrays
```

The CSV has columns: `droplet_id, x, y, vx, vy, theta` and is compatible with `analysis.trajectories.load_positions`.

---

## Running the simulation

```python
from main import DropletPipeline
p = DropletPipeline()
p.run_simulation()
```

Or directly:

```python
from simulation.main import run_simulation
run_simulation()
```

**Requirements:** `mpi4py`, `petsc4py`, `dolfinx`, `pyvista`, `tqdm`

The simulation runs on MPI. For a single-process run, no special command is needed. For multi-process (e.g. 4 cores):

```bash
mpirun -n 4 python main.py
```

---

## Analysing simulation output

```python
from main import DropletPipeline
p = DropletPipeline()
p.run_analysis(
    base_dir="simulation data",
    output_dir="./output_sim",
    dt=0.002,   # must match dt_val in simulation/main.py
)
```

The `dt` here must match `dt_val = 0.002` used in the simulation, otherwise all time-scaled quantities (τ_p, τ_r, PSD frequencies) will be wrong.

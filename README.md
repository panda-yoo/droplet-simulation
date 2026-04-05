# Chapter 2 — Background / Theory

## Key Concepts

### Vicsek angle $\theta$

```math
\theta_{t+1} = \theta_t + \chi_A \nabla A - \chi_B \nabla B
```

```math
\theta \sim f(\nabla A, \nabla B)
```

---

### About Vicsek Model [[1]](#ref1)

- The model consists of $N$ particles placed on a two-dimensional square box of size $L$.

- Each particle is characterized by:
  - Position: $x_n(t)$
  - Velocity:

```math
v_n(t) = v_0 e^{i\theta_n(t)}
```

- All particles move with the same speed $v_0$.  
  Only the direction $\theta_n(t)$ differs.

- Define:
  - $U_n(r_0; t)$ → circular neighborhood of radius $r_0$
  - $k_n(t)$ → number of particles in the neighborhood

- Average velocity:

```math
V_n(t) = \frac{1}{k_n(t)} \sum_{j : x_j \in U_n} v_j(t)
```

- This is the mean velocity of neighboring particles.

- Noise:
  - $\xi_n(t) \sim \text{Uniform}[-\pi, \pi]$
  - $\eta$ is the noise intensity

---

### Noise Models

#### Intrinsic Noise

```math
\theta_n(t + \Delta t) = \mathrm{Angle}(V_n(t)) + \eta \xi_n(t)
```

```math
v_n(t + \Delta t) = v_0 e^{i\theta_n(t+\Delta t)}
```

```math
x_n(t + \Delta t) = x_n(t) + v_n(t + \Delta t)\Delta t
```

---

#### Extrinsic Noise

```math
\theta_n(t + \Delta t) = \mathrm{Angle}\left(V_n(t) + \eta e^{i\xi_n(t)}\right)
```

```math
v_n(t + \Delta t) = v_0 e^{i\theta_n(t+\Delta t)}
```

```math
x_n(t + \Delta t) = x_n(t) + v_n(t + \Delta t)\Delta t
```

---

## Chemotaxis

Chemotaxis is the movement of cells or organisms in response to chemical gradients.

---

## Mathematical / Technical Details

### Field representing the medium

```math
A_{t+1} = A_t + D_A \nabla^2 A_t - \lambda A_t
```

---

### Field representing the droplet

```math
B_{t+1} = B_t + D_B \nabla^2 B_t - \lambda B_t
```

---

## Project Implementation

### 1. Vicsek Alignment Term

```math
\mathbf{V}_n^{(\text{align})} =
(\cos\theta_n, \sin\theta_n) \quad \text{if } k_n = 0
```

```math
\mathbf{V}_n^{(\text{align})} =
\frac{1}{k_n} \sum_{j \neq n,\; d_{nj}^2 < R^2}
v_j (\cos\theta_j, \sin\theta_j) \quad \text{otherwise}
```

---

### 2. Chemotaxis Term

```math
\mathbf{V}_n^{(\text{chem})} =
\left(
\partial_x A - \partial_x B,\;
\partial_y A - \partial_y B
\right)\Big|_{(x_n, y_n)}
```

---

### 3. Combined Velocity

```math
\mathbf{V}_n =
\mathbf{V}_n^{(\text{align})}
+
\mathbf{V}_n^{(\text{chem})}
```

---

### 4. Position Update (Periodic Boundary Conditions)

```math
x_n^{(\text{unwrapped})}(t+\Delta t) =
x_n^{(\text{unwrapped})}(t) + v_n \cos\theta_n(t)
```

```math
y_n^{(\text{unwrapped})}(t+\Delta t) =
y_n^{(\text{unwrapped})}(t) + v_n \sin\theta_n(t)
```

```math
x_n(t+\Delta t) =
x_n^{(\text{unwrapped})}(t+\Delta t) \bmod N_x
```

```math
y_n(t+\Delta t) =
y_n^{(\text{unwrapped})}(t+\Delta t) \bmod N_y
```

---

### 5. Field Update (ChemA and ChemB)

- Compute spatial average:

```math
\langle c \rangle
```

- Compute gradient:

```math
\nabla \langle c \rangle
```

- Field evolution (diffusion):

```math
c_{t+1} = \nabla^2 c_t
```

---

## Notes

- Use ```math blocks for GitHub rendering  
- Avoid $$ ... $$  
- Avoid aligned environments  

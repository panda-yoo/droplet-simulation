---
title: droplet-sim
updated: 2026-04-05 16:30:09Z
created: 2026-03-28 20:30:06Z
latitude: 13.08577790
longitude: 74.79425620
altitude: 0.0000
---

# Project Title

**Author:** Pranav Shinde
**Date:** 29 March 2026 
**Version:** v0.1  

---

# Chapter 1 — Overview

## Goal
- To simulate the droplet intraction dynamics
  
## Problem Statement

-

## Motivation

-

## Expected Outcome

-

---

# Chapter 2 — Background / Theory

## Key Concepts
### - Vicsek angle $\theta$ 
-
$$
	\theta _{t+1} = \theta _{t} + \chi_A \nabla A - \chi_B \nabla B
	$$
- 
  	$$ \theta \sim f(\nabla A,\nabla B)$$

- about vicsek model [[1]](#ref1)
	- The model consists of $N$ particles placed on a two-dimensional square box of sides $L$ 
	- Each particle is characterized by its position $x_n(t)$, and its velocity $v_n(t) = v_0 e^{i\theta_n(t)}$
	- All the particles move with the same speed $v_0$. What changes from one particle to another is the direction of motion given by the angle $\theta_n(t)$.
	-  to state mathematically the interaction rules between the particles, we can define $U_n(r_0;t)$ as the circular neighborhood with radius $r_0$ centered at $x_n(t)$, and $k_n(t)$ the number of particles within this neighborhood. Let $V_n(t)$ be the average velocity of the particles within each neighborhood $U_n(r_0;t)$:
		
       - $V_n(t) = \frac{1}{k_n(t)} \sum_{j : x_j \in U_n} v_j(t).$
       - The above formula is just the sum of the velocities of each particle in the neighborhood divided by the number of particles in that neighborhood. 
  - The interaction between the particles with noise is given by the simultaneous updating of all the velocity angles and all the particle positions according to the rules
  - Where $\xi_n(t)$ is a random variable uniformly distributed in the interval $[-\pi, \pi]$, and $\eta$ is the noise intensity, which is a positive parameter.
		
$$
\begin{aligned}
&\text{Intrinsic noise:} \\[10pt]
\theta_n(t + \Delta t) &= \mathrm{Angle}\!\left[ V_n(t) \right] + \eta \, \xi_n(t),\\[10pt]
v_n(t + \Delta t) &= v_0 e^{i\theta_n(t+\Delta t)},\\[10pt]
x_n(t + \Delta t) &= x_n(t) + v_n(t + \Delta t)\Delta t \\[15pt]

&\text{Extrinsic noise:} \\[10pt]
\theta_n(t + \Delta t) &= \mathrm{Angle}\!\left[ V_n(t) + \eta e^{i\xi_n(t)} \right],\\[10pt]
v_n(t + \Delta t) &= v_0 e^{i\theta_n(t+\Delta t)},\\[10pt]
x_n(t + \Delta t) &= x_n(t) + v_n(t + \Delta t)\Delta t
\end{aligned}
$$

### - Chemotaxis:
 - Chemotaxis : Chemotaxis is the phenomenon whereby somatic cells, bacteria, and other single-cell or multicellular organisms direct their movements in response to certain chemicals in their environment.

## Mathematical / Technical Details

- field representing the medium
$$
  A_{t+1}  = A_{t} + D_A \nabla^2A_{t} - \lambda A_{t}
$$

- field representing the droplet
$$
  B_{t+1}  = B_{t} + D_B \nabla^2B_{t} - \lambda B_{t}
$$ 


### - project implmentation 
1. Vicsek alignment term
$$
\begin{aligned}
\mathbf{V}_n^{(\text{align})}
&=
\begin{cases}
\left(\cos\theta_n,\; \sin\theta_n\right), & \text{if } k_n = 0, \\[8pt]
\displaystyle \frac{1}{k_n}
\sum_{\substack{j \neq n \\ d_{nj}^2 < R^2}}
v_j \left(\cos\theta_j,\; \sin\theta_j\right), & \text{otherwise , here R is a particle property}\\
\end{cases}
\end{aligned}
$$


2. Chemotaxis term    
$$
\begin{aligned}
\mathbf{V}_n^{(\text{chem})}
=
\left(
\partial_x A - \partial_x B,\;
\partial_y A - \partial_y B
\right)\Big|_{(x_n,\,y_n)}
\end{aligned}
$$

3. Combined interpretation
$$
\mathbf{V}_n
=
\mathbf{V}_n^{(\text{align})}
+
\mathbf{V}_n^{(\text{chem})}
$$
4. Position update (with periodic boundary conditions)
$$
\begin{aligned}
x_n^{\text{(unwrapped)}}(t+\Delta t) &= x_n^{\text{(unwrapped)}}(t) + v_n \cos\theta_n(t), \\
y_n^{\text{(unwrapped)}}(t+\Delta t) &= y_n^{\text{(unwrapped)}}(t) + v_n \sin\theta_n(t), \\[8pt]

x_n(t+\Delta t) &= x_n^{\text{(unwrapped)}}(t+\Delta t) \;\bmod \; N_x, \\
y_n(t+\Delta t) &= y_n^{\text{(unwrapped)}}(t+\Delta t) \;\bmod \; N_y.
\end{aligned}
$$
5 . For field chemA and ChemB
- We Take Average of fields $<c>$
- And then calculate Gradient  $\nabla$ $<c>$
- for evolution of the field we consider $c_{t+1}$ = $\nabla^2c_{t}$


## References

1. Paper / Article 


    [1]  <a id="ref1"></a>	[MIT Vicsek Model of self-propelled particles](https://web.mit.edu/8.334/www/grades/projects/projects10/Hernandez-Lopez-Rogelio/dynamics_2.html)

   
2. Book  
3. Documentation  

---

# Chapter 3 — Resources

## Papers

-

## Tutorials

-

## Useful Links

- [Link name](https://example.com)

---

# Chapter 4 — Plan / Approach

## Strategy
* Data Types
  	- we will have 3 main fields/grid
  		- `chemA`
		- `chemB`
		- `is_occupied` 
	-   `Particle` class
	```python
	class Particle:
    def __init__(self, i, x, y, theta, radius) -> None:
        self.id = i
        # wrapped (for simulation)
        self.x = x
        self.y = y
        # unwrapped (for tracking)
        self.x_unwrapped = x
        self.y_unwrapped = y
        self.theta = theta
        self.radius = radius
        self.v = 2.0
        self.influence = 500.0

	```
## Steps

- [ ] Step 1  
- [ ] Step 2  
- [ ] Step 3  

---

# Chapter 5 — Implementation

## Setup

- Language:
- Tools:
- Libraries:
- Hardware:

---

## Folder Structure

```text
project/
│── src/
│── data/
│── results/
│── notebooks/
│── images/
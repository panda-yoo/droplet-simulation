"""
simulation.main
===============
Full 2D active-droplet simulation.

Diffusion PDE for the chemical trail + stochastic SDE for droplet motion.
Uses FEniCSx (DOLFINx) with PETSc for the PDE and PyVista for visualisation.

This version is tuned to look more like the paper figures:
- clustered initial injection
- shorter time window
- stronger short-time randomness
- trail-driven self-avoidance
- no artificial heading-angle propulsion

Call ``run_simulation()`` to execute.
"""

from mpi4py import MPI
from petsc4py.PETSc import ScalarType

import os
from pathlib import Path

from dolfinx import mesh, fem, plot
from dolfinx.fem.petsc import LinearProblem

import matplotlib as mpl
import matplotlib.pyplot as plt
from analysis.plot_style import set_publication_style, get_droplet_color, StyleTokens, apply_grid, apply_legend
import numpy as np
import pyvista
import ufl
from tqdm import tqdm

# ── local imports ────────────────────────────────────────
from simulation.utils import (
    initial_condition,
    update_paper_parameters,
    save_trajectory_csv,
)


# =========================================================
# Main simulation entry point
# =========================================================

def run_simulation():
    """
    Run the complete active-droplet simulation.

    Workflow
    --------
    1. Build the FE mesh and function space
    2. Place droplets at clustered initial positions
    3. Time-step: solve diffusion PDE → compute gradient forces → SDE update
    4. Record an animated GIF and final-field screenshot
    5. Save trajectory data (.npz and .csv)
    """

    # =========================================================
    # Output location
    # =========================================================
    output = 'simulation data'
    os.makedirs(output, exist_ok=True)

    # =========================================================
    # Simulation settings
    # =========================================================
    Ndrop = 6
    seed = 4
    np.random.seed(seed)

    Nx = 400
    Ny = 400

    dt_val = 0.002          # closer to the paper's time scale
    frame_every = 1

    # =========================================================
    # Physical parameters (dimensionless)
    # =========================================================

    D_val = 1e-5

    beta = 400.0

    sigma = 5e-2

    T = 2.0                # shorter lifetime => less spaghetti

    R0 = 0.040
    Rshrink = 2.5e-4
    Rmin = 0.008

    delta = 0.060
    eta = 1.0


    # =========================================================
    # Mesh and finite-element space
    # =========================================================
    domain = mesh.create_unit_square(MPI.COMM_WORLD, Nx, Ny)
    comm = domain.comm

    V = fem.functionspace(domain, ("Lagrange", 1))

    un = fem.Function(V)
    un.interpolate(initial_condition)

    uh = fem.Function(V)

    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    # =========================================================
    # Time constants used by the PDE solve
    # =========================================================
    dt = fem.Constant(domain, ScalarType(dt_val))
    D = fem.Constant(domain, ScalarType(D_val))

    # =========================================================
    # Initial droplet positions
    # Central injection, like the paper's multiple-droplet setup
    # =========================================================
    x0 = np.clip(0.50 + 0.03 * np.random.randn(Ndrop), 0.30, 0.70)
    y0 = np.clip(0.50 + 0.03 * np.random.randn(Ndrop), 0.30, 0.70)

    X = [fem.Constant(domain, ScalarType(x0[i])) for i in range(Ndrop)]
    Y = [fem.Constant(domain, ScalarType(y0[i])) for i in range(Ndrop)]
    R = [fem.Constant(domain, ScalarType(R0)) for _ in range(Ndrop)]
    alpha = [fem.Constant(domain, ScalarType(0.0)) for _ in range(Ndrop)]
    gamma_c = [fem.Constant(domain, ScalarType(1.0)) for _ in range(Ndrop)]

    # =========================================================
    # Coordinates and PDE source terms
    # =========================================================
    x_coord = ufl.SpatialCoordinate(domain)
    source = 0
    Fx_forms = []
    Fy_forms = []

    for i in range(Ndrop):
        r2 = (x_coord[0] - X[i])**2 + (x_coord[1] - Y[i])**2

        Gi = (1.0 / (2.0 * np.pi * R[i]**2)) * ufl.exp(
            -r2 / (2.0 * R[i]**2)
        )

        source += alpha[i] * Gi

        Fx_forms.append(fem.form(Gi * ufl.grad(uh)[0] * ufl.dx))
        Fy_forms.append(fem.form(Gi * ufl.grad(uh)[1] * ufl.dx))

    # =========================================================
    # Weak form and linear solver
    # =========================================================
    a = (u * v + dt * D * ufl.dot(ufl.grad(u), ufl.grad(v))) * ufl.dx
    L = (un * v + dt * source * v) * ufl.dx
    problem = LinearProblem(
        a,
        L,
        u=uh,
        petsc_options_prefix="diffusion_",
        petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
    )

    # =========================================================
    # PyVista visualisation setup
    # =========================================================
    topology, cells, geometry = plot.vtk_mesh(V)
    grid = pyvista.UnstructuredGrid(topology, cells, geometry)
    grid.point_data["u"] = un.x.array.copy()

    plotter = pyvista.Plotter(off_screen=True, window_size=(700, 700))
    viridis = mpl.colormaps["viridis"]

    actor = plotter.add_mesh(
        grid,
        cmap=viridis,
        show_edges=False,
        lighting=False,
        clim=[0, 0.02],
    )

    xc = [float(X[i].value) for i in range(Ndrop)]
    yc = [float(Y[i].value) for i in range(Ndrop)]
    walker_points = pyvista.PolyData(np.c_[xc, yc, np.zeros(Ndrop)])

    plotter.add_mesh(
        walker_points,
        color="white",
        point_size=12,
        render_points_as_spheres=True,
    )

    plotter.view_xy()
    plotter.add_text(
        "Frame: 0  t = 0.00",
        position="upper_left",
        font_size=10,
        color="white",
        name="frame_text",
    )
    plotter.open_gif(os.path.join(output, "multi_droplet_active.gif"))

    # =========================================================
    # Trajectory storage
    # =========================================================
    trajectory_x = [[] for _ in range(Ndrop)]
    trajectory_y = [[] for _ in range(Ndrop)]

    # =========================================================
    # Time stepping
    # =========================================================
    num_steps = int(T / dt_val)
    t = 0.0

    for step in tqdm(range(num_steps)):
        t += dt_val

        # Update paper-style shrinking parameters
        Rval, alphaval, gammaval = update_paper_parameters(
            t=t, R0=R0, Rshrink=Rshrink, Rmin=Rmin, delta=delta, eta=eta
        )

        for i in range(Ndrop):
            R[i].value = ScalarType(Rval)
            alpha[i].value = ScalarType(alphaval)
            gamma_c[i].value = ScalarType(gammaval)

        # Solve diffusion PDE
        problem.solve()
        un.x.array[:] = uh.x.array

        # Assemble gradient force on each droplet
        Fx = np.zeros(Ndrop)
        Fy = np.zeros(Ndrop)

        for i in range(Ndrop):
            Fx[i] = comm.allreduce(fem.assemble_scalar(Fx_forms[i]), op=MPI.SUM)
            Fy[i] = comm.allreduce(fem.assemble_scalar(Fy_forms[i]), op=MPI.SUM)

        # SDE update: trail-driven drift + noise
        drift_prefactor = beta * (Rval**2) / (delta**2) / gammaval
        noise_scale = np.sqrt(sigma * dt_val) / gammaval

        for i in range(Ndrop):
            Xnew = float(X[i].value) - dt_val * drift_prefactor * Fx[i] + noise_scale * np.random.randn()
            Ynew = float(Y[i].value) - dt_val * drift_prefactor * Fy[i] + noise_scale * np.random.randn()

            # Reflecting boundaries
            if Xnew < 0:
                Xnew = -Xnew
            if Xnew > 1:
                Xnew = 2 - Xnew

            if Ynew < 0:
                Ynew = -Ynew
            if Ynew > 1:
                Ynew = 2 - Ynew

            Xnew = min(max(Xnew, 0.0), 1.0)
            Ynew = min(max(Ynew, 0.0), 1.0)

            X[i].value = ScalarType(Xnew)
            Y[i].value = ScalarType(Ynew)

            trajectory_x[i].append(Xnew)
            trajectory_y[i].append(Ynew)

        # Render frame
        if step % frame_every == 0:
            scalars = uh.x.array.copy()
            grid.point_data["u"] = scalars
            actor.mapper.scalar_range = (0.0, max(float(np.max(scalars)), 1e-12))

            xc = [float(X[i].value) for i in range(Ndrop)]
            yc = [float(Y[i].value) for i in range(Ndrop)]
            walker_points.points = np.c_[xc, yc, np.zeros(Ndrop)]
            walker_points.Modified()

            plotter.add_text(
                f"Frame: {step // frame_every}  t = {t:.2f}",
                position="upper_left",
                font_size=10,
                color="white",
                name="frame_text",
            )
            plotter.write_frame()

    plotter.screenshot(os.path.join(output, "final_field.png"))
    plotter.close()

    # =========================================================
    # Plot trajectories
    # =========================================================
    set_publication_style()
    fig, ax = plt.subplots(figsize=(8, 8))

    for i in range(Ndrop):
        c = get_droplet_color(i + 1)
        ax.plot(
            trajectory_x[i],
            trajectory_y[i],
            linewidth=1.5,
            color=c,
            label=f"Droplet {i+1}"
        )
        # Add start and end markers
        ax.plot(trajectory_x[i][0], trajectory_y[i][0], "o", color=c, ms=6)
        ax.plot(trajectory_x[i][-1], trajectory_y[i][-1], "s", color=c, ms=6)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    apply_grid(ax, alpha=0.25)
    ax.set_title(f"Simulation : Trajectories\nN = {Ndrop} droplets")
    apply_legend(ax)
    fig.savefig(os.path.join(output, f"droplet_{Ndrop}_trajectories.png"), dpi=150)
    plt.close(fig)

    save_trajectory_csv(
        trajectory_x=trajectory_x,
        trajectory_y=trajectory_y,
        dt_val=dt_val,
        Ndrop=Ndrop,
        output_path=os.path.join(output, f"droplet_{Ndrop}_tracking.csv"),
    )
    np.savez(
        os.path.join(output, f"data_droplet_{Ndrop}_trajectories.npz"),
        trajectory_x=np.array(trajectory_x, dtype=object),
        trajectory_y=np.array(trajectory_y, dtype=object),
    )

    print(f"Saving data_droplet_{Ndrop}_trajectories.npz and \ndroplet_{Ndrop}_tracking.csv in {output} folder")

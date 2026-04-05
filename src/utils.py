import numpy as np
from numpy.typing import NDArray
from typing import Tuple, Optional, List
from src.particle import Particle


def periodic_distance_sq(x, y, xc, yc, Nx, Ny):
    """Calculate squared distance with periodic boundary conditions."""
    dx = x - xc
    dy = y - yc
    dx = dx - Nx * np.round(dx / Nx)
    dy = dy - Ny * np.round(dy / Ny)
    return dx * dx + dy * dy


def vick_theta(particle: Particle, particles: List[Particle],
               Nx: int, Ny: int) -> tuple[float, float]:
    dx = 0.0
    dy = 0.0
    particle_inside = 0

    for i in range(len(particles)):
        if particles[i].id == particle.id:
            continue
        d_sqr = periodic_distance_sq(x=particle.x, y=particle.y,
                                           xc=particles[i].x, yc=particles[i].y,
                                           Nx=Nx, Ny=Ny)
        if d_sqr < particle.influence * particle.influence:
            #     do something
            dx = dx + particles[i].v * np.cos(particles[i].theta)
            dy = dy + particles[i].v * np.sin(particles[i].theta)
            particle_inside += 1

    if particle_inside == 0:
        return np.cos(particle.theta), np.sin(particle.theta)

    return dx / particle_inside, dy / particle_inside


def chemotaxis_theta(particle: Particle,
                     gradAx: NDArray, gradAy: NDArray,
                     gradBx: NDArray, gradBy: NDArray) -> Tuple[float, float]:
    Ny, Nx = gradAx.shape

    i = int(particle.y) % Ny
    j = int(particle.x) % Nx
    fx = gradAx[i, j] - gradBx[i, j]
    fy = gradAy[i, j] - gradBy[i, j]

    return fx, fy


def update_theta(particle: Particle, particles: List[Particle],
                 gradAx: NDArray, gradAy: NDArray,
                 gradBx: NDArray, gradBy: NDArray,
                 alpha=1.0, beta=1.0, dt=0.1) -> float:
    vax, vay = chemotaxis_theta(particle=particle, gradAx=gradAx, gradAy=gradAy,
                                gradBx=gradBx, gradBy=gradBy)

    vbx, vby = vick_theta(particle=particle, particles=particles,
                          Nx=gradAx.shape[1], Ny=gradAx.shape[0])
    vx = alpha * vax + beta * vbx
    vy = alpha * vay + beta * vby
    if vx * vx + vy * vy < 1e-8:
        theta = 0.1 * np.random.randn()
    else:
        theta = np.arctan2(vy, vx)

    return theta + 0.1 * np.random.randn()


def populate_is_occupied(radius: float, xc: float, yc: float, id: int,
                         is_occupied: NDArray, collisions: List[tuple]):
    """Mark cells as occupied and track collisions."""
    Ny, Nx = is_occupied.shape
    xc_int, yc_int = int(xc), int(yc)

    for x in range(xc_int - int(radius), xc_int + int(radius) + 1):
        for y in range(yc_int - int(radius), yc_int + int(radius) + 1):
            dist2 = periodic_distance_sq(x, y, xc, yc, Nx, Ny)
            if dist2 <= radius * radius:
                xp = x % Nx
                yp = y % Ny
                if is_occupied[yp, xp] != -1:
                    # Collision detected
                    other_id = is_occupied[yp, xp]
                    if (other_id, id) not in collisions and (id, other_id) not in collisions:
                        collisions.append((id, other_id))
                is_occupied[yp, xp] = id


def move(particle: Particle, Nx: int, Ny: int):
    """Update particle position based on velocity
     and angle."""
    particle.x_unwrapped += particle.v * np.cos(particle.theta)
    particle.y_unwrapped += particle.v * np.sin(particle.theta)
    particle.x = particle.x_unwrapped % Nx
    particle.y = particle.y_unwrapped % Ny


def local_field_avg(chem: NDArray, R: float, dx: float = 1.0) -> NDArray:
    """Compute circular neighborhood averages."""
    from scipy.ndimage import convolve

    radius_cells = max(1, int(np.ceil(R / dx)))
    offsets = np.arange(-radius_cells, radius_cells + 1, dtype=float)
    yy, xx = np.meshgrid(offsets, offsets, indexing="ij")
    kernel = ((xx * dx) ** 2 + (yy * dx) ** 2) <= R * R
    kernel = kernel.astype(np.float32)
    kernel_sum = kernel.sum()
    if kernel_sum == 0:
        kernel[radius_cells, radius_cells] = 1.0
        kernel_sum = 1.0
    kernel = kernel / kernel_sum

    avg = convolve(chem, kernel, mode='wrap')
    return avg.astype(chem.dtype)


def gradient(field: NDArray) -> Tuple[NDArray, NDArray]:
    """Compute gradient of a field."""
    Ny, Nx = field.shape
    grad_x = np.full(field.shape, np.nan, dtype=field.dtype)
    grad_y = np.full(field.shape, np.nan, dtype=field.dtype)
    for y in range(Ny):
        for x in range(Nx):
            grad_x[y, x] = (field[y, (x + 1) % Nx] - field[y, (x - 1) % Nx]) * 0.5
            grad_y[y, x] = (field[(y + 1) % Ny, x] - field[(y - 1) % Ny, x]) * 0.5
    return grad_x, grad_y


def laplacian(field: NDArray) -> NDArray:
    """Compute laplacian of a field."""
    lap_field = np.full(shape=field.shape, fill_value=np.nan, dtype=field.dtype)
    Ny, Nx = field.shape
    for y in range(Ny):
        for x in range(Nx):
            term = (field[(y + 1) % Ny, x] + field[(y - 1) % Ny, x] +
                    field[y, (x + 1) % Nx] + field[y, (x - 1) % Nx] -
                    4 * field[y, x])
            lap_field[y, x] = term
    return lap_field


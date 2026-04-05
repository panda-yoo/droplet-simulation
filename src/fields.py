import numpy as np
from numpy.typing import NDArray
from typing import Tuple, Optional
from src.particle import Particle
from src import utils


def deposit_chemA(
        chemA: NDArray,
        is_occupied: NDArray,
        modify_fn=None,
        **kwargs
):
    """
    Modify chemA locally (typically inside particles) without overwriting the field.

    Parameters:
    - chemA: field
    - is_occupied: occupancy grid
    - modify_fn: function(old_values, **kwargs) -> new_values
    """
    mask = (is_occupied != -1)

    if modify_fn is None:
        # default: simple depletion
        chemA[mask] = 0.0
    else:
        # apply custom modification ONLY where mask is true
        chemA[mask] = modify_fn(chemA[mask], **kwargs)

def noise(values, scale=0.1):
    return values + scale * np.random.randn(*values.shape)

def weaken(values, factor=0.5):
    return values * factor

def deposit_chemB(particle: Particle, chemB: NDArray, strength: float = 1.0):
    """
    Add contribution from a particle locally (no overwrite).
    """

    Ny, Nx = chemB.shape

    xc, yc = particle.x, particle.y
    r = particle.radius

    xc_int, yc_int = int(xc), int(yc)

    for x in range(xc_int - int(r), xc_int + int(r) + 1):
        for y in range(yc_int - int(r), yc_int + int(r) + 1):

            dist2 = utils.periodic_distance_sq(x, y, xc, yc, Nx, Ny)

            if dist2 <= r * r:
                xp = x % Nx
                yp = y % Ny

                chemB[yp, xp] += strength


def decay(field, rate=0.01):
    return field * (1 - rate)


def init_gaussian(chemB: NDArray, x0, y0, sigma=10.0, amplitude=1.0):
    Ny, Nx = chemB.shape
    for y in range(Ny):
        for x in range(Nx):
            dx = x - x0
            dy = y - y0
            chemB[y, x] += amplitude * np.exp(-(dx * dx + dy * dy) / (2 * sigma * sigma))


def make_circular_kernel(R: float, dx: float = 1.0) -> NDArray:
    radius_cells = max(1, int(np.ceil(R / dx)))
    offsets = np.arange(-radius_cells, radius_cells + 1, dtype=float)
    yy, xx = np.meshgrid(offsets, offsets, indexing="ij")
    footprint = ((xx * dx) ** 2 + (yy * dx) ** 2) <= R * R
    kernel = footprint.astype(np.float64)
    kernel_sum = kernel.sum()
    if kernel_sum == 0:
        kernel[radius_cells, radius_cells] = 1.0
        kernel_sum = 1.0
    return kernel / kernel_sum

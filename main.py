import matplotlib.pyplot as plt
import numpy as np
import os
from numpy.typing import NDArray
from tqdm import tqdm

from src.particle import Particle
from src.fields import *
from src import utils
from src.utils import populate_is_occupied


def folder():
    folder_name = 'snapshots'
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    return folder_name


def snap(is_occupied: NDArray,
         chemA: NDArray, chemB: NDArray,
         path: str, i: int) -> None:
    fig, ax = plt.subplots(1, 3)
    ax = ax.flatten()
    ax[0].set_title(f'at_{i}_is_occupied')
    ax[0].imshow(is_occupied)

    ax[1].set_title(f'at_{i}_chemA')
    ax[1].imshow(chemA)

    ax[2].set_title(f'at_{i}_chemB')
    ax[2].imshow(chemB)

    fig.savefig(os.path.join(path, 'at_{i}_snapshot.png'))
    plt.close(fig)


def main(is_occupied: NDArray,
         chemA: NDArray, chemB: NDArray,
         positions: NDArray, particles: list[Particle],
         collisions: list[tuple[int, int]],
         num_particles: int, num_steps: int
         ) -> None:
    for i in tqdm(range(num_steps)):
        is_occupied.fill(-1)
        collisions.clear()
        new_theta = []

        avgA = utils.local_field_avg(chemA, R=10)
        avgB = utils.local_field_avg(chemB, R=10)

        dAdx, dAdy = utils.gradient(avgA)
        dBdx, dBdy = utils.gradient(avgB)
        chemA.fill(1.0)

        # --- deposit ---
        for particle in particles:
            deposit_chemA(particle, chemA)
            deposit_chemB(particle, chemB)

        # --- particle update ---
        for particle in particles:
            utils.populate_is_occupied(radius=particle.radius, xc=particle.x, yc=particle.y,
                                       id=particle.id,
                                       is_occupied=is_occupied, collisions=collisions)
        # Compute theta's
        for particle in particles:
            theta = utils.update_theta(particle=particle, particles=particles,
                                       gradAx=dAdx, gradAy=dAdy, gradBx=dBdx, gradBy=dBdy)
            new_theta.append(theta)

        # assign theta's
        for idx, particle in enumerate(particles):
            particle.theta = new_theta[idx]

        for particle in particles:
            utils.move(particle, Nx, Ny)
        for particle in particles:
            positions[i, particle.id - 1, 0] = particle.x_unwrapped
            positions[i, particle.id - 1, 1] = particle.y_unwrapped

        # --- diffuse ---
        chemA = utils.laplacian(chemA)
        chemB = utils.laplacian(chemB)

        # --- decay ---
        chemB = decay(chemB, rate=0.01)
        if i % 10 == 0:
            path = folder()
            snap(is_occupied=is_occupied, chemA=chemA, chemB=chemB,
                 path=path, i=i)


pass

if __name__ == '__main__':

    L = 10
    BASE_NX = 100
    BASE_NY = 100
    Nx = L * BASE_NX
    Ny = L * BASE_NY
    num_particles = 4
    num_steps = 20

    is_occupied = np.full((Ny, Nx), fill_value=-1, dtype=np.int16)
    chemA = np.full((Ny, Nx), fill_value=0.0, dtype=np.float32)
    chemB = np.full((Ny, Nx), fill_value=0.0, dtype=np.float32)
    positions = np.full(shape=(num_steps, num_particles, 2), fill_value=np.nan)

    particles = []
    collisions = []

    #  Initialize particles with random positions, orientations, and radii
    for i in range(1, num_particles + 1):
        x = Nx * np.random.rand()
        y = Ny * np.random.rand()
        theta = 2 * np.pi * np.random.rand()
        radius = 20.0 + 5. * np.random.rand()
        particles.append(Particle(i, x, y, theta, radius))
        populate_is_occupied(radius=radius, xc=x, yc=y, id=i,
                             is_occupied=is_occupied, collisions=collisions)

        for particle in particles:
            deposit_chemA(particle=particle, chemA=chemA, is_occupied=is_occupied)
            deposit_chemB(particle=particle, chemB=chemB)

        main(is_occupied=is_occupied,
             chemA=chemA, chemB=chemB,
             positions=positions, particles=particles,
             collisions=collisions,
             num_particles=num_particles, num_steps=num_steps)

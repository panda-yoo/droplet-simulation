"""
canonical_correlation_analysis.py

Utilities for singular value analysis of cross-correlation matrices.
"""

from typing import Tuple

import numpy as np
from numpy.typing import NDArray


def compute_singular_values(
    correlation_matrices: NDArray[np.floating],
) -> NDArray[np.floating]:
    """
    Compute singular values of every correlation matrix.

    Parameters
    ----------
    correlation_matrices : ndarray
        Shape = (n_lags, n, n)

    Returns
    -------
    ndarray
        Shape = (n_lags, n)
    """

    n_lags = correlation_matrices.shape[0]
    n_modes = correlation_matrices.shape[1]

    singular_values = np.empty((n_lags, n_modes))

    for tau in range(n_lags):

        M = correlation_matrices[tau].copy()

        # Ignore trivial self-correlations
        np.fill_diagonal(M, 0)

        singular_values[tau] = np.linalg.svd(
            M,
            compute_uv=False
        )

    return singular_values


def compute_largest_singular_values(
    correlation_matrices: NDArray[np.floating],
) -> NDArray[np.floating]:
    """
    Largest singular value as a function of lag.
    """

    singular_values = compute_singular_values(correlation_matrices)

    return singular_values[:, 0]


def compute_singular_vectors(
    correlation_matrices: NDArray[np.floating],
) -> Tuple[
    NDArray[np.floating],
    NDArray[np.floating],
    NDArray[np.floating],
]:
    """
    Compute SVD of every correlation matrix.

    Returns
    -------
    U : (n_lags,n,n)
    S : (n_lags,n)
    Vt : (n_lags,n,n)
    """

    n_lags = correlation_matrices.shape[0]
    n = correlation_matrices.shape[1]

    U = np.empty((n_lags, n, n))
    S = np.empty((n_lags, n))
    Vt = np.empty((n_lags, n, n))

    for tau in range(n_lags):

        M = correlation_matrices[tau].copy()

        np.fill_diagonal(M, 0)

        U[tau], S[tau], Vt[tau] = np.linalg.svd(M)

    return U, S, Vt


def compute_mode_energy(
    singular_values: NDArray[np.floating],
) -> NDArray[np.floating]:
    """
    Fraction of total energy carried by each singular mode.

    E_i = sigma_i^2 / sum_j sigma_j^2
    """

    sigma2 = singular_values**2

    return sigma2 / sigma2.sum(axis=1, keepdims=True)


def compute_effective_rank(
    singular_values: NDArray[np.floating],
) -> NDArray[np.floating]:
    """
    Effective rank based on spectral entropy.

    r_eff = exp( -Σ p_i log(p_i) )
    """

    p = singular_values / singular_values.sum(
        axis=1,
        keepdims=True,
    )

    eps = np.finfo(float).eps

    entropy = -np.sum(
        p * np.log(p + eps),
        axis=1,
    )

    return np.exp(entropy)
from analysis.vacf import compute_velocity
from analysis.types import PositionsDict


from typing import Tuple,Dict
import numpy as np
from numpy.typing import NDArray
import matplotlib.pyplot as plt


def compute_cross_vcaf(v1:NDArray[np.floating],
                    v2:NDArray[np.floating],
                    lag_time:int = None)->tuple[NDArray[np.floating],int]:

    vi = compute_velocity(v1, dt=0.2)
    vj = compute_velocity(v2, dt=0.2)

    n = min(len(vi), len(vj))

    if lag_time is None:
        lag_time = n - 1
        
    cij = np.zeros(lag_time)

    norm_i = np.mean(vi[:,0]**2 + vi[:,1]**2)
    norm_j = np.mean(vj[:,0]**2 + vj[:,1]**2)

    norm = np.sqrt(norm_i * norm_j)

    for tau in range(lag_time):

        corr = (
            vi[:n-tau,0] * vj[tau:,0]
            +
            vi[:n-tau,1] * vj[tau:,1]
        )

        cij[tau] = np.mean(corr) / norm
        
    return cij,lag_time



def compute_cross_correlation_matrix(trajectories: PositionsDict,lagtime:int)->NDArray[np.floating]:
    """computes a cross-correlation matrix for all tau from 0 to lagtime

    Args:
        trajectories (PositionsDict): _description_
        lagtime (int): _description_
        is_max (bool, optional): _description_. Defaults to True.

    Returns:
        NDArray[np.floating]: _description_
    """
    nsize = len(trajectories.keys())
    m = np.zeros((lagtime,nsize,nsize))

    for i,key1 in enumerate(trajectories.keys()):
        for j,key2 in enumerate(trajectories.keys()):

            vi = trajectories[key1]
            vj = trajectories[key2]
            cij,_ = compute_cross_vcaf(v1 = vi,v2= vj,lag_time=lagtime)
    
            m[:,i,j] = cij.squeeze()
       
    return m

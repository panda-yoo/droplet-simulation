import numpy as np
from analysis.extract import load_positions   

from typing import Tuple
from analysis.types import PositionsDict,PairDistanceData


def compute_inter_distance_squared(trajectories: PositionsDict) -> PairDistanceData:
    """_summary_

    Args:
        trajectories (PositionsDict): _description_

    Returns:
        PairDistanceData: _description_
    """
    
    n = len(trajectories)
    trange = trajectories[list(trajectories.keys())[0]].shape[0]
    dro_ids = list(trajectories.keys())
    
    pair_map = []

    for i in range(n):
        for j in range(i+1,n):
            pair_map.append((dro_ids[i], dro_ids[j]))
    
    dij = np.zeros((trange, len(pair_map)))
    for t in range(trange):
        
        for k,(i,j) in enumerate(pair_map):
                    
            r1 = trajectories[i]
            r2 = trajectories[j]
            
            dr = r1[t] - r2[t]
            dij[t, k] = np.linalg.norm(dr)
     
                
    return PairDistanceData(distances= dij,
                            pair_map=pair_map)


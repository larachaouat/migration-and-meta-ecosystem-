
# This script
# - wraps repeated model runs of a 3-compartment ecosystem model (migration)
# - scans two parameters across a 2D grid and records persistence, mean stocks and timing indices
# - intended for use in figures / supplementary material for publication
#

import numpy as np
from scipy.integrate import solve_ivp
from class_meta_eco_mig_3C_lin import two_eco_s_3C_lin
import matplotlib.pyplot as plt


def calculate_tot_production_3C(I1, e1, a1, amp1, F, F1, F2, h1, h2, rc, rp, mp, mc1, mc2, nb_years, X0, f,  s, ext, B, Fmig=None, amp2=None, var1='I1', var2='a1', range1=None, range2=None):    # Set default ranges if none are provided
    """
    Calculate persistence and summary statistics over a 2D parameter grid for the 3-component ecosystem model
    with and without migration (migration scenario).
    This function sweeps two parameters (var1 and var2) over specified ranges (or defaults) to build a B x B grid,
    solves the underlying ecological model for each parameter combination both with migration and without migration,
    and computes persistence categories plus mean/min/max stock summaries for several state variables.
    Parameters
    ----------
    I1, e1, a1, amp1, F, F1, F2, h1, h2, rc, rp, mp, mc1, mc2 : float
        Model parameters passed directly to two_eco_s_3C_lin. Meaning depends on the model
        implementation (ingestion rates, efficiencies, attack rates, amplitudes, feeding rates, handling times,
        mortality/competition coefficients, etc.).
    nb_years : float or int
        Number of years to simulate. Used together with f to determine indexing into the returned time series.
    X0 : sequence of float
        Initial state vector passed to the model solvers. Must match the state ordering expected by two_eco_s_3C_lin.
    f : float
        Frequency parameter used to compute the time window for statistics (the function computes a time index
        window near the simulation end using nb_years and f). The code assumes model time-series are sampled at
        10 samples per unit time (see Notes).
    s : float
        A parameter forwarded into the model variable dictionary (purpose depends on model).
    ext : float
        Extinction threshold. Max values below ext are considered extinct in persistence classification.
    B : int
        Number of grid points per axis for the parameter sweep (grid is B x B).
    Fmig : float, optional
        Default migration fraction used when building the model variable dictionary (default None in signature,
        but the function sets a default of 0.5 internally).
    amp2 : float, optional
        Secondary amplitude parameter; if None the function sets amp2 = 0 internally.
    var1, var2 : str, optional
        Names of the two parameters to sweep in the variables dictionary. By default var1='I1' and var2='a1'.
    range1, range2 : array-like, optional
        Sequences of length B giving the values to sweep for var1 and var2. If None, they default to
        np.linspace(0, 2, B).

    Returns
    -------
    p_list, p1_list, p2_list : (B,B) arrays
        Persistence indicators:
            0 -> extinct in final window (below `ext`)
            3 -> persistent (above `ext`)
        p_list reports persistence of the combined consumer (C + Ca as used in model output),
        p1_list for consumer 1, p2_list for consumer 2.
    b_list_m : (8,B,B) array
        Mean final-window stocks for the 8 tracked quantities (order used in this code):
         0: N1 mean, 1: N2 mean, 2: P1 mean, 3: P2 mean,
         4: C mean, 5: Ca mean, 6: C2 mean, 7: C1 mean
    max_idx_list : (5,B,B) array
        Indices (relative to the sampled final window) of maxima for:
         [0] maximum uptake index,
         [1] unused (kept for potential extension),
         [2] maximum combined consumer (C) index,
         [3] unused,
         [4] maximum P1 index
        If no maximum was found (empty), index is set to 0.
    """
    
    
    if range1 is None:
        range1 = np.linspace(0, 2, B)
    if range2 is None:
        range2 = np.linspace(0, 2, B)
    it=0
    b_list = np.zeros((8, B, B))
    b_list_m = np.zeros((8, B, B))
    max_idx_list = np.zeros((5, B, B))
    t_end = int(nb_years * 1 / f)
    t_end_m = int(nb_years * 1 / f)
    p_list = np.zeros((B,B))
    p2_list = np.zeros((B,B))
    p1_list = np.zeros((B,B))

    X=X0
    X_m=X0
    for i, x_i in enumerate(range1):
        for j, x_j in enumerate(range2):
            variables = {
                'I1': I1, 'e1': e1, 
                'a1': a1, 
                'F': F, 'F1': F1, 'F2':F2,
                'h1': h1, 'h2': h2, 'rc': rc, 'rp': rp, 'mp': mp,
                'mc1': mc1, 'mc2': mc2, 'f': f, 'amp1': amp1, 
                'amp2':0, 's': s, 'Fmig': .5, 'per_s1': .5, 'per_s2': .5, 'mean_norm':False
            }

            # Update the specific variables to be varied
            variables[var1] = x_i
            variables[var2] = x_j
            variables['Fmig'] = x_j
            print(variables)

            model = two_eco_s_3C_lin(**variables)
            time, N1, P1, C1, C, N2, P2, C2, Ca = model.solve_model_nomig(nb_years=nb_years, X0=X0, eco=1)
            t_m, N1_m, P1_m, C1_m, C_m, N2_m, P2_m, C2_m, Ca_m = model.solve_model(nb_years=nb_years, X0=[X0[0],X0[1],X0[2], X0[7], X0[4],X0[5],X0[6], X0[3]], eco=2)

####persistence
            max_c = np.max(C[int(t_end - 1 / f) * 10:t_end * 10 - 2])
            max_c_m = np.max(C_m[int(t_end_m - 1 / f) * 10:t_end_m * 10 - 2] + Ca_m[int(t_end_m - 1 / f) * 10:t_end_m * 10 - 2])
            max_c_2 = np.max(C2[int(t_end - 1 / f) * 10:t_end * 10 - 2])
            max_c_2m = np.max(C2_m[int(t_end_m - 1 / f) * 10:t_end_m * 10 - 2])
            max_c_1 = np.max(C1[int(t_end - 1 / f) * 10:t_end * 10 - 2])
            max_c_1m = np.max(C1_m[int(t_end_m - 1 / f) * 10:t_end_m * 10 - 2])
            
            if max_c < ext and max_c_m < ext:  # non 
                p_list[i, j] = 0
            elif max_c > ext and max_c_m  < ext:  # C1 non migration survives
                p_list[i, j] = 1
            elif max_c < ext and max_c_m  > ext:  # C1 migrator survives
                p_list[i, j] = 2
            else:  # both case survival
                p_list[i, j] = 3
                
            if max_c_2 < ext and max_c_2m < ext:# non 
                p2_list[i,j] = 0
            elif max_c_2 > ext and max_c_2m  < ext: # C1 non migration survives
                p2_list[i,j] = 1
            elif max_c_2 < ext and max_c_2m  > ext:  # C1  migrator survives
                p2_list[i,j] = 2
            else: #both case survival
                p2_list[i,j] = 3

            if max_c_1 < ext and max_c_1m < ext:# non 
                p1_list[i,j] = 0
            elif max_c_1 > ext and max_c_1m  < ext: # C1 non migration survives
                p1_list[i,j] = 1
            elif max_c_1 < ext and max_c_1m  > ext:  # C1  migrator survives
                p1_list[i,j] = 2
            else: #both case survival
                p1_list[i,j] = 3

####stocks
# for N
            b_list[0,i,j] = np.mean(N1[int(t_end - 1 / f) * 10:t_end * 10 - 2])
            b_list[1,i,j] = np.mean(N2[int(t_end - 1 / f) * 10:t_end * 10 - 2])
            b_list_m[0,i,j]= np.mean(N1_m[int(t_end_m - 1 / f) * 10:t_end_m * 10 - 2])
            b_list_m[1,i,j] = np.mean(N2_m[int(t_end_m - 1 / f) * 10:t_end_m * 10 - 2])

#for P      
            b_list[2,i,j] = np.mean(P1[int(t_end - 1 / f) * 10:t_end * 10 - 2])
            b_list[3,i,j] = np.mean(P2[int(t_end - 1 / f) * 10:t_end * 10 - 2])
            b_list_m[2,i,j] = np.mean(P1_m[int(t_end_m - 1 / f) * 10:t_end_m * 10 - 2])
            b_list_m[3,i,j] = np.mean(P2_m[int(t_end_m - 1 / f) * 10:t_end_m * 10 - 2])
#for C
            b_list[4,i,j] = np.mean(C[int(t_end - 1 / f) * 10:t_end * 10 - 2])
            b_list[5,i,j] = np.mean(Ca[int(t_end - 1 / f) * 10:t_end * 10 - 2])
            b_list[6,i,j] = np.mean(C2[int(t_end - 1 / f) * 10:t_end * 10 - 2])
            b_list[7,i,j] = np.mean(C1[int(t_end - 1 / f) * 10:t_end * 10 - 2])
            b_list_m[4,i,j] = np.mean(C_m[int(t_end_m - 1 / f) * 10:t_end_m * 10 - 2])
            b_list_m[5,i,j] = np.mean(Ca_m[int(t_end_m - 1 / f) * 10:t_end_m * 10 - 2])
            b_list_m[6,i,j] = np.mean(C2_m[int(t_end_m - 1 / f) * 10:t_end_m * 10 - 2])
            b_list_m[7,i,j] = np.mean(C1_m[int(t_end_m - 1 / f) * 10:t_end_m * 10 - 2])
#timing max
            uptake = variables['a1']*(1+variables['amp1']*model.sin_seasonality(time,variables['per_s1'],model.fact_season_extrem( variables['per_s1'])))
            max_uptake_idx = np.where( uptake[int(t_end - 1 / f) * 10:int(t_end * 10 - 2)]== np.max(uptake[int(t_end - 1 / f) * 10:int(t_end * 10 - 2)]))
            max_P1_idx =  np.where( P1[int(t_end - 1 / f) * 10:int(t_end * 10 - 2)] == np.max(P1[int(t_end - 1 / f) * 10:int(t_end * 10 - 2)]))
            max_P1_m_idx =  np.where( P1_m[int(t_end - 1 / f) * 10:int(t_end * 10 - 2)] == np.max(P1_m[int(t_end - 1 / f) * 10:int(t_end * 10 - 2)]))
            max_C_idx =  np.where( C[int(t_end - 1 / f) * 10:int(t_end * 10 - 2)] == max_c)
            max_C_m_idx =  np.where( C_m[int(t_end - 1 / f) * 10:int(t_end * 10 - 2)] == np.max(C_m[int(t_end - 1 / f) * 10:int(t_end * 10 - 2)] ))

            print(max_c, max_c_m, max_C_idx,  max_C_m_idx)

            max_idx_list[0,i,j] = max_uptake_idx[0][0] if max_uptake_idx[0].size > 0 else  0 #maximum uptake index
            max_idx_list[1,i,j] = max_C_idx[0][0] if max_C_idx[0].size > 0 else  0  #maximum non mig conso index
            max_idx_list[2,i,j] = max_C_m_idx[0][0] if max_C_m_idx[0].size > 0  else  0    #maximum mig conso index
            max_idx_list[3,i,j] = max_P1_idx[0][0] if max_P1_idx[0].size > 0 else  0  #maximum non mig conso index
            max_idx_list[4,i,j] = max_P1_m_idx[0][0] if max_P1_m_idx[0].size > 0  else  0    #maximum mig conso index

            it += 1
            print(it)


    return p_list, p1_list, p2_list, b_list, b_list_m, max_idx_list



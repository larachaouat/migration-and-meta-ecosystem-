
# This script
# - wraps repeated model runs of a 3-compartment ecosystem model (no migration)
# - scans two parameters across a 2D grid and records persistence, mean stocks and timing indices
# - intended for use in figures / supplementary material for publication
#


import numpy as np
from scipy.integrate import solve_ivp
from class_meta_eco_mig_3C_lin import two_eco_s_3C_lin
import matplotlib.pyplot as plt


def calculate_tot_production_3C_nomig(
    I1, e1, a1, amp1, F, F1, F2, h1, h2, rc, rp, mp, mc1, mc2,
    nb_years, X0, f, s, ext, B, Fmig=None, amp2=None,
    var1='I1', var2='a1', range1=None, range2=None
):
    """
    Run the two_eco_s_3C_lin model across a 2D parameter grid and collect persistence,
    mean stocks and timing of maxima.

    Parameters
    ----------
    I1, e1, a1, amp1, F, F1, F2, h1, h2, rc, rp, mp, mc1, mc2 : floats
        Base parameter values passed to the model. Individual parameters may be
        overridden by the scanning variables var1/var2 and the arrays range1/range2.
    nb_years : float or int
        Number of simulated years (used to compute time windows for statistics).
    X0 : array-like
        Initial condition vector passed to model.solve_model_nomig.
    f : float
        Frequency (per year) used to compute the time indexing for the final-year window.
        The code assumes model outputs are sampled at a resolution such that multiplying
        times by 10 is appropriate for indexing (see note below).
    s : float
        Seasonality phase/parameter passed to the model (kept here for completeness).
    ext : float
        Extinction threshold. If consumer biomass < ext in the final window, species is
        considered extinct (p_list entry 0). Otherwise marked as persisting (3).
    B : int
        Grid resolution for the scan. If range1/range2 are None, default linspace(0,2,B)
        is used for each scanned variable.
    Fmig, amp2 : optional
        Unused in the current function but accepted for API compatibility.
    var1, var2 : str
        Names of the parameters that will be varied on the two axes (must match keys
        expected by two_eco_s_3C_lin).
    range1, range2 : array-like, optional
        Explicit arrays of values to scan for var1 and var2. If None, default ranges are used.

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

    # Set default ranges if none are provided
    if range1 is None:
        range1 = np.linspace(0, 2, B)
    if range2 is None:
        range2 = np.linspace(0, 2, B)

    # Preallocate output containers
    it = 0
    b_list = np.zeros((8, B, B))        # not filled in this function version, kept for API compatibility
    b_list_m = np.zeros((8, B, B))      # mean stocks for the "nomig" runs
    max_idx_list = np.zeros((5, B, B))  # indices of maxima recorded in the final-year window
    # t_end and t_end_m represent the integer number of sampling-blocks (years scaled by f)
    t_end = int(nb_years * 1 / f)
    t_end_m = int(nb_years * 1 / f)

    # Persistence matrices for combined consumer (p_list), consumer1 (p1_list), consumer2 (p2_list)
    p_list = np.zeros((B, B))
    p2_list = np.zeros((B, B))
    p1_list = np.zeros((B, B))

    # X/X_m are not used intensively in this function, but kept for potential future state updates
    X = X0
    X_m = X0

    # Loop over the 2D grid of parameter values
    for i, x_i in enumerate(range1):
        for j, x_j in enumerate(range2):
            # Build parameters dict for this run. Keep defaults for parameters we don't scan.
            variables = {
                'I1': I1, 'e1': e1,
                'a1': a1,
                'F': F, 'F1': F1, 'F2': F2,
                'h1': h1, 'h2': h2, 'rc': rc, 'rp': rp, 'mp': mp,
                'mc1': mc1, 'mc2': mc2, 'f': f, 'amp1': amp1,
                'amp2': 0, 's': s, 'Fmig': .5, 'per_s1': .5, 'per_s2': .5, 'mean_norm': False
            }

            # Update the two parameters being scanned for this grid cell
            variables[var1] = x_i
            variables[var2] = x_j

            # For debugging and reproducibility, print the parameter set being evaluated
            print(variables)

            # Instantiate the model with these parameters
            model = two_eco_s_3C_lin(**variables)

            # Solve the model for the "no migration" case.
            # The returned arrays are time series for the state variables used below.
            # Note: this call expects model.solve_model_nomig to accept nb_years, X0 and eco.
            t_m, N1_m, P1_m, C1_m, C_m, N2_m, P2_m, C2_m, Ca_m = model.solve_model_nomig(
                nb_years=nb_years, X0=X0, eco=1
            )

            # --- PERSISTENCE CHECKS ---
            # The code computes maxima of the final-year window and compares with extinction threshold.
            # The slicing uses multipliers of 10 (e.g. int(t_end_m - 1 / f) * 10) — this assumes that
            # model output is sampled at a resolution such that multiplying time indices by 10 aligns with samples.
            # If sampling resolution changes, these index computations must be adapted.

            max_c_m = np.max(
                C_m[int(t_end_m - 1 / f) * 10:t_end_m * 10 - 2]
                + Ca_m[int(t_end_m - 1 / f) * 10:t_end_m * 10 - 2]
            )
            max_c_2m = np.max(C2_m[int(t_end_m - 1 / f) * 10:t_end_m * 10 - 2])
            max_c_1m = np.max(C1_m[int(t_end_m - 1 / f) * 10:t_end_m * 10 - 2])

            # Mark persistence (0 = extinct, 3 = persisting)
            if max_c_m < ext:
                p_list[i, j] = 0
            else:
                p_list[i, j] = 3

            if max_c_2m < ext:
                p2_list[i, j] = 0
            else:
                p2_list[i, j] = 3

            if max_c_1m < ext:
                p1_list[i, j] = 0
            else:
                p1_list[i, j] = 3

            # --- MEAN STOCKS IN FINAL WINDOW ---
            # Compute mean values of each tracked variable within the final sampled window.
            # The indices are identical to those used for persistence above.
            b_list_m[0, i, j] = np.mean(N1_m[int(t_end_m - 1 / f) * 10:t_end_m * 10 - 2])
            b_list_m[1, i, j] = np.mean(N2_m[int(t_end_m - 1 / f) * 10:t_end_m * 10 - 2])

            b_list_m[2, i, j] = np.mean(P1_m[int(t_end_m - 1 / f) * 10:t_end_m * 10 - 2])
            b_list_m[3, i, j] = np.mean(P2_m[int(t_end_m - 1 / f) * 10:t_end_m * 10 - 2])

            b_list_m[4, i, j] = np.mean(C_m[int(t_end_m - 1 / f) * 10:t_end_m * 10 - 2])
            b_list_m[5, i, j] = np.mean(Ca_m[int(t_end_m - 1 / f) * 10:t_end_m * 10 - 2])
            b_list_m[6, i, j] = np.mean(C2_m[int(t_end_m - 1 / f) * 10:t_end_m * 10 - 2])
            b_list_m[7, i, j] = np.mean(C1_m[int(t_end_m - 1 / f) * 10:t_end_m * 10 - 2])

            # --- TIMING OF MAXIMA (indices within the final window) ---
            # Compute the uptake time series (using model's seasonality helper methods)
            uptake = variables['a1'] * (
                1 + variables['amp1'] * model.sin_seasonality(
                    t_m, variables['per_s1'], model.fact_season_extrem(variables['per_s1'])
                )
            )

            # Find index (relative to sliced final-window) where uptake is maximal
            max_uptake_idx = np.where(
                uptake[int(t_end - 1 / f) * 10:int(t_end * 10 - 2)]
                == np.max(uptake[int(t_end - 1 / f) * 10:int(t_end * 10 - 2)])
            )

            # Index of maximum P1 in the final window
            max_P1_m_idx = np.where(
                P1_m[int(t_end - 1 / f) * 10:int(t_end * 10 - 2)]
                == np.max(P1_m[int(t_end - 1 / f) * 10:int(t_end * 10 - 2)])
            )

            # Index of maximum combined consumer C in the final window
            max_C_m_idx = np.where(
                C_m[int(t_end - 1 / f) * 10:int(t_end * 10 - 2)]
                == np.max(C_m[int(t_end - 1 / f) * 10:int(t_end * 10 - 2)])
            )

            # Store the first index found (or 0 if none)
            max_idx_list[0, i, j] = max_uptake_idx[0][0] if max_uptake_idx[0].size > 0 else 0
            max_idx_list[2, i, j] = max_C_m_idx[0][0] if max_C_m_idx[0].size > 0 else 0
            max_idx_list[4, i, j] = max_P1_m_idx[0][0] if max_P1_m_idx[0].size > 0 else 0

            it += 1
            print(it)

    return p_list, p1_list, p2_list, b_list_m, max_idx_list

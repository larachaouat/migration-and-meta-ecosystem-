
# This script runs simulations of a meta-ecosystem model with three consummers for figures 5,6
# and writes existence / production results to CSV files 
# For each combination of two varying parameters, three model variants are
# evaluated:
#   - calculate_tot_production_3C_nomig : no migration
#   - calculate_tot_production_3C       : full model (with migration)
#   - calculate_tot_production_3C_migfix : migration fixed for some components
#
# Each call returns arrays for:
#   p_list    : overall species existence (or occupancy) matrix
#   p1_list   : existence / occupancy for consumer 1
#   p2_list   : existence / occupancy for consumer 2
#   b_list    : biomass / total production (multi-dimensional)
#   max_idx_list : index of the dominant state / configuration
#
# The script reshapes b_list and max_idx_list to 2-D and saves all outputs
# as CSV files named by the varying parameter pair and model info.
#
# Notes for reproducibility:
#   - threshold is used as an extinction threshold when deciding existence
#   - B is the grid resolution for the varying parameter(s)
#   - nb_years is the length of each simulation (used by the model functions)
#   - X0 is the initial state vector for the model
#   - range1 and range2 are numeric vectors defining the parameter ranges
#     (typically generated with np.linspace)
#
# Usage:
#   run this script from the directory that contains the imported functions
#   fun_biom_prod3c_f_m.py, fun_biom_prod3c_f_migfix_m.py, fun_biom_prod3c_f_nomig.py
#   or adjust import paths accordingly.

import matplotlib.pyplot as plt          # plotting library (not used directly in this script, but likely used in helper functions)
import matplotlib.patches as mpatches    # patch shapes for matplotlib (not used directly here)
import numpy as np                       # numerical arrays and file IO utilities
from fun_biom_prod3c_f_migfix import calculate_tot_production_3C_migfix

# Number of discrete steps used when creating parameter ranges (B x B grid)
B = 8

# Extinction threshold used by the model (populations below this considered extinct)
threshold = 1e-3

set = [
    ['F', 'Fmig', np.linspace(.3, .9, B), np.linspace(0.3, .7, B)] 

]

# ---------------------------
# Model parameter definitions
# ---------------------------
I1, e1, a1, amp1, amp2, F, F1, F2, h1, h2, rc, rp, mp, mc1, mc2, nb_years, X0, f, s, mean_norm = .8, .5, 1.5, .8, 0, .6, .6, .6, 5, 5, 0, 0, .5, .1, .1, 10, [1, 1, 0,
                                                                                                                                                                   0,
                                                                                                                                                                   1, 1, 0,
                                                                                                                                                                   1], 1/1000, 1, False


# ---------------------------
# Main loop: run each model variant for every parameter combination
# ---------------------------
for var1, var2, range1, range2 in set:
    # Print current initial condition vector to stdout for logging
    print(X0)

    # --- First scenario: "1C_amp8" ---
    # info string used to label output files for this scenario
    info = '1C_amp8'  # label for output files from this run

    # Call the simulation function that includes migration-fix behavior.
    # The function returns many arrays:
    # p_list, p1_list, p2_list: presence/existence metrics (per grid point)
    # b_list, b_list_m: total production arrays (resident and migrant)
    # b_list_min, b_list_min_m, b_list_max, b_list_max_m: production min/max variants
    p_list, p1_list, p2_list, b_list, b_list_m, b_list_min, b_list_min_m, b_list_max, b_list_max_m = calculate_tot_production_3C_migfix(
        I1, e1, a1, amp1, F, F1, F2, h1, h2, rc, rp, mp, mc1, mc2,
        nb_years, X0, f, s,
        amp2=amp2, ext=threshold, B=B, var1=var1, var2=var2, range1=range1, range2=range2
    )

    # Reshape multi-dimensional production arrays to 2D for saving to CSV.
    # reshape(dim0, -1) flattens all remaining axes into one column per time or species as appropriate.
    b_reshaped = b_list.reshape(b_list.shape[0], -1)
    b_m_reshaped = b_list_m.reshape(b_list.shape[0], -1)
    b_reshaped_min = b_list_min.reshape(b_list.shape[0], -1)
    # Potential mismatch: original code uses b_list_m.shape[0] for b_list_min_m reshape.
    # This preserves original behavior but may be a typo; if shapes differ, this will raise an error.
    b_m_reshaped_min = b_list_min_m.reshape(b_list_m.shape[0], -1)
    b_reshaped_max = b_list_max.reshape(b_list.shape[0], -1)
    b_m_reshaped_max = b_list_max_m.reshape(b_list_max_m.shape[0], -1)

    # Save existence and production results to CSV files using informative filenames.
    # Filenames include var1, var2, scenario info, and number of years simulated (nb_years).
    np.savetxt(f'{var1}vs{var2}_{info}_exist_{nb_years}.csv', p_list, delimiter=",")
    np.savetxt(f'{var1}vs{var2}_{info}_exist2_{nb_years}.csv', p2_list, delimiter=",")
    np.savetxt(f'{var1}vs{var2}_{info}_exist1_{nb_years}.csv', p1_list, delimiter=",")

    np.savetxt(f'{var1}vs{var2}_{info}_prod_{nb_years}.csv', b_reshaped, delimiter=",")
    np.savetxt(f'{var1}vs{var2}_{info}_prod_m_{nb_years}.csv', b_m_reshaped, delimiter=",")
    np.savetxt(f'{var1}vs{var2}_{info}_prod_min_{nb_years}.csv', b_reshaped_min, delimiter=",")
    np.savetxt(f'{var1}vs{var2}_{info}_prod_m_min_{nb_years}.csv', b_m_reshaped_min, delimiter=",")
    np.savetxt(f'{var1}vs{var2}_{info}_prod_max_{nb_years}.csv', b_reshaped_max, delimiter=",")
    np.savetxt(f'{var1}vs{var2}_{info}_prod_m_max_{nb_years}.csv', b_m_reshaped_max, delimiter=",")

    # --- Second scenario: "2C_F1.6" ---
    # Change initial condition vector X0 to represent a different community composition.
    X0 = [1, 1, 1,
          0,
          1, 1, 0,
          1]
    info = '2C_F1.6'  # label for this scenario

    # Call the simulation function again with the new initial conditions.
    # This call returns slightly different outputs in the original script: here max_idx_list is returned.
    p_list, p1_list, p2_list, b_list, b_list_m, max_idx_list = calculate_tot_production_3C_migfix(
        I1, e1, a1, amp1, F, F1, F2, h1, h2, rc, rp, mp, mc1, mc2,
        nb_years, X0, f, s,
        amp2=amp2, ext=threshold, B=B, var1=var1, var2=var2, range1=range1, range2=range2
    )

    # Reshape production arrays and the max index array for saving.
    b_reshaped = b_list.reshape(b_list.shape[0], -1)
    b_m_reshaped = b_list_m.reshape(b_list.shape[0], -1)
    max_idx_list_reshaped = max_idx_list.reshape(max_idx_list.shape[0], -1)

    # Save results to CSV with scenario-specific filenames.
    np.savetxt(f'{var1}vs{var2}_{info}_exist_{nb_years}.csv', p_list, delimiter=",")
    np.savetxt(f'{var1}vs{var2}_{info}_exist1_{nb_years}.csv', p1_list, delimiter=",")
    np.savetxt(f'{var1}vs{var2}_{info}_exist2_{nb_years}.csv', p2_list, delimiter=",")
    np.savetxt(f'{var1}vs{var2}_{info}_prod_{nb_years}.csv', b_reshaped, delimiter=",")
    np.savetxt(f'{var1}vs{var2}_{info}_prod_m_{nb_years}.csv', b_m_reshaped, delimiter=",")
    np.savetxt(f'{var1}vs{var2}_{info}_max_idx_{nb_years}.csv', max_idx_list_reshaped, delimiter=",")

    # --- Third scenario: "3C_F1.6_F2.6_amp8" ---
    info = '3C_F1.6_F2.6_amp8'  # label for this scenario

    # Run the simulation again (similar to first call) to produce min/max variants as well.
    p_list, p1_list, p2_list, b_list, b_list_m, b_list_min, b_list_min_m, b_list_max, b_list_max_m = calculate_tot_production_3C_migfix(
        I1, e1, a1, amp1, F, F1, F2, h1, h2, rc, rp, mp, mc1, mc2,
        nb_years, X0, f, s,
        amp2=amp2, ext=threshold, B=B, var1=var1, var2=var2, range1=range1, range2=range2
    )

    # Reshape arrays as before for CSV export.
    b_reshaped = b_list.reshape(b_list.shape[0], -1)
    b_m_reshaped = b_list_m.reshape(b_list.shape[0], -1)
    b_reshaped_min = b_list_min.reshape(b_list.shape[0], -1)
    b_m_reshaped_min = b_list_min_m.reshape(b_list_m.shape[0], -1)
    b_reshaped_max = b_list_max.reshape(b_list.shape[0], -1)
    b_m_reshaped_max = b_list_max_m.reshape(b_list_max_m.shape[0], -1)

    # Save outputs to CSV files with descriptive names for downstream analysis/figures.
    np.savetxt(f'{var1}vs{var2}_{info}_exist_{nb_years}.csv', p_list, delimiter=",")
    np.savetxt(f'{var1}vs{var2}_{info}_exist2_{nb_years}.csv', p2_list, delimiter=",")
    np.savetxt(f'{var1}vs{var2}_{info}_exist1_{nb_years}.csv', p1_list, delimiter=",")

    np.savetxt(f'{var1}vs{var2}_{info}_prod_{nb_years}.csv', b_reshaped, delimiter=",")
    np.savetxt(f'{var1}vs{var2}_{info}_prod_m_{nb_years}.csv', b_m_reshaped, delimiter=",")
    np.savetxt(f'{var1}vs{var2}_{info}_prod_min_{nb_years}.csv', b_reshaped_min, delimiter=",")
    np.savetxt(f'{var1}vs{var2}_{info}_prod_m_min_{nb_years}.csv', b_m_reshaped_min, delimiter=",")
    np.savetxt(f'{var1}vs{var2}_{info}_prod_max_{nb_years}.csv', b_reshaped_max, delimiter=",")
    np.savetxt(f'{var1}vs{var2}_{info}_prod_m_max_{nb_years}.csv', b_m_reshaped_max, delimiter=",")




# This script runs simulations of a meta-ecosystem model with one consummer for figures 2
# and writes existence / production results to CSV files 
#
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

import numpy as np
from fun_biom_prod3c_f import calculate_tot_production_3C
from fun_biom_prod3c_f_migfix import calculate_tot_production_3C_migfix


# Default grid resolution used when no B provided to the function call.
B = 5

# Extinction threshold used by the model (passed as 'ext' to calculation routines).
threshold = 1e-6



# ---------------------------
# Model parameter definitions
# ---------------------------
I1, e1, a1, amp1, amp2, F, F1, F2, h1, h2, rc, rp, mp, mc1, mc2, nb_years, X0, f, s, mean_norm= 1, .4, .6, .8 , 0,  .1, .1, .1, 5, 5, 0, 0, .3, .1, .1, 10, [1, 1, 0, 
                                                                                                                                                                 0, 
                                                                                                                                                                 1, 1, 0, 
                                                                                                                                                                 1], 1/1000, 1, False

set = [
    ['amp1', 'per_s1', np.linspace(0.6, 1, B), np.linspace(0.1, .9, B)],

]

# ---------------------------
# Main loop: run each model variant for every parameter combination
# ---------------------------
for var1, var2, range1, range2 in set:
    # ---------- Base model (1C) ----------
    info = '1C'

    # Call the simulation function that includes migration-fix behavior.
    # The function returns many arrays:
    # p_list, p1_list, p2_list: presence/existence metrics (per grid point)
    # b_list, b_list_m: total production arrays (resident and migrant)
    # b_list_min, b_list_min_m, b_list_max, b_list_max_m: production min/max variants
    p_list, p1_list, p2_list, b_list, b_list_m, b_list_min, b_list_min_m, b_list_max, b_list_max_m = calculate_tot_production_3C(
        I1, e1, a1, amp1, F, F1, F2, h1, h2, rc, rp, mp, mc1, mc2,
        nb_years, X0, f, s,
        amp2=amp2, ext=threshold, B=B, var1=var1, var2=var2, range1=range1, range2=range2
    )

    # Some returned arrays are higher-dimensional (e.g., [time, i, j]); reshape
    # them to 2D (rows x columns) so they can be saved as CSV.
    b_reshaped = b_list.reshape(b_list.shape[0], -1)
    b_m_reshaped = b_list_m.reshape(b_list_m.shape[0], -1)

    # Save outputs. Filenames include var names, model info, and nb_years.
    np.savetxt(f'{var1}vs{var2}_{info}_exist_{nb_years}.csv', p_list, delimiter=",")
    np.savetxt(f'{var1}vs{var2}_{info}_exist1_{nb_years}.csv', p1_list, delimiter=",")
    np.savetxt(f'{var1}vs{var2}_{info}_exist2_{nb_years}.csv', p2_list, delimiter=",")
    np.savetxt(f'{var1}vs{var2}_{info}_prod_{nb_years}.csv', b_reshaped, delimiter=",")
    np.savetxt(f'{var1}vs{var2}_{info}_prod_m_{nb_years}.csv', b_m_reshaped, delimiter=",")


    # ---------- Migration-fixed variant (1C_migfix) ----------
    info = '1C_migfix'

    p_list, p1_list, p2_list, b_list, b_list_m, b_list_min, b_list_min_m, b_list_max, b_list_max_m = calculate_tot_production_3C(
        I1, e1, a1, amp1, F, F1, F2, h1, h2, rc, rp, mp, mc1, mc2,
        nb_years, X0, f, s,
        amp2=amp2, ext=threshold, B=B, var1=var1, var2=var2, range1=range1, range2=range2
    )


    # Reshape multi-dimensional outputs before saving.
    b_reshaped = b_list.reshape(b_list.shape[0], -1)
    b_m_reshaped = b_list_m.reshape(b_list_m.shape[0], -1)


    # Save outputs for the migration-fixed variant.
    np.savetxt(f'{var1}vs{var2}_{info}_exist_{nb_years}.csv', p_list, delimiter=",")
    np.savetxt(f'{var1}vs{var2}_{info}_exist1_{nb_years}.csv', p1_list, delimiter=",")
    np.savetxt(f'{var1}vs{var2}_{info}_exist2_{nb_years}.csv', p2_list, delimiter=",")
    np.savetxt(f'{var1}vs{var2}_{info}_prod_{nb_years}.csv', b_reshaped, delimiter=",")
    np.savetxt(f'{var1}vs{var2}_{info}_prod_m_{nb_years}.csv', b_m_reshaped, delimiter=",")


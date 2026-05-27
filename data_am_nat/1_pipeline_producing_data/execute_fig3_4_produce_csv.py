
# Enhanced comments for publication (figures 3 and 4 simulations)
#
# Purpose:
#   Run parameter sweeps for a three-consumer meta-ecosystem model,
#   evaluate three model variants (no migration, full migration, migration fixed),
#   and write existence / production results to CSV for later plotting/analysis.
#
# How it works (high level):
#   For each pair of varying parameters (var1, var2) the script:
#     - calls a model-evaluation function that performs a B x B grid sweep
#       over the supplied numerical ranges (range1, range2).
#     - each call returns:
#         p_list        : overall species existence/occupancy matrix (2D: B x B)
#         p1_list       : existence / occupancy for consumer 1 (2D: B x B)
#         p2_list       : existence / occupancy for consumer 2 (2D: B x B)
#         b_list        : biomass / total production (multi-dimensional)
#         max_idx_list  : index of the dominant state/configuration (multi-dimensional)
#     - b_list and max_idx_list are reshaped to 2-D (rows correspond to grid points)
#       before being saved to CSV so they are easier to load into plotting scripts.
#
# Notes on reproducibility / important parameters:
#   - threshold : extinction threshold (biomass below this considered extinct)
#   - B         : grid resolution (number of points per varying parameter)
#   - nb_years  : length (time) of each simulation run used by the model functions
#   - X0        : initial state vector expected by the model functions (flattened)
#   - range1/range2 : numeric vectors (e.g., np.linspace) defining parameter ranges
#
# File naming convention for outputs:
#   <var1>vs<var2>_<info>_exist_<nb_years>.csv     -> overall occupancy matrix
#   <var1>vs<var2>_<info>_exist1_<nb_years>.csv    -> consumer 1 occupancy
#   <var1>vs<var2>_<info>_exist2_<nb_years>.csv    -> consumer 2 occupancy
#   <var1>vs<var2>_<info>_prod_<nb_years>.csv      -> reshaped production/biomass
#   <var1>vs<var2>_<info>_max_idx_<nb_years>.csv   -> reshaped dominant-state index
#
# Usage:
#   Run this script from the directory that contains the imported functions
#   fun_biom_prod3c_f_m.py, fun_biom_prod3c_f_migfix_m.py, fun_biom_prod3c_f_nomig.py
#   or adjust import paths accordingly.
#
# Citation guidance (for the Methods section in a publication):
#   - State the three model variants and the criteria used to declare a species "present"
#     (extinction threshold). Report B (grid resolution) and nb_years used for
#     time integration. Provide the full parameter vector and the parameter ranges
#     (range1, range2) used in the sweeps. Mention that outputs are saved as CSV
#     for reproducibility and downstream plotting.

import numpy as np
from fun_biom_prod3c_f_m import calculate_tot_production_3C
from fun_biom_prod3c_f_migfix_m import calculate_tot_production_3C_migfix
from fun_biom_prod3c_f_nomig import calculate_tot_production_3C_nomig

# Grid resolution for parameter sweeps (B x B grid)
B = 5

# Extinction threshold: biomass below this is considered extinct / absent
threshold = 1e-3

# ---------------------------
# Model parameter definitions
# ---------------------------
# Core parameters used by all model variants. When reporting methods,
# list these values explicitly and state which were varied vs held fixed.
I1, e1, a1, amp1, amp2, F, F1, F2, h1, h2, rc, rp, mp, mc1, mc2, nb_years, X0, f, s, mean_norm = \
        .8, .5, 1.5, .8, 0, .6, .6, .7, 5, 5, 0, 0, .5, .1, .1, 20, \
        [1, 1, 1,   # initial state vector flattened for model use
         0,
         1, 1, 1,
         1], 1/1000, 1, False

# Define parameter-pair sweeps:
# Each entry: [var1_name, var2_name, numeric_range1, numeric_range2]
# Use descriptive variable names so output filenames communicate what was varied.
parameter_sets = [
        ['amp1', 'F', np.linspace(0, 1, B), np.linspace(.3, .9, B)],
        ['amp1', 'per_s1', np.linspace(0.6, 1, B), np.linspace(0.1, .9, B)]
]

# Loop over each parameter-pair sweep and evaluate model variants.
for var1, var2, range1, range2 in parameter_sets:
        # Print initial state for logging / reproducibility
        print("Initial state X0:", X0)

        # 1) No migration baseline
        info = '3C_nomig_F.6_F2.9_F1.6'  # short descriptor for filenames; document meaning in the paper
        p_list, p1_list, p2_list, b_list, max_idx_list = calculate_tot_production_3C_nomig(
                I1, e1, a1, amp1, F, F1, F2, h1, h2, rc, rp, mp, mc1, mc2,
                nb_years, X0, f, s, amp2=amp2, ext=threshold, B=B,
                var1=var1, var2=var2, range1=range1, range2=range2
        )

        # Reshape multi-dimensional outputs to 2-D arrays so each grid point is a row
        b_reshaped = b_list.reshape(b_list.shape[0], -1)
        max_idx_list_reshaped = max_idx_list.reshape(max_idx_list.shape[0], -1)

        # Save outputs. CSV files are convenient for sharing and plotting.
        np.savetxt(f'{var1}vs{var2}_{info}_exist_{nb_years}.csv', p_list, delimiter=",")
        np.savetxt(f'{var1}vs{var2}_{info}_exist1_{nb_years}.csv', p1_list, delimiter=",")
        np.savetxt(f'{var1}vs{var2}_{info}_exist2_{nb_years}.csv', p2_list, delimiter=",")
        np.savetxt(f'{var1}vs{var2}_{info}_prod_{nb_years}.csv', b_reshaped, delimiter=",")
        np.savetxt(f'{var1}vs{var2}_{info}_max_idx_{nb_years}.csv', max_idx_list_reshaped, delimiter=",")

        # 2) Full model with migration fixed for some components (migfix)
        info = '3C_F.6_F2.9_F1.6'
        p_list, p1_list, p2_list, b_list, max_idx_list = calculate_tot_production_3C_migfix(
                I1, e1, a1, amp1, F, F1, F2, h1, h2, rc, rp, mp, mc1, mc2,
                nb_years, X0, f, s, amp2=amp2, ext=threshold, B=B,
                var1=var1, var2=var2, range1=range1, range2=range2
        )
        b_reshaped = b_list.reshape(b_list.shape[0], -1)
        max_idx_list_reshaped = max_idx_list.reshape(max_idx_list.shape[0], -1)
        np.savetxt(f'{var1}vs{var2}_{info}_exist_{nb_years}.csv', p_list, delimiter=",")
        np.savetxt(f'{var1}vs{var2}_{info}_exist1_{nb_years}.csv', p1_list, delimiter=",")
        np.savetxt(f'{var1}vs{var2}_{info}_exist2_{nb_years}.csv', p2_list, delimiter=",")
        np.savetxt(f'{var1}vs{var2}_{info}_prod_{nb_years}.csv', b_reshaped, delimiter=",")
        np.savetxt(f'{var1}vs{var2}_{info}_max_idx_{nb_years}.csv', max_idx_list_reshaped, delimiter=",")

        # 3) A 2-consumer initial condition variant (modify X0), again with migfix
        #    Document in the paper that this run probes the model sensitivity to initial occupancy.
        X0 = [
                1, 1, 1,  # consumer 1 and base states present
                0,        # some resource / state absent initially
                1, 1, 0,  # consumer 2 partially absent
                1
        ]
        info = '2C_F1.6'
        p_list, p1_list, p2_list, b_list, max_idx_list = calculate_tot_production_3C_migfix(
                I1, e1, a1, amp1, F, F1, F2, h1, h2, rc, rp, mp, mc1, mc2,
                nb_years, X0, f, s, amp2=amp2, ext=threshold, B=B,
                var1=var1, var2=var2, range1=range1, range2=range2
        )
        b_reshaped = b_list.reshape(b_list.shape[0], -1)
        max_idx_list_reshaped = max_idx_list.reshape(max_idx_list.shape[0], -1)
        np.savetxt(f'{var1}vs{var2}_{info}_exist_{nb_years}.csv', p_list, delimiter=",")
        np.savetxt(f'{var1}vs{var2}_{info}_exist1_{nb_years}.csv', p1_list, delimiter=",")
        np.savetxt(f'{var1}vs{var2}_{info}_exist2_{nb_years}.csv', p2_list, delimiter=",")
        np.savetxt(f'{var1}vs{var2}_{info}_prod_{nb_years}.csv', b_reshaped, delimiter=",")
        np.savetxt(f'{var1}vs{var2}_{info}_max_idx_{nb_years}.csv', max_idx_list_reshaped, delimiter=",")

        # 4) Repeat with a different descriptor/parameterization for comparison:
        info = '3C_nomig_F.55_F2.9_F1.6'
        p_list, p1_list, p2_list, b_list, max_idx_list = calculate_tot_production_3C_nomig(
                I1, e1, a1, amp1, F, F1, F2, h1, h2, rc, rp, mp, mc1, mc2,
                nb_years, X0, f, s, amp2=amp2, ext=threshold, B=B,
                var1=var1, var2=var2, range1=range1, range2=range2
        )
        b_reshaped = b_list.reshape(b_list.shape[0], -1)
        max_idx_list_reshaped = max_idx_list.reshape(max_idx_list.shape[0], -1)
        np.savetxt(f'{var1}vs{var2}_{info}_exist_{nb_years}.csv', p_list, delimiter=",")
        np.savetxt(f'{var1}vs{var2}_{info}_exist1_{nb_years}.csv', p1_list, delimiter=",")
        np.savetxt(f'{var1}vs{var2}_{info}_exist2_{nb_years}.csv', p2_list, delimiter=",")
        np.savetxt(f'{var1}vs{var2}_{info}_prod_{nb_years}.csv', b_reshaped, delimiter=",")
        np.savetxt(f'{var1}vs{var2}_{info}_max_idx_{nb_years}.csv', max_idx_list_reshaped, delimiter=",")

        # 5) Full migration model
        info = '3C_F.55_F2.9_F1.6'
        p_list, p1_list, p2_list, b_list, max_idx_list = calculate_tot_production_3C(
                I1, e1, a1, amp1, F, F1, F2, h1, h2, rc, rp, mp, mc1, mc2,
                nb_years, X0, f, s, amp2=amp2, ext=threshold, B=B,
                var1=var1, var2=var2, range1=range1, range2=range2
        )
        b_reshaped = b_list.reshape(b_list.shape[0], -1)
        max_idx_list_reshaped = max_idx_list.reshape(max_idx_list.shape[0], -1)
        np.savetxt(f'{var1}vs{var2}_{info}_exist_{nb_years}.csv', p_list, delimiter=",")
        np.savetxt(f'{var1}vs{var2}_{info}_exist1_{nb_years}.csv', p1_list, delimiter=",")
        np.savetxt(f'{var1}vs{var2}_{info}_exist2_{nb_years}.csv', p2_list, delimiter=",")
        np.savetxt(f'{var1}vs{var2}_{info}_prod_{nb_years}.csv', b_reshaped, delimiter=",")
        np.savetxt(f'{var1}vs{var2}_{info}_max_idx_{nb_years}.csv', max_idx_list_reshaped, delimiter=",")

        # 6) Migration-fixed variant (again, for comparison)
        info = '3C_migfix_F.55_F2.9_F1.6'
        p_list, p1_list, p2_list, b_list, max_idx_list = calculate_tot_production_3C_migfix(
                I1, e1, a1, amp1, F, F1, F2, h1, h2, rc, rp, mp, mc1, mc2,
                nb_years, X0, f, s, amp2=amp2, ext=threshold, B=B,
                var1=var1, var2=var2, range1=range1, range2=range2
        )
        b_reshaped = b_list.reshape(b_list.shape[0], -1)
        max_idx_list_reshaped = max_idx_list.reshape(max_idx_list.shape[0], -1)
        np.savetxt(f'{var1}vs{var2}_{info}_exist_{nb_years}.csv', p_list, delimiter=",")
        np.savetxt(f'{var1}vs{var2}_{info}_exist1_{nb_years}.csv', p1_list, delimiter=",")
        np.savetxt(f'{var1}vs{var2}_{info}_exist2_{nb_years}.csv', p2_list, delimiter=",")
        np.savetxt(f'{var1}vs{var2}_{info}_prod_{nb_years}.csv', b_reshaped, delimiter=",")
        np.savetxt(f'{var1}vs{var2}_{info}_max_idx_{nb_years}.csv', max_idx_list_reshaped, delimiter=",")

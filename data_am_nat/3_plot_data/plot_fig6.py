#####
# Short description:
# Plot Fig.6 — load coexistence and production CSVs for a B x B parameter grid,
# compute summary statistics comparing runs with and without migration, and
# produce two line-plot PNGs that summarize those differences across one slice
# of the parameter grid.
#
# Inputs:
# - CSV files expected in `path` with names formatted as:
#   "{var1}vs{var2}_{info}_exist_{yr}.csv",
#   "{var1}vs{var2}_{info}_prod_{yr}.csv",
#   "{var1}vs{var2}_{info}_prod_m_{yr}.csv"
#   where {info} is an entry from `set_info` and {yr} is the simulation length.
# - Each production CSV stores results for 8 state variables flattened over the
#   B x B parameter grid; these flat arrays are reshaped to shape (8, B, B).
#
# Output:
# - Two PNG files written to `path`:
#   - 2024-12-18_metaeco_timing.png: percent change in total production (migration vs no migration)
#   - 2024-12-18_mig_timing.png: absolute change in consumer biomass (migration vs no migration)
#
# How to run:
# 1. Install dependencies: numpy, matplotlib
#    e.g. pip install numpy matplotlib
# 2. Adjust `path`, `set_info`, `a`, `yr`, `B`, and var1/var2 ranges near the top if needed.
# 3. Run the script:
#    python /Users/chaouala/zenodo_manuscript/data_am_nat/3_plot_data/plot_fig6.py
#
# Notes:
# - The script assumes a B x B parameter grid; B is defined below.
# - Files are loaded with numpy.loadtxt; NaNs in production arrays are set to 0
#   before reshaping and further computation.
# - The script accumulates total production across all state variables and also
#   computes consumer biomass as C1 + C2 for the comparisons.

import matplotlib.pyplot as plt #matplotlib.__version__ '3.8.2'
import numpy as np #numpy.__version__ '1.26.3'
from matplotlib import ticker
from matplotlib.colors import TwoSlopeNorm, Normalize, LinearSegmentedColormap
import os
### layout

#plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 34

### find path to data folder
# compute a path relative to this script file so the script works regardless of CWD
try:
    _script_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    _script_dir = os.getcwd()
path = os.path.normpath(os.path.join(_script_dir, '..', '2_data_csv', 'fig6_csv'))

# grid size (B x B)
B=50
# simulation duration (used in filenames)
yr = 10


# variables defining the two axes and their numeric ranges:
# var1, var2 are the variable labels used in filenames; range1, range2 are the
# numeric tick values (percent/scales) for the plotted axes.
var1, var2, range1, range2=['F','Fmig',  (np.linspace(.3,.9,B)/.6-1)*100, np.linspace(0.3,.7,B)*100] 

# list of simulation parameter set identifiers (file suffixes)
set_info = ['1C',
    '2C',
        '3C_F1.6_F2.6'


]


# index along the first axis at which to take a horizontal slice for plotting
a= 24

for info in set_info: 
    print(info)
    ### Load Data ###

    # load existence grid (used to infer grid shape)
    p_list =np.loadtxt(f"{path}/{var1}vs{var2}_{info}_exist_{yr}.csv", delimiter=',')
    # load production arrays: no migration and with migration
    b_list_2D =np.loadtxt(f"{path}/{var1}vs{var2}_{info}_prod_{yr}.csv", delimiter=",")
    b_list_2D[np.isnan(b_list_2D)] = 0  # replace NaNs with zeros
    b_list_2D_m =np.loadtxt(f"{path}/{var1}vs{var2}_{info}_prod_m_{yr}.csv", delimiter=",")
    b_list_2D_m[np.isnan(b_list_2D_m)] = 0
    # reshape flattened production arrays into (8 state variables, B, B)
    b_list = b_list_2D.reshape((8, p_list.shape[0], p_list.shape[1]))
    b_list_m = b_list_2D_m.reshape((8, p_list.shape[0], p_list.shape[1]))


    ####plot###

    ## meta-ecological labels and axis descriptions ###
    names = ['N1', 'N2', 'P1', 'P2','C1', 'C2']
    parameters = ['amp1','mc1','Fmig','F','F1','a1','F2', 'per_s1', 'rc']
    axis = ['Amplitude of seasonality', 'Mortality rate of conumer 1', 'Time spent in the seasonal \n ecosystem($\%$)',
        'Competitive difference ($\%$) \n ($F_{mig}-F_1$)', 'Attack rate of the consumer1 ',  'Uptake rate of the primary producer',
        'F2-F', 'Percentage of summer1', 'Recycling rate of the consumer']


    x = range1 
    y = range2 
    step_x = abs(x[0]-x[-1])/(2*B)
    step_y = abs(y[0]-y[-1])/(2*B)
    ### Sum Total Prod, Ratio ###
    # sum production across all 8 state variables to get total production per grid cell
    data = np.sum(b_list[:, :, :], axis=0) # no migration
    data_m = np.sum(b_list_m[:, :, :], axis=0)  # with migration
    # compute percent change (migration vs no migration) and store by info label
    if info == '1C':
        dp1C = (data_m - data) / data * 100
    elif info.startswith('2C'):
        dp2C = (data_m - data) / data * 100
    else:
        dp3C = (data_m - data) / data * 100

   
    #     ### Sum Total consumer biomass ###
    # consumer biomass is C1 + C2 (indices 4 and 5 in the 8-state production array)
    C_list = b_list[5,:,:]+b_list[4,:,:]

    C_list_m = b_list_m[5,:,:]+b_list_m[4,:,:]


    datac =  C_list # consumer biomass without migration
    datac_m = C_list_m # consumer biomass with migration

    # store absolute difference (with - without) by info label
    if info == '1C':
        dpc_1C = datac_m - datac
    elif info.startswith('2C'):
        dpc_2C = datac_m - datac
    else:
        dpc_3C = datac_m - datac
# =======================
# Percent production difference plot (line slice across one grid row)
# =======================   



fig, axs = plt.subplots(1, 1, figsize=(8, 11))

# === Plot lines for percent change in total production ===
axs.plot(y, dp1C[a, :], color='#C36346', label='1 consumer')
axs.plot(y, dp2C[a, :], color='#4DAF4A', label='2 consumers')  # new green line
axs.plot(y, dp3C[a, :], color='#A997DF', label='3 consumers')

# === Mark maxima for each line ===
max_y1 = np.max(dp1C[a, :]); max_x1 = np.round(y[np.argmax(dp1C[a, :])])
max_y2 = np.max(dp2C[a, :]); max_x2 = np.round(y[np.argmax(dp2C[a, :])])
max_y3 = np.max(dp3C[a, :]); max_x3 = np.round(y[np.argmax(dp3C[a, :])])

axs.plot(max_x1, max_y1, 'o', color='#C36346', markersize=6)
axs.plot(max_x2, max_y2, 'o', color='#4DAF4A', markersize=6)
axs.plot(max_x3, max_y3, 'o', color='#A997DF', markersize=6)

# === Add labels for maxima ===
axs.text(max_x1, 31, f"{max_x1:.1f}", color='#C36346', ha='center')
axs.text(max_x2, 31, f"{max_x2:.1f}", color='#4DAF4A', ha='center')
axs.text(max_x3, 31, f"{max_x3:.1f}", color='#636363', ha='center')

# === Add vertical lines at maxima ===
axs.axvline(max_x1, color='#C36346', linestyle='dotted', linewidth=1.5)
axs.axvline(max_x2, color='#4DAF4A', linestyle='dotted', linewidth=1.5)
axs.axvline(max_x3, color='#636363', linestyle='dotted', linewidth=1.5)

# === Extra descriptive text ===
axs.text(40, 29, "Late migration", fontsize=22, color='#fb9d02ff', ha='center')
axs.text(60, 29, "Early migration", fontsize=22, color='#61bafaff', ha='center')

# === Axis formatting and labels ===
j = parameters.index(var2)
y_name = axis[j]
axs.margins(x=0)
axs.set_ylabel('Difference in biomass stock \nwith and without migration ($\%$)')
axs.set_xlabel(f'{y_name}')
axs.legend(loc='center left')
plt.xticks([30, 50, 70])
plt.tight_layout(pad=0)
plt.show()
fig.savefig(f'{path}/2024-12-18_metaeco_timing.png')



# =======================
# Absolute consumer biomass difference plot (line slice)
# =======================
fig, axs = plt.subplots(1, 1, figsize=(8, 11))

# === Plot lines for absolute difference in consumer biomass ===
axs.plot(y, dpc_1C[a, :], color='#C36346', label='1 consumer')
axs.plot(y, dpc_2C[a, :], color='#4DAF4A', label='2 consumers')
axs.plot(y, dpc_3C[a, :], color='#A997DF', label='3 consumers')

# === Mark maxima ===
max_y1 = np.max(dpc_1C[a, :]); max_x1 = np.round(y[np.argmax(dpc_1C[a, :])])
max_y2 = np.max(dpc_2C[a, :]); max_x2 = np.round(y[np.argmax(dpc_2C[a, :])])
max_y3 = np.max(dpc_3C[a, :]); max_x3 = np.round(y[np.argmax(dpc_3C[a, :])])

axs.plot(max_x1, max_y1, 'o', color='#C36346', markersize=6)
axs.plot(max_x2, max_y2, 'o', color='#4DAF4A', markersize=6)
axs.plot(max_x3, max_y3, 'o', color='#A997DF', markersize=6)

# === Add labels near maxima ===
axs.text(max_x1-4.5, 2.9, f"{max_x1:.1f}", color='#C36346', ha='center')
axs.text(max_x2+2.5, 2.9, f"{max_x2:.1f}", color='#4DAF4A', ha='center')
axs.text(max_x3+2, 2.9, f"{max_x3:.1f}", color='#636363', ha='center')

# === Vertical lines at maxima ===
axs.axvline(max_x1, color='#C36346', linestyle='dotted', linewidth=1.5)
axs.axvline(max_x2, color='#4DAF4A', linestyle='dotted', linewidth=1.5)
axs.axvline(max_x3, color='#636363', linestyle='dotted', linewidth=1.5)

# === Extra descriptive text ===
axs.text(40, 2.75, "Late migration", fontsize=22, color='#fb9d02ff', ha='center')
axs.text(60, 2.74, "Early migration", fontsize=22, color='#61bafaff', ha='center')

# === Axis formatting and labels ===
j = parameters.index(var2)
y_name = axis[j]
axs.set_ylabel('Difference in biomass stock \n with and without migration (a.u.)')
axs.set_xlabel(f'{y_name}')
axs.margins(x=0)
plt.xticks([30, 50, 70])
plt.tight_layout(pad=0)
plt.show()
fig.savefig(f'{path}/2024-12-18_mig_timing.png')


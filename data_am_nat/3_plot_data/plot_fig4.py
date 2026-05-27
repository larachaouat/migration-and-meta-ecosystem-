#####


# Short description:
# Plot Fig.4 — reads coexistence result CSVs for different parameter grids,
# creates a colored coexistence map and saves a PNG.
#
# Inputs:
# - CSV files expected in `path` (set below) with names like:
#   "{var1}vs{var2}_{info[a]}_exist_{yr}.csv",
#   "{var1}vs{var2}_{info[a]}_exist1_{yr}.csv",
#   "{var1}vs{var2}_{info[a]}_exist2_{yr}.csv"
# - The `info` list and index `a` select which simulation folder / parameter set to use.
#
# Output:
# - A PNG saved to the same `path`
#
# How to run (short & clear):
# 1. Make sure Python 3 and packages are installed: numpy, matplotlib.
#    e.g. pip install numpy matplotlib
# 2. Edit `path`, `a`, `yr`, or other parameters near the top if needed.
# 3. Run the script:
#    python /Users/chaouala/zenodo_manuscript/data_am_nat/3_plot_data/plot_fig3.py
#
# Notes:
# - This script expects a B x B grid of CSV results (B set below).
# - It overlays masks and places textual labels for different coexistence outcomes.

### libraries 
import matplotlib
import matplotlib.pyplot as plt #matplotlib.__version__ '3.8.2'
import numpy as np #numpy.__version__ '1.26.3'
from matplotlib.colors import ListedColormap 
import matplotlib.patches as mpatches
import os

### layout
matplotlib.use('TkAgg')  # Or 'Qt5Agg', 'GTK3Agg', etc.

#plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 32 

### find path file
# locate data folder relative to this script (robust if __file__ is unavailable)
try:
    _script_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    _script_dir = os.getcwd()
path = os.path.normpath(os.path.join(_script_dir, '..', '2_data_csv', 'fig4_csv'))

# size of simulation = B*B
B=100
# number of iteration per simulations yr*1000
yr = 30

# carecteristics of csv
info = [
    '2C_F.55_F1.6',
    '2C_F.55_F1.6_migfix',
    '3C_F.55_F2.9_F1.6',
    '3C_F.55_F2.9_F1.6_migfix'
    ]

a=2

# variables of interest: define parameter pairs and ranges to load the correct CSVs
set = [ 
    ['amp1', 'per_s1' ,np.linspace(0.6,1,B)*100, np.linspace(0.1,.9,B)*100],
]

for var1, var2, range1, range2 in set: 
    # Load precomputed coexistence matrices (main & two variants)
    p_list = np.loadtxt(f"{path}/{var1}vs{var2}_{info[a]}_exist_{yr}.csv", delimiter=',')
    p1_list = np.loadtxt(f"{path}/{var1}vs{var2}_{info[a]}_exist1_{yr}.csv", delimiter=',')
    p2_list = np.loadtxt(f"{path}/{var1}vs{var2}_{info[a]}_exist2_{yr}.csv", delimiter=',')

    # labels and axis names (used for plot annotation)
    names = ['N1', 'N2', 'P1', 'P2','C1', 'C2']
    parameters = ['amp1','Fmig',
          'F','F1',
          'F2', 'per_s1']
    axis = ['Amplitude of seasonality ($\%$)',  'Percentage of time in ecosystem1',
        'Competitive difference ($\%$) \n($F_{mig}-F_1$)', ' Attack rate (C1) ', 
        'Attack rate (C2)', 'Length of summer over a year ($\%$)']

    i= parameters.index(var1)
    j= parameters.index(var2)
    x_name = axis[i]
    y_name=axis[j]

    ####  build coexistence classification map
    conditions = [
    (np.isin(p_list, [0, 1]) & np.isin(p1_list, [2, 3]) &  np.isin(p2_list, [0, 1])),
    (np.isin(p_list, [0, 1]) & np.isin(p1_list, [0, 1]) &  np.isin(p2_list, [2, 3])), 
    (np.isin(p_list, [2, 3]) & np.isin(p1_list, [0, 1]) &  np.isin(p2_list, [0, 1])),
    (np.isin(p_list, [2, 3]) & np.isin(p1_list, [2, 3]) &  np.isin(p2_list, [0, 1])),
    (np.isin(p_list, [2, 3]) & np.isin(p1_list, [0, 1]) &  np.isin(p2_list, [2, 3])),
    (np.isin(p_list, [0, 1]) & np.isin(p1_list, [2, 3]) &  np.isin(p2_list, [2, 3])),
    (np.isin(p_list, [2, 3]) & np.isin(p1_list, [2, 3]) &  np.isin(p2_list, [2, 3])),
    ]
    outputs = [2, 3, 1, 12, 13, 23, 123]
    coex_tot_m = np.select(conditions, outputs, default=0)

    # Labels used for annotations/legend
    legend_labels = {
    0: '',
    1: '$C_{mig}$', 
    2: '$C_1$', 
    3: '$C_2$', 
    12: '$C_{mig}$ \n+ $C_1$', 
    13: '$C_{mig}$\n + $C_2$', 
    23: '$C_1$ \n+ $C_2$', 
    123: '$C_{mig}$\n + $C_1$ \n+ $C_2$', 
    }

    # Pastel colors and patches (visual choices for the figure)
    legend_patches = [
    mpatches.Patch(color='#F5F5F5', label=legend_labels[0]),  # Pastel whitesmoke
    mpatches.Patch(color='#FFB3B3', label=legend_labels[1]),  # Pastel red
    mpatches.Patch(color='#ADD8E6', label=legend_labels[2]),  # Pastel blue
    mpatches.Patch(color='#FFFACD', label=legend_labels[3]),  # Pastel yellow
    mpatches.Patch(color='#D8BFD8', label=legend_labels[12]),  # Pastel purple
    mpatches.Patch(color='#FFDAB9', label=legend_labels[13]),  # Pastel orange
    mpatches.Patch(color='#98FB98', label=legend_labels[23]),  # Pastel green
    mpatches.Patch(color='lightgrey', label=legend_labels[123])  # Pastel pink
    ]

    colorsList = ['#F5F5F5', '#FFB3B3', '#ADD8E6', '#FFFACD', '#D8BFD8', '#FFDAB9', '#98FB98', '#FF69B4']
    cmap = ListedColormap(colorsList)

    # mapping from coexistence value -> color index
    value_to_color = {
    0: 0,  
    1: 1,  
    2: 2,  
    3: 3, 
    12: 4, 
    13: 5, 
    23: 6, 
    123: 7 
    }
    color_indices = np.vectorize(value_to_color.get)(coex_tot_m)

    # plotting grid ranges and half-step for extent
    x = range1 
    y = range2  
    step_x = abs(x[0]-x[-1])/(2*B)
    step_y = abs(y[0]-y[-1])/(2*B)

    fig, ax = plt.subplots( figsize=(8,11 ))

    # Build RGBA masks to overlay presence/absence of species types (translucent)
    red_mask = np.zeros((*p_list.shape, 4), dtype=float)    # 4 channels for RGBA
    blue_mask = np.zeros((*p1_list.shape, 4), dtype=float)
    green_mask = np.zeros((*p2_list.shape, 4), dtype=float)
    alpha = 0.3

    # red mask from p_list (non 0/1 entries)
    red_mask[..., 0] =  ((p_list != 0) & (p_list != 1))*1                                    
    red_mask[..., 3] = ((p_list != 0) & (p_list != 1)) * alpha  

    # blue mask from p1_list
    blue_mask[..., 2] = ((p1_list != 0) & (p1_list != 1)) * 1                                   
    blue_mask[..., 3] = ((p1_list != 0) & (p1_list != 1)) * alpha               

    # green_mask is defined but not set (left blank in original code)
    # green_mask[..., 1] etc. can be enabled if p2_list should be visualized similarly.

    # overlay masks (imshow uses extent to map array indices to parameter values)
    ax.imshow(np.rot90(red_mask),  extent=[x.min()-step_x, x.max()+step_x, y.min()-step_y, y.max()+step_y], aspect='auto')
    ax.imshow(np.rot90(green_mask),  extent=[x.min()-step_x, x.max()+step_x, y.min()-step_y, y.max()+step_y], aspect='auto')
    ax.imshow(np.rot90(blue_mask),  extent=[x.min()-step_x, x.max()+step_x, y.min()-step_y, y.max()+step_y], aspect='auto')

    # Place text labels for each unique coexistence class by averaging coordinates
    unique_indices = np.unique(color_indices)
    for idx in unique_indices:
        y_pos, x_pos = np.where(np.flipud(np.rot90(color_indices))== idx)
        if len(y_pos) > 0 and len(x_pos) > 0:
            avg_x = x[x_pos].mean()-step_x*20 
            avg_y = y[y_pos].mean()+step_y*3
            ax.text(avg_x, avg_y, list(legend_labels.values())[idx], fontsize=32, color='black')

    # Axis labels and layout tweaks
    ax.set_xlabel(f'{x_name}')
    ax.set_ylabel(f'{y_name}')
    plt.locator_params(axis='y', nbins=5)
    plt.locator_params(axis='x', nbins=3)
    plt.tight_layout(pad=0)

    # Save figure
    fig.savefig(f'{path}/2025-01-07_ampvspers_coexistence_{info[a]}.png')

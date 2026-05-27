#####


# Short description:
# Plot a coexistence map — reads coexistence result CSV for specified parameter grids,
# creates a colored coexistence map with text labels for different outcomes and saves a PNG.
#
# Inputs:
# - CSV files expected in `path` (set below) with names like:
#   "{var1}vs{var2}_{info[1]}_exist_{yr}.csv"
#   (this script currently loads the second entry of `info` via info[1])
# - The `set` list defines pairs of variables and their numeric ranges to plot.
#
# Output:
# - A PNG saved to the same `path` (filename: persistence_mismatched.png)
#
# How to run (short & clear):
# 1. Make sure Python 3 and packages are installed: numpy, matplotlib.
#    e.g. pip install numpy matplotlib
# 2. Edit `path`, `yr`, `B`, or entries in `set` / `info` near the top if needed.
# 3. Run the script:
#    python /Users/chaouala/zenodo_manuscript/3_plot_data/plot_fig2.py
#
# Notes:
# - This script expects a B x B grid of CSV results (B set below).
# - It loops over variable pairs in `set`, loads the matching CSV (uses info[1]),
#   and overlays text labels at the average positions of each outcome class.


### libraries 
import matplotlib.pyplot as plt #matplotlib.__version__ '3.8.2'
import numpy as np #numpy.__version__ '1.26.3'
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
import os

### layout
#plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 34

### find path file
# locate data folder relative to this script (robust if __file__ is unavailable)
try:
    _script_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    _script_dir = os.getcwd()
path = os.path.normpath(os.path.join(_script_dir, '..', '2_data_csv', 'fig2_csv'))


# size of simulation = B*B
B=100
# number of iteration per simulations yr*1000
yr = 30

# carecteristics of csv

info = [ '1C' ,
    '1C_migfix' 

]

a=0
# variables of interset

set = [ ['amp1', 'per_s1' ,np.linspace(0,1,B)*100, np.linspace(0.1,.9,B)*100],
    ['amp1', 'F', np.linspace(0,1,B), np.linspace(.2,.4,B)-.3]
       ]

for var1, var2, range1, range2 in set: 
    ### Load Data ###
    p_list = np.loadtxt(f'{path}/{var1}vs{var2}_{info[a]}_exist_{yr}.csv', delimiter=',') #_{info[1]}

    ### Plot ###
    # variables 
    names = ['N1', 'N2', 'P1', 'P2','C1', 'C2']
    parameters = ['amp1','Fmig',
          'F','F1',
          'F2', 'per_s1']
    axis = ['Amplitude of seasonality ($\%$)',  'Percentage of time in ecosystem1',
        'F1-F', ' Attack rate (C1) ', 
        'Attack rate (C2)', 'Length of summer over a year ($\%$)']
    
    # Define labels for the legends
    legend_labels = {
    0: 'Extinct',
    1: 'Persistence in \n a seasonal meta-ecosystem', 
    2: 'Persistence \n with \n migration', 
    3: 'Persistence \nwith and without \n migration'
    }

    # Create custom colored patches for the legend
    legend_patches = [
    mpatches.Patch(color='whitesmoke', label=legend_labels[0]),
    mpatches.Patch(color='orange', label=legend_labels[1]),
    mpatches.Patch(color='purple', label=legend_labels[2]),
    mpatches.Patch(color='lightblue', label=legend_labels[3])
    ]

    fig, ax = plt.subplots( figsize=(8, 11))

    i= parameters.index(var1)
    j= parameters.index(var2)
    x_name = axis[i]
    y_name=axis[j]

    ### Create persistence map ###
    x = range1
    y = range2
    centers = [x[0],x[-1],y[0],y[-1]]
    dx, = np.diff(centers[:2])/(B)
    dy, = -np.diff(centers[2:])/(B)

    extent = [centers[0]-dx/2, centers[1]+dx/2, centers[2]+dy/2, centers[3]-dy/2]

    colorsList = ['whitesmoke', 'orange', 'plum', 'lightblue']
    cmap = ListedColormap(colorsList)
    d = p_list 

    cax = ax.imshow(np.flipud(np.rot90(d)),  cmap=cmap, origin='lower', 
                 extent=extent, aspect='auto', vmin=0, vmax=3 )
###uncomment for graph latex
    ax.set_xlabel(f'{x_name}')
    ax.set_ylabel(f'{y_name}')


    # Add text labels directly on the plot based on color_indices
    unique_indices = np.unique(d)

    for idx in unique_indices:
        y_pos, x_pos = np.where(np.flipud(np.rot90(d))== idx)

    # Check if positions are valid
    if len(y_pos) > 0 and len(x_pos) > 0:
        # Calculate the average position for placing the label
        avg_x = x[x_pos].mean()
        avg_y = y[y_pos].mean()+2*dy
        print(avg_y)
        ax.text(avg_x, avg_y, list(legend_labels.values())[int(idx)], color='black',fontsize=34, ha='center', va='center')


    plt.locator_params(axis='y', nbins=6)
    plt.tight_layout(pad=0)
    plt.show()
    fig.savefig(f'{path}/persistence_mismatched.png')




#####


# Short description:
# Plot Fig.5 — read coexistence / production CSVs for a grid of parameter values,
# compute summary statistics (total production change and consumer biomass change
# when migration is present vs absent), and save two heatmap PNGs.
#
# Inputs:
# - CSV files expected in `path` (set below) with names formatted as:
#   "{var1}vs{var2}_{info}_exist_{yr}.csv",
#   "{var1}vs{var2}_{info}_prod_{yr}.csv",
#   "{var1}vs{var2}_{info}_prod_m_{yr}.csv"
#   where {info} is selected from `set_info` and {yr} is the simulation length.
# - Each production CSV contains an 2D grid of results stacked for 8 state variables,
#   so the flat production array is reshaped into shape (8, B, B).
#
# Output:
# - Two PNG files saved to `path`:
#   - heatmap of percent change in total production with migration
#   - heatmap of absolute change in consumer biomass with migration
#
# How to run:
# 1. Install dependencies: numpy, matplotlib
#    e.g. pip install numpy matplotlib
# 2. Adjust `path`, `set_info`, `a`, `yr`, `B`, and var1/var2 ranges near the top if needed.
# 3. Run the script:
#    python /Users/chaouala/zenodo_manuscript/data_am_nat/3_plot_data/plot_fig5.py
#
# Notes:
# - The script assumes a B x B parameter grid (B given below).
# - Files are loaded using numpy.loadtxt and NaNs in production arrays are set to 0.
# - The script creates two figures per `info`: percent difference in total production
#   and absolute difference in consumer biomass (C1+C2) between migration and no migration.

import matplotlib.pyplot as plt #matplotlib.__version__ '3.8.2'
import numpy as np #numpy.__version__ '1.26.3'
from matplotlib import ticker
from matplotlib.colors import TwoSlopeNorm, Normalize, LinearSegmentedColormap
import os

### layout
# global plotting settings (font size)
#plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 32

### find path to data folder
# compute path relative to this script directory to be robust to different CWDs
try:
        _script_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
        _script_dir = os.getcwd()
path = os.path.normpath(os.path.join(_script_dir, '..', '2_data_csv', 'fig5_csv'))

# size of parameter grid: B x B
B=50
# number of time units in each simulation (used in filename)
yr = 30


# Define which two parameters are being compared (var1 vs var2) and their axis ranges.
# `range1` and `range2` are 1D arrays of length B describing the grid values for each axis.
var1, var2, range1, range2=['F','Fmig',  (np.linspace(.3,.9,B)/.6-1)*100, np.linspace(0.3,.7,B)*100] 



# List of simulation folders / parameter set identifiers to iterate over
set_info = [#'1C',
        #'2C_F1.6',
                '3C_F1.6_F2.6'


]


a= 0

for info in set_info: 
        print(info)
         ### Load Data ###
        # p_list: coexistence or existence map (not used below but loaded for consistency)
        p_list =np.loadtxt(f"{path}/{var1}vs{var2}_{info}_exist_{yr}.csv", delimiter=',')
        # b_list_2D: stacked production (no migration) as a 2D array that will be reshaped
        b_list_2D =np.loadtxt(f"{path}/{var1}vs{var2}_{info}_prod_{yr}.csv", delimiter=",")
        # replace NaNs (missing data) with 0 for production arrays
        b_list_2D[np.isnan(b_list_2D)] = 0
        # b_list_2D_m: stacked production with migration
        b_list_2D_m =np.loadtxt(f"{path}/{var1}vs{var2}_{info}_prod_m_{yr}.csv", delimiter=",")
        b_list_2D_m[np.isnan(b_list_2D_m)] = 0
        # reshape flat 2D production array into (8 variables, B, B) so we can index variables
        b_list = b_list_2D.reshape((8, p_list.shape[0], p_list.shape[1]))
        b_list_m = b_list_2D_m.reshape((8, p_list.shape[0], p_list.shape[1]))


        #### plot total production change (percent) ###

        # labels and parameter descriptions used for axis titles
        names = ['N1', 'N2', 'P1', 'P2','C1', 'C2']
        parameters = ['amp1','mc1','Fmig','F','F1','a1','F2', 'per_s1', 'rc']
        axis = ['Amplitude of seasonality', 'Mortality rate of conumer 1', 'Time spent in the seasonal ecosystem ($\%$)',
                        'Competitive difference ($\%$) \n ($F_{mig}-F_i$)', 'Attack rate of the consumer1 ',  'Uptake rate of the primary producer',
                        'F2-F', 'Percentage of summer1', 'Recycling rate of the consumer']

        # Choose a colormap for the heatmaps
        end_color = "#712f0aff"  # dark purple (not directly used below)
        start_color = "#f9a741ff"      # orange (not directly used below)
        custom_cmap = plt.colormaps['PuOr'].reversed()

        # x and y coordinates for the grid and half-step for image extent
        x = range1 
        y = range2 
        step_x = abs(x[0]-x[-1])/(2*B)
        step_y = abs(y[0]-y[-1])/(2*B)

        # Sum production across all 8 variables to get total production per grid cell
        data = np.sum(b_list[:, :, :], axis=0)       # total production without migration
        data_m = np.sum(b_list_m[:, :, :], axis=0)   # total production with migration
        # percent change in total production when migration is present
        dp3C = (data_m-data)/data *100

        # Plot percent change heatmap
        max = 12
        cmap =  custom_cmap
        fig, axs = plt.subplots(1, 1,figsize=(8, 11))    
        # Use a diverging normalization centered at zero to show increases vs decreases
        norm = TwoSlopeNorm(vmin=-max, vcenter=0, vmax=max)
        cmap.set_under('w')
        # flip/rotate matrix for correct orientation on axes and set extent so tick labels match parameter ranges
        cax = axs.imshow(np.flipud(np.rot90(dp3C)), aspect='auto', origin='lower', cmap=cmap, norm=norm,
                                                        extent=[x.min()-step_x, x.max()+step_x, y.min()-step_y, y.max()+step_y])
        fig.colorbar(cax, ax=axs, label='Biomass difference with and without migration(%)')
        # set axis labels using parameter descriptions
        i= parameters.index(var1)
        j= parameters.index(var2)
        x_name = axis[i]
        y_name=axis[j]

        axs.set_xlabel(f'{x_name}')
        axs.set_ylabel(f'{y_name}')
        plt.locator_params(axis='y', nbins=5)
        plt.tight_layout(pad=0) 
        # save percent-change heatmap
        fig.savefig(f'{path}/2024-12-18__heatmap_biom_meta.png')



        ### Sum Total consumer biomass change (absolute) ###
        # Consumers are assumed to be in indices 4 and 5 (C1, C2) in the b_list arrays
        C_list = b_list[5,:,:]+b_list[4,:,:]        # total consumer biomass without migration
        C_list_m = b_list_m[5,:,:]+b_list_m[4,:,:]  # total consumer biomass with migration

        datac =  C_list
        datac_m = C_list_m
        # absolute difference (with migration - without migration)
        dpc_3C = (datac_m-datac)
        
        ### Create the consumer biomass heatmap ###
        fig, axs = plt.subplots(1, 1,figsize=(8, 11)) 
        cmap.set_under('w')
        # center normalization around zero (symmetric vmin/vmax)
        norm = TwoSlopeNorm(vmin=-dpc_3C.max(), vcenter=0, vmax=dpc_3C.max())
        cax = axs.imshow(np.flipud(np.rot90(dpc_3C)), aspect='auto', origin='lower', cmap=custom_cmap, norm=norm,
                                                   extent=[range1.min(), range1.max(), range2.min(), range2.max()])
        cb= fig.colorbar(cax, ax=axs, label='Biomass difference with and without migration(a.u.)')
        # axis labels
        i= parameters.index(var1)
        j= parameters.index(var2)
        x_name = axis[i]
        y_name=axis[j]

        axs.set_xlabel(f'{x_name}')
        axs.set_ylabel(f'{y_name}')
        # set colorbar tick locator and update ticks for readability
        tick_locator = ticker.MaxNLocator(nbins=5)
        cb.locator = tick_locator
        cb.update_ticks()   
        plt.locator_params(axis='y', nbins=5)
        plt.tight_layout(pad=0)
        plt.show()
        # save consumer biomass heatmap
        #fig.savefig(f'{path}/2024-12-18__heatmap_biom_mig.png')


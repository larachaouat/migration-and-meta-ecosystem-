

GENERAL INFORMATION
This README.txt file was updated on 2026/06/15S
 
A. Paper associated with this archive
Citation: 
X

Brief abstract: Migration is a ubiquitous process that links ecosystems with distinct seasonal dynamics, transferring biomass and species interactions across space. Despite being widely altered by global change, studies commonly overlook the interaction of seasonal characteristics and bidirectional migration on species coexistence and biomass stocks across meta-ecosystems. We developed a mathematical model to study how migration interacts with key characteristics of seasonality—amplitude of seasonal variation and length of summer—to influence migrant persistence, coexistence of migratory and non-migratory species, and biomass stocks. Our findings demonstrate that seasonal characteristics mediate the effect of migration on coexistence and biomass stocks across meta-ecosystems. However, the effects strongly depend on migration timing: phenological mismatches can reduce biomass at local and meta-ecosystem scales and lead to the extinction of migratory and non-migratory consumers. Our study highlights how migration and seasonality interact to shape community structure and ecosystem function across scales, emphasizing the importance of system-level approaches for studying ecological outcomes of global change.

B. Originators
Lara Chaouat (1,2), Florian Altermatt (1,2), Tianna Peller (1,2,3)
 1. Department of Aquatic Ecology, Eawag: Swiss Federal Institute of Aquatic Science and Technology, Ueberlandstrasse 133, Dübendorf 8600, Switzerland. 

 2. Department of Evolutionary Biology and Environmental Studies, University of Zürich, Winterthurerstrasse 190, Zürich 8006, Switzerland. 

 3. Department of Ecology and Evolutionary Biology, University of Toronto, 25 Wilcocks St, Toronto, Ontario M5S 1A1, Canada.
    
C. Contact information
lara.chaouat@uzh.ch

D. Dates of data collection
January 2025

F. Funding Sources
Funding is from the Swiss National Science Foundation (grant 310030 197410)

ACCESS INFORMATION
1. Licenses/restrictions placed on the data or code: CC0 1.0 Universal (CC0 1.0) Public Domain Dedication
2. Data derived from other sources: No
3. Recommended citation for this data/code archive: X

DATA & CODE FILE OVERVIEW

Of 2_data_csv

This data repository consist of 62 data files, 14 code scripts, and this README document, with the following data and code filenames and variables
Data files and variables
In the data folder there is a folder per figure with the csv files corresponding to the figure (one can generate these csv files with the pipeline in folder 1 explained in the Readme in the folder).
Each csv file contains a biomass/existence/dominance matrix (2D: B x B).
File naming convention for csv files:
	<var1>vs<var2>_<info>_exist_<nb_years>.csv     -> overall occupancy matrix
	<var1>vs<var2>_<info>_exist1_<nb_years>.csv    -> consumer 1 occupancy
	<var1>vs<var2>_<info>_exist2_<nb_years>.csv    -> consumer 2 occupancy
	<var1>vs<var2>_<info>_prod_<nb_years>.csv      -> production/biomass
	<var1>vs<var2>_<info>_max_idx_<nb_years>.csv   -> dominant-state index

 Notes on reproducibility / important parameters:

   - B         : grid resolution (number of points per varying parameter)
   - nb_years  : length (time) of each simulation run used by the model functions
   - info	: important parameters of the model (number of consumers, migration...)
   - var1 and var2 are the x and y axis of the plots


Code scripts and workflow

Of the 1_pipeline_producing_data

	1. Class meta_eco
  		- Defines the meta-ecosystem object and includes the ODE solving.
   		- Output: one simulation for one parameter set.

	2. Function for simulation
   		- Runs the ODE solver for a defined period of time (number of years).
  		- Outputs CSV files saved in the data folder.
   		- Output: CSV files calculating persistence, biomass average over a year, maximum indexes, etc.
   		- File naming convention: var1_vs_var2_type_data_info_meta-eco_nb-years.csv

	3. Execution script
   		- Executes the function for different variables and meta-ecosystem structures.
   		- Output: a folder in data named by date (e.g., data/date_csv).

Of the 3_plot_data

	Plot Fig.2 — reads coexistence result CSV for specified parameter grids, creates a colored coexistence map and saves a PNG.
	Plot Fig.3 — reads coexistence result CSVs for different parameter grids, creates a colored coexistence map and saves a PNG.
	Plot Fig.4 — reads coexistence result CSV for specified parameter grids, creates a colored coexistence map and saves a PNG.
	Plot Fig.5 — reads coexistence result CSVs for different parameter grids, creates a heatmap of total production change and consumer biomass change and saves a PNG.
	Plot Fig.6 — reads coexistence result CSVs for different parameter grids, creates a line-plot that summarize those differences across one slice of the parameter grid and saves a PNG.



SOFTWARE VERSIONS
Python 3 and packages are installed: numpy (numpy.__version__ '1.26.3'), matplotlib (matplotlib.__version__ '3.8.2') and scipy (scipy.__version__ '1.13.0')






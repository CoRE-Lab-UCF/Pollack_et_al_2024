# Pollack_et_al_2024

# The statistical framework for generating boundary conditions for compound flood models

This directory contains MATLAB scripts used for estimating the return periods of compound flood drivers in coastal regions, focusing on the mixed populations (tropical cyclones (TCs) and non-tropical cyclones (non-TCs)). The methodology and results are detailed in the paper "[Funding rules that promote equity in climate adaptation outcomes](https://osf.io/preprints/osf/6ewmu)."

<br>

>Pollack, Adam & Santamaria-Aguilar, Sara & Maduwantha, Pravin & Keller, Klaus. (2024). Funding rules that promote equity in climate adaptation outcomes. 10.31219/osf.io/6ewmu.

<be>

## Overview of Scripts
**The input and output files are generated by the MATLAB scripts (except the raw data inputs). Under the "Input" files, the file name is given first following the format of the file**

### 1. Bias Correction
- **`Bias_Correction.m`**

**Description**: This script applies a bias correction to the rainfall data using quantile mapping. The method adjusts the hourly rainfall gauge data to match the basin-average rainfall values derived from the Analysis of Record for Calibration (AORC) data.-

**Input**:
- `Pr_Measured.mat`:- [date-time, hourly precipitation]
- `Pr_AORC.mat`:- [date-time, hourly precipitation]

**output**: 
- `Bias Corrected MS at {location}_With_AORC_of_GC.mat`:- [date-time, hourly precepitation]

 
### 2. Calculating accumulated rainfall 
- **`Creating_Accumulated_RF.m`**

**Description**: The script calculates the accumulated rainfall from 1 to 48 hours suing the hourly rainfall data.
    
**Input**: 
- `Bias Corrected MS at “location”_With_AORC_of_GC.mat`:- [date-time, hourly precipitation]

**output**:  
- `Hourly_accumulation_Bias_Corrected_RF_data_”location”.mat`:- [date-time, hourly accumulated Rf (nx48)]



   

### 3. 	Defining peak-over-threshold extreme events
- **`Scripts/Creating_POT_Extremes_NTR.m`**

**Description**: The script uses the Peak over threshold (POT) approach to define events over thresholds for non-tidal residuals (NTR).
Functions: ut_solv, ut_reconstr are needed and should be downloaded through U-tide package (https://www.mathworks.com/matlabcentral/fileexchange/46523-utide-unified-tidal-analysis-and-prediction-functions). The years with substantial data gaps are removed manually. 

**Input**:
- `Water level data` :- [date-time, hourly water levels]

**output**: 
- `POT_NTR_and_NTR_timeseries.mat` :- structure[NTR.POT.WL_raw.Tide.MovA_WL.Threshold]


- **`Scripts/Creating_POT_Extremes_RF.m`**

**Description**: The script uses the Peak over threshold (POT) approach to define events over thresholds for basin average rainfall. The code calculates the threshold exceedances for rainfall accumulation times from 1 to 48.

**Input**:
- `Hourly_accumulation_Bias_Corrected_RF_data_”location”.mat`:- [date-time, hourly accumulated Rf (nx48)]

**output**: 
- `POT_Events_for_each_RF_Accumulation_time.mat`: -structure[events[POT.Threshold]




### 4. Two-way sampling
4. 1 Conditioned on NTR

**Description**: The script finds the maximum accumulated rainfall (for all the hourly accumulations from 1 to 48) within a given time window around the selected POT event when conditioning NTR.
- **`Scripts/Sampling_Con_NTR.m`**

**Input**:
- `Hourly_accumulation_Bias_Corrected_RF_data_”location”.mat`:- [date-time, hourly accumulated Rf (nx48)]
- `POT_NTR_and_NTR_timeseries.mat`:- structure[NTR.POT.WL_raw.Tide.MovA_WL.Threshold]

**output**:  
- `Maximum_RF_events_for_each_POT_NTR.ma`:-structure[Acumulation[Time_NTR.POT_NTR.Time_RF.Max_RF]

4. 2 Conditioned on RF

**Description**: The script finds the maximum NTR (for all the selected accumulation time) within a given time window around the selected POT event when conditioning RF.
- **`Scripts/Sampling_Con_RF.m`**

**Input**:
- `POT_NTR_and_NTR_timeseries.mat`:- structure[NTR.POT.WL_raw.Tide.MovA_WL.Threshold]
- `POT_Events_for_each_RF_Accumulation_time.mat`: -structure[events[POT.Threshold]


**output**: 
- `Maximum_NTR_events_for_each_POT_RF_for_”Selected_accumulation_time”_RF_acc.mat`:-structure [Time_NTR.Max_NTR.Time_RF.POT_RF]

## 5. Stratification
**Description**: The two conditional samples are stratified into two sets as: 1. The events caused by Tropical cyclones 2. The events that were not caused by tropical cyclones. The Hurdat 2  data set is used to identify the events induced by tropical cyclones.
    
5.1 Conditioned on NTR
- **`Scripts/Stratification_Con_NTR.m`**

**Input**: 
- `Maximum_RF_events_for_each_POT_NTR.mat`:-structure[Acumulation[Time_NTR.POT_NTR.Time_RF.Max_RF]
- `Cyclone_Track_data_from_1850.mat`:- [date-time, lat, lon, speed, distance from the city]

**output**: 
- `TC_events_conditioning_POT_NTR.mat`:- structure[event.[Time_NTR,NTR,Time_RF,RF]
- `ETC_events_conditioning_POT_NTR.mat`:- structure[event.[Time_NTR,NTR,Time_RF,RF]

5.2 Conditioned on RF
- **`Scripts/Stratification_Con_RF.m`**

**Input**:
- `Maximum_NTR_events_for_each_POT_RF_for_”Selected_accumulation_time”_RF_acc.mat`:-structure [Time_NTR.Max_NTR.Time_RF.POT_RF]
- `Cyclone_Track_data_from_1850.mat`:- [date-time, lat, lon, speed, distance from the city]

**output**: 
- `TC_events_conditioning_POT_RF.mat`:- structure[event.[Time_NTR,NTR,Time_RF,RF]
- `ETC_events_conditioning_POT_RF.mat`:- structure[event.[Time_NTR,NTR,Time_RF,RF]


## 6. Return period estimation from most likely estimate

**Description**: The stratified conditional samples are fitted with Generalized Pareto Distributions to estimate the return periods: 1. The following "XXX_main.mat" program is used for generating the scaled events based on the return periods calculated through the function "Uni_Return_level_calc.mat". 
    
- **`Scripts/Univariate_event_generation_main.m`**
- **`Scripts/Uni_Return_level_calc.m`**


**Input**: 
- `n_sim` = number of Bootstrapping samples required [number]
- `ET_NTR` = vector of POT Non-TC events conditioned on NTR [time_NTR, NTR, Time_RF,RF]
- `TC_NTR` = vector of POT TC events conditioned on NTR [time_NTR, NTR,Time_RF,RF]
- `ET_RF` = vector of POT Non-TC events conditioned on RF [time_NTR, NTR,Time_RF,RF]
- `TC_RF` = vector of POT TC events conditioned on RF [time_NTR, NTR,Time_RF,RF]
- `Thres_NTR` = single numeric vector of NTR threshold
- `Thres_RF` = single numeric vector of RF threshold
- `n_years` = Total number of years
- `Q_RP` = Vector of Return periods of interest
- `l_b_NTR` = Lover bound of discretized NTR space for combining two populations;
- `U_b_NTR` = Upper bound of discretized NTR space for combining two populations;
- `l_b_R`F = Lover bound of discretized RF space for combining two populations;
- `U_b_RF` = Upper bound of discretized RF space for combining two populations;


**output**: 

- `RL_NTR.mat`:- combined return levels of NTR, for return periods from 0 to 500 in 0.1 interval
- `RL_RF.mat`:- Combined return levels of RF, for return periods from 0 to 500 in 0.1 interval

- `Design_events_Cropped_Irene_.mat`:- structure[WL,RF,RP,Tide,NTR]


## 7. Return period estimation for the confidence interval

**Description**: The stratified conditional samples are fitted with Generalized Pareto Distributions to estimate the return periods: 1. The following "XXX_main_bootstrap.mat" program is used for generating the scaled events based on the return periods calculated through the function "Uni_Return_level_calc_with_boostrap.mat". 
    
- **`Scripts/Univariate_event_generation_main.m`**
- **`Scripts/Uni_Return_level_calc_with_bootstrap.m`**


**Input**: 
- `n_sim` = number of Bootstrapping samples required [number]
- `ET_NTR` = vector of POT Non-TC events conditioned on NTR [time_NTR, NTR, Time_RF,RF]
- `TC_NTR` = vector of POT TC events conditioned on NTR [time_NTR, NTR,Time_RF,RF]
- `ET_RF` = vector of POT Non-TC events conditioned on RF [time_NTR, NTR,Time_RF,RF]
- `TC_RF` = vector of POT TC events conditioned on RF [time_NTR, NTR,Time_RF,RF]
- `Thres_NTR` = single numeric vector of NTR threshold
- `Thres_RF` = single numeric vector of RF threshold
- `n_years` = Total number of years
- `Q_RP` = Vector of Return periods of interest
- `l_b_NTR` = Lover bound of discretized NTR space for combining two populations;
- `U_b_NTR` = Upper bound of discretized NTR space for combining two populations;
- `l_b_R`F = Lover bound of discretized RF space for combining two populations;
- `U_b_RF` = Upper bound of discretized RF space for combining two populations;
- `RL_NTR`, - `RL_RF` = calculated combined return level vectors for a given return period in the vector Q_RP


**output**: 
- `RL_NTR_com.mat`:- combined return levels of NTR, for return periods from 0 to 500 in 0.1 interval
- `RL_RF_com.mat`:- combined return levels of RF, for return periods from 0 to 500 in 0.1 interval
- `RL_TC_NTR.mat`:- return levels of TC induced NTR, for return periods from 0 to 500 in 0.1 intervals, for each bootstrapped sample
- `RL_ETC_NTR.mat`:- return levels of non-TC induced NTR, for return periods from 0 to 500 in 0.1 intervals, for each bootstrapped sample
- `RL_TC_RF.mat`:- return levels of TC-induced RF, for return periods from 0 to 500 in 0.1 intervals, for each bootstrapped sample
- `RL_ETC_RF.mat`:- return levels of non-TC induced RF, for return periods from 0 to 500 in 0.1 intervals, for each bootstrapped sample
- `RP_NTR.mat`:- vector of return periods considered [1:0.1:500]
- `RP_RF.mat`:- vector of return periods considered [1:0.1:500]

- `Design_events_Cropped_Irene_.mat`:- structure[WL,RF,RP,Tide,NTR]
  
The Cropped events are used as boundary conditions for the compound flood model. 


# The Flood Hazard Modeling using SFINCS model

In this section, we describe how to perfom the flood hazard analysis used in the paper "[Funding rules that promote equity in climate adaptation outcomes](https://osf.io/preprints/osf/6ewmu)." using the rainfall and water levels derived in the previous section. We used the SFINCS flood model, which is an open access and open source flood model developed by Deltares and can be found (10.5281/zenodo.8038534, https://github.com/Deltares/SFINCS) 

### 1. SFINCS Set Up

Documentation about the SFINCS model and how to set it up can be found (https://sfincs.readthedocs.io/en/latest/index.html). 
Information about the raw input data used to set up SFINCS and the sources of the data can be found in the Supporting Material of the paper *LINK TO SUP MAT*
We used ArcGIS pro to preprocess the DEM (clip to catchments domain), the Delft Dashboard (Van Ormondt et al. 2020 https://doi.org/10.2166/hydro.2020.092), and the Open Earth Tools (https://svn.oss.deltares.nl/repos/openearthtools/trunk/matlab/applications/sfincs/) to generate some of the input files of SFINCS (DEM, Subgrid tables, and msk file). 
To generate the Roughness file, we converted the land cover classes to Manning coefficients using the "Mean" value of the conversion table "Unique_Land_Classes_CN.xls" in the Data_Flood_Modeling. 
The input model files of SFINCS for Gloucester City (NJ, US) used in this study can be downloaded from (https://doi.org/10.5281/zenodo.14251309).

### 2. Boundary Conditions and Rainfall Forcing

As explain in the methodlogy of "[Funding rules that promote equity in climate adaptation outcomes](https://osf.io/preprints/osf/6ewmu).", we simulated flooding from different return period events of extreme water levels and rainfall independtly. 
The boudary files of water levels and rainfall fields for SFINCS from the events derived from the Statistical Framework can be generated using the Matlab script "SFINCS_BC_generator.m" and the function "sfincs_write_netcdf_amprfile.m" from the Open Earth Tools (https://svn.oss.deltares.nl/repos/openearthtools/trunk/matlab/applications/sfincs/). 

**Input**: 
- `Design_events_Cropped_Irene_.mat`
- `Lat_lon_Time_Irene.mat`
- `Gloucester_street_light_utm.tif` in "Data_Flood_modeling"

**Output**:
- `Rain_XXX.nc` (Rainfall field of event of RP XXX)
- `sfincs.bzs`  (Water Level time series of RP XXX)

### 3. Model File & Running SFINCS

In "Data_Flood_modeling" can be found one of the model files ("sfincs.inp") used for one of the SFINCS simulations. This file can be open with any text editor and contains all the information about the configuration of SFINCS used for this study. In order to run the different events, the fields "bzsfile" and "netamprfile" need to be modified accordingly with the name given to the water levels boudary file and rainfall field files generated from the previous section.  
For this study, we used the Version 2.0.0 Alpe D'Huez Release 2022 (10.5281/zenodo.8038534) of SFINCS.

### 4. Downscaling floodplains to subgrid resolution

For this study, SFINCS is configurated to run using the subgrid approach and thus the outputs of the simulations ("sfincs_map.nc") can be downscaled to 1m resolution. For that, we use the Matlab script "SUBGRID_Downscale.m" in "Script_Flood" that was written with help of Marteen Van Ormondt. That script requires the DEM ("dem_subgrid_1m_nbd.mat" in "Data_Flood_Modeling") and a number of Matlab functions from the Open Earth Tools (https://svn.oss.deltares.nl/repos/openearthtools/trunk/matlab/applications/sfincs/ and https://svn.oss.deltares.nl/repos/openearthtools/trunk/matlab/io/netcdf/snctools/).
The output of this code is the maximum water depths from the simulation at 1m resolution ("'Max_WD_.mat'")

**Input**: 
- `sfincs_map.nc`
- `dem_subgrid_1m_nbd.mat`
- `sfincs.inp`

**Output**: 
- `'Max_WD_.mat'` 

**List of Functions Required from Open Earth Tools**: "deblank2.m", "downscale_zs.m", "mc_bute.m", "nc_char.m", "nc_double.m", "nc_float.m", "nc_getattsinfo_tmw.m", "nc_int.m", "nc_int64.m", "nc_nat.m", "nc_short.m", "nc_ubyte.m", "nc_unit.m", "nc_unit64.m", "nc_ushort.m", "nc_varget.m", "nc_varget_tmw.m", "nc_vargetr.m", "sfincs_get_values_for_dem.m", "sfincs_make_grid.m", sfincs_read_indices_for_dem.m", "sfincs_read_input.m", "snc_format.m", "snc_get_indexing.m", "snc_read_backend.m" 

### 5. Combination of Coastal and Pluvial Floodplains

In this study, the flood hazard is simulated separately for each driver and then their floodplains are combined using the Matlab script "

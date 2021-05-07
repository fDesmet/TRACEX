#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Creation: May 2021
Authors: Flora Desmet and Eike Koehn, ETH Zuerich

Description:
Define settings, variables and paths prior to running launch_TRACEX.py 
This includes:

1) Variables and thresholds for detection 
- the variable name (in model outputs) of the variable used for the detection and the type of threshold: 'var', 'thresh' + if compound extremes: 'var2', 'thresh2'

2) Settings 
- the different settings (depth limit, minimum duration, compound/single extremes, moving/fixed basline, type of connectivity (edge/corners), restart)

3) Model outputs 
- the time period: years vector 'years' and associated daily vector 't' 
- the path to the daily model outputs (one file per year): 'dir_in' 
- the setup used in the filenames of model output (format "setup_year_avg.nc")

4) Output folder 
- name of output folder: 'outpath'

5) Grid data
- grid file to get number of z levels (romsGrd.NZ), area
- 'mask3D' masking out the regions not used in the detection 
- field of grid cells height 'dz'
- field of grid cells depth 'zlevs_rho'

6) Thresholds fields (if not absolute thresholds)
- if thresh (or thresh2) are strings (= not absolute thresholds): the threshold fields 'thresh_main' + if compound 'thresh_main2' + if moving baseline 'slope_main' field (+'slope_main2')
- the threshold to be used for additional diagnostics such as intensity of O2 or temperature: 'thresh_omega', 'thresh_pH', 'thresh_O2', 'thresh_T'

7) Load the 365 and 366 days climatology for temperature to compute delta_temp in extremes 
- 'clm_file_output' (netCDF4 file)
- 'path_in' (path where the Temp_climatologies are located)
- 'filename_clm_366j' (name of a .npy file)

8) If needs additional properties for events such as propagation distance,  mean vertical occupation... This will run only after the detection has completed
- longitude 2D field
- latitude 2D field
- distance the coast 2D field

"""

# Load classic python modules
import numpy as np
from datetime import date
import sys
import netCDF4
import os
# Load package functions
from TRACEX_functions import *


###################################################################################################
##############   1) Choose variables and thresholds for detection   ###############################
###################################################################################################
var = 'omega' # the detection variable
var2 = 'pH'   # the second detection variable if compound event
thresh = 'stat_1th_omega' # the variable threshold: used to pick the right variable for detection
thresh2 = 'stat_1th_pH'   # the second variable threshold: used to pick the right variable for detection if compound 

# if statistical threhsold give the percentile used: used to load the proper file
perc_choice = '1'
# give the suffixe naming of the file which contain the threshold
suffixe_thresh_file = '1984_2019_hc003_daily_pactcs30'

# the threshold value or the threshold filename if defined in an array for each variable used in additional
# diagnostics (not used as a criteria for detection). If 0 then no threshold is taken into account. 
thresh_omega = 'stat_1th_omega' # 1 for example if absolute threhsold of 1; 0 if no detection on this variable
thresh_pH = 'stat_1th_pH' 
thresh_O2 = 0
thresh_T = 0

###################################################################################################
#############################     2)   Choose settings       ######################################
###################################################################################################

# how do extremes connect?
edges_or_corners = 1 
edges_only = np.abs(1-edges_or_corners)

# do you want two variables to be extreme simultaneously (compound extremes)? Set to True.
two_threshs = True
# do you want the threshold to be based on a moving baseline? Set to True.
moving_baseline = False
# is your model outputs split into two time periods kept in different folders? Set to True.
use_extended_hindcast = True
# what is your lower depth limit? Detection is from surface to "depth"
depth = 150       # in meters
# what is your lower duration limit? Duration>= minDuration
minDuration = 5   # in days
# are you restarting the detection from a year in the middle of your period? Set to 1. 
restart = 0
# do you want additional diagnostics to be output such as oxgen and temperature informations? Set to True.
additional_diags = True
# do you want to add additional properties such as propagation distance, mean vertical occupation...? Set to True. This will run only after the detection has completed
additional_properties = True

###################################################################################################
#############################  3)  Choose model outputs      ######################################
###################################################################################################

# time period
years = np.arange(1984,2020,1)
# get daily vector
t = np.arange(date(years[0],1,1).toordinal(),date(years[-1]+1,1,1).toordinal()+1)

# directory of your model outputs
dir_in = '/net/kryo/work/fdesmet/roms/output/pactcs30/hc003_daily_pactcs30/avg'
# if model outputs split into two folders, second directory
if use_extended_hindcast:
    dir_in2 = '/net/kryo/work/fdesmet/roms/output/pactcs30/hc003_daily_pactcs30_extended_2013_2019/avg'

# setup for filenames of model output (format "setup_year_avg.nc")
setup = 'pactcs30'

###################################################################################################
#############################  4) Choose your output folder   #####################################
###################################################################################################
outpath = '/net/kryo/work/fdesmet/figures/acidification_diagnostic/detection_output/{}_{}/detection001/'.format(years[0],years[-1])

# If the outpath does not exist, create the directory     
try:
    os.makedirs(outpath)
except FileExistsError:
    # directory already exists
    pass

###################################################################################################
########### 5) Load grid, mask and compute z arrays if needed using ROMS_tools functions ##########
###################################################################################################  
# Your model grid file
grdfile = '/net/kryo/work/martinfr/Roms/Inputs/pactcs30/pactcs30_grd.nc'
# Your file where grid informations are stored such as size of grid cells
clmfile = '/net/kryo/work/martinfr/Roms/Inputs/pactcs30/N64ts10tb4hc250_grd_merged_SiO3_PO4_fix/pactcs30_clm.nc'
# Using the getGrid function to get Lat/lon of the grid and area
romsGrd = getGrid(grdfile)
romsGrd.getAttrs(clmfile)
romsGrd.getLatLon()
romsGrd.getArea()

# Which domain do you want to run the detection on?
fmask = '/net/kryo/work/fdesmet/Data/Masks/north_pacific_CCS_option1.nc'
ncfmasks = netCDF4.Dataset(fmask, 'r')
areamask = ncfmasks.variables['mask'][:]
#reduce domain size: from the mask, reduce the domain of the detection 
id1 = np.where(areamask==1)  
idx_eta_min = np.min(id1[0])
idx_eta_max = np.max(id1[0])
idx_xi_min = np.min(id1[1])
idx_xi_max = np.max(id1[1])
indices = [idx_eta_min,idx_eta_max,idx_xi_min,idx_xi_max]
#reduce mask array
areamask = areamask[idx_eta_min:idx_eta_max,idx_xi_min:idx_xi_max]
# convert mask from float32 array to boolean array 
mask2D = areamask.astype(bool)
mask3D = np.repeat(mask2D[np.newaxis, :, :], romsGrd.NZ, axis=0)
del areamask

# Load field of grid cells height and depths 
path_in = '/net/kryo/work/fdesmet/Data/pactcs30/'
dz = load_npy_file_3D(path_in,'z_array_for_detection/dz.npy',indices,pickle=True)
zlevs_rho = load_npy_file_3D(path_in,'z_array_for_detection/zlevs_rho.npy',indices,pickle=True)
# Load area (in square meters)
area = romsGrd.area[idx_eta_min:idx_eta_max,idx_xi_min:idx_xi_max]

###################################################################################################
################ 6) If statistical threshold: load arrays of thresholds ###########################
###################################################################################################

# the threshold for the detect_var is already chosen above -> just load the corresponding dataset

if moving_baseline:
    # choose here threshold files used for the extremes detection
    print('Choosing moving baseline...')
    if thresh == 'stat_1th_omega': # variable for detection is omega aragonite
        # Load threshold computed on detrended timeserie and trend (slope) of the timeserie
        print('Loading the 3D array of the {}th percentile of omega values...'.format(perc_choice))
        thresh_main = netCDF4.Dataset('/net/kryo/work/fdesmet/Data/pactcs30/OA_thresh/linear_reg_omega_arag_offl.nc','r').variables['thresh_on_res'][:,idx_eta_min:idx_eta_max,idx_xi_min:idx_xi_max]
        slope_main = netCDF4.Dataset('/net/kryo/work/fdesmet/Data/pactcs30/OA_thresh/linear_reg_omega_arag_offl.nc','r').variables['slope'][:,idx_eta_min:idx_eta_max,idx_xi_min:idx_xi_max]     
        print('Shape of main threhsold is {}'.format(thresh_main.shape))   
        # if compound extreme: second variable is pH
        if two_threshs:
            # Load threshold computed on detrended timeserie and trend (slope) of the timeserie
            print('Also loading the 3D array of the {}th percentile of pH values...'.format(perc_choice))
            thresh_main2 = netCDF4.Dataset('/net/kryo/work/fdesmet/Data/pactcs30/OA_thresh/linear_reg_pH_offl.nc','r').variables['thresh_on_res'][:,idx_eta_min:idx_eta_max,idx_xi_min:idx_xi_max]
            slope_main2 = netCDF4.Dataset('/net/kryo/work/fdesmet/Data/pactcs30/OA_thresh/linear_reg_pH_offl.nc','r').variables['slope'][:,idx_eta_min:idx_eta_max,idx_xi_min:idx_xi_max]   
            print('Shape of secondary threhsold is {}'.format(thresh_main2.shape))     

    elif thresh == 'stat_1th_pH': # variable for detection is pH
        # Load threshold computed on detrended timeserie and trend (slope) of the timeserie
        print('Loading the 3D array of the {}th percentile of pH values...'.format(perc_choice))
        thresh_main = netCDF4.Dataset('/net/kryo/work/fdesmet/Data/pactcs30/OA_thresh/linear_reg_pH_offl.nc','r').variables['thresh_on_res'][:,idx_eta_min:idx_eta_max,idx_xi_min:idx_xi_max]
        slope_main = netCDF4.Dataset('/net/kryo/work/fdesmet/Data/pactcs30/OA_thresh/linear_reg_pH_offl.nc','r').variables['slope'][:,idx_eta_min:idx_eta_max,idx_xi_min:idx_xi_max]      
        print('Shape of main threhsold is {}'.format(thresh_main.shape))  
        # if compound extreme: second variable is omega aragonite
        if two_threshs:
            # Load threshold computed on detrended timeserie and trend (slope) of the timeserie
            print('Also loading the 3D array of the {}th percentile of omega values...'.format(perc_choice))
            thresh_main2 = netCDF4.Dataset('/net/kryo/work/fdesmet/Data/pactcs30/OA_thresh/linear_reg_omega_arag_offl.nc','r').variables['thresh_on_res'][:,idx_eta_min:idx_eta_max,idx_xi_min:idx_xi_max]
            slope_main2 = netCDF4.Dataset('/net/kryo/work/fdesmet/Data/pactcs30/OA_thresh/linear_reg_omega_arag_offl.nc','r').variables['slope'][:,idx_eta_min:idx_eta_max,idx_xi_min:idx_xi_max]    
            print('Shape of secondary threhsold is {}'.format(thresh_main2.shape))    

    else:
        raise NotImplemented("Sorry, this threshold {} is not implemented".format(thresh)) 
        
    
    # choose here threshold files for additional diagnostics computation: omega,pH,O2,T
    # if 0 then no threhsold is taken into account.
    if thresh_omega == 'stat_1th_omega': 
        # Load threshold computed on detrended timeserie and trend (slope) of the timeserie
        thresh_omega = netCDF4.Dataset('/net/kryo/work/fdesmet/Data/pactcs30/OA_thresh/linear_reg_omega_arag_offl.nc','r').variables['thresh_on_res'][:,idx_eta_min:idx_eta_max,idx_xi_min:idx_xi_max]
        slope_omega = netCDF4.Dataset('/net/kryo/work/fdesmet/Data/pactcs30/OA_thresh/linear_reg_omega_arag_offl.nc','r').variables['slope'][:,idx_eta_min:idx_eta_max,idx_xi_min:idx_xi_max]        
    if thresh_pH == 'stat_1th_pH':
        # Load threshold computed on detrended timeserie and trend (slope) of the timeserie
        thresh_pH = netCDF4.Dataset('/net/kryo/work/fdesmet/Data/pactcs30/OA_thresh/linear_reg_pH_offl.nc','r').variables['thresh_on_res'][:,idx_eta_min:idx_eta_max,idx_xi_min:idx_xi_max]
        slope_pH = netCDF4.Dataset('/net/kryo/work/fdesmet/Data/pactcs30/OA_thresh/linear_reg_pH_offl.nc','r').variables['slope'][:,idx_eta_min:idx_eta_max,idx_xi_min:idx_xi_max]        
    if thresh_O2 == 'stat_1th_O2':
        thresh_O2 = load_npy_file_3D(path_in,'OA_thresh/O2_thresh_1thperc_{}.npy'.format(suffixe_thresh_file),indices)
    if thresh_T == 'stat_1th_T':
        thresh_T = load_npy_file_3D(path_in,'OA_thresh/T_thresh_1thperc_{}.npy'.format(suffixe_thresh_file),indices)

# if using thresholds based on a fixed baseline
else:
    # choose here threshold files used for the extremes detection
    if thresh == 'stat_1th_omega':  # variable for detection is omega aragonite
        print('Loading the 3D array of the {}th percentile of omega values...'.format(perc_choice))
        thresh_main = load_npy_file_3D(path_in,'OA_thresh/omega_thresh_{}thperc_{}.npy'.format(perc_choice,suffixe_thresh_file),indices)
        # if compound extreme: second variable is pH
        if two_threshs:
            print('Also loading the 3D array of the {}th percentile of pH values...'.format(perc_choice))
            thresh_main2 = load_npy_file_3D(path_in,'OA_thresh/pH_thresh_{}thperc_{}.npy'.format(perc_choice,suffixe_thresh_file),indices)
    elif thresh == 'stat_1th_pH':  # variable for detection is pH
        print('Loading the 3D array of the {}th percentile of pH values...'.format(perc_choice))
        thresh_main = load_npy_file_3D(path_in,'OA_thresh/pH_thresh_{}thperc_{}.npy'.format(perc_choice,suffixe_thresh_file),indices)
        # if compound extreme: second variable is omega aragonite
        if two_threshs:
            print('Also loading the 3D array of the {}th percentile of omega values...'.format(perc_choice))
            thresh_main2 = load_npy_file_3D(path_in,'OA_thresh/omega_thresh_{}thperc_{}.npy'.format(perc_choice,suffixe_thresh_file),indices)
    else:
        raise NotImplemented("Sorry, this threshold {} is not implemented".format(thresh)) 
        
        
    # choose here threshold files for additional diagnostics computation: omega,pH,O2,T
    # if 0 then no threhsold is taken into account.
    if thresh_omega == 'stat_1th_omega':
        thresh_omega = load_npy_file_3D(path_in,'OA_thresh/omega_thresh_{}thperc_{}.npy'.format(perc_choice,suffixe_thresh_file),indices)
    if thresh_pH == 'stat_1th_pH':
        thresh_pH = load_npy_file_3D(path_in,'OA_thresh/pH_thresh_{}thperc_{}.npy'.format(perc_choice,suffixe_thresh_file),indices)
    if thresh_O2 == 'stat_1th_O2':
        thresh_O2 = load_npy_file_3D(path_in,'OA_thresh/O2_thresh_1thperc_{}.npy'.format(suffixe_thresh_file),indices)
    if thresh_T == 'stat_1th_T':
        thresh_T = load_npy_file_3D(path_in,'OA_thresh/T_thresh_1thperc_{}.npy'.format(suffixe_thresh_file),indices)
        


###################################################################################################
#### 7) Load the 365 days climatology for temperature to compute delta_temp in extremes ###########
###################################################################################################        
clm_file_output = '/net/kryo/work/fdesmet/roms/output/pactcs30/hc003_daily_pactcs30_extended_2013_2019/clim/daily_clim_1984_2019_temp_smooth.nc'
filename_clm_366j = 'Temp_climatologies/daily_clim_1984_2019_temp_smooth_366j.npy'


###################################################################################################
#### 8) Load the latitude, longitude and distance to the coast arrays to compute ##################
####    additional properties of the extremes such as propagation distance...    ##################
################################################################################################### 
# This will run only after the detection has completed
if additional_properties:
    flonlat = '/net/kryo/work/martinfr/Roms/Inputs/pactcs30/pactcs30_grd.nc'
    lon_varname = 'lon_rho'
    lat_varname = 'lat_rho'
    fdcoast = '/net/kryo/work/fdesmet/Data/pactcs30/pactcs30_dist2coast.nc'
    dcoast_varname = 'dcoast'

   
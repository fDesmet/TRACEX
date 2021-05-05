#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Creation: May 2021
Authors: Flora Desmet and Eike Koehn, ETH Zuerich

Description:
Main script for the detection of extremes in space and time.
An step-by-step example is given for a 2D field in the jupyter notebook: TRACEX_unit_test/2D_example_TRACEX.ipynb
"""

#########################
## Load needed modules ##
#########################
import numpy as np
from datetime import date
import netCDF4
import time
from collections import Counter
import scipy.ndimage as ndimage
import pandas as pd
import os
import pickle
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pdb
## TRACEX functions 
from TRACEX_functions import *

########################
## Load all variables ##
########################
from define_TRACEX_variables import *

# define chunck length (can vary depending on available memory of the machine)
chunklen = int(31) 
numchunk = int(np.ceil(366./chunklen)+1)

#mask values below chosen depth (default is 100m) 
# Check the sign of zrho levs file
if np.nansum(zlevs_rho)<0:
    zlevs_rho = zlevs_rho*-1
mask3D[zlevs_rho>depth] = False

#####################################
## definition of the core rountine ##
#####################################
def tracex():

    #### Define restart filenames where dictionnaries will be saved at the end of each year to be able to restart the script at any year within the study period
    frestartname_charac = outpath + 'extreme_detection_restart_file.pck'
    frestartname_coords= outpath + 'extreme_detection_extremes_dict_restart_file.pck'

    #### Load filelist for detection and data of restart if it is a restart
    if not restart:
        filelist_name = ['{}_{}_avg.nc'.format(setup, yr) for yr in years]
        t_ind=0
        yr = 0        
        if use_extended_hindcast:
            # Add path to files
            filelist = []
            for yr_i in years:
                if yr_i<2013:
                    filelist.append(dir_in+'/{}_{}_avg.nc'.format(setup, yr_i))
                else:
                    print('For year {} use {}'.format(yr_i,dir_in2))
                    filelist.append(dir_in2+'/{}_{}_avg.nc'.format(setup, yr_i))
    else:
        num = np.load(outpath + 't_ind_yr_n_events_prev_restart_file.npy')
        yr = num[1] + 1 # The year saved in num is the last saved step --> years[num[1]] has already been processed and saved
        t_ind = num[0]
        print('t_ind {} date: {}'.format(t_ind,date.fromordinal(t[t_ind])))
        restart_yr = years[yr]
        print('Year {}'.format(restart_yr))
        years_restart = np.arange(restart_yr,years[-1]+1,1)
        filelist_name  = ['{}_{}_avg.nc'.format(setup, yr_i) for yr_i in years_restart]    
        if use_extended_hindcast:   
            # Add path to files
            filelist = []
            for yr_i in years_restart:
                if yr_i<2013:
                    filelist.append(dir_in+'/{}_{}_avg.nc'.format(setup, yr_i))
                else:
                    print('For year {} use {}'.format(yr_i,dir_in2))
                    filelist.append(dir_in2+'/{}_{}_avg.nc'.format(setup, yr_i))

    if not use_extended_hindcast:
        filelist = [dir_in+'/'+f for f in filelist_name]

    
    
    ### Write detection settings into a text file   
    file_settings = open(outpath+'detection_settings.txt','w') 
    file_settings.writelines(['Detection with DeteXt branch Flora \n','Detection on variable {} \n'.format(var), 'Type of threshold used :{} \n'.format(thresh), \
                              'Detection on two variables : {} \n'.format(two_threshs), \
                              'Moving baseline: {} \n'.format(moving_baseline), \
                              'Using extended hindcast: {} \n'.format(use_extended_hindcast), \
                              'Threshold computed over {} \n'.format(suffixe_thresh_file), 'Mask used : {} \n'.format(fmask), \
                              'Outputs from directory : {} \n'.format(dir_in), 'Writting detection output in : {} \n'.format(outpath), \
                              'Checking additional diags : {} \n'.format(additional_diags), 'Lower depth limit : {} \n'.format(depth), \
			      'Minimum duration limit (days) : {} \n'.format(minDuration), 'Link by edges or corners : {} \n'.format(edges_or_corners), \
			      'Chunk length used to load data (days) : {} \n'.format(chunklen)])
    file_settings.close()
    
    
    
    ### Define properties of event (keys of dictionary)
    if not restart:
        # create the dictionary with the global properties
        mhw = init_property_extreme_dict(additional_diags)
            
        # Create the dictionnary to keep coordinates of grid cells of each events
        extremes = init_extremes_dict()
                
        ## Write detection steps into a text file   
        file_log = open(outpath+'detection_log.txt','w')

    else:
        ## Write detectionsteps into a text file   
        file_log = open(outpath+'detection_log.txt','a')
        file_log.writelines(['Restart \n','Year {} \n'.format(restart_yr),'Date of t_ind : {} {}\n'.format(t_ind,date.fromordinal(t[t_ind]))])


    ###############################################################################################################################   
    #                                                                                                                             #
    #      *****    *********    ******   *******    *********       ******     *******  *********  *******     *****  *********  #
    #    ***   ***  *  ***  *   **  ***   ***   ***  *  ***  *       ***  ***   *******  *  ***  *  *******   *******  *  ***  *  #
    #    ***           ***     **   ***   ***   ***     ***          ***   ***  ***         ***     ***      ***          ***     #
    #      *****       ***    ***   ***   ***   **      ***          ***   ***  *******     ***     *******  ***          ***     #
    #          ***     ***    *********   ******        ***          ***   ***  ***         ***     ***      ***          ***     #
    #    ***   ***     ***    ***   ***   ***  **       ***          ***  ***   *******     ***     *******   *******     ***     #
    #      *****       ***    ***   ***   ***   ***     ***          ******     *******     ***     *******     *****     ***     #
    #                                                                                                                             #
    ###############################################################################################################################
          
    #Loop over years
    for f in filelist:     
        print('Processing file ' + f + '...')
        file_log.writelines(['Processing file {} \n'.format(f)])
        #retrieve data if restart and first year of restart
        if restart and f==filelist[0]:
            #import from dictionnary and events_prev and detect_var_prev and n_events_prev
            mhw = pickle.load(open(frestartname_charac,'rb'),encoding='ascii')
            extremes = pickle.load(open(frestartname_coords,'rb'),encoding='ascii')
            events_prev = np.load(outpath + 'events_prev_restart_file.npy')
            detect_var_prev = np.load(outpath + 'detect_var_prev_restart_file.npy')
            if two_threshs:
                detect_var2_prev = np.load(outpath + 'detect_var2_prev_restart_file.npy')
            n_events_prev = num[2]

        ###########################################################################################
        ##  Loop over the chunks to speed up the script (by avoiding to load everything at once) ##
        ###########################################################################################      
        for m in range(1,numchunk):
            file_log.writelines(['Processing m = {}'.format(m)])
            #if we are processing the first file not in restart mode   
            if not restart and f==filelist[0] and m==1:          
                nc_in = netCDF4.Dataset(f,'r')
                detect_var,mask4D,T_file,omega_arag,pH,O2,T = load_var_datachunks_optimized(nc_in,var,mask3D,m,chunklen,numchunk,indices,additional_diags)
                if two_threshs:
                    detect_var2 = load_var_2_datachunks_optimized(nc_in,var2,mask3D,m,chunklen,numchunk,indices)
                print(T_file)
                print(t_ind)
                if additional_diags:
                    if len(nc_in.variables['temp'][:,0,0,0]) == 365:
                        temp_clm_ident=netCDF4.Dataset(clm_file_output, 'r')
                        T = T - load_datachunk(temp_clm_ident,'temp_smooth',m,chunklen,numchunk,indices) 
                    else:
                        T = T - load_npy_file(path_in,filename_clm_366j,m,chunklen,numchunk,indices) 
                    
                if thresh == 'stat_1th_omega' or thresh == 'stat_1th_pH' or thresh == 'stat_01th_fitted_O2':
                    # generate 4D array of thresh
                    if not moving_baseline:
                        thresh_4D = np.repeat(thresh_main[np.newaxis,: , :, :], detect_var.shape[0], axis=0)
                    else:
                        # Load the 4D moving baseline threhsold
                        thresh_4D = load_threshold_linear_trend_removed(thresh_main,slope_main,range(t_ind,t_ind+detect_var.shape[0])) 
                        
                    # Detect events through space and time: if exceed_bool is positive it means 
                    # that detect_var is below the 1st percentile and therefore it is an extreme : True
                    # !! To change if extreme is defined as above a certain threshold !!
                    exceed_bool = thresh_4D - detect_var
                    # Intensity
                    I = thresh_4D - detect_var
                    del thresh_4D  
                    
                    if two_threshs:
                        if not moving_baseline:
                            thresh_4D = np.repeat(thresh_main2[np.newaxis,: , :, :], detect_var2.shape[0], axis=0)
                        else:
                            # Load the 4D moving baseline threhsold
                            thresh_4D = load_threshold_linear_trend_removed(thresh_main2,slope_main2,range(t_ind,t_ind+detect_var2.shape[0])) 
                            
                        # Detect events through space and time: if exceed_bool is positive it means 
                        # that detect_var is below the 1st percentile and therefore it is an extreme : True
                        # !! To change if extreme is defined as above a certain threshold !!
                        exceed_bool2 = thresh_4D - detect_var2
                        # Intensity
                        I2 = thresh_4D - detect_var2
                        del thresh_4D  
                    
                else:
                    ###############################################################################
                    ##################### Absolute detect_var threshold  ##########################
                    ###############################################################################
                    # Detect events through space and time: if exceed_bool is positive it means 
                    # that detect_var is smaller than thresh and therefore it is an extreme : True
                    # !! To change if extreme is defined as above a certain threshold !!
                    print(thresh)
                    exceed_bool = thresh - detect_var
                    # Intensity
                    I = thresh - detect_var 
                    
                
                exceed_bool[exceed_bool<=0] = False
                exceed_bool[exceed_bool>0] = True
                exceed_bool[np.isnan(exceed_bool)]=False     # Discard land mass points
                # if second threshold (compound extremes), process the second exceed_bool array and multiply the two
                if two_threshs:
                    exceed_bool2[exceed_bool2<=0] = False
                    exceed_bool2[exceed_bool2>0] = True
                    exceed_bool2[np.isnan(exceed_bool2)]=False  # Discard land mass points
                    # Multiply the two boolean arrays
                    exceed_bool = np.multiply(exceed_bool,exceed_bool2)
                    
                # Find contiguous regions of exceed_bool = True
                if edges_only:
                    events, n_events = ndimage.label(exceed_bool)
                elif edges_or_corners:
                    dims = len(np.shape(exceed_bool))
                    s = ndimage.generate_binary_structure(dims,dims)  # Make diagonal connections
                    events, n_events = ndimage.label(exceed_bool,s)

                # Save for next chunk 
                detect_var_prev = np.reshape(detect_var[-1,:,:,:],(1,detect_var.shape[1],detect_var.shape[2],detect_var.shape[3]))
                if two_threshs:
                    detect_var2_prev = np.reshape(detect_var2[-1,:,:,:],(1,detect_var2.shape[1],detect_var2.shape[2],detect_var2.shape[3]))
               
                del exceed_bool
                
                mhw['n_events'] = n_events
                print('There are {} events in chunk {} of year {} \n'.format(n_events,m,years[yr]))
                file_log.writelines(['There are {} events in chunk {} of year {} \n'.format(n_events,m,years[yr])])

                ##################
                # Loop over events
                ##################
                for ev in range(1,n_events+1):  
                    temp = np.where(events == ev)
                    # temp is a tuple of 4 arrays of int with the index of the location 
                    # of the events: temp[0] is an array of all the time index where ev is present
                    
                    ##################################################
                    # Check duration
                    ##################################################
                    # Condition on duration: if not al least "minDuration" different times then not long enough and if the last time of the event is not the last time in the chunk (otherwise could continue in the next chunk and be long enough) 
                    if len(Counter(temp[0]).keys())<minDuration and max(temp[0],default=0)!=I.shape[0]-1:
                        # Discard this event
                        events[events == ev] = 0   
                        continue
                    
                    ##################################################
                    # Otherwise, fill the dictionnary with the event
                    ##################################################
                    extremes,mhw,v_gridded,v = fill_extreme_dictionaries_with_core_chars(mhw,extremes,ev,temp,I,t_ind,t,area,dz,zlevs_rho,detect_var)
                    
                    if additional_diags:
                        # Store properties for additional diagnostics
                        if not moving_baseline:
                            mhw = fill_extreme_dictionaries_with_add_chars(mhw,temp,omega_arag,pH,O2,T,v_gridded,v,thresh,thresh_omega=thresh_omega,thresh_pH=thresh_pH,thresh_O2=thresh_O2)  
                        else:
                            # Load the 4D moving baseline threshold
                            thresh_omega_tmp = load_threshold_linear_trend_removed(thresh_omega,slope_omega,range(t_ind,t_ind+omega_arag.shape[0])) 
                            thresh_pH_tmp = load_threshold_linear_trend_removed(thresh_pH,slope_pH,range(t_ind,t_ind+pH.shape[0])) 
                            mhw = fill_extreme_dictionaries_with_add_chars(mhw,temp,omega_arag,pH,O2,T,v_gridded,v,thresh,thresh_omega=thresh_omega_tmp,thresh_pH=thresh_pH_tmp,thresh_O2=thresh_O2)    
                            del thresh_omega_tmp,thresh_pH_tmp

                # Keep track of the global time indice over the entire detection 
                t_ind += T_file
                
            ###################################################################              
            #### If not the first chunk of first year or if it is a restart ###
            ###################################################################            
            else:
                # If processing the first file in restart mode   
                if restart and f==filelist[0] and m==1: 
                    print('Data already loaded from restart...')
                else:
                    # Save last time step of previous events in last chunk if not restart time
                    events_prev = events[-1,:,:,:]
                    # Save the previous total number of events
                    n_events_prev = mhw['n_events']    
                    
                # Load new data 
                nc_in_curr = netCDF4.Dataset(f,'r')                
                detect_var_curr,mask4D,T_file,omega_arag,pH,O2,T = load_var_datachunks_optimized(nc_in_curr,var,mask3D,m,chunklen,numchunk,indices,additional_diags)
                if two_threshs:
                    detect_var2_curr = load_var_2_datachunks_optimized(nc_in_curr,var2,mask3D,m,chunklen,numchunk,indices)
                    
                print(T_file)
                print(t_ind)

                if additional_diags:
                    temp_clm_ident=netCDF4.Dataset(clm_file_output, 'r')
                    if len(nc_in_curr.variables['temp'][:,0,0,0]) == 365:
                        T = T - load_datachunk(temp_clm_ident,'temp_smooth',m,chunklen,numchunk,indices) # Temperature variable name in climatology file
                    else:
                        T = T - load_npy_file(path_in,filename_clm_366j,m,chunklen,numchunk,indices)                  
       
                if thresh == 'stat_1th_omega' or thresh == 'stat_1th_pH'  or thresh == 'stat_01th_fitted_O2':
                    # generate 4D array of thresh
                    if not moving_baseline:
                        thresh_4D = np.repeat(thresh_main[np.newaxis,: , :, :], detect_var_curr.shape[0], axis=0)
                    else:
                        # Load the 4D moving baseline threhsold
                        thresh_4D = load_threshold_linear_trend_removed(thresh_main,slope_main,range(t_ind,t_ind+detect_var_curr.shape[0])) 
          
                    # Intensity
                    I = thresh_4D - detect_var_curr
                    #detect events through space and time: if exceed_bool is positive it means 
                    #that detect_var is below the 1st percentile and therefore it is an extreme : True
                    # !! To change if extreme is defined as above a certain threshold !!
                    detect_var = np.concatenate((detect_var_prev,detect_var_curr),axis=0)
                    thresh_4D = np.concatenate((np.reshape(thresh_4D[-1,:,:,:],(1,thresh_4D.shape[1],thresh_4D.shape[2],thresh_4D.shape[3])),thresh_4D),axis=0)  
                    exceed_bool = thresh_4D - detect_var
                    del thresh_4D
                    if two_threshs:
                        if not moving_baseline:
                            thresh_4D = np.repeat(thresh_main2[np.newaxis,: , :, :], detect_var2_curr.shape[0], axis=0)
                        else:
                            # Load the 4D moving baseline threhsold
                            thresh_4D = load_threshold_linear_trend_removed(thresh_main2,slope_main2,range(t_ind,t_ind+detect_var2_curr.shape[0])) 
                            
                        #Strorage of intensity
                        I2 = thresh_4D - detect_var2_curr
                        #detect events through space and time: if exceed_bool is positive it means 
                        #that detect_var is below the 5th percentile and therefore it is an extreme : True
                        # !! To change if extreme is defined as above a certain threshold !!
                        detect_var2 = np.concatenate((detect_var2_prev,detect_var2_curr),axis=0)
                        thresh_4D = np.concatenate((np.reshape(thresh_4D[-1,:,:,:],(1,thresh_4D.shape[1],thresh_4D.shape[2],thresh_4D.shape[3])),thresh_4D),axis=0)  
                        exceed_bool2 = thresh_4D - detect_var2
                        del thresh_4D                      
    
                else:
                    ###################################################################
                    ### Absolute detect_var threshold
                    ###################################################################
                    #  Intensity
                    I = thresh - detect_var_curr
                    #detect events through space and time: if exceed_bool is positive it means 
                    #that detect_var is below 1 and therefore it is an extreme : True
                    # !! To change if extreme is defined as above a certain threshold !!
                    detect_var = np.concatenate((detect_var_prev,detect_var_curr),axis=0)
                    exceed_bool = thresh - detect_var
    
                
                exceed_bool[exceed_bool<=0] = False
                exceed_bool[exceed_bool>0] = True
                exceed_bool[np.isnan(exceed_bool)]=False   # Discard land mass points
                # if second threshold (compound extremes), process the second exceed_bool array and multiply the two
                if two_threshs:
                    exceed_bool2[exceed_bool2<=0] = False
                    exceed_bool2[exceed_bool2>0] = True
                    exceed_bool2[np.isnan(exceed_bool2)]=False    # Discard land mass points
                    # Multiply the two boolean arrays
                    exceed_bool = np.multiply(exceed_bool,exceed_bool2)
                    
                # Find contiguous regions of exceed_bool = True
                if edges_only:    
                    events_tmp, n_events = ndimage.label(exceed_bool)
                elif edges_or_corners:
                    dims = len(np.shape(exceed_bool))
                    s = ndimage.generate_binary_structure(dims,dims)  # Make diagonal connections
                    events_tmp, n_events = ndimage.label(exceed_bool,s)                
                    
                events_tmp[events_tmp != 0] += n_events_prev
                
                # Remove unecessary variables
                detect_var_prev = np.reshape(detect_var[-1,:,:,:],(1,detect_var.shape[1],detect_var.shape[2],detect_var.shape[3]))
                del detect_var,exceed_bool
                if two_threshs:
                    detect_var2_prev = np.reshape(detect_var2[-1,:,:,:],(1,detect_var2.shape[1],detect_var2.shape[2],detect_var2.shape[3]))
                    del detect_var2                    
                
                mhw['n_events'] += n_events
                print('There are {} events in chunk {} of year {} \n'.format(n_events,m,years[yr]))
                file_log.writelines(['There are {} events in chunk {} of year {} \n'.format(n_events,m,years[yr])])

                ###################################################################              
                #### Check for continuation of events to merge them
                ###################################################################                
                #global min/max of the new chunk without zeros
                if len(np.nonzero(events_tmp)[0])>0:
                    #Returns a tuple of arrays, one for each dimension of a, 
                    #containing the indices of the non-zero elements in that dimension.
                    m_global = np.min(events_tmp[np.nonzero(events_tmp)])
                    M_global = np.max(events_tmp[np.nonzero(events_tmp)])
                    file_log.writelines(['m global {} \n'.format(m_global),'M global {} \n'.format(M_global)])
                    # Events without last time step of previous chunk                    
                    events = np.delete(events_tmp,0,0)                    
                    ###########################################################
                    ### Check if some events continue on the next chunk
                    ###########################################################
                    # Events on last time step of previous chunk
                    tmp_prev = events_tmp[0,:,:,:]
                    del events_tmp
                    # if there are some events in the last time step of previous chunk
                    if len(np.nonzero(tmp_prev)[0])>0:
                        # take the min and max to loop over them
                        m_local = np.min(tmp_prev[np.nonzero(tmp_prev)])
                        M = np.max(tmp_prev[np.nonzero(tmp_prev)])
                        file_log.writelines(['m local {} \n'.format(m_local),'M local {} \n'.format(M)])
                        # for each number see if the event continues
                        for j in range(m_local,M+1):
                            # check if j is in previous year events detection
                            if len(tmp_prev[tmp_prev ==j])>0:
                                # check if this event is also in the new chunk 
                                if len(events[0,:,:,:][events[0,:,:,:] == j])>0:
                                    # if yes, then we have a common event 
                                    # find the corresponding event id in previous chunk
                                    # CAUTION it can correspond to several id! 
                                    #take each of these id in evs
                                    evs = np.unique(events_prev[tmp_prev == j])
                                    # in dictionnaries, merge all the corresponding events to the first id 
                                    idx_ev = np.where(mhw['ev_number']== evs[0])
                                    
                                    ############################################################################
                                    # First gather previous events that are, in the end, only one event  
                                    ############################################################################
                                    # if is indeed in dictionnary and if there are at least two events corresponding
                                    if len(idx_ev[0])!=0 and len(evs)>1:
                                        # Merge the other events in this index and remove them from dictionaaries 
                                        for ev_dict in range(1,len(evs)):
                                            idx_tmp_ev = np.where(mhw['ev_number']== evs[ev_dict]) # evs[ev_dict] is the number of the events to be merged # mhw['ev_number'][idx_ev[0][0]] is the number of the event that will remain in dict
                                            
                                            # security if check
                                            if len(idx_tmp_ev[0])!=0:
                                                # merge extremes in coordinates dictionary
                                                extremes = merge_extremes(mhw,idx_ev,evs,ev_dict,extremes)
                                                # merge extremes in characteristics dictionary
                                                mhw = merge_mhws(mhw,idx_tmp_ev,idx_ev,additional_diags)

                                                # replace the number of the merged event by the number of the mother event (the one kept in dictionnary) to be able to forward future event linked to that event to the mother event 
                                                events_prev[events_prev == evs[ev_dict]] = evs[0] 
                                                # if addition had already occurred for this event in previous j values, replace the event number in event array by the mother event number to keep track for future links to the next chunck
                                                if len (events[events == evs[ev_dict]]) != 0:
                                                    print('Changing numbers because additon already occured on event {}'.format(evs[ev_dict]))
                                                    file_log.writelines(['Changing numbers because additon already occured on event {} \n'.format(evs[ev_dict])])
                                                    # when you merge the event, check if this event number was in event array (if an event of the new chunck had already been added to that event, needs to change its number to mother event number)
                                                    events[events == evs[ev_dict]] = evs[0]

                                            else:
                                                # should not happen, write a line in the log file
                                                file_log.writelines(['Event {} to gather in {} was removed from dict already (already merged) \n'.format(evs[ev_dict],evs[0])])

    
                                    ############################################################################
                                    # Second, add the new parts of the event from the new chunk 
                                    ############################################################################   

                                    # security if check  
                                    if len(idx_ev[0])!=0:
                                        ev = evs[0]
                                        # replace the number of the event in the new chunk to not process it again later and keep track for next chunk 
                                        temp = np.where(events == j) 
                                        events[events == j] =ev
                                        # Append to coordinates dictionnary to the existing key ev
                                        extremes = add_to_extremes(ev,temp,t_ind,extremes)
                                        # Add in characteristics dictionnary
                                        print(j,ev,min(temp[0]),max(temp[0]))
                                        if not moving_baseline:
                                            mhw = add_to_mhws(mhw,additional_diags,temp,idx_ev,area,dz,zlevs_rho,I,t,t_ind,detect_var_curr,omega_arag,pH,T,O2,thresh,thresh_omega=thresh_omega,thresh_pH=thresh_pH,thresh_O2=thresh_O2)
                                        else:
                                            # Load the 4D moving baseline threhsold
                                            thresh_omega_tmp = load_threshold_linear_trend_removed(thresh_omega,slope_omega,range(t_ind,t_ind+omega_arag.shape[0])) 
                                            thresh_pH_tmp = load_threshold_linear_trend_removed(thresh_pH,slope_pH,range(t_ind,t_ind+pH.shape[0])) 
                                            mhw = add_to_mhws(mhw,additional_diags,temp,idx_ev,area,dz,zlevs_rho,I,t,t_ind,detect_var_curr,omega_arag,pH,T,O2,thresh,thresh_omega=thresh_omega_tmp,thresh_pH=thresh_pH_tmp,thresh_O2=thresh_O2)    
                                            del thresh_omega_tmp,thresh_pH_tmp 

                                    else:
                                        # should not happen, write a line in the log file
                                        file_log.writelines(['Event {} in which we want to add new part of the event was removed from dict already (already merged) \n'.format(evs[0]),'j={} and evs are {}: and idx_ev is {} \n'.format(j,evs,idx_ev)])                                        

                    ##############################################################################
                    # Once continuous events have been taken care of, deal with all the new events
                    ##############################################################################
                   
                    # Loop over events
                    for ev in range(m_global,M_global+1):
                        temp = np.where(events == ev)
                        # if it was a continuous event, its number will not be between m_global and M_global+1
                        # so the counter keys will be less than minDuration so not processed
                        if len(Counter(temp[0]).keys())<minDuration and max(temp[0], default=0)!=I.shape[0]-1:
                            # Remove this event from events
                            events[events == ev] = 0   
                            continue
                        
                        # Fill the coordinates and characteristics dictionnaries with the new event
                        extremes,mhw,v_gridded,v = fill_extreme_dictionaries_with_core_chars(mhw,extremes,ev,temp,I,t_ind,t,area,dz,zlevs_rho,detect_var_curr)             
                        ### Other metrics
                        if additional_diags:
                            if not moving_baseline:
                                mhw = fill_extreme_dictionaries_with_add_chars(mhw,temp,omega_arag,pH,O2,T,v_gridded,v,thresh,thresh_omega=thresh_omega,thresh_pH=thresh_pH,thresh_O2=thresh_O2)
                            else:
                                # Load the 4D moving baseline threhsold
                                thresh_omega_tmp = load_threshold_linear_trend_removed(thresh_omega,slope_omega,range(t_ind,t_ind+omega_arag.shape[0])) 
                                thresh_pH_tmp = load_threshold_linear_trend_removed(thresh_pH,slope_pH,range(t_ind,t_ind+pH.shape[0])) 
                                mhw = fill_extreme_dictionaries_with_add_chars(mhw,temp,omega_arag,pH,O2,T,v_gridded,v,thresh,thresh_omega=thresh_omega_tmp,thresh_pH=thresh_pH_tmp,thresh_O2=thresh_O2)  
                                del thresh_omega_tmp,thresh_pH_tmp                                             
                
                # If no events in this chunk, keep events array anyway for the script to continue
                else:
                    events = np.delete(events_tmp,0,0)
                # Keep track of the global time indice over the entire detection 
                t_ind += T_file
            
        ###########################################################################
        ### Save the data from this year, to be able to perform a restart if needed
        ###########################################################################        
        mhw_prev = mhw.copy() # save the dictionnary before processing this year if restart is needed
        extremes_prev = extremes.copy()
        events_prev = events[-1,:,:,:]
        n_events_prev = mhw['n_events']
        ## save those + t_ind, yr etc for restart
        pickle.dump(mhw_prev,open(frestartname_charac,'wb'))
        pickle.dump(extremes_prev,open(frestartname_coords,'wb'),protocol=pickle.HIGHEST_PROTOCOL)
        np.save(outpath + 'events_prev_restart_file.npy',events_prev)
        np.save(outpath + 'detect_var_prev_restart_file.npy',detect_var_prev)
        if two_threshs:
            np.save(outpath + 'detect_var2_prev_restart_file.npy',detect_var2_prev)
        np.save(outpath + 't_ind_yr_n_events_prev_restart_file.npy',np.asarray([t_ind,yr,n_events_prev]))
        print('....  data from previous year saved for restart')   
        # Go to the next year
        yr += 1
    
    ###########################################################################
    ### Global processing prior to finishing 
    ###########################################################################
    # remove all events shorter than minDuration from the dictionnaries (ie. events that did not continue on the next chunk)
    # Remove from coordinates dictionnay 
    for k in list(extremes):
        if np.asarray(mhw['duration'])[np.asarray(mhw['ev_number'])==int(k)]<minDuration:
            del extremes[k]
    # Remove from characteristics dictionnay                   
    for k in mhw.keys():
        if k == 'duration' or k == 'n_events':
            continue
        track = 0
        for loop in range(len(mhw['duration'])):
            if mhw['duration'][loop]<minDuration:
                mhw[k].pop(track)
            else:
                track += 1       

    mhw['duration'] = np.asarray(mhw['duration'])[np.asarray(mhw['duration'])>=minDuration]
                
    # Replace the number of events by the length of the dictionnary to have the number of events fulfilling conditions
    mhw['n_events'] = len(mhw['ev_number'])  
        
    ###########################################################################
    ####################         Saving output       ##########################
    ########################################################################### 
    file_log.close() 
    # as dictionnary files
    mhw_detected = mhw.copy()
    extreme_detected_loc = extremes.copy()
    # save characteristics
    fname = outpath + 'extreme_detection_{}m_{}_{}_minD_{}_3.pck'.format(depth,years[0],years[-1],minDuration)
    pickle.dump(mhw_detected,open(fname,'wb'))   
    # convert the dictionary into a DataFrame
    data = pd.DataFrame.from_dict(mhw)
    # write the dataFrame into csv file
    data.to_csv(outpath + 'extreme_detection_{}m_{}_{}_minD_{}_3.csv'.format(depth,years[0],years[-1],minDuration))
    
    # save coordinates
    fname = outpath + 'extreme_detection_loc_{}m_{}_{}_minD_{}_3.pck'.format(depth,years[0],years[-1],minDuration)
    pickle.dump(extreme_detected_loc,open(fname,'wb'),protocol=pickle.HIGHEST_PROTOCOL)   

    print('Detection done')
                                      

###################
## Launch script ##
###################
tracex()
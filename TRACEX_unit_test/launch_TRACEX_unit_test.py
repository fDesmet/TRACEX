#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Creation: May 2021
Authors: Flora Desmet, ETH Zuerich

Description:
Main script for the detection of extremes in a small 2D array as an example for the package.

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
import xarray as xr
## TRACEX functions 
from TRACEX_functions_unit_test import *

#######################################################
## Load input (normally in define_TRACEX_variables.py) 
#######################################################

# Chunk length is 4 days
chunklen = int(4)
numchunk = int(np.ceil(8./chunklen)+1)
# One variable only but could do the testing with two booleans for two thresholds as well
two_threshs = False 
var = 'extremes'
# No additional diagnostics wanted
additional_diags = False
outpath  = '/Users/floradesmet/Documents/ETHz_flora/TRACEX/TRACEX_unit_test/output_TRACEX_unit_test/'
t = np.arange(20,28,1)
# Not a restart. Restart would not really work in this example because only one single file is used. 
restart = 0
# Select edge or corners connection
edges_or_corners = 1
edges_only = np.abs(1-edges_or_corners)
# Select fixed baseline
moving_baseline = False

#####################################
## definition of the core rountine ##
#####################################
def tracex():

    #### Define restart filenames where dictionnaries will be saved at the end of each year to be able to restart the script at any year within the study period
    frestartname_charac = outpath + 'extreme_detection_restart_file.pck'
    frestartname_coords= outpath + 'extreme_detection_extremes_dict_restart_file.pck'

    #### Load filelist for detection and data of restart if it is a restart
    if not restart:
        filelist = ['/Users/floradesmet/Documents/ETHz_flora/TRACEX/TRACEX_unit_test/boolean_array_for_TRACEX_unit_test.nc']
        t_ind=0       
    else:
        num = np.load(outpath + 't_ind_yr_n_events_prev_restart_file.npy') 
        t_ind = num[0]
        filelist = ['/Users/floradesmet/Documents/ETHz_flora/TRACEX/TRACEX_unit_test/boolean_array_for_TRACEX_unit_test.nc']
    
    ### Define properties of event (keys of dictionary)
    if not restart:
        # create the dictionary with the global properties
        mhw = init_property_extreme_dict(additional_diags)
            
        # Create the dictionnary to keep gridcells of each events
        extremes = init_extremes_dict()
                
        ## Write detection steps into a text file   
        file_log = open(outpath+'detection_log.txt','w')

    else:
        ## Write detectionsteps into a text file   
        file_log = open(outpath+'detection_log.txt','a')


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

        ##############################################################################################
        ###     loop over the chunks to speed up the script  (by avoiding to load everything at once)
        ##############################################################################################
        for m in range(1,numchunk):
            file_log.writelines(['Processing m = {}'.format(m)])
            #if we are processing the first file not in restart mode   
            if not restart and f==filelist[0] and m==1:          
                nc_in = netCDF4.Dataset(f,'r')
                exceed_bool,T_file = load_var_datachunks_optimized(nc_in,var,m,chunklen,numchunk,additional_diags)
                if two_threshs:
                    exceed_bool2 = load_var_2_datachunks_optimized(nc_in,var2,mask3D,m,chunklen,numchunk,indices)
                print(T_file)
                print(t_ind)
                if additional_diags:
                    if len(nc_in.variables['temp'][:,0,0,0]) == 365:
                        temp_clm_ident=netCDF4.Dataset(clm_file_output, 'r')
                        T = T - load_datachunk(temp_clm_ident,'temp_smooth',m,chunklen,numchunk,indices) 
                    else:
                        T = T - load_npy_file(path_in,filename_clm_366j,m,chunklen,numchunk,indices)                  
                                   
                exceed_bool[exceed_bool<=0] = False
                exceed_bool[exceed_bool>0] = True
                exceed_bool[np.isnan(exceed_bool)]=False     # Discard land mass points
                # if second threshold (compound extremes), process the second exceed_bool array and multiply the two
                if two_threshs:
                    exceed_bool2[exceed_bool2<=0] = False
                    exceed_bool2[exceed_bool2>0] = True
                    exceed_bool2[np.isnan(exceed_bool2)]=False     # Discard land mass pointspoints
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
                detect_var_prev = np.reshape(exceed_bool[-1,:],(1,exceed_bool.shape[1]))
                if two_threshs:
                    detect_var2_prev = np.reshape(detect_var2[-1,:,:,:],(1,detect_var2.shape[1],detect_var2.shape[2],detect_var2.shape[3]))
                
                mhw['n_events'] = n_events
                print('There are {} events in chunk {} \n'.format(n_events,m))
                file_log.writelines(['There are {} events in chunk {} \n'.format(n_events,m)])

                # Only for the step-by-step explanation: save the array numbers after first step 
                ds = xr.Dataset(data_vars=dict(events=(["time", "x"], events)),
                                coords=dict(time=np.arange(0,events.shape[0],1),x=np.arange(0,events.shape[1],1)),
                                attrs=dict(description="Event number in chunk 1"))
                ds.to_netcdf(outpath+'events_chunk_1_for_TRACEX_unit_test.nc')

                ##################
                # Loop over events
                ##################
                for ev in range(1,n_events+1):  
                    temp = np.where(events == ev)
                    # temp is a tuple of 4 arrays of int with the index of the location 
                    # of the events: temp[0] is an array of all the time index where ev is present

                    # no minimum duration requirement in this exmple 

                    ##################################################
                    # Fill the dictionnary with the event
                    ##################################################
                    extremes,mhw,v_gridded,v = fill_extreme_dictionaries_with_core_chars(mhw,extremes,ev,temp,t_ind,t,exceed_bool)
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
                #if we are processing the first file in restart mode   
                if restart and f==filelist[0] and m==1: 
                    print('Data already loaded from restart...')
                else:
                    # Save last time step of previous events if not restart time
                    events_prev = events[-1,:]
                    #save the previous total number of events
                    n_events_prev = mhw['n_events']    
                    
                #open new file 
                nc_in_curr = netCDF4.Dataset(f,'r')                
                exceed_bool_curr,T_file = load_var_datachunks_optimized(nc_in_curr,var,m,chunklen,numchunk,additional_diags)
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
                       
                exceed_bool = np.concatenate((detect_var_prev,exceed_bool_curr),axis=0)

                exceed_bool[exceed_bool<=0] = False
                exceed_bool[exceed_bool>0] = True
                exceed_bool[np.isnan(exceed_bool)]=False     # Discard land mass points
                # if second threshold (compound extremes), process the second exceed_bool array and multiply the two
                if two_threshs:
                    exceed_bool2[exceed_bool2<=0] = False
                    exceed_bool2[exceed_bool2>0] = True
                    exceed_bool2[np.isnan(exceed_bool2)]=False     # Discard land mass points points
                    # Multiply the two boolean arrays
                    exceed_bool = np.multiply(exceed_bool,exceed_bool2)
                    
                # Find contiguous regions of exceed_bool = True
                if edges_only:    
                    events_tmp, n_events = ndimage.label(exceed_bool)
                elif edges_or_corners:
                    dims = len(np.shape(exceed_bool))
                    s = ndimage.generate_binary_structure(dims,dims)  # Make diagonal connections
                    events_tmp, n_events = ndimage.label(exceed_bool,s)                


                # Only for the step-by-step explanation: save the array numbers after detection in second chunk
                ds = xr.Dataset(data_vars=dict(events=(["time", "x"], events_tmp)),
                                coords=dict(time=np.arange(0,events_tmp.shape[0],1),x=np.arange(0,events_tmp.shape[1],1)),
                                attrs=dict(description="Event number in chunk 2 raw"))
                ds.to_netcdf(outpath+'events_chunk_2_raw_TRACEX_unit_test.nc')

                events_tmp[events_tmp != 0] += n_events_prev

                # Only for the step-by-step explanation: save the array numbers after adding the maximum ID number of first chunk to the second chunck detection
                ds = xr.Dataset(data_vars=dict(events=(["time", "x"], events_tmp)),
                                coords=dict(time=np.arange(0,events_tmp.shape[0],1),x=np.arange(0,events_tmp.shape[1],1)),
                                attrs=dict(description="Event number in chunk 2 added numbers"))
                ds.to_netcdf(outpath+'events_chunk_2_added_numbers_for_TRACEX_unit_test.nc')
                
                # Remove unecessary variables
                detect_var_prev = np.reshape(exceed_bool[-1,:],(1,exceed_bool.shape[1]))
                if two_threshs:
                    detect_var2_prev = np.reshape(detect_var2[-1,:,:,:],(1,detect_var2.shape[1],detect_var2.shape[2],detect_var2.shape[3]))
                    del detect_var2                    
                
                mhw['n_events'] += n_events
                print('There are {} events in chunk {}  \n'.format(n_events,m))
                file_log.writelines(['There are {} events in chunk {} \n'.format(n_events,m,)])

                loop_log = 1
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
                    #event without previous chunk                    
                    events = np.delete(events_tmp,0,0)                    
                    ###########################################################
                    ### Check if some events continue on the next chunk
                    ###########################################################
                    # events on last time of previous chunk
                    tmp_prev = events_tmp[0,:]
                    del events_tmp
                    # if there are some events in last time of previous chunk
                    if len(np.nonzero(tmp_prev)[0])>0:
                        # take the min and max to loop over them
                        m_local = np.min(tmp_prev[np.nonzero(tmp_prev)])
                        M = np.max(tmp_prev[np.nonzero(tmp_prev)])
                        file_log.writelines(['m local {} \n'.format(m_local),'M local {} \n'.format(M)])
                        #for each number see if the event continues
                        for j in range(m_local,M+1):
                            #check if j is in previous year events detection
                            if len(tmp_prev[tmp_prev ==j])>0:
                                # check if this event is also in the new chunk 
                                if len(events[0,:][events[0,:] == j])>0:
                                    # if yes, then we have a common event 
                                    #find the corresponding event id in previous chunk
                                    # CAUTION it can correspond to several id! 
                                    #take each of these id in evs
                                    evs = np.unique(events_prev[tmp_prev == j])
                                    # in dictionnary we gather all the corresponding events to the first id 
                                    idx_ev = np.where(mhw['ev_number']== evs[0])
                                    
                                    ############################################################################
                                    # First gather previous events that are, in the end, only one event  
                                    ############################################################################
                                    # if is indeed in dictionnary and if there are at least two events corresponding
                                    if len(idx_ev[0])!=0 and len(evs)>1:
                                        #gather the other events from dict in this index and remove it from dict 
                                        for ev_dict in range(1,len(evs)):
                                            idx_tmp_ev = np.where(mhw['ev_number']== evs[ev_dict]) # evs[ev_dict] is the number of the events to be merged # mhw['ev_number'][idx_ev[0][0]] is the number of the event that will remain in dict
                                            
                                            # security if check
                                            if len(idx_tmp_ev[0])!=0:
                                                # merge extremes in coordinates dictionary
                                                extremes = merge_extremes(mhw,idx_ev,evs,ev_dict,extremes)
                                                # merge extremes in characteristics dictionary
                                                mhw = merge_mhws(mhw,idx_tmp_ev,idx_ev,additional_diags)

                                                # replace the number of the merged event by the number of the mother event (so one kept in dictionnary) to be able to forward future event linked to that event to the mother event 
                                                events_prev[events_prev == evs[ev_dict]] = evs[0] 
                                                # Only for the step-by-step explanation
                                                ds = xr.Dataset(data_vars=dict(events=(["x"], events_prev)),
                                                                coords=dict(x=np.arange(0,len(events_prev),1)),
                                                                attrs=dict(description="Event number in last time of chunk 1 "))
                                                ds.to_netcdf(outpath+'events_chunk_2_loop_log_{}_for_TRACEX_unit_test.nc'.format(loop_log))
                                                loop_log += 1

                                                # if addition had already occurred for this event in previous j values, replace the event number in event array by the mother event number to keep track for future links to the next chunck
                                                if len (events[events == evs[ev_dict]]) != 0:
                                                    print('Changing numbers because additon already occured on event {}'.format(evs[ev_dict]))
                                                    file_log.writelines(['Changing numbers because additon already occured on event {} \n'.format(evs[ev_dict])])
                                                    # when you merge the event, check if this event number was in event array (if an event of the new chunck had already been added to that event, needs to change its number to mother event number)
                                                    events[events == evs[ev_dict]] = evs[0]
                                                    # Only for the step-by-step explanation
                                                    ds = xr.Dataset(data_vars=dict(events=(["time", "x"], events)),
                                                                coords=dict(time=np.arange(0,events.shape[0],1),x=np.arange(0,events.shape[1],1)),
                                                                attrs=dict(description="Event number in chunk 2 "))
                                                    ds.to_netcdf(outpath+'events_chunk_2_loop_log_{}_for_TRACEX_unit_test.nc'.format(loop_log))
                                                    loop_log += 1


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
                                        # Only for the step-by-step explanation
                                        ds = xr.Dataset(data_vars=dict(events=(["time", "x"], events)),
                                                                coords=dict(time=np.arange(0,events.shape[0],1),x=np.arange(0,events.shape[1],1)),
                                                                attrs=dict(description="Event number in chunk 2 "))
                                        ds.to_netcdf(outpath+'events_chunk_2_loop_log_{}_for_TRACEX_unit_test.nc'.format(loop_log))
                                        loop_log += 1

                                        # Append to coordinates dictionnary to the existing key ev
                                        extremes = add_to_extremes(ev,temp,t_ind,extremes)
                                        # Add in characteristics dictionnary
                                        print(j,ev,min(temp[0]),max(temp[0]))
                                        if not moving_baseline:
                                            mhw = add_to_mhws(mhw,additional_diags,temp,idx_ev,t,t_ind,exceed_bool_curr)
                                        else:
                                            # Load the 4D moving baseline threhsold
                                            thresh_omega_tmp = load_threshold_linear_trend_removed(thresh_omega,slope_omega,range(t_ind,t_ind+omega_arag.shape[0])) 
                                            thresh_pH_tmp = load_threshold_linear_trend_removed(thresh_pH,slope_pH,range(t_ind,t_ind+pH.shape[0])) 
                                            mhw = add_to_mhws(mhw,additional_diags,temp,idx_ev,area,dz,zlevs_rho,I,t,t_ind,detect_var_curr,omega_arag,pH,T,O2,thresh,thresh_omega=thresh_omega_tmp,thresh_pH=thresh_pH_tmp,thresh_O2=thresh_O2)    
                                            del thresh_omega_tmp,thresh_pH_tmp
                                               

                                    else:
                                        # should not happen, write a line in the log file
                                        file_log.writelines(['Event {} in which we want to add new part of the event was removed from dict already (already merged) \n'.format(evs[0]),'j={} and evs are {}: and idx_ev is {} \n'.format(j,evs,idx_ev)])                                        

                    ###############################################################################
                    # Once continuous events have been taken care of, deal with all the new events
                    ###############################################################################
                   
                    # Loop over events
                    for ev in range(m_global,M_global+1):
                        temp = np.where(events == ev)
                        # if it was a continuous event, its number will not be between m_global and M_global+1
                        # so the counter keys will be less than minDuration so not processed
                        if len(Counter(temp[0]).keys())<1:
                            # Remove this event from events
                            events[events == ev] = 0   
                            continue
                        
                        # Fill the coordinates and characteristics dictionnaries with the new event
                        extremes,mhw,v_gridded,v = fill_extreme_dictionaries_with_core_chars(mhw,extremes,ev,temp,t_ind,t,exceed_bool_curr)             
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
        events_prev = events[-1,:]
        n_events_prev = mhw['n_events']
        ## save those + t_ind, yr etc 
        pickle.dump(mhw_prev,open(frestartname_charac,'wb'))
        pickle.dump(extremes_prev,open(frestartname_coords,'wb'),protocol=pickle.HIGHEST_PROTOCOL)
        np.save(outpath + 'events_prev_restart_file.npy',events_prev)
        np.save(outpath + 'detect_var_prev_restart_file.npy',np.ma.filled(detect_var_prev))
        if two_threshs:
            np.save(outpath + 'detect_var2_prev_restart_file.npy',detect_var2_prev)
        np.save(outpath + 't_ind_n_events_prev_restart_file.npy',np.asarray([t_ind,n_events_prev]))
        print('....  data from previous year saved for restart')   

    
    ###########################################################################
    ### Global processing prior to finishing 
    ###########################################################################
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
    fname = outpath + 'extreme_detection.pck'
    pickle.dump(mhw_detected,open(fname,'wb'))   
    # save coordinates
    fname = outpath + 'extreme_detection_loc.pck'
    pickle.dump(extreme_detected_loc,open(fname,'wb'),protocol=pickle.HIGHEST_PROTOCOL)   
   
    print('Detection done')                      


###################
## Launch script ##
###################
tracex()




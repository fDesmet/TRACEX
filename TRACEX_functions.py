#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Creation: May 2021
Authors: Flora Desmet and Eike Koehn, ETH Zuerich

Description:
Define functions needed for running launch_TRACEX.py 

"""
import numpy as np
import netCDF4

###### for characteristics dictionnary

def init_property_extreme_dict(additional_diags):
    """
    The properties stored to characterize the extremes are the following:
    
    ev_number: number of the event (tag common to both the characteristics and coordinates dictionnaries)

    ## Time properties:
    time_I_max    (datetime format) time of maximum intensity
    time_z_max    (datetime format) time when the event reaches its shallowest depth
    index_start   (ordinal date format) ordinal time when the event initiates
    index_end     (ordinal date format) ordinal time when the event ends
    duration      (index_end - index_start) [days]

    ## Magnitude properties:
    intensity_max                 maximum [Thresh-Value] value (if extremes are defined as below the x percentile) 
    intensityxvolume_detect_var   severity [km3.d.var unit]: cumulative intensity in space and time
    intensity_mean                cumulative intensity (severity) divided by the cumulative volume over time
    detect_var_min                minimum value of the variable used for detection reached within the extreme (relevant if extremes are defined as below the x percentile)
    detect_varxvolume_detect_var  volume weighted space and time sum of values in the event of the variable used for detection 
    detect_var_mean               time and volume weighted mean value of the variable used for detection in the event (detect_varxvolume_detect_var/cumulative volume in time)
        
 
    !!!! Note that for all depth-related properties computation, the ROMS grid has the shallowest depth associated to the highest s-level and the deepest depth to the 0 s-level !!!!!

    ## Size properties
    volume_tot          [km3.d] cumulative volume over time
    volume_mean         [km3] daily volume: cumulative volume / duration 
    vertical_extent     [m] difference between shallowest and deepest point

    ## Loaction properties
    z_max               [m] shallowest depth reached by the event 
    z_min_at_z_max      [m] deepest depth where the event reached the shallowest level
    loc_z_max_eta       ETA index (on reduced grid) of the shallowest point
    loc_z_max_xi        XI index (on reduced grid) of the shallowest point
    loc_Imax_eta        ETA index (on reduced grid) of the maximum intensity point
    loc_Imax_xi         XI index (on reduced grid) of the maximum intensity point

    ## Magnitude properties for addtional diagnostics (needed if compound extremes to get magnitude in the secondary variable used in the detection)  

    For now, implemented for omega aragonite, pH, oxygen:

    omega_min               minimum value of omega aragonite 
    omegaxvolume_omega      volume weighted space and time sum of omega aragonite in the event 
    omega_mean              time and volume weighted mean omega aragonite in the event (omegaxvolume_omega/cumulative volume in time)
    intensityxvolume_omega  severity [km3.d.var unit]: cumulative intensity for omega aragonite in space and time if a threshold is given; otherwise, equal to omegaxvolume_omega
    omega_I_mean            severity divided by the cumulative volume over time if a threshold is given; otherwise, equal to omega_mean
    omega_I_max             maximum [Thresh-Value] value (if extremes are defined as below the x percentile) if a threshold is given; otherwise only the maximum value of omega aragonite in the event

    For temperature: 

    delta_temp_max                  [degC] the maximum temperature anomaly (positive or negative) 
    intensityxvolume_delta_temp     [km3.d.degC] volume weighted space and time sum of temperature anomalies in the event 
    delta_temp_mean                 [degC] time and volume weighted mean temperature anomaly in the event (intensityxvolume_delta_temp/cumulative volume in time)

    """


    mhw = {}
    mhw['ev_number'] = []

    ## Time properties:
    mhw['time_I_max'] = []  # datetime format: time of maximum intensity
    mhw['time_z_max'] = []  # datetime format: time when the event reaches its shallowest depth
    mhw['index_start'] = [] # ordinal date format: ordinal time when the event initiates
    mhw['index_end'] = []   # ordinal date format: ordinal time when the event ends
    mhw['duration'] = []    # index_end - index_start [days]

    ## Magnitude properties:
    mhw['intensity_max'] = []               # maximum [Thresh-Value] value (if extremes are defined as below the x percentile)
    mhw['intensityxvolume_detect_var'] = [] # severity [km3.d.var unit]: cumulative intensity in space and time
    mhw['intensity_mean'] = []              # cumulative intensity divided by the cumulative volume over time
    mhw['detect_var_min'] = []              # minimum value of the variable used for detection reached within the extreme (relevant if extremes are defined as below the x percentile; should be turned into maximum otherwise, ie for temperature extremes)
    mhw['detect_varxvolume_detect_var'] =[] # volume weighted space and time sum of values in the event of the variable used for detection 
    mhw['detect_var_mean'] = []             # time and volume weighted mean value of the variable used for detection in the event (detect_varxvolume_detect_var/cumulative volume in time)

    ## Size properties:
    mhw['volume_tot'] = []      # [km3.d] cumulative volume over time
    mhw['volume_mean'] = []     # [km3] daily volume: cumulative volume / duration 
    mhw['vertical_extent'] =[]  # [m] difference between shallowest and deepest point 

    #  !!!! Note that for all depth-related properties computation, the ROMS grid has the shallowest depth associated to the highest s-level and the deepest depth to the 0 s-level !!!!! #

    # Loaction properties:
    mhw['z_max'] = []           # [m] shallowest depth reached by the event 
    mhw['z_min_at_z_max'] = []  # [m] deepest depth where the event reached the shallowest level
    mhw['loc_z_max_eta'] = []   # ETA index (on reduced grid) of the shallowest point
    mhw['loc_z_max_xi'] = []    # XI index (on reduced grid) of the shallowest point
    mhw['loc_Imax_eta'] = []    # ETA index (on reduced grid) of the maximum intensity point
    mhw['loc_Imax_xi'] = []     # XI index (on reduced grid) of the maximum intensity point

    if additional_diags:  
        ## Magnitude properties for addtional diagnostics (needed if compound extremes to get magnitude in the secondary variable used in the detection)  
        mhw['omega_min'] = []               #  minimum value of omega aragonite
        mhw['omegaxvolume_omega'] = []      # volume weighted space and time sum of omega aragonite in the event 
        mhw['omega_mean'] = []              # time and volume weighted mean omega aragonite in the event (omegaxvolume_omega/cumulative volume in time)
        mhw['intensityxvolume_omega'] = []  # severity [km3.d.omega aragonite unit]: cumulative intensity for omega aragonite in space and time if a threshold is given; otherwise, equal to omegaxvolume_omega
        mhw['omega_I_mean'] = []            # severity divided by the cumulative volume over time if a threshold is given; otherwise, equal to omega_mean
        mhw['omega_I_max'] = []             # maximum [Thresh-Value] value (if extremes are defined as below the x percentile) if a threshold is given; otherwise only the maximum value of omega aragonite in the event

        mhw['pH_min'] = [] 
        mhw['pHxvolume_pH'] = []
        mhw['pH_mean'] = [] 
        mhw['intensityxvolume_pH'] = [] 
        mhw['pH_I_mean'] = [] 
        mhw['pH_I_max'] = [] 

        mhw['oxygen_min'] = [] 
        mhw['oxygenxvolume_oxygen'] = []
        mhw['oxygen_mean'] = []
        mhw['intensityxvolume_oxygen'] = []      
        mhw['oxygen_I_mean'] = []    
        
        mhw['delta_temp_max'] = []              # [degC] the maximum temperature anomaly (positive or negative)
        mhw['intensityxvolume_delta_temp'] = [] # [km3.d.degC] volume weighted space and time sum of temperature anomalies in the event 
        mhw['delta_temp_mean'] = []             # [degC] time and volume weighted mean temperature anomaly in the event (intensityxvolume_delta_temp/cumulative volume in time)

    return mhw


###### for coordinates dictionnary

def init_extremes_dict():
    extremes = {}
    return extremes


def init_indiv_extreme_dict(extremes,number):
    extremes[str(number)] = dict.fromkeys(['time', 's_rho', 'eta_rho', 'xi_rho'])
    return extremes


def patch_to_dict(number,time,s_rho,eta_rho,xi_rho,extremes):
    # store it as a long 1D array instead of list 
    # check if it is a new event or if we are merging
    if str(number) not in extremes.keys():
        init_indiv_extreme_dict(extremes,number)
        extremes[str(number)]['time'] = np.asarray(time,dtype='int16')
        extremes[str(number)]['s_rho'] = np.asarray(s_rho,dtype='int8')
        extremes[str(number)]['eta_rho'] = np.asarray(eta_rho,dtype='int16')
        extremes[str(number)]['xi_rho'] = np.asarray(xi_rho,dtype='int16')
    else:
        extremes[str(number)]['time'] = np.hstack((extremes[str(number)]['time'],np.asarray(time,dtype='int16')))
        extremes[str(number)]['s_rho'] = np.hstack((extremes[str(number)]['s_rho'],np.asarray(s_rho,dtype='int8')))
        extremes[str(number)]['eta_rho'] = np.hstack((extremes[str(number)]['eta_rho'],np.asarray(eta_rho,dtype='int16')))
        extremes[str(number)]['xi_rho'] = np.hstack((extremes[str(number)]['xi_rho'],np.asarray(xi_rho,dtype='int16')))   

    return extremes



def load_datachunk(ncf,varname,m,chunklen,numchunk,indices): 
    idx_eta_min = indices[0]
    idx_eta_max = indices[1]
    idx_xi_min = indices[2]
    idx_xi_max = indices[3]
    if m == (numchunk-1):
        datachunk = ncf.variables[varname][chunklen*(m-1):,:,idx_eta_min:idx_eta_max,idx_xi_min:idx_xi_max]
    else:
        datachunk = ncf.variables[varname][chunklen*(m-1):chunklen*m,:,idx_eta_min:idx_eta_max,idx_xi_min:idx_xi_max]
        
    return datachunk


def load_var_datachunks_optimized(ncf,var,mask3D,m,chunklen,numchunk,indices,additional_diags):
    ''' 
    Load main variable for detection and 
    additional variables for diagnostics if additional_diags is set to true. 
    '''

    if additional_diags:
        omega = load_datachunk(ncf,'omega_arag_offl',m,chunklen,numchunk,indices)
        pH = load_datachunk(ncf,'pH_offl',m,chunklen,numchunk,indices)
        O2 = load_datachunk(ncf,'O2',m,chunklen,numchunk,indices)
        T = load_datachunk(ncf,'temp',m,chunklen,numchunk,indices)
        
        if var == 'omega':
            detect_var = omega
        elif var == 'pH':
            detect_var = pH
        elif var == 'O2':
            detect_var = O2
        elif var == 'temp':
            detect_var = T
            
        mask4D = np.repeat(mask3D[np.newaxis,: , :, :], detect_var.shape[0], axis=0)
        # Mask values outside the mask (the depth limit is included in the mask)
        detect_var = np.ma.masked_where(mask4D==False,detect_var)
        detect_var = np.ma.filled(detect_var,1e33) 
        T_file =  detect_var.shape[0]
    
    else:
        if var == 'omega':
            var_name = 'omega_arag_offl'
        elif var == 'pH':
            var_name = 'pH_offl'
        elif var == 'O2':
            var_name = 'O2'
        elif var == 'temp':
            var_name = 'temp'
            
        detect_var = load_datachunk(ncf,var_name,m,chunklen,numchunk,indices)
        mask4D = np.repeat(mask3D[np.newaxis,: , :, :], detect_var.shape[0], axis=0)
        # Mask values outside the mask (the depth limit is included in the mask)
        detect_var = np.ma.masked_where(mask4D==False,detect_var)
        detect_var = np.ma.filled(detect_var,1e33) 
        T_file =  detect_var.shape[0]
        
        omega = np.zeros_like(detect_var)
        pH = np.zeros_like(detect_var)
        O2 = np.zeros_like(detect_var)
        T = np.zeros_like(detect_var)
        
        print(var_name)
        print(detect_var.shape)
        
    return detect_var,mask4D,T_file,omega,pH,O2,T


def load_var_2_datachunks_optimized(ncf,var2,mask3D,m,chunklen,numchunk,indices):
    ''' 
    Load a second variable for detection on two variables simultaneously 
    (compound extremes)
    '''

    if var2 == 'omega':
        var_name = 'omega_arag_offl'
    elif var2 == 'pH':
        var_name = 'pH_offl'
    elif var2 == 'O2':
        var_name = 'O2'
    elif var2 == 'temp':
        var_name = 'temp'
        
    detect_var2 = load_datachunk(ncf,var_name,m,chunklen,numchunk,indices)
    mask4D = np.repeat(mask3D[np.newaxis,: , :, :], detect_var2.shape[0], axis=0)
    # Mask values outside the mask
    detect_var2 = np.ma.masked_where(mask4D==False,detect_var2)
    detect_var2 = np.ma.filled(detect_var2,1e33) 
       
    print(var_name)
    print(detect_var2.shape)
        
    return detect_var2


def load_npy_file(path_in,filename,m,chunklen,numchunk,indices):
    '''
    Load 4D field for numpy npy file for a particular chunk and grid reduction
    '''
    idx_eta_min = indices[0]
    idx_eta_max = indices[1]
    idx_xi_min = indices[2]
    idx_xi_max = indices[3]
    if m == (numchunk-1):
        variable = np.load(path_in+filename,allow_pickle=True)[chunklen*(m-1):,:,idx_eta_min:idx_eta_max,idx_xi_min:idx_xi_max]        
    else:
        variable = np.load(path_in+filename,allow_pickle=True)[chunklen*(m-1):chunklen*m,:,idx_eta_min:idx_eta_max,idx_xi_min:idx_xi_max]
    return variable


def load_npy_file_3D(path_in,filename,indices,pickle=False):
    '''
    Load 3D field for numpy npy file for a particular grid reduction
    '''
    idx_eta_min = indices[0]
    idx_eta_max = indices[1]
    idx_xi_min = indices[2]
    idx_xi_max = indices[3]
    if pickle==True:
        variable = np.load(path_in+filename,allow_pickle=True)[:,idx_eta_min:idx_eta_max,idx_xi_min:idx_xi_max]
    else:
        variable = np.load(path_in+filename)[:,idx_eta_min:idx_eta_max,idx_xi_min:idx_xi_max]        
    return variable



def fill_extreme_dictionaries_with_core_chars(mhw,extremes,ev,temp,I,t_ind,t,area,dz,zlevs_rho,detect_var):
    '''
    Fill all properties for a newly detected event
    '''
    #### Fill in the coordinates dictionnary with a new entry ev
    extremes = patch_to_dict(ev,temp[0]+t_ind,temp[1],temp[2],temp[3],extremes)

    #### Fill in the characteristics dictionnary with a new entry ev
    ## Event number
    mhw['ev_number'].append(ev)

    ## Time properties
    t_min = min(temp[0])
    t_max = max(temp[0])
    duration = t_max - t_min + 1 # [days] if you have daily outputs
    mhw['duration'].append(duration)
    mhw['index_start'].append(t[t_ind] + t_min)
    mhw['index_end'].append(t[t_ind] + t_max)    

    ## Size properties
    # area is 2D array so take 2 index arrays (eta,xi) ; dz is 3D so take 3 index arrays (depth, eta,xi)
    v_gridded = np.multiply(area[temp[2:4]],dz[temp[1:4]])*10**-9 # km3
    v = np.sum(v_gridded)  # km3.d (sum of volume over each time step)
    v_mean = v/duration 
    mhw['volume_tot'].append(v)
    mhw['volume_mean'].append(v_mean)
    
    ## Magnitude properties
    # Maximum intensity
    mhw['intensity_max'].append(max(I[temp[0:4]]))
    mhw['detect_var_min'].append(min(detect_var[temp[0:4]]))
    mhw['loc_Imax_eta'].append(temp[2][np.argmax(I[temp[0:4]])])
    mhw['loc_Imax_xi'].append(temp[3][np.argmax(I[temp[0:4]])])
    mhw['time_I_max'].append(t[t_ind] + temp[0][np.argmax(I[temp[0:4]])])
    # Mean intensity
    IxV = np.sum(np.multiply(I[temp[0:4]],v_gridded))
    mhw['intensityxvolume_detect_var'].append(IxV)
    mhw['intensity_mean'].append(IxV/v) 
    # Mean value of variable
    IxV = np.sum(np.multiply(detect_var[temp[0:4]],v_gridded))
    mhw['detect_varxvolume_detect_var'].append(IxV)
    mhw['detect_var_mean'].append(IxV/v)   
    
    ## Location and vertical extent properties
    l = np.argmin(zlevs_rho[temp[1:4]])
    l_m_zma =  np.argmax(zlevs_rho[temp[1],temp[2][l],temp[3][l]])
    lm = np.argmax(zlevs_rho[temp[1:4]])
    mhw['z_max'].append(zlevs_rho[temp[1][l],temp[2][l],temp[3][l]])
    mhw['z_min_at_z_max'].append(zlevs_rho[temp[1],temp[2][l],temp[3][l]][l_m_zma]) 
    mhw['z_min'].append(zlevs_rho[temp[1][lm],temp[2][lm],temp[3][lm]])
    mhw['loc_z_max_eta'].append(temp[2][l])
    mhw['loc_z_max_xi'].append(temp[3][l])   
    mhw['vertical_extent'].append(zlevs_rho[temp[1][lm],temp[2][lm],temp[3][lm]]-zlevs_rho[temp[1][l],temp[2][l],temp[3][l]])
    mhw['time_z_max'].append(t[t_ind] + temp[0][l])              
    
    return extremes,mhw,v_gridded,v
    


def fill_extreme_dictionaries_with_add_chars(mhw,temp,omega,pH,O2,T,v_gridded,v,thresh,thresh_omega=0,thresh_pH=0,thresh_O2=0):
    '''
    Fill properties for additional diagnostics for a newly detected event
    omega aragonite, pH, oxygen and temperature
    '''
    
    ################
    #### omega  ####
    ################    
    mhw['omega_min'].append(min(omega[temp[0:4]]))
    IxV = np.sum(np.multiply(omega[temp[0:4]],v_gridded))
    mhw['omegaxvolume_omega'].append(IxV)
    mhw['omega_mean'].append(IxV/v)
    if len(np.shape(thresh_omega))==0:  
        # absolute threshold or none (float, not array)
        thresh_4D_omega = np.zeros_like(omega)+thresh_omega
    elif thresh_omega.ndim == 3:  
        # thresh is a 3D array (fixed baseline)
        thresh_4D_omega = np.repeat(thresh_omega[np.newaxis,: , :, :], omega.shape[0], axis=0)
    else: 
        # thresh is a 4D array (moving baseline)
        thresh_4D_omega = thresh_omega[:]
    IxV = np.sum(np.multiply(thresh_4D_omega[temp[0:4]]-omega[temp[0:4]],v_gridded))
    mhw['intensityxvolume_omega'].append(IxV)
    mhw['omega_I_mean'].append(IxV/v)    
    mhw['omega_I_max'].append(max(thresh_4D_omega[temp[0:4]]-omega[temp[0:4]])) 
    del thresh_4D_omega
    
    ################
    ####   pH   ####
    ################
    mhw['pH_min'].append(min(pH[temp[0:4]]))
    IxV = np.sum(np.multiply(pH[temp[0:4]],v_gridded))
    mhw['pHxvolume_pH'].append(IxV)
    mhw['pH_mean'].append(IxV/v)
    if len(np.shape(thresh_pH))==0:  
        # absolute threshold or none (float, not array)
        thresh_4D_pH = np.zeros_like(pH)+thresh_pH
    elif thresh_pH.ndim == 3:  
        # thresh is a 3D array (fixed baseline)
        thresh_4D_pH = np.repeat(thresh_pH[np.newaxis,: , :, :], pH.shape[0], axis=0)
    else: 
        # thresh is a 4D array (moving baseline)
        thresh_4D_pH = thresh_pH[:]
    IxV = np.sum(np.multiply(thresh_4D_pH[temp[0:4]]-pH[temp[0:4]],v_gridded))
    mhw['intensityxvolume_pH'].append(IxV)
    mhw['pH_I_mean'].append(IxV/v)
    mhw['pH_I_max'].append(max(thresh_4D_pH[temp[0:4]]-pH[temp[0:4]])) 
    del thresh_4D_pH
    
    ################
    ####   O2   ####
    ################    
    mhw['oxygen_min'].append(min(O2[temp[0:4]]))
    IxV = np.sum(np.multiply(O2[temp[0:4]],v_gridded))
    mhw['oxygenxvolume_oxygen'].append(IxV)    
    mhw['oxygen_mean'].append(IxV/v)
    if len(np.shape(thresh_O2))==0:  
        # absolute threshold or none (float, not array)
        thresh_4D_O2 = np.zeros_like(O2)+thresh_O2
    else:  
        # thresh is a 3D array (fixed baseline)
        thresh_4D_O2 = np.repeat(thresh_O2[np.newaxis,: , :, :], O2.shape[0], axis=0)
    IxV = np.sum(np.multiply(thresh_4D_O2[temp[0:4]]-O2[temp[0:4]],v_gridded))
    del thresh_4D_O2
    mhw['intensityxvolume_oxygen'].append(IxV)
    mhw['oxygen_I_mean'].append(IxV/v)
    
   
    ################
    ######  T  #####
    ################    
    T_ex = [min(T[temp[0:4]]), max(T[temp[0:4]])]
    mhw['delta_temp_max'].append(T_ex[np.argmax([abs(np.asarray(T_ex))])])
    IxV = np.sum(np.multiply(T[temp[0:4]],v_gridded))
    mhw['intensityxvolume_delta_temp'].append(IxV)
    mhw['delta_temp_mean'].append(IxV/v)
    
    return mhw
    

def merge_extremes(mhw,idx_ev,evs,ev_dict,extremes):
    '''
    Merge two existing events in the coordinates dictionnary
    '''
    extremes = patch_to_dict(mhw['ev_number'][idx_ev[0][0]],extremes[str(evs[ev_dict])]['time'],extremes[str(evs[ev_dict])]['s_rho'],extremes[str(evs[ev_dict])]['eta_rho'],extremes[str(evs[ev_dict])]['xi_rho'],extremes)
    # delete the merged extremes from coordinates dictionary
    extremes.pop(str(evs[ev_dict]))
    return extremes


def merge_mhws(mhw,idx_tmp_ev,idx_ev,additional_diags):
    '''
    Merge two existing events in the characteristics dictionnary
    '''
    if mhw['index_start'][idx_tmp_ev[0][0]]<mhw['index_start'][idx_ev[0][0]]:
        mhw['index_start'][idx_ev[0][0]] = mhw['index_start'][idx_tmp_ev[0][0]]
    
    if mhw['index_end'][idx_tmp_ev[0][0]]>mhw['index_end'][idx_ev[0][0]]:
        # needed in the case of an addition occurring in the event being merged prior to its merging because its final time can be later than the previous chunk end time
        mhw['index_end'][idx_ev[0][0]] = mhw['index_end'][idx_tmp_ev[0][0]]
    mhw['duration'][idx_ev[0][0]] = mhw['index_end'][idx_ev[0][0]] - mhw['index_start'][idx_ev[0][0]] + 1
        
    # Both volumes are added
    mhw['volume_tot'][idx_ev[0][0]] += mhw['volume_tot'][idx_tmp_ev[0][0]]
    
    # Maximum intensity between the two
    if mhw['intensity_max'][idx_tmp_ev[0][0]]>mhw['intensity_max'][idx_ev[0][0]]:     
        #change intensity max
        mhw['intensity_max'][idx_ev[0][0]] =  mhw['intensity_max'][idx_tmp_ev[0][0]]
        #Change location of Imax
        mhw['loc_Imax_eta'][idx_ev[0][0]] = mhw['loc_Imax_eta'][idx_tmp_ev[0][0]] 
        mhw['loc_Imax_xi'][idx_ev[0][0]] = mhw['loc_Imax_xi'][idx_tmp_ev[0][0]]        
        mhw['time_I_max'][idx_ev[0][0]] = mhw['time_I_max'][idx_tmp_ev[0][0]]
   
    # Both severities are added (cumulative intensities)
    mhw['intensityxvolume_detect_var'][idx_ev[0][0]] += mhw['intensityxvolume_detect_var'][idx_tmp_ev[0][0]]
        
    # Minimum variable value between the two
    if mhw['detect_var_min'][idx_tmp_ev[0][0]]<mhw['detect_var_min'][idx_ev[0][0]]: 
        mhw['detect_var_min'][idx_ev[0][0]] = mhw['detect_var_min'][idx_tmp_ev[0][0]]
    # Both cumulative variable fields are added
    mhw['detect_varxvolume_detect_var'][idx_ev[0][0]] += mhw['detect_varxvolume_detect_var'][idx_tmp_ev[0][0]]
          
    # Spatial and vertical extent
    if mhw['z_min'][idx_tmp_ev[0][0]]>mhw['z_min'][idx_ev[0][0]]:
        mhw['z_min'][idx_ev[0][0]] = mhw['z_min'][idx_tmp_ev[0][0]]
    if mhw['z_max'][idx_tmp_ev[0][0]]<mhw['z_max'][idx_ev[0][0]]:
        mhw['z_max'][idx_ev[0][0]] = mhw['z_max'][idx_tmp_ev[0][0]]
        mhw['z_min_at_z_max'][idx_ev[0][0]] = mhw['z_min_at_z_max'][idx_tmp_ev[0][0]]
        mhw['vertical_extent'][idx_ev[0][0]] = mhw['z_min'][idx_ev[0][0]] - mhw['z_max'][idx_ev[0][0]] 
        mhw['time_z_max'][idx_ev[0][0]] = mhw['time_z_max'][idx_tmp_ev[0][0]] 
        mhw['loc_z_max_eta'][idx_ev[0][0]] = mhw['loc_z_max_eta'][idx_tmp_ev[0][0]] 
        mhw['loc_z_max_xi'][idx_ev[0][0]] =   mhw['loc_z_max_xi'][idx_tmp_ev[0][0]]                                                  

        
    if additional_diags:    
        
        if mhw['omega_min'][idx_tmp_ev[0][0]]<mhw['omega_min'][idx_ev[0][0]]: 
            mhw['omega_min'][idx_ev[0][0]] = mhw['omega_min'][idx_tmp_ev[0][0]]
        mhw['omegaxvolume_omega'][idx_ev[0][0]] += mhw['omegaxvolume_omega'][idx_tmp_ev[0][0]]
        mhw['intensityxvolume_omega'][idx_ev[0][0]] += mhw['intensityxvolume_omega'][idx_tmp_ev[0][0]]
        if mhw['omega_I_max'][idx_tmp_ev[0][0]]>mhw['omega_I_max'][idx_ev[0][0]]:     
            # Update maximum intensity
            mhw['omega_I_max'][idx_ev[0][0]] =  mhw['omega_I_max'][idx_tmp_ev[0][0]]    
        
        if mhw['pH_min'][idx_tmp_ev[0][0]]<mhw['pH_min'][idx_ev[0][0]]: 
            mhw['pH_min'][idx_ev[0][0]] = mhw['pH_min'][idx_tmp_ev[0][0]]
        mhw['pHxvolume_pH'][idx_ev[0][0]] += mhw['pHxvolume_pH'][idx_tmp_ev[0][0]]
        mhw['intensityxvolume_pH'][idx_ev[0][0]] += mhw['intensityxvolume_pH'][idx_tmp_ev[0][0]]
        if mhw['pH_I_max'][idx_tmp_ev[0][0]]>mhw['pH_I_max'][idx_ev[0][0]]:     
            # Update maximum intensity
            mhw['pH_I_max'][idx_ev[0][0]] =  mhw['pH_I_max'][idx_tmp_ev[0][0]]   
                           
        if mhw['oxygen_min'][idx_tmp_ev[0][0]]<mhw['oxygen_min'][idx_ev[0][0]]:
            mhw['oxygen_min'][idx_ev[0][0]] = mhw['oxygen_min'][idx_tmp_ev[0][0]]
        mhw['intensityxvolume_oxygen'][idx_ev[0][0]] += mhw['intensityxvolume_oxygen'][idx_tmp_ev[0][0]]
        mhw['oxygenxvolume_oxygen'][idx_ev[0][0]] += mhw['oxygenxvolume_oxygen'][idx_tmp_ev[0][0]]
    
        if abs(mhw['delta_temp_max'][idx_tmp_ev[0][0]])>abs(mhw['delta_temp_max'][idx_ev[0][0]]):
            mhw['delta_temp_max'][idx_ev[0][0]] = mhw['delta_temp_max'][idx_tmp_ev[0][0]]
        mhw['intensityxvolume_delta_temp'][idx_ev[0][0]] += mhw['intensityxvolume_delta_temp'][idx_tmp_ev[0][0]] 
        
   # remove the event that has been merged from dictionnary
    for k in mhw.keys():
        if k == 'n_events':
            continue 
        mhw[k].pop(idx_tmp_ev[0][0])
        
    return mhw



def add_to_extremes(ev,temp,t_ind,extremes):
    '''
    Add the coordinates of a persisting event to the already existing entry from the previous chunk
    '''
    extremes = patch_to_dict(ev,temp[0]+t_ind,temp[1],temp[2],temp[3],extremes)
    return extremes



def add_to_mhws(mhw,additional_diags,temp,idx_ev,area,dz,zlevs_rho,I,t,t_ind,detect_var,omega,pH,T,O2,thresh,thresh_omega=0,thresh_pH=0,thresh_O2=0):
    '''
    Adjust the characteristics of an already existing entry from the previous chunk with the characteristics of the persisting event in the new chunk 
    '''                 
    t_min = mhw['index_start'][idx_ev[0][0]]
    t_max = max(temp[0])+t[t_ind]
    duration = t_max - t_min + 1
    print(t_min,t_max,duration)
    if duration > mhw['duration'][idx_ev[0][0]]:
        mhw['duration'][idx_ev[0][0]] = duration
        mhw['index_end'][idx_ev[0][0]] = t_max 
    
    ## Size properties                               
    v_gridded = np.multiply(area[temp[2:4]],dz[temp[1:4]])*10**-9 # km3.d
    v = np.sum(v_gridded)
    mhw['volume_tot'][idx_ev[0][0]] += v
    mhw['volume_mean'][idx_ev[0][0]] = mhw['volume_tot'][idx_ev[0][0]]/mhw['duration'][idx_ev[0][0]]
    
    ## Maximum intensity
    I_max = max(I[temp[0:4]])
    if I_max>mhw['intensity_max'][idx_ev[0][0]]:     
        # Change maximum intensity
        mhw['intensity_max'][idx_ev[0][0]] = I_max
        # Change location of maximum intensity 
        mhw['loc_Imax_eta'][idx_ev[0][0]] = temp[2][np.argmax(I[temp[0:4]])]
        mhw['loc_Imax_xi'][idx_ev[0][0]] = temp[3][np.argmax(I[temp[0:4]])]              
        mhw['time_I_max'][idx_ev[0][0]] = t[t_ind] + temp[0][np.argmax(I[temp[0:4]])]

    # Add to cumulative intensity
    IxV = np.sum(np.multiply(I[temp[0:4]],v_gridded))
    mhw['intensityxvolume_detect_var'][idx_ev[0][0]] += IxV
    I_mean = mhw['intensityxvolume_detect_var'][idx_ev[0][0]]/mhw['volume_tot'][idx_ev[0][0]]
    mhw['intensity_mean'][idx_ev[0][0]] = I_mean
    
    # Adjust minimum of variable
    O_min = min(detect_var[temp[0:4]])
    if O_min<mhw['detect_var_min'][idx_ev[0][0]]:     
        # If minimum of variable is lower in the new part of the event (new chunk):
        mhw['detect_var_min'][idx_ev[0][0]] = O_min         
        
    # Add to cumulative variable field
    IxV = np.sum(np.multiply(detect_var[temp[0:4]],v_gridded))
    mhw['detect_varxvolume_detect_var'][idx_ev[0][0]] += IxV
    I_mean = mhw['detect_varxvolume_detect_var'][idx_ev[0][0]]/mhw['volume_tot'][idx_ev[0][0]]
    mhw['detect_var_mean'][idx_ev[0][0]] = I_mean      
    
    # Spatial and vertical extent
    l = np.argmin(zlevs_rho[temp[1:4]])
    l_m_zma =  np.argmax(zlevs_rho[temp[1],temp[2][l],temp[3][l]])
    lm = np.argmax(zlevs_rho[temp[1:4]])
    z_max = zlevs_rho[temp[1][l],temp[2][l],temp[3][l]]
    z_min_max = zlevs_rho[temp[1],temp[2][l],temp[3][l]][l_m_zma] 
    z_min = zlevs_rho[temp[1][lm],temp[2][lm],temp[3][lm]]
    if z_min>mhw['z_min'][idx_ev[0][0]]:
        mhw['z_min'][idx_ev[0][0]] = z_min
    if z_max<mhw['z_max'][idx_ev[0][0]]:
        mhw['z_max'][idx_ev[0][0]] = z_max
        mhw['z_min_at_z_max'][idx_ev[0][0]] = z_min_max
        mhw['vertical_extent'][idx_ev[0][0]] = mhw['z_min'][idx_ev[0][0]] - z_max 
        mhw['time_z_max'][idx_ev[0][0]] = t[t_ind] + temp[0][l]
        mhw['loc_z_max_eta'][idx_ev[0][0]] = temp[2][l]
        mhw['loc_z_max_xi'][idx_ev[0][0]] = temp[3][l]
                           
                         
    if additional_diags:
  
        ################
        #### omega  ####
        ################  
        omega_min = min(omega[temp[0:4]])    
        if omega_min<mhw['omega_min'][idx_ev[0][0]]: 
            mhw['omega_min'][idx_ev[0][0]] = omega_min    
        IxV = np.sum(np.multiply(omega[temp[0:4]],v_gridded))
        mhw['omegaxvolume_omega'][idx_ev[0][0]] += IxV
        mhw['omega_mean'][idx_ev[0][0]] = mhw['omegaxvolume_omega'][idx_ev[0][0]]/mhw['volume_tot'][idx_ev[0][0]]        
        if len(np.shape(thresh_omega))==0:  
            # absolute threshold or none (float, not array)
            thresh_4D_omega = np.zeros_like(omega)+thresh_omega
        elif thresh_omega.ndim == 3:  
            # thresh is a 3D array (fixed baseline)
            thresh_4D_omega = np.repeat(thresh_omega[np.newaxis,: , :, :], omega.shape[0], axis=0)    
        else: 
            # thresh is a 4D array (moving baseline)
            thresh_4D_omega = thresh_omega[:]
        IxV = np.sum(np.multiply(thresh_4D_omega[temp[0:4]]-omega[temp[0:4]],v_gridded))
        mhw['intensityxvolume_omega'][idx_ev[0][0]] += IxV
        mhw['omega_I_mean'][idx_ev[0][0]] = mhw['intensityxvolume_omega'][idx_ev[0][0]]/mhw['volume_tot'][idx_ev[0][0]]  
        omega_I_max = max(thresh_4D_omega[temp[0:4]]-omega[temp[0:4]])
        if omega_I_max>mhw['omega_I_max'][idx_ev[0][0]]:     
            # Update maximum intensity
            mhw['omega_I_max'][idx_ev[0][0]] =  omega_I_max 
        del thresh_4D_omega
    
        
        ################
        ####  pH    ####
        ################    
        pH_min = min(pH[temp[0:4]])
        if pH_min<mhw['pH_min'][idx_ev[0][0]]: 
            mhw['pH_min'][idx_ev[0][0]] = pH_min
        # Add to cumulative intensity
        IxV = np.sum(np.multiply(pH[temp[0:4]],v_gridded))
        mhw['pHxvolume_pH'][idx_ev[0][0]] += IxV
        mhw['pH_mean'][idx_ev[0][0]] = mhw['pHxvolume_pH'][idx_ev[0][0]]/mhw['volume_tot'][idx_ev[0][0]]
        if len(np.shape(thresh_pH))==0:  
            # absolute threshold or none (float, not array)
            thresh_4D_pH = np.zeros_like(pH)+thresh_pH
        elif thresh_pH.ndim == 3:  
            # thresh is a 3D array (fixed baseline)
            thresh_4D_pH = np.repeat(thresh_pH[np.newaxis,: , :, :], pH.shape[0], axis=0) 
        else: 
            # thresh is a 4D array (moving baseline)
            thresh_4D_pH = thresh_pH[:]
        IxV = np.sum(np.multiply(thresh_4D_pH[temp[0:4]]-pH[temp[0:4]],v_gridded))
        mhw['intensityxvolume_pH'][idx_ev[0][0]] += IxV
        mhw['pH_I_mean'][idx_ev[0][0]] = mhw['intensityxvolume_pH'][idx_ev[0][0]]/mhw['volume_tot'][idx_ev[0][0]]
        pH_I_max = max(thresh_4D_pH[temp[0:4]]-pH[temp[0:4]])
        if pH_I_max>mhw['pH_I_max'][idx_ev[0][0]]:     
            # Update maximum intensity
            mhw['pH_I_max'][idx_ev[0][0]] =  pH_I_max
        del thresh_4D_pH 


        ################
        ####  O2    ####
        ################  
        O2_min = min(O2[temp[0:4]])
        if O2_min<mhw['oxygen_min'][idx_ev[0][0]]:
            mhw['oxygen_min'][idx_ev[0][0]] = O2_min
        IxV = np.sum(np.multiply(O2[temp[0:4]],v_gridded))
        mhw['oxygenxvolume_oxygen'][idx_ev[0][0]] += IxV
        mhw['oxygen_mean'][idx_ev[0][0]] = mhw['oxygenxvolume_oxygen'][idx_ev[0][0]]/mhw['volume_tot'][idx_ev[0][0]]
        if len(np.shape(thresh_O2))==0:  
            # absolute threshold or none (float, not array)
            thresh_4D_O2 = np.zeros_like(O2)+thresh_O2
        else:  
            # thresh is a 3D array (fixed baseline)
            thresh_4D_O2 = np.repeat(thresh_O2[np.newaxis,: , :, :], O2.shape[0], axis=0) 
        IxV = np.sum(np.multiply(thresh_4D_O2[temp[0:4]]-O2[temp[0:4]],v_gridded))
        del thresh_4D_O2
        mhw['intensityxvolume_oxygen'][idx_ev[0][0]] += IxV
        mhw['oxygen_I_mean'][idx_ev[0][0]] = mhw['intensityxvolume_oxygen'][idx_ev[0][0]]/mhw['volume_tot'][idx_ev[0][0]]    
        
        ######################
        ####  Temperature ####
        ######################      
        T_ex = [min(T[temp[0:4]]), max(T[temp[0:4]])]
        T_max = T_ex[np.argmax([abs(np.asarray(T_ex))])]
        if abs(T_max)>abs(mhw['delta_temp_max'][idx_ev[0][0]]):
            mhw['delta_temp_max'][idx_ev[0][0]] = T_max
        IxV = np.sum(np.multiply(T[temp[0:4]],v_gridded))
        mhw['intensityxvolume_delta_temp'][idx_ev[0][0]] += IxV
        mhw['delta_temp_mean'][idx_ev[0][0]] = mhw['intensityxvolume_delta_temp'][idx_ev[0][0]]/mhw['volume_tot'][idx_ev[0][0]] 

    return mhw

def load_threshold_linear_trend_removed(thresh_residuals,slope,tdata):
    '''
    Load the 4D threshold field in the case of a moving baseline. 
    '''

    tdata_4D = np.repeat(np.repeat(np.repeat(np.asarray(tdata)[:,np.newaxis,np.newaxis,np.newaxis],slope.shape[0], axis=1),slope.shape[1],axis=2),slope.shape[2],axis=3)

    thresh_4D = np.repeat(slope[np.newaxis,: , :, :], len(tdata), axis=0)*tdata_4D + np.repeat(thresh_residuals[np.newaxis,: , :, :], len(tdata), axis=0)

    del tdata_4D

    return thresh_4D
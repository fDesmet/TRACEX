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
    vertical_extent     [m] difference between shallowest and deepest point over the whole lifetime

    ## Loaction properties
    z_max               [m] shallowest depth reached by the event 
    z_min               [m] deepest depth reached by the event 
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

    delta_temp_max                  [degC] maximum temperature anomaly (positive or negative) 
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
    mhw['vertical_extent'] =[]  # [m] difference between shallowest and deepest point over the whole lifetime

    #  !!!! Note that for all depth-related properties computation, the ROMS grid has the shallowest depth associated to the highest s-level and the deepest depth to the 0 s-level !!!!! #

    # Loaction properties:
    mhw['z_max'] = []           # [m] shallowest depth reached by the event 
    mhw['z_min'] = []           # [m] deepest depth reached by the event
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
    # Update mean volume
    mhw['volume_mean'][idx_ev[0][0]] = mhw['volume_tot'][idx_ev[0][0]]/mhw['duration'][idx_ev[0][0]]

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
    # Update mean intensity
    mhw['intensity_mean'][idx_ev[0][0]] = mhw['intensityxvolume_detect_var'][idx_ev[0][0]]/mhw['volume_tot'][idx_ev[0][0]]
        
    # Minimum variable value between the two
    if mhw['detect_var_min'][idx_tmp_ev[0][0]]<mhw['detect_var_min'][idx_ev[0][0]]: 
        mhw['detect_var_min'][idx_ev[0][0]] = mhw['detect_var_min'][idx_tmp_ev[0][0]]
    # Both cumulative variable fields are added
    mhw['detect_varxvolume_detect_var'][idx_ev[0][0]] += mhw['detect_varxvolume_detect_var'][idx_tmp_ev[0][0]]
    # Update mean value
    mhw['detect_var_mean'][idx_ev[0][0]] = mhw['detect_varxvolume_detect_var'][idx_ev[0][0]]/mhw['volume_tot'][idx_ev[0][0]]
          
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
        mhw['omega_mean'][idx_ev[0][0]] = mhw['omegaxvolume_omega'][idx_ev[0][0]]/mhw['volume_tot'][idx_ev[0][0]]
        mhw['intensityxvolume_omega'][idx_ev[0][0]] += mhw['intensityxvolume_omega'][idx_tmp_ev[0][0]]
        mhw['omega_I_mean'][idx_ev[0][0]] = mhw['intensityxvolume_omega'][idx_ev[0][0]]/mhw['volume_tot'][idx_ev[0][0]]
        if mhw['omega_I_max'][idx_tmp_ev[0][0]]>mhw['omega_I_max'][idx_ev[0][0]]:     
            # Update maximum intensity
            mhw['omega_I_max'][idx_ev[0][0]] =  mhw['omega_I_max'][idx_tmp_ev[0][0]]    
        
        if mhw['pH_min'][idx_tmp_ev[0][0]]<mhw['pH_min'][idx_ev[0][0]]: 
            mhw['pH_min'][idx_ev[0][0]] = mhw['pH_min'][idx_tmp_ev[0][0]]
        mhw['pHxvolume_pH'][idx_ev[0][0]] += mhw['pHxvolume_pH'][idx_tmp_ev[0][0]]
        mhw['pH_mean'][idx_ev[0][0]] = mhw['pHxvolume_pH'][idx_ev[0][0]]/mhw['volume_tot'][idx_ev[0][0]]
        mhw['intensityxvolume_pH'][idx_ev[0][0]] += mhw['intensityxvolume_pH'][idx_tmp_ev[0][0]]
        mhw['pH_I_mean'][idx_ev[0][0]] = mhw['intensityxvolume_pH'][idx_ev[0][0]]/mhw['volume_tot'][idx_ev[0][0]]
        if mhw['pH_I_max'][idx_tmp_ev[0][0]]>mhw['pH_I_max'][idx_ev[0][0]]:     
            # Update maximum intensity
            mhw['pH_I_max'][idx_ev[0][0]] =  mhw['pH_I_max'][idx_tmp_ev[0][0]]   
                           
        if mhw['oxygen_min'][idx_tmp_ev[0][0]]<mhw['oxygen_min'][idx_ev[0][0]]:
            mhw['oxygen_min'][idx_ev[0][0]] = mhw['oxygen_min'][idx_tmp_ev[0][0]]
        mhw['intensityxvolume_oxygen'][idx_ev[0][0]] += mhw['intensityxvolume_oxygen'][idx_tmp_ev[0][0]]
        mhw['oxygen_I_mean'][idx_ev[0][0]] = mhw['intensityxvolume_oxygen'][idx_ev[0][0]]/mhw['volume_tot'][idx_ev[0][0]]
        mhw['oxygenxvolume_oxygen'][idx_ev[0][0]] += mhw['oxygenxvolume_oxygen'][idx_tmp_ev[0][0]]
        mhw['oxygen_mean'][idx_ev[0][0]] = mhw['oxygenxvolume_oxygen'][idx_ev[0][0]]/mhw['volume_tot'][idx_ev[0][0]]

        if abs(mhw['delta_temp_max'][idx_tmp_ev[0][0]])>abs(mhw['delta_temp_max'][idx_ev[0][0]]):
            mhw['delta_temp_max'][idx_ev[0][0]] = mhw['delta_temp_max'][idx_tmp_ev[0][0]]
        mhw['intensityxvolume_delta_temp'][idx_ev[0][0]] += mhw['intensityxvolume_delta_temp'][idx_tmp_ev[0][0]] 
        mhw['delta_temp_mean'][idx_ev[0][0]] =  mhw['intensityxvolume_delta_temp'][idx_ev[0][0]]/mhw['volume_tot'][idx_ev[0][0]]
        
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


### To compute additional properties
import math

def compute_distance_earth(lat1,lon1,lat2,lon2):
    '''
    Compute the distance in kilometers on the Earth based on two lat/lon coordinates point.
    Uses the Haversine formula considering the Earth as a sphere of R=6371km.
    This formula is valid for latitudes close to 40°N or S and gets biased otherwise due to the ellipsoidal shape of the Earth.
    '''
    # change 0-360 °E lon to -180:180 with negative for °W
    lon1 -= 180
    lon2 -= 180
    dlat = (lat2 - lat1)*np.pi/180
    dlon = (lon2 - lon1)*np.pi/180
    phy1 = lat1 *np.pi/180
    phy2 = lat2 *np.pi/180
    R = 6371 # km
    
    # Haversine formula
    a = math.sin(dlat/2)*math.sin(dlat/2) + math.cos(phy1)*math.cos(phy2)*math.sin(dlon/2)*math.sin(dlon/2)
    c = 2*math.atan2(np.sqrt(a),np.sqrt(1-a))
    d = R*c
    return d

def implement_additonal_characteristics(mhw,fdcoast,dcoast_varname,flonlat,lon_varname,lat_varname,zlevs_rho,area,dz,extremes,indices):
    '''
    Implement addtional characteristics to the dictionnary for each event:

    - depth                     [m] time average of volume weighted depth
    - vertOcc                   [m] time average of area weighted vertical occupation of an event
    - dcoast                    [km] time average of volume weighted distance to the coast
    - lat                       [degN] time average of volume weighted latitude
    - lon                       [degE] time average of volume weighted longitute
    - horProp                   [km] distance between the center of gravity of an extreme at initiation and at ending
    - vertProp                  [m] height difference between the center of gravity of an extreme at initiation and at ending: positive value means a deepening; negative value means a shoaling of the events over its lifetime
    - propVel                   [km/d] horizontal propagation divided by duration of an event
    - area_mean                 [km2] area average in time and depth
    - duration_at_cell_mean     [days] average duration at each grid cell involved in the extreme
    - duration_at_cell_max      [days] maximum duration at one grid cell
    - severity_omega_D_I        [days. omega unit] duration times mean intensity with regard to aragonite saturation state
    - severity_pH_D_I           [days. pH unit] duration times mean intensity with regard to pH
    
    '''

    # import sparse matrix module
    import sparse as spa 

    add_car = ['depth','vertOcc','dcoast','lat','lon','horProp','vertProp','propVel','area_mean',
    'duration_at_cell_mean','duration_at_cell_max','severity_omega_D_I','severity_pH_D_I']

    print('Implementing addtional characteristics to the dictionnary: {}'.format(add_car))

    ev_nu = np.asarray(mhw['ev_number'][:])
    # Get number of s-levels
    NZ = dz.shape[0]

    # initialize
    for var in add_car:
        mhw[var] = []
 
    # expand to 3D fields
    lat = netCDF4.Dataset(flonlat, 'r').variables[lat_varname][indices[0]:indices[1],indices[2]:indices[3]]   
    lat3D = np.expand_dims(lat, axis=0)
    lat3D = np.repeat(lat3D,NZ, axis=0)
    del lat
    lon = netCDF4.Dataset(flonlat, 'r').variables[lon_varname][indices[0]:indices[1],indices[2]:indices[3]]   
    lon3D = np.expand_dims(lon, axis=0)
    lon3D = np.repeat(lon3D,NZ, axis=0)
    del lon
    dcoast = netCDF4.Dataset(fdcoast, 'r').variables[dcoast_varname][indices[0]:indices[1],indices[2]:indices[3]]   
    dcoast3D = np.expand_dims(dcoast, axis=0)
    dcoast3D = np.repeat(dcoast3D,NZ, axis=0)
    del dcoast

    # Fill the new characteristics
    for ev,i in zip(ev_nu,range(len(ev_nu))):
        key = str(ev)
        print('Event {}'.format(key))

        # Get first and last time of occurrence
        t_min = np.nanmin(extremes[key]['time'])
        t_max = np.nanmax(extremes[key]['time'])

        # Initialize 
        v_tot = 0
        depth_w = 0
        dcoast_w = 0
        lat_w = 0
        lon_w = 0 
        a_tot_vo = 0
        vo_w = 0
        area_ev_per_time = []

        # severities
        mhw['severity_omega_D_I'].append(np.multiply(mhw['duration'][i],mhw['omega_I_mean'][i]))
        mhw['severity_pH_D_I'].append(np.multiply(mhw['duration'][i],mhw['pH_I_mean'][i]))

        # Loop over lifetime of the event
        for ti in range(t_min,t_max+1):
            # select grid cells hit by the event at that time
            ss = np.asarray(extremes[key]['s_rho'])[np.asarray(extremes[key]['time'])==ti]
            etas = np.asarray(extremes[key]['eta_rho'])[np.asarray(extremes[key]['time'])==ti]
            xis = np.asarray(extremes[key]['xi_rho'])[np.asarray(extremes[key]['time'])==ti]

            
            v_gridded = np.multiply(area[etas,xis],dz[ss,etas,xis])*10**-9  # volume at t=ti [km3] 
            v_tot += np.sum(v_gridded) # add to cumulative volume over lifetime

            # Add the volume-integrated values for depth, dcoast, latitude and longitude for that timestep
            tmp = zlevs_rho[ss,etas,xis]
            depth_w += np.sum(np.multiply(tmp,v_gridded))
            tmp = dcoast3D[ss,etas,xis]
            dcoast_w += np.sum(np.multiply(tmp,v_gridded))
            tmp = lat3D[ss,etas,xis]
            lat_w += np.sum(np.multiply(tmp,v_gridded))
            tmp = lon3D[ss,etas,xis]
            lon_w += np.sum(np.multiply(tmp,v_gridded))
            del tmp

            # Compute area-integrated vertical occupation for that timestep
            data = np.ones(len(ss))
            coords = np.asarray([ss,etas,xis]).astype('int64') 
            # use sparse matrix
            sparse_arr = spa.COO(coords, data,has_duplicates=False)
            #print(sparse_arr.shape)
            # multiply sparse array with dz values
            array_tmp = sparse_arr*dz[:sparse_arr.shape[0],:sparse_arr.shape[1],:sparse_arr.shape[2]]
            vo_t_tmp = array_tmp.sum(axis=0) # do the sum of the extreme cells in depth --> return a 2D spatial array of vertical occupation
            # multiply this vertical occupation by area
            vo_t_a_w = np.multiply(vo_t_tmp,area[:sparse_arr.shape[1],:sparse_arr.shape[2]])
            # Add to area-integrated vertical occupation
            vo_w += np.sum(vo_t_a_w)
            # Compute area involved 
            area_tmp = np.where(vo_t_tmp!=0,area[:sparse_arr.shape[1],:sparse_arr.shape[2]],0) # Count the area only where the vertical occupation is not zero
            a_tot_vo += np.sum(area_tmp) # [m2]

            # Append area hit by event at each s-level       
            ss_min = np.nanmin(ss)
            ss_max = np.nanmax(ss)
            # Loop over depths
            for s in range(ss_min,ss_max+1):
                etas_ss = etas[ss==s]
                xis_ss = xis[ss==s]
                area_ev_per_time.append(np.nansum(area[etas_ss,xis_ss])*10**-6) # [km2]

            if ti ==t_min:
                # If ti is the first timestep, compute center of gravity of the event for propagation
                tmp = lat3D[ss,etas,xis]
                lat_cdg_ini = np.sum(np.multiply(tmp,v_gridded))/np.sum(v_gridded) 
                tmp = lon3D[ss,etas,xis]
                lon_cdg_ini = np.sum(np.multiply(tmp,v_gridded))/np.sum(v_gridded) 
                tmp = zlevs_rho[ss,etas,xis]
                depth_cdg_ini = np.sum(np.multiply(tmp,v_gridded))/np.sum(v_gridded) 
                del tmp
            if ti==t_max:
                # If ti is the last timestep, compute center of gravity of the event for propagation
                tmp = lat3D[ss,etas,xis]
                lat_cdg_fin = np.sum(np.multiply(tmp,v_gridded))/np.sum(v_gridded) 
                tmp = lon3D[ss,etas,xis]
                lon_cdg_fin = np.sum(np.multiply(tmp,v_gridded))/np.sum(v_gridded)                
                tmp = zlevs_rho[ss,etas,xis]
                depth_cdg_fin = np.sum(np.multiply(tmp,v_gridded))/np.sum(v_gridded) 
                del tmp

        # Fill in dictionnary for that event
        mhw['depth'].append(depth_w/v_tot)
        mhw['dcoast'].append(dcoast_w/v_tot)
        mhw['lat'].append(lat_w/v_tot)
        mhw['lon'].append(lon_w/v_tot)

        mhw['vertProp'].append(depth_cdg_fin-depth_cdg_ini) # [m] positive means a deepening; negative means a shoaling of the events over its lifetime
        HP = compute_distance_earth(lat_cdg_ini,lon_cdg_ini,lat_cdg_fin,lon_cdg_fin)
        mhw['horProp'].append(HP)

        mhw['propVel'].append(HP/mhw['duration'][i]) # [km/day]

        mhw['area_mean'].append(np.nanmean(area_ev_per_time)) 
        
        mhw['vertOcc'].append(vo_w/a_tot_vo) # [m]
        del data,coords,sparse_arr,array_tmp,vo_t_tmp,vo_t_a_w,vo_w,a_tot_vo
        del ss,etas,xis,t_min,t_max,v_tot,depth_w,dcoast_w,lat_w,lon_w,v_gridded,HP,lat_cdg_ini,lon_cdg_ini,lon_cdg_fin,lat_cdg_fin,area_ev_per_time,ss_min,ss_max,etas_ss,xis_ss

        # local duration per grid cell 
        times = extremes[key]['time']
        ss = extremes[key]['s_rho']
        etas= extremes[key]['eta_rho']
        xis = extremes[key]['xi_rho'] 
        coords = np.asarray([times,ss,etas,xis]).astype('int64') 
        del ss,times,etas,xis
        # create an 2D array from Z,Y,X coordinates only
        array_2D_3coords = np.vstack((coords[1,:],coords[2,:],coords[3,:]))
        coord_unique,duration_cell = np.unique(array_2D_3coords,axis=1,return_counts=True) # return the cell coordinates and the number of days this cell was extreme  
        mhw['duration_at_cell_mean'].append(np.nanmean(duration_cell))
        mhw['duration_at_cell_max'].append(np.nanmax(duration_cell))
        del duration_cell,coord_unique,array_2D_3coords

    return mhw

def select_events_reaching_top_depth_limit(mhw,zlevs_rho,extremes,depth_lim=100):
    '''
    Remove all events that did not reach a certain depth. 
    '''
    
    print('Removing all extremes that do not reach the top {} m...'.format(depth_lim))

    ev_to_remove = []
    idx_to_remove = []
    for ev,idx in zip(mhw['ev_number'],range(mhw['n_events'])):
        key = str(ev)
        # see if there is a zlevs shallower than depth_limit for this extreme 
        ss = np.asarray(extremes[key]['s_rho'])
        etas = np.asarray(extremes[key]['eta_rho'])
        xis = np.asarray(extremes[key]['xi_rho'])

        shallowest_depth_hit = np.nanmin(zlevs_rho[ss,etas,xis])
        if shallowest_depth_hit > depth_lim:
            # if the shallowest depth hit by this extreme is larger than (thus below) the depth limit of detection 
            ev_to_remove.append(ev)
            idx_to_remove.append(idx)     

    if len(ev_to_remove)>0:
        print('{} events did not reach the {} m depth'.format(len(ev_to_remove),depth_lim))
        # Remove all extremes not reaching the depth limit from coordinates dictionnay 
        for k in ev_to_remove:
            del extremes[str(k)]
        
        # Remove all extremes not reaching the depth limit from characteristics dictionnay   
        # convert the array into list 
        mhw['duration'] = mhw['duration'].tolist()               
        for k in mhw.keys():
            if k == 'n_events':
                continue
            for track,n in zip(idx_to_remove,range(len(idx_to_remove))):
                mhw[k].pop(int(track-n))          
        # Replace the number of events by the length of the dictionnary to have the number of events fulfilling conditions
        mhw['n_events'] = len(mhw['ev_number'])  
    else:
        print('All events reach the {} m depth. Nothing to remove.'.format(depth_lim))

    return mhw,extremes
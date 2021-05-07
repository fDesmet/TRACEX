# TRACE-X
***TRACking EXtremes***<br>
Detecting marine extremes in space (3D) and time.<br>
Extremes with regards to pH, aragonite saturation state, oxygen, temperature or any other variable are detected in space (3D) and time from large model output. <br>

![](images/TRACEX_gif.gif)

## Content of the package
This package contains:
1. A namelist python file (define_TRACEX_variables.py) setting ALL necessary parameters and path-names
2. A module containing all necessary functions (TRACEX_functions.py)
3. The main routine to detect the extremes (launch_TRACEX.py) calling both module and namelist 
4. A simple example on a 2D array with a step-by-step description of the algorithm (TRACEX_unit_test folder)

## Goal of the package
The package aims at tracking in time spatially coherent structures of grid cells fulfilling a given criterion in large datasets, such as earth system model output, that cannot be treated as one block due to memory exceedance and thus need to be chunked. <br>
This algortithm takes care of the detection of coherent structures in each given chunk, of the merging across successive chunks and of the derivation of a few characterisitcs related to extreme events. <br>
The output of the main routine are two Python dictionaries:
1. A dictionnary with all coordinates (T,Z,Y,X, for ROMS time,s-level,eta,xi) of all detected extremes
2. A dictionary of global statistical properties summarizing and characerizing the extremes  

## Description of the package
At the core of the package is the ndimage.label() function, which is used to detect "neighbouring" 1s in a boolean array and to group them together as patches. In this package, the boolean array is the result of the question "Is the variable X 'extreme' at the analyzed locations?", comparing a value matrix to a threshold matrix. The patches returned by ndimage.label() can hence be interpreted as stand-alone extremes. The threshold array needs to be provided as input, based on absolute or relative threshold approaches. <br>
**Note**: The algorithm is written for extremes defined as **below** a certain threshold (useful for oxygen, pH, aragonite saturation state for instance). It needs slight modification if the extremes are defined as above a certain threshold (heat waves for instance). 

"Neighboring" can be chosen to mean, whether pixels need to share edges, corners, or both of them. It can also be adapted to connect the cells only along certain dimensions. <br>
As model output (ROMS in this case) can be memory heavy, this package allows to chunk the data in time, i.e. calculate the patches over a limited number of output time steps. When transitioning from one chunk to the next chunk, the package ensures, that extremes that were separated by chunking the data, are rejoined. <br>

<p align="center">
  <img width="100%" src="/TRACEX_unit_test/boolean_array_TRACEX_unit_test.png" alt="Boolean array 2D example"/>
</p>
<p align="center">
  <img width="100%" src="/TRACEX_unit_test/identified_extremes_TRACEX_unit_test.png" alt="Extremes 2D example"/>
</p>

The package:
1. is compatible with different domain set-up of ROMS and should also be easily transferable to other model output<br>
**Note**: ROMS run on terrain following s-coordinates for which the highest s-level represents the surface and the s=0 level represents the deepest point at a given location. 
2. can process different time periods, output frequency with different chunk lengthstimesteps
3. can ignore extremes in certain parts of the model output (i.e. below a certain depth or outside a certain domain, given as a mask)
4. can require the extremes to have a minimum duration 
5. can detect compound extremes, ie patches simultaneously extreme with regard to two different variables. 

## Input requirement
A comprehensive list of input needed to run the package can be found here:
1. Variables and thresholds for detection: the variable name (in model outputs) of the variable used for the detection and the type of threshold: 'var', 'thresh' *+ if compound extremes: 'var2', 'thresh2'*
2. Settings: the different settings (depth limit, minimum duration, compound/single extremes, moving/fixed basline, type of connectivity (edge/corners), restart)
3. Model outputs:
    - the time period: years vector and associated daily (or other frequency depending on the model output frequency) vector  
    - the path to the daily model outputs (one file per year)
    - the setup used in the filenames of model output (format "setup_year_avg.nc")
4. Output: name of output folder 
5. Grid data:
    - number of z levels 
    - 2D array of the area of each model grid cell ($m^2$)
    - 3D mask array of 0 and 1 masking out the regions to discard from the study 
    - 3D field of grid cells height 
    - 3D field of grid cells depth 
6. Threshold field (if relative threshold):
    - if thresh *(or thresh2)* is a string (= not absolute threshold): the threshold field 'thresh_main' *+ if compound 'thresh_main2' + if moving baseline 'slope_main' field (+'slope_main2')*
    - the threshold to be used for additional diagnostics such as intensity of oxygen or temperature; if 0, then the mean value will be returned instead of intensity. 
7. The 365 and 366 days climatologies (or other frequency depending on model output frequency) for temperature to compute delta_temp in extremes 

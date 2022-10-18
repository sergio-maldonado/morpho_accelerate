Acceleration of morphodynamic simulations based on local trends
in the bed evolution
Ellie Newell

- Introduction

This project aims to improve the method used at accelerate the numerical simulation of morphodynamic change. Sample files have been included in this folder to simulate a Linear profile with a Mild slope and Fine sediment for 24 hours under energetic wave conditions.

- Contents of this folder

MATLAB scripts:
AccelerationAlgorithm.m - Matlab script for the acceleration algorithm
accelerate_simulation.m - Matlab function for the extrapolation of an XBeach simulation
extract_jons.m - Matlab funstion for the extraction of variables from jonswap.txt file
extract_params.m - Matlab function for the extraction of variables from params.txt file

XBeach files:
bed.dep - sample input bathymetry file
x.grd - sample x coordinate file
y.grd - sample y coordinate file
params.txt - sample input parameters file
jonswap.txt - sample jonswap parameters file

- Other requirements

The acceleration algorithm requires the XBeach MATLAB toolbox (which can be downloaded from https://oss.deltares.nl/web/xbeach/tools) and the XBeach executable (which can be downloaded from  https://download.deltares.nl/en/download/xbeach-open-source/)

- Using the acceleration algorithm

Running the AcccelerationAlgorithm.m script (ensuring that the XBeach executable and input files are all in the current folder) will start a series of XBeach simulations, the progress of which will appear in the MATLAB command window. When the full simulation is complete, the final predicted beach profile will be plotted in a figure alongside the initial profile for comparison. 

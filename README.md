# SCM responses
This repository contains the data and scripts required to reproduce the results of the paper "Characterizing Convection Schemes Using Their Responses to Convective Tendency Perturbations" submitted to the *Journal of Advances of Modeling Earth Systems (JAMES)*. 

## Brief description of project

An SCM intercomparison project based on the Linear Response Function framework of Kuang (2010), where we examine the temperature and moisture responses to small convective tendency perturbations.

Participating Single-Column Models (SCMs): WRF, LMDZ, CNRM, UM, SCAM

Convection schemes tested:
1. WRF Kain-Fritsch
2. WRF New-Tiedtke
3. WRF New-Simplified-Arakawa-Schubert
4. WRF Betts-Miller-Janjic
5. WRF Zhang-McFarlane 
6. CNRM PCMT scheme
7. UM Simplified-Betts-Miller
8. UM Mass-Flux scheme (Gregory & Rowntree)
9. SCAM Zhang-McFarlane
10. LMDZ modified-Emanuel scheme + Cold pool formulation

## Repository structure

1. ```/data``` directory: contains the data from the SCMs, in ```csv``` format
   - ```/data/[model_name]/``` directories: contain the data of the individual models. There are four sub-directories (five for WRF) for each model:
     - ```REF``` folder: contains data for RCE mean state, including Temperature (T) and Relative Humidity (RH)
     - ```matrix_X_raw```folder: contains raw data for the *T* and *q* responses to *dT/dt* and *dq/dt* perturbations
     - ```matrix_M_inv``` folder: contains the post-processed (normalized and standardized) **M**<sup>-1</sup> matrix data 
     - ```response_profiles``` folder: contains the post-processed response profiles (vertical column) for perturbation at two levels (850 and 650 hPa)
     - (for WRF only) ```pbl_mp_sensitivity``` folder: data for WRF PBL and MP sensitivity tests. There are three sub-directories in this folder:
       - ```mean_states``` folder: contains data for the RCE mean state sensitivity to PBL and MP schemes
       - ```response_profiles``` folder: contains data for the sensitivity of *T* and *q* responses to PBL and MP schemes
       - ```response_profiles_non_idealized``` folder: contains data to compare sensitivity of *T* and *q* responses to PBL and MP schemes between idealized and non-idealized setups
3. ```/scripts``` directory: contains the python scripts to postprocess and plot the figures in the paper

   
 
  

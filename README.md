# SCM responses
This repository contains the data and scripts required to reproduce the results of the paper "Characterizing Convection Schemes Using Their Responses to Convective Tendency Perturbations" published in the *Journal of Advances of Modeling Earth Systems (JAMES)*. 

**The published manuscript about this work should be cited as:**  
Hwong, Y. L., Song, S., Sherwood, S. C., Stirling, A. J., Rio, C., Roehrig, R., ... & Touzé‐Peiffer, L. (2021). Characterizing convection schemes using their responses to imposed tendency perturbations. Journal of Advances in Modeling Earth Systems, 13(5), e2021MS002461, https://doi.org/10.1029/2021MS002461 

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
     - ```REF``` directory: contains data for RCE mean state, including Temperature (T) and Relative Humidity (RH)
     - ```matrix_X_raw```directory: contains raw data for the *T* and *q* responses to *dT/dt* and *dq/dt* perturbations
     - ```matrix_M_inv``` directory: contains the post-processed (normalized and standardized) **M**<sup>-1</sup> matrix data 
     - ```response_profiles``` directory: contains the post-processed response profiles (vertical column) for perturbation at two levels (850 and 650 hPa)
     - (for WRF only) ```pbl_mp_sensitivity``` directory: data for WRF PBL and MP sensitivity tests. There are three sub-directories in this folder:
       - ```mean_states``` directory: contains data for the RCE mean state sensitivity to PBL and MP schemes
       - ```response_profiles``` directory: contains data for the sensitivity of *T* and *q* responses to PBL and MP schemes
       - ```response_profiles_non_idealized``` directory: contains data to compare sensitivity of *T* and *q* responses to PBL and MP schemes between idealized and non-idealized setups
3. ```/scripts``` directory: contains the python scripts to post-process and plot the figures in the paper
   - ```/scripts/plot_inidividual_matrix/``` directory: contains scripts to plot the **M**<sup>-1</sup> matrix and the 2-levels response profiles for individual SCMs
     - ```plot_matrix.py``` : script to post-process and plot the matrix and response profiles for selected SCM using the raw data in the ```matrix_X_raw``` folder. Option available to save the (post-processed) outputs as ```csv``` files in the ```matrix_M_inv``` and ```response_profiles``` folders of the selected SCM
   - ```/scripts/plot_figures/``` directory: contains scripts to plot the Figures in the paper
     - ```plot_rce_mean_states_all_models.py``` : script to plot **Figure 1** (RCE mean state T and RH of all SCMs)
     - ```plot_anomaly_profiles.py``` : script to plot **Figure 2 - 3** (vertical profiles of *T* and *q* responses for selected SCMs)
     - ```plot_matrix_all_models.py``` : script to plot **Figures 4 - 7** (**M**<sup>-1</sup> matrices of all SCMs)
     - ```plot_rh_q_correlation_matrix.py``` : script to plot **Figure 8** (correlation matrix of RH vs. q')
     - ```plot_pbl_mp_sensitivity_mean_states.py``` : script to plot **Figures 9 - 10** (sensitivity of RCE mean state to PBL and MP schemes)
     - ```plot_pbl_mp_sensitivity_responses.py``` : script to plot **Figures 11 - 12** (sensitivity of *T* and *q* responses to PBL and MP schemes)
     - ```plot_pbl_mp_sensitivity_responses_non_idealized.py``` : script to plot **Figures A1 - A2** (sensitivity of *T* and *q* responses to PBL and MP schemes comparison between idealized and non-idealized setups)

   
 
  

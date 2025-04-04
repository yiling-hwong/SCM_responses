#!/usr/bin/env python3

"""
Created on 13 August 2020

@author: Yi-Ling HWONG <yiling.hwong@gmail.com>

This script plots the four quadrats of the M^-1 matrix (Kuang, 2010)
in either matrix form, or vertical profiles.

"""

from sam import SAM
from wrf import WRF
from umbm import UMBM
from ummf import UMMF
from lmdz import LMDZ
from cnrm import CNRM
from scam import SCAM



"""
1. SET object parameters here
"""

# Set model. Options: sam,wrf,umbm,ummf,lmdz,cnrm,scam
model = "sam" # sam,wrf,umbm,ummf,lmdz,cnrm,scam

# For WRF only. Options: kfeta,ntiedtke,nsas,camzm,bmj
scheme = "camzm"

# For LMDZ only. Options: 5a, 6a, 6ab
lmdz_version = "6ab"

# Set to True if standardise to Kuang's (SAM) CRM power input
standardise_kuang = True

# Select functions to run
plot_M_inv_matrix = False
plot_anomaly_profiles_2_levels = True # plot anomaly profiles for 850 and 650 hPa

# Set one of the following to True (perturb dT/dt or dq/dt)
perturb_t = False
perturb_q = True

# Set state anomaly to either "T" or "q"
state_anomaly = "T"

# Perturbation amplitude (0.5 K/d & 0.2 g/kg/d, or 0.2 K/d & 0.1 g/kg/d )
t_amplitude = 0.5
q_amplitude = 0.2

"""
2. Create object for selected Model
"""

# SAM: 8 for 850hpa, 10 for 730hpa, 11 for 650hpa, 12 for 560hpa
# WRF: 12 for 850hpa, 21 for 730hpa, 27 for 650hpa, 33 for 560hpa
# UMBM & UMMF: 12 for 850hpa, 18 for 730hpa, 21 for 650hpa, 24 for 560hpa
# LMDZ: 21 for 850hpa, 26 for 730hpa, 29 for 650hpa, 31 for ~560hpa (579hpa)
# CNRM: 15 for 850hpa, 20 for 730hpa, 23 for 650hpa, 26 for 560hpa
# SCAM: 11 for 850hpa, 16 for 730hpa, 18 for 650hpa, 21 for ~560hpa (546hpa)

label_level_1 = "850 hPa"
label_level_2 = "650 hPa"

if model == "sam":
    Model = SAM(standardise_kuang=standardise_kuang, perturb_t=perturb_t, perturb_q=perturb_q,state_anomaly=state_anomaly,t_amplitude=t_amplitude, q_amplitude=q_amplitude)
    target_level_1 = 8
    target_level_2 = 11
elif model == "wrf":
    Model = WRF(scheme=scheme,standardise_kuang=standardise_kuang,perturb_t=perturb_t,perturb_q=perturb_q,state_anomaly=state_anomaly,t_amplitude=t_amplitude,q_amplitude=q_amplitude)
    target_level_1 = 12
    target_level_2 = 27
elif model == "umbm":
    Model = UMBM(standardise_kuang=standardise_kuang,perturb_t=perturb_t,perturb_q=perturb_q,state_anomaly=state_anomaly,t_amplitude=t_amplitude,q_amplitude=q_amplitude)
    target_level_1 = 12
    target_level_2 = 21
elif model == "ummf":
    Model = UMMF(standardise_kuang=standardise_kuang,perturb_t=perturb_t,perturb_q=perturb_q,state_anomaly=state_anomaly,t_amplitude=t_amplitude,q_amplitude=q_amplitude)
    target_level_1 = 12
    target_level_2 = 21
elif model == "lmdz":
    Model = LMDZ(lmdz_version=lmdz_version,standardise_kuang=standardise_kuang,perturb_t=perturb_t,perturb_q=perturb_q,state_anomaly=state_anomaly,t_amplitude=t_amplitude,q_amplitude=q_amplitude)
    target_level_1 = 21
    target_level_2 = 29
elif model == "cnrm":
    Model = CNRM(standardise_kuang=standardise_kuang,perturb_t=perturb_t,perturb_q=perturb_q,state_anomaly=state_anomaly,t_amplitude=t_amplitude,q_amplitude=q_amplitude)
    target_level_1 = 15
    target_level_2 = 23
elif model == "scam":
    Model = SCAM(standardise_kuang=standardise_kuang,perturb_t=perturb_t,perturb_q=perturb_q,state_anomaly=state_anomaly,t_amplitude=t_amplitude,q_amplitude=q_amplitude)
    target_level_1 = 11
    target_level_2 = 18

"""
3. Plot selected 
"""

if plot_M_inv_matrix == True:
    Model.plot_matrix_M_inv(write_m_inv_to_file=False,vmax_kuang_t=0.5,vmax_kuang_q=0.3,vmax_power_t=0.1,vmax_power_q=0.05)

if plot_anomaly_profiles_2_levels == True:
    Model.plot_anomaly_profile(write_anomaly_to_file=False,target_level_1=target_level_1,target_level_2=target_level_2,label_level_1=label_level_1,label_level_2=label_level_2)



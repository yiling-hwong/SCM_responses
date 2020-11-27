#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on 13 August 2020

@author: Yi-Ling HWONG <yiling.hwong@gmail.com>

This script plots the RCE mean states of all the models (SAM and SCMs)
Temperature (theta_es) and RH are plotted in the same figure

"""

import numpy as np
import matplotlib.pyplot as plt


def get_temperature_data():

    """
    Get RCE Temperature profile of all models
    :return: 2 lists: T and pressures of all models
    """

    lines0 = open("../../data/SAM/REF/T.csv","r").readlines()
    lines0b = open("../../data/SAM/REF/p.csv", "r").readlines()
    lines1 = open("../../data/WRF/REF/kfeta_pbl_ysu_T.csv", "r").readlines()
    lines2 = open("../../data/WRF/REF/ntiedtke_pbl_ysu_T.csv", "r").readlines()
    lines3 = open("../../data/WRF/REF/nsas_pbl_ysu_T.csv", "r").readlines()
    lines4 = open("../../data/WRF/REF/bmj_pbl_ysu_T.csv", "r").readlines()
    lines5 = open("../../data/WRF/REF/camzm_pbl_mp_uw_T.csv", "r").readlines()
    lines6 = open("../../data/CNRM/REF/RCE_mean_T.csv", "r").readlines()
    lines7 = open("../../data/UMMF/REF/RCE_mean_T.csv", "r").readlines()
    lines8 = open("../../data/UMBM/REF/RCE_mean_T.csv", "r").readlines()
    lines9 = open("../../data/SCAM/REF/RCE_mean_T.csv", "r").readlines()
    lines10 = open("../../data/LMDZ/REF/lmdz6a/RCE_mean_T.csv","r").readlines()
    lines11 = open("../../data/LMDZ/REF/lmdz6ab/RCE_mean_T.csv", "r").readlines()

    sam_ref = []
    sam_pres = []
    lmdz6a_ref = []
    lmdz6a_pres = []
    lmdz6ab_ref = []
    lmdz6ab_pres = []
    scam_ref = []
    scam_pres = []
    wrf_zm_ref = []
    wrf_zm_pres = []
    wrf_kf_ref = []
    wrf_kf_pres = []
    wrf_nt_ref = []
    wrf_nt_pres = []
    wrf_nsas_ref = []
    wrf_nsas_pres = []
    wrf_bmj_ref = []
    wrf_bmj_pres = []
    umbm_ref = []
    umbm_pres = []
    ummf_ref = []
    ummf_pres = []
    cnrm_ref = []
    cnrm_pres = []

    ref_all = []
    pressures_all = []


    # 0. SAM (Kuang CRM)
    for line in lines0:
        spline = line.rstrip("\n")
        sam_ref.append(float(spline))

    for line in lines0b:
        spline = line.rstrip("\n")
        sam_pres.append(float(spline))

    if sam_pres[0] > sam_pres[1]:
        sam_pres.reverse()
        sam_ref.reverse()

    ref_all.append(sam_ref)
    pressures_all.append(sam_pres)

    # 1. WRF kfeta
    for line in lines1:
        spline = line.rstrip("\n").split(",")
        wrf_kf_pres.append(float(spline[0]))
        wrf_kf_ref.append(float(spline[1]))

    if wrf_kf_pres[0] > wrf_kf_pres[1]:
        wrf_kf_pres.reverse()
        wrf_kf_ref.reverse()

    ref_all.append(wrf_kf_ref)
    pressures_all.append(wrf_kf_pres)

    # 2. WRF ntiedtke
    for line in lines2:
        spline = line.rstrip("\n").split(",")
        wrf_nt_pres.append(float(spline[0]))
        wrf_nt_ref.append(float(spline[1]))

    if wrf_nt_pres[0] > wrf_nt_pres[1]:
        wrf_nt_pres.reverse()
        wrf_nt_ref.reverse()

    ref_all.append(wrf_nt_ref)
    pressures_all.append(wrf_nt_pres)

    # 3. WRF nsas
    for line in lines3:
        spline = line.rstrip("\n").split(",")
        wrf_nsas_pres.append(float(spline[0]))
        wrf_nsas_ref.append(float(spline[1]))

    if wrf_nsas_pres[0] > wrf_nsas_pres[1]:
        wrf_nsas_pres.reverse()
        wrf_nsas_ref.reverse()

    ref_all.append(wrf_nsas_ref)
    pressures_all.append(wrf_nsas_pres)

    # 4. WRF bmj
    for line in lines4:
        spline = line.rstrip("\n").split(",")
        wrf_bmj_pres.append(float(spline[0]))
        wrf_bmj_ref.append(float(spline[1]))

    if wrf_bmj_pres[0] > wrf_bmj_pres[1]:
        wrf_bmj_pres.reverse()
        wrf_bmj_ref.reverse()

    ref_all.append(wrf_bmj_ref)
    pressures_all.append(wrf_bmj_pres)

    # 5. WRF camzm
    for line in lines5:
        spline = line.rstrip("\n").split(",")
        wrf_zm_pres.append(float(spline[0]))
        wrf_zm_ref.append(float(spline[1]))

    if wrf_zm_pres[0] > wrf_zm_pres[1]:
        wrf_zm_pres.reverse()
        wrf_zm_ref.reverse()

    ref_all.append(wrf_zm_ref)
    pressures_all.append(wrf_zm_pres)

    # 6. CNRM
    for line in lines6:
        spline = line.rstrip("\n").split(",")
        cnrm_pres.append(float(spline[0]))
        cnrm_ref.append(float(spline[1]))

    if cnrm_pres[0] > cnrm_pres[1]:
        cnrm_pres.reverse()
        cnrm_ref.reverse()

    ref_all.append(cnrm_ref)
    pressures_all.append(cnrm_pres)


    # 7. UMMF
    for line in lines7:
        spline = line.rstrip("\n").split(",")
        ummf_pres.append(float(spline[0]))
        ummf_ref.append(float(spline[1]))

    if ummf_pres[0] > ummf_pres[1]:
        ummf_pres.reverse()
        ummf_ref.reverse()

    ref_all.append(ummf_ref)
    pressures_all.append(ummf_pres)

    # 8. UMBM
    for line in lines8:
        spline = line.rstrip("\n").split(",")
        umbm_pres.append(float(spline[0]))
        umbm_ref.append(float(spline[1]))

    if umbm_pres[0] > umbm_pres[1]:
        umbm_pres.reverse()
        umbm_ref.reverse()

    ref_all.append(umbm_ref)
    pressures_all.append(umbm_pres)

    # 9. SCAM
    for line in lines9:
        spline = line.rstrip("\n").split(",")
        scam_pres.append(float(spline[0]))
        scam_ref.append(float(spline[1]))

    if scam_pres[0] > scam_pres[1]:
        scam_pres.reverse()
        scam_ref.reverse()

    ref_all.append(scam_ref)
    pressures_all.append(scam_pres)

    # 10. LMDZ6A
    for line in lines10:
        spline = line.rstrip("\n").split(",")
        lmdz6a_pres.append(float(spline[0]))
        lmdz6a_ref.append(float(spline[1]))

    if lmdz6a_pres[0] > lmdz6a_pres[1]:
        lmdz6a_pres.reverse()
        lmdz6a_ref.reverse()

    ref_all.append(lmdz6a_ref)
    pressures_all.append(lmdz6a_pres)

    # 11. LMDZ6Ab
    for line in lines11:
        spline = line.rstrip("\n").split(",")
        lmdz6ab_pres.append(float(spline[0]))
        lmdz6ab_ref.append(float(spline[1]))

    if lmdz6ab_pres[0] > lmdz6ab_pres[1]:
        lmdz6ab_pres.reverse()
        lmdz6ab_ref.reverse()

    ref_all.append(lmdz6ab_ref)
    pressures_all.append(lmdz6ab_pres)

    print ()
    print ("Length of RCE T list (all models):", len(ref_all),len(pressures_all))
    print ()

    return ref_all,pressures_all


def get_RH_data():

    """
    Get RCE RH profile of all models
    :return: 2 lists: RH and pressures of all models
    """

    lines0 = open("../../data/SAM/REF/rh.csv", "r").readlines()
    lines0b = open("../../data/SAM/REF/p.csv", "r").readlines()
    lines1 = open("../../data/WRF/REF/kfeta_pbl_ysu_rh.csv", "r").readlines()
    lines2 = open("../../data/WRF/REF/ntiedtke_pbl_ysu_rh.csv", "r").readlines()
    lines3 = open("../../data/WRF/REF/nsas_pbl_ysu_rh.csv", "r").readlines()
    lines4 = open("../../data/WRF/REF/bmj_pbl_ysu_rh.csv", "r").readlines()
    lines5 = open("../../data/WRF/REF/camzm_pbl_mp_uw_rh.csv", "r").readlines()
    lines6 = open("../../data/CNRM/REF/RCE_mean_RH.csv", "r").readlines()
    lines7 = open("../../data/UMMF/REF/RCE_mean_RH.csv", "r").readlines()
    lines8 = open("../../data/UMBM/REF/RCE_mean_RH.csv", "r").readlines()
    lines9 = open("../../data/SCAM/REF/RCE_mean_RH.csv", "r").readlines()
    lines10 = open("../../data/LMDZ/REF/lmdz6a/RCE_mean_RH.csv", "r").readlines()
    lines11 = open("../../data/LMDZ/REF/lmdz6ab/RCE_mean_RH.csv", "r").readlines()

    sam_ref = []
    sam_pres = []
    lmdz6a_ref = []
    lmdz6a_pres = []
    lmdz6ab_ref = []
    lmdz6ab_pres = []
    scam_ref = []
    scam_pres = []
    wrf_zm_ref = []
    wrf_zm_pres = []
    wrf_kf_ref = []
    wrf_kf_pres = []
    wrf_nt_ref = []
    wrf_nt_pres = []
    wrf_nsas_ref = []
    wrf_nsas_pres = []
    wrf_bmj_ref = []
    wrf_bmj_pres = []
    umbm_ref = []
    umbm_pres = []
    ummf_ref = []
    ummf_pres = []
    cnrm_ref = []
    cnrm_pres = []

    ref_all = []
    pressures_all = []

    # 0. SAM (Kuang CRM)
    for line in lines0:
        spline = line.rstrip("\n")
        sam_ref.append(float(spline))

    for line in lines0b:
        spline = line.rstrip("\n")
        sam_pres.append(float(spline))

    if sam_pres[0] > sam_pres[1]:
        sam_pres.reverse()
        sam_ref.reverse()

    ref_all.append(sam_ref)
    pressures_all.append(sam_pres)

    # 1. WRF kfeta
    for line in lines1:
        spline = line.rstrip("\n").split(",")
        wrf_kf_pres.append(float(spline[0]))
        wrf_kf_ref.append(float(spline[1]))

    if wrf_kf_pres[0] > wrf_kf_pres[1]:
        wrf_kf_pres.reverse()
        wrf_kf_ref.reverse()

    ref_all.append(wrf_kf_ref)
    pressures_all.append(wrf_kf_pres)

    # 2. WRF ntiedtke
    for line in lines2:
        spline = line.rstrip("\n").split(",")
        wrf_nt_pres.append(float(spline[0]))
        wrf_nt_ref.append(float(spline[1]))

    if wrf_nt_pres[0] > wrf_nt_pres[1]:
        wrf_nt_pres.reverse()
        wrf_nt_ref.reverse()

    ref_all.append(wrf_nt_ref)
    pressures_all.append(wrf_nt_pres)

    # 3. WRF nsas
    for line in lines3:
        spline = line.rstrip("\n").split(",")
        wrf_nsas_pres.append(float(spline[0]))
        wrf_nsas_ref.append(float(spline[1]))

    if wrf_nsas_pres[0] > wrf_nsas_pres[1]:
        wrf_nsas_pres.reverse()
        wrf_nsas_ref.reverse()

    ref_all.append(wrf_nsas_ref)
    pressures_all.append(wrf_nsas_pres)

    # 4. WRF bmj
    for line in lines4:
        spline = line.rstrip("\n").split(",")
        wrf_bmj_pres.append(float(spline[0]))
        wrf_bmj_ref.append(float(spline[1]))

    if wrf_bmj_pres[0] > wrf_bmj_pres[1]:
        wrf_bmj_pres.reverse()
        wrf_bmj_ref.reverse()

    ref_all.append(wrf_bmj_ref)
    pressures_all.append(wrf_bmj_pres)

    # 5. WRF camzm
    for line in lines5:
        spline = line.rstrip("\n").split(",")
        wrf_zm_pres.append(float(spline[0]))
        wrf_zm_ref.append(float(spline[1]))

    if wrf_zm_pres[0] > wrf_zm_pres[1]:
        wrf_zm_pres.reverse()
        wrf_zm_ref.reverse()

    ref_all.append(wrf_zm_ref)
    pressures_all.append(wrf_zm_pres)

    # 6. CNRM
    for line in lines6:
        spline = line.rstrip("\n").split(",")
        cnrm_pres.append(float(spline[0]))
        cnrm_ref.append(float(spline[1]))

    if cnrm_pres[0] > cnrm_pres[1]:
        cnrm_pres.reverse()
        cnrm_ref.reverse()

    ref_all.append(cnrm_ref)
    pressures_all.append(cnrm_pres)

    # 7. UMMF
    for line in lines7:
        spline = line.rstrip("\n").split(",")
        ummf_pres.append(float(spline[0]))
        ummf_ref.append(float(spline[1]))

    if ummf_pres[0] > ummf_pres[1]:
        ummf_pres.reverse()
        ummf_ref.reverse()

    ref_all.append(ummf_ref)
    pressures_all.append(ummf_pres)

    # 8. UMBM
    for line in lines8:
        spline = line.rstrip("\n").split(",")
        umbm_pres.append(float(spline[0]))
        umbm_ref.append(float(spline[1]))

    if umbm_pres[0] > umbm_pres[1]:
        umbm_pres.reverse()
        umbm_ref.reverse()

    ref_all.append(umbm_ref)
    pressures_all.append(umbm_pres)

    # 9. SCAM
    for line in lines9:
        spline = line.rstrip("\n").split(",")
        scam_pres.append(float(spline[0]))
        scam_ref.append(float(spline[1]))

    if scam_pres[0] > scam_pres[1]:
        scam_pres.reverse()
        scam_ref.reverse()

    ref_all.append(scam_ref)
    pressures_all.append(scam_pres)

    # 10. LMDZ6A
    for line in lines10:
        spline = line.rstrip("\n").split(",")
        lmdz6a_pres.append(float(spline[0]))
        lmdz6a_ref.append(float(spline[1]))

    if lmdz6a_pres[0] > lmdz6a_pres[1]:
        lmdz6a_pres.reverse()
        lmdz6a_ref.reverse()

    ref_all.append(lmdz6a_ref)
    pressures_all.append(lmdz6a_pres)

    # 11. LMDZ6Ab
    for line in lines11:
        spline = line.rstrip("\n").split(",")
        lmdz6ab_pres.append(float(spline[0]))
        lmdz6ab_ref.append(float(spline[1]))

    if lmdz6ab_pres[0] > lmdz6ab_pres[1]:
        lmdz6ab_pres.reverse()
        lmdz6ab_ref.reverse()

    ref_all.append(lmdz6ab_ref)
    pressures_all.append(lmdz6ab_pres)

    print()
    print("Length of RCE RH list (all models):", len(ref_all), len(pressures_all))
    print()

    return ref_all, pressures_all

def get_theta_es():

    """
    Convert T to saturation equivalent potential temperature (theta_es)
    :return: list of theta_es for all models
    """

    print()
    print("########################")
    print("Getting theta_es...")
    print()

    import metpy.calc as mpcalc
    from metpy.units import units

    temperature_all_models = get_temperature_data()[0]
    pressure_all_models = get_temperature_data()[1]

    theta_es_all = []

    for i in range(len(temperature_all_models)):

        t_units = temperature_all_models[i] * units.kelvin
        p_units = pressure_all_models[i] * units.mbar

        theta_es = mpcalc.saturation_equivalent_potential_temperature(p_units, t_units)
        theta_es_all.append(theta_es.magnitude)

    print ("Length of theta es list (all models):",len(theta_es_all))

    return theta_es_all

def get_diff(ref_list,pressure_list):

    """
    Get difference between SCM RCE profiles and ensemble mean of SCM
    :return: 2 lists: diffs for all SCMs and common pressure list
    """

    # construct common pressure list

    pressure_start = 0
    pressure_end = 1000

    common_pressure_list = []

    for n in range(pressure_start, pressure_end + 1):

        if n % 50 == 0:
            common_pressure_list.append(n)

    lmdz6a_common = []
    lmdz6ab_common = []
    scam_common = []
    wrf_kf_common = []
    wrf_zm_common = []
    wrf_nt_common = []
    wrf_bmj_common = []
    wrf_nsas_common = []
    ummf_common = []
    umbm_common = []
    cnrm_common = []


    for n in range(len(common_pressure_list)):

        wrf_kf_common.append(np.interp(common_pressure_list[n], pressure_list[1], ref_list[1]))
        wrf_nt_common.append(np.interp(common_pressure_list[n], pressure_list[2], ref_list[2]))
        wrf_nsas_common.append(np.interp(common_pressure_list[n], pressure_list[3], ref_list[3]))
        wrf_bmj_common.append(np.interp(common_pressure_list[n], pressure_list[4], ref_list[4]))
        wrf_zm_common.append(np.interp(common_pressure_list[n], pressure_list[5], ref_list[5]))
        cnrm_common.append(np.interp(common_pressure_list[n], pressure_list[6], ref_list[6]))
        ummf_common.append(np.interp(common_pressure_list[n], pressure_list[7], ref_list[7]))
        umbm_common.append(np.interp(common_pressure_list[n], pressure_list[8], ref_list[8]))
        scam_common.append(np.interp(common_pressure_list[n], pressure_list[9], ref_list[9]))
        lmdz6a_common.append(np.interp(common_pressure_list[n], pressure_list[10], ref_list[10]))
        lmdz6ab_common.append(np.interp(common_pressure_list[n], pressure_list[11], ref_list[11]))

    zipped = zip(lmdz6a_common,lmdz6ab_common,scam_common,wrf_zm_common,wrf_kf_common,wrf_nt_common,wrf_nsas_common,wrf_bmj_common,umbm_common,ummf_common,cnrm_common)

    combined = []

    for z in zipped:
        z = list(z)
        combined.append(z)

    combined = np.array(combined)
    scm_average = np.mean(combined,axis=1)

    diff_all = []

    wrf_kf_diff = list(np.array(wrf_kf_common) - np.array(scm_average))
    diff_all.append(wrf_kf_diff)
    wrf_nt_diff = list(np.array(wrf_nt_common) - np.array(scm_average))
    diff_all.append(wrf_nt_diff)
    wrf_nsas_diff = list(np.array(wrf_nsas_common) - np.array(scm_average))
    diff_all.append(wrf_nsas_diff)
    wrf_bmj_diff = list(np.array(wrf_bmj_common) - np.array(scm_average))
    diff_all.append(wrf_bmj_diff)
    wrf_zm_diff = list(np.array(wrf_zm_common) - np.array(scm_average))
    diff_all.append(wrf_zm_diff)
    cnrm_diff = list(np.array(cnrm_common) - np.array(scm_average))
    diff_all.append(cnrm_diff)
    ummf_diff = list(np.array(ummf_common) - np.array(scm_average))
    diff_all.append(ummf_diff)
    umbm_diff = list(np.array(umbm_common) - np.array(scm_average))
    diff_all.append(umbm_diff)
    scam_diff = list(np.array(scam_common) - np.array(scm_average))
    diff_all.append(scam_diff)
    lmdz6a_diff = list(np.array(lmdz6a_common) - np.array(scm_average))
    diff_all.append(lmdz6a_diff)
    lmdz6ab_diff = list(np.array(lmdz6ab_common) - np.array(scm_average))
    diff_all.append(lmdz6ab_diff)

    print ()
    print ("Length of diff list (all models):",len(diff_all))
    print ()

    return diff_all,common_pressure_list


def plot_mean_states():

    """
    Main function to plot RCE T and RH profiles for all models
    """


    pressures_all = get_RH_data()[1]
    temperature_all = get_theta_es()
    rh_all = get_RH_data()[0]

    common_pressure_list = get_diff(temperature_all,pressures_all)[1]
    temperature_diff = get_diff(temperature_all,pressures_all)[0]
    rh_diff = get_diff(rh_all,pressures_all)[0]

    # For plotting zero line

    x_zero = []

    for n in range(0, len(common_pressure_list)):
        x_zero.append(0.0)



    #############################
    # PLOT TEMPERATURE
    #############################

    # plot parameters
    tick_fontsize = 16
    label_fontsize = 16
    legend_fontsize = 15

    scm_labels = ["WRF-KF", "WRF-NT", "WRF-NSAS", "WRF-BMJ", "WRF-ZM", "CNRM", "UM-MF", "UM-SBM", "SCAM", "LMDZ6A", "LMDZ6Ab"]

    y_label = "P [hPa]"

    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(12, 12), facecolor='w', edgecolor='k')
    fig.subplots_adjust(hspace=.5, wspace=.001)

    axs = axs.ravel()

    # ---------------------
    # PLOT SAM T

    axs[0].plot(temperature_all[0], pressures_all[0], label="SAM (CRM)", color="black", linewidth=2.5, zorder=10)
    axs[0].legend(loc="lower left", fontsize=legend_fontsize)

    axs[0].set_ylim([50, 1000])
    axs[0].set_yticks([100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
    axs[0].set_ylim(axs[0].get_ylim()[::-1])  # invert y axis when y axis is pressure


    axs[0].set_xlim([314, 366])
    axs[0].set_xticks([320,340,360])
    x_label = "$\it{θ_{es}}$ [K]"

    axs[0].tick_params(labelsize=tick_fontsize)

    axs[0].set_xlabel(x_label, fontsize=label_fontsize)
    axs[0].set_ylabel(y_label, fontsize=label_fontsize)

    axs[0].annotate("(a)", (296, 110), fontsize=label_fontsize, annotation_clip=False)

    plt.tight_layout()

    # ---------------------
    # PLOT SCM T

    for i in range(len(scm_labels)):

        if i == 9:
            axs[1].plot(temperature_all[i+1], pressures_all[i+1], label=scm_labels[i], color="#ffda33")

        else:
            axs[1].plot(temperature_all[i+1], pressures_all[i+1], label=scm_labels[i])


    axs[1].set_ylim([50, 1000])
    axs[1].set_yticks([100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
    axs[1].set_ylim(axs[1].get_ylim()[::-1])  # invert y axis when y axis is pressure


    axs[1].set_xlim([314, 366])  # for T'
    axs[1].set_xticks([320,340,360])
    x_label = "$\it{θ_{es}}$ [K]"

    axs[1].tick_params(labelsize=tick_fontsize)

    axs[1].set_xlabel(x_label, fontsize=label_fontsize)
    axs[1].set_ylabel(y_label, fontsize=label_fontsize)

    axs[1].annotate("(b)", (295, 110), fontsize=label_fontsize, annotation_clip=False)

    plt.tight_layout()

    # ---------------------
    # PLOT SCM DELTA T

    print()

    for i in range(len(scm_labels)):

        if i == 9:
            axs[2].plot(temperature_diff[i], common_pressure_list, label=scm_labels[i],color="#ffda33" )

        else:
            axs[2].plot(temperature_diff[i], common_pressure_list, label=scm_labels[i])

    axs[2].plot(x_zero, common_pressure_list, "black", linestyle="dashed", linewidth="0.8")

    axs[2].set_ylim([50, 1000])
    axs[2].set_yticks([100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
    axs[2].set_ylim(axs[2].get_ylim()[::-1])  # invert y axis when y axis is pressure

    axs[2].set_xlim([-24, 24])  # for T'
    axs[2].set_xticks([-20, -10, 0, 10, 20])
    x_label = "Δ$\it{θ_{es}}$ [K]"

    axs[2].tick_params(labelsize=tick_fontsize)

    axs[2].set_xlabel(x_label, fontsize=label_fontsize)
    axs[2].set_ylabel(y_label, fontsize=label_fontsize)

    axs[2].annotate("(c)", (-41, 110), fontsize=label_fontsize, annotation_clip=False)

    plt.tight_layout()


    #############################
    # PLOT RH
    #############################

    #---------------------
    # PLOT SAM RH

    axs[3].plot(rh_all[0], pressures_all[0], label="SAM (CRM)", color="black", linewidth=2.5, zorder=10)
    axs[3].legend(loc="lower left", fontsize=legend_fontsize)

    axs[3].set_ylim([50, 1000])
    axs[3].set_yticks([100,200, 300, 400, 500, 600, 700, 800, 900, 1000])
    axs[3].set_ylim(axs[3].get_ylim()[::-1])  # invert y axis when y axis is pressure


    axs[3].set_xlim([0, 100])
    axs[3].set_xticks([0, 25, 50, 75, 100])
    x_label = "RH [%]"

    axs[3].tick_params(labelsize=tick_fontsize)

    axs[3].set_xlabel(x_label,fontsize=label_fontsize)
    axs[3].set_ylabel(y_label,fontsize=label_fontsize)

    axs[3].annotate("(d)", (-35, 110), fontsize=label_fontsize, annotation_clip=False)

    plt.tight_layout()

    #---------------------
    # PLOT SCM RH

    for i in range(len(scm_labels)):

        if i == 9:
            axs[4].plot(rh_all[i + 1], pressures_all[i + 1], label=scm_labels[i], color="#ffda33")

        else:
            axs[4].plot(rh_all[i + 1], pressures_all[i + 1], label=scm_labels[i])

    axs[4].set_ylim([50, 1000])
    axs[4].set_yticks([100,200, 300, 400, 500, 600, 700, 800, 900, 1000])
    axs[4].set_ylim(axs[4].get_ylim()[::-1])  # invert y axis when y axis is pressure


    axs[4].set_xlim([0, 100])
    axs[4].set_xticks([0, 25, 50, 75, 100])
    x_label = "RH [%]"

    axs[4].tick_params(labelsize=tick_fontsize)

    axs[4].set_xlabel(x_label,fontsize=label_fontsize)
    axs[4].set_ylabel(y_label, fontsize=label_fontsize)

    axs[4].annotate("(e)", (-35, 110), fontsize=label_fontsize, annotation_clip=False)

    plt.tight_layout()

    #---------------------
    # PLOT SCM DELTA RH

    for i in range(len(scm_labels)):

        if i == 9:
            axs[5].plot(rh_diff[i], common_pressure_list, label=scm_labels[i], color="#ffda33")

        else:
            axs[5].plot(rh_diff[i], common_pressure_list, label=scm_labels[i])

    axs[5].plot(x_zero, common_pressure_list, "black", linestyle="dashed", linewidth="0.8")

    axs[5] = plt.gca()

    axs[5].set_ylim([50, 1000])
    axs[5].set_yticks([100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
    axs[5].set_ylim(axs[5].get_ylim()[::-1])  # invert y axis when y axis is pressure

    axs[5].set_xlim([-60, 60])  # for T'
    axs[5].set_xticks([-50, -25, 0, 25, 50])
    x_label = "ΔRH [%]"

    axs[5].tick_params(labelsize=tick_fontsize)

    axs[5].set_xlabel(x_label, fontsize=label_fontsize)
    axs[5].set_ylabel(y_label, fontsize=label_fontsize)

    axs[5].annotate("(f)", (-100, 110), fontsize=label_fontsize, annotation_clip=False)

    plt.tight_layout()

    # COMMON legend
    fig.subplots_adjust(top=0.86, hspace=0.2, wspace=0.4)
    axs.flatten()[2].legend(bbox_to_anchor=(0.4, 1.34), ncol=4, fontsize=legend_fontsize)

    plt.show()


if __name__ == "__main__":

    plot_mean_states()



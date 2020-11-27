#!/usr/bin/env python3

"""
Created on 13 August 2020

@author: Yi-Ling HWONG <yiling.hwong@gmail.com>

This script plots the correlation matrix for RCE RH vs. q response to dT/dt perturbation

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

"""
Select parameters to plot here
"""

perturbation = "q_dtdt"
use_only_magnitude_of_responses = True # set to True if ignore signs of responses and use only magnitude


class MidpointNormalize(colors.Normalize):
    """
    Normalise the colorbar so that diverging bars work their way either side from a prescribed midpoint value)

    e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
    """

    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))


def get_common_pressure_list():

    # construct common pressure list (250,350, ... 950 hPa)

    pressure_start = 200
    pressure_end = 950

    pressure_list = []

    for n in range(pressure_start,pressure_end+1):

        if n%100 == 0:
            pressure_list.append(n+50)

    return pressure_list

def  get_rh_all_models():

    lines0 = open("../../data/SAM/REF/rh.csv","r").readlines()
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
    lines10 = open("../../data/LMDZ/REF/lmdz6ab/RCE_mean_RH.csv", "r").readlines()
    lines11 = open("../../data/LMDZ/REF/lmdz6a/RCE_mean_RH.csv", "r").readlines()
    lines12 = open("../../data/LMDZ/REF/lmdz5a/RCE_mean_RH.csv", "r").readlines()

    sam_ref = []
    sam_pres = []
    lmdz5a_ref = []
    lmdz5a_pres = []
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

    # 10. LMDZ6Ab
    for line in lines10:
        spline = line.rstrip("\n").split(",")
        lmdz6ab_pres.append(float(spline[0]))
        lmdz6ab_ref.append(float(spline[1]))

    if lmdz6ab_pres[0] > lmdz6ab_pres[1]:
        lmdz6ab_pres.reverse()
        lmdz6ab_ref.reverse()

    ref_all.append(lmdz6ab_ref)
    pressures_all.append(lmdz6ab_pres)

    # 11. LMDZ6A
    for line in lines11:
        spline = line.rstrip("\n").split(",")
        lmdz6a_pres.append(float(spline[0]))
        lmdz6a_ref.append(float(spline[1]))

    if lmdz6a_pres[0] > lmdz6a_pres[1]:
        lmdz6a_pres.reverse()
        lmdz6a_ref.reverse()

    ref_all.append(lmdz6a_ref)
    pressures_all.append(lmdz6a_pres)

    # 12. LMDZ5A
    for line in lines12:
        spline = line.rstrip("\n").split(",")
        lmdz5a_pres.append(float(spline[0]))
        lmdz5a_ref.append(float(spline[1]))

    if lmdz5a_pres[0] > lmdz5a_pres[1]:
        lmdz5a_pres.reverse()
        lmdz5a_ref.reverse()

    # ref_all.append(lmdz5a_ref)
    # pressures_all.append(lmdz5a_pres)

    # print()
    # print("Length of RCE RH list (all models):", len(ref_all), len(pressures_all))
    # print()

    return ref_all, pressures_all

def get_mean_RH_between_two_pressure_levels(pressure_list, rh_list):


    if len(rh_list) == 1:

        mean_rh = rh_list[0]

    else:

        middle_values = [] # midpoint of all individual lines y=mx

        for index,value in enumerate(pressure_list):

            if index == len(pressure_list) - 1:

                break

            else:

                m = abs((rh_list[index] - rh_list[index+1]) / (pressure_list[index] - pressure_list[index+1])) # gradient of slope

                half_value_pressure = abs(((pressure_list[index] - pressure_list[index+1]) / 2))
                mid_value_rh = min(rh_list[index], rh_list[index+1]) + (m*(half_value_pressure))

                middle_values.append(mid_value_rh)


        mean_rh = np.average(middle_values)

    return mean_rh

def get_mean_anomalies():

    """
    Get mean anomalies across a horizontal stripe in M^-1 matrices
    """

    #perturbation = "q_dtdt"

    folder_list = ["SAM","WRF","WRF","WRF","WRF","WRF","CNRM","UMMF","UMBM","SCAM","LMDZ","LMDZ"]
    model_list = ["sam","wrf_kfeta","wrf_ntiedtke","wrf_nsas","wrf_bmj","wrf_camzm","cnrm","ummf","umbm","scam","lmdz6ab","lmdz6a"]

    matrix_all = []

    for i in range(len(folder_list)):

        lines = open("../../data/"+folder_list[i]+"/matrix_M_inv/M_inv_"+model_list[i]+"_"+perturbation+"_norm_kuang.csv","r").readlines()

        matrix = []

        for line in lines:
            spline = line.rstrip("\n").split(",")
            temp = [float(x) for x in spline]

            if use_only_magnitude_of_responses == True:
                temp = [abs(float(x)) for x in spline]

            matrix.append(temp)

        matrix_all.append(matrix)

    matrix_mean_across_horizontal_stripes = []

    for n in range(len(matrix_all)):

        matrix_mean = (np.mean(matrix_all[n],axis=1)).tolist()
        matrix_mean.reverse()
        matrix_mean_across_horizontal_stripes.append(matrix_mean)

    print ("Length of matrix mean list (all models):",len(matrix_mean_across_horizontal_stripes),len(matrix_mean_across_horizontal_stripes[0]))

    return matrix_mean_across_horizontal_stripes


def get_rh_at_common_pressure_levels():

    common_pressure_levels = get_common_pressure_list()
    common_pressure_levels.reverse()

    ref_all,pressures_all = get_rh_all_models()

    rh_means_at_common_pressure_levels = []

    for i in range(len(common_pressure_levels)):

        p2 = common_pressure_levels[i] + 50
        p1 = common_pressure_levels[i] - 50

        ref_mean_all = []

        for n in range(len(pressures_all)):

            model_ref = ref_all[n]
            model_pres = pressures_all[n]

            model_ref_list = []
            model_pres_list = []

            for index, value in enumerate(model_pres):

                if value >= p1:
                    if value <= p2:
                        model_ref_list.append(model_ref[index])
                        model_pres_list.append(model_pres[index])

            model_rh_mean = get_mean_RH_between_two_pressure_levels(model_pres_list, model_ref_list)
            ref_mean_all.append(model_rh_mean)

        rh_means_at_common_pressure_levels.append(ref_mean_all)

    print ()
    print ("Length of RH means at common pressure levels:",len(rh_means_at_common_pressure_levels))
    print ("Number of models at each common pressure level:",len(rh_means_at_common_pressure_levels[0]))

    return rh_means_at_common_pressure_levels

def get_pressure_levels_kuang_all_models():

    """
    Get the 18 pressure levels of M^-1 matrix for all models
    """

    folder_list = ["SAM","WRF","WRF","WRF","WRF","WRF","CNRM","UMMF","UMBM","SCAM","LMDZ","LMDZ"]
    model_list = ["sam","kfeta","ntiedtke","nsas","bmj","camzm","cnrm","ummf","umbm","scam","lmdz","lmdz"]

    pressures_kuang_all = []

    for n in range(len(folder_list)):

        if folder_list[n] == "WRF":
            lines = open("../../data/"+folder_list[n]+"/matrix_X_raw/"+model_list[n]+"/pressures_kuang","r").readlines()
        else:
            lines = open("../../data/"+folder_list[n]+"/matrix_X_raw/pressures_kuang","r").readlines()

        pressures = []

        for line in lines:
            spline = line.rstrip("\n")
            pressures.append(float(spline) / 100)

        pressures.reverse()
        pressures_kuang_all.append(pressures)

    print ("Length of kuang pressure list (all models):",len(pressures_kuang_all),len(pressures_kuang_all[0]))

    return pressures_kuang_all


def get_anomalies_at_common_pressure_levels():

    common_pressure_levels = get_common_pressure_list()
    anom_means_all = get_mean_anomalies() # mean across horizontal stripes for all models
    pressures_kuang_all = get_pressure_levels_kuang_all_models()

    anom0 = []
    anom1 = []
    anom2 = []
    anom3 = []
    anom4 = []
    anom5 = []
    anom6 = []
    anom7 = []
    anom8 = []
    anom9 = []
    anom10 = []
    anom11 = []
    anom12 = []


    for n in range(len(common_pressure_levels)):

        anom0.append(np.interp(common_pressure_levels[n], pressures_kuang_all[0], anom_means_all[0]))
        anom1.append(np.interp(common_pressure_levels[n], pressures_kuang_all[1], anom_means_all[1]))
        anom2.append(np.interp(common_pressure_levels[n], pressures_kuang_all[2], anom_means_all[2]))
        anom3.append(np.interp(common_pressure_levels[n], pressures_kuang_all[3], anom_means_all[3]))
        anom4.append(np.interp(common_pressure_levels[n], pressures_kuang_all[4], anom_means_all[4]))
        anom5.append(np.interp(common_pressure_levels[n], pressures_kuang_all[5], anom_means_all[5]))
        anom6.append(np.interp(common_pressure_levels[n], pressures_kuang_all[6], anom_means_all[6]))
        anom7.append(np.interp(common_pressure_levels[n], pressures_kuang_all[7], anom_means_all[7]))
        anom8.append(np.interp(common_pressure_levels[n], pressures_kuang_all[8], anom_means_all[8]))
        anom9.append(np.interp(common_pressure_levels[n], pressures_kuang_all[9], anom_means_all[9]))
        anom10.append(np.interp(common_pressure_levels[n], pressures_kuang_all[10], anom_means_all[10]))
        anom11.append(np.interp(common_pressure_levels[n], pressures_kuang_all[11], anom_means_all[11]))
        #anom12.append(np.interp(common_pressure_levels[n], pressures_kuang_all[12], anom_means_all[12]))

    print ()
    print ("COMMON PRESSURE LEVELS (hPa):",common_pressure_levels)
    print()
    print("ANOM MEANS AT COMMON PRESSURE LEVELS:")
    print("SAM = ", anom0)
    print("WRF_KF = ", anom1)
    print("WRF_NT = ", anom2)
    print("WRF_NSAS = ", anom3)
    print("WRF_BMJ = ", anom4)
    print("WRF_ZM = ", anom5)
    print("CNRM = ", anom6)
    print("UM_MF = ", anom7)
    print("UM_SBM = ", anom8)
    print("SCAM = ", anom9)
    print("LMDZ6Ab = ", anom10)
    print("LMDZ6A = ", anom11)
    #print("LMDZ5A = ", anom12)
    print()

    return anom0,anom1,anom2,anom3,anom4,anom5,anom6,anom7,anom8,anom9,anom10,anom11

def get_correlations():

    common_pressure_levels = get_common_pressure_list()
    rh_at_common_pressure_levels = get_rh_at_common_pressure_levels()
    anom0,anom1,anom2,anom3,anom4,anom5,anom6,anom7,anom8,anom9,anom10,anom11 = get_anomalies_at_common_pressure_levels()

    from scipy.stats import pearsonr, spearmanr

    zipped = zip(anom0,anom1,anom2,anom3,anom4,anom5,anom6,anom7,anom8,anom9,anom10,anom11)

    anom_zipped = []

    for z in zipped:
        z = list(z)
        anom_zipped.append(z)

    corr_list_all = []
    p_values_all = []

    for i in range(len(common_pressure_levels)):

        ref_mean_all = rh_at_common_pressure_levels[i]

        correlation_list = []
        p_values = []

        for n in range(len(common_pressure_levels)):

            #correlation,p_value = pearsonr(ref_mean_all,anom_zipped[n])
            correlation, p_value = spearmanr(ref_mean_all, anom_zipped[n])

            p = str(round(p_value,3))

            p_values.append(p)

            #-----------------

            x = str(round(correlation,2))
            correlation_list.append(x)

        correlation_list = [float(x) for x in correlation_list]

        corr_list_all.append(correlation_list)
        p_values_all.append(p_values)

    # Transpose array
    corr_list_all = np.array(corr_list_all)
    corr_list_all = corr_list_all.transpose()

    p_values_all = np.array(p_values_all)
    p_values_all = p_values_all.transpose()

    return corr_list_all,p_values_all

def plot_corr_matrix():

    import pandas as pd
    import seaborn as sn

    corr_list,p_values = get_correlations()

    df_corr = pd.DataFrame(corr_list)
    df_pvals = pd.DataFrame(p_values)

    print ()
    print ("######################################")
    print ()
    print ("CORRELATION MATRIX:")
    print(df_corr)
    print ()
    print ("P VALUES:")
    print (df_pvals)

    #############################
    # PLOT
    #############################

    sn.set_style("whitegrid")  # set background to white

    x_labels = ["950", "850", "750", "650", "550", "450", "350", "250"]
    y_labels = ["250", "350", "450", "550", "650", "750", "850", "950"]

    cbar_ticks = [-1, -0.75, -0.5, -0.25, 0.00, 0.25, 0.5, 0.75, 1.0]

    g = sn.heatmap(df_corr, annot=True, cmap="RdBu_r", square=True, vmin=-1, vmax=1,
                   norm=MidpointNormalize(midpoint=0.0, vmin=-1, vmax=1), xticklabels=x_labels, yticklabels=y_labels,
                   cbar_kws=dict(ticks=cbar_ticks))  # camp = RdBu_r

    g.yaxis.set_ticks_position("left")  # set ticks on left (otherwise will have ticks on right as well)
    g.xaxis.set_ticks_position("bottom")

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    plt.xlabel("Pressure level of RH [hPa]", fontsize=16, labelpad=12)
    plt.ylabel("Pressure level of q' [hPa]", fontsize=16, labelpad=12)

    plt.tight_layout()

    #plt.savefig('/Users/yi-linghwong/Documents/postdoc/publications/paracon_scm_paper/_PLOTS/_FINAL/corr_RH_vs_Q_DTDT_magnitude_only.png',dpi=300)
    #plt.savefig('/Users/yi-linghwong/Documents/postdoc/publications/paracon_scm_paper/_PLOTS/_FINAL/corr_RH_vs_Q_DTDT.png',dpi=300)

    plt.show()


if __name__ == "__main__":

    plot_corr_matrix()
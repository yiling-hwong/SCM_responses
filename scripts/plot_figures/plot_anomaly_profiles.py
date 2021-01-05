#!/usr/bin/env python3

"""
Created on 13 August 2020

@author: Yi-Ling HWONG <yiling.hwong@gmail.com>

This script plots the T and q response vertical profiles for all models or selected models.
These profiles are extracted columns from the M_inv matrix

"""

"""
Set parameters here
"""

plot_all_profiles = False
plot_selected_profiles = True

# set one of the following to True (perturb dT/dt or dq/dt)
# Only used if plot_all_profiles = True
perturb_t = True
perturb_q = False

# set state anomaly to either "T" or "q"
state_anomaly = "T"

# select index of model to plot (listed below)
# used only if plot_selected_profiles = True
if state_anomaly == "T":
    model_index_to_plot = [0,1,10,7] # SAM, CNRM + UM-MF for T'
elif state_anomaly == "q":
    model_index_to_plot = [0,1,4,6] # SAM, WRF-BMJ and WRF-NSAS for q'

# 0=perturbation shape
# 1=SAM
# 2=WRF-KF
# 3=WRF-NT
# 4=WRF-BMJ
# 5=WRF-ZM
# 6=WRF-NSAS
# 7=UM-MF
# 8=UM-SBM
# 9=SCAM
# 10=CNRM
# 11=LMDZ5A
# 12=LMDZ6A
# 13=LMDZ6Ab

#-------------------------

import matplotlib.pyplot as plt

folder_list = ["SAM", "WRF", "WRF", "WRF", "WRF", "WRF", "UMMF", "UMBM", "SCAM", "CNRM", "LMDZ", "LMDZ", "LMDZ"]
model_list = ["SAM", "kfeta", "ntiedtke", "bmj", "camzm", "nsas", "UMMF", "UMBM", "SCAM", "CNRM", "LMDZ5A", "LMDZ6A",
              "LMDZ6Ab"]

def get_response_data_per_quadrat():

    """
    Get responses for all model for either one of the four quadrats
    i.e. T_DTDT or Q_DTDT or T_DQDT or Q_DQDT
    """

    if perturb_t == True:

        if state_anomaly == "T":
            perturbation = "T_DTDT"
        elif state_anomaly == "q":
            perturbation = "Q_DTDT"

    elif perturb_q == True:

        if state_anomaly == "T":
            perturbation = "T_DQDT"
        elif state_anomaly == "q":
            perturbation = "Q_DQDT"

    pressures_all = []
    anom_850_all = []
    anom_650_all = []


    if perturbation == "T_DTDT" or perturbation == "Q_DTDT":
        lines = open("../../data/perturbation_profiles/DTDT_profiles.csv","r").readlines()
    if perturbation == "T_DQDT" or perturbation == "Q_DQDT":
        lines = open("../../data/perturbation_profiles/DQDT_profiles.csv","r").readlines()

    pressures_perturb = []
    perturb_850 = []
    perturb_650 = []

    for line in lines[1:]:
        spline = line.rstrip("\n").split(",")
        pressures_perturb.append(float(spline[0]))
        perturb_850.append(float(spline[1]))
        perturb_650.append(float(spline[2]))

    pressures_all.append(pressures_perturb)
    anom_850_all.append(perturb_850)
    anom_650_all.append(perturb_650)

    for n in range(len(folder_list)):

        if perturb_t == True:
            if model_list[n] == "SCAM" or model_list[n] == "UMMF" or model_list[n] == "LMDZ6Ab":
                perturbation_amplitude = "02"
            else:
                perturbation_amplitude = "05"
        elif perturb_q == True:
            if model_list[n] == "SCAM" or model_list[n] == "UMMF" or model_list[n] == "LMDZ6Ab":
                perturbation_amplitude = "01"
            else:
                perturbation_amplitude = "02"

        lines = open("../../data/"+folder_list[n]+"/anomalies/normalised/"+model_list[n]+"_"+perturbation+"_"+perturbation_amplitude+"_norm_kuang.csv","r").readlines()

        pressures = []
        anom_850 = []
        anom_650 = []

        for line in lines[1:]:
            spline = line.rstrip("\n").split(",")
            pressures.append(float(spline[0]))
            anom_850.append(float(spline[1]))
            anom_650.append(float(spline[2]))

        pressures_all.append(pressures)
        anom_850_all.append(anom_850)
        anom_650_all.append(anom_650)

    return pressures_all,anom_850_all,anom_650_all

def get_responses_data_for_state_anomaly():

    """
    Get responses for all model for either one of the state anomalies
    i.e. either "T" or "q" for both perturbations (dT/dt and dq/dt)
    """

    pressures_dtdt = []
    pressures_dqdt = []
    anom_850_dtdt = []
    anom_650_dtdt = []
    anom_850_dqdt = []
    anom_650_dqdt = []

    lines0a = open("../../data/perturbation_profiles/DTDT_profiles.csv", "r").readlines()
    lines0b = open("../../data/perturbation_profiles/DQDT_profiles.csv", "r").readlines()

    pressures_perturb_dtdt = []
    pressures_perturb_dqdt = []
    perturb_850_dtdt = []
    perturb_650_dtdt = []
    perturb_850_dqdt = []
    perturb_650_dqdt = []

    for line in lines0a[1:]:
        spline = line.rstrip("\n").split(",")
        pressures_perturb_dtdt.append(float(spline[0]))
        perturb_850_dtdt.append(float(spline[1]))
        perturb_650_dtdt.append(float(spline[2]))

    for line in lines0b[1:]:
        spline = line.rstrip("\n").split(",")
        pressures_perturb_dqdt.append(float(spline[0]))
        perturb_850_dqdt.append(float(spline[1]))
        perturb_650_dqdt.append(float(spline[2]))

    pressures_dtdt.append(pressures_perturb_dtdt)
    pressures_dqdt.append(pressures_perturb_dqdt)
    anom_850_dtdt.append(perturb_850_dtdt)
    anom_850_dqdt.append(perturb_850_dqdt)
    anom_650_dtdt.append(perturb_650_dtdt)
    anom_650_dqdt.append(perturb_650_dqdt)

    # DTDT
    for n in range(len(folder_list)):

        if model_list[n] == "SCAM" or model_list[n] == "UMMF" or model_list[n] == "LMDZ6Ab":
            lines = open("../../data/" + folder_list[n] + "/anomalies/normalised/" + model_list[n] + "_" + state_anomaly + "_DTDT_02_norm_kuang.csv", "r").readlines()
        else:
            lines = open("../../data/" + folder_list[n] + "/anomalies/normalised/" + model_list[n] + "_" + state_anomaly + "_DTDT_05_norm_kuang.csv", "r").readlines()

        pressures = []
        anom_850 = []
        anom_650 = []

        for line in lines[1:]:
            spline = line.rstrip("\n").split(",")
            pressures.append(float(spline[0]))
            anom_850.append(float(spline[1]))
            anom_650.append(float(spline[2]))

        pressures_dtdt.append(pressures)
        anom_850_dtdt.append(anom_850)
        anom_650_dtdt.append(anom_650)

        # DQDT
        for n in range(len(folder_list)):

            if model_list[n] == "SCAM" or model_list[n] == "UMMF" or model_list[n] == "LMDZ6Ab":
                lines = open("../../data/" + folder_list[n] + "/anomalies/normalised/" + model_list[n] + "_" + state_anomaly + "_DQDT_01_norm_kuang.csv","r").readlines()
            else:
                lines = open("../../data/" + folder_list[n] + "/anomalies/normalised/" + model_list[n] + "_" + state_anomaly + "_DQDT_02_norm_kuang.csv", "r").readlines()

            pressures = []
            anom_850 = []
            anom_650 = []

            for line in lines[1:]:
                spline = line.rstrip("\n").split(",")
                pressures.append(float(spline[0]))
                anom_850.append(float(spline[1]))
                anom_650.append(float(spline[2]))

            pressures_dqdt.append(pressures)
            anom_850_dqdt.append(anom_850)
            anom_650_dqdt.append(anom_650)

    return pressures_dtdt,anom_850_dtdt,anom_650_dtdt,pressures_dqdt,anom_850_dqdt,anom_650_dqdt


def plot_responses_all():

    pressures_all,anom_850_all,anom_650_all = get_response_data_per_quadrat()

    # plot parameters
    tick_fontsize = 12
    label_fontsize = 12
    legend_fontsize = 12
    title_fontsize = 13

    # For zero line
    x_zero = []

    for n in range(0,len(pressures_all[0])):
        x_zero.append(0.0)

    #############################
    # PLOT
    #############################

    print ()
    print ("Plotting responses for all models...")

    plot_titles = ["(a) Perturbation","(b) SAM (CRM)","(c) WRF-KF","(d) WRF-NT","(e) WRF-BMJ",
                   "(f) WRF-ZM","(g) WRF-NSAS","(h) UM-MF","(i) UM-SBM","(j) SCAM","(k) CNRM","(l) LMDZ5A","(m) LMDZ6A","(n) LMDZ6Ab"]

    fig, axs = plt.subplots(nrows=4, ncols=4, figsize=(12, 16), facecolor='w', edgecolor='k')
    fig.subplots_adjust(hspace=.5, wspace=.001)
    fig.delaxes(axs[3, 2])
    fig.delaxes(axs[3, 3])

    axs = axs.ravel()

    y_label = "P [hPa]"

    for i in range(len(pressures_all)):

        axs[i].plot(anom_850_all[i], pressures_all[i], "r", label="850 hPa")
        axs[i].plot(anom_650_all[i], pressures_all[i], "b", label="650 hPa")
        axs[i].plot(x_zero, pressures_all[0], "black", linestyle="dashed", linewidth="0.8")

        axs[i].set_ylim([200, 1000])
        axs[i].set_yticks([200,300,400,500,600,700,800,900,1000])
        axs[i].set_ylim(axs[i].get_ylim()[::-1])  # invert y axis when y axis is pressure

        if i == 0:

            axs[i].legend(loc="upper right", fontsize=legend_fontsize)

            if perturb_t == True:
                axs[i].set_xlim([-0.04, 0.54])  # for T'
                axs[i].set_xticks([0.0,0.1,0.2,0.3,0.4,0.5])
                x_label = "dT/dt [K d$\mathregular{^{-1}}$]"

            if perturb_q == True:
                axs[i].set_xlim([-0.015, 0.215])  # for T'
                axs[i].set_xticks([0.0,0.1,0.2])
                x_label = "dq/dt [g kg$\mathregular{^{-1}}$ d$\mathregular{^{-1}}$]"

        else:

            if state_anomaly == "T":
                axs[i].set_xlim([-0.35, 1.05])  # for T'
                axs[i].set_xticks([-0.3, 0.0, 0.3, 0.6, 0.9])
                x_label = "T' [K]"

            if state_anomaly == "q":
                axs[i].set_xlim([-0.3, 0.5])  # for q'
                axs[i].set_xticks([-0.2, 0.0, 0.2, 0.4])
                x_label = "q' [g kg$\mathregular{^{-1}}$]"

        axs[i].tick_params(labelsize=tick_fontsize)

        axs[i].set_xlabel(x_label, fontsize=label_fontsize)
        axs[i].set_ylabel(y_label, fontsize=label_fontsize)

        axs[i].set_title(plot_titles[i], fontsize=title_fontsize)

    plt.tight_layout()

    plt.show()

def plot_responses_selected():

    pressures_dtdt,anom_850_dtdt,anom_650_dtdt,pressures_dqdt,anom_850_dqdt,anom_650_dqdt = get_responses_data_for_state_anomaly()

    # plot parameters
    tick_fontsize = 19
    label_fontsize = 19
    legend_fontsize = 19
    title_fontsize = 22

    # For zero line
    x_zero = []

    for n in range(0,len(pressures_dtdt[0])):
        x_zero.append(0.0)

    #############################
    # PLOT
    #############################

    print ()
    print ("Plotting responses for selected models...")

    plot_titles = ["Perturbation","SAM (CRM)","WRF-KF","WRF-NT","WRF-BMJ",
                   "WRF-ZM","WRF-NSAS","UM-MF","UM-SBM","SCAM","CNRM","LMDZ5A","LMDZ6A","LMDZ6Ab"]

    fig, axs = plt.subplots(nrows=2, ncols=4, figsize=(18, 13), facecolor='w', edgecolor='k')
    fig.subplots_adjust(hspace=0.5, wspace=.001)

    axs = axs.ravel()

    # DTDT
    label_list = ["(a) dT/dt ","(b) ", "(c) ", "(d) "]
    for i in range(len(model_index_to_plot)):

        axs[i].plot(anom_850_dtdt[model_index_to_plot[i]], pressures_dtdt[model_index_to_plot[i]], "r", label="850 hPa")
        axs[i].plot(anom_650_dtdt[model_index_to_plot[i]], pressures_dtdt[model_index_to_plot[i]], "b", label="650 hPa")
        axs[i].plot(x_zero, pressures_dtdt[0], "black", linestyle="dashed", linewidth="0.8")

        y_label = "P [hPa]"

        axs[i].set_ylim([200, 1000])
        axs[i].set_yticks([200,300,400,500,600,700,800,900,1000])
        axs[i].set_ylim(axs[i].get_ylim()[::-1])  # invert y axis when y axis is pressure

        if i == 0:
            axs[i].set_xlim([-0.04, 0.54])  # for T'
            axs[i].set_xticks([0.0,0.1,0.2,0.3,0.4,0.5])
            x_label = "K per day"

        else:

            if state_anomaly == "T":
                axs[i].set_xlim([-0.35, 1.05])  # for T'
                axs[i].set_xticks([-0.3, 0.0, 0.3, 0.6, 0.9])
                x_label = "T' [K]"

            if state_anomaly == "q":
                axs[i].set_xlim([-0.3, 0.5])  # for q'
                axs[i].set_xticks([-0.2, 0.0, 0.2, 0.4])
                x_label = "q' [g kg$\mathregular{^{-1}}$]"

        axs[i].tick_params(labelsize=tick_fontsize)

        axs[i].set_xlabel(x_label, fontsize=label_fontsize)

        plot_title = label_list[i]+plot_titles[model_index_to_plot[i]]

        axs[i].set_title(plot_title, fontsize=title_fontsize)

        if i != 0:
            empty_string_labels = [''] * 9
            axs[i].set_yticklabels(empty_string_labels)

        if i == 0:
            axs[i].legend(loc="upper right",fontsize=legend_fontsize)
            axs[i].set_ylabel(y_label, fontsize=label_fontsize)


    # DQDT
    label_list = ["(e) dq/dt ", "(f) ", "(g) ", "(h) "]
    for i in range(len(model_index_to_plot)):

        axs[i+4].plot(anom_850_dqdt[model_index_to_plot[i]], pressures_dqdt[model_index_to_plot[i]], "r",label="850 hPa")
        axs[i+4].plot(anom_650_dqdt[model_index_to_plot[i]], pressures_dqdt[model_index_to_plot[i]], "b",label="650 hPa")
        axs[i+4].plot(x_zero, pressures_dtdt[0], "black", linestyle="dashed", linewidth="0.8")

        y_label = "P [hPa]"

        axs[i+4].set_ylim([200, 1000])
        axs[i+4].set_yticks([200, 300, 400, 500, 600, 700, 800, 900, 1000])
        axs[i+4].set_ylim(axs[i+4].get_ylim()[::-1])  # invert y axis when y axis is pressure

        if i == 0:
            axs[i+4].set_xlim([-0.015, 0.215])  # for T'
            axs[i+4].set_xticks([0.0, 0.1, 0.2])
            x_label = "g kg$\mathregular{^{-1}}$ per day"

        else:

            if state_anomaly == "T":
                axs[i+4].set_xlim([-0.35, 1.05])  # for T'
                axs[i+4].set_xticks([-0.3, 0.0, 0.3, 0.6, 0.9])
                x_label = "T' [K]"

            if state_anomaly == "q":
                axs[i+4].set_xlim([-0.3, 0.5])  # for q'
                axs[i+4].set_xticks([-0.2, 0.0, 0.2, 0.4])
                x_label = "q' [g kg$\mathregular{^{-1}}$]"


        axs[i+4].tick_params(labelsize=tick_fontsize)

        axs[i+4].set_xlabel(x_label, fontsize=label_fontsize)

        plot_title = label_list[i]+plot_titles[model_index_to_plot[i]]

        axs[i+4].set_title(plot_title, fontsize=title_fontsize)

        if i != 0:
            empty_string_labels = [''] * 9
            axs[i+4].set_yticklabels(empty_string_labels)

        if i == 0:
            axs[i+4].legend(loc="upper right",fontsize=legend_fontsize)
            axs[i+4].set_ylabel(y_label, fontsize=label_fontsize)

    plt.tight_layout()

    plt.show()

def main():

    if plot_all_profiles == True:
        plot_responses_all()

    if plot_selected_profiles == True:
        plot_responses_selected()


if __name__ == "__main__":

    main()

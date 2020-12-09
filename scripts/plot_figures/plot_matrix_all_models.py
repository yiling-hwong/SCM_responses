#!/usr/bin/env python3

"""
Created on 13 August 2020

@author: Yi-Ling HWONG <yiling.hwong@gmail.com>

This script plots the M^-1 matrix for all models

"""

"""
Select quadrat to plot here
i.e. T'_DTDT, Q'_DTDT, T'_DQDT, or Q'_DQDT
"""

# set one of the following to True (perturb dT/dt or dq/dt)
perturb_t = True
perturb_q = False

# set state anomaly to either "T" or "q"
state_anomaly = "q"

#-------------------------


import matplotlib.pyplot as plt
import numpy as np

def get_matrix_data_all():

    if perturb_t == True:
        if state_anomaly == "T":
            perturbation = "t_dtdt"
        elif state_anomaly == "q":
            perturbation = "q_dtdt"

    elif perturb_q == True:
        if state_anomaly == "T":
            perturbation = "t_dqdt"
        elif state_anomaly == "q":
            perturbation = "q_dqdt"

    folder_list = ["SAM","WRF","WRF","WRF","WRF","WRF","CNRM","UMMF","UMBM","SCAM","LMDZ","LMDZ"]
    model_list = ["sam","kfeta","ntiedtke","nsas","bmj","camzm","cnrm","ummf","umbm","scam","lmdz6a","lmdz6ab"]

    pressures_all = []
    matrix_all = []

    for n in range(len(folder_list)):

        if folder_list[n] == "WRF":
            lines1 = open("../../data/"+folder_list[n]+"/matrix_X_raw/"+model_list[n]+"/pressures_kuang","r").readlines()
            lines2 = open("../../data/"+folder_list[n]+"/matrix_M_inv/M_inv_wrf_"+model_list[n]+"_"+perturbation+"_norm_kuang.csv","r").readlines()

        else:
            lines1 = open("../../data/"+folder_list[n]+"/matrix_X_raw/pressures_kuang","r").readlines()
            lines2 = open("../../data/"+folder_list[n]+"/matrix_M_inv/M_inv_"+model_list[n]+"_"+perturbation+"_norm_kuang.csv","r").readlines()

        pressures = []
        anomalies = []

        for line in lines1:
            spline = line.rstrip("\n")
            pressures.append(float(spline))

        pressures_all.append(pressures)

        for line in lines2:
            spline = line.rstrip("\n").split(",")
            spline = [float(s) for s in spline]
            anomalies.append(spline)

        matrix_all.append(anomalies)

    print ()
    print ("Length of pressures_all:",len(pressures_all))
    print ("Length of matrix_all:",len(matrix_all))

    return pressures_all,matrix_all

def plot_matrix_all():

    pressures_all,matrix_all = get_matrix_data_all()

    # Plot parameters
    tick_fontsize = 14
    label_fontsize = 14
    title_fontsize = 16

    if state_anomaly == "T":
        vmax_ = 0.5
        vmin_ = -vmax_
    if state_anomaly == "q":
        vmax_ = 0.3
        vmin_ = -vmax_

    #############################
    # PLOT
    #############################

    print ()
    print ("Plotting all matrices...")

    plot_titles = ["(a) SAM (CRM)","(b) WRF-KF","(c) WRF-NT","(d) WRF-NSAS","(e) WRF-BMJ",
                   "(f) WRF-ZM","(g) CNRM","(h) UM-MF","(i) UM-SBM","(j) SCAM","(k) LMDZ6A", "(l) LMDZ6Ab"]

    fig, axs = plt.subplots(nrows=3, ncols=4, figsize=(18, 13), facecolor='w', edgecolor='k')

    axs = axs.ravel()

    for i in range(len(pressures_all)):

        matrix = matrix_all[i]
        matrix = np.array(matrix)
        pressures = pressures_all[i]

        X1, Y1 = np.meshgrid(pressures, pressures)

        im1 = axs[i].pcolor(X1, Y1, matrix, vmin=vmin_, vmax=vmax_, cmap="rainbow")

        axs[i].set_ylim(axs[i].get_ylim()[::-1])  # invert y axis when y axis is pressure
        axs[i].set_xlim(axs[i].get_xlim()[::-1])

        xy_label = "P [hPa]"
        axs[i].set_xlabel(xy_label, fontsize=label_fontsize)
        axs[i].set_ylabel(xy_label, fontsize=label_fontsize)

        axs[i].set_xlim([1000, 200])
        axs[i].set_xticks([1000, 800, 600, 400, 200])

        axs[i].set_ylim([1000, 200])
        axs[i].set_yticks([1000, 800, 600, 400, 200])

        axs[i].tick_params(labelsize=tick_fontsize)
        axs[i].set_title(plot_titles[i], fontsize=title_fontsize)

        plt.tight_layout()  # so that the two plots won't overlap

    # COLORBAR 2
    m1 = vmin_
    m5 = vmax_
    m2 = float(1 * (m5 - m1) / 4 + m1)
    m3 = float(2 * (m5 - m1) / 4 + m1)
    m4 = float(3 * (m5 - m1) / 4 + m1)
    colorbar_labels = [m1, m2, m3, m4, m5]

    colorbar_labels_str = ["{:.2f}".format(round(x, 2)) for x in colorbar_labels]

    fig.subplots_adjust(bottom=0.12, hspace=0.32, wspace=0.30)
    cbaxes = fig.add_axes([0.3, 0.05, 0.438, 0.013])  # [left, bottom, width, height]
    cbar = fig.colorbar(im1, cax=cbaxes, orientation="horizontal")
    cbar.set_ticks(colorbar_labels)
    cbar.set_ticklabels(colorbar_labels_str)
    cbar.ax.tick_params(labelsize=tick_fontsize)

    if state_anomaly == "T":
        cbar.ax.set_xlabel("[K]", rotation=0, fontsize=title_fontsize)
    elif state_anomaly == "q":
        cbar.ax.set_xlabel("[g kg$\mathregular{^{-1}}$]", rotation=0, fontsize=title_fontsize)


    plt.show()

    # fig, axs = plt.subplots(nrows=3, ncols=4, figsize=(18, 12), facecolor='w', edgecolor='k')
    # fig.subplots_adjust(hspace=.5, wspace=.001)
    #
    # fig.delaxes(axs[2, 3])
    #
    # axs = axs.ravel()
    #
    # for i in range(len(pressures_all)):
    #
    #     matrix = matrix_all[i]
    #     matrix = np.array(matrix)
    #     pressures = pressures_all[i]
    #
    #     X1, Y1 = np.meshgrid(pressures,pressures)
    #
    #     im1 = axs[i].pcolor(X1,Y1,matrix,vmin=vmin_,vmax=vmax_,cmap="rainbow")
    #
    #     axs[i].set_ylim(axs[i].get_ylim()[::-1])  # invert y axis when y axis is pressure
    #     axs[i].set_xlim(axs[i].get_xlim()[::-1])
    #
    #     xy_label = "P [hPa]"
    #     axs[i].set_xlabel(xy_label, fontsize=label_fontsize)
    #     axs[i].set_ylabel(xy_label, fontsize=label_fontsize)
    #
    #     axs[i].set_xlim([1000, 200])
    #     axs[i].set_xticks([1000, 800, 600, 400, 200])
    #
    #     axs[i].set_ylim([1000, 200])
    #     axs[i].set_yticks([1000, 800, 600, 400, 200])
    #
    #     axs[i].tick_params(labelsize=tick_fontsize)
    #     axs[i].set_title(plot_titles[i], fontsize=title_fontsize)
    #
    #     plt.tight_layout()
    #
    # # COLORBAR
    # m1 = vmin_
    # m5 = vmax_
    # m2 = float(1 * (m5 - m1) / 4 + m1)
    # m3 = float(2 * (m5 - m1) / 4 + m1)
    # m4 = float(3 * (m5 - m1) / 4 + m1)
    # colorbar_labels = [m1, m2, m3, m4, m5]
    #
    # colorbar_labels_str = ["{:.2f}".format(round(x, 2)) for x in colorbar_labels]
    #
    # cbaxes = fig.add_axes([0.796, 0.052, 0.013, 0.258]) # [left, bottom, width, height]
    # cbar = fig.colorbar(im1, cax=cbaxes)
    # cbar.set_ticks(colorbar_labels)
    # cbar.set_ticklabels(colorbar_labels_str)
    # cbar.ax.tick_params(labelsize=tick_fontsize)
    #
    # if state_anomaly == "T":
    #     cbar.ax.set_ylabel("[K]", rotation=0, fontsize=title_fontsize, labelpad=20)
    # elif state_anomaly == "q":
    #     cbar.ax.set_ylabel("[g kg$\mathregular{^{-1}}$]", rotation=0, fontsize=title_fontsize, labelpad=40)

    #-----


if __name__ == "__main__":

    plot_matrix_all()

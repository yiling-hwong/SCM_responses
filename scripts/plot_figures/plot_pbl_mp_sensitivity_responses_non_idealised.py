#!/usr/bin/env python3

"""
Created on 13 August 2020

@author: Yi-Ling HWONG <yiling.hwong@gmail.com>

This script plots the comparison of the sensitivity of KF or BMJ responses (T' and q') to PBL and MP schemes
with and without idealised settings (ideal radiation and surface fluxes)

"""

"""
Select scheme to plot (KF = Kain-Fritsch; BMJ = Betts-Miller-Janjic)
"""

plot_kf_sensitivity = True
plot_bmj_sensitivity = True

#---------------------------

import numpy as np
import matplotlib.pyplot as plt

def get_data(scheme):

    """
    Get sensitivity data for scheme
    """

    configurations = ["yes_both","yes_both","no_rad","no_rad","no_sfc","no_sfc"] * 4
    sensitivities = (["pbl_sensitivity"] * 12) + (["mp_sensitivity"] * 12)
    vars = ((["t"] * 6) + (["q"] * 6)) * 2
    perturbations = ["dtdt","dqdt"] * 12

    pressures = []
    responses_all = []

    for i in range(24):

        responses = []

        lines = open("../../data/WRF/pbl_mp_sensitivity/responses_normalised_non_idealised/"+scheme+"_"
                     +configurations[i]+"_"+sensitivities[i]+"_"+vars[i]
                     +"_"+perturbations[i]+"_norm.csv","r").readlines()

        if i == 0:
            for line in lines[1:]:
                spline = line.rstrip("\n").split(",")
                pressures.append(float(spline[0]))

        for line in lines[1:]:
            spline = line.rstrip("\n").split(",")
            spline = [float(s) for s in spline]
            responses.append(spline[1:])

        responses_all.append(responses)


    """
    AVERAGE DTDT and DQDT
    """

    yes_both_combined_t_pbl = np.array(list(zip(responses_all[0],responses_all[1])))
    no_rad_combined_t_pbl = np.array(list(zip(responses_all[2],responses_all[3])))
    no_sfc_combined_t_pbl = np.array(list(zip(responses_all[4], responses_all[5])))
    yes_both_combined_q_pbl = np.array(list(zip(responses_all[6],responses_all[7])))
    no_rad_combined_q_pbl = np.array(list(zip(responses_all[8],responses_all[9])))
    no_sfc_combined_q_pbl = np.array(list(zip(responses_all[10], responses_all[11])))

    yes_both_combined_t_mp = np.array(list(zip(responses_all[12],responses_all[13])))
    no_rad_combined_t_mp = np.array(list(zip(responses_all[14],responses_all[15])))
    no_sfc_combined_t_mp = np.array(list(zip(responses_all[16], responses_all[17])))
    yes_both_combined_q_mp = np.array(list(zip(responses_all[18],responses_all[19])))
    no_rad_combined_q_mp = np.array(list(zip(responses_all[20],responses_all[21])))
    no_sfc_combined_q_mp = np.array(list(zip(responses_all[22], responses_all[23])))

    yes_both_t_pbl = np.mean(yes_both_combined_t_pbl, axis=1)
    no_rad_t_pbl = np.mean(no_rad_combined_t_pbl, axis=1)
    no_sfc_t_pbl = np.mean(no_sfc_combined_t_pbl, axis=1)
    yes_both_q_pbl = np.mean(yes_both_combined_q_pbl, axis=1)
    no_rad_q_pbl = np.mean(no_rad_combined_q_pbl, axis=1)
    no_sfc_q_pbl = np.mean(no_sfc_combined_q_pbl, axis=1)

    yes_both_t_mp = np.mean(yes_both_combined_t_mp, axis=1)
    no_rad_t_mp = np.mean(no_rad_combined_t_mp, axis=1)
    no_sfc_t_mp = np.mean(no_sfc_combined_t_mp, axis=1)
    yes_both_q_mp = np.mean(yes_both_combined_q_mp, axis=1)
    no_rad_q_mp = np.mean(no_rad_combined_q_mp, axis=1)
    no_sfc_q_mp = np.mean(no_sfc_combined_q_mp, axis=1)

    return pressures,yes_both_t_pbl,no_rad_t_pbl,no_sfc_t_pbl,yes_both_q_pbl,no_rad_q_pbl,no_sfc_q_pbl,yes_both_t_mp,no_rad_t_mp,no_sfc_t_mp,yes_both_q_mp,no_rad_q_mp,no_sfc_q_mp

def plot_sensitivity(scheme):

    if scheme == "kfeta":
        print ()
        print ("Plot sensitivity comparison for WRF Kain Fritsch...")
        pressures,yes_both_t_pbl,no_rad_t_pbl,no_sfc_t_pbl,yes_both_q_pbl,no_rad_q_pbl,no_sfc_q_pbl,yes_both_t_mp,no_rad_t_mp,no_sfc_t_mp,yes_both_q_mp,no_rad_q_mp,no_sfc_q_mp = get_data("kfeta")

    if scheme == "bmj":
        print ()
        print ("Plot sensitivity comparison for WRF Betts Miller Janjic...")
        pressures,yes_both_t_pbl,no_rad_t_pbl,no_sfc_t_pbl,yes_both_q_pbl,no_rad_q_pbl,no_sfc_q_pbl,yes_both_t_mp,no_rad_t_mp,no_sfc_t_mp,yes_both_q_mp,no_rad_q_mp,no_sfc_q_mp = get_data("bmj")


    x_zero = []

    for n in range(len(pressures)):
        x_zero.append(0.0)

    #############################
    # PLOT
    #############################

    # plot parameters
    legend_fontsize = 16
    tick_fontsize = 19
    label_fontsize = 19
    plot_title_fontsize = 22

    plot_titles = ["(a) Ideal rad & sfcflx on", "(b) Ideal rad. off", "(c) Ideal sfcflx off"]

    fig, ax = plt.subplots(nrows=4, ncols=3, figsize=(14, 22), facecolor='w', edgecolor='k')

    ax = ax.ravel()

    for i in range(12):

        if i == 0:
            ref_to_plot_dtdt = yes_both_t_pbl
        elif i == 1:
            ref_to_plot_dtdt = no_rad_t_pbl
        elif i == 2:
            ref_to_plot_dtdt = no_sfc_t_pbl
        elif i == 3:
            ref_to_plot_dtdt = yes_both_q_pbl
        elif i == 4:
            ref_to_plot_dtdt = no_rad_q_pbl
        elif i == 5:
            ref_to_plot_dtdt = no_sfc_q_pbl
        elif i == 6:
            ref_to_plot_dtdt = yes_both_t_mp
        elif i == 7:
            ref_to_plot_dtdt = no_rad_t_mp
        elif i == 8:
            ref_to_plot_dtdt = no_sfc_t_mp
        elif i == 9:
            ref_to_plot_dtdt = yes_both_q_mp
        elif i == 10:
            ref_to_plot_dtdt = no_rad_q_mp
        elif i == 11:
            ref_to_plot_dtdt = no_sfc_q_mp


        if i <= 5:
            ax[i].plot(ref_to_plot_dtdt[:, 0], pressures, "r", label="YSU")
            ax[i].plot(ref_to_plot_dtdt[:, 1], pressures, "b", label="MYNN2")
            ax[i].plot(ref_to_plot_dtdt[:, 2], pressures, "g", label="ACM2")
            ax[i].plot(ref_to_plot_dtdt[:, 3], pressures, "y", label="GBM")

        else:
            ax[i].plot(ref_to_plot_dtdt[:, 0], pressures, "m", label="WSM6")
            ax[i].plot(ref_to_plot_dtdt[:, 1], pressures, "c", label="Kessler")
            ax[i].plot(ref_to_plot_dtdt[:, 2], pressures, "olive", label="Thompson")
            ax[i].plot(ref_to_plot_dtdt[:, 3], pressures, "orange", label="Morrison")

        ax[i].plot(x_zero, pressures, "black", linestyle="dashed", linewidth="0.8")

        ax[i].set_ylim([200, 1000])
        ax[i].set_yticks([200, 300, 400, 500, 600, 700, 800, 900, 1000])

        ax[i].set_ylim(ax[i].get_ylim()[::-1])  # invert y axis when y axis is pressure

        ax[i].tick_params(labelsize=tick_fontsize)

        if i <= 2:
            ax[i].set_title(plot_titles[i], fontsize=plot_title_fontsize, pad = 20)

        if i == 0 or i == 1 or i == 2 or i == 6 or i == 7 or i == 8:

            ax[i].set_xlim([-0.35, 1.05])  # for T'
            ax[i].set_xticks([-0.3, 0.0, 0.3, 0.6, 0.9])

            x_label = "T' [K]"

        else:

            ax[i].set_xlim([-0.3, 0.5])  # for q'
            ax[i].set_xticks([-0.2, 0.0, 0.2, 0.4])

            x_label = "q' [g kg$\mathregular{^{-1}}$]"

        ax[i].set_xlabel(x_label, fontsize=label_fontsize)

        if i == 0 or i == 3 or i == 6 or i == 9:
            ax[i].set_ylabel("P [hPa]", fontsize=label_fontsize)
            ax[i].legend(loc="center right", fontsize=legend_fontsize)

    plt.tight_layout()

    plt.show()


def main():

    if plot_kf_sensitivity == True:
        plot_sensitivity("kfeta")

    if plot_bmj_sensitivity == True:
        plot_sensitivity("bmj")


if __name__ == "__main__":

    main()

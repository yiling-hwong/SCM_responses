#!/usr/bin/env python3

"""
Created on 13 August 2020

@author: Yi-Ling HWONG <yiling.hwong@gmail.com>

This script plots the sensitivity of responses (T' and q') to PBL and MP schemes

"""

"""
Select sensitivity to plot here (PBL or MP)
"""

plot_pbl_sensitivity = True
plot_mp_sensitivity = True

#---------------------------

import numpy as np
import matplotlib.pyplot as plt


def get_data(sensitivity):

    lines1 = open("../../data/WRF/pbl_mp_sensitivity/responses_normalised/kfeta_" + sensitivity + "_sensitivity_t_dtdt_norm.csv","r").readlines()
    lines2 = open("../../data/WRF/pbl_mp_sensitivity/responses_normalised/ntiedtke_" + sensitivity + "_sensitivity_t_dtdt_norm.csv","r").readlines()
    lines3 = open("../../data/WRF/pbl_mp_sensitivity/responses_normalised/nsas_" + sensitivity + "_sensitivity_t_dtdt_norm.csv","r").readlines()
    lines4 = open("../../data/WRF/pbl_mp_sensitivity/responses_normalised/bmj_" + sensitivity + "_sensitivity_t_dtdt_norm.csv","r").readlines()

    lines5 = open("../../data/WRF/pbl_mp_sensitivity/responses_normalised/kfeta_" + sensitivity + "_sensitivity_t_dqdt_norm.csv","r").readlines()
    lines6 = open("../../data/WRF/pbl_mp_sensitivity/responses_normalised/ntiedtke_" + sensitivity + "_sensitivity_t_dqdt_norm.csv","r").readlines()
    lines7 = open("../../data/WRF/pbl_mp_sensitivity/responses_normalised/nsas_" + sensitivity + "_sensitivity_t_dqdt_norm.csv","r").readlines()
    lines8 = open("../../data/WRF/pbl_mp_sensitivity/responses_normalised/bmj_" + sensitivity + "_sensitivity_t_dqdt_norm.csv","r").readlines()

    lines9 = open("../../data/WRF/pbl_mp_sensitivity/responses_normalised/kfeta_" + sensitivity + "_sensitivity_q_dtdt_norm.csv","r").readlines()
    lines10 = open("../../data/WRF/pbl_mp_sensitivity/responses_normalised/ntiedtke_" + sensitivity + "_sensitivity_q_dtdt_norm.csv","r").readlines()
    lines11 = open("../../data/WRF/pbl_mp_sensitivity/responses_normalised/nsas_" + sensitivity + "_sensitivity_q_dtdt_norm.csv","r").readlines()
    lines12 = open("../../data/WRF/pbl_mp_sensitivity/responses_normalised/bmj_" + sensitivity + "_sensitivity_q_dtdt_norm.csv","r").readlines()

    lines13 = open("../../data/WRF/pbl_mp_sensitivity/responses_normalised/kfeta_" + sensitivity + "_sensitivity_q_dqdt_norm.csv","r").readlines()
    lines14 = open("../../data/WRF/pbl_mp_sensitivity/responses_normalised/ntiedtke_" + sensitivity + "_sensitivity_q_dqdt_norm.csv","r").readlines()
    lines15 = open("../../data/WRF/pbl_mp_sensitivity/responses_normalised/nsas_" + sensitivity + "_sensitivity_q_dqdt_norm.csv","r").readlines()
    lines16 = open("../../data/WRF/pbl_mp_sensitivity/responses_normalised/bmj_" + sensitivity + "_sensitivity_q_dqdt_norm.csv","r").readlines()

    kfeta_t_dtdt = []
    ntiedtke_t_dtdt = []
    nsas_t_dtdt = []
    bmj_t_dtdt = []

    kfeta_t_dqdt = []
    ntiedtke_t_dqdt = []
    nsas_t_dqdt = []
    bmj_t_dqdt = []

    kfeta_q_dtdt = []
    ntiedtke_q_dtdt = []
    nsas_q_dtdt = []
    bmj_q_dtdt = []

    kfeta_q_dqdt = []
    ntiedtke_q_dqdt = []
    nsas_q_dqdt = []
    bmj_q_dqdt = []

    # ----
    # T_DTDT

    for line in lines1[1:]:
        spline = line.rstrip("\n").split(",")
        spline = [float(s) for s in spline]
        kfeta_t_dtdt.append(spline)

    for line in lines2[1:]:
        spline = line.rstrip("\n").split(",")
        spline = [float(s) for s in spline]
        ntiedtke_t_dtdt.append(spline)

    for line in lines3[1:]:
        spline = line.rstrip("\n").split(",")
        spline = [float(s) for s in spline]
        nsas_t_dtdt.append(spline)

    for line in lines4[1:]:
        spline = line.rstrip("\n").split(",")
        spline = [float(s) for s in spline]
        bmj_t_dtdt.append(spline)

    # ----
    # T_DQDT

    for line in lines5[1:]:
        spline = line.rstrip("\n").split(",")
        spline = [float(s) for s in spline]
        kfeta_t_dqdt.append(spline)

    for line in lines6[1:]:
        spline = line.rstrip("\n").split(",")
        spline = [float(s) for s in spline]
        ntiedtke_t_dqdt.append(spline)

    for line in lines7[1:]:
        spline = line.rstrip("\n").split(",")
        spline = [float(s) for s in spline]
        nsas_t_dqdt.append(spline)

    for line in lines8[1:]:
        spline = line.rstrip("\n").split(",")
        spline = [float(s) for s in spline]
        bmj_t_dqdt.append(spline)

    # ----
    # Q_DTDT

    for line in lines9[1:]:
        spline = line.rstrip("\n").split(",")
        spline = [float(s) for s in spline]
        kfeta_q_dtdt.append(spline)

    for line in lines10[1:]:
        spline = line.rstrip("\n").split(",")
        spline = [float(s) for s in spline]
        ntiedtke_q_dtdt.append(spline)

    for line in lines11[1:]:
        spline = line.rstrip("\n").split(",")
        spline = [float(s) for s in spline]
        nsas_q_dtdt.append(spline)

    for line in lines12[1:]:
        spline = line.rstrip("\n").split(",")
        spline = [float(s) for s in spline]
        bmj_q_dtdt.append(spline)

    # ----
    # Q_DQDT

    for line in lines13[1:]:
        spline = line.rstrip("\n").split(",")
        spline = [float(s) for s in spline]
        kfeta_q_dqdt.append(spline)

    for line in lines14[1:]:
        spline = line.rstrip("\n").split(",")
        spline = [float(s) for s in spline]
        ntiedtke_q_dqdt.append(spline)

    for line in lines15[1:]:
        spline = line.rstrip("\n").split(",")
        spline = [float(s) for s in spline]
        nsas_q_dqdt.append(spline)

    for line in lines16[1:]:
        spline = line.rstrip("\n").split(",")
        spline = [float(s) for s in spline]
        bmj_q_dqdt.append(spline)


    """
    AVERAGE DTDT and DQDT
    """

    kfeta_t_combined = np.array(list(zip(kfeta_t_dtdt,kfeta_t_dqdt)))
    ntiedtke_t_combined = np.array(list(zip(ntiedtke_t_dtdt, ntiedtke_t_dqdt)))
    nsas_t_combined = np.array(list(zip(nsas_t_dtdt, nsas_t_dqdt)))
    bmj_t_combined = np.array(list(zip(bmj_t_dtdt, bmj_t_dqdt)))

    kfeta_q_combined = np.array(list(zip(kfeta_q_dtdt,kfeta_q_dqdt)))
    ntiedtke_q_combined = np.array(list(zip(ntiedtke_q_dtdt, ntiedtke_q_dqdt)))
    nsas_q_combined = np.array(list(zip(nsas_q_dtdt, nsas_q_dqdt)))
    bmj_q_combined = np.array(list(zip(bmj_q_dtdt, bmj_q_dqdt)))

    kfeta_t = np.mean(kfeta_t_combined,axis=1)
    ntiedtke_t = np.mean(ntiedtke_t_combined, axis=1)
    nsas_t = np.mean(nsas_t_combined, axis=1)
    bmj_t = np.mean(bmj_t_combined, axis=1)

    kfeta_q = np.mean(kfeta_q_combined,axis=1)
    ntiedtke_q = np.mean(ntiedtke_q_combined, axis=1)
    nsas_q = np.mean(nsas_q_combined, axis=1)
    bmj_q = np.mean(bmj_q_combined, axis=1)

    return kfeta_t,ntiedtke_t,nsas_t,bmj_t,kfeta_q,ntiedtke_q,nsas_q,bmj_q


def plot_sensitivity(sensitivity):

    if sensitivity == "pbl":
        print()
        print ("Plotting PBL sensitivity ...")
        kfeta_t, ntiedtke_t, nsas_t, bmj_t, kfeta_q, ntiedtke_q, nsas_q, bmj_q = get_data("pbl")

    if sensitivity == "mp":
        print()
        print ("Plotting MP sensitivity ...")
        kfeta_t, ntiedtke_t, nsas_t, bmj_t, kfeta_q, ntiedtke_q, nsas_q, bmj_q = get_data("mp")


    x_zero = []
    for n in range(0, len(kfeta_t[:, 0])):
        x_zero.append(0.0)

    #############################
    # PLOT
    #############################

    # plot parameters
    mp_legend_fontsize = 16
    pbl_legend_fontsize = 19
    tick_fontsize = 19
    label_fontsize = 19
    plot_title_fontsize = 22
    figure_title_fontsize = 24

    plot_titles_top = ["(a) WRF-KF", "(b) WRF-NT", "(c) WRF-NSAS", "(d) WRF-BMJ"]
    plot_titles_bottom = ["(e) WRF-KF", "(f) WRF-NT", "(g) WRF-NSAS", "(h) WRF-BMJ"]

    fig, big_axes = plt.subplots(figsize=(18.0, 15.0), nrows=2, ncols=1, sharey=True)

    for row, big_ax in enumerate(big_axes, start=1):

        if row == 1:
            big_ax.set_title("Temperature responses \n", fontsize=figure_title_fontsize, pad=20)
        if row == 2:
            big_ax.set_title("Moisture responses \n", fontsize=figure_title_fontsize, pad=20)

        # Turn off axis lines and ticks of the big subplot
        # obs alpha is 0 in RGBA string!
        # big_ax.tick_params(labelcolor=(1., 1., 1., 0.0), top='off', bottom='off', left='off', right='off')
        big_ax.set_xticks([])
        big_ax.set_yticks([])
        # removes the white frame
        big_ax._frameon = False

    # T'
    for i in range(1, 5):

        if i == 1:
            ref_to_plot_t = kfeta_t
        elif i == 2:
            ref_to_plot_t = ntiedtke_t
        elif i == 3:
            ref_to_plot_t = nsas_t
        elif i == 4:
            ref_to_plot_t = bmj_t

        ax = fig.add_subplot(2, 4, i)

        ax.set_title(plot_titles_top[i - 1], fontsize=plot_title_fontsize)

        if sensitivity == "pbl":

            ax.plot(ref_to_plot_t[:, 1], ref_to_plot_t[:, 0], "r", label="YSU")
            ax.plot(ref_to_plot_t[:, 2], ref_to_plot_t[:, 0], "b", label="MYNN2")
            ax.plot(ref_to_plot_t[:, 3], ref_to_plot_t[:, 0], "g", label="ACM2")
            ax.plot(ref_to_plot_t[:, 4], ref_to_plot_t[:, 0], "y", label="GBM")


        if sensitivity == "mp":

            ax.plot(ref_to_plot_t[:, 1], ref_to_plot_t[:, 0], "m", label="WSM6")
            ax.plot(ref_to_plot_t[:, 2], ref_to_plot_t[:, 0], "c", label="Kessler")
            ax.plot(ref_to_plot_t[:, 3], ref_to_plot_t[:, 0], "olive", label="Thompson")
            ax.plot(ref_to_plot_t[:, 4], ref_to_plot_t[:, 0], "orange", label="Morrison")

        ax.plot(x_zero, ref_to_plot_t[:, 0], "black", linestyle="dashed", linewidth="0.8")

        ax = plt.gca()

        ax.set_ylim([200, 1000])
        ax.set_yticks([200, 300, 400, 500, 600, 700, 800, 900, 1000])

        ax.set_ylim(ax.get_ylim()[::-1])  # invert y axis when y axis is pressure

        ax.tick_params(labelsize=tick_fontsize)


        ax.set_xlim([-0.35, 1.05])  # for T'
        ax.set_xticks([-0.3, 0.0, 0.3, 0.6, 0.9])

        x_label = "T' [K]"

        ax.set_xlabel(x_label, fontsize=label_fontsize)

        if i != 1:
            empty_string_labels = [''] * 9
            ax.set_yticklabels(empty_string_labels)

        if i == 1:
            ax.set_ylabel("P [hPa]", fontsize=label_fontsize)

            if sensitivity == "pbl":
                ax.legend(loc="center right", fontsize=pbl_legend_fontsize)
            if sensitivity == "mp":
                ax.legend(loc="center right", fontsize=mp_legend_fontsize)

    # Q'
    for i in range(5, 9):

        if i == 5:
            ref_to_plot_q = kfeta_q
        elif i == 6:
            ref_to_plot_q = ntiedtke_q
        elif i == 7:
            ref_to_plot_q = nsas_q
        elif i == 8:
            ref_to_plot_q = bmj_q

        ax = fig.add_subplot(2, 4, i)

        ax.set_title(plot_titles_bottom[i - 5], fontsize=plot_title_fontsize)

        if sensitivity == "pbl":

            ax.plot(ref_to_plot_q[:, 1], ref_to_plot_q[:, 0], "r", label="YSU")
            ax.plot(ref_to_plot_q[:, 2], ref_to_plot_q[:, 0], "b", label="MYNN2")
            ax.plot(ref_to_plot_q[:, 3], ref_to_plot_q[:, 0], "g", label="ACM2")
            ax.plot(ref_to_plot_q[:, 4], ref_to_plot_q[:, 0], "y", label="GBM")

        if sensitivity == "mp":

            ax.plot(ref_to_plot_q[:, 1], ref_to_plot_q[:, 0], "m", label="WSM6")
            ax.plot(ref_to_plot_q[:, 2], ref_to_plot_q[:, 0], "c", label="Kessler")
            ax.plot(ref_to_plot_q[:, 3], ref_to_plot_q[:, 0], "olive", label="Thompson")
            ax.plot(ref_to_plot_q[:, 4], ref_to_plot_q[:, 0], "orange", label="Morrison")

        ax.plot(x_zero, ref_to_plot_q[:, 0], "black", linestyle="dashed", linewidth="0.8")

        ax = plt.gca()

        ax.set_ylim([200, 1000])
        ax.set_yticks([200, 300, 400, 500, 600, 700, 800, 900, 1000])

        ax.set_ylim(ax.get_ylim()[::-1])  # invert y axis when y axis is pressure

        ax.tick_params(labelsize=tick_fontsize)

        ax.set_xlim([-0.3, 0.5])  # for q'
        ax.set_xticks([-0.2, 0.0, 0.2, 0.4])

        x_label = "q' [g kg$\mathregular{^{-1}}$]"

        ax.set_xlabel(x_label, fontsize=label_fontsize)

        if i != 5:
            empty_string_labels = [''] * 9
            ax.set_yticklabels(empty_string_labels)

        if i == 5:
            ax.set_ylabel("P [hPa]", fontsize=label_fontsize)

            if sensitivity == "pbl":
                ax.legend(loc="center right", fontsize=pbl_legend_fontsize)
            if sensitivity == "mp":
                ax.legend(loc="center right", fontsize=mp_legend_fontsize)


    fig.set_facecolor('w')
    plt.tight_layout()

    plt.show()

def main():

    if plot_pbl_sensitivity == True:
        plot_sensitivity("pbl")

    if plot_mp_sensitivity == True:
        plot_sensitivity("mp")


if __name__ == "__main__":

    main()

#!/usr/bin/env python3

"""
Created on 13 August 2020

@author: Yi-Ling HWONG <yiling.hwong@gmail.com>

This script plots the RCE mean states (T and RH) sensitivity to PBL and MP schemes

"""

"""
Select mean state to plot here (T or RH)
"""

plot_T_sensitivity = True # plot theta_es sensitivity
plot_rh_sensitivity = True

# plot parameters
cape_fontsize = 17
tick_fontsize = 19
label_fontsize = 19
legend_fontsize = 19
plot_title_fontsize = 22
figure_title_fontsize = 24

#-------------------------

import numpy as np
import matplotlib.pyplot as plt

def get_theta_e_sfc_and_cape():

    """
    Get theta_e_sfc and cape for plotting theta_es sensitivity
    """

    lines1 = open("../../data/WRF/pbl_mp_sensitivity/mean_states/kfeta_pbl_sensitivity_theta_es.csv","r").readlines()
    lines2 = open("../../data/WRF/pbl_mp_sensitivity/mean_states/ntiedtke_pbl_sensitivity_theta_es.csv","r").readlines()
    lines3 = open("../../data/WRF/pbl_mp_sensitivity/mean_states/nsas_pbl_sensitivity_theta_es.csv","r").readlines()
    lines4 = open("../../data/WRF/pbl_mp_sensitivity/mean_states/bmj_pbl_sensitivity_theta_es.csv","r").readlines()

    lines5 = open("../../data/WRF/pbl_mp_sensitivity/mean_states/kfeta_mp_sensitivity_theta_es.csv","r").readlines()
    lines6 = open("../../data/WRF/pbl_mp_sensitivity/mean_states/ntiedtke_mp_sensitivity_theta_es.csv","r").readlines()
    lines7 = open("../../data/WRF/pbl_mp_sensitivity/mean_states/nsas_mp_sensitivity_theta_es.csv","r").readlines()
    lines8 = open("../../data/WRF/pbl_mp_sensitivity/mean_states/bmj_mp_sensitivity_theta_es.csv","r").readlines()

    theta_es_pbl = []
    cape_pbl = []
    theta_es_mp = []
    cape_mp = []

    # PBL
    for line in lines1[:1]:
        spline = line.rstrip("\n").split(",")
        theta_es_pbl.append(float(spline[0]))
        cape_pbl.append(int(spline[1]))

    for line in lines2[:1]:
        spline = line.rstrip("\n").split(",")
        theta_es_pbl.append(float(spline[0]))
        cape_pbl.append(int(spline[1]))

    for line in lines3[:1]:
        spline = line.rstrip("\n").split(",")
        theta_es_pbl.append(float(spline[0]))
        cape_pbl.append(int(spline[1]))

    for line in lines4[:1]:
        spline = line.rstrip("\n").split(",")
        theta_es_pbl.append(float(spline[0]))
        cape_pbl.append(int(spline[1]))

    # MP
    for line in lines5[:1]:
        spline = line.rstrip("\n").split(",")
        theta_es_mp.append(float(spline[0]))
        cape_mp.append(int(spline[1]))

    for line in lines6[:1]:
        spline = line.rstrip("\n").split(",")
        theta_es_mp.append(float(spline[0]))
        cape_mp.append(int(spline[1]))

    for line in lines7[:1]:
        spline = line.rstrip("\n").split(",")
        theta_es_mp.append(float(spline[0]))
        cape_mp.append(int(spline[1]))

    for line in lines8[:1]:
        spline = line.rstrip("\n").split(",")
        theta_es_mp.append(float(spline[0]))
        cape_mp.append(int(spline[1]))

    return theta_es_pbl,cape_pbl,theta_es_mp,cape_mp


def get_sensitivity_data(var):

    """
    Get PBL and MP sensitivity data
    """

    lines1 = open("../../data/WRF/pbl_mp_sensitivity/mean_states/kfeta_pbl_sensitivity_"+var+".csv","r").readlines()
    lines2 = open("../../data/WRF/pbl_mp_sensitivity/mean_states/ntiedtke_pbl_sensitivity_"+var+".csv","r").readlines()
    lines3 = open("../../data/WRF/pbl_mp_sensitivity/mean_states/nsas_pbl_sensitivity_"+var+".csv","r").readlines()
    lines4 = open("../../data/WRF/pbl_mp_sensitivity/mean_states/bmj_pbl_sensitivity_"+var+".csv","r").readlines()

    lines5 = open("../../data/WRF/pbl_mp_sensitivity/mean_states/kfeta_mp_sensitivity_"+var+".csv","r").readlines()
    lines6 = open("../../data/WRF/pbl_mp_sensitivity/mean_states/ntiedtke_mp_sensitivity_"+var+".csv","r").readlines()
    lines7 = open("../../data/WRF/pbl_mp_sensitivity/mean_states/nsas_mp_sensitivity_"+var+".csv","r").readlines()
    lines8 = open("../../data/WRF/pbl_mp_sensitivity/mean_states/bmj_mp_sensitivity_"+var+".csv","r").readlines()


    kfeta_pbl = []
    ntiedtke_pbl = []
    nsas_pbl = []
    bmj_pbl = []

    kfeta_mp = []
    ntiedtke_mp = []
    nsas_mp = []
    bmj_mp = []

    # PBL
    for line in lines1[2:]:
        spline = line.rstrip("\n").split(",")
        spline = [float(s) for s in spline]
        kfeta_pbl.append(spline)

    for line in lines2[2:]:
        spline = line.rstrip("\n").split(",")
        spline = [float(s) for s in spline]
        ntiedtke_pbl.append(spline)

    for line in lines3[2:]:
        spline = line.rstrip("\n").split(",")
        spline = [float(s) for s in spline]
        nsas_pbl.append(spline)

    for line in lines4[2:]:
        spline = line.rstrip("\n").split(",")
        spline = [float(s) for s in spline]
        bmj_pbl.append(spline)


    # MP
    for line in lines5[2:]:
        spline = line.rstrip("\n").split(",")
        spline = [float(s) for s in spline]
        kfeta_mp.append(spline)

    for line in lines6[2:]:
        spline = line.rstrip("\n").split(",")
        spline = [float(s) for s in spline]
        ntiedtke_mp.append(spline)

    for line in lines7[2:]:
        spline = line.rstrip("\n").split(",")
        spline = [float(s) for s in spline]
        nsas_mp.append(spline)

    for line in lines8[2:]:
        spline = line.rstrip("\n").split(",")
        spline = [float(s) for s in spline]
        bmj_mp.append(spline)

    kfeta_pbl = np.array(kfeta_pbl)
    ntiedtke_pbl = np.array(ntiedtke_pbl)
    nsas_pbl = np.array(nsas_pbl)
    bmj_pbl = np.array(bmj_pbl)

    kfeta_mp = np.array(kfeta_mp)
    ntiedtke_mp = np.array(ntiedtke_mp)
    nsas_mp = np.array(nsas_mp)
    bmj_mp = np.array(bmj_mp)

    return kfeta_pbl,ntiedtke_pbl,nsas_pbl,bmj_pbl,kfeta_mp,ntiedtke_mp,nsas_mp,bmj_mp


def plot_sensitivity(var):

    if var == "T":
        print ()
        print ("Plot T sensitivity ...")
        kfeta_pbl, ntiedtke_pbl, nsas_pbl, bmj_pbl, kfeta_mp, ntiedtke_mp, nsas_mp, bmj_mp = get_sensitivity_data("theta_es")

    if var == "rh":
        print ()
        print ("Plot rh sensitivity ...")
        kfeta_pbl, ntiedtke_pbl, nsas_pbl, bmj_pbl, kfeta_mp, ntiedtke_mp, nsas_mp, bmj_mp = get_sensitivity_data("rh")


    #############################
    # PLOT
    #############################

    plot_titles_top = ["(a) WRF-KF", "(b) WRF-NT", "(c) WRF-NSAS", "(d) WRF-BMJ"]
    plot_titles_bottom = ["(e) WRF-KF", "(f) WRF-NT", "(g) WRF-NSAS", "(h) WRF-BMJ"]

    fig, big_axes = plt.subplots(figsize=(18.0, 16.0), nrows=2, ncols=1, sharey=True)

    for row, big_ax in enumerate(big_axes, start=1):

        if row == 1:
            big_ax.set_title("PBL sensitivity \n", fontsize=figure_title_fontsize, pad=20)
        if row == 2:
            big_ax.set_title("MP sensitivity \n", fontsize=figure_title_fontsize, pad=20)

        # Turn off axis lines and ticks of the big subplot
        # obs alpha is 0 in RGBA string!
        big_ax.set_xticks([])
        big_ax.set_yticks([])
        # removes the white frame
        big_ax._frameon = False


    for i in range(1, 5):

        if i == 1:
            ref_to_plot = kfeta_pbl
        elif i == 2:
            ref_to_plot = ntiedtke_pbl
        elif i == 3:
            ref_to_plot = nsas_pbl
        elif i == 4:
            ref_to_plot = bmj_pbl

        ax = fig.add_subplot(2, 4, i)
        ax.set_title(plot_titles_top[i-1], fontsize=plot_title_fontsize)
        ax.plot(ref_to_plot[:,4], ref_to_plot[:,0], "r", label="YSU")
        ax.plot(ref_to_plot[:,5], ref_to_plot[:,1], "b", label="MYNN2")
        ax.plot(ref_to_plot[:,6], ref_to_plot[:,2], "g", label="ACM2")
        ax.plot(ref_to_plot[:,7], ref_to_plot[:,3], "y", label="GBM")

        ax = plt.gca()

        ax.set_ylim([50, 1000])
        ax.set_yticks([100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])

        ax.set_ylim(ax.get_ylim()[::-1])  # invert y axis when y axis is pressure

        ax.tick_params(labelsize=tick_fontsize)

        if var == "T":

            theta_es_pbl, cape_pbl, theta_es_mp, cape_mp = get_theta_e_sfc_and_cape()

            ax.set_xlim([314, 366])
            ax.set_xticks([320,340,360])
            x_label = "$\it{θ_{es}}$ [K]"
            ax.set_xlabel(x_label, fontsize=label_fontsize)

            theta_x_position = theta_es_pbl[i-1] - 2

            ax.annotate("$\overline{θ_{e}^{b}}$", (theta_x_position, 975), fontsize=label_fontsize, color="red")

            ax.annotate("$\overline{CAPE}$ = "+str(cape_pbl[i-1])+" J kg$\mathregular{^{-1}}$", (333.5, 750), fontsize=cape_fontsize)

            # Plot 2nd overlapping x axis for theta_e_sfc

            ax2 = ax.twiny()

            # Move twinned axis ticks and label from top to bottom
            ax2.xaxis.set_ticks_position("bottom")
            ax2.xaxis.set_label_position("bottom")

            # Offset the twin axis below the host
            ax2.spines["bottom"].set_position(("axes", 0))

            # Turn on the frame for the twin axis, but then hide all
            # but the bottom spine
            ax2.set_frame_on(True)
            ax2.patch.set_visible(False)

            ax2.spines["bottom"].set_visible(True)

            ax2.set_xlim([314, 366])  # for T'
            new_tick_locations = np.array([theta_es_pbl[i-1]])

            ax2.set_xticks(new_tick_locations)

            invisible_xticks = [" "]
            ax2.set_xticklabels(invisible_xticks,fontsize=label_fontsize)

            ax2.tick_params(axis='x', colors='red',size=6)

        if var == "rh":

            ax.set_xlim([0, 100])  # for T'
            ax.set_xticks([0, 25, 50, 75, 100])
            x_label = "RH [%]"

            ax.set_xlabel(x_label, fontsize=label_fontsize)


        if i != 1:
            empty_string_labels = [''] * 4
            ax.set_yticklabels(empty_string_labels)

        if i == 1:
            ax.set_ylabel("P [hPa]", fontsize=label_fontsize)
            ax.legend(loc="center right", fontsize=legend_fontsize)


    for i in range(5, 9):

        if i == 5:
            ref_to_plot = kfeta_mp
        elif i == 6:
            ref_to_plot = ntiedtke_mp
        elif i == 7:
            ref_to_plot = nsas_mp
        elif i == 8:
            ref_to_plot = bmj_mp

        ax = fig.add_subplot(2, 4, i)
        ax.set_title(plot_titles_bottom[i-5], fontsize=plot_title_fontsize)
        ax.plot(ref_to_plot[:, 4], ref_to_plot[:, 0], "m", label="WSM6")
        ax.plot(ref_to_plot[:, 5], ref_to_plot[:, 1], "c", label="Kessler")
        ax.plot(ref_to_plot[:, 6], ref_to_plot[:, 2], "olive", label="Thompson")
        ax.plot(ref_to_plot[:, 7], ref_to_plot[:, 3], "orange", label="Morrison")

        ax = plt.gca()

        ax.set_ylim([50, 1000])
        ax.set_yticks([100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])

        ax.set_ylim(ax.get_ylim()[::-1])  # invert y axis when y axis is pressure

        ax.tick_params(labelsize=tick_fontsize)

        if var == "T":

            theta_es_pbl, cape_pbl, theta_es_mp, cape_mp = get_theta_e_sfc_and_cape()

            ax.set_xlim([314, 366])  # for T'
            ax.set_xticks([320, 340, 360])
            x_label = "$\it{θ_{es}}$ [K]"
            ax.set_xlabel(x_label, fontsize=label_fontsize)

            theta_x_position = theta_es_mp[i - 5] - 2

            ax.annotate("$\overline{θ_{e}^{b}}$", (theta_x_position, 975), fontsize=label_fontsize, color="red")

            ax.annotate("$\overline{CAPE}$ = " + str(cape_mp[i - 5]) + " J kg$\mathregular{^{-1}}$", (333.5, 750),
                        fontsize=cape_fontsize)

            # Plot 2nd overlapping x axis for theta_e_sfc

            ax2 = ax.twiny()

            # Move twinned axis ticks and label from top to bottom
            ax2.xaxis.set_ticks_position("bottom")
            ax2.xaxis.set_label_position("bottom")

            # Offset the twin axis below the host
            ax2.spines["bottom"].set_position(("axes", 0))

            # Turn on the frame for the twin axis, but then hide all
            # but the bottom spine
            ax2.set_frame_on(True)
            ax2.patch.set_visible(False)

            ax2.spines["bottom"].set_visible(True)

            ax2.set_xlim([314, 366])  # for T'
            new_tick_locations = np.array([theta_es_mp[i - 5]])

            ax2.set_xticks(new_tick_locations)

            invisible_xticks = [" "]
            ax2.set_xticklabels(invisible_xticks, fontsize=label_fontsize)

            ax2.tick_params(axis='x', colors='red', size=6)

        if var == "rh":

            ax.set_xlim([0, 100])  # for T'
            ax.set_xticks([0, 25, 50, 75, 100])
            x_label = "RH [%]"

            ax.set_xlabel(x_label, fontsize=label_fontsize)

        if i != 5:
            empty_string_labels = [''] * 4
            ax.set_yticklabels(empty_string_labels)

        if i == 5:
            ax.set_ylabel("P [hPa]", fontsize=label_fontsize)
            ax.legend(loc="center right", fontsize=legend_fontsize)

    fig.set_facecolor('w')
    plt.tight_layout()

    plt.show()


def main():

    if plot_T_sensitivity == True:
        plot_sensitivity("T")

    if plot_rh_sensitivity == True:
        plot_sensitivity("rh")


if __name__ == "__main__":

    main()

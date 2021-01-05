#!/usr/bin/env python3

"""
Created on 13 August 2020
@author: Yi-Ling HWONG
"""

import numpy as np
import matplotlib.pyplot as plt
import math


def get_perturbations_T(pressures):

    num_levels = len(pressures)

    TtendAmp = 0.5

    tend_list = []

    for n in range(1,num_levels+1):

        tten = []

        for i in range(1,num_levels+1):

            if i == n:
                deltaf = 1.0
            else:
                deltaf = 0.0

            tten_temp = 0.5 * (TtendAmp / 86400 * (deltaf + math.exp(-((pressures[i - 1] - pressures[n - 1]) / 7500.) ** 2.)))

            tten.append(tten_temp*86400) # in units of K/d

        tend_list.append(tten)

    Y_tendency_matrix = np.array(tend_list)
    Y_tendency_matrix = np.array(Y_tendency_matrix).transpose()

    return (Y_tendency_matrix)


def get_perturbations_q(pressures):

    num_levels = len(pressures)

    QtendAmp = 0.2

    tend_list = []

    for n in range(1, num_levels + 1):

        qten = []

        for i in range(1, num_levels + 1):

            if i == n:
                deltaf = 1.0
            else:
                deltaf = 0.0

            qten_temp = 0.5 * (
                        QtendAmp / 86400 * (deltaf + math.exp(-((pressures[i - 1] - pressures[n - 1]) / 7500.) ** 2.)))

            qten.append(qten_temp * 86400)  # in units of g/kg/d

        tend_list.append(qten)

    Y_tendency_matrix = np.array(tend_list)
    Y_tendency_matrix = np.array(Y_tendency_matrix).transpose()

    return (Y_tendency_matrix)


def get_power_input_T(pressures,tten,psfc):

    c_p = 1004.6 # units: J/K/kg
    g = 9.81

    delta_p = []

    delta_p.append(psfc*100-pressures[0])

    for index,value in enumerate(pressures[:-1]):
        delt_p = pressures[index] - pressures[index+1]
        delta_p.append(delt_p)

    print ()
    print ("------------------------------")
    print ("GET POWER INPUT T")
    print ("Delta pressure length:",len(delta_p))
    print (delta_p)

    """
    Calculate total heat
    """

    total_heat = []

    for n in range(len(pressures)):

        heat_profile = []
        tten_profile = tten[:,n].tolist()

        for i in range(len(tten_profile)):
            dTdt = tten_profile[i] / 86400
            del_p = delta_p[i]
            power = (c_p / g) * del_p * dTdt
            heat_profile.append(power)

        total_heat_per_column = np.sum(heat_profile)
        total_heat.append(total_heat_per_column)

    print ()
    print ("TOTAL HEAT length (in W m^-2):", len(total_heat))
    print (total_heat)

    # for th in total_heat:
    #     print (th)

    return (total_heat)

def get_power_input_q(pressures,qten,psfc):

    L_v = 2.5e6  # units: J/kg
    g = 9.81

    delta_p = []

    delta_p.append(psfc * 100 - pressures[0])

    for index, value in enumerate(pressures[:-1]):
        delt_p = pressures[index] - pressures[index + 1]
        delta_p.append(delt_p)

    print()
    print("------------------------------")
    print("GET POWER INPUT Q")
    print("Delta pressure length:", len(delta_p))
    print(delta_p)

    """
    Calculate total heat
    """

    total_heat = []

    for n in range(len(pressures)):

        heat_profile = []
        qten_profile = qten[:, n].tolist()

        for i in range(len(qten_profile)):
            dqdt = (qten_profile[i] / 86400) / 1000
            del_p = delta_p[i]
            power = (L_v / g) * del_p * dqdt
            heat_profile.append(power)

        total_heat_per_column = np.sum(heat_profile)
        total_heat.append(total_heat_per_column)

    print()
    print("TOTAL HEAT length (in W m^-2):", len(total_heat))
    print(total_heat)

    # for th in total_heat:
    #     print (th)

    return (total_heat)

def get_power_input_kuang(perturb_t,perturb_q):

    print ()
    print ("#########################")
    print ("GET KUANG POWER INPUT")
    print ()

    start_level_t = 1
    end_level_t = 20
    start_level_q = 1
    end_level_q = 18
    num_levels = 18

    lines1 = open("../../data/SAM/matrix_X_raw/p.csv","r").readlines()
    lines3 = open("../../data/SAM/matrix_X_raw/dX.csv","r").readlines()

    pressures = []
    y_tendencies = []

    for line in lines1:
        spline = line.rstrip("\n")
        pressures.append(float(spline) * 100)

    if pressures[0] < pressures[1]:
        pressures.reverse()

    pressures = pressures[:num_levels]

    """
    Get tendency matrix
    """

    for line in lines3:
        spline = line.rstrip("\n").split(",")
        spline = [-float(s) for s in spline]
        y_tendencies.append(spline)

    Y_tendency_matrix_complete = np.array(y_tendencies)

    if perturb_t == True:
        Y_tendency_matrix_extracted = Y_tendency_matrix_complete[:(end_level_q-start_level_q+1),:(end_level_q-start_level_q+1)]
    elif perturb_q == True:
        Y_tendency_matrix_extracted = Y_tendency_matrix_complete[(end_level_t-start_level_t+1):,(end_level_t-start_level_t+1):]

    """
    Get power input
    """

    psfc = 1014.8

    if perturb_t == 1:
        power_input_list = get_power_input_T(pressures,Y_tendency_matrix_extracted,psfc)
    elif perturb_q == 1:
        power_input_list = get_power_input_q(pressures,Y_tendency_matrix_extracted,psfc)

    print ()
    print ("KUANG power input:",len(power_input_list))
    print (power_input_list)
    print ()
    print ("#########################")

    return (power_input_list)


def get_m_inverse(X,Y):

    print ()
    print ("---------------------------")
    print ("GET M INVERSE")
    print ()

    Y_diag = np.diag(Y)
    Y_inv = np.linalg.inv(Y_diag)

    print ("Y_diag and Y_inv shapes:", Y_diag.shape,Y_inv.shape)
    # print (Y_diag)
    # print ()
    # print (Y_inv)

    M_inv = np.matmul(X,Y_inv)

    return (M_inv)

def plot_matrix(matrix,pressures_x,pressures_y,plot_title,vmin_,vmax_):

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)

    to_plot = matrix

    X1, Y1 = np.meshgrid(pressures_x, pressures_y)

    im1 = ax1.pcolor(X1, Y1, to_plot, vmin=vmin_, vmax=vmax_, cmap="rainbow")

    ax1.set_ylim(ax1.get_ylim()[::-1])  # invert y axis when y axis is pressure
    ax1.set_xlim(ax1.get_xlim()[::-1])

    plt.xlabel("P (hPa)", fontsize=12)
    plt.ylabel("P (hPa)", fontsize=12)

    #------
    # colorbar

    m1 = vmin_
    m5 = vmax_
    m2 = float(1 * (m5 - m1) / 4 + m1)
    m3 = float(2 * (m5 - m1) / 4 + m1)
    m4 = float(3 * (m5 - m1) / 4 + m1)
    colorbar_labels = [m1, m2, m3, m4, m5]

    colorbar_labels_str = ["{:.2f}".format(round(x, 2)) for x in colorbar_labels]

    cbar = plt.colorbar(im1)
    cbar.set_ticks(colorbar_labels)
    cbar.set_ticklabels(colorbar_labels_str)
    cbar.ax.tick_params(labelsize=12)

    plt.title(plot_title, fontsize=12)
    plt.gca().xaxis.tick_bottom()

    ax1.tick_params(labelsize=12)

    ax1.set_xlim([1000, 200])
    plt.xticks([1000, 800, 600, 400, 200])

    ax1.set_ylim([1000, 200])
    plt.yticks([1000, 800, 600, 400, 200])

    plt.show()

def plot_profile(standardise_kuang,pressures,x_profile_1,x_profile_2,m_profile_1,m_profile_2,plot_title,label_level_1,label_level_2,state_anomaly):

    #--------------------
    # For plotting zero line

    x_zero = []

    for n in range(0,len(pressures)):
        x_zero.append(0.0)

    y_label = "P (hPa)"

    #---------------
    # plot

    if standardise_kuang == True:

        fig1 = plt.figure(figsize=(4, 5))
        ax1 = fig1.add_subplot(111)

        ax1.plot(m_profile_1, pressures, "r", label=label_level_1)
        ax1.plot(m_profile_2, pressures, "b", label=label_level_2)
        ax1.plot(x_zero, pressures, "black", linestyle="dashed", linewidth="0.8")

        ax1 = plt.gca()

        ax1.set_ylim([200, 1000])
        ax1.set_ylim(ax1.get_ylim()[::-1])  # invert y axis when y axis is pressure

        if state_anomaly == "T":

            ax1.set_xlim([-0.35, 1.05])  # for T'
            plt.xticks([-0.3, 0.0, 0.3, 0.6, 0.9])
            x_label = "T' (K)"

        elif state_anomaly == "q":

            ax1.set_xlim([-0.3, 0.5])  # for q'
            plt.xticks([-0.2, 0.0, 0.2, 0.4])
            x_label = "q' (g/kg)"


        ax1.tick_params(labelsize=12)

        plt.xlabel(x_label,fontsize=12)
        plt.ylabel(y_label,fontsize=12)

        plt.title(plot_title, fontsize=12)

        plt.tight_layout()  # so that the two plots won't overlap

        plt.show()

    elif standardise_kuang == False:

        fig1 = plt.figure(figsize=(4, 5))
        ax1 = fig1.add_subplot(111)

        ax1.plot(m_profile_1, pressures, "r", label=label_level_1)
        ax1.plot(m_profile_2, pressures, "b", label=label_level_2)
        ax1.plot(x_zero, pressures, "black", linestyle="dashed", linewidth="0.8")

        ax1 = plt.gca()

        ax1.set_ylim([200, 1000])  # for q'
        ax1.set_ylim(ax1.get_ylim()[::-1])  # invert y axis when y axis is pressure

        if state_anomaly == "T":

            ax1.set_xlim([-0.055, 0.18])  # for T'
            plt.xticks([-0.05, 0.0, 0.05, 0.10, 0.15])
            x_label = "T' (K/(W m$\mathregular{^{-2}}$)"

        elif state_anomaly == "q":

            ax1.set_xlim([-0.035, 0.095])  # for q'
            plt.xticks([-0.03, 0.00, 0.03, 0.06, 0.09])
            x_label = "q' (g/kg/(W m$\mathregular{^{-2}}$)"

        ax1.tick_params(labelsize=12)

        leg1 = ax1.legend(loc="upper left",fontsize=12)

        plt.xlabel(x_label,fontsize=12)
        plt.ylabel(y_label,fontsize=12)

        plt.title(plot_title+" M$\mathregular{^{-1}}$ (norm.power)", fontsize=9)

        plt.tight_layout()  # so that the two plots won't overlap

        plt.show()

        #-------
        # Plot X raw

        fig1 = plt.figure(figsize=(4, 5))
        ax1 = fig1.add_subplot(111)

        ax1.plot(x_profile_1, pressures, "r", label=label_level_1)
        ax1.plot(x_profile_2, pressures, "b", label=label_level_2)
        ax1.plot(x_zero, pressures, "black", linestyle="dashed", linewidth="0.8")

        ax1 = plt.gca()

        ax1.set_ylim([200, 1000])  # for q'
        ax1.set_ylim(ax1.get_ylim()[::-1])  # invert y axis when y axis is pressure

        if state_anomaly == "T":

            ax1.set_xlim([-0.35, 1.05])  # for T'
            plt.xticks([-0.3, 0.0, 0.3, 0.6, 0.9])
            x_label = "T' (K)"

        elif state_anomaly == "q":

            ax1.set_xlim([-0.3, 0.5]) # for q'
            plt.xticks([-0.2,0.0,0.2,0.4])
            x_label = "q' (g/kg)"

        ax1.tick_params(labelsize=8)

        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(plot_title+" X (raw)", fontsize=9)

        plt.tight_layout()  # so that the two plots won't overlap

        plt.show()

def write_profiles_to_file(standardise_kuang,pressures, x_profile_1, x_profile_2, m_profile_1, m_profile_2,model_folder,file_name,
                  label_level_1,label_level_2):

    print ()
    print ("WRITING to file ...")

    m_zipped = zip(pressures, m_profile_1, m_profile_2)
    x_zipped = zip(pressures, x_profile_1, x_profile_2)

    m_profiles_str = []
    x_profiles_str = []

    for mz in m_zipped:
        mz_list = list(mz)
        mz_list_str = [str(m) for m in mz_list]
        m_profiles_str.append(mz_list_str)

    for xz in x_zipped:
        xz_list = list(xz)
        xz_list_str = [str(x) for x in xz_list]
        x_profiles_str.append(xz_list_str)

    m_profiles_str.insert(0, ["pressures(hPa)," + label_level_1 + "," + label_level_2])
    x_profiles_str.insert(0, ["pressures(hPa)," + label_level_1 + "," + label_level_2])

    if standardise_kuang == True:

        path1 = "../../data/"+model_path+"/response_profiles/"+file_name+"_norm_kuang.csv"
        f_norm_kuang = open(path1,"w")
        print ()
        print (path1)

        for ms in m_profiles_str:
            f_norm_kuang.write(",".join(ms) + "\n")

        f_norm_kuang.close()

    elif standardise_kuang == False:

        path2 = "../../data/"+model_folder+"/response_profiles/"+file_name+"_norm_power.csv"
        path3 = "../../data/"+model_folder+"/response_profiles/"+file_name+"_raw.csv"

        print ()
        print (path2)
        print (path3)

        f_norm_power = open(path2,"w")

        for ms in m_profiles_str:
            f_norm_power.write(",".join(ms) + "\n")

        f_norm_power.close()

        f_raw = open(path3,"w")

        for ms in x_profiles_str:
            f_raw.write(",".join(ms) + "\n")

        f_raw.close()

def write_matrix_to_file(M_inv,model_folder,file_name):

    print ()
    print ("##########################")
    print ("Writing M_INV to file...")
    print ()

    M_inv_list = M_inv.tolist()

    M_inv_str = []

    for mi in M_inv_list:
        temp = [str(x) for x in mi]
        M_inv_str.append(temp)

    print ("Length M_INV to write:",len(M_inv_str),len(M_inv_str[0]))

    f = open("/../../data/"+model_folder+"/matrix_M_inv/"+"M_inv_"+file_name+".csv","w")

    for mi in M_inv_str:
        f.write(",".join(mi)+"\n")

    f.close()




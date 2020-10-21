#!/usr/bin/env python3

"""
Created on 13 August 2020
@author: Yi-Ling HWONG
"""

import sys
import numpy as np
from utils import *

class SAM():

    def __init__(self,standardise_kuang,perturb_t,perturb_q,var_x,t_amplitude,q_amplitude):

        self.standardise_kuang = standardise_kuang
        self.perturb_t = perturb_t
        self.perturb_q = perturb_q
        self.var_x = var_x
        self.t_amplitude = t_amplitude
        self.q_amplitude = q_amplitude

        self.start_level_t = 1
        self.end_level_t = 20
        self.start_level_q = 1
        self.end_level_q = 18
        self.num_levels = 18

        print ()
        print ("SAM")
        print ()

        lines1 = open("/Users/yi-linghwong/Documents/postdoc/publications/paracon_scm_paper/_DATA/_SAM/matrix_X/p.csv","r").readlines()

        pressures = []

        for line in lines1:
            spline = line.rstrip("\n")
            pressures.append(float(spline) * 100)

        if pressures[0] < pressures[1]:
            pressures.reverse()

        self.pressures = pressures[:self.num_levels]
        pressures_hpa = [ p/100 for p in pressures]
        self.pressures_t = pressures_hpa[:self.end_level_t]
        self.pressures_q = pressures_hpa[:self.end_level_q]

        self.pressures_x = self.pressures_q
        self.pressures_y = self.pressures_q

        # if perturb_t == 1 and var_x == "T_PHY":
        #     self.pressures_x = self.pressures_t
        #     self.pressures_y = self.pressures_t
        # elif perturb_t == 1 and var_x == "QVAPOR":
        #     self.pressures_x = self.pressures_t
        #     self.pressures_y = self.pressures_q
        # elif perturb_q == 1 and var_x == "T_PHY":
        #     self.pressures_x = self.pressures_q
        #     self.pressures_y = self.pressures_t
        # elif perturb_q == 1 and var_x == "QVAPOR":
        #     self.pressures_x = self.pressures_q
        #     self.pressures_y = self.pressures_q


    def get_X_anomaly(self):

        """
        Get anomaly matrix
        """

        perturb_t = self.perturb_t
        perturb_q = self.perturb_q
        var_x = self.var_x

        start_level_t = self.start_level_t
        end_level_t = self.end_level_t
        start_level_q = self.start_level_q
        end_level_q = self.end_level_q

        lines2 = open("/Users/yi-linghwong/Documents/postdoc/publications/paracon_scm_paper/_DATA/_SAM/matrix_X/X.csv","r").readlines()

        x_anomalies = []

        for line in lines2:
            spline = line.rstrip("\n").split(",")
            spline = [float(s) for s in spline]
            x_anomalies.append(spline)

        X_matrix_complete = np.array(x_anomalies)

        if perturb_t == 1 and var_x == "T_PHY":
            X_anomaly_matrix = X_matrix_complete[:(end_level_q - start_level_q + 1), :(end_level_q - start_level_q + 1)]
        elif perturb_t == 1 and var_x == "QVAPOR":
            X_anomaly_matrix = X_matrix_complete[(end_level_t - start_level_t + 1):, :(end_level_q - start_level_q + 1)]
        elif perturb_q == 1 and var_x == "T_PHY":
            X_anomaly_matrix = X_matrix_complete[:(end_level_q - start_level_q + 1), (end_level_t - start_level_t + 1):]
        elif perturb_q == 1 and var_x == "QVAPOR":
            X_anomaly_matrix = X_matrix_complete[(end_level_t - start_level_t + 1):, (end_level_t - start_level_t + 1):]

        print ("Shape X_anomaly_matrix:",X_anomaly_matrix.shape)

        return (X_anomaly_matrix)


    def get_Y_tendency(self):

        """
        Get tendency matrix
        """

        perturb_t = self.perturb_t
        perturb_q = self.perturb_q
        start_level_t = self.start_level_t
        end_level_t = self.end_level_t
        start_level_q = self.start_level_q
        end_level_q = self.end_level_q


        lines3 = open("/Users/yi-linghwong/Documents/postdoc/publications/paracon_scm_paper/_DATA/_SAM/matrix_X/dX.csv","r").readlines()

        y_tendencies = []

        for line in lines3:
            spline = line.rstrip("\n").split(",")
            spline = [-float(s) for s in spline]
            y_tendencies.append(spline)

        Y_tendency_matrix_complete = np.array(y_tendencies)

        if perturb_t == 1:
            Y_tendency_matrix_extracted = Y_tendency_matrix_complete[:(end_level_q-start_level_q+1),:(end_level_q-start_level_q+1)]
        elif perturb_q == 1:
            Y_tendency_matrix_extracted = Y_tendency_matrix_complete[(end_level_t-start_level_t+1):,(end_level_t-start_level_t+1):]


        print ()
        print ("Y_tendency_matrix_complete shape:", Y_tendency_matrix_complete.shape)
        print ("Y_tendency_matrix_extracted shape:", Y_tendency_matrix_extracted.shape)

        return (Y_tendency_matrix_complete,Y_tendency_matrix_extracted)

    def get_power_input(self):

        """
        Get power input
        """

        pressures = self.pressures
        perturb_t = self.perturb_t
        perturb_q = self.perturb_q
        standardise_kuang = self.standardise_kuang

        Y_tendency_matrix_extracted = self.get_Y_tendency()[1]

        psfc = 1014.8

        if perturb_t == 1:
            power_input_list = get_power_input_T(pressures,Y_tendency_matrix_extracted,psfc)
        elif perturb_q == 1:
            power_input_list = get_power_input_q(pressures,Y_tendency_matrix_extracted,psfc)

        if standardise_kuang == 1:
            heat_kuang = get_power_input_kuang(perturb_t,perturb_q)

            if len(heat_kuang) == len(power_input_list):
                standardised_heat = [ p / h for p,h in zip(power_input_list,heat_kuang)]
                power_input_list = standardised_heat

            else:
                print ("LENGTH OF HEAT KUANG AND HEAT MODEL NOT EQUAL, EXITING...")
                sys.exit()

        print ()
        print ("Heat input extracted length:",len(power_input_list))
        print (power_input_list)

        return power_input_list


    def get_M_inv(self):

        """
        Get M^-1
        M_inv = X * Y^-1
        """

        power_input_list = self.get_power_input()
        power_input_array = np.array(power_input_list)
        X_anomaly_matrix = self.get_X_anomaly()

        power_input_array = np.array(power_input_list)

        # get number of columns of X matrix, which will be the size of Y_inv

        num_col_X = np.size(X_anomaly_matrix,1)
        power_input_array = power_input_array[:num_col_X]

        M_inv = get_m_inverse(X_anomaly_matrix,power_input_array)

        print ()
        print ("M_inv shape:",M_inv.shape)

        return (M_inv)

    def plot_matrix_M_inv(self,write_m_inv_to_file,vmax_kuang_t,vmax_kuang_q,vmax_power_t,vmax_power_q):

        """
        PLOT MATRIX M_INV (NORMALISED)
        1. Normalised by standard Kuang's heat input (in K or /gkg) (standardised_kuang = 1)
        2. Normlised by model power input (in K/(Wm^-2) or g/kg/(Wm^2)) (standardised_kuang = 0)
        """

        perturb_t = self.perturb_t
        perturb_q = self.perturb_q
        var_x = self.var_x
        pressures_x = self.pressures_x
        pressures_y = self.pressures_y
        M_inv = self.get_M_inv()
        standardise_kuang = self.standardise_kuang

        if standardise_kuang == 1:

            vmax_t = vmax_kuang_t
            vmax_q = vmax_kuang_q

            if perturb_t == 1 and var_x == "T_PHY":
                plt_title = "T' to dT/dt perturb. [K]"
                file_name = "t_dtdt_norm_kuang"
                vmax = vmax_t
                vmin = -vmax
            elif perturb_t == 1 and var_x == "QVAPOR":
                plt_title = "q' to dT/dt perturb. [g/kg]"
                file_name = "q_dtdt_norm_kuang"
                vmax = vmax_q
                vmin = -vmax
            elif perturb_q == 1 and var_x == "T_PHY":
                plt_title = "T' to dq/dt perturb. [K]"
                file_name = "t_dqdt_norm_kuang"
                vmax = vmax_t
                vmin = -vmax
            elif perturb_q == 1 and var_x == "QVAPOR":
                plt_title = "q' to dq/dt perturb. [g/kg]"
                file_name = "q_dqdt_norm_kuang"
                vmax = vmax_q
                vmin = -vmax

        elif standardise_kuang == 0:

            vmax_t = vmax_power_t
            vmax_q = vmax_power_q

            if perturb_t == 1 and var_x == "T_PHY":
                plt_title = "T' to dT/dt perturb. [K/(W m$\mathregular{^{-2}}$)]"
                file_name = "t_dtdt_norm_power"
                vmax = vmax_t
                vmin = -vmax
            elif perturb_t == 1 and var_x == "QVAPOR":
                plt_title = "q' to dT/dt perturb. [g/kg/(W m$\mathregular{^{-2}}$)]"
                file_name = "q_dtdt_norm_power"
                vmax = vmax_q
                vmin = -vmax
            elif perturb_q == 1 and var_x == "T_PHY":
                plt_title = "T' to dq/dt perturb. [K/(W m$\mathregular{^{-2}}$)]"
                file_name = "t_dqdt_norm_power"
                vmax = vmax_t
                vmin = -vmax
            elif perturb_q == 1 and var_x == "QVAPOR":
                plt_title = "q' to dq/dt perturb. [g/kg/(W m$\mathregular{^{-2}}$)]"
                file_name = "q_dqdt_norm_power"
                vmax = vmax_q
                vmin = -vmax

        plot_title = "M$\mathregular{^{-1}}$, SAM " + plt_title

        print()
        print("M INVERSE SHAPE:", M_inv.shape)
        print("MAX value of M-1:", np.amax(M_inv))
        print("MIN value of M-1:", np.amin(M_inv))

        plot_matrix(M_inv,pressures_x,pressures_y,plot_title,vmin,vmax)

        if write_m_inv_to_file == 1:

            file_name = "sam_" + file_name

            write_matrix_to_file(M_inv,"_SAM",file_name)



    def plot_matrix_X_raw(self,vmax_t,vmax_q):

        """
        PLOT MATRIX X RAW (UNNORMALISED)
        """

        print ()
        print ("#######################")
        print ("PLOT X RAW")
        print ()

        perturb_t = self.perturb_t
        perturb_q = self.perturb_q
        var_x = self.var_x
        pressures_x = self.pressures_x
        pressures_y = self.pressures_y

        start_level_t = self.start_level_t
        end_level_t = self.end_level_t
        start_level_q = self.start_level_q
        end_level_q = self.end_level_q

        #----------------------------------------------
        # GET raw X

        lines2 = open("/Users/yi-linghwong/Documents/postdoc/publications/paracon_scm_paper/_DATA/_SAM/X.csv","r").readlines()

        x_anomalies = []

        for line in lines2:
            spline = line.rstrip("\n").split(",")
            spline = [float(s) for s in spline]
            x_anomalies.append(spline)

        X_matrix_complete = np.array(x_anomalies)

        if perturb_t == 1 and var_x == "T_PHY":
            X_anomaly_matrix = X_matrix_complete[:(end_level_q - start_level_q + 1), :(end_level_q - start_level_q + 1)]
        elif perturb_t == 1 and var_x == "QVAPOR":
            X_anomaly_matrix = X_matrix_complete[(end_level_t - start_level_t + 1):, :(end_level_q - start_level_q + 1)]
        elif perturb_q == 1 and var_x == "T_PHY":
            X_anomaly_matrix = X_matrix_complete[:(end_level_q - start_level_q + 1), (end_level_t - start_level_t + 1):]
        elif perturb_q == 1 and var_x == "QVAPOR":
            X_anomaly_matrix = X_matrix_complete[(end_level_t - start_level_t + 1):, (end_level_t - start_level_t + 1):]

        print ("Shape X_anomaly_matrix:",X_anomaly_matrix.shape)

        #--------------------------------------------------------

        vmax_t = vmax_t
        vmax_q = vmax_q

        if perturb_t == 1 and var_x == "T_PHY":
            plt_title = "T' to dT/dt perturb. [K]"
            vmax = vmax_t
            vmin = -vmax
        elif perturb_t == 1 and var_x == "QVAPOR":
            plt_title = "q' to dT/dt perturb. [g/kg]"
            vmax = vmax_q
            vmin = -vmax
        elif perturb_q == 1 and var_x == "T_PHY":
            plt_title = "T' to dq/dt perturb. [K]"
            vmax = vmax_t
            vmin = -vmax
        elif perturb_q == 1 and var_x == "QVAPOR":
            plt_title = "q' to dq/dt perturb. [g/kg]"
            vmax = vmax_q
            vmin = -vmax

        plot_title = "X, SAM " + plt_title

        print()
        print("X RAW anomaly SHAPE:", X_anomaly_matrix.shape)
        print("MAX value of X:", np.amax(X_anomaly_matrix))
        print("MIN value of X:", np.amin(X_anomaly_matrix))

        plot_matrix(X_anomaly_matrix,pressures_x,pressures_y,plot_title,vmin,vmax)


    def plot_anomaly_profile(self,write_anomaly_to_file,target_level_1,target_level_2,label_level_1,label_level_2):

        """
        PLOT anomaly profiles (2 levels)
        Optionally write profiles to file
        """

        standardise_kuang = self.standardise_kuang
        perturb_t = self.perturb_t
        perturb_q = self.perturb_q
        var_x = self.var_x
        t_amplitude = self.t_amplitude
        q_amplitude = self.q_amplitude
        M_inv = self.get_M_inv()
        X_anomaly_matrix = self.get_X_anomaly()
        pressures_y = self.pressures_y


        m_profile_1 = M_inv[:,target_level_1-1]
        x_profile_1 = X_anomaly_matrix[:,target_level_1-1]
        m_profile_2 = M_inv[:,target_level_2-1]
        x_profile_2 = X_anomaly_matrix[:,target_level_2-1]

        print ("M^-1 profile length:", len(m_profile_1), len(m_profile_2))
        print ("X profile length:", len(x_profile_1), len(x_profile_2))


        plot_title = "SAM"

        plot_profile(standardise_kuang,pressures_y,x_profile_1,x_profile_2,m_profile_1,m_profile_2,plot_title,label_level_1,label_level_2,var_x)

        if write_anomaly_to_file == 1:

            if t_amplitude == 0.5:
                t_amp = "05"
            elif t_amplitude == 0.2:
                t_amp = "02"

            if q_amplitude == 0.2:
                q_amp = "02"
            elif q_amplitude == 0.1:
                q_amp = "01"

            if perturb_t == 1 and var_x == "T_PHY":
                file_name = "SAM_T_DTDT_" + t_amp
            elif perturb_t == 1 and var_x == "QVAPOR":
                file_name = "SAM_Q_DTDT_" + t_amp
            elif perturb_q == 1 and var_x == "T_PHY":
                file_name = "SAM_T_DQDT_" + q_amp
            elif perturb_q == 1 and var_x == "QVAPOR":
                file_name = "SAM_Q_DQDT_" + q_amp

            write_profiles_to_file(standardise_kuang,pressures_y,x_profile_1,x_profile_2,m_profile_1,m_profile_2,"_SAM",file_name,
                          label_level_1,label_level_2)

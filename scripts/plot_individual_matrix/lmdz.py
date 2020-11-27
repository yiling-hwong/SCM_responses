#!/usr/bin/env python3

"""
Created on 13 August 2020
@author: Yi-Ling HWONG
"""

import sys
import numpy as np
from utils import *

class LMDZ():

    def __init__(self,lmdz_version,standardise_kuang,perturb_t,perturb_q,state_anomaly,t_amplitude,q_amplitude):

        self.lmdz_version = lmdz_version
        self.standardise_kuang = standardise_kuang
        self.perturb_t = perturb_t
        self.perturb_q = perturb_q
        self.state_anomaly = state_anomaly
        self.t_amplitude = t_amplitude
        self.q_amplitude = q_amplitude

        self.start_level_t = 1
        self.start_level_q = 1
        self.num_levels = 42 # num level up to 200 hPa

        print ()
        print ("LMDZ")
        print ()

        lines1 = open("../../data/LMDZ/matrix_X_raw/pressures","r").readlines()

        if standardise_kuang == True:
            lines3 = open("../../data/LMDZ/matrix_X_raw/level_index_kuang.csv","r").readlines()

            self.end_level_t = 18
            self.end_level_q = 18

        elif standardise_kuang == False:
            lines3 = open("../../data/LMDZ/matrix_X_raw/level_index.csv","r").readlines()

            self.end_level_t = 42
            self.end_level_q = 42

        pressures = []
        level_index = []

        for line in lines3:
            spline = line.rstrip("\n")
            level_index.append(int(spline))

        print("Length level index:", len(level_index))
        print(level_index)
        print()

        self.level_idx = [x - 1 for x in level_index]
        print("Level index - 1:", self.level_idx)
        print()

        for line in lines1:
            spline = line.rstrip("\n")
            pressures.append(float(spline) * 100)

        if pressures[0] < pressures[1]:
            pressures.reverse()

        self.pressures = pressures[:self.num_levels]

        print("Length pressures (up to 200 hPa):", len(self.pressures))
        print(self.pressures)

        pressures_extracted = []
        for index, value in enumerate(self.pressures):
            if index in self.level_idx:
                pressures_extracted.append(value / 100)

        self.pressures_t = pressures_extracted[:self.end_level_t]
        self.pressures_q = pressures_extracted[:self.end_level_q]

        self.pressures_x = self.pressures_t
        self.pressures_y = self.pressures_q


    def get_X_anomaly(self):

        """
        Get anomaly matrix
        """

        lmdz_version = self.lmdz_version
        perturb_t = self.perturb_t
        perturb_q = self.perturb_q
        state_anomaly = self.state_anomaly
        t_amplitude = self.t_amplitude
        q_amplitude = self.q_amplitude

        # ----------------------

        if perturb_t == True and self.t_amplitude == 0.5:
            perturbation = "DTDT"
            perturbation_amplitude = "05"
        elif perturb_t == True and self.t_amplitude == 0.2:
            perturbation = "DTDT"
            perturbation_amplitude = "02"
        elif perturb_q == True and self.q_amplitude == 0.2:
            perturbation = "DQDT"
            perturbation_amplitude = "02"
        elif perturb_q == True and self.q_amplitude == 0.1:
            perturbation = "DQDT"
            perturbation_amplitude = "01"

        if state_anomaly == "T":
            var = "T"
        elif state_anomaly == "q":
            var = "Q"

        if lmdz_version == "5a":
            lines = open("../../data/LMDZ/matrix_X_raw/lmdz5a/" + var + "_" + perturbation + "_" + perturbation_amplitude + ".csv","r").readlines()
        elif lmdz_version == "6a":
            lines = open("../../data/LMDZ/matrix_X_raw/lmdz6a/" + var + "_" + perturbation + "_" + perturbation_amplitude + ".csv","r").readlines()
        elif lmdz_version == "6ab":
            lines = open("../../data/LMDZ/matrix_X_raw/lmdz6ab/" + var + "_" + perturbation + "_" + perturbation_amplitude + ".csv","r").readlines()

        x_anomalies = []

        for line in lines:
            spline = line.rstrip("\n").split(",")

            if perturb_t == True:
                spline = [float(x) * (0.5 / t_amplitude) for x in spline]

            if perturb_q == True:
                spline = [float(x) * (0.2 / q_amplitude) for x in spline]

            x_anomalies.append(spline)

        X_anomaly_matrix = np.array(x_anomalies)

        # get matrix at selected index

        row_idx = np.array(self.level_idx)
        col_idx = np.array(self.level_idx)

        X_anomaly_extracted = X_anomaly_matrix[row_idx[:,None], col_idx] # extract row and column by index

        print ()
        print ("Shape X_anomaly_matrix:", X_anomaly_extracted.shape)

        return (X_anomaly_extracted)


    def get_Y_tendency(self):

        """
        Get tendency matrix
        """

        perturb_t = self.perturb_t
        perturb_q = self.perturb_q

        if perturb_t == True:
            Y_tendency_matrix_complete = get_perturbations_T(self.pressures)
        elif perturb_q == True:
            Y_tendency_matrix_complete = get_perturbations_q(self.pressures)

        row_idx = np.array(self.level_idx)
        col_idx = np.array(self.level_idx)

        Y_tendency_matrix_extracted = Y_tendency_matrix_complete[row_idx[:, None], col_idx]  # extract row and column by index

        print ()
        print("Y_tendency_matrix_complete shape:", Y_tendency_matrix_complete.shape)
        print("Y_tendency_matrix_extracted shape:", Y_tendency_matrix_extracted.shape)

        return (Y_tendency_matrix_complete, Y_tendency_matrix_extracted)

    def get_power_input(self):

        """
        Get power input
        """

        pressures = self.pressures
        perturb_t = self.perturb_t
        perturb_q = self.perturb_q
        standardise_kuang = self.standardise_kuang
        level_idx = self.level_idx


        Y_tendency_matrix_complete = self.get_Y_tendency()[0]

        psfc = 1014.8

        if perturb_t == True:
            power_input_list = get_power_input_T(pressures,Y_tendency_matrix_complete,psfc)
        elif perturb_q == True:
            power_input_list = get_power_input_q(pressures,Y_tendency_matrix_complete,psfc)


        power_input_extracted = []
        for index,value in enumerate(power_input_list):
            if index in level_idx:
                power_input_extracted.append(value)

        if standardise_kuang == True:
            heat_kuang = get_power_input_kuang(perturb_t,perturb_q)

            if len(heat_kuang) == len(power_input_extracted):
                standardised_heat = [ p / h for p,h in zip(power_input_extracted,heat_kuang)]
                power_input_extracted = standardised_heat

            else:
                print ("LENGTH OF HEAT KUANG AND HEAT MODEL NOT EQUAL, EXITING...")
                sys.exit()

        print ()
        print ("FINAL Heat input extracted length:",len(power_input_extracted))
        print (power_input_extracted)

        return power_input_extracted


    def get_M_inv(self):

        """
        Get M^-1
        M_inv = X * Y^-1
        """

        power_input_extracted = self.get_power_input()
        power_input_array = np.array(power_input_extracted)
        X_anomaly_matrix = self.get_X_anomaly()

        # get number of columns of X matrix, which will be the size of Y_inv

        num_col_X = np.size(X_anomaly_matrix, 1)
        power_input_array = power_input_array[:num_col_X]

        M_inv = get_m_inverse(X_anomaly_matrix, power_input_array)

        print()
        print("M_inv shape:", M_inv.shape)

        return (M_inv)

    def plot_matrix_M_inv(self,write_m_inv_to_file,vmax_kuang_t,vmax_kuang_q,vmax_power_t,vmax_power_q):

        """
        PLOT MATRIX M_INV (NORMALISED)
        1. Normalised by standard Kuang's heat input (in K or g/kg) (standardised_kuang = True)
        2. Normlised by model power input (in K/(Wm^-2) or g/kg/(Wm^2)) (standardised_kuang = False)
        """

        lmdz_version = self.lmdz_version
        perturb_t = self.perturb_t
        perturb_q = self.perturb_q
        state_anomaly = self.state_anomaly
        pressures_x = self.pressures_x
        pressures_y = self.pressures_y
        M_inv = self.get_M_inv()
        standardise_kuang = self.standardise_kuang

        if standardise_kuang == True:

            vmax_t = vmax_kuang_t
            vmax_q = vmax_kuang_q

            if perturb_t == True and state_anomaly == "T":
                plt_title = "T' to dT/dt perturb. [K]"
                file_name = "t_dtdt_norm_kuang"
                vmax = vmax_t
                vmin = -vmax
            elif perturb_t == True and state_anomaly == "q":
                plt_title = "q' to dT/dt perturb. [g/kg]"
                file_name = "q_dtdt_norm_kuang"
                vmax = vmax_q
                vmin = -vmax
            elif perturb_q == True and state_anomaly == "T":
                plt_title = "T' to dq/dt perturb. [K]"
                file_name = "t_dqdt_norm_kuang"
                vmax = vmax_t
                vmin = -vmax
            elif perturb_q == True and state_anomaly == "q":
                plt_title = "q' to dq/dt perturb. [g/kg]"
                file_name = "q_dqdt_norm_kuang"
                vmax = vmax_q
                vmin = -vmax

        elif standardise_kuang == False:

            vmax_t = vmax_power_t
            vmax_q = vmax_power_q

            if perturb_t == True and state_anomaly == "T":
                plt_title = "T' to dT/dt perturb. [K/(W m$\mathregular{^{-2}}$)]"
                file_name = "t_dtdt_norm_power"
                vmax = vmax_t
                vmin = -vmax
            elif perturb_t == True and state_anomaly == "q":
                plt_title = "q' to dT/dt perturb. [g/kg/(W m$\mathregular{^{-2}}$)]"
                file_name = "q_dtdt_norm_power"
                vmax = vmax_q
                vmin = -vmax
            elif perturb_q == True and state_anomaly == "T":
                plt_title = "T' to dq/dt perturb. [K/(W m$\mathregular{^{-2}}$)]"
                file_name = "t_dqdt_norm_power"
                vmax = vmax_t
                vmin = -vmax
            elif perturb_q == True and state_anomaly == "q":
                plt_title = "q' to dq/dt perturb. [g/kg/(W m$\mathregular{^{-2}}$)]"
                file_name = "q_dqdt_norm_power"
                vmax = vmax_q
                vmin = -vmax

        if lmdz_version == "5a":
            plot_title = "M$\mathregular{^{-1}}$, LMDZ5A " + plt_title
        elif lmdz_version == "6a":
            plot_title = "M$\mathregular{^{-1}}$, LMDZ6A " + plt_title
        elif lmdz_version == "6ab":
            plot_title = "M$\mathregular{^{-1}}$, LMDZ6Ab " + plt_title

        print()
        print("M INVERSE SHAPE:", M_inv.shape)
        print("MAX value of M-1:", np.amax(M_inv))
        print("MIN value of M-1:", np.amin(M_inv))

        plot_matrix(M_inv,pressures_x,pressures_y,plot_title,vmin,vmax)

        if write_m_inv_to_file == 1:

            if lmdz_version == "5a":
                file_name = "lmdz5a_" + file_name
                write_matrix_to_file(M_inv,"LMDZ",file_name)
            elif lmdz_version == "6a":
                file_name = "lmdz6a_" + file_name
                write_matrix_to_file(M_inv,"LMDZ",file_name)
            elif lmdz_version == "6ab":
                file_name = "lmdz6ab_" + file_name
                write_matrix_to_file(M_inv,"LMDZ",file_name)

    def plot_anomaly_profile(self,write_anomaly_to_file,target_level_1,target_level_2,label_level_1,label_level_2):

        """
        PLOT anomaly profiles (2 levels)
        Optionally write profiles to file
        """

        lmdz_version = self.lmdz_version
        standardise_kuang = self.standardise_kuang
        perturb_t = self.perturb_t
        perturb_q = self.perturb_q
        state_anomaly = self.state_anomaly
        t_amplitude = self.t_amplitude
        q_amplitude = self.q_amplitude
        level_idx = self.level_idx
        M_inv = self.get_M_inv()
        X_anomaly_matrix = self.get_X_anomaly()
        pressures_y = self.pressures_y

        print ()
        print ("########################")
        print ("Plotting X profiles ...")
        print ("Level index:",level_idx)

        for index,value in enumerate(level_idx):
            if value == target_level_1-1:
                target_level1 = index
            if value == target_level_2-1:
                target_level2 = index

        print ("TARGET LEVEL indices:",target_level1,target_level2)
        print ("TARGET LEVELS:",level_idx[target_level1]+1,level_idx[target_level2]+1)
        print ()

        m_profile_1 = M_inv[:,target_level1]
        x_profile_1 = X_anomaly_matrix[:,target_level1]
        m_profile_2 = M_inv[:,target_level2]
        x_profile_2 = X_anomaly_matrix[:,target_level2]

        print ("M^-1 profile length:", len(m_profile_1), len(m_profile_2))
        print ("X profile length:", len(x_profile_1), len(x_profile_2))

        if lmdz_version == "5a":
            plot_title = "LMDZ5A"
        elif lmdz_version == "6a":
            plot_title = "LMDZ6A"
        elif lmdz_version == "6ab":
            plot_title = "LMDZ6Ab"

        plot_profile(standardise_kuang,pressures_y,x_profile_1,x_profile_2,m_profile_1,m_profile_2,plot_title,label_level_1,label_level_2,state_anomaly)

        if write_anomaly_to_file == True:

            if t_amplitude == 0.5:
                t_amp = "05"
            elif t_amplitude == 0.2:
                t_amp = "02"

            if q_amplitude == 0.2:
                q_amp = "02"
            elif q_amplitude == 0.1:
                q_amp = "01"

            if perturb_t == True and state_anomaly == "T":
                file_name = "LMDZ_T_DTDT_" + t_amp
            elif perturb_t == True and state_anomaly == "q":
                file_name = "LMDZ_Q_DTDT_" + t_amp
            elif perturb_q == True and state_anomaly == "T":
                file_name = "LMDZ_T_DQDT_" + q_amp
            elif perturb_q == True and state_anomaly == "q":
                file_name = "LMDZ_Q_DQDT_" + q_amp

            write_profiles_to_file(standardise_kuang,pressures_y,x_profile_1,x_profile_2,m_profile_1,m_profile_2,"LMDZ",file_name,
                          label_level_1,label_level_2)

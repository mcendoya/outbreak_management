# -*- coding: utf-8 -*-
"""
Code to perform the simulations

@author: Martina Cendoya
"""

import model_functions as fn
import pandas as pd
import numpy as np
from numba import set_num_threads
set_num_threads(8)

# Individuals data
xy = pd.read_csv("coords.csv")
coords = xy[['x', 'y']].values

##################################################################################

# Initial outbreak (t=0) at 5 foci, 10 infected per focus, each focus 1km radious

# The dataframe already contains the initial outbreak 
# s_t0 = fn.aggregate_introd(nI = 10, nT = len(xy), coords = coords, f = 5, f_rad = 1000)
# xy['t0'] = s_t0

##################################################################################

xy['t0'][xy['t0']==1] = 0.015 # We consider the initial infected as asymptomatic infected (Ia)
alm_cells = xy['cell_ha'] # cell in which each individual is located

# Folder to save results
res_v = "results"

##################################################################################

# Simulations of disease spread without intervention (BASELINE SCENARIO)

# Spread time, time interval to save
t_max, t_save = 30*12, 6

# Initial status (t=0)
s_init = np.array(xy['t0'])
# Spread parameters
beta, rg, nu, sigma = 0.015, 800, 1, 8
# cell in which each individual is located (to numpy.array)
alm_cells = np.array(alm_cells)

for s in range(1,51): # 50 simulations
    s1 = fn.stat(s_init, alm_cells, coords, beta, rg, nu, sigma, t_max, t_save)
    np.savetxt("./"+str(res_v)+"/no_intervention/spread"+str(s)+".txt", np.transpose(s1), fmt='%f')

##################################################################################

# Simulations with the implementation of outbreak management plans

# Folders must first be created in the directory with the following structure: 
    # res_v: name of the folder to store the results
        # surv: name of the folder to store the results with the surveillance approach 'surv'
            # ztXX_eYY_bZZ: Name of the folder inside 'surv' with buffer zone XX, 
            # eradication radius YY and reduction of the transmission rate ZZ 
            # (replace XX, YY, ZZ by the corresponding values)
            
# Extracts the results of disease spread (each 6t) and surveillance (each 12t)

##################################################################################

# Spread time
t_max=30*12

# Buffer zone sizes
r_zt = [2500, 5000]

# Eradication radii
r_er = [50, 100]
##################################################################################

# One-step suirveillance approach, CL_low
surv = "one_step_00"

# no reduction in the transmission rate in the BZ
b = 0 
for zt in r_zt:
    for er in r_er:
        for s in range(1,51): # 50 simulations
            parameters = {'spread': {'beta':0.015, 'rg':800, 'nu':1, 'sigma':8},
                          'survey_d': {'cl':0.8, 'dp':0.01, 'ms':0.8},
                          'survey_zt': {'cl':0.9, 'dp':0.01, 'ms':0.8},
                          'control': {'r_zt':zt+50, 'r_er':er}}
            s_final, s_survey = fn.strategies(xy, parameters, t_max)
            np.savetxt("./"+str(res_v)+"/"+str(surv)+"/zt"+str(zt)+"_e"+str(er)+"_b"+str(b)+"/spread"+str(s)+".txt", np.transpose(s_final), fmt='%f')
            s_survey1 = s_survey.astype(str)
            np.savetxt("./"+str(res_v)+"/"+str(surv)+"/zt"+str(zt)+"_e"+str(er)+"_b"+str(b)+"/surveys"+str(s)+".txt", np.transpose(s_survey1), fmt='%s')
            
# 50% reduction in the transmission rate in the BZ 
b = 50 
for zt in r_zt:
    for er in r_er:
        for s in range(1,51): # 50 simulations
            parameters = {'spread': {'beta':0.015, 'beta_zt':0.015*0.5, 'rg':800, 'nu':1, 'sigma':8},
                          'survey_d': {'cl':0.8, 'dp':0.01, 'ms':0.8},
                          'survey_zt': {'cl':0.9, 'dp':0.01, 'ms':0.8},
                          'control': {'r_zt':zt+50, 'r_er':er}}
            s_final, s_survey = fn.strategies(xy, parameters, t_max, additional_control=True)
            np.savetxt("./"+str(res_v)+"/"+str(surv)+"/zt"+str(zt)+"_e"+str(er)+"_b"+str(b)+"/spread"+str(s)+".txt", np.transpose(s_final), fmt='%f')
            s_survey1 = s_survey.astype(str)
            np.savetxt("./"+str(res_v)+"/"+str(surv)+"/zt"+str(zt)+"_e"+str(er)+"_b"+str(b)+"/surveys"+str(s)+".txt", np.transpose(s_survey1), fmt='%s')

# 90% reduction in the transmission rate in the BZ  
b = 90 
for zt in r_zt:
    for er in r_er:
        for s in range(1,51): # 50 simulations
            parameters = {'spread': {'beta':0.015, 'beta_zt':0.015*0.1, 'rg':800, 'nu':1, 'sigma':8},
                          'survey_d': {'cl':0.8, 'dp':0.01, 'ms':0.8},
                          'survey_zt': {'cl':0.9, 'dp':0.01, 'ms':0.8},
                          'control': {'r_zt':zt+50, 'r_er':er}}
            s_final, s_survey = fn.strategies(xy, parameters, t_max, additional_control=True)
            np.savetxt("./"+str(res_v)+"/"+str(surv)+"/zt"+str(zt)+"_e"+str(er)+"_b"+str(b)+"/spread"+str(s)+".txt", np.transpose(s_final), fmt='%f')
            s_survey1 = s_survey.astype(str)
            np.savetxt("./"+str(res_v)+"/"+str(surv)+"/zt"+str(zt)+"_e"+str(er)+"_b"+str(b)+"/surveys"+str(s)+".txt", np.transpose(s_survey1), fmt='%s')
##################################################################################

# One-step suirveillance approach, CL_high  
surv = "one_step_01"

# no reduction in the transmission rate in the BZ   
b = 0
for zt in r_zt:
    for er in r_er:
        for s in range(1,51): # 50 simulations
            parameters = {'spread': {'beta':0.015, 'rg':800, 'nu':1, 'sigma':8},
                          'survey_d': {'cl':0.99, 'dp':0.01, 'ms':0.8},
                          'survey_zt': {'cl':0.99, 'dp':0.01, 'ms':0.8},
                          'control': {'r_zt':zt+50, 'r_er':er}}
            s_final, s_survey = fn.strategies(xy, parameters, t_max)
            np.savetxt("./"+str(res_v)+"/"+str(surv)+"/zt"+str(zt)+"_e"+str(er)+"_b"+str(b)+"/spread"+str(s)+".txt", np.transpose(s_final), fmt='%f')
            s_survey1 = s_survey.astype(str)
            np.savetxt("./"+str(res_v)+"/"+str(surv)+"/zt"+str(zt)+"_e"+str(er)+"_b"+str(b)+"/surveys"+str(s)+".txt", np.transpose(s_survey1), fmt='%s')

# 50% reduction in the transmission rate in the BZ 
b = 50
for zt in r_zt:
    for er in r_er:
        for s in range(1,51): # 50 simulations
            parameters = {'spread': {'beta':0.015,'beta_zt':0.015*0.5,  'rg':800, 'nu':1, 'sigma':8},
                          'survey_d': {'cl':0.99, 'dp':0.01, 'ms':0.8},
                          'survey_zt': {'cl':0.99, 'dp':0.01, 'ms':0.8},
                          'control': {'r_zt':zt+50, 'r_er':er}}
            s_final, s_survey = fn.strategies(xy, parameters, t_max, additional_control=True)
            np.savetxt("./"+str(res_v)+"/"+str(surv)+"/zt"+str(zt)+"_e"+str(er)+"_b"+str(b)+"/spread"+str(s)+".txt", np.transpose(s_final), fmt='%f')
            s_survey1 = s_survey.astype(str)
            np.savetxt("./"+str(res_v)+"/"+str(surv)+"/zt"+str(zt)+"_e"+str(er)+"_b"+str(b)+"/surveys"+str(s)+".txt", np.transpose(s_survey1), fmt='%s')

# 90% reduction in the transmission rate in the BZ 
b = 90
for zt in r_zt:
    for er in r_er:
        for s in range(1,51): # 50 simulations
            parameters = {'spread': {'beta':0.015,'beta_zt':0.015*0.1,  'rg':800, 'nu':1, 'sigma':8},
                          'survey_d': {'cl':0.99, 'dp':0.01, 'ms':0.8},
                          'survey_zt': {'cl':0.99, 'dp':0.01, 'ms':0.8},
                          'control': {'r_zt':zt+50, 'r_er':er}}
            s_final, s_survey = fn.strategies(xy, parameters, t_max, additional_control=True)
            np.savetxt("./"+str(res_v)+"/"+str(surv)+"/zt"+str(zt)+"_e"+str(er)+"_b"+str(b)+"/spread"+str(s)+".txt", np.transpose(s_final), fmt='%f')
            s_survey1 = s_survey.astype(str)
            np.savetxt("./"+str(res_v)+"/"+str(surv)+"/zt"+str(zt)+"_e"+str(er)+"_b"+str(b)+"/surveys"+str(s)+".txt", np.transpose(s_survey1), fmt='%s')

##################################################################################

# Two-step suirveillance approach, CL_high
surv = "two_step_00"

# no reduction in the transmission rate in the BZ   
b = 0
for zt in r_zt:
    for er in r_er:
        for s in range(1,51): # 50 simulations
            parameters = {'spread': {'beta':0.015, 'rg':800, 'nu':1, 'sigma':8},
                          'survey_d': {'step1': {'dp':0.01, 'ms':0.8},
                                       'step2': {'cl':0.8, 'dp':0.01}},
                          'survey_zt': {'step1': {'dp':0.01, 'ms':0.8},
                                       'step2': {'cl':0.9, 'dp':0.01}},
                          'control': {'r_zt':zt+50, 'r_er':er}}
            s_final, s_survey = fn.two_steps_strategy(xy, parameters, t_max, additional_control=False)
            np.savetxt("./"+str(res_v)+"/"+str(surv)+"/zt"+str(zt)+"_e"+str(er)+"_b"+str(b)+"/spread"+str(s)+".txt", np.transpose(s_final), fmt='%f')
            s_survey1 = s_survey.astype(str)
            np.savetxt("./"+str(res_v)+"/"+str(surv)+"/zt"+str(zt)+"_e"+str(er)+"_b"+str(b)+"/surveys"+str(s)+".txt", np.transpose(s_survey1), fmt='%s')

# 50% reduction in the transmission rate in the BZ 
b = 50
for zt in r_zt:
    for er in r_er:
        for s in range(1,51): # 50 simulations
            parameters = {'spread': {'beta':0.015, 'beta_zt':0.015*0.5, 'rg':800, 'nu':1, 'sigma':8},
                          'survey_d': {'step1': {'dp':0.01, 'ms':0.8},
                                       'step2': {'cl':0.8, 'dp':0.01}},
                          'survey_zt': {'step1': {'dp':0.01, 'ms':0.8},
                                       'step2': {'cl':0.9, 'dp':0.01}},
                          'control': {'r_zt':zt+50, 'r_er':er}}
            s_final, s_survey = fn.two_steps_strategy(xy, parameters, t_max, additional_control=True)
            np.savetxt("./"+str(res_v)+"/"+str(surv)+"/zt"+str(zt)+"_e"+str(er)+"_b"+str(b)+"/spread"+str(s)+".txt", np.transpose(s_final), fmt='%f')
            s_survey1 = s_survey.astype(str)
            np.savetxt("./"+str(res_v)+"/"+str(surv)+"/zt"+str(zt)+"_e"+str(er)+"_b"+str(b)+"/surveys"+str(s)+".txt", np.transpose(s_survey1), fmt='%s')

# 90% reduction in the transmission rate in the BZ 
b = 90
for zt in r_zt:
    for er in r_er:
        for s in range(1,51): # 50 simulations
            parameters = {'spread': {'beta':0.015, 'beta_zt':0.015*0.1, 'rg':800, 'nu':1, 'sigma':8},
                          'survey_d': {'step1': {'dp':0.01, 'ms':0.8},
                                       'step2': {'cl':0.8, 'dp':0.01}},
                          'survey_zt': {'step1': {'dp':0.01, 'ms':0.8},
                                       'step2': {'cl':0.9, 'dp':0.01}},
                          'control': {'r_zt':zt+50, 'r_er':er}}
            s_final, s_survey = fn.two_steps_strategy(xy, parameters, t_max, additional_control=True)
            np.savetxt("./"+str(res_v)+"/"+str(surv)+"/zt"+str(zt)+"_e"+str(er)+"_b"+str(b)+"/spread"+str(s)+".txt", np.transpose(s_final), fmt='%f')
            s_survey1 = s_survey.astype(str)
            np.savetxt("./"+str(res_v)+"/"+str(surv)+"/zt"+str(zt)+"_e"+str(er)+"_b"+str(b)+"/surveys"+str(s)+".txt", np.transpose(s_survey1), fmt='%s')
            

##################################################################################

# Two-steps suirveillance approach, CL_low
surv = "two_step_01"

# no reduction in the transmission rate in the BZ   
b = 0
for zt in r_zt:
    for er in r_er:
        for s in range(1,51): # 50 simulations
            parameters = {'spread': {'beta':0.015, 'rg':800, 'nu':1, 'sigma':8},
                          'survey_d': {'step1': {'dp':0.01, 'ms':0.8},
                                       'step2': {'cl':0.6, 'dp':0.01}},
                          'survey_zt': {'step1': {'dp':0.01, 'ms':0.8},
                                       'step2': {'cl':0.7, 'dp':0.01}},
                          'control': {'r_zt':zt+50, 'r_er':er}}
            s_final, s_survey = fn.two_steps_strategy(xy, parameters, t_max, additional_control=False)
            np.savetxt("./"+str(res_v)+"/"+str(surv)+"/zt"+str(zt)+"_e"+str(er)+"_b"+str(b)+"/spread"+str(s)+".txt", np.transpose(s_final), fmt='%f')
            s_survey1 = s_survey.astype(str)
            np.savetxt("./"+str(res_v)+"/"+str(surv)+"/zt"+str(zt)+"_e"+str(er)+"_b"+str(b)+"/surveys"+str(s)+".txt", np.transpose(s_survey1), fmt='%s')

# 50% reduction in the transmission rate in the BZ 
b = 50
for zt in r_zt:
    for er in r_er:
        for s in range(1,51): # 50 simulations
            parameters = {'spread': {'beta':0.015, 'beta_zt':0.015*0.5, 'rg':800, 'nu':1, 'sigma':8},
                          'survey_d': {'step1': {'dp':0.01, 'ms':0.8},
                                       'step2': {'cl':0.6, 'dp':0.01}},
                          'survey_zt': {'step1': {'dp':0.01, 'ms':0.8},
                                       'step2': {'cl':0.7, 'dp':0.01}},
                          'control': {'r_zt':zt+50, 'r_er':er}}
            s_final, s_survey = fn.two_steps_strategy(xy, parameters, t_max, additional_control=True)
            np.savetxt("./"+str(res_v)+"/"+str(surv)+"/zt"+str(zt)+"_e"+str(er)+"_b"+str(b)+"/spread"+str(s)+".txt", np.transpose(s_final), fmt='%f')
            s_survey1 = s_survey.astype(str)
            np.savetxt("./"+str(res_v)+"/"+str(surv)+"/zt"+str(zt)+"_e"+str(er)+"_b"+str(b)+"/surveys"+str(s)+".txt", np.transpose(s_survey1), fmt='%s')
            
# 90% reduction in the transmission rate in the BZ 
b = 90
for zt in r_zt:
    for er in r_er:
        for s in range(1,51): # 50 simulations
            parameters = {'spread': {'beta':0.015, 'beta_zt':0.015*0.1, 'rg':800, 'nu':1, 'sigma':8},
                          'survey_d': {'step1': {'dp':0.01, 'ms':0.8},
                                       'step2': {'cl':0.6, 'dp':0.01}},
                          'survey_zt': {'step1': {'dp':0.01, 'ms':0.8},
                                       'step2': {'cl':0.7, 'dp':0.01}},
                          'control': {'r_zt':zt+50, 'r_er':er}}
            s_final, s_survey = fn.two_steps_strategy(xy, parameters, t_max, additional_control=True)
            np.savetxt("./"+str(res_v)+"/"+str(surv)+"/zt"+str(zt)+"_e"+str(er)+"_b"+str(b)+"/spread"+str(s)+".txt", np.transpose(s_final), fmt='%f')
            s_survey1 = s_survey.astype(str)
            np.savetxt("./"+str(res_v)+"/"+str(surv)+"/zt"+str(zt)+"_e"+str(er)+"_b"+str(b)+"/surveys"+str(s)+".txt", np.transpose(s_survey1), fmt='%s')

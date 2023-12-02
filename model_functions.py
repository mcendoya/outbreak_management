# -*- coding: utf-8 -*-
"""
Functions for the simulation of disease spread and for the implementation of 
outbreak management plans.

@author: Martina Cendoya
"""
import numpy as np
import numba as nb
import itertools
import scipy.special as sc
from random import sample
from shapely.geometry import Point
from shapely.ops import unary_union
import geopandas as gpd
# enable shapely.speedups which makes some of the spatial queries running faster
import shapely.speedups
shapely.speedups.enable()
from scipy.spatial.distance import cdist

@nb.njit(fastmath=True)
def dist_eu(coords_i, coords_j):
    """
    Euclidean distance function between two individuals i and j

    Parameters
    ----------
    coords_i : array with coordinates [x, y] of individual i.
    coords_j : array with coordinates [x, y] of individual j.

    Returns
    -------
    d: distance.

    """
    s = 0.0
    for i in range(2):
        tmp = coords_i[i] - coords_j[i]
        s += tmp * tmp
    d = s**0.5
    return float(d)

@nb.njit(fastmath=True)
def cor_matern(dist, rg, nu):
    """
    Matérn correlation

    Parameters
    ----------
    dist : Euclidean distance.
    rg : range parameter.
    nu : smoothness parameter.

    Returns
    -------
    cor : correlation.

    """
    nu = float(nu)
    phi = (np.sqrt(8*nu)/rg)
    a = (2**(nu - 1) * sc.gamma(nu))**(-1)
    b = (dist * phi)**nu * sc.kv(nu, dist * phi)
    cor = a * b
    return cor

def nb_dist(D_MAX, coords, indI, indS):
    """
    Neighbors within a distance radius D_MAX

    Parameters
    ----------
    D_MAX : Maximum distance at which neighbors are considered.
    coords : array with coordinates [x, y] of individuals.
    indI : Index of the individuals from whom their neighbors are to be collected.
    indS : Index of the individuals that are evaluated if they are neighbors of indI 
        (i.e., if they are at a distance <= D_MAX).

    Returns
    -------
    nb_0 : List with the indS individuals that are at a distance <= D_MAX for each indI.

    """
    nb_0 = []
    for i in indI:
        nb_i = []
        for j in indS:
            dist = dist_eu(coords[i], coords[j])
            if dist <= D_MAX:
                nb_i.append(j)
        nb_0.append(nb_i)
    return nb_0

def aggregate_introd(nI, nT, coords, f, f_rad):
    """
    Initial outbreak aggregated in foci

    Parameters
    ----------
    nI : Number of infected in each focus.
    nT : Total number of infected.
    coords : coordinates [x,y] of the individuals.
    f : Number of initial foci.
    f_rad : Radius of each initial focus (m).

    Returns
    -------
    Array with values 0, 1 of length equal to the number of individuals. 
    0 = susceptible (not infected) individuals.
    1 = infected individuals (initial outbreak).

    """
    nI = int(nI)
    nT = int(nT)
    f = int(f)
    s_t0 = np.zeros(nT)
    idx = sample(range(0,len(s_t0)), f)
    s_t0[idx] = 1
    indI = idx
    indS = np.where(s_t0==0)[0]
    nb = nb_dist(f_rad, coords, indI, indS)
    nb_I = []
    for i in range(len(nb)):
        idxI = sample(nb[i],(nI-1))
        nb_I.append(idxI)
    nb_I_flat = [j for i in nb_I for j in i]
    s_t0[nb_I_flat]=1
    return(s_t0)


@nb.njit(fastmath=True)
def f_inf(coords, indS_i, indI, s_t, beta_arr, rg, nu, dist_max):
    """
    Calculation of the probability of infection of the susceptible individual i 
    (indS_i) based on the force of infection.

    Parameters
    ----------
    coords : coordinates [x,y] of the individuals.
    indS_i : susceptible individual index i.
    indI : indexes of infected individuals.
    s_t : array with the states of individuals based on the parameter of 
        transmission rate reduction as a function of symptom expression.
        0 = susceptible
        0.015 = asymptomatic infected (from infection to symptom expression)
        [0.03,1] = symptomatic infected (+0.015 each time t from symptom expression)
    beta_arr : array with the transmission rate of infected individuals.
        (an array with the value for each individual is used so that the transmission 
         rate reduction can be applied later.)
    rg : range parameter of the Matérn correlation function.
    nu : smoothness parameter of the Matérn correlation function.
    dist_max : Maximum distance to evaluate spatial correlation.

    Returns
    -------
    dx_p : Probability of infection if individual i.

    """
    dx=0
    for j in indI:
        dist = dist_eu(coords[indS_i], coords[j])
        if dist < dist_max:
            cor_m = cor_matern(dist, rg, nu)
            lmbda = s_t[j]
            dx += (cor_m * beta_arr[j] * lmbda)
    dx_p = 1 - np.exp(-(dx))
    return dx_p


@nb.njit(parallel=True, fastmath=True)
def stat_t(s_t, alm_cells, coords, dist_max, rg, nu, beta_arr, sigma, sint_expr):
    """
    Disease spread for time t

    Parameters
    ----------
    s_t : array with the states of individuals based on the parameter of 
        transmission rate reduction as a function of symptom expression.
    alm_cells : cell in which each individual is located.
    coords : coordinates [x,y] of the individuals.
    dist_max : Maximum distance to evaluate spatial correlation.
    rg : range parameter of the Matérn correlation function.
    nu : smoothness parameter of the Matérn correlation function.
    beta_arr : array with the transmission rate of infected individuals.
        (an array with the value for each individual is used so that the transmission
         rate reduction can be applied later.)
    sigma : Time from infection to symptom expression.
    sint_expr : Counter since infection for symptom expression.

    Returns
    -------
    s_t : array with the states of individuals;
    sint_expr : Counter updated since infection for symptom expression for each
        asymptomatic infected individual.

    """
    indI = np.where(s_t > 0)[0]
    indS = np.where(s_t == 0)[0]
    indIA = np.where((s_t > 0) & (s_t < 1))[0]
    # cell in which each susceptible individual is located
    indS_cell = alm_cells[indS]
    cells = np.unique(indS_cell)
    for cell in nb.prange(len(cells)):
        # Selection of one susceptible individual per cell
        i_cell = np.random.choice(indS[indS_cell == cells[cell]], 1)[0]
        # Probability of infection of each susceptible selected in each cell
        prob = f_inf(coords, i_cell, indI, s_t, beta_arr, rg, nu, dist_max)
        for i in indS[indS_cell == cells[cell]]:
            # Random variable from a Bernoulli distribution with probability 'prob'
            u = np.random.binomial(1, prob)
            # If the Bernoulli random variable (u) is equal to 1, the individual becomes Ia
            s_t[i] = 0.015 if u==1 else 0
    for j in nb.prange(len(indIA)):
        if sint_expr[indIA[j]] < sigma:
            sint_expr[indIA[j]] += 1
        elif s_t[indIA[j]] < 0.99:
            s_t[indIA[j]] += 0.015
        else:
            s_t[indIA[j]] += 0.01
    return(s_t, sint_expr)

#############################################################################

# Dispersión sin intervenciones    
@nb.njit(parallel=True, fastmath=True)
def stat(s_init, alm_cells, coords, beta, rg, nu, sigma, t_max, t_save):
    """
    Disease spread without interventions.

    Parameters
    ----------
    s_init : Initial disease state of each individual (t=0)
    alm_cells : cell in which each individual is located.
    coords : coordinates [x,y] of the individuals.
    beta : transmission rate.
    rg : range parameter of the Matérn correlation function.
    nu : smoothness parameter of the Matérn correlation function.
    sigma : Time from infection to symptom expression.
    t_max : Number of time periods to simulate disease spread (t final).
    t_save : Interval of time periods t to be saved in the output file
        (e.g., t_save = 6, will be saved every 6 t).

    Returns
    -------
    s_final : status of each individual at each time t.

    """
    t_s = np.arange(0, t_max+1, t_save)
    t_s = t_s.astype(nb.int32)
    s_final = np.empty((len(t_s), len(s_init)), dtype=np.float64)
    s_final[0] = s_init
    s_t = np.copy(s_init)
    beta_arr = np.array([beta]*len(s_t))
    dist_max = rg + (rg/2)
    sint_expr = np.array([0]*len(coords))
    sint_expr[s_init>0] = 1
    for t in range(1, t_max+1):
        s_t, sint_expr = stat_t(s_t, alm_cells, coords, dist_max, rg, nu, beta_arr, sigma, sint_expr)
        if t in t_s:
            idx = int(np.where(t_s == t)[0][0])
            s_final[idx] = s_t
    return s_final


##################################################################

@nb.njit(fastmath=True)
def sample_size(n, cl, dp, ms):
    """
    Sample size based on hypergeometric distribution.

    Parameters
    ----------
    n : Population size
    cl : Confidence level
    dp : Design prevalence
    ms : Method sensitivity

    Returns
    -------
    c : Sample size

    """
    a = 1-(1-cl)**(1/(n*dp))
    b = n-(0.5*(n*dp*ms-1))
    c = (a*b)/ms
    return (c)

def bound_delim(centroid, coords, r):
    """
    Delimitation perimeter

    Parameters
    ----------
    centroid : index of individuals to be the center of the delimitation circumference.
    coords : coordinates [x,y] of the individuals.
    r : delimitation radius.

    Returns
    -------
    bound : Polygon with perimeter of delimitation.

    """
    centroid_coords = coords[centroid]
    c = [Point(x, y).buffer(r) for x, y in centroid_coords]
    bound = gpd.GeoSeries(unary_union(c))
    return(bound)

def nb_delim(all_individuals, centroid, coords, r):
    """
    Individuals within the delimitation perimeter
    
    Parameters
    ----------
    all_individuals : index of individuals to be assessed if they are within the delimitation perimeter.
    centroid : index of individuals to be the center of the delimitation circumference.
    coords : coordinates [x,y] of the individuals.
    r : delimitation radius.

    Returns
    -------
    nb_0 : index of individuals within the delimitation perimeter.

    """
    d = cdist(coords[centroid], coords[all_individuals])
    m = np.unique(np.concatenate([np.where(d[i]<=r)[0] for i in range(len(d))]))
    nb_0 = all_individuals[m]
    return nb_0

###########################################

@nb.njit(parallel=True, fastmath=True)
def stat_6t(t_end, s_t, alm_cells, coords, dist_max, rg, nu, beta_arr, sigma, sint_expr):
    """
    Disease spread without intervention during 6 time periods t

    Parameters
    ----------
    t_end : time t in which the dispersion must be finished
    s_t : Disease state of each individual
    alm_cells : cell in which each individual is located.
    coords : coordinates [x,y] of the individuals.
    dist_max : Maximum distance to evaluate spatial correlation.
    rg : range parameter of the Matérn correlation function.
    nu : smoothness parameter of the Matérn correlation function.
    beta_arr : array with the transmission rate of infected individuals.
        (an array with the value for each individual is used so that the transmission
         rate reduction can be applied later.)
    sigma : Time from infection to symptom expression.
    sint_expr : Counter since infection for symptom expression.

    Returns
    -------
    s_t : status of each individual at t = t_end.
    sint_expr: Counter updated since infection for symptom expression for each
        asymptomatic infected individual.
 
    """
    for i in range(t_end-5, t_end+1):
        indI = np.where(s_t > 0)[0]
        indS = np.where(s_t == 0)[0]
        indIA = np.where((s_t > 0) & (s_t < 1))[0]
        # ceda a la que pertenece cada indS
        indS_cell = alm_cells[indS]
        cells = np.unique(indS_cell)
        for cell in nb.prange(len(cells)):
            # Selección de un individuo susceptible por celda
            i_cell = np.random.choice(indS[indS_cell == cells[cell]], 1)[0]
            # Probabilidad de infectarse de cada susceptible seleccionado en cada celda
            prob = f_inf(coords, i_cell, indI, s_t, beta_arr, rg, nu, dist_max)
            for i in indS[indS_cell == cells[cell]]:
                # Random variable from a Bernoulli distribution with probability 'prob'
                u = np.random.binomial(1, prob)
                # If the Bernoulli random variable (u) is equal to 1, the individual becomes Ia
                s_t[i] = 0.015 if u==1 else 0
        for j in nb.prange(len(indIA)):
            if sint_expr[indIA[j]] < sigma:
                sint_expr[indIA[j]] += 1
            elif s_t[indIA[j]] < 0.99:
                s_t[indIA[j]] += 0.015
            else:
                s_t[indIA[j]] += 0.01
    return(s_t, sint_expr)

################################################################################

def survey(ind, params):
    """
    Function for the selection of individuals to be sampled with the one-step 
    surveillance approach.

    Parameters
    ----------
    ind : indices of all individuals in the area to be surveyed
    params : dictionary with parameters for sample size calculation 
        'cl' : conficence level 
        'dp' : design prevalence
        'ms' : method sensitivity

    Returns
    -------
    ind_sample : indexes of selected individuals.

    """
    n_sample = round(sample_size(n=len(ind), cl=params['cl'], 
                                 dp=params['dp'], ms=params['ms']))
    if n_sample <= len(ind):
        ind_sample = np.array(sample(list(ind), n_sample))
    else:
        ind_sample = np.empty(0, int)
    return ind_sample

def strategies(xy, parameters, t_max, additional_control=False):
    """
    Implementation of outbreak management plans with the one-step surveillance approach.

    Parameters
    ----------
    xy : DataFrame of the individuals in the population, one row per individual, with columns:
        x : x coordinates
        y : y coordinates
        cell_ha : cell in which each individual is located.
        t0 : disease status (numeric)
    parameters : nested dictionary with the parameters
        'spread': disease spread parameters
            'beta': transmission rate
            'beta_zt': optional, transmission rate reduction in the buffer zone.
                If specified, aditional_control should be True.
            'rg': range parameter of the Matérn correlation function.
            'nu': smoothness parameter of the Matérn correlation function.
            'sigma': time from infection to symptom expression.
        'survey_d': detection survey parameters
            'cl' : conficence level 
            'dp' : design prevalence
            'ms' : method sensitivity
        'survey_zt': survey in buffer zone parameters
            'cl' : conficence level 
            'dp' : design prevalence
            'ms' : method sensitivity
        'control': control strategies parameters
            'r_zt': radius of delimitation of the buffer zone
            'r_er': eradication radious
    t_max : time to finish the simulation
    additional_control : True or False, optional
        Disease transmission rate reduction in the buffer zone
        The default is False. If True, 'beta_zt' should be specified.

    Returns
    -------
    s_final : result with the status of each individual at each time t.
    s_survey : results of the surveys.

    """
    coords = xy[['x', 'y']].values
    alm_cells = np.array(xy['cell_ha'])
    s_t = np.copy(np.array(xy['t0'])) # copy of the initial status
    s_final = np.copy(np.array(xy['t0']))
    s_survey = np.empty(len(s_t), dtype=np.float64)
    s_survey[:] = np.nan
    sint_expr = np.array([0]*len(coords))
    sint_expr[s_t>0] = 1
    ind_outZT = np.array(range(len(s_t)))
    zt, ze = [], []
    t_eradication = 0
    ind_sample = np.empty(0, int)
    s_er = []
    beta_arr = np.array([parameters['spread']['beta']]*len(s_t))
    rg = parameters['spread']['rg']
    nu, sigma = parameters['spread']['nu'], parameters['spread']['sigma']
    dist_max = rg + (rg/2)
    for t in np.arange(6, (t_max)+6, 6):
        # Spread 6 t
        s_t, sint_expr = stat_6t(t, s_t, alm_cells, coords, dist_max, rg, nu, beta_arr, sigma, sint_expr)
        if t % 12 == 0:
            # Detection survey each 12 t
            if len(ind_outZT)>0:
                ind_sample = survey(ind_outZT, parameters['survey_d'])
            # Annual survey in Bufffer zone (BZ)
            if len(zt) > 0:
                if len(ind_outZT)==0: # If that t there has been no detection survey
                    ind_sample = np.empty(0, int)
                ind_sample_zt = survey(zt, parameters['survey_zt'])
                ind_sample = np.append(ind_sample, ind_sample_zt)
            # Saves the results of the surveys
            res_survey = np.empty(len(s_t), dtype=np.float64)
            res_survey[:] = np.nan
            res_survey[ind_sample] = s_t[ind_sample]
            s_survey = np.vstack([s_survey, res_survey])
        # BZ Delimitation/Update
            if any(s_t[ind_sample]>0):
                positivos_detectados = ind_sample[s_t[ind_sample]>0] 
                zt_i = nb_delim(all_individuals=np.array(ind_outZT), 
                                centroid=positivos_detectados,
                                coords=coords, r=parameters['control']['r_zt'])
                # Added to the previous BZ
                zt.extend(zt_i)
                if additional_control:
                    beta_arr[zt] = parameters['spread']['beta_zt']
                ind_outZT = list(np.setdiff1d(ind_outZT, zt))
                zt = list(np.setdiff1d(zt, s_er))
                # Calculation of those to be removed
                ze = nb_delim(all_individuals=np.array(zt), centroid=positivos_detectados, 
                              coords=coords, r=parameters['control']['r_er'])
                t_eradication = t + 6
        # Eradication
        elif t == t_eradication:
            s_t[ze] = -1
            s_er = list(np.where(s_t == -1)[0]) # Eradicated
            zt = list(np.setdiff1d(zt, s_er))
        # Store the spread result for time t (every 6 t).
        s_final = np.vstack([s_final, s_t])
    return(s_final, s_survey)
    
#################################################################################

def two_steps(cells, alm_cells, params):
    """
    Function for the selection of individuals to be sampled with the two-steps 
    surveillance approach.

    Parameters
    ----------
    cells : cells in the area to be surveyed.
    alm_cells : cell in which each individual is located.
    params : nested dictionary with parameters for sample size calculation 
        'step1' : parameters for step 1 (individuals to be sampled in each cell)
            'dp' : design prevalence
            'ms' : method sensitivity
        'step2' : parameters for step 2 (cells to be inspected)
            'cl' : conficence level 
            'dp' : design prevalence

    Returns
    -------
    ind_sample : indexes of selected individuals.

    """
    mask = np.isin(alm_cells, cells)
    alm_cells_mask = alm_cells[mask]
    n_alms = np.unique(alm_cells_mask, return_counts=True)[1]
    n_mean = round(np.mean(n_alms))
    cl_test = np.around(np.arange(0.05, 1, 0.05), 2)
    ni_test = np.array([round(sample_size(n=n_mean, cl=cl_i, 
                                      dp=params['step1']['dp'], 
                                      ms=params['step1']['ms'])) 
                    for cl_i in cl_test])
    ms_test = np.around(np.arange(0.05, 1, 0.05), 2)
    nh_test = np.array([round(sample_size(n=len(cells), 
                                          cl=params['step2']['cl'], 
                                          dp=params['step2']['dp'],
                                          ms=ms_i))
                        for ms_i in ms_test])
    cl_uplim = max(cl_test[ni_test<=n_mean])
    ms_dlim = min(ms_test[nh_test<=len(cells)])
   
    if cl_uplim < ms_dlim:
        cl1 = cl_uplim
        n_h = len(cells)
    else:
        if cl_uplim <= 0.5:
            cl1 = cl_uplim
        else:
            if ms_dlim <= 0.5:
                cl1 = 0.5
            else:
                cl1 = ms_dlim
        n_h = nh_test[ms_test==cl1][0]
    n_i = ni_test[cl_test==cl1][0]    
    n_total = n_i * n_h
    arr = np.arange(len(cells))
    sample_h= sample(list(arr), n_h) # índices
    while sum(n_alms[sample_h])<n_total:
        arr = np.setdiff1d(arr, sample_h)
        if len(arr)>0:
            sample_h = np.delete(sample_h, np.where(n_alms[sample_h]==min(n_alms[sample_h]))[0][0])
            sample_h = np.append(sample_h, sample(list(arr), 1))
        else:
            n_total = sum(n_alms[sample_h])
    prop_sample = np.round((n_alms[sample_h]/sum(n_alms[sample_h]))*n_total)
    ind_sample = [sample(list(np.where(alm_cells==cells[sample_h[i]])[0]), int(prop_sample[i]))
                for i in range(len(sample_h))]
    ind_sample = np.array(list(itertools.chain(*ind_sample)))
    return(ind_sample)

def two_steps_strategy(xy, parameters, t_max, additional_control=False):
    """
    Implementation of outbreak management plans with the two-steps surveillance approach.

    Parameters
    ----------
    xy : DataFrame of the individuals in the population, one row per individual, with columns:
        x : x coordinates
        y : y coordinates
        cell_ha : cell in which each individual is located.
        t0 : disease status (numeric)
    parameters : nested dictionary with the parameters
        'spread': disease spread parameters
            'beta': transmission rate
            'beta_zt': optional, transmission rate reduction in the buffer zone.
                If specified, aditional_control should be True.
            'rg': range parameter of the Matérn correlation function.
            'nu': smoothness parameter of the Matérn correlation function.
            'sigma': time from infection to symptom expression.
        'survey_d': detection survey parameters
            'step1' : parameters for step 1 (individuals to be sampled in each cell)
                'dp' : design prevalence
                'ms' : method sensitivity
            'step2' : parameters for step 2 (cells to be inspected)
                'cl' : conficence level 
                'dp' : design prevalence
        'survey_zt': survey in buffer zone parameters
            'step1' : parameters for step 1 (individuals to be sampled in each cell)
                'dp' : design prevalence
                'ms' : method sensitivity
            'step2' : parameters for step 2 (cells to be inspected)
                'cl' : conficence level 
                'dp' : design prevalence
        'control': control strategies parameters
            'r_zt': radius of delimitation of the buffer zone
            'r_er': eradication radious
    t_max : time to finish the simulation
    additional_control : True or False, optional
        Disease transmission rate reduction in the buffer zone
        The default is False. If True, 'beta_zt' should be specified.

    Returns
    -------
    s_final : result with the status of each individual at each time t.
    s_survey : results of the surveys.

    """
    coords = xy[['x', 'y']].values
    s_t = np.copy(np.array(xy['t0'])) 
    s_final = np.copy(np.array(xy['t0']))
    s_survey = np.empty(len(s_t), dtype=np.float64)
    s_survey[:] = np.nan
    sint_expr = np.array([0]*len(coords))
    sint_expr[s_t>0] = 1
    zt, ze = [], []
    t_eradication = 0
    s_er = []
    ind_outZT = np.array(range(len(s_t)))
    beta_arr = np.array([parameters['spread']['beta']]*len(s_t))
    rg = parameters['spread']['rg']
    nu, sigma = parameters['spread']['nu'], parameters['spread']['sigma']
    dist_max = rg + (rg/2)
    # cells:
    alm_cells = np.copy(np.array(xy['cell_ha']))
    cells_outZT = np.unique(alm_cells)
    cells_ZT = []
    
    for t in np.arange(6, (t_max)+6, 6):
        # Spread  6t
        s_t, sint_expr = stat_6t(t, s_t, alm_cells, coords, dist_max, rg, nu, beta_arr, sigma, sint_expr)
        if t % 12 == 0:
            # Detection survey each 12 t
            if len(cells_outZT) > 0:
                ind_sample = two_steps(cells_outZT, alm_cells, parameters['survey_d'])
               
            # Annual surveys in BZ
            if len(zt) > 0:
                if len(cells_outZT) == 0: # If that t there has been no detection survey
                    ind_sample = np.empty(0, int)
                ind_sample_zt = two_steps(cells_ZT, alm_cells, parameters['survey_zt'])
                ind_sample = np.append(ind_sample, ind_sample_zt)
            # Save survey result
            res_survey = np.empty(len(s_t), dtype=np.float64)
            res_survey[:] = np.nan
            res_survey[ind_sample] = s_t[ind_sample]
            s_survey = np.vstack([s_survey, res_survey])
            # BZ Delimitation/Update
            if any(s_t[ind_sample]>0):
                positivos_detectados = ind_sample[s_t[ind_sample]>0] 
                zt_i = nb_delim(all_individuals=np.array(ind_outZT), 
                                centroid=positivos_detectados,
                                coords=coords, r=parameters['control']['r_zt'])
                # Added to the previous BZ
                zt.extend(zt_i)
                zt = list(np.unique(zt))
                if additional_control:
                    beta_arr[zt] = parameters['spread']['beta_zt']
                ind_outZT = list(np.setdiff1d(ind_outZT, zt))
                zt = list(np.setdiff1d(zt, s_er))
                # Calculation of those to be removed
                ze = nb_delim(all_individuals=np.array(zt), centroid=positivos_detectados, 
                              coords=coords, r=parameters['control']['r_er'])
                # Cells inside and outside the BZ
                cells_ZT = np.unique(alm_cells[zt])
                cells_outZT = np.setdiff1d(cells_outZT, cells_ZT)
                t_eradication = t + 6
            # Eradication
        elif t == t_eradication:
            s_t[ze] = -1
            s_er = list(np.where(s_t == -1)[0]) # Eradicated
            zt = list(np.setdiff1d(zt, s_er))
            cells_ZT = np.unique(alm_cells[zt])
        # Store the spread result for time t (every 6 t).
        s_final = np.vstack([s_final, s_t])
    return(s_final, s_survey)
    
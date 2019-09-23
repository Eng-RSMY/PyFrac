# -*- coding: utf-8 -*-
"""
This file is part of PyFrac.

Created by Haseeb Zia on Fri Sep 6 16:53:19 2019.
Copyright (c) "ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory", 2016-2019. All rights
reserved. See the LICENSE.TXT file for more details.
"""

from elastohydrodynamic_solver import finiteDiff_operator_laminar
import numpy as np
from math import ceil
# from elasticity import calculate_fluid_flow_characteristics_laminar
import copy

s_max = 1000
a = np.zeros(s_max)
b = np.zeros(s_max)
mu = np.zeros(s_max)
nu = np.zeros(s_max)

b[:2] = 1 / 3
a[:2] = 1 - b[:2]
for j in range(2, s_max):
    b[j] = (j * j + j - 2) / (2 * j * (j + 1))
    a[j] = 1 - b[j]
    mu[j] = (2 * j - 1) * b[j] / (j * b[j - 1])
    nu[j] = - (j - 1) * b[j] / (j * b[j - 2])



def solve_width_pressure_RKL2(Fr_lstTmStp, sim_properties, fluid_properties, mat_properties, EltTip, partlyFilledTip, C,
                         FillFrac, EltCrack, InCrack, LkOff, wTip, timeStep, Qin, perfNode, Vel, corr_ribbon):

    C_EltTip = np.copy(C[np.ix_(EltTip[partlyFilledTip],
                                EltTip[partlyFilledTip])])  # keeping the tip element entries to restore current
    #  tip correction. This is done to avoid copying the full elasticity matrix.

    # filling fraction correction for element in the tip region
    FillF = FillFrac[partlyFilledTip]
    for e in range(0, len(partlyFilledTip)):
        r = FillF[e] - .25
        if r < 0.1:
            r = 0.1
        ac = (1 - r) / r
        C[EltTip[partlyFilledTip[e]], EltTip[partlyFilledTip[e]]] *= (1. + ac * np.pi / 4.)


    dt_CFL = 5 * fluid_properties.viscosity * min(Fr_lstTmStp.mesh.hx, Fr_lstTmStp.mesh.hy) ** 3 / \
             (mat_properties.Eprime * np.max(Fr_lstTmStp.w) ** 3)
    s = ceil(-0.5 + (8 + 16 * timeStep / dt_CFL) ** 0.5 / 2)
    # s = 60
    delt_wTip = wTip - Fr_lstTmStp.w[EltTip]
    tip_delw_step = delt_wTip / s
    n_channel = len(Fr_lstTmStp.EltChannel)

    EltCrack_k = np.concatenate((Fr_lstTmStp.EltChannel, EltTip))
    corr_nei = find_neighbor_in_crack(EltCrack_k, Fr_lstTmStp.mesh)

    cond_0 = finiteDiff_operator_laminar(Fr_lstTmStp.w,
                                         EltCrack_k,
                                         Fr_lstTmStp.muPrime,
                                         Fr_lstTmStp.mesh,
                                         InCrack,
                                         corr_nei)

    # pf_0 = Fr_lstTmStp.pFluid[EltCrack_k]
    pf_0 = np.dot(C[np.ix_(EltCrack_k, EltCrack_k)], Fr_lstTmStp.w[EltCrack_k])
    M_0 = np.dot(cond_0[:n_channel, :-1], pf_0) + Qin[Fr_lstTmStp.EltChannel] / Fr_lstTmStp.mesh.EltArea

    W_0 = Fr_lstTmStp.w[EltCrack_k]
    W_1 = np.zeros(len(EltCrack_k))
    mu_t_1 = 4 / (3 * (s * s + s - 2))
    W_1[:n_channel] = Fr_lstTmStp.w[Fr_lstTmStp.EltChannel] + timeStep * mu_t_1 * M_0
    W_1[n_channel:] = Fr_lstTmStp.w[EltTip] + tip_delw_step

    W_jm1 = np.copy(W_1)
    W_jm2 = np.copy(W_0)
    tau_M0 = timeStep * M_0
    param_pack = (Fr_lstTmStp.muPrime, Fr_lstTmStp.mesh, InCrack, corr_nei)
    sum_mut = copy.copy(mu_t_1)
    sum_gamma = 0

    for j in range(2, s + 1):

        # C_EltTip = np.copy(C[np.ix_(EltTip[partlyFilledTip],
        #                             EltTip[partlyFilledTip])])  # keeping the tip element entries to restore current
        # #  tip correction. This is done to avoid copying the full elasticity matrix.
        #
        # # filling fraction correction for element in the tip region
        # FillF = j / s * FillFrac[partlyFilledTip]
        # for e in range(0, len(partlyFilledTip)):
        #     r = FillF[e] - .25
        #     if r < 0.1:
        #         r = 0.1
        #     ac = (1 - r) / r
        #     C[EltTip[partlyFilledTip[e]], EltTip[partlyFilledTip[e]]] *= (1. + ac * np.pi / 4.)


        W_j, sum_mut, sum_gamma = RKL_substep(j, s, W_jm1, W_jm2, W_0, EltCrack_k, n_channel, tip_delw_step, param_pack, C, timeStep, tau_M0, Qin[Fr_lstTmStp.EltChannel], wTip, sum_mut, sum_gamma)
        W_jm2 = W_jm1
        W_jm1 = W_j

        # C[np.ix_(EltTip[partlyFilledTip], EltTip[partlyFilledTip])] = C_EltTip

    w_s = np.zeros(Fr_lstTmStp.mesh.NumberOfElts)
    pf_s = np.zeros(Fr_lstTmStp.mesh.NumberOfElts)

    w_s[EltCrack_k] = W_j
    pf_s[EltCrack_k] = np.dot(C[np.ix_(EltCrack_k, EltCrack_k)], W_j)

    C[np.ix_(EltTip[partlyFilledTip], EltTip[partlyFilledTip])] = C_EltTip

    return w_s, pf_s, None


def RKL_substep(j, s, W_jm1, W_jm2, W_0, crack, n_channel, tip_delw_step, param_pack, C, tau, tau_M0, Qin, W0_tip, sum_mut, sum_gamma):

    pf = np.dot(C[np.ix_(crack, crack)], W_jm1)
    muPrime, Mesh, InCrack, neiInCrack = param_pack
    w_jm1 = np.zeros(Mesh.NumberOfElts)
    w_jm1[crack] = W_jm1
    cond = finiteDiff_operator_laminar(w_jm1, crack, muPrime, Mesh, InCrack, neiInCrack)
    M_jm1 = np.dot(cond[:n_channel, :-1], pf) + Qin / Mesh.EltArea
    mu_t = 4 * (2 * j - 1) * b[j] / (j * (s * s + s - 2) * b[j - 1])
    gamma_t = -a[j - 1] * mu_t
    W_j = np.zeros(len(crack))
    W_j[:n_channel] = mu[j] * W_jm1[:n_channel] + nu[j] * W_jm2[:n_channel] + (1 - mu[j] - nu[j]) * W_0[:n_channel] + mu_t * tau * M_jm1 + gamma_t * tau_M0
    sum_mut += mu_t
    sum_gamma += gamma_t
    W_j[n_channel:] = W0_tip + j * tip_delw_step
    # W_j[n_channel:] = s * tip_delw_step

    return W_j, sum_mut, sum_gamma


def find_neighbor_in_crack(EltCrack, mesh):
    # The code below finds the indices(in the EltCrack list) of the neighbours of all the cells in the crack.
    # This is done to avoid costly slicing of the large numpy arrays while making the linear system during the
    # fixed point iterations. For neighbors that are outside the fracture, len(EltCrack) + 1 is returned.
    corr_nei = np.full((len(EltCrack), 4), len(EltCrack), dtype=np.int)
    for i, elem in enumerate(EltCrack):
        corresponding = np.where(EltCrack == mesh.NeiElements[elem, 0])[0]
        if len(corresponding) > 0:
            corr_nei[i, 0] = corresponding
        corresponding = np.where(EltCrack == mesh.NeiElements[elem, 1])[0]
        if len(corresponding) > 0:
            corr_nei[i, 1] = corresponding
        corresponding = np.where(EltCrack == mesh.NeiElements[elem, 2])[0]
        if len(corresponding) > 0:
            corr_nei[i, 2] = corresponding
        corresponding = np.where(EltCrack == mesh.NeiElements[elem, 3])[0]
        if len(corresponding) > 0:
            corr_nei[i, 3] = corresponding

    return corr_nei

def tip_edges(EltRibbon, EltTip, mesh):
    corr_tElts = []
    for elt in EltRibbon:
        t_edges = []
        for neighbor in mesh.NeiElements[elt]:
            tip_elt = np.where(neighbor == EltTip)[0]
            if len(tip_elt) > 0:
                t_edges.append(tip_elt[0])
        corr_tElts.append(t_edges)
    return corr_tElts


def solve_width_pressure_RKL2_2(Fr_lstTmStp, sim_properties, fluid_properties, mat_properties, EltTip, partlyFilledTip, C,
                         FillFrac, EltCrack, InCrack, LkOff, wTip, timeStep, Qin, perfNode, Vel, corr_ribbon):

    C_EltTip = np.copy(C[np.ix_(EltTip[partlyFilledTip],
                                EltTip[partlyFilledTip])])  # keeping the tip element entries to restore current
    #  tip correction. This is done to avoid copying the full elasticity matrix.

    # filling fraction correction for element in the tip region
    FillF = FillFrac[partlyFilledTip]
    for e in range(0, len(partlyFilledTip)):
        r = FillF[e] - .25
        if r < 0.1:
            r = 0.1
        ac = (1 - r) / r
        C[EltTip[partlyFilledTip[e]], EltTip[partlyFilledTip[e]]] *= (1. + ac * np.pi / 4.)


    dt_CFL = 5 * fluid_properties.viscosity * min(Fr_lstTmStp.mesh.hx, Fr_lstTmStp.mesh.hy) ** 3 / \
             (mat_properties.Eprime * np.max(Fr_lstTmStp.w) ** 3)
    s = ceil(-0.5 + (8 + 16 * timeStep / dt_CFL) ** 0.5 / 2)
    # s = 599
    delt_wTip = wTip - Fr_lstTmStp.w[EltTip]
    tip_delw_step = delt_wTip / s
    n_channel = len(Fr_lstTmStp.EltChannel)

    EltCrack_k = np.concatenate((Fr_lstTmStp.EltChannel, EltTip))
    corr_nei = find_neighbor_in_crack(EltCrack_k, Fr_lstTmStp.mesh)

    corr_rbn_crk = np.zeros(len(EltTip), dtype=int)
    for i in range(len(EltTip)):
        corr_rbn_crk[i] = np.where(EltCrack_k == corr_ribbon[i])[0]

    cond_0 = finiteDiff_operator_laminar(Fr_lstTmStp.w,
                                         EltCrack_k,
                                         Fr_lstTmStp.muPrime,
                                         Fr_lstTmStp.mesh,
                                         InCrack,
                                         corr_nei)

    # pf_0 = Fr_lstTmStp.pFluid[EltCrack_k]
    pf_0 = np.dot(C[np.ix_(EltCrack_k, EltCrack_k)], Fr_lstTmStp.w[EltCrack_k])
    M_0 = np.dot(cond_0[:, :-1], pf_0) + Qin[EltCrack_k] / Fr_lstTmStp.mesh.EltArea

    W_0 = Fr_lstTmStp.w[EltCrack_k]
    # W_1 = np.zeros(len(EltCrack_k))
    mu_t_1 = 4 / (3 * (s * s + s - 2))
    W_1 = Fr_lstTmStp.w[EltCrack_k] + timeStep * mu_t_1 * M_0
    # W_1[n_channel:] = Fr_lstTmStp.w[EltTip] + tip_delw_step

    W_1_tip = Fr_lstTmStp.w[EltTip] + tip_delw_step
    dwc = W_1[n_channel:] - W_1_tip
    W_1[corr_rbn_crk] += dwc
    W_1[n_channel:] = W_1_tip

    W_jm1 = np.copy(W_1)
    W_jm2 = np.copy(W_0)
    tau_M0 = timeStep * M_0
    param_pack = (Fr_lstTmStp.muPrime, Fr_lstTmStp.mesh, InCrack, corr_nei)
    sum_mut = copy.copy(mu_t_1)
    sum_gamma = 0

    for j in range(2, s + 1):

        # C_EltTip = np.copy(C[np.ix_(EltTip[partlyFilledTip],
        #                             EltTip[partlyFilledTip])])  # keeping the tip element entries to restore current
        # #  tip correction. This is done to avoid copying the full elasticity matrix.
        #
        # # filling fraction correction for element in the tip region
        # FillF = j / s * FillFrac[partlyFilledTip]
        # for e in range(0, len(partlyFilledTip)):
        #     r = FillF[e] - .25
        #     if r < 0.1:
        #         r = 0.1
        #     ac = (1 - r) / r
        #     C[EltTip[partlyFilledTip[e]], EltTip[partlyFilledTip[e]]] *= (1. + ac * np.pi / 4.)


        W_j, sum_mut, sum_gamma = RKL_substep_2(j, s, W_jm1, W_jm2, W_0, EltCrack_k, n_channel, tip_delw_step, param_pack, C, timeStep, tau_M0, Qin, wTip, sum_mut, sum_gamma, corr_rbn_crk)
        W_jm2 = W_jm1
        W_jm1 = W_j

        # C[np.ix_(EltTip[partlyFilledTip], EltTip[partlyFilledTip])] = C_EltTip

    w_s = np.zeros(Fr_lstTmStp.mesh.NumberOfElts)
    pf_s = np.zeros(Fr_lstTmStp.mesh.NumberOfElts)

    w_s[EltCrack_k] = W_j
    pf_s[EltCrack_k] = np.dot(C[np.ix_(EltCrack_k, EltCrack_k)], W_j)

    C[np.ix_(EltTip[partlyFilledTip], EltTip[partlyFilledTip])] = C_EltTip

    return w_s, pf_s, None


def RKL_substep_2(j, s, W_jm1, W_jm2, W_0, crack, n_channel, tip_delw_step, param_pack, C, tau, tau_M0, Qin, W0_tip, sum_mut, sum_gamma, corr_rbn_crk):

    pf = np.dot(C[np.ix_(crack, crack)], W_jm1)
    muPrime, Mesh, InCrack, neiInCrack = param_pack
    w_jm1 = np.zeros(Mesh.NumberOfElts)
    w_jm1[crack] = W_jm1
    cond = finiteDiff_operator_laminar(w_jm1, crack, muPrime, Mesh, InCrack, neiInCrack)
    M_jm1 = np.dot(cond[:, :-1], pf) + Qin[crack] / Mesh.EltArea
    mu_t = 4 * (2 * j - 1) * b[j] / (j * (s * s + s - 2) * b[j - 1])
    gamma_t = -a[j - 1] * mu_t
    W_j = mu[j] * W_jm1 + nu[j] * W_jm2 + (1 - mu[j] - nu[j]) * W_0 + mu_t * tau * M_jm1 + gamma_t * tau_M0
    # sum_mut += mu_t
    # sum_gamma += gamma_t
    W_j_t = W0_tip + j * tip_delw_step

    dwc = W_j[n_channel:] - W_j_t
    W_j[corr_rbn_crk] += dwc
    W_j[n_channel:] = W_j_t

    # W_j[n_channel:] = s * tip_delw_step

    return W_j, sum_mut, sum_gamma

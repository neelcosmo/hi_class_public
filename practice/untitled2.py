#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 10:03:31 2023

@author: user
"""

import numpy as np
import matplotlib.pyplot as plt
from classy import Class
plt.rcParams["figure.dpi"] = 200

cosmo = Class()

params = {
    "lensing":"yes",
    #"non linear":"halofit",
    "output":"tCl, pCl, lCl, mPk, dTk, vTk",
    "k_output_values":"0.0001, 0.01, 0.1",
    "z_pk":"0., 1., 10",
    "input_verbose":"1",
    "background_verbose":"1",
    "thermodynamics_verbose":"1",
    "perturbations_verbose":"1",
    "transfer_verbose":"1",
    "primordial_verbose":"1",
    "spectra_verbose":"1",
    "nonlinear_verbose":"1",
    "lensing_verbose":"1",
    "output_verbose":"1",
    "write warnings":"y"}

m_B = 0.5

params_smg = {
    "Omega_Lambda":"0",
    "Omega_fld":"0",
    "Omega_smg":"-1",
    "output_background_smg":"10",
    "gravity_model":"propto_hubble_n",
    "parameters_smg":f"1., 0.5, 0., 0., 1., 1., {m_B}, 1., 1.",
    "expansion_model":"lcdm",
    "expansion_smg":"0.5",
    # "cs2_safe_smg":"1e-10",
    "kineticity_safe_smg":"1e-4",
    } #default is 0, otherwise the default non-zero value is 1e-4 in hi_class.ini.

cosmo.set(params)
cosmo.set(params_smg)
cosmo.compute()
# errorline
alphaBarr = cosmo.get_background()["braiding_smg"]
# practicearr = cosmo.get_background()["practice_braidingplusone_smg"]
alphaMarr = cosmo.get_background()["Mpl_running_smg"]
alphaTarr = cosmo.get_background()["tensor_excess_smg"]
alphaKarr = cosmo.get_background()["kineticity_smg"]

#print(practicearr - alphaBarr)

M2arr = cosmo.get_background()["M*^2_smg"]
cs2arr = cosmo.get_background()["c_s^2"]
Darr = cosmo.get_background()["kin (D)"]

muarr = cosmo.get_background()["mu_smg"]
Sigmaarr = cosmo.get_background()["Sigma_smg"]

beta2_test = alphaBarr*(1+alphaTarr) + 2*(alphaMarr-alphaTarr)
beta4_test = alphaBarr*(alphaTarr-alphaMarr) - 0.5*alphaBarr*alphaBarr*(1+alphaTarr)
beta1_test = Darr*cs2arr - beta2_test - beta4_test
beta3_test = beta1_test*(1+alphaTarr) + beta2_test*(1+alphaMarr)

mu_test = 2*beta3_test/(M2arr*(2*beta1_test + beta2_test*(2-alphaBarr)))
Sigma_test = (beta1_test+beta2_test+beta3_test)/(M2arr*(2*beta1_test + beta2_test*(2-alphaBarr)))

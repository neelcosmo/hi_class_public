#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  5 16:49:21 2023

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
    "temperature contributions":"lisw",
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

# params_smg = {
#     "Omega_Lambda":"0",
#     "Omega_fld":"0",
#     "Omega_smg":"-1",
#     "gravity_model":"propto_omega",
#     "parameters_smg":"1., 0., 0., 0., 1.",
#     "expansion_model":"lcdm",
#     "expansion_smg":"0.5"}

m_B = 0.25

params_smg = {
    "Omega_Lambda":"0",
    "Omega_fld":"0",
    "Omega_smg":"-1",
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

lnzarr = np.log(cosmo.get_background()["z"])
lnHarr = np.log(cosmo.get_background()["H [1/Mpc]"])
lnalphaBarr = np.log(cosmo.get_background()["braiding_smg"])
Hgradient = np.gradient(lnHarr, lnzarr)
alphaBgradient = np.gradient(lnalphaBarr, lnzarr)
inferredm_B = -4*(Hgradient/alphaBgradient)
diff = inferredm_B - m_B

fig, ax = plt.subplots()
ax.set_xlabel("z")
ax.set_ylabel("loglog slope")
fig.suptitle(f"$m_B =${m_B}")
ax.set_xscale("log")
ax.plot(lnzarr, Hgradient, label="H")
ax.plot(lnzarr, alphaBgradient, label="$\\alpha_B$")
# ax.plot(lnzarr, diff, label=("inferred $m_B$ - input $m_B$"))
ax.legend()

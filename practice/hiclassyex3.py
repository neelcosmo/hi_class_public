#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 16:17:39 2023

@author: user
"""

import numpy as np
import matplotlib.pyplot as plt
from classy import Class
plt.rcParams["figure.dpi"] = 200

cosmo = Class()

params = {
    "lensing":"yes",
    "non linear":"halofit",
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

params_mnu = {
    "N_ur":"2.0328",
    "N_ncdm":"1",
    "m_ncdm":"0.06"}

mnu_arr = [0, 0.06, 0.1, 1]

fig, ax = plt.subplots()
ax.set_xlabel("z")
ax.set_ylabel("$\\rho_{\\nu}$")

fig2, ax2 = plt.subplots()
ax2.set_xlabel("$\ell$")
ax2.set_ylabel("$[\ell(\ell+1)/2\pi]C_\ell^{TT}$")
ax2.set_xscale("log")

karr = np.logspace(-4, 0.75)

fig3, ax3 = plt.subplots()
ax3.set_xlabel("k [h/Mpc]")
ax3.set_ylabel("$P_{nl}(k)$")

# full_pkarr = []

for mnu_ind in range(len(mnu_arr)):
    mnu = mnu_arr[mnu_ind]
    
    cosmo.set(params)
    cosmo.set(params_mnu)
    cosmo.set({"m_ncdm":f"{mnu}"})
    cosmo.compute()
    h = cosmo.get_current_derived_parameters(["h"])[("h")]
    
    zarr = cosmo.get_background()["z"]
    rhonuarr = cosmo.get_background()["(.)rho_ncdm[0]"]
    ax.loglog(zarr, rhonuarr, label="$m_{\\nu}=%1.2f eV" % mnu)
    
    ell = cosmo.lensed_cl()["ell"][2:]
    cl = cosmo.lensed_cl()["tt"][2:]
    dl = cl*ell*(ell+1)/(2*np.pi)
    ax2.plot(ell, dl, label="m_{\\nu}=%1.2f eV" % mnu)
    
    pkarr = []
    # print(cosmo.get_pk_and_k_and_z(nonlinear=False))
    for k_ind in range(len(karr)):
        pkarr.append(cosmo.pk(karr[k_ind]*h, 0)*(h**3))
        # full_pkarr.append(cosmo.pk_lin(karr[k_ind]*h, 0)*(h**3))
        # pkarr.append(cosmo.pk_lin(karr[k_ind], 0))
        # full_pkarr.append(cosmo.pk_lin(karr[k_ind], 0))
    pkarr = np.array(pkarr)
    # print(pkarr)

    ax3.loglog(karr, pkarr, label="$m_{\\nu}$=%1.2f eV" % mnu)
    
    cosmo.struct_cleanup()
    cosmo.empty()

ax.legend()
ax2.legend()
ax3.legend()

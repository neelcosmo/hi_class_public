#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 16:08:21 2023

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

params_smg = {
    "Omega_Lambda":"0",
    "Omega_fld":"0",
    "Omega_smg":"-1",
    "gravity_model":"propto_omega",
    "parameters_smg":"1., 0., 0., 0., 1.",
    "expansion_model":"lcdm",
    "expansion_smg":"0.5",
    "kineticity_safe_smg":"1e-4"} #default is 0, otherwise the default non-zero value is 1e-4 in hi_class.ini.

cBarr = [0, 0.625, 1.25, 1.875, 2.5]

cosmo.set(params)
# cosmo.set(params_smg)

cosmo.compute()

cl_ref = cosmo.lensed_cl()["tt"][2:2501]
ell = cosmo.lensed_cl()["ell"][2:2501]
dl_ref = cl_ref*ell*(ell+1)/(2*np.pi)

fig, axs = plt.subplots(2, 1, sharex = True, figsize = (8.0, 7.0))
fig.suptitle("Only lisw, lensed, kineticity_safe=10$^{-4}$")
axs[1].set_xlabel("$\ell$")
axs[0].set_ylabel("$[\ell(\ell+1)/2\pi]C_\ell^{TT}$")
axs[1].set_ylabel("rel. dev. (%)")
axs[1].set_xscale("log")
axs[0].plot(ell, dl_ref, "k--", linewidth=3, label="$\Lambda$CDM")
axs[1].plot(ell, np.zeros(dl_ref.shape), "k--", linewidth=3)

h = cosmo.get_current_derived_parameters(["h"])["h"]
karr = np.logspace(np.log10(1e-5/h), np.log10(1/h))
pk_ref = []
for k_ind in range(len(karr)):
    k = karr[k_ind]
    pk_ref.append(cosmo.pk(k*h, 0)*(h**3))
pk_ref = np.array(pk_ref)

fig2, axs2 = plt.subplots(2, 1, sharex=True, figsize=(8.0, 7.0))
fig2.suptitle("Only lisw")
axs2[1].set_xlabel("k [h/Mpc]")
axs2[0].set_ylabel("$P(k)$ [Mpc/h]^3")
axs2[1].set_ylabel("rel. dev. (%)")
axs2[0].loglog(karr, pk_ref, "k--", linewidth=3, label="$\Lambda$CDM")
axs2[1].set_xscale("log")
axs2[1].plot(karr, np.zeros(pk_ref.shape), "k--", linewidth=3)

phiphi_ref = cosmo.raw_cl()["pp"][2:2501] #same as that from raw_cl
kappakappa_ref = 0.25*phiphi_ref*(ell**4)
dlkk_ref = kappakappa_ref*ell*(ell+1)/(2*np.pi)

fig3, axs3 = plt.subplots(2, 1, sharex=True, figsize=(8, 7))
axs3[1].set_xscale("log")
axs3[1].set_xlabel("$\ell$")
axs3[0].set_ylabel("$(\ell^4/4)[\ell(\ell+1)/2\pi]C_\ell^{\phi\phi}$")
axs3[0].loglog(ell, dlkk_ref, "k--", linewidth=3, label="$\Lambda$CDM")
axs3[1].set_ylabel("rel dev. (%)")
axs3[1].plot(ell, np.zeros(dlkk_ref.shape), "k--", linewidth=3)


for cB_ind in range(len(cBarr)):
    cosmo.struct_cleanup()
    cosmo.empty()
    cosmo.set(params)
    cosmo.set(params_smg)
    cosmo.set({"parameters_smg":f"1., {cBarr[cB_ind]}, 0., 0., 1."})
    
    cosmo.compute()
    cl = cosmo.lensed_cl()["tt"][2:2501]
    dl = cl*ell*(ell+1)/(2*np.pi)
    
    axs[0].plot(ell, dl, label=f"$c_B$ = {cBarr[cB_ind]}")
    axs[1].plot(ell, 100*(dl-dl_ref)/dl_ref)
    
    pk = []
    for k_ind in range(len(karr)):
        k = karr[k_ind]
        pk.append(cosmo.pk(k*h, 0)*(h**3))
    pk = np.array(pk)
    axs2[0].loglog(karr, pk, label=f"$c_B$ = {cBarr[cB_ind]}")  
    axs2[1].plot(karr, 100*(pk-pk_ref)/pk)
    
    phiphi = cosmo.raw_cl()["pp"][2:2501]
    kappakappa = 0.25*phiphi*(ell**4)
    dlkk = kappakappa*ell*(ell+1)/(2*np.pi)
    
    axs3[0].plot(ell, kappakappa*ell*(ell+1)/(2*np.pi), label=f"$c_B$ = {cBarr[cB_ind]}")
    axs3[1].plot(ell, 100*(dlkk-dlkk_ref)/dlkk_ref, alpha=1)
    
axs[0].legend()
axs2[0].legend()
axs3[0].legend()

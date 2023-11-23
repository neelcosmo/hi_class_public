#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 22:13:13 2023

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

params_smg = {
    "Omega_Lambda":"0",
    "Omega_fld":"0",
    "Omega_smg":"-1",
    "gravity_model":"propto_omega",
    "parameters_smg":"1., 0.625, 0., 0., 1.",
    "expansion_model":"lcdm",
    "expansion_smg":"0.5"}

cosmo.set(params)
cosmo.set(params_smg)
cosmo.compute()

cl_ref = cosmo.lensed_cl()["tt"][2:]
ell = cosmo.lensed_cl()["ell"][2:]
dl_ref = cl_ref*ell*(ell+1)/(2*np.pi)

fig, axs = plt.subplots(2, 1, sharex=True)
fig.suptitle("$c_B = 0.625$")
axs[0].set_xscale("log")
axs[1].set_xlabel("$\ell$")
axs[1].set_ylabel("rel. dev. (%)")
axs[0].set_ylabel("[$\ell(\ell+1)/2\pi]C_\ell^{TT}$")
axs[0].plot(ell, dl_ref, label="$\propto \Omega$ (ref.)")
axs[1].plot(ell, np.zeros(dl_ref.shape))

# fig2, ax2 = plt.subplots()
# fig2.suptitle("$c_B = 0.625$")
# ax2.set_xscale("log")
# ax2.set_xlabel("$\ell$")
# ax2.set_ylabel("rel. diff. (%)")
# ax2.plot(ell, np.zeros(dl_ref.shape), label="$\propto \Omega$ (ref.)")

karr = karr = np.logspace(np.log10(1e-4), np.log10(4))
h = cosmo.get_current_derived_parameters(["h"])["h"]
pk_ref = []
for k_ind in range(len(karr)):
    pk_ref.append(cosmo.pk(karr[k_ind]*h, 0)*h**3)
pk_ref = np.array(pk_ref)

fig2, axs2 = plt.subplots(2, 1, sharex=True)
fig2.suptitle("$c_B = 0.625$")
axs2[1].set_xlabel("k [h/Mpc]")
axs2[1].set_ylabel("rel. dev. (%)")
axs2[0].set_ylabel("P(k) [Mpc/h]$^3$")
axs2[0].loglog(karr, pk_ref, label="$\propto \Omega$ (ref.)")
axs2[1].set_xscale("log")
axs2[1].plot(karr, np.zeros(pk_ref.shape))

n_arr = [0.5, 1, 2, 3, 4]

for n_ind in range(len(n_arr)):
    n = n_arr[n_ind]
    cosmo.struct_cleanup()
    cosmo.empty()
    cosmo.set(params)
    cosmo.set(params_smg)
    cosmo.set({"gravity_model":"propto_omega_n"})
    cosmo.set({"parameters_smg":f"1, .625, 0, 0, 1, {n}"})
    cosmo.compute()
    
    cl = cosmo.lensed_cl()["tt"][2:]
    dl = cl*ell*(ell+1)/(2*np.pi)   
    axs[0].plot(ell, dl, label="$\propto \Omega^{%1.1f}$" % n)
    axs[1].plot(ell, 100*(dl-dl_ref)/dl_ref, label="$\propto \Omega^{%1.1f}$" % n)
    
    pk = []
    for k_ind in range(len(karr)):
        pk.append(cosmo.pk(karr[k_ind]*h, 0)*h**3)
    pk = np.array(pk)
    
    axs2[0].loglog(karr, pk, label="$\propto \Omega^{%1.1f}$" % n)
    axs2[1].plot(karr, 100*(pk-pk_ref)/pk_ref, label="$\propto \Omega^{%1.1f}$" % n)

axs[0].legend()
axs2[0].legend()

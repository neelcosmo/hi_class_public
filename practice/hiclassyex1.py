#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 16:44:43 2023

@author: user
"""

import numpy as np
import matplotlib.pyplot as plt
from classy import Class
plt.rcParams["figure.dpi"] = 200

cosmo = Class()

#The write options don't work because there isn't permission to write on files
# params = {
#     "lensing":"yes",
#     "non linear":"halofit",
#     "output":"tCl, pCl, lCl, mPk, dTk, vTk",
#     "write background":"yes",
#     "write thermodynamics":"yes",
#     "write primordial":"y",
#     "k_output_values":"0.0001, 0.01, 0.1",
#     "z_pk":"0., 1., 10",
#     "input_verbose":"1",
#     "background_verbose":"1",
#     "thermodynamics_verbose":"1",
#     "perturbations_verbose":"1",
#     "transfer_verbose":"1",
#     "primordial_verbose":"1",
#     "spectra_verbose":"1",
#     "nonlinear_verbose":"1",
#     "lensing_verbose":"1",
#     "output_verbose":"1",
#     "write warnings":"y",
#     "write parameters":"yup",
#     "alternative facts":"fake"}

#Fake input only gives a warning on the terminal, but produces an error in
#python and prevents it from running.
# params = {
#     "lensing":"yes",
#     "non linear":"halofit",
#     "output":"tCl, pCl, lCl, mPk, dTk, vTk",
#     "k_output_values":"0.0001, 0.01, 0.1",
#     "z_pk":"0., 1., 10",
#     "input_verbose":"1",
#     "background_verbose":"1",
#     "thermodynamics_verbose":"1",
#     "perturbations_verbose":"1",
#     "transfer_verbose":"1",
#     "primordial_verbose":"1",
#     "spectra_verbose":"1",
#     "nonlinear_verbose":"1",
#     "lensing_verbose":"1",
#     "output_verbose":"1",
#     "write warnings":"y",
#     "alternative facts":"fake"}

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

cosmo.set(params)

cosmo.compute()

b = cosmo.get_background()

fig, ax = plt.subplots()
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("z")
ax.set_ylabel("$\\rho_\gamma$")
ax.plot(b["z"], b["(.)rho_g"])

fig, ax = plt.subplots()
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("z")
ax.set_ylabel("$\\rho_b$")
ax.plot(b["z"], b["(.)rho_b"])

fig, ax = plt.subplots()
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("z")
ax.set_ylabel("$\\rho_{ur}$")
ax.plot(b["z"], b["(.)rho_ur"])

fig, ax = plt.subplots()
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("z")
ax.set_ylabel("$\\rho_c$")
ax.plot(b["z"], b["(.)rho_cdm"])

fig, ax = plt.subplots()
ax.set_xscale("log")
ax.set_xlabel("z")
ax.set_ylabel("$\Omega_\gamma$")
ax.plot(b["z"], b["(.)rho_g"]/b["(.)rho_crit"]) #g is photons, neutrinos are in rho_ur

fig, ax = plt.subplots()
ax.set_xscale("log")
ax.set_xlabel("z")
ax.set_ylabel("$\Omega_b$")
ax.plot(b["z"], b["(.)rho_b"]/b["(.)rho_crit"])

fig, ax = plt.subplots()
ax.set_xscale("log")
ax.set_xlabel("z")
ax.set_ylabel("$\Omega_ur$")
ax.plot(b["z"], b["(.)rho_ur"]/b["(.)rho_crit"])

fig, ax = plt.subplots()
ax.set_xscale("log")
ax.set_xlabel("z")
ax.set_ylabel("$\Omega_c$")
ax.plot(b["z"], b["(.)rho_cdm"]/b["(.)rho_crit"])

fig, ax = plt.subplots()
ax.set_xscale("log")
ax.set_xlabel("z")
ax.set_ylabel("$\Sigma\Omega_i$")
ax.plot(b["z"], (b["(.)rho_g"]+b["(.)rho_b"]+b["(.)rho_ur"]+b["(.)rho_cdm"]+b["(.)rho_lambda"])/b["(.)rho_crit"])

fig, ax = plt.subplots()
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("z")
ax.set_ylabel("$D_L$")
ax.plot(b["z"], b["lum. dist."])

fig, ax = plt.subplots()
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("z")
ax.set_ylabel("$D_A$")
ax.plot(b["z"], b["ang.diam.dist."])

th = cosmo.get_thermodynamics()

tt = cosmo.raw_cl()["tt"][2:]
te = cosmo.raw_cl()["te"][2:]
ee = cosmo.raw_cl()["ee"][2:]
ell = cosmo.raw_cl()["ell"][2:]

fig, ax = plt.subplots()
ax.set_xscale("log")
ax.set_xlabel("$\ell$")
ax.set_ylabel("$[\ell(\ell+1)/2\pi]C_\ell^{TT}$")
ax.plot(ell, tt*ell*(ell+1)/(2*np.pi))

fig, ax = plt.subplots()
ax.set_xscale("log")
ax.set_xlabel("$\ell$")
ax.set_ylabel("$[\ell(\ell+1)/2\pi]C_\ell^{TE}$")
ax.plot(ell, te*ell*(ell+1)/(2*np.pi))

fig, ax = plt.subplots()
ax.set_xscale("log")
ax.set_xlabel("$\ell$")
ax.set_ylabel("$[\ell(\ell+1)/2\pi]C_\ell^{EE}$")
ax.plot(ell, ee*ell*(ell+1)/(2*np.pi))

ttl = cosmo.lensed_cl()["tt"][2:]
elll = cosmo.lensed_cl()["ell"][2:]
tt25 = cosmo.raw_cl(lmax=2500)["tt"][2:]

fig, ax = plt.subplots()
ax.set_xscale("log")
ax.set_xlabel("$\ell$")
ax.set_ylabel("$[\ell(\ell+1)/2\pi]C_\ell^{TT}$")
ax.plot(elll, tt25*elll*(elll+1)/(2*np.pi), label="unlensed")
ax.plot(elll, ttl*elll*(elll+1)/(2*np.pi), label="lensed")
ax.legend()

tt25Dl = tt25*elll*(elll+1)/(2*np.pi)
ttlDl = ttl*elll*(elll+1)/(2*np.pi)

fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True)
ax2.set_xscale("log")
ax2.set_xlabel("$\ell$")
ax1.plot(elll, tt25Dl, label="unlensed")
ax1.plot(elll, ttlDl, label="lensed")
ax1.set_ylabel("$[\ell(\ell+1)/2\pi]C_\ell^{TT}$")
ax1.legend()
ax2.plot(elll, (ttlDl-tt25Dl)/tt25Dl, "g")
ax2.set_ylabel("Residuals")

phiphi = cosmo.raw_cl()["pp"][2:]
kappakappa = 0.25*phiphi*(ell**4)

fig, ax = plt.subplots()
ax.set_yscale("log")
ax.set_xscale("log")
ax.set_xlabel("$\ell$")
ax.set_ylabel("$(\ell^4/4)C_\ell^{\phi\phi}$")
ax.plot(ell, kappakappa*ell*(ell+1)/(2*np.pi))

# phiphi = cosmo.lensed_cl()["pp"][2:] #same as the phiphi from raw_cl
# kappakappa = 0.25*phiphi*(elll**4)

# fig, ax = plt.subplots()
# ax.set_yscale("log")
# ax.set_xscale("log")
# ax.set_xlabel("$\ell$")
# ax.set_ylabel("$(\ell^4/4)C_\ell^{\phi\phi}$")
# ax.plot(elll, kappakappa*elll*(elll+1)/(2*np.pi))

#cosmo.get_pk_and_k_and_z(nonlinear=False)

karr = np.logspace(np.log10(1e-4), np.log10(5), 50)
zarr = np.array([0, 1, 10])
pkarr = np.zeros((len(karr), len(zarr)))
h = cosmo.get_current_derived_parameters(["h"])["h"]

for z_ind in range(len(zarr)):
    for k_ind in range(len(karr)):
        pkarr[k_ind][z_ind] = cosmo.pk_lin(karr[k_ind]*h, zarr[z_ind])*(h**3)
pkarr = pkarr.transpose()

fig, ax = plt.subplots()
ax.set_xlabel("k (h/Mpc)")
ax.set_ylabel("$P_{lin}(k)$ (Mpc/h)$^3$")
for z_ind in range(len(zarr)):
    ax.loglog(karr, pkarr[z_ind], label=f"z={zarr[z_ind]}")
ax.legend()

pkarrnl = np.zeros((len(karr), len(zarr)))

for z_ind in range(len(zarr)):
    for k_ind in range(len(karr)):
        pkarrnl[k_ind][z_ind] = cosmo.pk(karr[k_ind]*h, zarr[z_ind])*(h**3)
pkarrnl = pkarrnl.transpose()

fig, axs = plt.subplots(1, 3, sharey=True)
axs[0].set_ylabel("$P_{nl}(k)$ (Mpc/h)$^3$")
for z_ind in range(len(zarr)):
    axs[z_ind].set_xlabel("k (h/Mpc)")
    axs[z_ind].set_title(f"z={zarr[z_ind]}")
    axs[z_ind].loglog(karr, pkarr[z_ind], label="linear")
    axs[z_ind].loglog(karr, pkarrnl[z_ind], label="nonlinear")
    axs[z_ind].legend()

karrt = cosmo.get_transfer()["k (h/Mpc)"]
transb0 = cosmo.get_transfer()["d_b"]
transb1 = cosmo.get_transfer(1)["d_b"]
transb10 = cosmo.get_transfer(10)["d_b"]

fig, ax = plt.subplots()
ax.loglog(karrt, np.abs(transb0), label="z=0")
ax.loglog(karrt, np.abs(transb1), label="z=1")
ax.loglog(karrt, np.abs(transb10), label="z=10")
ax.set_xlabel("k (h/Mpc)")
ax.set_ylabel("|Baryon density transfer function|")
ax.legend()

tauarr1 = cosmo.get_perturbations()["scalar"][0]["tau [Mpc]"]
tauarr2 = cosmo.get_perturbations()["scalar"][1]["tau [Mpc]"]
tauarr3 = cosmo.get_perturbations()["scalar"][2]["tau [Mpc]"]
deltabk1 = cosmo.get_perturbations()["scalar"][0]["delta_b"]
deltabk2 = cosmo.get_perturbations()["scalar"][1]["delta_b"]
deltabk3 = cosmo.get_perturbations()["scalar"][2]["delta_b"]

fig, ax = plt.subplots()
ax.set_xlabel("Conformal time (Mpc)")
ax.set_ylabel("$|\delta_k|$")
ax.loglog(tauarr1, np.abs(deltabk1), label="k=0.0001")
ax.loglog(tauarr2, np.abs(deltabk2), label="k=0.01")
ax.loglog(tauarr3, np.abs(deltabk3), label="k=0.1")
ax.legend()

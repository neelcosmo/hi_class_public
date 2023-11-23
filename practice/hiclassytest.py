#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 15:32:28 2023

@author: user
"""

from classy import Class

params = {
    "output":"tCl, lCl",
    "l_max_scalars":"2000",
    "lensing":"yes",
    "A_s":"2.3e-9",
    "n_s":"0.9624",
    "h":"0.6711",
    "omega_b":"0.022068",
    "omega_cdm":"0.12029"}

cosmo = Class()

#cosmo.set_default()

cosmo.set(params)

cosmo.compute()

clas = cosmo.lensed_cl(2000)

print(clas)

cosmo.struct_cleanup()

cosmo.empty()
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 08:42:07 2023

@author: user
"""

import numpy as np
import matplotlib.pyplot as plt
from classy import Class
plt.rcParams["figure.dpi"] = 200

data = np.loadtxt("/home/user/hi_class/output/hi_class_practice_background.dat").transpose()

alphaBarr = data[20]
practicearr = data[25]

print(practicearr - alphaBarr)
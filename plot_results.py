# -*- coding: utf-8 -*-
"""
Program to plot PPXF results.

"""

from __future__ import print_function, division

import os 

import glob
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, hstack

from ppxf.ppxf import ppxf
import ppxf.ppxf_util as util


import context
from muse_resolution import get_muse_fwhm, broad2res

import pickle

image = np.zeros((85, 85))

for xpix in range(1, 85):
    for ypix in range(1, 85):

        filename = os.path.join(context.plots_dir, "ppxf_values/ssps_ppxf_x{}_y{}.fits".format(xpix, ypix))

        hdu_list = fits.open(filename, memmap=True)
	
        evt_data = Table(hdu_list[1].data)

        print(evt_data['col0'][0])

        print(xpix, ypix)

        image[ypix][xpix] = evt_data['col0'][0]


vmin, vmax = 9999999999, 0

for i in range(2, len(image[0])):
	for j in range(2, len(image[0])):
		vmin = min(vmin, image[i][j])
		vmax = max(vmax, image[i][j])


print(vmin)
print(vmax)

plt.imshow(image, cmap='gist_earth', vmin=vmin, vmax=vmax, origin="bottom")

plt.colorbar()
plt.show()
        

        
        

# -*- coding: utf-8 -*-
"""
Program to plot PPXF results.

"""

from __future__ import print_function, division

import os 

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

import context

image_vel = np.zeros((85, 85))
image_sigma = np.zeros((85, 85))
image_h3 = np.zeros((85, 85))
image_h4 = np.zeros((85, 85))

#b = []

for xpix in range(1, 85):
    for ypix in range(1, 85):

        filename = os.path.join(context.data_dir, "ppxf_values/ssps_ppxf_x{}_y{"
                                               "}.fits".format(xpix, ypix))

        hdu_list = fits.open(filename, memmap=True)
	
        evt_data = Table(hdu_list[1].data)


        #image_vel[ypix][xpix] = evt_data['col0'][0]

        #image_sigma[ypix][xpix] = evt_data['col1'][0]

        #image_h3[ypix][xpix] = evt_data['col2'][0]

        image_h4[ypix][xpix] = evt_data['col3'][0]


nbins = 7225

hist_vel = image_vel.flatten()
hist_sigma = image_sigma.flatten()
hist_h3 = image_h3.flatten()
hist_h4 = image_h4.flatten()


#plt.imshow(image_vel, cmap='Spectral', vmin=10148, vmax=10153, origin="bottom")
#plt.colorbar()
#plt.show()

#bordas_vel = (9950, 10250)

#histogram_vel = plt.hist(hist_vel, nbins, bordas_vel)

#plt.show()       


#plt.imshow(image_sigma, cmap='Spectral', vmin=240, vmax=270, origin="bottom")
#plt.colorbar()
#plt.show()


#bordas_sigma = (200, 400)

#histogram_sigma = plt.hist(hist_sigma, nbins, bordas_sigma)

#plt.show()


#plt.imshow(image_h3, cmap='Spectral', vmin=-0.01, vmax=0.01, origin="bottom")
#plt.colorbar()
#plt.show()

#bordas_h3 = (-0.3, 0.3)

#histogram_h3 = plt.hist(hist_h3, nbins, bordas_h3)

#plt.show()


plt.imshow(image_h4, cmap='Spectral', vmin=0, vmax=0.02, origin="bottom")
plt.colorbar()
plt.show()

bordas_h4 = (-1, 1)

histogram_h4 = plt.hist(hist_h4, nbins, bordas_h4)

plt.show()
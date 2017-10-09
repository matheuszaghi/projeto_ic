# -*- coding: utf-8 -*-
"""

Created on 09/10/2017

@Author: Carlos Eduardo Barbosa

Test imports and first pPXF run.
"""
import os

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

from ppxf.ppxf import ppxf

import context

def test_image():
    """ First test to import modules and read fits file. """
    filename = os.path.join(context.data_dir,
                            "ADP.2017-06-16T13:59:19.245.fits")
    image = fits.getdata(filename)
    image = np.ma.array(image, mask=np.isnan(image))
    vmin = np.percentile(image[~np.isnan(image)], 1)
    vmax = np.percentile(image[~np.isnan(image)], 99)
    im = plt.imshow(image, origin="bottom", vmin=vmin, vmax=vmax)
    plt.minorticks_on()
    plt.colorbar(im)
    plt.show()


if __name__ == "__main__":
    test_image()
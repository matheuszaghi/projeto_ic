# -*- coding: utf-8 -*-
"""

Created on 09/10/2017

@Author: Carlos Eduardo Barbosa

Test imports and first pPXF run.
"""
import os

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from mpdaf.obj import Spectrum, WaveCoord

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

def test_ppxf():
    """ Test loading data """
    filename = os.path.join(context.data_dir, "test_slice.fits")
    # Loading data and header from cube
    data = fits.getdata(filename, 1)
    header = fits.getheader(filename, 1)
    # Picking one spectrum
    xpix = 180
    ypix = 157
    specdata = data[:,ypix-1,xpix-1]
    # Wavelenght data
    wave = WaveCoord(cdelt=header["cd3_3"], crval=header["crval3"],
                     cunit= u.angstrom)
    wave1 = (np.arange(header["naxis3"])) * header["cd3_3"] + header["crval3"]
    spec = Spectrum(wave=wave, data=specdata)
    ax = plt.subplot(111)
    spec.plot()
    ax.plot(wave1, specdata, "-")
    plt.show()
    # print(header["naxis1"], header["naxis2"], header["naxis3"], header["naxis"])
    # zdim, ydim, xdim = data.shape

if __name__ == "__main__":
    # test_image()
    test_ppxf()
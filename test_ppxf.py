# -*- coding: utf-8 -*-
"""

Created on 09/10/2017

@Author: Carlos Eduardo Barbosa

Test imports and first pPXF run.
"""
from __future__ import print_function, division

import os


import glob
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

import ppxf
import ppxf.ppxf_util as util

import context
from muse_resolution import get_muse_fwhm, broad2res

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

def test_ppxf(filename):
    """ Test pPXF run in a cube sampling the central part of the galaxy.  """

    # Loading data and header from cube
    data = fits.getdata(filename, 1)
    header = fits.getheader(filename, 1)
    # Picking one spectrum for this test
    xpix = 1
    ypix = 1
    specdata = data[:,ypix-1,xpix-1]
    # Wavelenght data
    wave = (np.arange(header["naxis3"])) * header["cd3_3"] + header["crval3"]
    # Before rebin to log-scale, it is necessary to homogeneize the
    # MUSE spectrum to the same FWHM
    f = get_muse_fwhm()
    fwhm_data = f(wave)
    fwhm_max = fwhm_data.max()
    print ("Broadening spectra to homogeneize FHWM to {}".format(fwhm_max))
    spec = broad2res(wave, specdata, fwhm_data, fwhm_max) # Broad to 2.9 AA
    # Rebin the data to logarithm scale
    galaxy, logwave1, velscale = util.log_rebin([wave[0], wave[-1]], specdata)
    galaxy = galaxy/np.median(galaxy)  # Normalize spectrum to avoid numerical issues (??)

    file_dir = os.path.dirname(os.path.realpath(__file__))  # path of this procedure


    estrelas = glob.glob(file_dir + '/projeto_ic/models/Mbi1.30Z*.fits')


    hdu = fits.open(estrelas[0])
    data2 = hdu[0].data
    header2 = hdu[0].header
    wave2 = header2['CRVAL1'] + np.array([0., header2['CDELT1']*(header2['NAXIS1'] - 1)])


    fwhm_data2 = broad2res(wave2, estrelas, 2.51)

    data2New, logwave2, velscale_temp = util.log_rebin(wave2, data2, velscale=velscale/velscale_ratio)
    templates = np.empty((data2New.size, len(estrelas)))





    #FWHM_dif = np.sqrt(fwhm_max**2 - FWHM_tem**2)     # NAO ENTENDI
    #sigma = FWHM_dif/2.355/header2['CDELT1']          # ESSA PARTE

    # Ajusta o espectro dos templates e da galaxia para o mesmo comprimento
    # de onda inicial

    c = 299792.458
    if velscale_ratio > 1:
        dv = (np.mean(logwave2[:velscale_ratio]) - logwave1[0])*c  # km/s
    else:
        dv = (logwave2[0] - logwave1[0])*c  # km/s   


    z = 0.03 # Redshift of the galaxy

    # Exclude the emission lines of the gas
    goodPixels = util.determine_goodpixels(logwave1, wave2, z)

    vel = c*np.log(1 + z)
    start = [vel, 10000.]   
    t = clock()

    pp = ppxf(templates, galaxy, noise, velscale, start,
    	      goodpixels=goodPixels, plot=True, moments=4,
    	      degree=4, vsyst=dv, velscale_ratio=velscale_ratio)

    print("Formal errors:")
    print("     dV    dsigma   dh3      dh4")
    print("".join("%8.2g" % f for f in pp.error*np.sqrt(pp.chi2)))

    print('Elapsed time in pPXF: %.2f s' % (clock() - t))

if __name__ == "__main__":
    # test_image()
    filename = os.path.join(context.data_dir, "cube_10x10.fits")
    test_ppxf(filename)
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

from ppxf import ppxf
import ppxf_util as util

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
    xpix = 1
    ypix = 1
    specdata = data[:,ypix-1,xpix-1]
    # Wavelenght data
    wave = WaveCoord(cdelt=header["cd3_3"], crval=header["crval3"],
                     cunit= u.angstrom)
    wave1 = (np.arange(header["naxis3"])) * header["cd3_3"] + header["crval3"]
    spec = Spectrum(wave=wave, data=specdata)
    #ax = plt.subplot(111)
    #spec.plot()
    #ax.plot(wave1, specdata, "-")
    #plt.show()
    # print(header["naxis1"], header["naxis2"], header["naxis3"], header["naxis"])
    # zdim, ydim, xdim = data.shape
    ####################

    FWHM_gal = 4.2 ###CHUTE### VERIFICAR

    galaxy, logwave1, velscale = util.log_rebin(wave1, specdata)
    galaxy = galaxy/np.median(galaxy)  # Normalize spectrum to avoid numerical issues (??)

    estrelas = os.path.join(context.data_dir, "/models/Mbi1.30Z*.fits")


    FWHM_tem = 2.51 ###CHUTE### VERIFICAR


    velscale_ratio = 2  # ????


    #

    hdu = fits.open(estrelas[0])
    data2 = hdu[0].data
    header2 = hdu[0].header
    wave2 = h2['CRVAL1'] + np.array([0., h2['CDELT1']*(h2['NAXIS1'] - 1)])
    data2New, logwave2, velscale_temp = util.log_rebin(wave2, data2, velscale=velscale/velscale_ratio)
    templates = np.empty((sspNew.size, len(estrelas)))



    FWHM_dif = np.sqrt(FWHM_gal**2 - FWHM_tem**2)     # NAO ENTENDI
    sigma = FWHM_dif/2.355/h2['CDELT1']               # ESSA PARTE


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
    test_ppxf()
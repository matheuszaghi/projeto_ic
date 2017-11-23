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

from ppxf.ppxf import ppxf
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
    fwhm_max = fwhm_data.max() + 0.01
    print ("Broadening spectra to homogeneize FHWM to {}".format(fwhm_max))
    specdata = broad2res(wave, specdata, fwhm_data, fwhm_max)[0] # Broad to 2.91 AA
    # Rebin the data to logarithm scale
    velscale = 30.
    galaxy, logwave1, velscale = util.log_rebin([wave[0], wave[-1]], specdata, velscale=velscale)
    galaxy = galaxy/np.median(galaxy)  # Normalize spectrum to avoid numerical issues (??)
    # Preparing templates: TODO: make this part a function outside of the fitting routine
    # Start checking if file exists and other conditions
    redo = False
    file_dir = os.path.dirname(os.path.realpath(__file__))  # path of this procedure
    outdir = os.path.join(file_dir, "templates")
    template_file = os.path.join(outdir, "templates_velscale{}.fits".format(int(velscale)))
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    if os.path.exists(template_file) and not redo:
        templates = fits.getdata(template_file, 0)
        logwave2 = fits.getdata(template_file, 1)
        ntemp = templates.shape[1]
    else:
        if context.username == "kadu":
            ssps = glob.glob(os.path.join(file_dir, "ppxf/miles_models/Mun*.fits"))
        else:
            ssps = glob.glob(os.path.join(file_dir, "models/Mbi1.30Z*.fits"))
        ntemp = len(ssps)
        header = fits.getheader(ssps[0], 0)
        wave2 = header['CRVAL1'] + np.arange(header['NAXIS1']) * header['CDELT1']
        for i, ssp in enumerate(ssps):
            print("Template {} / {}".format(i+1, ntemp))
            data = fits.getdata(ssp, 0)
            newssp = broad2res(wave2, data, 2.51, fwhm_max)[0]
            newssp, logwave2, velscale = util.log_rebin([wave2[0], wave2[-1]],
                                                        newssp, velscale=velscale)
            if i == 0:
                templates = np.zeros((len(logwave2), ntemp))
            templates[:,i] = newssp / np.median(newssp)
        hdu = fits.PrimaryHDU(templates)
        hdu2 = fits.ImageHDU(logwave2)
        hdulist = fits.HDUList([hdu, hdu2])
        hdulist.writeto(template_file, overwrite=True)
    # Calculating velocity shift between templates and galaxy spectrum
    c = 299792.458
    dv = (logwave2[0] - logwave1[0])*c  # km/s
    z = 0.034 # Redshift of the galaxy used for initial guess

    # Exclude the emission lines of the gas
    # goodPixels = util.determine_goodpixels(logwave1, wave2, z)
    goodPixels = np.ones_like(galaxy)
    vel = c*np.log(1 + z)
    start = [vel, 100., 0., 0.]
    noise = np.ones_like(galaxy)
    pp = ppxf(templates, galaxy, noise, velscale, start,
    	      goodpixels=None, plot=True, moments=4,
    	      degree=20, vsyst=dv, clean=True)
    # print("Formal errors:")
    # print("     dV    dsigma   dh3      dh4")
    # print("".join("%8.2g" % f for f in pp.error*np.sqrt(pp.chi2)))
    plt.show()

if __name__ == "__main__":
    # test_image()
    velscale = 30.
    filename = os.path.join(context.data_dir, "cube_10x10.fits")
    test_ppxf(filename)
# -*- coding: utf-8 -*-
"""
Program to run pPXF over a data cube.

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

def load_templates(velscale, fwhm, redo=False):
    """ Load templates to be used in the fitting process. """
    file_dir = os.path.dirname(os.path.realpath(__file__))
    template_file = os.path.join(file_dir, "templates",
             "templates_velscale{}_fwhm{}.fits".format(int(velscale), fwhm))
    if os.path.exists(template_file) and not redo:
        return template_file
    # TODO: use the same set of models for both users
    if context.username == "kadu":
        ssps = glob.glob(os.path.join(file_dir, "ppxf/miles_models/Mun*.fits"))
    else:
        ssps = glob.glob(os.path.join(file_dir, "models/Mbi1.30Z*.fits"))
    ntemp = len(ssps)
    header = fits.getheader(ssps[0], 0)
    wave2 = header['CRVAL1'] + np.arange(header['NAXIS1']) * header['CDELT1']
    for i, ssp in enumerate(ssps):
        print("Template {} / {}".format(i + 1, ntemp))
        data = fits.getdata(ssp, 0)
        newssp = broad2res(wave2, data, 2.51, fwhm)[0]
        newssp, logwave2, velscale = util.log_rebin([wave2[0], wave2[-1]],
                                                    newssp, velscale=velscale)
        if i == 0:
            templates = np.zeros((len(logwave2), ntemp))
        templates[:, i] = newssp / np.median(newssp)
    hdu = fits.PrimaryHDU(templates)
    hdu2 = fits.ImageHDU(logwave2)
    hdulist = fits.HDUList([hdu, hdu2])
    hdulist.writeto(template_file, overwrite=True)
    return template_file

def run_ppxf (filename):
    """ Run pPXF for all spectra in a given filename. """
    # Constants to be used in routine if not optional
    velscale = 30. # km/s
    c = 299792.458
    z = 0.034  # Redshift of the galaxy used for initial guess
    vel = c * np.log(1 + z)
    start = [vel, 100., 0., 0.]
    # Setting the FWHM of the fitting
    f = get_muse_fwhm()
    fwhm_data = f(np.linspace(4500, 10000, 1000))
    fwhm_max = np.round(fwhm_data.max() + 0.01) # This will use FWHM=3.0
    ###########################################################################
    # Templates are loaded only once, they should not be included in the loop
    template_file = load_templates(velscale, fwhm_max)
    templates = fits.getdata(template_file, 0)
    logwave2 = fits.getdata(template_file, 1)
    ntemp = templates.shape[1]
    # Loading data and header from cube 
    data = fits.getdata(filename, 1)
    header = fits.getheader(filename, 1)
    ###########################################################################
    # How to iterate over complete array using only one loop
    zdim, ydim, xdim = data.shape
    pixels = np.array(np.meshgrid(np.arange(xdim)+1,
                      np.arange(ydim)+1)).reshape((2,-1)).T
    ###########################################################################
    # Setting output lists
    log_dir = os.path.join(context.plots_dir, "ppxf_results")
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)
    finalList = []
    tempList = []
    for xpix,ypix in pixels:
        # Picking one spectrum for this test
        specdata = data[:,ypix-1,xpix-1]
        # Wavelenght data
        wave = (np.arange(header["naxis3"])) * header["cd3_3"] + \
               header["crval3"]
        print ("Broadening spectra to homogeneize FHWM to {}".format(
               fwhm_max))
        specdata = broad2res(wave, specdata, fwhm_data, fwhm_max)[0] # Broad
        # to 3.0 AA
        # Rebin the data to logarithm scale
        galaxy, logwave1, velscale = util.log_rebin([wave[0], wave[-1]],
                                               specdata, velscale=velscale)
        galaxy = galaxy/np.median(galaxy)  # Normalize spectrum to avoid numerical issues (??)
        dv = (logwave2[0] - logwave1[0])*c  # km/s
        # Exclude the emission lines of the gas
        wavetemp = np.exp(logwave2)
        # goodPixels = np.ones_like(galaxy)
        noise = np.ones_like(galaxy)
        pp = ppxf(templates, galaxy, noise, velscale, start,
                  goodpixels=None, plot=True, moments=4,
                  degree=8, vsyst=dv, clean=True)
        #print("Formal errors:")
        #print("     dV    dsigma   dh3      dh4")
        #print("".join("%8.2g" % f for f in pp.error*np.sqrt(pp.chi2)))

        solutionList = ["%.2f" % x for x in pp.sol]
        print(solutionList)
        plt.savefig(os.path.join(log_dir, 'ppxf_x{}_y{}.png'.format(xpix,
                                                                    ypix)))
        plt.clf()


        tempList.append(solutionList)

    finalList.append(tempList)
    with open ('resultadosMatia.txt', 'w') as output:
        for x in range(len(finalList)):
            for y in range(len(finalList[x])):
                print(x + 1, y + 1, file=output)
                print(finalList[x][y], file=output)




if __name__ == "__main__":
    velscale = 30.
    filename = os.path.join(context.data_dir, "cube_10x10.fits")
    run_ppxf(filename)
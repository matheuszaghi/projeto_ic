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
from astropy.table import Table, hstack

from ppxf.ppxf import ppxf
import ppxf.ppxf_util as util


import context
from muse_resolution import get_muse_fwhm, broad2res

def load_templates(velscale, fwhm, redo=False):
    """ Load templates to be used in the fitting process. """
    file_dir = os.path.dirname(os.path.realpath(__file__))
    template_dir = os.path.join(context.data_dir, "templates")
    if not os.path.exists(template_dir):
        os.mkdir(template_dir)
    output = os.path.join(template_dir,
             "templates_velscale{}_fwhm{}.fits".format(int(velscale), fwhm))
    if os.path.exists(output) and not redo:
        return output

    ssps = glob.glob(os.path.join(file_dir, "ppxf/miles_models/Mun*.fits")) 

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
    hdulist.writeto(output, overwrite=True)
    return output

def run_ppxf (filename):
    """ Run pPXF for all spectra in a given filename. """
    # Constants to be used in routine if not optional
    velscale = 30. # km/s
    c = 299792.458
    z = 0.034  # Redshift of the galaxy used for initial guess
    vel = c * np.log(1 + z)
    start = [vel, 100., 0., 0.]
    start = [start, start]
    # Setting the FWHM of the fitting
    f = get_muse_fwhm()
    fwhm_data = f(np.linspace(4500, 10000, 1000))
    fwhm_max = np.round(fwhm_data.max() + 0.01) # This will use FWHM=3.0
    ###########################################################################
    # Templates are loaded only once, they should not be included in the loop
    template_file = load_templates(velscale, fwhm_max)
    templates = fits.getdata(template_file, 0)
    logwave2 = fits.getdata(template_file, 1)
    n_ssps = templates.shape[1]
    # Loading data and header from cube 
    data = fits.getdata(filename, 1)
    header = fits.getheader(filename, 1)
    wave = (np.arange(header["naxis3"])) * header["cd3_3"] + \
             header["crval3"]
    # Loading emission line template
    spec0 = data[:,0,0]
    logwave = util.log_rebin([wave[0], wave[-1]], spec0, velscale=velscale)[1]
    lam = np.exp(logwave)
    gas_templates, gas_names, line_wave = \
        util.emission_lines(logwave2, [lam[0], lam[-1]], fwhm_max)
    n_gas = gas_templates.shape[1]
    # Adding emission templates to the stellar templates
    templates = np.column_stack([templates, gas_templates])
    # Specifying components
    components = np.ones(n_ssps + n_gas, dtype=np.int)
    components[:n_ssps] = 0
    ###########################################################################
    # How to iterate over complete array using only one loop
    zdim, ydim, xdim = data.shape
    pixels = np.array(np.meshgrid(np.arange(xdim)+1,
                      np.arange(ydim)+1)).reshape((2,-1)).T
    pixels = Table(pixels, names=["x", "y"])
    ###########################################################################
    # Sky lines
    skylines = np.array([5577])
    # Setting output lists
    log_dir = os.path.join(context.plots_dir, "ppxf_results")
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)
    test_dir = os.path.join(context.plots_dir, "ppxf_values")
    if not os.path.exists(test_dir):
        os.mkdir(test_dir)
    tempList = []
    values = []
    sols = []
    sol_gas = []
    tmpssps = []
    tmpgas = []
    for i, (xpix,ypix) in enumerate(pixels):
        print(xpix, ypix)
        output = os.path.join(test_dir, 'gas_ppxf_x{}_y{}.fits'.format(xpix, ypix))
        if os.path.exists(output):
            print('Already calculated values...')
            continue
        # Picking one spectrum for this test
        specdata = data[:,ypix-1,xpix-1]
        # Wavelenght data
        print ("Broadening spectra to homogeneize FHWM to {}".format(
               fwhm_max))
        specdata = broad2res(wave, specdata, fwhm_data, fwhm_max)[0] # Broad
        # to 3.0 AA
        # Rebin the data to logarithm scale
        galaxy, logwave1, velscale = util.log_rebin([wave[0], wave[-1]],
                                               specdata, velscale=velscale)
        lam = np.exp(logwave1)
        badpixels = []
        for skyline in skylines:
            idx1 = np.where(lam < skyline + 15)[0]
            idx2 = np.where(lam > skyline - 15)[0]
            idx = np.intersect1d(idx1, idx2)
            badpixels.append(idx)
        badpixels = np.unique(np.hstack(badpixels))

        goodpixels = np.arange(len(galaxy))
        goodpixels = np.delete(goodpixels, badpixels)

        #selecting non finite numbers
        index_galaxy_nan = np.where(np.isnan(galaxy))
        index_galaxy_inf = np.where(np.isinf(galaxy))
        index_galaxy_ninf = np.where(np.isneginf(galaxy))

        #setting nans to zero
        galaxy[np.isnan(galaxy)] = 0
        galaxy[np.isinf(galaxy)] = 0
        galaxy[np.isneginf(galaxy)] = 0
        
        #nans cant be used in goodpixels
        goodpixels = np.delete(goodpixels, index_galaxy_nan)
        goodpixels = np.delete(goodpixels, index_galaxy_inf)
        goodpixels = np.delete(goodpixels, index_galaxy_ninf)


        galaxy = galaxy/np.median(galaxy)  # Normalize spectrum to avoid numerical issues (??)
        dv = (logwave2[0] - logwave1[0])*c  # km/s
        # Exclude the emission lines of the gas
        wavetemp = np.exp(logwave2)
        noise = np.ones_like(galaxy)
        pp = ppxf(templates, galaxy, noise, velscale, start,
                  goodpixels=goodpixels, plot=True, moments=4,
                  degree=8, vsyst=dv, clean=True, lam=lam,
                  component=components)

        tmpssps = Table(pp.sol[0])
        tmpgas = Table(pp.sol[1])
        tmpssps.write(os.path.join(test_dir, 'ssps_ppxf_x{}_y{}.fits'.format(xpix, ypix)))
        tmpgas.write(os.path.join(test_dir, 'gas_ppxf_x{}_y{}.fits'.format(xpix, ypix)))

        sols.append(pp.sol[0]) # Only SSPs
        sol_gas.append(pp.sol[1]) #Only gas

        plt.savefig(os.path.join(log_dir, 'ppxf_x{}_y{}.png'.format(xpix, ypix)))
        plt.clf()

    #saving gas results in a fits table
    sol_gas = Table(np.array(sol_gas), names=['vel', 'sigma', 'h3', 'h4'])
    table_gas = hstack([pixels, sol_gas])
    table_gas.write('Resultados/Resultados_gas.fits', format='fits', overwrite=True)
    #saving ssps results in a fits table
    sols = Table(np.array(sols), names=['vel', 'sigma', 'h3', 'h4'])
    table = hstack([pixels, sols])
    table.write('Resultados/Resultados_ssps.fits', format='fits', overwrite=True)
    print(table_gas)
    print(table)

if __name__ == "__main__":
    velscale = 30.
    filename = os.path.join(context.data_dir, "cube_85x85.fits")
    run_ppxf(filename)
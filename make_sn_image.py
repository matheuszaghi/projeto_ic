# -*- coding: utf-8 -*-
""" 

Created on 11/05/18

Author : Carlos Eduardo Barbosa

"""
from __future__ import print_function, division

import os
import itertools

import numpy as np
from astropy.io import fits

import context
from der_snr import DER_SNR

def make_sn_image():
    """ Put short description here. """
    cubo = os.path.join(context.data_dir, "cube_10x10.fits")
    flux = fits.getdata(cubo, 1)
    ############################################################################
    # Solução proposta
    zdim, ydim, xdim = flux.shape
    snimage = np.zeros((ydim, xdim))
    for j, i in itertools.product(range(ydim), range(xdim)):
        snimage[j,i] = DER_SNR(flux[:, j, i])
    # Na prática, vamos usar a função collapse cube abaixo
    ############################################################################

    #criando vetor 2d
    snr = []
    for i in range(10):
        snr.append([])
        for j in range(10):
            snr[i].append(0)



    #as dimensoes do cubo estao na forma (z,y,x)
    #loop da variavel y
    for j in range(0,10):
        #loop da variavel x
        for i in range(0,10):
            snr[i][j] = DER_SNR(flux[:, j, i]) #calculo do sinal-ruido
    #salvando imagem fits
    final_image = fits.PrimaryHDU(snr)
    final_image.writeto('Image_SN.fits', overwrite=True)

def collapse_cube(cubename, outfile, redo=False):
    """ Collapse a MUSE data cube to produce a white-light image and a
    noise image.

    The noise is estimated with the same definition of the DER_SNR algorithm.

    Input Parameters
    ----------------
    cubename : str
        Name of the MUSE data cube

    outfile : str
        Name of the output file

    redo : bool (optional)
        Redo calculations in case the oufile exists.
    """
    if os.path.exists(outfile) and not redo:
        return
    data = fits.getdata(cubename, 1)
    h = fits.getheader(cubename, 1)
    h2 = fits.getheader(cubename, 2)
    h["NAXIS"] = 2
    del h["NAXIS3"]
    h2["NAXIS"] = 2
    del h2["NAXIS3"]
    print("Starting collapsing process...")
    newdata = np.nanmedian(data, axis=0)
    noise = 1.482602 / np.sqrt(6.) * np.nanmedian(np.abs(2.* data - \
           np.roll(data, 2, axis=0) - np.roll(data, -2, axis=0)), \
           axis=0)
    hdu = fits.PrimaryHDU(newdata, h)
    hdu2 = fits.ImageHDU(noise, h2)
    hdulist = fits.HDUList([hdu, hdu2])
    hdulist.writeto(outfile, overwrite=True)
    return


    
if __name__ == "__main__":
    make_sn_image()
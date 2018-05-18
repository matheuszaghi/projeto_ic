# -*- coding: utf-8 -*-
""" 

Created on 11/05/18

Author : Carlos Eduardo Barbosa

"""

import context

import os

import numpy as np

from astropy.io import fits

from der_snr import DER_SNR

def make_sn_image():
    """ Put short description here. """
    cubo = os.path.join(context.data_dir, "cube_10x10.fits")
    flux = fits.getdata(cubo, 1)

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

    
if __name__ == "__main__":
    make_sn_image()
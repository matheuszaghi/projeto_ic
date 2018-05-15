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
    #header = fits.getheader(cubo, 1)
    #wave = (np.arange(header["naxis3"])) * header["cd3_3"] + \
    #         header["crval3"]
    

    print(flux)
    for arrayList in flux:

        for array in arrayList:
            snr = DER_SNR(array)
            print(snr)

    
if __name__ == "__main__":
    make_sn_image()
# -*- coding: utf-8 -*-
"""
Created on 01/11/2017
@Author: Carlos Eduardo Barbosa
Test routine for cube slicing
"""
from __future__ import print_function, division

import os

from mpdaf.obj import Cube

import context

def test_slicing(filename, size=10, ext=1):
    """ Test slicing of the cube with MPDAF. """
    cube = Cube(cubename, ext=ext)
    center = (7.02162, 229.186)
    lbda = (cube.wave.get_range()[0], 5900)
    cube = cube.subcube(center, size, lbda)
    newcubename = os.path.join(context.data_dir, "cube_{0}arcsec.fits".format(
        size))
    cube.write(newcubename)
    return

if __name__ == "__main__":
    cubename = os.path.join(context.data_dir,
                            "ADP.2017-06-16T13:59:19.244.fits")
    test_slicing(cubename, size=16)
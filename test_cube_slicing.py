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

def test_slicing():
    """ Test slicing of the cube with MPDAF. """
    cubename = os.path.join(context.data_dir,
                            "ADP.2017-06-16T13:59:19.244.fits")
    cube = Cube(cubename, ext=1)
    zdim, ydim, xdim = cube.shape
    center = (ydim // 2, xdim // 2)
    newcube = cube.select_lambda(cube.wave.get_range()[0], 7000)
    newcubename = os.path.join(context.data_dir, "test_slice.fits")
    newcube.write(newcubename)
    return

if __name__ == "__main__":
    test_slicing()

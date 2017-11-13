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
                            "ADP.2017-06-16T13_59_19.244.fits")


    cube = Cube(cubename, ext=1)

    center = (7.02162, 229.186)
    lbda = (cube.wave.get_range()[0], 5900)

    #gera um cubo 10x10 pixeis para teste do ppxf
   
    cube = cube.subcube(center, 2, lbda)

    newcubename = os.path.join(context.data_dir, "test_slice.fits")
    cube.write(newcubename)
    return

if __name__ == "__main__":
    test_slicing()
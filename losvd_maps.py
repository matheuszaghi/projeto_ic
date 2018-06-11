# -*- coding: utf-8 -*-
""" 

Created on 11/06/18

Author : Carlos Eduardo Barbosa

"""
from __future__ import print_function, division

import os

import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

import context

def read_results(redo=False):
    """ Read tables from pPXF and produces an array to be plotted"""
    data_dir = os.path.join(context.data_dir, "ppxf_values")
    moments = ["vel", "sigma", "h3", "h4"]
    outputs = [os.path.join(data_dir, "{}.fits".format(_)) for _ in
               moments]
    if all([os.path.exists(_) for _ in outputs]):
        results = []
        for out in outputs:
            data = fits.getdata(out)
            results.append(data)
        return results
    ssps = [_ for _ in os.listdir(data_dir) if _.startswith("ssps")]
    xs = np.unique([int(_.split("_")[2][1:]) for _ in ssps])
    ys = np.unique([int(_.split("_")[3][1:-5]) for _ in ssps])
    vel = np.zeros((ys[-1], xs[-1])) * np.nan
    sigma = np.zeros_like(vel) * np.nan
    h3 = np.zeros_like(vel) * np.nan
    h4 = np.zeros_like(vel) * np.nan
    for x in xs:
        for y in ys:
            filename = os.path.join(data_dir,
                                    "ssps_ppxf_x{}_y{}.fits".format(x, y))
            if not(os.path.exists(filename)):
                continue
            data = Table.read(filename)
            vel[y-1,x-1] = data["col0"]
            sigma[y-1,x-1] = data["col1"]
            h3[y-1, x-1] = data["col2"]
            h4[y-1,x-1] = data["col3"]
    results = [vel, sigma, h3, h4]
    for res, out in zip(results, outputs):
        hdu = fits.PrimaryHDU(res)
        hdu.writeto(out, overwrite=True)
    return [vel, sigma, h3, h4]

def make_maps(moments):
    """ Produces the plots. """
    lims = [[10140, 10160], [210, 320], [-0.05, 0.05], [-0.04, 0.08]]
    labels = ["$\Delta V$ (km/s)", "$\sigma$ (km/s)", "$h_3$", "$h_4$"]
    names = ["vel", "sigma", "h3", "h4"]
    meanV = np.mean(lims[0])
    moments[0] -= meanV
    lims[0] = [lims[0][0] - meanV, lims[0][1] - meanV]
    for i, mom in enumerate(moments):
        out = os.path.join(context.plots_dir, "{}.png".format(names[i]))
        fig = plt.figure(i+1, figsize=(5,4))
        ax = plt.subplot(111)
        ax.minorticks_on()
        im = ax.imshow(gaussian_filter(mom, sigma=1.2),
                      origin="bottom",
                   cmap="Spectral_r", vmin=lims[i][0], vmax=lims[i][1])
        ax.set_xlabel("X (pixel)")
        ax.set_ylabel("Y (pixel)")
        cbar = plt.colorbar(im)
        cbar.set_label(labels[i])
        plt.subplots_adjust(left=0.1, bottom=0.1, top=0.98, right=0.96)
        plt.savefig(out, dpi=200)
    # plt.show()


if __name__ == "__main__":
    results = read_results()
    make_maps(results)
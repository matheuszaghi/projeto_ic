



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


def run_ppxf (filename):
    # Loading data and header from cube
    data = fits.getdata(filename, 1)
    header = fits.getheader(filename, 1)

    finalList = []
    tempList = []
    
    for xpix in range(1, 3):
        for ypix in range(1, 3):
    
            # Picking one spectrum for this test
            specdata = data[:,ypix-1,xpix-1]
            # Wavelenght data
            wave = (np.arange(header["naxis3"])) * header["cd3_3"] + header["crval3"]
            f = get_muse_fwhm()
            fwhm_data = f(wave)
            fwhm_max = fwhm_data.max() + 0.01
            print ("Broadening spectra to homogeneize FHWM to {}".format(fwhm_max))
            specdata = broad2res(wave, specdata, fwhm_data, fwhm_max)[0] # Broad to 2.91 AA
            # Rebin the data to logarithm scale
            velscale = 30.
            galaxy, logwave1, velscale = util.log_rebin([wave[0], wave[-1]], specdata, velscale=velscale)
            galaxy = galaxy/np.median(galaxy)  # Normalize spectrum to avoid numerical issues (??)


            #### VERIFCAR ####
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


            #### VERIFICAR ####





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
            

            #print("Formal errors:")
            #print("     dV    dsigma   dh3      dh4")
            #print("".join("%8.2g" % f for f in pp.error*np.sqrt(pp.chi2)))

            solutionList = ["%.2f" % x for x in pp.sol]
            print(solutionList)
            #plt.show()


            plt.savefig('plot_figures/Plot_Pixel({},{})'.format(xpix,ypix))
            plt.clf()


            tempList.append(solutionList)

        finalList.append(tempList)
        tempList = []

    with open ('resultadosMatia.txt', 'w') as output:
        for x in range(len(finalList)):
            for y in range(len(finalList[x])):
                print(x + 1, y + 1, file=output)
                print(finalList[x][y], file=output)




if __name__ == "__main__":
    velscale = 30.
    filename = os.path.join(context.data_dir, "cube_10x10.fits")
    run_ppxf(filename)
#!/usr/bin/env python

import os
import sys
#import numpy as np
#import scipy
import warnings
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
#import photutils
#from skimage import measure, morphology
#from scipy.cluster import vq
#import gaps
#import matplotlib.pyplot as plt
#import pylab
#from scipy.signal import argrelmin
#import scipy.ndimage.filters as filters
import time


if __name__ == "__main__":
    #import argparse
    #parser = argparse.ArgumentParser(description="Run DAOPHOT PSF photometry on FITS image.")
    #
    #parser.add_argument('--file', '-f', action="store", help="The FITS file to process", default=None)
    ##parser.add_argument('--datarepo', '-d', action="store", help="The data repository directory", default="/data/lsst/decam/redux/cp/cosmos/")
    ##parser.add_argument('--outfile', '-o', action="store", help="The output filename for the metrics.", default="qametrics.csv")
    #parser.add_argument('--verbose', '-v', action="store_true", help="Print out the data as it is gathered.", default=False)
    #parser.add_argument('--clobber', '-c', action="store_true", help="Overwrite the output file if it already exists.", default=False)
    #
    #args = parser.parse_args()
    #file = args.file
    
    # I need to update the FITS header with the fpack settings
    # For the flux images *_??.fits and *_??_comb.fits use these settings:
    # sxaddpar,head,'FZALGOR','RICE_1'
    # sxaddpar,head,'FZQMETHD','SUBTRACTIVE_DITHER_1'
    # sxaddpar,head,'FZQVALUE',4
    # sxaddpar,head,'FZDTHRSD','CHECKSUM'
    # For the bad pixel mask *_??_comb.bpm.fits don't add anything
    # since these are integers and will be automatically compressed
    # losslessly with PLIO.

    print "Adding fpack keywords to SMASH FITS headers"
    t0 = time.time()

    # Find all of the fits files
    #basedir = "/datalab/users/dnidever/smash/dr1/"
    #basedir = "/data/smash/dr1/"
    basedir = "/data/smash/cp/red/photred/test/"

    warnings.simplefilter('ignore', category=AstropyWarning)

    for root, dirs, files in os.walk(basedir):
        for file in files:
            if file.endswith(".fits"):
                if not file.endswith("_comb.bpm.fits"):
                    fullfile = os.path.join(root, file)
                    print fullfile, " adding fpack keyword values"

                    # Apply changes to header (MODIFY IN PLACE)
                    #hdulist = fits.open(fullfile, mode='update') # modify IN PLACE
                    #fitshdr = hdulist[0].header # use only first in list
		    data, fitshdr = fits.getdata(fullfile,header=True)
                    # add the keywords
                    fitshdr['FZALGOR'] = 'RICE_1'
                    fitshdr['FZQMETHD'] = 'SUBTRACTIVE_DITHER_1'
                    fitshdr['FZQVALUE'] = 8
                    fitshdr['FZDTHRSD'] = 'CHECKSUM'
                    # remove "BEGIN MAIN ...." and "BEGIN EXTENSION ..."
                    for val in fitshdr:
                        if val[0:5] == "BEGIN":
                            fitshdr.remove(val)
                    #hdulist.close(output_verify='ignore')         # now FITS header is MODIFIED
		    fits.writeto(fullfile,data,header=fitshdr,clobber=True)
                    sys.exit()
                else:
                    print fullfile, "do nothing"

    print time.time()-t0

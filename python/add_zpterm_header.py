#!/usr/bin/env python

import os
import sys
import numpy as np
import warnings
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
from astropy.table import Table
import time
import shutil
import re
#import subprocess
import glob
#import logging
#import socket
#from scipy.signal import convolve2d
#from scipy.ndimage.filters import convolve

# Add the ZEROPOINT to the FITS header using the "chips" file

if __name__ == "__main__":

    t0 = time.time()

    #print sys.argv

    # Not enough inputs
    n = len(sys.argv)
    if n < 2:
        print "Syntax - add_zpterm_header.py fitsfile chipfile"
        sys.exit()

    # File names
    fitsfile = sys.argv[1]
    chipfile = sys.argv[2]
    # Check that the files exist
    if os.path.exists(fitsfile) == False:
        print fitsfile, "file NOT FOUND"
        sys.exit()
    if os.path.exists(chipfile) == False:
        print chipfile, "file NOT FOUND"
        sys.exit()
    
    #
    print("Updating ZERO-POINT in header for "+fitsfile)

    # Load the information on the chips, make it an astropy table
    allchips = Table(fits.getdata(chipfile,1))
    # Find the record for this chip
    base = os.path.basename(fitsfile)
    base = os.path.splitext(os.path.splitext(base)[0])[0]
    ind, = np.where(allchips['BASE'] == base)
    nind = len(ind)
    if nind == 0:
        print(base+" not found in "+chipfile)
        sys.exit()
    chip = allchips[ind[0]]

    # Load the FITS header
    data, head = fits.getdata(fitsfile, ext=0, header=True)

    # Calculate the zero-point for this chip file
    #  The transformation equation used from SMASHRED_APPLY_PHOTTRANSEQN.PRO
    #   V = mV - v1 - v2 * XV - v3 * (B-V) - v4 * XV*(B-V) - v5 * (B-V) * (B-V)
    #         +(aperture correction) + (time correction) 
    #   The aperture corrections need to be POSITIVE
    #  outmag = inmag - t.zpterm - t.amterm*t.airmass - t.colterm*clr - t.amcolterm*t.airmass*clr - t.colsqterm*clr*clr - t.apcor + 2.5*alog10(t.exptime)
    # Leave out the aperture correction.  This is only for PSF photometry not aperture photometry.

    # We want to "remove" the color term, so use a "fiducial" color for a given band
    coltab = { 'u': 1.36, 'g': 0.6672, 'r': 0.5861, 'i': 0.3179, 'z': 0.3465}
    color = coltab[chip['BAND']]
    zpterm = chip['ZPTERM'] + chip['AIRMASS']*chip['AMTERM'] + color*chip['COLTERM']
    zpterm = -zpterm   # convert to ADDITIVE term
    # DAOPHOT adds 25.0 for the instrumental magnitudes, add that now
    zpterm += 25.0
    print('ZPTERM = %f' % zpterm)
    # The LAMBDA program will correct for the exposure time

    # Update the header
    head['MAGZERO'] = (zpterm,'SMASH ZEROPOINT')
    # Write out the file
    print('Updating '+fitsfile)
    fits.writeto(fitsfile, data, head, overwrite=True)

# Written in/for python 3.x.x
# This script uses photometry to calculate flux density and spectral indices on single epoch images using the positions of sources in the quick look catalogue
# run the script by invoking the main inputs in the format:
# python photometry_elliptical.py <input catalogue> <zeroth order FITS image> <first order FITS image>
# or run it without invoking inputs and it will prompt you to enter each input separately

import os
import sys
from termcolor import colored
import numpy as np
from astropy.stats import SigmaClip
import astropy.units as u
from astropy.io import fits
from radio_beam import Beam
from astropy.wcs import WCS
from astropy.table import Table
from astropy.coordinates import SkyCoord
from photutils import Background2D, MedianBackground
from photutils import aperture_photometry
from photutils.aperture import SkyEllipticalAperture

def beam_info(fits_image_file):
    header = fits.getheader(fits_image_file)
    beam_temp = Beam.from_fits_header(header)
    beam = beam_temp.to_header_keywords()
    
    major = float(f"{beam['BMAJ']*3600:.2f}")
    minor = float(f"{beam['BMIN']*3600:.2f}")
    pa = float(f"{beam['BPA']*3600:.2f}")
    pixel_size = float('{:.2f}'.format(np.abs(header['CDELT1'])*3600)) # CDELT1 is the change in value, in the x dirction, in beam units. Therefore each step (pixel) is equal to |CDELT1| beam units. Multiply 3600 to convert degs to arcsecs
    pixels_per_beam = beam['BMAJ']/np.abs(header['CDELT1']) # number of pixels across the beam
    return major, minor, pa, pixel_size;

def do_photometry(image, table):
    hdu = fits.open(image)
    if hdu[0].data.shape[0] == 1 or hdu[0].data.shape[1] == 1:
        hdu[0].data = hdu[0].data.reshape(hdu[0].data.shape[-2], hdu[0].data.shape[-1]) # change shape to (n,n) from (1,1,n,n) or (1,n,n)
    data = u.Quantity(hdu[0].data, unit=hdu[0].header['BUNIT'])
    w = WCS(hdu[0].header).celestial

    a, b, pa, pixel_size = beam_info(image)
    positions = SkyCoord(table['RA_Source'], table['DEC_Source'], frame='fk5', unit='deg')
    
    sigma_clip = SigmaClip(sigma=3)
    bkg_estimator = MedianBackground()
    bkg = Background2D(data, (114,114), filter_size=(19,19), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    data -= bkg.background

    apertures = SkyEllipticalAperture(positions, a*u.arcsec, b*u.arcsec, theta=pa*u.arcsec)
    photometry_table = aperture_photometry(data, apertures, wcs=w)
    factor = 1.133*(a*b)/pow(pixel_size,2) # I've used this to define the beam for photutils
    photometry_table['aperture_sum'] /= factor
    return photometry_table;

if (len(sys.argv) > 1 and len(sys.argv) < 4) or len(sys.argv) > 4:
    print(colored("\nWrong format used! Input format required:", "red"), colored('\n\tIn terminal','blue'),'\n\t\t$ python photometry_elliptical.py <input catalogue> <zeroth order FITS image> <first order FITS image>', colored('\n\n\tOr in iPython','blue'), colored('\n\t\t In[1]:','green'), ' %run photometry_elliptical.py <input catalogue> <zeroth order FITS image> <first order FITS image>\n')
    print(colored("Alternatively, run the script without invoking inputs and it will prompt you to enter each input separately!\n", "red"))
    print('\t\t $python photometry_elliptical.py\n')
    sys.exit()
try:
    catalogue_SE = sys.argv[1]
except:
    print('\n\n')
    os.system('ls *.csv')
    catalogue_SE = input("\nEnter a filename for the input catalogue:\n>")
try:
    fits_tt0 = sys.argv[2]
except:
    print('\n\n')
    os.system('ls *tt0**fits')
    fits_tt0 = input("\nEnter a filename for the zeroth order single ephoch image:\n>")
try:
    fits_tt1 = sys.argv[3]
except:
    print('\n\n')
    os.system('ls *tt1**fits')
    fits_tt1 = input("\nEnter a filename for the first order single epoch image:\n>")

target = fits_tt0.split('_')[0]
table = Table.read(catalogue_SE)

photometry_on_tt0 = do_photometry(fits_tt0, table)
photometry_on_tt1 = do_photometry(fits_tt1, table)

print('\nPhotometry on tt0 image\n', photometry_on_tt0)
print('\nPhotometry on tt1 image\n', photometry_on_tt1)

table_tt0_tt1 = photometry_on_tt0[:]
table_tt0_tt1.add_column(table['Source_name'], index=0)
table_tt0_tt1.add_column(photometry_on_tt1['aperture_sum'], rename_duplicate=True)
table_tt0_tt1.add_column(table_tt0_tt1['aperture_sum_1']/table_tt0_tt1['aperture_sum'], name='Spectral Index')

print('\nCombined table with spectral indices\n', table_tt0_tt1)
table_tt0_tt1.write(f'{target}_measured_tt0-tt1_Sv_and_alpha.csv', overwrite=True)


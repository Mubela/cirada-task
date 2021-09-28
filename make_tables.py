# Written in/for python 3.x.x
# The script accepts a catalogue and a FITS image and identifies which sources in the catalogue are present in
# the image and which are not. Either set of sources is written out to a separate table.
# Execution format: python make_tables.py <input catalogue> <zeroth order FITS image>

from astropy.coordinates  import SkyCoord
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
from termcolor import colored
from tqdm import tqdm
import sys
import os

if (len(sys.argv) > 1 and len(sys.argv) < 3) or len(sys.argv) > 3:
    print(colored("\nWrong format used! Input format required:", "red"), colored('\n\tIn terminal','blue'),'\n\t\t$ python make_tables.py <input catalogue> <zeroth order FITS image> <first order FITS image>', colored('\n\n\tOr in iPython','blue'), colored('\n\t\t In[1]:','green'), ' %run make_tables.py <input catalogue> <zeroth order FITS image>\n')
    print(colored("Alternatively, run the script without invoking inputs and it will prompt you to enter each input separately!\n", "red"))
    print('\t\t $python make_tables.py\n')
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

target = fits_tt0.split('_')[0] # used in output file names
table = Table.read(catalogue_SE)
hdu0 = fits.open(fits_tt0)

if hdu0[0].data.shape[0] == 1 or hdu0[0].data.shape[1] == 1:
    hdu0[0].data = hdu0[0].data.reshape(hdu0[0].data.shape[-2], hdu0[0].data.shape[-1])
    
w0 = WCS(hdu0[0].header).celestial

# table has >700k sources, this will take a while - needs optimizing
in_table = [] # list of sources in single epoch image
for i in tqdm(range(len(table))): # tqdm monitors the progress
    ra = table['RA_Source'][i]
    dec = table['DEC_Source'][i]
    position = SkyCoord(ra, dec, unit='deg', frame='fk5')
    in_image = position.contained_by(w0)
    if in_image:
        in_table.append(i)

# table of sources in QL catalogue but not in SE image
sources_not_in_se_image = table[:]
sources_not_in_se_image.remove_rows(in_table)
sources_not_in_se_image.write('sources_not_in_{target}.csv', overwrite=True)

# table of sources that are present in the SE image (to be used as input for photometry)
table_se_image = table[in_table]
table_se_image.write('sources_in_{target}.csv', format='csv', overwrite=True)

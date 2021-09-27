from astropy.coordinates  import SkyCoord
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
from tqdm import tqdm

catalogue = 'CIRADA_VLASS1QL_table2_hosts.csv'
fits_tt0 = 'J100200+023000_.image.pbcor.tt0.subim.fits'
fits_tt1 = 'J100200+023000_.image.pbcor.tt1.subim.fits'

table = Table.read(catalogue)
hdu0 = fits.open(fits_tt0)
hdu1 = fits.open(fits_tt1)

if hdu0[0].data.shape[0] == 1 or hdu0[0].data.shape[1] == 1:
    hdu0[0].data = hdu0[0].data.reshape(hdu0[0].data.shape[-2], hdu0[0].data.shape[-1])
if hdu1[0].data.shape[0] == 1 or hdu1[0].data.shape[1] == 1:
    hdu1[0].data = hdu1[0].data.reshape(hdu1[0].data.shape[-2], hdu1[0].data.shape[-1])
    
w0 = WCS(hdu0[0].header).celestial
w1 = WCS(hdu1[0].header).celestial

# table has >700k sources, this will take a while - needs optimizing
in_table = [] # list of sources in single epoch image
for i in range(len(table)):
    ra = table['RA_Source'][i]
    dec = table['DEC_Source'][i]
    position = SkyCoord(ra, dec, unit='deg', frame='fk5')
    in_image = position.contained_by(w0)
    if in_image:
        in_table.append(i)

# table of sources in QL catalogue but not in SE image
sources_not_in_se_image = table[:]
sources_not_in_se_image.remove_rows(in_table)
sources_not_in_se_image.write('sources_not_in_J100200+02300.csv', overwrite=True)

table_se_image = table[in_table]
table_se_image.write('sources_in_J100200+02300.csv', overwrite=True)

# flux density comparison
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib.ticker import AutoMinorLocator

def best_fit(x_data, y_data):

    #x_data  = list(x_data)
    #x_data.remove(max(x_data))
    #y_data = list(y_data)
    #y_data.remove(max(y_data))
    
    xbar = sum(x_data)/len(x_data)
    ybar = sum(y_data)/len(y_data)
    n = len(x_data) # or len(y_data)

    numer = sum([xi*yi for xi,yi in zip(x_data, y_data)]) - n * xbar * ybar
    denum = sum([xi**2 for xi in x_data]) - n * xbar**2

    b = numer / denum
    a = ybar - b * xbar

    print('best fit line:\ny = {:.2f} + {:.2f}x'.format(a, b))

    return a, b;



if (len(sys.argv) > 1 and len(sys.argv) < 3) or len(sys.argv) > 3:
    print(colored("\nWrong format used! Input format required:", "red"), colored('\n\tIn terminal','blue'),'\n\t\t$ python flux_comparison.py <input catalogue> <zeroth order FITS image> <first order FITS image>', colored('\n\n\tOr in iPython','blue'), colored('\n\t\t In[1]:','green'), ' %run flux_comparison.py <input catalogue> <zeroth order FITS image>\n')
    print(colored("Alternatively, run the script without invoking inputs and it will prompt you to enter each input separately!\n", "red"))
    print('\t\t $python make_tables.py\n')
    sys.exit()
try:
    quick_look_catalogue = sys.argv[1]
except:
    print('\n\n')
    os.system('ls *.csv')
    quick_look_catalogue = input("\nEnter a filename for the input catalogue:\n>")
try:
    photometry_catalogue = sys.argv[2]
except:
    print('\n\n')
    os.system('ls *tt0**fits')
    photometry_catalogue = input("\nEnter a filename for the zeroth order single ephoch image:\n>")
#quick_look_catalogue = 'sources_in_J100200+02300.csv'
#photometry_catalogue = 'J100200+023000_measured_tt0-tt1_Sv_and_alpha.csv'
target = photometry_catalogue.split('_')[0]

table_ql = Table.read(quick_look_catalogue)
table_phot = Table.read(photometry_catalogue)
table_ql['aperture_sum'] = table_phot['aperture_sum']
table_ql.sort('Peak_flux_source')

flux_ql = table_ql['Peak_flux_source']
flux_phot = table_ql['aperture_sum']

flux_phot *= 1000

a, b = best_fit(flux_ql, flux_phot)
yfit = [a+b*x for x in flux_ql]

fig, ax = plt.subplots(figsize=(8,6))
ax.plot(flux_ql, flux_phot, marker = '.', linestyle='', label='data')
ax.plot([flux_ql[0]-2, flux_ql[-1]+1], [flux_ql[0]-2 ,flux_ql[-1]+1], color='purple', label='y=x')
ax.plot(flux_ql, yfit, linestyle='--', color='grey', label='best fit')
#ax.set_title('$l=321.5$', loc='right')
ax.set_ylabel('Photometry Measured Flux Density ($mJy.beam^{-1}$)')
ax.set_xlabel('Quick Look Flux Density ($mJy.beam^{-1}$)')
ax.tick_params(axis='both', direction='in', which='both', bottom=True, left=True)
#ax.ticklabel_format(axis='both', style='sci',scilimits=(0,0))
ax.legend()
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
fig.show()
fig.savefig(f'quick_look_vs_photometry_flux_on_{target}.png', dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format=None, transparent=False, bbox_inches=None, pad_inches=0.1, metadata=None)

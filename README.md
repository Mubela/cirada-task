## TASKS

1. Accept a catalogue and two single epoch image file names (taylor coefficientis: tt0, tt1 images) as inputs
2. Measure the flux densities and spectral indices of the sources in the images at the positions of the sources listed in the catalogue
3. Write out the results to a CSV file
4. Identify sources appearing in the catalogue but not present in the image

5. Auxiliary step: Plot the fluxes in the quick look catalogue and those measure from the photometry on the single epoch image against each other

### This respository contains four scripts
```
make_tables.py
```
> Identifies the sources in the single epoch images (as each pair of single epoch images covers the same area of sky, it's run on only one). The script writes out two catalogues: the first contains the sources present in the single epoch image; the second contains the remaining sources, seen in the QL catalogue but in the single epoch images.

```
photometry_elliptical.py
```

> performs photometry on the two single epoch images (tt0 and tt1) using the output from make_tables.py that contains only the sources present in the single epoch images. It measures the flux densities in the zeroth order and first order images, and measures the spectral index using the ratio: I_alpha = I_1/I_0. These are written out to a CSV format file.

```
cirada.py
```

> performs the tasks implemented in make_tables.py and photometry_elliptical.py in a single script.

```
flux_comparison.py
```

> plots the fluxes in the quick look catalogue against the fluxes measured in the single epoch image

### Order of execution

Either:
```
> make_tables.py -> photometry_elliptical.py -> flux_comparison.py
```
Or:
```
> cirada.py -> flux_comparison.py
```

### Output table columns

'Source_name'      | Name as in the quick look catalogue  
'id'               | ID number as per the sources contained in the single epoch image used in the photometry  
'xcenter'          | x coordinate in pixels  
'ycenter'          | y coordinate in pixels  
'sky_center'       | sky coordinates in degrees (ra,dec)  
'aperture_sum'     | photometry measurement on the zeroth order image (tt0)  
'aperture_sum_1'   | photometry measurement on the first order image (tt1)  
'Spectral Index'   | Spectral index measurement, I_alpha = I_tt1/I_tt0  

### Dependencies

os  
sys  
tqdm  
termcolor  
numpy  
astropy  
radio_beam  
photutils  
matplotlib  

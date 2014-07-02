#! /usr/bin/env python

from astropy.table import Table, Column

photo_table = Table.read('photometry.fits')
M31_table = Table.read('tables/MasseyXL_M31_photo_err.fit')

for star in photo_table:
    if star['ID'] in M31_table['LGGS']:
        idx = M31_table['LGGS'].index(star['ID'])

        

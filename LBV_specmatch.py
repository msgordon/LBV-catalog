#! /usr/bin/env python
import numpy as np
import pyfits
from astropy.table import Table
import glob
import os

direc = 'M312013/M31B_Red'
outfile = 'M31B_Red.list'

table = Table.read('photometry_corr.tsv',format='ascii.tab')

RAphot = [x.split('+')[0][1:].split('.')[0] for x in table['ID']]

filelist = glob.glob(os.path.join(direc,'*.fits'))
RAspec = [os.path.basename(x).split('.')[0].split('-')[-1] for x in filelist]

f = open(outfile,'w')
for idx, phot in enumerate(RAphot):
    if phot in RAspec:
        idy = RAspec.index(phot)
        f.write('%s\t%s\n' % (table['ID'][idx], filelist[idy]))

f.close()


        



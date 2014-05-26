#! /usr/bin/env python
import argparse
from astropy.table import Table,Column
from star import Star
import numpy as np


def star_photometry(starList):

    # Sort list
    starList.sort(key=lambda x: x['ID'])

    # star table
    tList = Table([dict(star.get_photometry()) for star in starList])

    #t.add_column(c,index=0)

    t = Table()
    c = Column([star.ID for star in starList],name='ID')
    t.add_column(c,index=0)

    
    for col in tList.colnames:
        c = Column([float(row) if row else None for row in tList[col]],name=col)
        t.add_column(c)
    
    # get phot from colors
    Bcol = np.around([row['B_V'] + row['V'] if not row['B'] else row['B'] for row in t], decimals=2)
    t['B'] = Column(Bcol,name='B')

    Rcol = np.around([-(row['V_R'] - row['V']) if not row['R'] else row['R'] for row in t], decimals=2)
    t['R'] = Column(Rcol,name='R')

    Icol = [-(row['R_I'] - row['R']) if (row['R'] is not None and row['R_I'] is not None) else None for row in t]
    Icol = [np.around(row, decimals=2) if row is not None else None for row in Icol]
    t['I'] = Column(Icol,name='I')

    Ucol = [row['U_B'] + row['B'] if (row['B'] is not None and row['U_B'] is not None) else None for row in t]
    Ucol = [np.around(row, decimals=2) if row is not None else None for row in Ucol]
    t['U'] = Column(Ucol,name='U')


    # Convert to flux (Jy)
    ## http://ssc.spitzer.caltech.edu/warmmission/propkit/pet/magtojy/ref.html
    zp = {'U':1823,'B':4130,'V':3781,'R':2941,'I':2635}  #Jy
    

    ## http://www.stsci.edu/hst/nicmos/documents/handbooks/current_NEW/Appendix_B.14.3.html#329940
    for col in [t['U'],t['B'],t['V'],t['R'],t['I']]:
        c = [np.power(10.0,-val/2.5)*zp[col.name] if val is not None else None for val in col]
        t[col.name] = Column(c,name=col.name,description='Zeropoint: %i Jy' % zp[col.name],unit='Jy')


    # convert to F_lambda, erg/s/cm^2/micron
    wave = {'U':0.36,'B':0.44,'V':0.55,'R':0.71,'I':0.97}  #microns
    for col in [t['U'],t['B'],t['V'],t['R'],t['I']]:
        c = [3.0e-9 * val / (wave[col.name]**2) if val is not None else None for val in col]
        t[col.name] = Column(c,name=col.name,unit='erg*s^-1*cm^-2*micron^-1',description='Zeropoint: %i Jy, Eff_wave: %f' %(zp[col.name],wave[col.name]))

    print t

def main():
    parser = argparse.ArgumentParser(description="Returns photometric tables of catalog stars")

    parser.add_argument('catalog',type=str,help='JSON catalog of sources.')
    parser.add_argument('-WISE',type=str,required=True,help='Directory of WISE tables')
    parser.add_argument('-2MASS',type=str,required=True,help='Directory of 2MASS tables')

    args = parser.parse_args()

    # Get starlist
    print 'Loading JSON data from: %s' % args.catalog
    theList = Star.load(args.catalog)
    print '\tLoaded %i sources.' % len(theList)
    
    star_photometry(theList)




if __name__ == '__main__':
    main()

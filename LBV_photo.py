#! /usr/bin/env python
import argparse
from astropy.table import Table,Column,Row
from star import Star
import numpy as np
import atpy
import os

## http://ssc.spitzer.caltech.edu/warmmission/propkit/pet/magtojy/ref.html
zp = {'U':1823,'B':4130,'V':3781,'R':2941,'I':2635,'3.6':277.5,'4.5':179.5,'8.0':63.1,'J':1594,'H':1024,'K':666.7,'W1':309.540,'W2':171.787,'W3':31.674,'W4':8.363}  #Jy

wave = {'U':0.36,'B':0.44,'V':0.55,'R':0.71,'I':0.97,'3.6':3.55,'4.5':4.439,'8.0':7.872,'J':1.235,'H':1.662,'K':2.159,'W1':3.3526,'W2':4.6028,'W3':11.5608,'W4':22.0883}  #microns


###################
#  Read in table data
###################
def read_table_xml(filename,verbose=False):
    try:
        #print 'Reading table:  ' + filename
        t = atpy.Table(filename,verbose=False)
    except:
        #print 'No sources found in ' + filename
        return None

    if verbose:
        num = len(t)
        if num == 1:
            print 'Found ' + `num` + ' source in ' + filename
        else:
            print 'Found ' + `num` + ' sources in ' + filename
    return t


##############
# Find closest source to 'star' in 'table'
##############
def find_closest_source(star,table):
    minDist = np.Infinity
    Mindex = 0
    for i in xrange(0,len(table)):
        cDist = (star.RAd**2 + star.DECd**2) - (table['ra'][i]**2 - table['dec'][i]**2)
        if cDist < minDist:
            minDist = cDist
            Mindex = i

    return table[i]



def add_2MASS(table, starList, directory):

    # Sort list
    starList.sort(key=lambda x: x['ID'])

    JHK = []
    for star in starList:        
        # Get filename of source
        catfile = os.path.join(directory,star.ID+'.tbl')

        # Read table
        catTable = read_table_xml(catfile,verbose=False)

        # if found, find closest source
        if catTable is not None:
            source = find_closest_source(star,catTable)
            phot = {x:source[x] for x in ('j_m','h_m','k_m')}

        else:
            phot = {x:None for x in ('j_m','h_m','k_m')}

        JHK.append(phot)

    jCol = Column([x['j_m'] for x in JHK],name='J',unit='mag')
    hCol = Column([x['h_m'] for x in JHK],name='H',unit='mag')
    kCol = Column([x['k_m'] for x in JHK],name='K',unit='mag')

    # Convert to Jy
    #zp = {'J':1594,'H':1024,'K':666.7} # Jy
    #wave = {'J':1.235,'H':1.662,'K':2.159}  #microns

    print 'Appending J,H,K photometry from %s' % directory
    for x in (jCol,hCol,kCol):
        '''
        for idx,old in enumerate(table[x.name]):
            if old is None:
                table[x.name][idx] = x[idx]
            else:
                continue
        '''
        table[x.name] = x

        c = [np.power(10.0,-val/2.5)*zp[x.name] if val is not None else None for val in x]
        table.add_column(Column(c,name='F_%s_Jy' % x.name,description='Zeropoint: %i Jy' % zp[x.name],unit='Jy'))


    # convert to F_lambda, erg/s/cm^2/micron
        c = [3.0e-9 * val / (wave[x.name]**2) if val is not None else None for val in table['F_%s_Jy'%x.name]]
        table.add_column(Column(c,name='F_%.2f_um' % wave[x.name],unit='erg*s^-1*cm^-2*micron^-1',description='Zeropoint: %i Jy, Eff_wave: %f' %(zp[x.name],wave[x.name])))

        lc = [val * wave[x.name] if val is not None else None for val in c]
        table.add_column(Column(lc,name='lam_F_%.2f_um' % wave[x.name],unit='erg*s^-1*cm^-2',description='Zeropoint: %i Jy, Eff_wave: %f' %(zp[x.name],wave[x.name])))


    return table
    
def add_WISE(table, starList, directory):

    #http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html#WISEZMA
    # Sort list
    starList.sort(key=lambda x: x['ID'])

    Wbands = []
    for star in starList:        
        # Get filename of source
        catfile = os.path.join(directory,star.ID+'.tbl')

        # Read table
        catTable = read_table_xml(catfile,verbose=False)

        # if found, find closest source
        if catTable is not None:
            source = find_closest_source(star,catTable)
            phot = {x:source[x] for x in ('w1mag','w2mag','w3mag','w4mag')}

        else:
            phot = {x:None for x in ('w1mag','w2mag','w3mag','w4mag')}

        Wbands.append(phot)

    w1Col = Column([x['w1mag'] for x in Wbands],name='W1',unit='mag')
    w2Col = Column([x['w2mag'] for x in Wbands],name='W2',unit='mag')
    w3Col = Column([x['w3mag'] for x in Wbands],name='W3',unit='mag')
    w4Col = Column([x['w4mag'] for x in Wbands],name='W4',unit='mag')

    # Convert to Jy
    #zp = {'W1':309.540,'W2':171.787,'W3':31.674,'W4':8.363}  #Jy
    #wave = {'W1':3.3526,'W2':4.6028,'W3':11.5608,'W4':22.0883}  #microns
    
    print 'Appending 3-22 um photometry from %s' % directory
    for x in (w1Col,w2Col,w3Col,w4Col):
        table[x.name] = x

        c = [np.power(10.0,-val/2.5)*zp[x.name] if val is not None else None for val in x]
        table.add_column(Column(c,name='F_%s_Jy' % x.name,description='Zeropoint: %i Jy' % zp[x.name],unit='Jy'))


    # convert to F_lambda, erg/s/cm^2/micron
        c = [3.0e-9 * val / (wave[x.name]**2) if val is not None else None for val in table['F_%s_Jy'%x.name]]
        table.add_column(Column(c,name='F_%.2f_um' % wave[x.name],unit='erg*s^-1*cm^-2*micron^-1',description='Zeropoint: %i Jy, Eff_wave: %f' %(zp[x.name],wave[x.name])))

        lc = [val * wave[x.name] if val is not None else None for val in c]
        table.add_column(Column(lc,name='lam_F_%.2f_um' % wave[x.name],unit='erg*s^-1*cm^-2',description='Zeropoint: %i Jy, Eff_wave: %f' %(zp[x.name],wave[x.name])))


    return table

            
    


def star_photometry(starList):

    # Sort list
    starList.sort(key=lambda x: x['ID'])

    # star table
    tList = Table([dict(star.get_photometry()) for star in starList])

    #t.add_column(c,index=0)

    t = Table()
    c = Column([star.ID for star in starList],name='ID')
    t.add_column(c,index=0)
    c = Column([star.get_gal_name() for star in starList],name='Gal')
    t.add_column(c,index=1)

    for col in tList.colnames:
        if col in ['__3_6_','__4_5_','__8_0_']:
            name = '.'.join(col.strip('_').split('_'))
        else:
            name = col
        c = Column([float(row) if row else None for row in tList[col]],name=name)
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
    #zp = {'U':1823,'B':4130,'V':3781,'R':2941,'I':2635,'3.6':277.5,'4.5':179.5,'8.0':63.1}  #Jy
    

    ## http://www.stsci.edu/hst/nicmos/documents/handbooks/current_NEW/Appendix_B.14.3.html#329940
    for col in ['U','B','V','R','I','3.6','4.5','8.0']:
        c = [np.power(10.0,-val/2.5)*zp[col] if val is not None else None for val in t[col]]
        t.add_column(Column(c,name='F_%s_Jy' % col,description='Zeropoint: %i Jy' % zp[col],unit='Jy'))


    # convert to F_lambda, erg/s/cm^2/micron
    #wave = {'U':0.36,'B':0.44,'V':0.55,'R':0.71,'I':0.97,'3.6':3.55,'4.5':4.439,'8.0':7.872}  #microns
    for col in ['U','B','V','R','I','3.6','4.5','8.0']:
        c = [3.0e-9 * val / (wave[col]**2) if val is not None else None for val in t['F_%s_Jy'%col]]
        t.add_column(Column(c,name='F_%.2f_um' % wave[col],unit='erg*s^-1*cm^-2*micron^-1',description='Zeropoint: %i Jy, Eff_wave: %f' %(zp[col],wave[col])))

        lc = [val * wave[col] if val is not None else None for val in c]
        t.add_column(Column(lc,name='lam_F_%.2f_um' % wave[col],unit='erg*s^-1*cm^-2',description='Zeropoint: %i Jy, Eff_wave: %f' %(zp[col],wave[col])))


    return t




def main():
    parser = argparse.ArgumentParser(description="Returns photometric tables of catalog stars")

    parser.add_argument('catalog',type=str,help='JSON catalog of sources.')
    parser.add_argument('-WISE',type=str,required=True,help='Directory of WISE tables')
    parser.add_argument('-2MASS',type=str,dest='MASS',required=True,help='Directory of 2MASS tables')

    args = parser.parse_args()

    # Get starlist
    print 'Loading JSON data from: %s' % args.catalog
    theList = Star.load(args.catalog)
    print '\tLoaded %i sources.' % len(theList)


    print
    print 'Getting photometry from catalogs...'
    t = star_photometry(theList)
    
    print 'Locating sources in %s' % args.MASS
    t = add_2MASS(t,theList,args.MASS)

    print 'Locating sources in %s' % args.WISE
    t = add_WISE(t,theList,args.WISE)
    
    print

    #t.write('photometry.tsv',format='ascii.tab')
    #exit()

    outfile = 'photometry.fits'
    
    print
    

    colnames = [x for x in t.colnames if 'lam' in x]
    for col in colnames:
        t[col] = [99.99 if ((x is None) or (x is 'None')) else x for x in t[col]]


    # Output tab file for Roberta
    t.write('photometry.tsv', format='ascii.tab')
    exit()
        
    newCols = []
    for col in colnames:
        c = Column([np.float(x) for x in t[col]],name=col,dtype=np.float)
        newCols.append(c)
    #print colnames
    eTable = Table()
    eTable.add_column(t['ID'])
    #eTable.add_column(t['lam_F_0.55_um'])
    #c = Column([str(x) for x in t['ID']],name='ID',dtype=str)
    #eTable.add_column(c)
    #eTable.add_columns([t[x] for x in colnames])
    eTable.add_columns(newCols)

    #print eTable['ID'].dtype
    #exit()
               
    #print eTable.colnames
    #for col in eTable.colnames:
    #    print eTable[col]

    #print eTable
    print 'Writing table to %s' % outfile
    eTable.write(outfile)




if __name__ == '__main__':
    main()

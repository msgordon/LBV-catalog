#! /usr/bin/env python
import argparse
from astropy.table import Table,Column,Row,vstack
from star import Star
import numpy as np
import atpy
import os
from astropy.coordinates import SkyCoord,match_coordinates_sky
import astropy.units as u

## http://ssc.spitzer.caltech.edu/warmmission/propkit/pet/magtojy/ref.html
zp = {'U':1823,'B':4130,'V':3781,'R':2941,'I':2635,'3.6':277.5,'4.5':179.5,'5.8':116.6,'8.0':63.1,'J':1594,'H':1024,'K':666.7,'W1':309.540,'W2':171.787,'W3':31.674,'W4':8.363}  #Jy

wave = {'U':0.36,'B':0.44,'V':0.55,'R':0.71,'I':0.97,'3.6':3.55,'4.5':4.439,'5.8':5.731,'8.0':7.872,'J':1.235,'H':1.662,'K':2.159,'W1':3.3526,'W2':4.6028,'W3':11.5608,'W4':22.0883}  #microns


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
def find_closest_source2(star,table,rad=3):
    star_c = SkyCoord(ra=star.RAd*u.deg,dec=star.DECd*u.deg,frame='icrs')

    table_c = []
    for row in table:
        table_c.append((row['ra'],row['dec']))
    table_ra,table_dec = zip(*table_c)
    table_c = SkyCoord(ra=table_ra*u.deg,dec=table_dec*u.deg,frame='icrs')

    idx,sep2d,dist3d = match_coordinates_sky(star_c,table_c,
                                             storekdtree=u'_kdtree_sky')

    return table[int(idx)]
    if sep2d.is_within_bounds(upper=rad*u.arcsec):
        return table[int(idx)]
    else:
        return None

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
    Mblank = []
    for star in starList:        
        # Get filename of source
        catfile = os.path.join(directory,star.ID+'.tbl')

        # Read table
        catTable = read_table_xml(catfile,verbose=False)

        # if found, find closest source
        if catTable is not None:
            source = find_closest_source2(star,catTable,rad=4)
            phot = {x:source[x] for x in ('j_m','h_m','k_m')}

        else:
            phot = {x:None for x in ('j_m','h_m','k_m')}
            Mblank.append(star.ID)

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

    print 'Failed to find 2MASS photometry for %i sources' % len(Mblank)

    return table
    
def add_WISE(table, starList, directory):

    #http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html#WISEZMA
    # Sort list
    starList.sort(key=lambda x: x['ID'])

    Wbands = []
    Wblank = []
    for star in starList:        
        # Get filename of source
        catfile = os.path.join(directory,star.ID+'.tbl')

        # Read table
        catTable = read_table_xml(catfile,verbose=False)

        # if found, find closest source
        if catTable is not None:
            source = find_closest_source2(star,catTable,rad=5)

            try:
                phot = {x:source[x] for x in ('w1mag','w2mag','w3mag','w4mag')}
            except:
                phot = {x:source[y] for x,y in zip(('w1mag','w2mag','w3mag','w4mag'),('w1mpro','w2mpro','w3mpro','w4mpro'))}

        else:
            phot = {x:None for x in ('w1mag','w2mag','w3mag','w4mag')}
            Wblank.append(star.ID)

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


    print 'Failed to find WISE photometry for %i sources' % len(Wblank)
        
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
        if col in ['__3_6_','__4_5_','__5_8_','__8_0_']:
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
    for col in ['U','B','V','R','I','3.6','4.5','5.8','8.0']:
        c = [np.power(10.0,-val/2.5)*zp[col] if val is not None else None for val in t[col]]
        t.add_column(Column(c,name='F_%s_Jy' % col,description='Zeropoint: %i Jy' % zp[col],unit='Jy'))


    # convert to F_lambda, erg/s/cm^2/micron
    #wave = {'U':0.36,'B':0.44,'V':0.55,'R':0.71,'I':0.97,'3.6':3.55,'4.5':4.439,'8.0':7.872}  #microns
    for col in ['U','B','V','R','I','3.6','4.5','5.8','8.0']:
        c = [3.0e-9 * val / (wave[col]**2) if val is not None else None for val in t['F_%s_Jy'%col]]
        t.add_column(Column(c,name='F_%.2f_um' % wave[col],unit='erg*s^-1*cm^-2*micron^-1',description='Zeropoint: %i Jy, Eff_wave: %f' %(zp[col],wave[col])))

        lc = [val * wave[col] if val is not None else None for val in c]
        t.add_column(Column(lc,name='lam_F_%.2f_um' % wave[col],unit='erg*s^-1*cm^-2',description='Zeropoint: %i Jy, Eff_wave: %f' %(zp[col],wave[col])))


    return t


def photo_corr(phot_table):
    # replace broke ass values
    # J013337.00+303637.5 -> J013337.04+303637.6
    # J013415.38+302816.3 -> J013415.42+302816.4 
    id1 = np.where(phot_table['ID'] == 'J013337.00+303637.5')
    id2 = np.where(phot_table['ID'] == 'J013415.38+302816.3')
    phot_table['ID'][id1] = 'J013337.04+303637.6'
    phot_table['ID'][id2] = 'J013415.42+302816.4'

    M31_table = Table.read('tables/MasseyXL_M31_photo_err.fit')
    M33_table = Table.read('tables/MasseyXL_M33_photo_err.fit')

    idx = [np.where(M31_table['LGGS'] == star['ID']) for star in phot_table if star['ID'] in M31_table['LGGS']]
    rows31 = M31_table[idx]
    idx = [np.where(M33_table['LGGS'] == star['ID']) for star in phot_table if star['ID'] in M33_table['LGGS']]
    rows33 = M33_table[idx]

    rows = vstack([rows31, rows33])

    e_Vmag = Column(data=np.ndarray.flatten(rows['e_Vmag']),name='e_Vmag',dtype=np.float)#**2-rows['e_Vmag']**2))
    e_Bmag = Column(data=np.ndarray.flatten(rows['e_B-V']),name='e_Bmag',dtype=np.float)
    e_Umag = Column(data=np.ndarray.flatten(rows['e_U-B']),name='e_Umag',dtype=np.float)
    e_Rmag = Column(data=np.ndarray.flatten(rows['e_V-R']),name='e_Rmag',dtype=np.float)
    e_Imag = Column(data=np.ndarray.flatten(rows['e_R-I']),name='e_Imag',dtype=np.float)

    e_Vflux = Column(data=np.log(10.)/2.5*np.ndarray.flatten(np.array([x*y for x,y in zip(phot_table['F_V_Jy'],e_Vmag)])),name='e_F_V_Jy')
    e_Bflux = Column(data=np.log(10.)/2.5*np.ndarray.flatten(np.array([x*y for x,y in zip(phot_table['F_B_Jy'],e_Bmag)])),name='e_F_B_Jy')
    e_Uflux = Column(data=np.log(10.)/2.5*np.ndarray.flatten(np.array([x*y for x,y in zip(phot_table['F_U_Jy'],e_Umag)])),name='e_F_U_Jy')
    e_Rflux = Column(data=np.log(10.)/2.5*np.ndarray.flatten(np.array([x*y for x,y in zip(phot_table['F_R_Jy'],e_Rmag)])),name='e_F_R_Jy')
    e_Iflux = Column(data=np.log(10.)/2.5*np.ndarray.flatten(np.array([x*y for x,y in zip(phot_table['F_I_Jy'],e_Imag)])),name='e_F_I_Jy')


    e_Vlamflux = Column(data=np.ndarray.flatten(np.array([x for x in e_Vflux]))*3.0e-9/0.55**2,name='e_F_0.55_um')
    e_Blamflux = Column(data=np.ndarray.flatten(np.array([x for x in e_Bflux]))*3.0e-9/0.44**2,name='e_F_0.44_um')
    e_Ulamflux = Column(data=np.ndarray.flatten(np.array([x for x in e_Uflux]))*3.0e-9/0.36**2,name='e_F_0.36_um')
    e_Rlamflux = Column(data=np.ndarray.flatten(np.array([x for x in e_Rflux]))*3.0e-9/0.71**2,name='e_F_0.71_um')
    e_Ilamflux = Column(data=np.ndarray.flatten(np.array([x for x in e_Iflux]))*3.0e-9/0.97**2,name='e_F_0.97_um')

    e_Ve = Column(data=np.ndarray.flatten(np.array([x for x in e_Vlamflux]))*0.55,name='e_lam_F_0.55_um')
    e_Be = Column(data=np.ndarray.flatten(np.array([x for x in e_Blamflux]))*0.44,name='e_lam_F_0.44_um')
    e_Ue = Column(data=np.ndarray.flatten(np.array([x for x in e_Ulamflux]))*0.36,name='e_lam_F_0.36_um')
    e_Re = Column(data=np.ndarray.flatten(np.array([x for x in e_Rlamflux]))*0.71,name='e_lam_F_0.71_um')
    e_Ie = Column(data=np.ndarray.flatten(np.array([x for x in e_Ilamflux]))*0.97,name='e_lam_F_0.97_um')


    phot_table.add_columns([e_Umag,e_Bmag,e_Vmag,e_Rmag,e_Imag,e_Uflux,e_Bflux,e_Vflux,e_Rflux,e_Iflux,e_Ulamflux,e_Blamflux,e_Vlamflux,e_Rlamflux,e_Ilamflux,e_Ue,e_Be,e_Ve,e_Re,e_Ie])

    return phot_table




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

    outfile = 'photometry_ZOMG'

    colnames = [x for x in t.colnames if 'lam' in x]
    for col in colnames:
        t[col] = [99.99 if ((x is None) or (x is 'None')) else x for x in t[col]]

    #for Roberta
    rTable = Table()
    rTable.add_columns([t[x] for x in ['ID','Gal']])
    for col in ['U','B','V','R','I','J','H','K',
                '3.6','4.5','5.8','8.0',
                'W1','W2','W3','W4',
                'F_U_Jy','F_B_Jy','F_V_Jy','F_R_Jy','F_I_Jy',
                'F_J_Jy','F_H_Jy','F_K_Jy',
                'F_3.6_Jy','F_4.5_Jy','F_5.8_Jy','F_8.0_Jy',
                'F_W1_Jy','F_W2_Jy','F_W3_Jy','F_W4_Jy',
                'F_0.36_um','lam_F_0.36_um','F_0.44_um','lam_F_0.44_um','F_0.55_um','lam_F_0.55_um','F_0.71_um','lam_F_0.71_um','F_0.97_um','lam_F_0.97_um',
                'F_1.24_um','lam_F_1.24_um','F_1.66_um','lam_F_1.66_um','F_2.16_um','lam_F_2.16_um',
                'F_3.55_um','lam_F_3.55_um','F_4.44_um','lam_F_4.44_um','F_5.73_um','lam_F_5.73_um','F_7.87_um','lam_F_7.87_um',
            'F_3.35_um','lam_F_3.35_um','F_4.60_um','lam_F_4.60_um','F_11.56_um','lam_F_11.56_um','F_22.09_um','lam_F_22.09_um']:
        c = Column([np.float(x) if x else 99.99 for x in t[col]],name=col,dtype=np.float)
        rTable.add_column(c)

    rTable = photo_corr(rTable)

    rTable.write(outfile+'.tsv',format='ascii.tab')
    rTable.write(outfile+'.fits')
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
    outfile = 'photometry.fits'
    print 'Writing table to %s' % outfile
    eTable.write(outfile)




if __name__ == '__main__':
    main()

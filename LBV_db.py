#! /usr/bin/env python

import numpy as np
import logging
import collections
from collections import OrderedDict as dict
from star import Star
import sys
#### Requires python2.7


##################
# LBV_DB class holds database of objects
#  Includes data read routines
#  Includes matching scheme

class LBV_DB:

    # Class members
    FILE_dict = dict()
    INCLUDE_dict = dict()
    DATA_dict = dict()  # Each entry will be a dictionary mapped to a set
    LOOKUP_dict = dict()
    NUM_STARS = 0
    
    def __call__(self,*args):
        return self


    # Read datasets, and return as dictionary
    # 'include' flags this catalog for inclusion
    def __init__(self):        
        WeisFile = "./tables/Weis.dat"
        WeisDat = np.genfromtxt(WeisFile,names=['Name','RA','DEC','U','B','V','R','I','priority'],dtype=['a10','f8','f8','f8','f8','f8','f8','f8','i4'],autostrip=True)
        self.DATA_dict['w'] = WeisDat
        self.FILE_dict['w'] = WeisFile
        self.INCLUDE_dict['w'] = False
        self.NUM_STARS += len(WeisDat)

        ValeevFile = "./tables/Valeev_altered2.dat"
        ValeevDat = np.genfromtxt(ValeevFile,skip_header=1,delimiter=';',names=['Name','RA','DEC','V','B_V','s','R','comment'],dtype=['a6','f8','f8','f8','f8','f8','f8','a100'],autostrip=True)
        self.DATA_dict['v'] = ValeevDat
        self.FILE_dict['v'] = ValeevFile
        self.INCLUDE_dict['v'] = False
        self.NUM_STARS += len(ValeevDat)
        
        HumphreysFile = "./tables/Humphreys.dat"
        HumphreysDat = np.genfromtxt(HumphreysFile,delimiter=';',names=['Gal','Name','RA','DEC','Spec','comment'],dtype=['a3','a10','f8','f8','a20','a100'],autostrip=True)
        self.DATA_dict['h'] = HumphreysDat
        self.FILE_dict['h'] = HumphreysFile
        self.INCLUDE_dict['h'] = False
        self.NUM_STARS += len(HumphreysDat)

        
        MasseyFile = "./tables/Massey17_final.tsv"
        MasseyDat = np.genfromtxt(MasseyFile,delimiter=';',names=['LGGS','Gal','Name','V','B_V','U_B','SpType','CrossID','Ref','RA','DEC'],dtype=['a20','a4','a20','a20','a20','a20','a20','a20','a20','f8','f8'],autostrip=True)
        self.DATA_dict['m'] = MasseyDat
        self.FILE_dict['m'] = MasseyFile
        self.INCLUDE_dict['m'] = False
        self.NUM_STARS += len(MasseyDat)

        MasseyLongFile = "./tables/MasseyLong.tsv"
        MasseyLongDat = np.genfromtxt(MasseyLongFile,skip_header=1,delimiter=';',names=['Gal','LGGS','n_LGGS','HAmag','HA','O3','V','B_V','U_B','Q','Type','Name','Ref','RA','DEC'],dtype=['a8','a20','a1','f8','f8','f8','f8','f8','f8','f8','a20','a20','a5','f8','f8'],autostrip=True)
        self.DATA_dict['mL'] = MasseyLongDat
        self.FILE_dict['mL'] = MasseyLongFile
        self.INCLUDE_dict['mL'] = False
        self.NUM_STARS += len(MasseyLongDat)
        

        McQuinnNV = "./tables/McQuinn_nonvariable.txt"
        McQuinnNVDat = np.genfromtxt(McQuinnNV,names=['Gal','ID','RA','DEC'],dtype=['a20','a20','f8','f8'],usecols=(0,1,2,3),delimiter=';',autostrip=True)
        self.DATA_dict['mqnv'] = McQuinnNVDat
        self.FILE_dict['mqnv'] = McQuinnNV
        self.INCLUDE_dict['mqnv'] = False
        for idx,val in enumerate(McQuinnNVDat['ID']):
            McQuinnNVDat['RA'][idx] = np.float(val[1:10])
            McQuinnNVDat['DEC'][idx] = np.float(val[11:19])
        self.NUM_STARS += len(McQuinnNVDat)

        McQuinnV = "./tables/McQuinn_variable.txt"
        McQuinnVDat = np.genfromtxt(McQuinnV,names=['Gal','ID','RA','DEC'],dtype=['a20','a20','f8','f8'],usecols=(0,1,2,3),delimiter=';',autostrip=True)
        self.DATA_dict['mqv'] = McQuinnVDat
        self.FILE_dict['mqv'] = McQuinnV
        self.INCLUDE_dict['mqv'] = False
        for idx,val in enumerate(McQuinnVDat['ID']):
            McQuinnVDat['RA'][idx] = np.float(val[1:10])
            McQuinnVDat['DEC'][idx] = np.float(val[11:19])
        self.NUM_STARS += len(McQuinnVDat)

        '''
        Drout31 = "./tables/Drout_M31.txt"
        Drout31data = np.genfromtxt(Drout31,delimiter=';',names=['RAJ','DEJ','recno','Name','rank','RA','DEC','Vobs','fVobs','r','Vexp','VomVe','V','B_V','srank','LGGS','SimbadName'],dtype=['a20','a20','a20','a20','a20','f8','f8','f8','a1','f8','f8','f8','f8','f8','a2','a4','a30'],autostrip=True)
        self.DATA_dict['d31'] = Drout31data
        self.FILE_dict['d31'] = Drout31
        self.INCLUDE_dict['d31'] = True
        self.NUM_STARS += len(Drout31data)

        
        Drout33_r = "./tables/Drout_M33_red.txt"
        Drout33_r_data = np.genfromtxt(Drout33_r,delimiter=';',names=['Name','RA','DEC','Vobs','fVobs','r','Vexp','VomVe','V','B_V','fBband','VmR','K','rank'],dtype=['a19','f8','f8','f8','a1','f8','f8','f8','f8','f8','a1','f8','f8','a1'],autostrip=True)
        self.DATA_dict['d33_r'] = Drout33_r_data
        self.FILE_dict['d33_r'] = Drout33_r
        self.INCLUDE_dict['d33_r'] = True
        self.NUM_STARS += len(Drout33_r_data)

        
        Drout33_y = "./tables/Drout_M33_yellow.txt"
        Drout33_y_data = np.genfromtxt(Drout33_y,delimiter=';',names=['Name','RA','DEC','Vobs','fVobs','r','Vexp','VomVe','V','B_V','ew','rank'],dtype=['a19','f8','f8','f8','a1','f8','f8','f8','f8','f8','f8','a2'],autostrip=True)
        self.DATA_dict['d33_y'] = Drout33_y_data
        self.FILE_dict['d33_y'] = Drout33_y
        self.INCLUDE_dict['d33_y'] = True
        self.NUM_STARS += len(Drout33_y_data)
        '''

        Drout31 = "./tables/Drout_M31_SG.tsv"
        Drout31data = np.genfromtxt(Drout31,delimiter=';',names=['recno','Name','RA','DEC','rank','V','fV','logT','flogT','logL','flogL','comment'],dtype=['i8','a20','f8','f8','i8','f8','a1','f8','a1','f8','a1','a31'],autostrip=True)
        self.DATA_dict['d31'] = Drout31data
        self.FILE_dict['d31'] = Drout31
        self.INCLUDE_dict['d31'] = True
        self.NUM_STARS += len(Drout31data)
        
        Drout33_r = "./tables/Drout_M33_redSG.txt"
        Drout33_r_data = np.genfromtxt(Drout33_r,delimiter=';',names=['Name','RA','DEC','rank','Teff','Lum','Var'],dtype=['a19','f8','f8','a1','f8','f8','a5'],autostrip=True)
        self.DATA_dict['d33_r'] = Drout33_r_data
        self.FILE_dict['d33_r'] = Drout33_r
        self.INCLUDE_dict['d33_r'] = True
        self.NUM_STARS += len(Drout33_r_data)
        
        Drout33_y = "./tables/Drout_M33_yellowSG.txt"
        Drout33_y_data = np.genfromtxt(Drout33_y,delimiter=';',names=['Name','RA','DEC','rank','Teff','Lum','Var'],dtype=['a19','f8','f8','a1','f8','f8','a5'],autostrip=True)
        self.DATA_dict['d33_y'] = Drout33_y_data
        self.FILE_dict['d33_y'] = Drout33_y
        self.INCLUDE_dict['d33_y'] = True
        self.NUM_STARS += len(Drout33_y_data)

        
        # Build lookup dictionary
        print 'Building lookup dictionary...'
        
        MasseyXLongFile = "./tables/MASSEYXLFINAL.dat"
        MasseyXLongDat = np.genfromtxt(MasseyXLongFile,skip_header=1,delimiter=';',names=['Gal','LGGS','n_LGGS','Name','V','B_V','U_B','V_R','R_I','Type','rType','RA','DEC'],dtype=['a5','a20','a1','a20','f8','f8','f8','f8','f8','a10','a5','f8','f8'],autostrip=True)
        for row in MasseyXLongDat:
            self.LOOKUP_dict[row['LGGS']] = row
        print '\t%i sources available for search' % len(self.LOOKUP_dict)
        

    # Match objects on 'RA' and 'DEC'.  Return match or None.
    @staticmethod
    def find_match(this, matchSet):
        for that in matchSet:
            if ((np.abs(this['RA'] - that['RA']) < 0.3) and (np.abs(this['DEC'] - that['DEC']) < 0.3)):
                return that

        return None


    # Update progress
    @staticmethod
    def update_progress(progress,toPrint=''):
        num_dash = int(progress/10)
        num_spaces = 10-num_dash
        print '\r    [{0}{1}] {2}%     {3}'.format('#'*num_dash,' '*num_spaces, int(progress),toPrint),
        sys.stdout.flush()

    # Match wrapper.  Returns list of stars.
    def get_matches(self):
        theList = set()   # Master list to hold all stars (matched and unique)
       
        # Loop through catalogs
        for this_cat, this_cat_dict in self.DATA_dict.iteritems():
            # Skip if this catalog is not meant to be included
            if self.INCLUDE_dict[this_cat] == False:
                continue

            logging.info('Matching from %s' % (self.FILE_dict[this_cat]))
            total_matches = 0

            # Loop through each star
            for idx, this_star in enumerate(this_cat_dict):
                these_matches = []  # Holds all matches to this_star
                these_matches.append((this_cat,this_star))
                 # Update progress
                LBV_DB.update_progress(idx*100./len(this_cat_dict),this_star['Name'])

                # Loop through other catalogs
                for that_cat, that_cat_dict in self.DATA_dict.iteritems():
                    # Skip this catalog
                    if this_cat == that_cat:
                        continue
                     
                    # Find match from each catalog
                    that_star = LBV_DB.find_match(this_star,that_cat_dict)
                    if that_star is not None:
                        these_matches.append((that_cat,that_star))
                                            
                if len(these_matches) > 1:  # Matches were found
                    total_matches += 1


                # For each catalog, assign a dictionary to each data item
                starDict = dict(these_matches)
                for catalog,data in starDict.iteritems():
                    starDict[catalog] = dict(zip(data.dtype.names,data))
                    
                    # add color to Weis
                    if catalog == 'w':
                        starDict['w']['B_V'] = starDict['w']['B']-starDict['w']['V']
                        starDict['w']['U_B'] = starDict['w']['U']-starDict['w']['B']
                    
                    # add fields for filename, include
                    starDict[catalog]['file'] = self.FILE_dict[catalog]
                    starDict[catalog]['include'] = self.INCLUDE_dict[catalog]

    
                # Convert to Star and add to master list
                star = Star(sdict=starDict)
                theList.add(star)  #If star is already in this Set, no action
            

                                
            logging.info('\n\tFound %i/%i matches.' % (total_matches, len(this_cat_dict)))
        logging.info('Found %i total sources.' % (len(theList)))

        logging.info('Checking against lookup dictionary...')
        match_count = 0
        for star in theList:
            # Check if match in mXL
            if star.ID in self.LOOKUP_dict:
                data = self.LOOKUP_dict[star.ID]
                star.__dict__['mXL'] = dict(zip(data.dtype.names,data))
                star['mXL']['file'] = 'tables/MASSEYXLFINAL.dat'
                star['mXL']['include'] = False
                match_count += 1
                             
        logging.info('\t%i/%i sources found in lookup dictionary.'%(match_count,len(theList)))
        return list(theList)


#! /usr/bin/env python
import numpy as np
#import logging
#import collections
from collections import OrderedDict as dict
from star import Star
import sys
import argparse
from astropy.table import Table,Column
from astropy.coordinates import match_coordinates_sky,SkyCoord
from astropy import units as u
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
    COORD_dict = dict() # SkyCoords of each object
    LOOKUP_dict = dict()
    NUM_STARS = 0
    
    def __call__(self,*args):
        return self


    # Read datasets, and return as dictionary
    # 'include' flags this catalog for inclusion
    def __init__(self,catList,masterFile):
        for filename in catList:
            print 'Reading table from %s' % filename
            table = Table.read(filename)
            title = table.meta['TITLE']

            if 'RAd' not in table.colnames:
                RAd,DECd = zip(*[Star.sex2deg(ra,dec) for ra,dec in zip(table['RA'],table['DEC'])])
                c = Column(RAd,name='RAd',dtype=np.float)
                table.add_column(c)
                c = Column(DECd,name='DECd',dtype=np.float)
                table.add_column(c)
            
            self.masterFile = masterFile
            self.DATA_dict[title] = table
            self.FILE_dict[title] = filename
            self.INCLUDE_dict[title] = table.meta['INCLUDE']
            self.NUM_STARS += len(table)

            # Add skycoords to each dictionary
            for cat,cat_dict in self.DATA_dict.iteritems():
                self.COORD_dict[cat] = SkyCoord(ra=cat_dict['RAd']*u.deg,
                                                dec=cat_dict['DECd']*u.deg,
                                                frame='icrs')
        


        print

        # Build lookup dictionary
        print 'Building lookup dictionary from %s' % masterFile
        masterTable = Table.read(masterFile)
        '''
        if 'RAd' not in masterTable.colnames:
                RAd,DECd = zip(*[Star.sex2deg(ra,dec) for ra,dec in zip(masterTable['RA'],masterTable['DEC'])])
                c = Column(RAd,name='RAd',dtype=np.float)
                masterTable.add_column(c)
                c = Column(DECd,name='DECd',dtype=np.float)
                masterTable.add_column(c)
        '''     
        self.LOOKUP_dict = dict(zip(masterTable['LGGS'],masterTable))

        #for row in MasseyXLongDat:
        #    self.LOOKUP_dict[row['LGGS']] = row
        print '\t%i sources available for search' % len(self.LOOKUP_dict)
        print


    @staticmethod
    def find_match2(this,matchCoords,matchSet,rad=3):
        # First, check if IDs match
        id_star = this['ID'] if 'ID' in this.columns else this['Name']
        ids_matchSet = matchSet['ID'] if 'ID' in matchSet.columns else matchSet['Name']
        try:
            idx = list(ids_matchSet).index(id_star)
        except:
            pass #just continue
        else:
            return matchSet[idx]
            
        
        c = SkyCoord(ra=this['RAd']*u.deg,
                     dec=this['DECd']*u.deg,
                     frame='icrs')
        #catalog = SkyCoord(ra=matchSet['RAd']*u.degree,dec=matchSet['DECd']*u.degree,frame='icrs')

        idx,sep2d,dist3d = match_coordinates_sky(c,matchCoords,
                                                 storekdtree=u'_kdtree_sky')

        if sep2d.is_within_bounds(upper=rad*u.arcsec):
            return matchSet[int(idx)]
        else:
            return None



        
    # Match objects on 'RA' and 'DEC'.  Return match or None.
    @staticmethod
    def find_match(this, matchSet,degree=False):
        if degree:
            if 'RAd' not in this.colnames:
                RAt,DECt = Star.sex2deg(this['RA'],this['DEC'])
            else:
                RAt = this['RAd']
                DECt = this['DECd']

            if 'RAd' not in matchSet.colnames:
                matchcoords = [Star.sex2deg(ra,dec) for ra,dec in zip(matchSet['RA'],matchSet['DEC'])]
                dist = [np.sqrt((RAt-ra)**2 + (DECt-dec)**2) for ra,dec in matchcoords]
            else:
                dist = [np.sqrt((RAt - that['RAd'])**2 + (DECt - that['DECd'])**2) for that in matchSet]
            if min(dist) < 0.3:
                return matchSet[np.argmin(dist)]

        else:
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

        sys.stdout.write("\033[K") # Clear to the end of line
        sys.stdout.flush()

    # Match wrapper.  Returns list of stars.
    def get_matches(self):
        theList = set()   # Master list to hold all stars (matched and unique)
       
        # Loop through catalogs
        for this_cat, this_cat_dict in self.DATA_dict.iteritems():
            # Skip if this catalog is not meant to be included
            if self.INCLUDE_dict[this_cat] == False:
                continue

            print 'Matching from %s' % (self.FILE_dict[this_cat])
            total_matches = 0

            # Loop through each star
            for idx, this_star in enumerate(this_cat_dict):
                these_matches = []  # Holds all matches to this_star
                these_matches.append((this_cat,this_star))
                 # Update progress
                LBV_DB.update_progress(idx*100./len(this_cat_dict),this_star['Name'])
                if idx == len(this_cat_dict)-1:
                    LBV_DB.update_progress(100,'Complete')

                # Loop through other catalogs
                for that_cat, that_cat_dict in self.DATA_dict.iteritems():
                    # Skip this catalog
                    if this_cat == that_cat:
                        continue
                     
                    # Find match from each catalog
                    ###########################
                    that_star = LBV_DB.find_match2(this_star,
                                                   self.COORD_dict[that_cat],
                                                   that_cat_dict)
                    #if that_cat in ['jm','mqnv','mqv','mcq']:
                    #    that_star = LBV_DB.find_match(this_star,that_cat_dict,degree=True)
                    #else:
                    #    that_star = LBV_DB.find_match(this_star,that_cat_dict)
                    #that_star = LBV_DB.find_match(this_star,that_cat_dict)
                    if that_star is not None:
                        these_matches.append((that_cat,that_star))
                                            
                if len(these_matches) > 1:  # Matches were found
                    total_matches += 1


                # For each catalog, assign a dictionary to each data item
                starDict = dict(these_matches)
                for catalog,data in starDict.iteritems():
                    starDict[catalog] = dict(zip(data.dtype.names,data))
                    
                                       
                    # add fields for filename, include
                    starDict[catalog]['file'] = self.FILE_dict[catalog]
                    starDict[catalog]['include'] = self.INCLUDE_dict[catalog]

    
                # Convert to Star and add to master list
                star = Star(sdict=starDict)
                theList.add(star)  #If star is already in this Set, no action
            

                                
            print'\n\tMatched %i/%i sources.' % (total_matches, len(this_cat_dict))
        print 'Found %i total sources.' % (len(theList))

        print 'Checking against lookup dictionary...'

        match_count = 0
        for star in theList:
            # Check if match in mXL
            if star.ID in self.LOOKUP_dict:
                data = self.LOOKUP_dict[star.ID]
                star.__dict__['mXL'] = dict(zip(data.dtype.names,data))
                star['mXL']['file'] = self.masterFile
                star['mXL']['include'] = False
                match_count += 1
                             
        print '\t%i/%i sources found in lookup dictionary.'%(match_count,len(theList))
        return list(theList)


def main():
    parser = argparse.ArgumentParser(description='Compile database from star catalogs.')
    parser.add_argument('catalogs',nargs='+',help='Filenames of catalogs in IPAC format.')
    parser.add_argument('--master',required=True,type=str,help='Master catalog file.')
    parser.add_argument('-o',default='catalog.json',type=str,dest='outfile',help='Output file in JSON format (Default: catalog.json)')

    args = parser.parse_args()

    LBV_data = LBV_DB(args.catalogs,args.master)
    theList = LBV_data.get_matches()

    print 'Saving JSON data to: %s' % args.outfile
    Star.save(theList, args.outfile)
    

    
if __name__ == '__main__':
    main()

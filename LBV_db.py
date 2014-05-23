#! /usr/bin/env python
import numpy as np
#import logging
#import collections
from collections import OrderedDict as dict
from star import Star
import sys
import argparse
from astropy.table import Table
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
    def __init__(self,catList,masterFile):
        for filename in catList:
            print 'Reading table from %s' % filename
            table = Table.read(filename)
            title = table.meta['TITLE']

            self.masterFile = masterFile
            self.DATA_dict[title] = table
            self.FILE_dict[title] = filename
            self.INCLUDE_dict[title] = table.meta['INCLUDE']
            self.NUM_STARS += len(table)

        print

        # Build lookup dictionary
        print 'Building lookup dictionary from %s' % masterFile
        masterTable = Table.read(masterFile)
        self.LOOKUP_dict = dict(zip(masterTable['LGGS'],masterTable))
        #for row in MasseyXLongDat:
        #    self.LOOKUP_dict[row['LGGS']] = row
        print '\t%i sources available for search' % len(self.LOOKUP_dict)
        print

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
                    that_star = LBV_DB.find_match(this_star,that_cat_dict)
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

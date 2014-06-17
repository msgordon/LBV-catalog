import numpy as np
import simplejson as json
from jsonDecoder import ConcatJSONDecoder
import collections
from collections import OrderedDict as dict
from copy import deepcopy
import HTML
#### Requires python2.7

#################
# Star class holds object references from each dataset
class Star:
    def __call__(self,*args):
        return self

    # Constructor builds Star from dictionary, generates unique ID
    def __init__(self,sdict):
        self.__dict__.update(sdict)

        # Generate ID and coordinates
        self.RAs, self.DECs, self.ID = self.gen_coord()

        # Add degree form coordinates
        self.RAd, self.DECd = self.get_degrees()

        # Add colon form coordinates
        self.RA, self.DEC = self.colon_coord()

    def __getitem__(self,val):
        return self.__dict__[val]

    # Return string in JSON format
    def __str__(self):
        return json.dumps(self.toJSON(),sort_keys=True)
       
    def __repr__(self):
        return self.__str__()

    def __eq__(self,other):
        return self.__hash__() == other.__hash__()

    # Make class hashable for Set
    def __hash__(self):
        return hash(self.ID)

    '''
    # Display ID should match 2MASS or WISE, if they exist
    #  Else return self.ID (with padded 0, if necessary)
    def get_dID(self):
        if   getattr(self,'MASSsource',None) is not None:
            return self['MASSsource']
        elif getattr(self,'WISEsource',None) is not None:
            return self['WISEsource']
        else:
            # Pad zero if necessary
            if self.ID[-2] == '.':
                return self.ID+'0'
            '''
            
    def get_dID(self):
        # Pad zero if neceesary
        #if self.ID[-2] == '.':
        #    return self.ID+'0'
        return self.ID
        

    # Display RA should match 2MASS or WISE, if they exist
    def get_dRA(self):
        starID = self.get_dID()
        return starID[1:3]+':'+starID[3:5]+':'+starID[5:10]
        
    # Display DEC should match 2MASS or WISE, if they exist
    def get_dDEC(self):
        starID = self.get_dID()
        return starID[10:13]+':'+starID[13:15]+':'+starID[15:]

    # Get table class for Javascript
    def get_class(self):
        Class = ' '.join([x.upper() for x in self.get_catalogs(get_all=True)])
        if len(self.get_catalogs()) > 1:  # duplicates are present
            Class += ' dupl'
        if self.MASStable is not None:
            Class += ' 2SS'
        if self.WISEtable is not None:
            Class += ' ISE'
        if self.get_dID()[2] == '1':
            Class += ' m33'
        else:   #elif star.get_dID()[3] == '4' or star.get_dID()[3] =='3':
            Class += ' m31'
        Class += ' data'  # include class for all rows

        return Class
        
    # Get names from star as dictionary {'w':'W(name)','m':'M(name)'} etc
    def get_names(self):
        names = dict()
        catList = self.get_catalogs(get_all=False)
        for cat in catList:
            if 'Name' in self[cat]:
                if ((self[cat]['Name'] != '') and (self[cat]['Name'] != u'')):
                    names[cat] = '%s(%s)' % (cat.upper(),self[cat]['Name'])

        return names
        
    
    # Get catalogs from star ['m','h','w'] etc
    #  If get_all == True, get long catalogs as well (['mL','mqv'])
    def get_catalogs(self,get_all=False):
        catList = []
        attrList = [x for x in dir(self) if (not callable(getattr(self,x))) and (not x.startswith('__')) and (not x.endswith('table'))]
        for attr in attrList:
            a = getattr(self,attr)
            if isinstance(a,dict):  #indicates catalog
                if get_all:
                    catList.append(attr)
                else:
                    if a['include'] == True:   #indicates not MasseyLong etc
                        catList.append(attr)

        catList.sort()
        return catList

    # Get photometry, return as dictionary {'V':'19.0','B_V':'0.32'} etc
    def get_photometry(self):
        # Photometric tags
        phot = ['V','J','H','K','B_V','U_B','U','B','R','I','V_R','R_I','__3_6_','__4_5_','__5_8_','__8_0_']

        photDict = dict.fromkeys(phot,'')  # Holds phot values
        
        catList = self.get_catalogs(get_all=True)

        # For each tag, find value in star (if it exists)
        for p in phot:
            for cat in catList:
                # Check if p is a key in this catalog, else continue
                if p not in self[cat]:
                    continue

                # Prioritize Massey
                if (cat == 'mL') or (cat == 'm') or (cat == 'mXL'):
                    photDict[p] = self[cat][p]
                    break
                else:
                    photDict[p] = self[cat][p]
                    continue

        # Special case, read J,H,K from 2MASS data
        if getattr(self,'MASStable',None) is not None:
            photDict['J'] = self['MASStable']['j_m']
            photDict['H'] = self['MASStable']['h_m']
            photDict['K'] = self['MASStable']['k_m']

        # Truncate all floats to strings
        for k,v in photDict.iteritems():
            if v == '' or v == u'':
                photDict[k] = ''
                continue
            if not isinstance(v,str):
                photDict[k] = np.around(v,decimals=2).astype('|S5')
            if (photDict[k] == '99.99') or (photDict[k] == 'nan'):
                photDict[k] = ''
        

        return photDict
    
    # Generate coordinates in form:
    #   (00:11:22.22, +33:44:55.5)
    def colon_coord(self):
        RA = []
        for x in self.RAs:
            if x.isdigit():
                RA.append(x)
            if x == '+':
                RA.append(':')
            if x == '.':
                RA.append('.')
        DEC = []
        for x in self.DECs:
            if x.isdigit():
                DEC.append(x)
            if x == '+':
                DEC.append(':')
            if x == '.':
                DEC.append('.')

        RA = ''.join([x for x in RA])
        DEC = '+' + ''.join([x for x in DEC])

        # RA leading 0
        if RA[1] == ':':
            RA = '0'+RA
        return (RA, DEC)
        
    # Generate coordinates and ID in form:
    #   ('00h11m22.2s','33d44m55.5s','J001122.22+334455.5')
    def gen_coord(self):
        RA = ''
        DE = ''
        
        # Prioritize coordinates from Massey and Humphreys
        if getattr(self,'m',None) is not None:
            RA = self.m['RA']
            DE = self.m['DEC']
        elif getattr(self,'mL',None) is not None:
            RA = self.mL['RA']
            DE = self.mL['DEC']
        elif getattr(self,'h',None) is not None:
            RA = self.h['RA']
            DE = self.h['DEC']
        else:  # choose from available catalogs
            attrList = [x for x in dir(self) if (not callable(getattr(self,x))) and (not x.startswith('__')) and (not x.endswith('table'))]
            attrList.sort()
            for attr in attrList:
                a = getattr(self,attr)
                if isinstance(a,dict): #indicates catalog
                    if (a.get('RA') is not None) and (a.get('DEC') is not None):
                        RA = a['RA']
                        DE = a['DEC']
                        break

        if (RA == '') or (DE == ''):
            # If Star was written without any catalogs?
            raise AttributeError('No coordinates found')

        RAt = str(np.fix(RA))
        RAstr =[RAt[:-6],'h+',RAt[-6:-4],'m+',RAt[-4:-2],str(np.around(RA-np.fix(RA),decimals=2))[1:],'s']
        if not RAstr[0]:
            RAstr[0] = '00'
        RAfin = ''.join([thing for thing in RAstr])

        DEt = str(np.fix(DE))
        DEstr =[DEt[:-6],'d+',DEt[-6:-4],'m+',DEt[-4:-2],str(np.around(DE-np.fix(DE),decimals=2))[1:],'s']
        if not DEstr[0]:
            DEstr[0]= '0'
        DEfin = ''.join([thing for thing in DEstr])

        RAlast = str(np.around(RA-np.fix(RA),decimals=2))[1:]
        if len(RAlast) < 3:
            RAlast = RAlast + '0'
        RAID = ''.join([x for x in [RAstr[0],RAstr[2],RAstr[4],RAlast]])
        # Leading zero
        if (float(RAID) <= 100000) and (float(RAID) >=10000):
            RAID = '0'+RAID

        DEID = ''.join([x for x in [DEstr[0],DEstr[2],DEstr[4],str(np.around(DE-np.fix(DE),decimals=1))[1:]]])

        ID = ''.join([x for x in ['J',RAID,'+',DEID]])

        # If star in Massey, ID should just be LGGS
        if getattr(self,'m',None) is not None:
            ID = self.m['LGGS']
        elif getattr(self,'mL',None) is not None:
            ID = self.mL['LGGS']
        return (RAfin,DEfin,ID)

    # Generate coordinates in degree decimal form
    def get_degrees(self):
        RAh = ''
        RAm = ''
        RAs = ''
        flag = 'h'
        for x in self.RAs:
            if flag == 'h':
                if x is not 'h':
                    RAh += x
                else:
                    flag = 'm'
                    continue
            if flag == 'm':
                if x is not 'm':
                    RAm += x
                else:
                    flag = 's'
                    continue
            if flag == 's':
                if x is not 's':
                    RAs += x
                    
        RAh = float(RAh)*15.0
        RAm = float(RAm)*0.25
        RAs = float(RAs)/240.0

        RAd = RAh+RAm+RAs

        DECd = ''
        DECm = ''
        DECs = ''
        flag = 'd'
        for x in self.DECs:
            if flag == 'd':
                if x is not 'd':
                    DECd += x
                else:
                    flag = 'm'
                    continue
            if flag == 'm':
                if x is not 'm':
                    DECm += x
                else:
                    flag = 's'
                    continue
            if flag == 's':
                if x is not 's':
                    DECs += x
                    
        DECd = float(DECd)
        DECm = float(DECm)/60.0
        DECs = float(DECs)/3600.0
        
        DECd = DECd+DECm+DECs
        return (RAd,DECd)


    # Generate ID from coords 
    @staticmethod
    def gen_ID(RA,DEC):
        RAo = []
        for x in RA:
            if x.isdigit():
                RAo.append(x)
            if x == '.':
                RAo.append('.')

        DECo = []
        for x in DEC:
            if x.isdigit():
                DECo.append(x)
            if x == '.':
                DECo.append('.')

        RAo = ''.join([x for x in RAo])
        DECo = '+'+''.join([x for x in DECo])
        return 'J'+RAo+DECo
    

    # Convert Star to JSON
    def toJSON(self):
        attrList = [x for x in dir(self) if (not callable(getattr(self,x))) and (not x.startswith('__')) and (not x.endswith('table'))]

        attrList.sort()
        attr = []
        for x in attrList:
            a = getattr(self,x)
            if a is not None:
                #  Numpy dtypes not jsonable?
                if isinstance(a,np.generic):
                    a = np.asscalar(a)

                if isinstance(a,dict):
                    for k,v in a.items():
                        if isinstance(v,np.generic):
                            a[k]= np.asscalar(v)
                attr.append(a)
            else:
                attr.append('')
            

        return dict(zip(attrList, attr))

    
    # Save starList to JSON
    @staticmethod
    def save(starList,filename):
        starList = list(starList)
        f = open(filename, 'w')
        for star in starList:
            stardict = star.toJSON()
            json.dump(stardict,f, sort_keys=True,indent=4)
            f.write('\n')
        f.close()

    # Rad starList from JSON
    @staticmethod
    def load(filename):
        f = open(filename, 'r')
        stars = json.load(f, cls=ConcatJSONDecoder,object_pairs_hook=collections.OrderedDict)
        f.close()

        starList = []
        for star in stars:
            for k,v in star.iteritems():
                if v is u'':
                    star[k] = None
            starList.append(Star(sdict=star))
        return starList


    # Get hoverdata as dictionary, keyed by catalog
    def get_hover(self,get_all=False):
        hover = dict()
        # For each catalog, zip attributes into string
        for cat in self.get_catalogs(get_all):
            pairs = deepcopy(self[cat])   # Do not want to edit dictionary
            # Truncate number attributes/delete file,include attribs
            for k,v in pairs.iteritems():
                if (not isinstance(v,str) and (not isinstance(v,unicode))):
                    pairs[k] = np.around(v,decimals=2).astype('|S5')
                if (v == '') or (v == u'') or (k == 'file') or (k == 'include'):
                    pairs[k] = ''
            pairs = zip(pairs.keys(),pairs.values())

            # Remove blank or unicode (or 'include')
            hover[cat] = ' '.join([`x`.strip('"\'')+':'+`y`.strip('"\'') for x,y in pairs if y != ''])
        return hover
                
        
    # Get name in HTML-ready form with hoverdata
    def get_namesHTML(self):
        html = []
        names = self.get_names()
        if len(names) == 0:
            return ''  # Return empty string if no names
        hover = self.get_hover()

        # For each name, build html hover text
        for k,v in names.iteritems():
            html.append('<span style="color:#DF0101" title="%s">%s</span>' % (hover[k],v))

        return ', '.join(html)

    # Catalogs in HTML-ready form with hoverdata
    def get_catalogsHTML(self):
        html = []
        cats = self.get_catalogs(get_all=True)
        hover = self.get_hover(get_all=True)

        # Remove MasseyLong from list if M already there
        if ('m' in cats) or ('mL' in cats):
            if ('m' in cats) and ('mL' in cats):
                cats.pop(cats.index('mL'))
            if 'mXL' in cats:
                cats.pop(cats.index('mXL'))
        # For each catalog, build html hover text
        for k in cats:
            html.append('<span style="color:#DF0101" title="%s">%s</span>' % (hover[k],k.upper()))

        return ', '.join(html)        

    # Render cells for 2MASS/WISE matches.  Return '' if no match
    def get_table_link(self,filename):
        return HTML.TableCell('<a href="%s" style="color:#0000FF">Y</a>' % (filename),align="center")

    # Get priority from file
    def get_priorityHTML(self,priorityDict):
        prior = priorityDict.get(self.get_dID(),'')
        if prior == '':
            attrib = 99
        else:
            attrib = prior
        return HTML.TableCell(prior,align="center",attribs={'sorttable_customkey':attrib})


    def get_gal_name(self):
        if self.get_dID()[2] == '1':
            return 'M33'
        else:
            return 'M31'


    ### THIS ONE IS THE SHIT
    @staticmethod
    def sex2deg(RA,DE):
        RAt = str(np.fix(RA))
        RAstr =[RAt[:-6],'h+',RAt[-6:-4],'m+',RAt[-4:-2],str(np.around(RA-np.fix(RA),decimals=2))[1:],'s']
        if not RAstr[0]:
            RAstr[0] = '00'
        RAfin = ''.join([thing for thing in RAstr])

        DEt = str(np.fix(DE))
        DEstr =[DEt[:-6],'d+',DEt[-6:-4],'m+',DEt[-4:-2],str(np.around(DE-np.fix(DE),decimals=2))[1:],'s']
        if not DEstr[0]:
            DEstr[0]= '0'
        DEfin = ''.join([thing for thing in DEstr])


        RAh = ''
        RAm = ''
        RAs = ''
        flag = 'h'
        for x in RAfin:
            if flag == 'h':
                if x is not 'h':
                    RAh += x
                else:
                    flag = 'm'
                    continue
            if flag == 'm':
                if x is not 'm':
                    RAm += x
                else:
                    flag = 's'
                    continue
            if flag == 's':
                if x is not 's':
                    RAs += x
                    
        RAh = float(RAh)*15.0
        RAm = float(RAm)*0.25
        RAs = float(RAs)/240.0

        RAd = RAh+RAm+RAs

        DECd = ''
        DECm = ''
        DECs = ''
        flag = 'd'
        for x in DEfin:
            if flag == 'd':
                if x is not 'd':
                    DECd += x
                else:
                    flag = 'm'
                    continue
            if flag == 'm':
                if x is not 'm':
                    DECm += x
                else:
                    flag = 's'
                    continue
            if flag == 's':
                if x is not 's':
                    DECs += x
                    
        DECd = float(DECd)
        DECm = float(DECm)/60.0
        DECs = float(DECs)/3600.0
        
        DECd = DECd+DECm+DECs
        return (RAd,DECd)

#! /usr/bin/env python

#### Requires virtualenv, "activateVE"
#### Requires HTML.py http://www.decalage.info/python/html
#### Requires jsonDecoder.py http://stackoverflow.com/questions/8730119/
#### Requires python2.7

import sys
import os
import numpy as np
from urllib import urlretrieve
import xml.etree.ElementTree as ET
import atpy
import logging
import argparse
from star import Star
from LBV_db import LBV_DB
import HTML
import collections
from collections import OrderedDict as dict


############
# Get priority dictionary from file
############
def get_priority_dict(filename):
    f = open(filename,'r')
    priorityDict = dict()
    for line in f:
        line = line.split()
        priorityDict[line[0]] = line[1]

    f.close()
    return priorityDict
    
##############
# Query URL.  If 'retrieve=False', build object but don't actually query
##############
def query(URL,star,radius,catalog,output,retrieve=True,verbose=True,overwrite=False,forcefind=False):
        
    URL[1] = star.RAs
    URL[3] = star.DECs
    URL[5] = `radius`
    URL[7] = catalog
    query = ''.join([thing for thing in URL])
    
    if retrieve:
        if overwrite == False:
            if os.access(output,os.F_OK):
                print('Skipping  '+ star.ID)
                return query
        if forcefind:
            if read_table(output) == None:
                URL[5] = `radius + 1`
                if catalog == 'fp_psc':
                    URL[7] = 'pt_src_6x2'
                else:
                    URL[7] = 'wise_allwise_p3as_psd'
                query = ''.join([thing for thing in URL])
                print 'Retrieving:  ' + query
                urlretrieve(query,output)
                print 'Writing to:  ' + output
                return query
                
            else:
                print('Skipping  '+ star.ID)
                return query
        else:
            urlretrieve(query,output)
            print 'Writing to:  ' + output
 
    return query

###################
#  Read in table data
###################
def read_table(filename,verbose=False):
    try:
        print 'Reading table:  ' + filename
        t = atpy.Table(filename,verbose=False)
    except:
        print 'No sources found in ' + filename
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

###########
# Query images.  If 'retrieve=False', build object but don't actually query
############
def query_images(star,url,xmloutput,verbose=False,retrieve=True,overwrite=False):
    log = logging.getLogger(__file__+':query_images ')
    if verbose:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.INFO)

    if retrieve:
        if overwrite == False:
            if os.access(xmloutput,os.F_OK):
                log.debug('Skipping '+ star.ID)
                return url
        log.debug('Retrieving:  ' + url)
        urlretrieve(url,xmloutput)
        log.debug('Writing to:  ' + xmloutput)

        # Retrieve images
        try:
            tree = ET.parse(xmloutput)
            root = tree.getroot()
        except:
            #weird bug where xml file isn't always written
            #so, do a wget
            os.system('wget "'+url+'" --output-document='+xmloutput)
            tree = ET.parse(xmloutput)
            root = tree.getroot()

        # Parse cutouts and thumbnails
        star.numImages = int(root.find('.//result/totalimages').text.strip())
        star.shrunk = []
        star.big = []
        for tag in root.findall('.//result/image'):
            star.shrunk.append(tag.find('shrunkjpgurl').text.strip())
            star.big.append(tag.find('jpgurl').text.strip())

        # If folder doesn't exist, create it
        assets = './dMASSxml/'+star.ID
        if overwrite == False:
            if os.access(assets,os.F_OK):
                log.debug('Skipping '+assets)
                return
            
        if not os.access(assets,os.F_OK):
            os.mkdir(assets)
        for site in star.shrunk:
            log.debug('Retrieving:  ' + site)
            output = assets+'/'+os.path.basename(site)
            urlretrieve(site,output)
            log.debug('Writing to:  ' + output)
        for site in star.big:
            log.debug('Retrieving:  ' + site)
            output = assets+'/'+os.path.basename(site)
            urlretrieve(site,output)
            log.debug('Writing to:  ' + output)


#########
# Build subpage
#########
def build_subpage(star,xmlfile,build=False):
    tree = ET.parse(xmlfile)
    root = tree.getroot()

    rows = []
    rows.append(['<h2>'+root.find('.//result/equCoord').text.strip()+'</h2>'])
    rows[0].append(root.find('.//result/eclCoord').text.strip() +'&nbsp;&nbsp;' +root.find('.//result/galCoord').text.strip())
    rows[0].append('Input Location: '+root.find('.//input/locstr').text.strip()+'&nbsp;&nbsp;&nbsp;&nbsp;'+root.find('.//input/subsetsize').text.strip()+'&nbsp;&nbsp;&nbsp;North up, East to the left')
    rows[0].append('Data tag: '+root.find('.//result/datatag').text.strip())
    rows[0] = '<br>'.join([x.strip('"\'') for x in rows[0]])

    base = ''
    rows.append([])
    #  DSS1
    for dss in root.findall('.//result/image'):
        if dss.find('.//surveyname').text.strip() != 'DSS':
            continue
        if dss.find('.//band').text.strip()[3] != '1':
            continue
        rows[1].append(HTML.TableCell('<a class="box" href="%s"> <img src="%s" border="1" alt="%s" />\n</a><br>\nDSS %s band<br>\nObs Date: %s' % (base+os.path.basename(dss.find('jpgurl').text.strip()),base+os.path.basename(dss.find('shrunkjpgurl').text.strip()),base+os.path.basename(dss.find('shrunkjpgurl').text.strip()),dss.find('band').text.strip(),dss.find('obsdate').text.strip())))

    rows.append([])
    # DSS2
    for dss in root.findall('.//result/image'):
        if dss.find('.//surveyname').text.strip() != 'DSS':
            continue
        if dss.find('.//band').text.strip()[3] != '2':
            continue
        rows[2].append(HTML.TableCell('<a class="box" href="%s"> <img src="%s" border="1" alt="%s" />\n</a><br>\nDSS %s band<br>\nObs Date: %s' % (base+os.path.basename(dss.find('jpgurl').text.strip()),base+os.path.basename(dss.find('shrunkjpgurl').text.strip()),base+os.path.basename(dss.find('shrunkjpgurl').text.strip()),dss.find('band').text.strip(),dss.find('obsdate').text.strip())))

    rows.append([])
    #2MASS
    for dss in root.findall('.//result/image'):
        if dss.find('.//surveyname').text.strip() != '2MASS':
            continue
        rows[3].append(HTML.TableCell('<a class="box" href="%s"> <img src="%s" border="1" alt="%s" />\n</a><br>\n2MASS %s band<br>\nObs Date: %s' % (base+os.path.basename(dss.find('jpgurl').text),base+os.path.basename(dss.find('shrunkjpgurl').text),base+os.path.basename(dss.find('shrunkjpgurl').text),dss.find('band').text.strip(),dss.find('obsdate').text.strip())))
        
    htmlpage = HTML.Table(rows[1:])
    buildpage = './dMASSxml/'+star.ID+'/'+star.ID+'.html'

    if build:
        f = open(buildpage,'w')
        f.write(rows[0])
        f.write(str(htmlpage))
        f.close()
    return


############
# Prepare row for each star.  Return as list.
############
def prepare_row(star,columns):
    row = []
    # For each column, pull appropriate data
    for col in columns:
        if '-' in col:   # Numpy recarrays do not allow '-' in keys
            col = col.replace('-','_')
            
        current = None   # This will add to row
        # If top-level attribute (ID,RA,DEC), just use that
        if getattr(star,col,None) is not None:
            current = star[col]
                

#########
#  Build main html page to htmlfile
############
def publish_html(starList,htmlfile,priorityDict,sortby='ID',write=True):
    log = logging.getLogger(__file__+':publish_html ')
    log.info('Building HTML Catalog')
    
    # Sort list
    starList.sort(key=lambda x: x[sortby])

    # htmlcode holds the table
    #htmlcode = HTML.Table(header_row=['ID','RA','DEC','Name','V','B-V','U-B','J','H','K','Catalog','2MASS','WISE'],attribs={'class':'sortable','id':'stardata'})

    #m31Rhtmlcode = HTML.Table(header_row=['ID','RA','DEC','Priority','V','B-V','U-B','J','H','K','Catalog','2MASS','WISE'],attribs={'class':'sortable m31R','id':'stardata31R'})



    # Original Header
#    m33Rhtmlcode = HTML.Table(header_row=['ID','RA','DEC','Priority','V','B-V','U-B','J','H','K','Catalog','2MASS','WISE'],attribs={'class':'sortable m33R','id':'stardata33R'})
#    m31Yhtmlcode = HTML.Table(header_row=['ID','RA','DEC','Priority','V','B-V','U-B','J','H','K','Catalog','2MASS','WISE'],attribs={'class':'sortable m31Y','id':'stardata31Y'})
#    m33Yhtmlcode = HTML.Table(header_row=['ID','RA','DEC','Priority','V','B-V','U-B','J','H','K','Catalog','2MASS','WISE'],attribs={'class':'sortable m33Y','id':'stardata33Y'})

    m33Rhtmlcode = HTML.Table(header_row=['ID','RA','DEC','V','B-V','U-B','V-R','R-I','J','H','K','Catalog','2MASS','WISE'],attribs={'class':'sortable m33R','id':'stardata33R'})
    m31Yhtmlcode = HTML.Table(header_row=['ID','RA','DEC','V','B-V','U-B','V-R','R-I','J','H','K','Catalog','2MASS','WISE'],attribs={'class':'sortable m31Y','id':'stardata31Y'})
    m33Yhtmlcode = HTML.Table(header_row=['ID','RA','DEC','V','B-V','U-B','V-R','R-I','J','H','K','Catalog','2MASS','WISE'],attribs={'class':'sortable m33Y','id':'stardata33Y'})


    # For use in 'button'
    hidden = 0   #number of hidden rows at the beginning
    
    # Append table rows for each star
    for star in starList:
        ID = star.get_dID()
        RA = star.get_dRA()
        DEC = star.get_dDEC()
        Name = star.get_namesHTML()
        priority = star.get_priorityHTML(priorityDict)


        # Get photometry
        phot = star.get_photometry()
        
        RI = phot['R_I']
        VR = phot['V_R']
        
        V = phot['V']
        BV = phot['B_V']
        UB = phot['U_B']
        J = phot['J']
        H = phot['H']
        K = phot['K']

        BV = HTML.TableCell(BV,attribs={'class':'BVdata'})
        UB = HTML.TableCell(UB,attribs={'class':'UBdata'})
        V  = HTML.TableCell(V,attribs={'class':'Vdata'})
        
        Catalog = star.get_catalogsHTML()

        MASS = ''
        WISE = ''
        if star.MASStable is not None:
            MASS = star.get_table_link(star.MASSfile)
            # Link 2MASS to ID
            ID = HTML.link(ID,'./dMASSxml/%s/%s.html'%(star.ID,star.ID))
        if star.WISEtable is not None:
            WISE = star.get_table_link(star.WISEfile)

        # Set table class for javascript
        Class = star.get_class()

        # Set id class
        ID = HTML.TableCell(ID,attribs={'class':'starID'})

        # Set H spectra to display none
        style = ''
        rowID = ''
        #if 'H' in Class:
        #    style = 'display: none;'
        #    rowID = 'H'
        #    hidden += 1

        #row = [ID,RA,DEC,Name,V,BV,UB,J,H,K,Catalog,MASS,WISE]
        row = [ID,RA,DEC,V,BV,UB,VR,RI,J,H,K,Catalog,MASS,WISE]
        #row = [ID,RA,DEC,priority,V,BV,UB,J,H,K,Catalog,MASS,WISE]
        if 'D31' in Class:
            m31Yhtmlcode.rows.append(HTML.TableRow(row,attribs={'class':Class,'style':style,'id':rowID}))
        elif 'D33_R' in Class:
            m33Rhtmlcode.rows.append(HTML.TableRow(row,attribs={'class':Class,'style':style,'id':rowID}))
#        elif 'm31Y' in Class:
#            m31Yhtmlcode.rows.append(HTML.TableRow(row,attribs={'class':Class,'style':style,'id':rowID}))
        elif 'D33_Y' in Class:
            m33Yhtmlcode.rows.append(HTML.TableRow(row,attribs={'class':Class,'style':style,'id':rowID}))


    # Add header information
    headerstuff = '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"\n   "http://www.w3.org/TR/html4/strict.dtd">\n<HTML>\n   <HEAD>\n      <TITLE>Drout Supergiants</TITLE>\n  <link href="./assets/tabs.css" rel="stylesheet" type="text/css">\n <link href="./assets/LBV.css" rel="stylesheet" type="text/css">\n     <style>table.sortable thead {\n    background-color:#eee;\n    font-weight: bold;\n    cursor: default;\n}</style>\n'

    # Add scripting
    scripting = '<script language="javascript" type="text/javascript" src="./assets/jquery.min.js"></script>\n<script language="javascript" type="text/javascript" src="./assets/jquery-ui.min.js"></script>\n<script language="javascript" type="text/javascript" src="./assets/jquery-ui.js"></script>\n<script src="sorttable.js"></script>\n<script src="html2tsv.js"></script>\n<script src="toggleMult.js"></script>\n<script src="contentSelect.js"></script>\n<script type="text/javascript" src="spawn.js"></script>\n'
    headerstuff += scripting +'</HEAD>\n   <BODY>\n'

    # Notes/body
    notes = '<h1>Drout Supergiants</h1>\n<h2><a href="./LBV.html">LBV Candidates</a></h2>\n<p>Click on the ID column for IPAC archive images centered on the coordinates.<br>Click on the 2MASS/WISE columns to go to the data set for the closest matching point source(s) in each survey.<br>2MASS matches are the closest sources within 3 arcseconds of the target coordinates.<br>WISE matches are the closest sources within 10 arcseconds of the target coordinates.<br>Hover over the Name (or Catalog) of each entry to see star data from the corresponding paper.<br>J, H, K magnitudes are from the 2MASS catalogs.<br>The table can be sorted by clicking on the header of each column.<br>The displayed table can be exported to a tab-delimited CSV file.</p>\n'

    # Buttons
    #button = """<form>\n<label><input type="checkbox" checked="checked" onClick="toggle('dupl');"/>Duplicates</label>\n"""
    button = """<form>\n<label><input type="checkbox" checked="checked" onClick="toggle('H');"/>H Spectra</label>\n"""
    button += """&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\n<input type="button" value="Export to CSV" id="linktotsv" onclick="html2tsv('stardata');"/>\n"""
    button +="""<input type="button" value="Plot CMD" onclick="spawn('CMD','B-V','V','BVdata','Vdata');"/>\n"""
    button +="""<input type="button" value="Plot Two Color" onclick="spawn('TwoColor','B-V','U-B','BVdata','UBdata');"/>\n"""
    button += """&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"""
    #button += '<span id="rowCount">%i/%i sources displayed.</span>' % (len(starList)-hidden,len(starList))
    button +='</form><br>\n'

    #  Tabs
    tabs = """<div id="tabs" class="tabs">\n <ul>\n  <li><a href="#tabs-m31">M31</a></li>\n  <li><a href="#tabs-m33">M33</a></li>\n </ul>\n"""
    tabs += """<div id="tabs-m31">\n"""
    tabs += """  <div id="subtabs31">\n  <ul>\n<li><a href="#subtabs31Y">Yellow</a></li>\n</ul>\n"""
#    tabs += """  <div id="subtabs31">\n  <ul>\n<li><a href="#subtabs31R">Red</a></li>\n<li><a href="#subtabs31Y">Yellow</a></li>\n</ul>\n"""
#    tabs += """<div id="subtabs31R">\n"""
#    tabs += '<span id="rowCount31R">%i/%i sources displayed.</span><br>\n' % (len(m31Rhtmlcode.rows),len(m31Rhtmlcode.rows))
#    tabs += str(m31Rhtmlcode)
#    tabs += """</div>\n"""
    tabs +="""<div id="subtabs31Y">\n"""
    tabs += '<span id="rowCount31Y">%i/%i sources displayed.</span><br>\n' % (len(m31Yhtmlcode.rows),len(m31Yhtmlcode.rows))
    tabs += str(m31Yhtmlcode)
    tabs += """</div>\n</div>\n</div>\n"""
    
    tabs += """<div id="tabs-m33">\n"""
    tabs += """  <div id="subtabs33">\n  <ul>\n<li><a href="#subtabs33R">Red</a></li>\n<li><a href="#subtabs33Y">Yellow</a></li>\n</ul>\n"""
    tabs += """<div id="subtabs33R">\n"""
    tabs += '<span id="rowCount33R">%i/%i sources displayed.</span><br>\n' % (len(m33Rhtmlcode.rows),len(m33Rhtmlcode.rows))
    tabs += str(m33Rhtmlcode)
    tabs += """</div>\n"""
    tabs +="""<div id="subtabs33Y">\n"""
    tabs += '<span id="rowCount33Y">%i/%i sources displayed.</span><br>\n' % (len(m33Yhtmlcode.rows),len(m33Yhtmlcode.rows))
    tabs += str(m33Yhtmlcode)
    tabs += """</div>\n</div>\n</div></div>\n"""
    
    tabs += """<script>\n$("#tabs,#subtabs31,#subtabs33").tabs();\n</script>"""
    
    website = ''.join([headerstuff,notes,button,tabs,'\n  </BODY>\n</HTML>'])
    # Write to file, or return string
    if write:
        log.info('\tWriting catalog to %s' % (htmlfile))
        f = open(htmlfile,'w')
        f.write(website)
        f.close()
        log.info('\t%i sources written to HTML' % (len(starList)))
        return
    
    return website




##### ICRS correction
def icrs_correction(RAs,DECs):
    RA = str(np.float(RAs.replace(':',''))-0.03)
    idx = RA.find('.')
    if idx == 4:  #4347.98, add 00
        RA = ''.join(['00:',RA[idx-4:idx-2],':',RA[idx-2:]])
    if idx == 5:
        RA = ''.join(['01:',RA[idx-4:idx-2],':',RA[idx-2:]])
    

    DEC = str(np.float(DECs[1:].replace(':',''))-0.14)
    DEC = ''.join([DEC[0:2],':',DEC[2:4],':',DEC[4:]])

    return (RA,DEC)

    
#########
## output hectospec catalog
###########
def hecto_catalog(starList, filename, priorityDict):
    #header = 'ra\tdec\tobject\tmag\trank\ttype\n--\t---\t------\t----\t---\t----\n'
    print 'Writing catalog to %s...' % filename
    f = open(filename,'w')
    #f.write(header)
    for star in starList:
        prior = priorityDict.get(star.get_dID(),'')
        if prior != '':
            # if in Massey, perform correction
            RA = star.get_dRA()
            DEC = star.get_dDEC()
            if [i for i in ['m','mXL','mL'] if i in star.get_catalogs(get_all=True)]:
                RA,DEC = icrs_correction(RA,DEC)
            else:
                DEC = DEC[1:] # gets rid of '+' sign

            #print star.get_dID(),RA,DEC
            nameDict = star.get_names()
            if 'h' in nameDict:
                name = nameDict['h']
            elif 'w' in nameDict:
                name = nameDict['w']
            elif 'v' in nameDict:
                name = nameDict['v']
            else:
                name = star.get_dID().split('+')[0][1:]   # get name as RA only

            name = 'D(%s)' % name

            line = '\t'.join([RA,DEC,name,star.get_photometry()['V'],prior,'target'])
            f.write(line)
            f.write('\n')
    f.close()
    return

    

def main():

    parser = argparse.ArgumentParser(description='Find matches of input JSON catalog in WISE and 2MASS')

    parser.add_argument('catalog',type=str,help='JSON catalog of sources.')

    args = parser.parse_args()
    
    # Master list to hold Stars
    print('Loading JSON data from: %s' % (args.catalog))
    theList = Star.load(args.catalog)
    print('\tLoaded %i sources.' % (len(theList)))


    ###############
    # 2MASS query
    ###############
    URL = ['http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?outfmt=1&objstr=','','+','','&spatial=Cone&radius=','','&radunits=arcsec&catalog=','']
    htmlbase='http://irsa.ipac.caltech.edu/cgi-bin/FinderChart/nph-finder?locstr='
    htmlend = '&markervis_orig=true&markervis_shrunk=true&mode=prog'
    
    for star in theList:
        star.MASSfile = './dMASStables/%s.tbl' % (star.ID)

        # Query URL to find 2MASS sources
        star.MASSquery = query(URL,star,3,'fp_psc',star.MASSfile,retrieve=True,verbose=True,overwrite=True,forcefind=True)
        #star.MASStable = read_table(star.MASSfile,verbose=True)

        # If found, find closest source
        #if star.MASStable is not None:
            #star.MASStable = find_closest_source(star,star.MASStable)
            #star.MASSsource = Star.gen_ID(star.MASStable['clon'],star.MASStable['clat'])

            # Parse xml files and pull images
            #star.MASShtml = htmlbase+star.MASStable['clon']+'+'+star.MASStable['clat']+htmlend
            #star.MASSxml = './dMASSxml/'+star.ID+'.xml'
            #query_images(star,star.MASShtml,star.MASSxml,retrieve=True,verbose=False,overwrite=False)
            #build_subpage(star,star.MASSxml,build=True)


    #################
    # WISE query
    ##################
    for star in theList:
        star.WISEfile = './dWISEtables/' + star.ID+'.tbl'
        # Query URL to find WISE sources
        star.WISEquery = query(URL,star,3,'wise_allsky_4band_p3as_psd',star.WISEfile,retrieve=True,verbose=True,overwrite=True,forcefind=True)
        #star.WISEtable = read_table(star.WISEfile)

        # If found, get closest source
        #if star.WISEtable is not None:
        #    star.WISEtable = find_closest_source(star,star.WISEtable)
        #    star.WISEsource = Star.gen_ID(star.WISEtable['clon'],star.WISEtable['clat'])


    exit()

    #########
    # Make HTML catalog
    ##########
    htmlfile = 'test.html'
    priorityDict = get_priority_dict('priority.tsv')
    publish_html(theList,htmlfile,priorityDict)
    ###########
    # Output fiber assignment table
    ##########
    m31 = []
    m33 = []
    for star in theList:
        if 'm31' in star.get_class():
            m31.append(star)
        elif 'm33' in star.get_class():
            m33.append(star)
    hecto_catalog(m31,'hecto_LBV_m31.tsv',priorityDict)
    hecto_catalog(m33,'hecto_LBV_m33.tsv',priorityDict)        

if __name__ == "__main__":
    sys.exit(main())

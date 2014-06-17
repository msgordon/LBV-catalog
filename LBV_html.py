#! /usr/bin/env python
import argparse
import HTML
from astropy.table import Table



def publish_html(tableFile,outfile):
    table = Table.read(tableFile)

    writecols = ['ID','Gal','U','B','V','R','I','J','H','K','3.6','4.5','5.8','8.0']

    htmlcode = HTML.Table(header_row=writecols,attribs={'class':'sortable','id':'datatable'})

    for star in table:
        if star['Gal'] == 'M33':
            continue

        cells = [HTML.TableCell(star[cell],attribs={'class':'cell_%s'%cell}) for cell in writecols]

        htmlcode.rows.append(HTML.TableRow(cells))

    # Add header information
    headerstuff = '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"\n   "http://www.w3.org/TR/html4/strict.dtd">\n<HTML>\n   <HEAD>\n      <TITLE>Drout Supergiants, M31</TITLE>\n   <link href="./assets/css/LBV.css" rel="stylesheet" type="text/css">\n     <style>table.sortable thead {\n    background-color:#eee;\n    font-weight: bold;\n    cursor: default;\n}</style>\n'

    # Add scripting
    scripting = '<script src="./assets/js/sorttable.js"></script>\n'
    headerstuff += scripting +'</HEAD>\n   <BODY>\n'

    # Notes/body
    notes = '<h1>Drout Supergiants, M31</h1>\n\n'

    # Build table
    website = ''.join([headerstuff,notes,str(htmlcode),'\n  </BODY>\n</HTML>'])

    f = open(outfile,'w')
    f.write(website)
    f.close()
    return

def main():
    parser = argparse.ArgumentParser(description='Generate HTML document from input star list')

    parser.add_argument('table',type=str,help='FITS table of data')
    parser.add_argument('-o',dest='outfile',type=str,default='photometry.html',help='Output HTML file (default="photometry.html")')

    args = parser.parse_args()

    publish_html(args.table,args.outfile)



if __name__ == '__main__':
    main()

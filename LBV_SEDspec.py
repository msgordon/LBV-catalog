#! /usr/bin/env python
from astropy.table import Table,Column
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import numpy as np
from scipy.optimize import curve_fit
import pyfits

h = 6.626e-27  # erg*s
c = 3.0e10     # cm/s
k = 1.381e-16  # erg/K


def get_spectrum(filename):
    spectrum,header = pyfits.getdata(filename, 0, header=True)

    if len(spectrum) != int(header['naxis1']):
        spectrum = spectrum[0]
          
    #lambda1 = header['crval1']
    #lambda2 = np.round(len(spectrum)*np.round(header['cdelt1'],2)+header['crval1'])
    #wave = np.arange(lambda1, lambda2, np.round(header['cdelt1'],2))
    waves = np.arange(0,len(spectrum))*header['CDELT1'] + header['CRVAL1']

    return (waves,spectrum)


def blackbody(waves, T, scale=1.0):
    # waves in microns

    # convert to cm
    wavesCM = waves*1e-4

    # flux in ergs/s/cm^2/cm/sr
    flux = 2*h*c**2/wavesCM**5 * 1.0/(np.exp(h*c/(k*wavesCM*T)) - 1.0)

    # rescale for fitting
    flux = scale*flux

    # return as intensity, ergs/s/cm^2
    return flux * wavesCM





def main():
    parser = argparse.ArgumentParser(description='Plot and fit of SEDs of input photometry')

    parser.add_argument('table',type=str,help='FITS table of data')
    parser.add_argument('idlist',type=str,help='.list file of ID/spec matches')
    parser.add_argument('-pdf',type=str,required=True,help='Output pdf file')

    args = parser.parse_args()

    phot_table = Table.read(args.table,format='ascii.tab')
    spec_dict = {}

    pp = PdfPages(args.pdf)

    data = np.genfromtxt(args.idlist,delimiter='\t',comments='#',names = ['ID','Bfilename','Rfilename'],dtype=['S100','S100','S100'],autostrip=True)

    for row in data:
        spec_dict[row['ID']] = [row['Bfilename'],row['Rfilename']]

    plotCols = [x for x in phot_table.columns if (x[0:3] == 'lam')]
    eCols = [x for x in phot_table.columns if (x[0:5] == 'e_lam')]
    wavesZ = np.array([float(x[6:10]) for x in plotCols])
    wavesDISP = np.linspace(wavesZ[0],wavesZ[-1],num=1000)

    for star in phot_table:
        if star['ID'] in spec_dict:
            waves = []
            phot = []
            err = []
        else:
            continue

        for x,y,e in zip(wavesZ,[star[z] for z in plotCols],[star[z] for z in eCols]):
            if y != 99.99:
                waves.append(x)
                phot.append(y)
                err.append(e)
                
        waves = np.array(waves)
        phot = np.array(phot)
        err = np.array(err)
        try:
            popt,pcov = curve_fit(blackbody,waves[0:6],phot[0:6],p0=[14000,phot[0]])
        except:
            continue

        T = popt[0]
        print 'Fitting %s, T = %.1f' % (star['ID'],T)
        fit = blackbody(wavesDISP,*popt)

        fig = plt.figure(figsize=(20,10), dpi=72)
        plt.subplot(3,1,1)
        plt.plot(np.log10(waves),np.log10(phot),'ro')
        #plt.errorbar(np.log10(waves),np.log10(phot),yerr=err/(2.303*np.log10(phot)),fmt='ro')#err/(2.303*phot)
        plt.plot(np.log10(wavesDISP),np.log10(fit),'b-',linewidth=2)
        plt.ylim([-15,-10.5])
        plt.xlabel(r'$log[\lambda]\,[\mu m]$')
        plt.ylabel(r'$log[\lambda\,F_{\lambda}]\,[erg\,sec^{-1}\,cm^{-2}]$')
        plt.title(star['ID'])
        plt.legend(['Photometry','BB(T = %.1f)'%T])


        spec_x1,spec_y1 = get_spectrum(spec_dict[star['ID']][0])
        spec_x2,spec_y2 = get_spectrum(spec_dict[star['ID']][1])
        if min(spec_x1) < min(spec_x2):
            blue_x,blue_y = (spec_x1, spec_y1)
            red_x,red_y = (spec_x2, spec_y2)
        else:
            blue_x,blue_y = (spec_x2, spec_y2)
            red_x,red_y = (spec_x1, spec_y1)
            
        
        plt.subplot(3,1,2)
        plt.plot(blue_x,blue_y,'b')
        plt.xlim([min(blue_x),max(blue_x)])
        plt.ylim([min(blue_y),max(blue_y)])
        plt.xlabel(r'$\lambda\,[\AA]$')

        plt.subplot(3,1,3)
        plt.plot(red_x,red_y,'r')
        plt.xlim([min(red_x),max(red_x)])
        plt.ylim([min(red_y),max(red_y)])
        plt.xlabel(r'$\lambda\,[\AA]$')
        

        #plt.show()
        pp.savefig(fig)
    pp.close()
        
        


if __name__ == '__main__':
    main()

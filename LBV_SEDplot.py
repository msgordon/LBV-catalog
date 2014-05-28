#! /usr/bin/env python
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import numpy as np
from scipy.optimize import curve_fit

h = 4.135668e-15  # eV*s
c = 3.0e10        # cm/s
k = 8.617e-5      # eV/K

bb = lambda wave,A,B: A/np.power(wave,5) * 1.0/(np.exp(B/wave)-1.0)
bb2= lambda wave,A,B,C,D: bb(wave,A,B) + bb(wave,C,D)

def plot_SED(tableFile):
    table = Table.read(tableFile)

    pp = PdfPages('BB2.pdf')
    
    waves = np.array([float(x[6:10]) for x in table.columns[1:]])
    for star in table:
        phot = [x if x != 99.99 else np.NaN for x in star.data][1:]

        try:
            #popt, pcov = curve_fit(bb,waves,phot,p0=(phot[0],1e-2))
            popt, pcov = curve_fit(bb2,waves,phot,p0=(phot[0],1e-2,phot[4],1e-3))
        except:
            continue

        print 'Fitting %s' % star['ID']
        T = (h*c/k)/popt[1]*10000.
        fit = bb2(waves,*popt)
        #fit = bb(waves,*popt)
        fig = plt.figure()
        plt.plot(waves,np.log10(phot),'ro')
        plt.plot(waves,np.log10(fit),'b-')
        #plt.plot(waves,bb(waves,1e-7,7000)/np.max(bb(waves,1e-7,7000)),linewidth=2)
        plt.xlabel(r'$\lambda\,[\mu m]$')
        plt.ylabel(r'$log[\lambda\,F_{\lambda}]\,[erg\,sec^{-1}\,cm^{-2}]$')
        plt.title(star['ID'])
        plt.legend(['Photometry','BB(T = %f)'%T])

        pp.savefig(fig)
        #plt.show()
        #exit()

    pp.close()
    #plt.show()
    exit()
        


def main():
    parser = argparse.ArgumentParser(description='Plot and fit of SEDs of input photometry')

    parser.add_argument('table',type=str,help='FITS table of data')

    args = parser.parse_args()

    plot_SED(args.table)



if __name__ == '__main__':
    main()

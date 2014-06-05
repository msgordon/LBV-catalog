#! /usr/bin/env python
from astropy.table import Table,Column
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import numpy as np
from scipy.optimize import curve_fit

h = 6.626e-27  # erg*s
c = 3.0e10     # cm/s
k = 1.381e-16  # erg/K

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


bbtemp = lambda wave,T: 2*h*c**2/np.power(wave,5) * 1.0/(np.exp(h*c/(k*wave*T)) - 1.0)

bb = lambda wave,A,B: A/np.power(wave,5) * 1.0/(np.exp(B/wave)-1.0)
bblog = lambda wave,A,B: A - 5.0*np.log10(wave) - np.log10(np.exp(B/wave)-1.0)

bb2= lambda wave,A,B,C,D: bb(wave,A,B) + bb(wave,C,D)
bblog2 = lambda wave,A,B,C,D: np.log10(10.0**bblog(wave,A,B) + 10.0**bblog(wave,C,D))

def plot_SED(tableFile):
    table = Table.read(tableFile)

    waves = np.array([float(x[6:10]) for x in table.columns[1:]])
    wavesDISP = np.linspace(waves[0],waves[-1],num=1000)

    pp = PdfPages('BB.pdf')
    
    for star in table:
        phot = [x if x != 99.99 else np.NaN for x in star.data][1:]
        if np.NaN in phot:
            continue

        try:
            popt,pcov = curve_fit(blackbody,waves[0:7],phot[0:7],p0=[14000,phot[0]])
        except:
            continue

        T = popt[0]
        print 'Fitting %s, T = %.1f' % (star['ID'],T)
        fit = blackbody(wavesDISP,*popt)

        '''
        f = open('thing.txt','w')
        for x,y in zip(np.log10(waves),np.log10(phot)):
            f.write('%g\t%g\n'%(x,y))
        f.close()

        exit()
        '''
        
        fig = plt.figure()
        plt.plot(np.log10(waves),np.log10(phot),'ro')
        plt.plot(np.log10(wavesDISP),np.log10(fit),'b-',linewidth=2)
        plt.ylim([-15,-10.5])
        plt.xlabel(r'$log[\lambda]\,[\mu m]$')
        plt.ylabel(r'$log[\lambda\,F_{\lambda}]\,[erg\,sec^{-1}\,cm^{-2}]$')
        plt.title(star['ID'])
        plt.legend(['Photometry','BB(T = %.1f)'%T])

        pp.savefig(fig)
    pp.close()



    exit()
    

    wavesCM = waves*1e-4

    wavesDISP = np.linspace(waves[0],waves[-1],num=100)
    #pp = PdfPages('BB2.pdf')
    
    for star in table:
        phot = [x if x != 99.99 else np.NaN for x in star.data][1:]

        print waves
        print phot
        exit()

        popt,pcov = curve_fit(bb2,waves,phot,p0=[phot[0],0.02,phot[0],0.2])
        
        #popt,pcov = curve_fit(bbtemp,wavesCM[0:2],phot[0:2],p0=14000)

        print popt
        fit = bb2(waves,*popt)
        
        plt.plot(waves,np.log10(phot),'ro')
        plt.plot(waves,np.log10(fit),'b-')
        plt.plot(waves,np.log10(bb(waves,*[popt[0],popt[1]])))
        plt.plot(waves,np.log10(bb(waves,*[popt[2],popt[3]])))
        plt.show()
        exit()


    


        '''
        try:
            #popt, pcov = curve_fit(bb,waves,phot,p0=(phot[0],1e-2))
            #popt, pcov = curve_fit(bb2,waves,phot,p0=(phot[0],1e-2,phot[4],1e-3))
            #popt, pcov = curve_fit(bblog2,waves,np.log10(phot),p0=(np.log10(phot[0]),-2.0,np.log10(phot[4]),-3.0))
            idx1 = range(0,5)
            idx2 = range(5,len(waves))
            popt1,pcov1 = curve_fit(bb,waves[idx1],phot[idx1])
            popt2,pcov2 = curve_fit(bb,waves[idx2],phot[idx2])
            
        except:
            continue
        '''
        print 'Fitting %s' % star['ID']

        popt1,pcov1 = curve_fit(bb,waves[0:5],phot[0:5])
        popt2,pcov2 = curve_fit(bb,waves[7:],phot[7:])
        #
        print popt1,popt2
        #fit = bblog2(waves,*popt)
        fit1 = bb(waves,*popt1)
        fit2 = bb(waves,*popt2)
        fig = plt.figure()
        plt.plot(waves,np.log10(phot),'ro')
        plt.plot(waves,np.log10(fit1+fit2),'b--',linewidth=2)
        plt.plot(waves,np.log10(fit1),'b-',linewidth=2)
        plt.plot(waves,np.log10(fit2),'m-',linewidth=2)
        #plt.plot(waves,np.log10(bb(waves,popt[0],popt[1])),'b-',linewidth=2)
        #plt.plot(waves,np.log10(bb(waves,popt[2],popt[3])),'b-',linewidth=2)
        plt.show()
        exit()
        #

        
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

    #pp.close()
    #plt.show()
    exit()
        


def main():
    parser = argparse.ArgumentParser(description='Plot and fit of SEDs of input photometry')

    parser.add_argument('table',type=str,help='FITS table of data')

    args = parser.parse_args()

    plot_SED(args.table)



if __name__ == '__main__':
    main()

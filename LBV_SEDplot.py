#! /usr/bin/env python
from astropy.table import Table
import matplotlib.pyplot as plt
import argparse

def plot_SED(tableFile):
    table = Table.read(tableFile)

    waves = [float(x[6:10]) for x in table.columns[1:]]
    for star in table:
        phot = [x for x in star.data][1:]
        plt.plot(waves,phot,'ro')

        plt.xlabel(r'$\lambda [\micron]$')
        plt.ylabel(r'$\lambda\,F_{\lambda} [erg sec^{-1} cm^{-2}]$')
        plt.show()
        exit()
        


def main():
    parser = argparse.ArgumentParser(description='Plot and fit of SEDs of input photometry')

    parser.add_argument('table',type=str,help='FITS table of data')

    args = parser.parse_args()

    plot_SED(args.table)



if __name__ == '__main__':
    main()

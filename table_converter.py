#! /usr/bin/env python
import numpy as np
from astropy.table import Table,Column


WeisFile = "Weis.dat"
WeisDat = np.genfromtxt(WeisFile,names=['Name','RA','DEC','U','B','V','R','I','priority'],dtype=['a10','f8','f8','f8','f8','f8','f8','f8','i4'],autostrip=True)
Wtab = Table(WeisDat,meta={'title':'w','include':False,'num':len(WeisDat)})
colorBV = Column(name='B_V',data=Wtab['B']-Wtab['V'])
colorUB = Column(name='U_B',data=Wtab['U']-Wtab['B'])
Wtab.add_column(colorBV)
Wtab.add_column(colorUB)
Wtab.write("Weis.fits")

ValeevFile = "Valeev_altered2.dat"
ValeevDat = np.genfromtxt(ValeevFile,skip_header=1,delimiter=';',names=['Name','RA','DEC','V','B_V','s','R','comment'],dtype=['a6','f8','f8','f8','f8','f8','f8','a100'],autostrip=True)
Vtab = Table(ValeevDat,meta={'title':'v','include':False,'num':len(ValeevDat)})
Vtab.write("Valeev.fits")

HumphreysFile = "Humphreys.dat"
HumphreysDat = np.genfromtxt(HumphreysFile,delimiter=';',names=['Gal','Name','RA','DEC','Spec','comment'],dtype=['a3','a10','f8','f8','a20','a100'],autostrip=True)
Htab = Table(HumphreysDat,meta={'title':'h','include':False,'num':len(HumphreysDat)})
Htab.write("Humphreys.fits")

MasseyFile = "Massey17_final.tsv"
MasseyDat = np.genfromtxt(MasseyFile,delimiter=';',names=['LGGS','Gal','Name','V','B_V','U_B','SpType','CrossID','Ref','RA','DEC'],dtype=['a20','a4','a20','a20','a20','a20','a20','a20','a20','f8','f8'],autostrip=True)
Mtab = Table(MasseyDat,meta={'title':'m','include':False,'num':len(MasseyDat)})
Mtab.write("Massey.fits")

MasseyLongFile = "MasseyLong.tsv"
MasseyLongDat = np.genfromtxt(MasseyLongFile,skip_header=1,delimiter=';',names=['Gal','LGGS','n_LGGS','HAmag','HA','O3','V','B_V','U_B','Q','Type','Name','Ref','RA','DEC'],dtype=['a8','a20','a1','f8','f8','f8','f8','f8','f8','f8','a20','a20','a5','f8','f8'],autostrip=True)
MLtab = Table(MasseyLongDat,meta={'title':'mL','include':False,'num':len(MasseyLongDat)})
MLtab.write("MasseyLong.fits")

McQuinnNV = "McQuinn_nonvariable.txt"
McQuinnNVDat = np.genfromtxt(McQuinnNV,names=['Gal','ID','RA','DEC'],dtype=['a20','a20','f8','f8'],usecols=(0,1,2,3),delimiter=';',autostrip=True)
for idx,val in enumerate(McQuinnNVDat['ID']):
    McQuinnNVDat['RA'][idx] = np.float(val[1:10])
    McQuinnNVDat['DEC'][idx] = np.float(val[11:19])
McQNVtab = Table(McQuinnNVDat,meta={'title':'mqnv','include':False,'num':len(McQuinnNVDat)})
McQNVtab.write("McQuinn_NV.fits")


McQuinnV = "McQuinn_variable.txt"
McQuinnVDat = np.genfromtxt(McQuinnV,names=['Gal','ID','RA','DEC'],dtype=['a20','a20','f8','f8'],usecols=(0,1,2,3),delimiter=';',autostrip=True)
for idx,val in enumerate(McQuinnVDat['ID']):
    McQuinnVDat['RA'][idx] = np.float(val[1:10])
    McQuinnVDat['DEC'][idx] = np.float(val[11:19])
McQVtab = Table(McQuinnVDat,meta={'title':'mqv','include':False,'num':len(McQuinnVDat)})
McQVtab.write("McQuinn_V.fits")

Drout31 = "Drout_M31_SG.tsv"
Drout31data = np.genfromtxt(Drout31,delimiter=';',names=['recno','Name','RA','DEC','rank','V','fV','logT','flogT','logL','flogL','comment'],dtype=['i8','a20','f8','f8','i8','f8','a1','f8','a1','f8','a1','a31'],autostrip=True)
D31tab = Table(Drout31data,meta={'title':'d31','include':True,'num':len(Drout31data)})
D31tab.write("Drout_M31_SG.fits")

Drout33_r = "Drout_M33_redSG.txt"
Drout33_r_data = np.genfromtxt(Drout33_r,delimiter=';',names=['Name','RA','DEC','rank','Teff','Lum','Var'],dtype=['a19','f8','f8','a1','f8','f8','a5'],autostrip=True)
D33tab_r = Table(Drout33_r_data,meta={'title':'d33_r','include':True,'num':len(Drout33_r_data)})
D33tab_r.write("Drout_M33_redSG.fits")

Drout33_y = "Drout_M33_yellowSG.txt"
Drout33_y_data = np.genfromtxt(Drout33_y,delimiter=';',names=['Name','RA','DEC','rank','Teff','Lum','Var'],dtype=['a19','f8','f8','a1','f8','f8','a5'],autostrip=True)
D33tab_y = Table(Drout33_y_data,meta={'title':'d33_y','include':True,'num':len(Drout33_y_data)})
D33tab_y.write("Drout_M33_yellowSG.fits")

MasseyXLongFile = "MASSEYXLFINAL.dat"
MasseyXLongDat = np.genfromtxt(MasseyXLongFile,skip_header=1,delimiter=';',names=['Gal','LGGS','n_LGGS','Name','V','B_V','U_B','V_R','R_I','Type','rType','RA','DEC'],dtype=['a5','a20','a1','a20','f8','f8','f8','f8','f8','a10','a5','f8','f8'],autostrip=True)
MXLtab = Table(MasseyXLongDat,meta={'title':'mXL','include':False,'num':len(MasseyXLongDat)})
MXLtab.write("MasseyXL.fit")

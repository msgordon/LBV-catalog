#! /usr/bin/env python
import numpy as np
from astropy.table import Table,Column
from astropy.coordinates import ICRS
from astropy import units as u
from star import Star


'''
MouldFile = 'tables/table1a'
MouldDat = np.genfromtxt(MouldFile,dtype=['a1','a30','a20','a20','a20','a20','a20','a20','a20','a20','a20','a20','a2'],autostrip=True)
newMould = []
for idx,row in enumerate(MouldDat):
    new = [x[0:-1] if '&' in x else x for x in row]
    new = new[1:-1]
    new.insert(1,new[0].split('+')[1])
    new[0] = new[0].split('+')[0]
    new.insert(0,''.join(['J00',new[0],'+',new[1]]))
    new[1:] = [float(x) for x in new[1:]]
    
    newMould.append(new)


newMould = np.core.records.array(newMould,dtype=[('ID','S19'),('RA',float),('DEC',float),('RAd',float),('DECd',float),('__3_6_',float),('e36',int),('__4_5_',float),('e45',int),('__5_8_',float),('e58',int),('__8_0_',float),('e80',int)])
#MouldDat=np.array(MouldDat,dtype=[('RA',float),('DEC',float),('RAd',float),('DECd',float),('__3_6_',float),('e36',int),('__4_5_',float),('e45',int),('__5_8_',float),('e58',int),('__8_0_',float),('e80',int)])

MTab = Table(newMould,meta={'title':'jm','include':False,'num':len(newMould)})
#print MTab['__3_6_']
MTab.write('tables/Mould.fits')
exit()
'''

'''
MouldFile = 'tables/m31.bigrot2'
MouldDat = np.genfromtxt(MouldFile,names=['Num','RAd','DECd','J','H','K','__3_6_','__4_5_','__5_8','__8_0_'])

MTab = Table(MouldDat,meta={'title':'jm','include':False,'num':len(MouldDat)})

mRA = []
mDEC = []
for decimal in zip(MouldDat['RAd'],MouldDat['DECd']):
    c = ICRS(ra=decimal[0],dec=decimal[1],unit=(u.degree,u.degree))
    RA,DEC = c.to_string(sep='').split()
    mRA.append(np.float(RA))
    mDEC.append(np.float(DEC))
cRA = Column(name='RA',data=mRA)
cDEC = Column(name='DEC',data=mDEC)
MTab.add_columns([cRA,cDEC])
MTab.write('tables/Mould.fits')

exit()
'''
'''
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

MasseyLongFile = "tables/MasseyLong.tsv"
MasseyLongDat = np.genfromtxt(MasseyLongFile,skip_header=1,delimiter=';',names=['Gal','LGGS','n_LGGS','HAmag','HA','O3','V','B_V','U_B','Q','Type','Name','Ref','RA','DEC'],dtype=['a8','a20','a1','f8','f8','f8','f8','f8','f8','f8','a20','a20','a5','f8','f8'],autostrip=True)
ra,dec = zip(*[Star.sex2deg(ra,dec) for ra,dec in zip(MasseyLongDat['RA'],MasseyLongDat['DEC'])])

MLtab = Table(MasseyLongDat,meta={'title':'mL','include':False,'num':len(MasseyLongDat)})
MLtab.add_column(Column(ra,name='RAd'))
MLtab.add_column(Column(dec,name='DECd'))

MLtab.write("tables/MasseyLong.fits")

exit()
'''
'''
McQuinnNV = "tables/McQuinn_nonvariable.txt"
McQuinnNVDat = np.genfromtxt(McQuinnNV,names=['Gal','ID','RA','DEC'],dtype=['a20','a20','f8','f8'],usecols=(0,1,2,3),delimiter=';',autostrip=True)
for idx,val in enumerate(McQuinnNVDat['ID']):
    McQuinnNVDat['RA'][idx] = np.float(val[1:10])
    McQuinnNVDat['DEC'][idx] = np.float(val[11:19])

ra,dec = zip(*[Star.sex2deg(ra,dec) for ra,dec in zip(McQuinnNVDat['RA'],McQuinnNVDat['DEC'])])
McQNVtab = Table(McQuinnNVDat,meta={'title':'mqnv','include':False,'num':len(McQuinnNVDat)})
McQNVtab.add_column(Column(ra,name='RAd'))
McQNVtab.add_column(Column(dec,name='DECd'))


McQNVtab.write("tables/McQuinn_NV.fits")


McQuinnV = "tables/McQuinn_variable.txt"
McQuinnVDat = np.genfromtxt(McQuinnV,names=['Gal','ID','RA','DEC'],dtype=['a20','a20','f8','f8'],usecols=(0,1,2,3),delimiter=';',autostrip=True)
for idx,val in enumerate(McQuinnVDat['ID']):
    McQuinnVDat['RA'][idx] = np.float(val[1:10])
    McQuinnVDat['DEC'][idx] = np.float(val[11:19])
McQVtab = Table(McQuinnVDat,meta={'title':'mqv','include':False,'num':len(McQuinnVDat)})
ra,dec = zip(*[Star.sex2deg(ra,dec) for ra,dec in zip(McQuinnVDat['RA'],McQuinnVDat['DEC'])])
McQVtab.add_column(Column(ra,name='RAd'))
McQVtab.add_column(Column(dec,name='DECd'))
McQVtab.write("tables/McQuinn_V.fits")
'''
'''
Drout31 = "tables/Drout_M31_SG.tsv"
Drout31data = np.genfromtxt(Drout31,delimiter=';',names=['recno','Name','RA','DEC','rank','V','fV','logT','flogT','logL','flogL','comment'],dtype=['i8','a20','f8','f8','i8','f8','a1','f8','a1','f8','a1','a31'],autostrip=True)

D31tab = Table(Drout31data,meta={'title':'d31','include':True,'num':len(Drout31data)})
ra,dec = zip(*[Star.sex2deg(ra,dec) for ra,dec in zip(Drout31data['RA'],Drout31data['DEC'])])
D31tab.add_column(Column(ra,name='RAd'))
D31tab.add_column(Column(dec,name='DECd'))
D31tab.write("tables/Drout_M31_SG.fits")
exit()

Drout33_r = "Drout_M33_redSG.txt"
Drout33_r_data = np.genfromtxt(Drout33_r,delimiter=';',names=['Name','RA','DEC','rank','Teff','Lum','Var'],dtype=['a19','f8','f8','a1','f8','f8','a5'],autostrip=True)
D33tab_r = Table(Drout33_r_data,meta={'title':'d33_r','include':True,'num':len(Drout33_r_data)})
D33tab_r.write("Drout_M33_redSG.fits")

Drout33_y = "Drout_M33_yellowSG.txt"
Drout33_y_data = np.genfromtxt(Drout33_y,delimiter=';',names=['Name','RA','DEC','rank','Teff','Lum','Var'],dtype=['a19','f8','f8','a1','f8','f8','a5'],autostrip=True)
D33tab_y = Table(Drout33_y_data,meta={'title':'d33_y','include':True,'num':len(Drout33_y_data)})
D33tab_y.write("Drout_M33_yellowSG.fits")
'''
'''
MasseyXLongFile = "tables/MASSEYXLFINAL.dat"
MasseyXLongDat = np.genfromtxt(MasseyXLongFile,skip_header=1,delimiter=';',names=['Gal','LGGS','n_LGGS','Name','V','B_V','U_B','V_R','R_I','Type','rType','RA','DEC'],dtype=['a5','a20','a1','a20','f8','f8','f8','f8','f8','a10','a5','f8','f8'],autostrip=True)
MXLtab = Table(MasseyXLongDat,meta={'title':'mXL','include':False,'num':len(MasseyXLongDat)})
ra,dec = zip(*[Star.sex2deg(ra,dec) for ra,dec in zip(MXLtab['RA'],MXLtab['DEC'])])
MXLtab.add_column(Column(ra,name='RAd'))
MXLtab.add_column(Column(dec,name='DECd'))
MXLtab.write("tables/MasseyXL.fit")
'''
McQuinn = 'tables/McQuinn_dl.fit'
McQuinnDat = Table.read(McQuinn)
RA = [np.float(x.split('+')[0][1:]) for x in McQuinnDat['SSTM3307']]
DEC = [np.float(x.split('+')[1]) for x in McQuinnDat['SSTM3307']]
RAc = Column(RA,name='RA',dtype=np.float)
DECc = Column(DEC,name='DEC',dtype=np.float)

McQuinnDat.add_column(RAc)
McQuinnDat.add_column(DECc)

McQuinnDat.meta['title'] = 'mcqS'
McQuinnDat.meta['include'] = False
McQuinnDat.meta['num'] = len(McQuinnDat)
ra,dec = zip(*[Star.sex2deg(ra,dec) for ra,dec in zip(McQuinnDat['RA'],McQuinnDat['DEC'])])
McQuinnDat.add_column(Column(ra,name='RAd'))
McQuinnDat.add_column(Column(dec,name='DECd'))
McQuinnDat.write('tables/McQuinn.fits')




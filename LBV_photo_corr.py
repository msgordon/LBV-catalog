#! /usr/bin/env python
import numpy as np
from astropy.table import Table, Column,vstack

phot_table = Table.read('photometry.fits')
# replace broke ass values
# J013337.00+303637.5 -> J013337.04+303637.6
# J013415.38+302816.3 -> J013415.42+302816.4 
id1 = np.where(phot_table['ID'] == 'J013337.00+303637.5')
id2 = np.where(phot_table['ID'] == 'J013415.38+302816.3')
phot_table['ID'][id1] = 'J013337.04+303637.6'
phot_table['ID'][id2] = 'J013415.42+302816.4'

M31_table = Table.read('tables/MasseyXL_M31_photo_err.fit')
M33_table = Table.read('tables/MasseyXL_M33_photo_err.fit')

idx = [np.where(M31_table['LGGS'] == star['ID']) for star in phot_table if star['ID'] in M31_table['LGGS']]
rows31 = M31_table[idx]
idx = [np.where(M33_table['LGGS'] == star['ID']) for star in phot_table if star['ID'] in M33_table['LGGS']]
rows33 = M33_table[idx]

rows = vstack([rows31, rows33])

e_Vmag = Column(data=np.ndarray.flatten(rows['e_Vmag']),name='e_Vmag',dtype=np.float)#**2-rows['e_Vmag']**2))
e_Bmag = Column(data=np.ndarray.flatten(rows['e_B-V']),name='e_Bmag',dtype=np.float)
e_Umag = Column(data=np.ndarray.flatten(rows['e_U-B']),name='e_Umag',dtype=np.float)
e_Rmag = Column(data=np.ndarray.flatten(rows['e_V-R']),name='e_Rmag',dtype=np.float)
e_Imag = Column(data=np.ndarray.flatten(rows['e_R-I']),name='e_Imag',dtype=np.float)

e_Vflux = Column(data=np.log(10.)/2.5*np.ndarray.flatten(np.array([x*y for x,y in zip(phot_table['F_V_Jy'],e_Vmag)])),name='e_F_V_Jy')
e_Bflux = Column(data=np.log(10.)/2.5*np.ndarray.flatten(np.array([x*y for x,y in zip(phot_table['F_B_Jy'],e_Bmag)])),name='e_F_B_Jy')
e_Uflux = Column(data=np.log(10.)/2.5*np.ndarray.flatten(np.array([x*y for x,y in zip(phot_table['F_U_Jy'],e_Umag)])),name='e_F_U_Jy')
e_Rflux = Column(data=np.log(10.)/2.5*np.ndarray.flatten(np.array([x*y for x,y in zip(phot_table['F_R_Jy'],e_Rmag)])),name='e_F_R_Jy')
e_Iflux = Column(data=np.log(10.)/2.5*np.ndarray.flatten(np.array([x*y for x,y in zip(phot_table['F_I_Jy'],e_Imag)])),name='e_F_I_Jy')


e_Vlamflux = Column(data=np.ndarray.flatten(np.array([x for x in e_Vflux]))*3.0e-9/0.55**2,name='e_F_0.55_um')
e_Blamflux = Column(data=np.ndarray.flatten(np.array([x for x in e_Bflux]))*3.0e-9/0.44**2,name='e_F_0.44_um')
e_Ulamflux = Column(data=np.ndarray.flatten(np.array([x for x in e_Uflux]))*3.0e-9/0.36**2,name='e_F_0.36_um')
e_Rlamflux = Column(data=np.ndarray.flatten(np.array([x for x in e_Rflux]))*3.0e-9/0.71**2,name='e_F_0.71_um')
e_Ilamflux = Column(data=np.ndarray.flatten(np.array([x for x in e_Iflux]))*3.0e-9/0.97**2,name='e_F_0.97_um')

e_Ve = Column(data=np.ndarray.flatten(np.array([x for x in e_Vlamflux]))*0.55,name='e_lam_F_0.55_um')
e_Be = Column(data=np.ndarray.flatten(np.array([x for x in e_Blamflux]))*0.44,name='e_lam_F_0.44_um')
e_Ue = Column(data=np.ndarray.flatten(np.array([x for x in e_Ulamflux]))*0.36,name='e_lam_F_0.36_um')
e_Re = Column(data=np.ndarray.flatten(np.array([x for x in e_Rlamflux]))*0.71,name='e_lam_F_0.71_um')
e_Ie = Column(data=np.ndarray.flatten(np.array([x for x in e_Ilamflux]))*0.97,name='e_lam_F_0.97_um')


phot_table.add_columns([e_Umag,e_Bmag,e_Vmag,e_Rmag,e_Imag,e_Uflux,e_Bflux,e_Vflux,e_Rflux,e_Iflux,e_Ulamflux,e_Blamflux,e_Vlamflux,e_Rlamflux,e_Ilamflux,e_Ue,e_Be,e_Ve,e_Re,e_Ie])

phot_table.write('photometry_corr.fits')
phot_table.write('photometry_corr.tsv',format='ascii.tab')


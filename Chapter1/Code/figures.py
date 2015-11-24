# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 12:00:30 2015

@author: nabobalis
"""

import numpy as np
import matplotlib.pyplot as plt
from sunkitsst import read_cubes
import astropy
from astropy.table import Table
from astropy import units as u
from astropy.constants import k_B, m_p

"""
Figure of the photopshere/QS
"""

imfile = '/data/Mounted/SWAT/monster/series/crispex.6302.fullstokes_aligned.icube'
spfile = '/data/Mounted/SWAT/monster/series/crispex.6302.fullstokes_aligned.sp.icube'
data = read_cubes.read_cubes(imfile,spfile)[1]
res = 0.059#*725
i = 30
extent = [0,data.shape[2]*res,0,data.shape[3]*res]
plt.figure()
plt.imshow(data[i,0],origin='lower',cmap=plt.get_cmap('gray'),
           extent = extent, vmin=3.5*np.std(data[i,0]))
plt.xlabel('Distance (arcsecond)')
plt.ylabel('Distance (arcsecond)')

circle1=plt.Circle((48.1,27),2,color='white',fill=False)
plt.annotate('A', (47,21), fontsize=24, color='w')

circle2=plt.Circle((31.4,26),2,color='white',fill=False)
plt.annotate('B', (30.5,20), fontsize=24, color='w')

fig = plt.gcf()
fig.gca().add_artist(circle1)
fig.gca().add_artist(circle2)

plt.savefig('/home/nabobalis/GitRepos/PhDThesis/Chapter1/Figs/QS.pdf',dpi=300,bbox_inches='tight')

"""
Figure of the VAL Model Atmosphere.

"""

#def read_VAL3c_MTW(VAL_file=None, MTW_file=None, mu=0.602):
#    """
#    Read in the data from Table 12 in Vernazza (1981) and combine with
#    McWhirter (1975).
#    Parameters
#    ----------
#    VAL_file : string
#        The data file for the VAL3c atmosphere, defaults to
#        `pysac.mhs_atmosphere.hs_model.VALIIIc_data`
#    MTW_file : string
#        The data file for the McWhirter atmosphere, defaults to
#        `pysac.mhs_atmosphere.hs_model.MTWcorona_data`
#    mu : float
#        The mean molecular weight ratio for the corona. defaults to 0.6.
#    Returns
#    -------
#    data : `astropy.table.Table`
#        The combined data, sorted by Z.
#    """
#    VAL3c = Table.read(VAL_file, format='ascii', comment='#')
#    VAL3c['Z'].unit = u.km
#    VAL3c['rho'].unit = u.Unit('g cm-3')
#    VAL3c['p'].unit = u.Unit('dyne/cm^2')
#    VAL3c['T'].unit = u.K
#    VAL3c['n_i'].unit = u.one/u.cm**3
#    VAL3c['n_e'].unit = u.one/u.cm**3
#    # Calculate the mean molecular weight ratio
#    VAL3c['mu'] = 4.0/(3*0.74+1+VAL3c['n_e']/VAL3c['n_i'])
#    VAL3c.sort('Z')
#
#    return VAL3c
#
#data = read_VAL3c_MTW(VAL_file='/home/nabobalis/GitRepos/PhDThesis/Chapter1/Code/VALIIIC.dat')
#
#fig, ax = plt.subplots()
#axes = [ax, ax.twinx()]
#lns1 = axes[0].semilogy(data['Z'],data['rho'],label=r'$\rho$', color='blue')
#lns2 = axes[1].semilogy(data['Z'],data['T'],label=r'$T$', color='red')
#lns = lns1+lns2
#labs = [l.get_label() for l in lns]
#
#axes[1].legend(lns, labs)
#axes[0].set_xlabel(r'Height (km)')
#axes[0].set_ylabel(r'Density (g / cm$^3$)', color='blue')
#axes[0].set_xlim(0,2350)
#axes[0].tick_params(axis='y', colors='blue')
#axes[1].set_ylabel(r'Temperture (K)', color='red')
#axes[1].tick_params(axis='y', colors='red')
#
#plt.savefig('/home/nabobalis/GitRepos/PhDThesis/Chapter1/Figs/val3.pdf',dpi=300,bbox_inches='tight')

"""
Figure of the AR
"""

#imfile = '/data/Mounted/SWAT/arlimb/crispex.6563.icube'
#spfile = '/data/Mounted/SWAT/arlimb/crispex.6563.sp.icube'
#data = read_cubes.read_cubes(imfile,spfile)[1]
#res = 0.059#*725
#i = -1
#extent = [0,data.shape[2]*res,0,data.shape[3]*res]
#plt.figure()
#plt.imshow(data[i,0],origin='lower',cmap=plt.get_cmap('gray'),
#           extent = extent)
#plt.xlabel('Distance (arcsecond)')
#plt.ylabel('Distance (arcsecond)')
#
#circle1=plt.Circle((5.7,46.65),10,color='white',fill=False)
#plt.annotate('A', (5,33), fontsize=24, color='w')
#
#circle2=plt.Circle((32.8,11.9),7,color='white',fill=False)
#plt.annotate('D', (31.8,1), fontsize=24, color='w')
#
#circle3=plt.Circle((34,37.7),2,color='white',fill=False)
#plt.annotate('B', (32.7,32), fontsize=24, color='w')
#
#circle4=plt.Circle((49.5,31.9),2,color='white',fill=False)
#plt.annotate('C', (48.5,26), fontsize=24, color='w')
#
#fig = plt.gcf()
#fig.gca().add_artist(circle1)
#fig.gca().add_artist(circle2)
#fig.gca().add_artist(circle3)
#fig.gca().add_artist(circle4)
#
#plt.savefig('/home/nabobalis/GitRepos/PhDThesis/Chapter1/Figs/AR.pdf',dpi=300,bbox_inches='tight')

"""
Figure of the SunSpot Number (Monthly)
"""

#mon = np.loadtxt('ISSN_M_tot.csv', delimiter=',')
#plt.figure()
#plt.plot(mon[:,0],mon[:,3], 'b', label='Monthly')
#plt.plot(mon[:,0],mon[:,5], 'r', label='Running Average')
#plt.xlabel('Time (years)')
#plt.ylabel('Number')
#plt.xlim(mon[:,0].min(),mon[:,0].max())
#plt.ylim(0,280)
#plt.legend()
#plt.savefig('/home/nabobalis/GitRepos/PhDThesis/Chapter1/Figs/sunspot_number.pdf',dpi=300,bbox_inches='tight')

"""
Figure of the Chromosphere
"""

#imfile = '/data/Mounted/SWAT/fastrbe/crisp.6563.icube'
#spfile = '/data/Mounted/SWAT/fastrbe/crisp.6563.sp.icube'
#data = read_cubes.read_cubes(imfile,spfile)[1]
#res = 0.059#*725
#i = 350
#extent = [0,data.shape[2]*res,0,data.shape[3]*res]
#plt.figure()
#plt.imshow(data[i,2],origin='lower',cmap=plt.get_cmap('gray'),
#           extent = extent)
#plt.xlabel('Distance (arcsecond)')
#plt.ylabel('Distance (arcsecond)')
#
#circle1=plt.Circle((12,33),10,color='white',fill=False)
#plt.annotate('A', (11,18), fontsize=24, color='w')
#
#circle2=plt.Circle((36,38.5),6,color='white',fill=False)
#plt.annotate('B', (34.8,28.5), fontsize=24, color='w')
#
#circle3=plt.Circle((42.3,29.9),4,color='white',fill=False)
#plt.annotate('C', (41.3,22), fontsize=24, color='w')
#
#fig = plt.gcf()
#fig.gca().add_artist(circle1)
#fig.gca().add_artist(circle2)
#fig.gca().add_artist(circle3)
#
#plt.savefig('/home/nabobalis/GitRepos/PhDThesis/Chapter1/Figs/Chromo.pdf',dpi=300,bbox_inches='tight')
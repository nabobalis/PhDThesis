# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 12:00:30 2015

@author: nabobalis
"""

from __future__ import division
import pycwt as wavelet
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sunpy.io.fits as fit
from scipy.io import readsav
import scipy.fftpack as fftpack
import matplotlib.animation as anim
import glob
from sunkitsst.read_cubes import read_cubes
from matplotlib.image import NonUniformImage
import sunpy.map as mapc
from scipy.stats import pearsonr
from scipy.signal import detrend

def area_inten(bound,lim_inten):

    area = np.zeros([bound.shape[0]])
    inten = np.zeros([bound.shape[0]])
    pore = np.zeros(bound.shape,dtype=np.int)

    for i in range(0,bound.shape[0]):
        pore[i] = (bound[i] <= lim_inten[i])
        area[i] = len(pore[i].nonzero()[0])
        inten[i] = np.sum(bound[i][pore[i].nonzero()])

    return area, inten

def find_closest(array, target):

    idx = array.searchsorted(target) #array must be sorted
    idx = np.clip(idx, 1, len(array)-1)
    left = array[idx-1]
    right = array[idx]
    idx -= target - left < right - target
    return idx

def wavel(signal,cadence):
    mother='morlet'
    sig_level = 0.95
    #/ signal.std()
    t1 = np.linspace(0,cadence*signal.size,signal.size)
    wave, scales, freqs, coi = wavelet.cwt((signal - signal.mean()),cadence,
                                           wavelet=mother, dj=1/100.)

    power = (np.abs(wave)) ** 2 # / scales[:,None]
    period = 1/freqs
    alpha = 0.0
    ## (variance=1 for the normalized SST)
    signif, fft_theor = wavelet.significance(signal, period, scales, 0, alpha,
                            significance_level=sig_level, wavelet=mother)
    sig95 = np.ones([1, signal.size]) * signif[:, None]
    sig95 = power / sig95

    ## indices for stuff
    idx = find_closest(period,coi.max())

    ## Into minutes
    t1 /= 60
    period /= 60
    coi /= 60

    return wave,scales,sig95,idx,t1,coi,period,power

def cross_wavelet(signal_1, signal_2, period, mother='morlet', plot=True):

    signal_1 = (signal_1 - signal_1.mean()) / signal_1.std()    # Normalizing
    signal_2 = (signal_2 - signal_2.mean()) / signal_2.std()    # Normalizing

    W12, cross_coi, freq, signif = wavelet.xwt(signal_1, signal_2, period, dj=1/100, s0=-1, J=-1,
                                         significance_level=0.95, wavelet=mother,
                                         normalize=True)

    cross_power = np.abs(W12)**2
    cross_sig = np.ones([1, signal_1.size]) * signif[:, None]
    cross_sig = cross_power / cross_sig
    cross_period = 1/freq

    WCT, aWCT, corr_coi, freq, sig = wavelet.wct(signal_1, signal_2, period, dj=1/100, s0=-1, J=-1,
                                            sig=False,significance_level=0.95, wavelet=mother,
                                            normalize=True)

    cor_sig = np.ones([1, signal_1.size]) * sig[:, None]
    cor_sig = np.abs(WCT) / cor_sig
    cor_period = 1/freq


    t1 = np.linspace(0,period*signal_1.size,signal_1.size)

    ## indices for stuff
    idx = find_closest(cor_period,corr_coi.max())

    ## Into minutes
    t1 /= 60
    cross_period /= 60
    cor_period /= 60
    cross_coi /= 60
    corr_coi /= 60

    extent_corr =  [t1.min(),t1.max(),0,max(cor_period)]

    plt.figure(figsize=(12,12))
    im3= plt.imshow(np.rad2deg(aWCT), origin='lower',interpolation='nearest', cmap='RdBu', extent=extent_corr)
    plt.fill(np.concatenate([t1, t1[-1:]+period, t1[-1:]+period,t1[:1]-period, t1[:1]-period]),
            (np.concatenate([corr_coi,[1e-9], cor_period[-1:], cor_period[-1:], [1e-9]])),
            'k', alpha=0.3,hatch='x')
    plt.ylim(([min(cor_period), cor_period[idx]]))
    plt.xlim(t1.min(),t1.max())
    cbar = plt.colorbar(im3)
    cbar.solids.set_edgecolor("face")
    plt.show()

"""
Figure of data analysis example
"""

### Creation of Signal and Saved
#time = np.linspace(0,1000,1000)
#data = np.sin(2*np.pi*time/5) + np.cos(2*np.pi*time/10) + 5*np.random.rand(1000)
#np.savetxt('/home/nabobalis/GitRepos/PhDThesis/Chapter2/Code/test_sig.txt',[time,data])

load = np.loadtxt('/home/nabobalis/GitRepos/PhDThesis/Chapter2/Code/test_sig.txt')
time = load[0]
data = load[1]
dt = time[1] - time[0]

# FFT
fft = fftpack.fft(data-np.mean(data))/time.shape[0]

freq = fftpack.fftfreq(data.shape[0],1)
fft_power = np.abs(fft)**2

# Wavelet
wave,scales,sig95,idx,t1,coi,period,power = wavel(data-data.mean(),dt)

# EMD

emd_data = np.loadtxt('/home/nabobalis/GitRepos/PhDThesis/Chapter2/Code/emd_sig.csv', delimiter=',')

# Plotting

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(figsize=(12,8),nrows=2, ncols=2)

ax1.plot(time/60, data, 'k', linewidth=1.5)
ax1.set_xlabel('Time (minutes)')
ax1.set_ylabel('Amplitude (No Units)')
ax1.set_xlim(time.min()/60,time.max()/60)

ax2.plot(freq[0:len(time)/2],fft_power[0:len(time)/2])
ax2.set_xlabel('Frequency (Hertz)')
ax2.set_ylabel('Power (No Units)')
ax2.axis((0.0, 0.5, 0.0, 0.30000000000000004))

extent = [time.min(),time.max(),0,max(period)]
im = NonUniformImage(ax3, interpolation='nearest', extent=extent)
im.set_cmap('cubehelix_r')
im.set_data(t1, period[:idx], power[:idx,:])
ax3.images.append(im)
ax3.set_ylabel('Period (minutes)')
ax3.set_xlabel('Time (minutes)')
ax3.contour(t1, period[:idx], sig95[:idx,:], [-99,1], colors='k', linewidths=2, extent=extent)
ax3.fill(np.concatenate([t1, t1[-1:]+dt, t1[-1:]+dt,t1[:1]-dt, t1[:1]-dt]),
        (np.concatenate([coi,[1e-9], period[-1:], period[-1:], [1e-9]])),
        'k', alpha=0.3,hatch='/', zorder=100, antialiased=True, rasterized=True)
ax3.set_xlim(t1.min(),t1.max())
ax3.set_ylim(([min(period), period[idx]]))
ax3.axis((0.0, 16.673398306304815, 0.033366700033366704, 1.2766635538646571))

ax4.plot(time/60, emd_data[:,3], 'k', linewidth=1.5)
ax4.set_xlabel('Time (minutes)')
ax4.set_ylabel('Amplitude (No Units)')
ax4.set_xlim(time.min()/60,time.max()/60)

fig.tight_layout()
plt.savefig('/home/nabobalis/GitRepos/PhDThesis/Chapter2/Figs/signal_overview.pdf',dpi=300,bbox_inches='tight')

"""
Data Sources
"""

## Sunspot
#sun_data = readsav('/home/nabobalis/Data/Reduced_data_blue.sav')['blue_spk_cor']
#sun_dt = 6.8
#sun_pixel_arc = 0.097
#sun_bounding_box = [318,636,386,694]
#sun_background_box = [300,700,100,400]
#sun__full_extent = [0,sun_pixel_arc*sun_data.shape[1],0,sun_pixel_arc*sun_data.shape[2]]
#sun_cut_extent = [sun_bounding_box[0]*sun_pixel_arc,sun_bounding_box[1]*sun_pixel_arc,
#                  sun_bounding_box[2]*sun_pixel_arc,sun_bounding_box[3]*sun_pixel_arc]
#
## Pore
#pore_data = fits.getdata('/home/nabobalis/Data/gband_pore.fits')
#pore_bounding_box = [450,700,550,750]
#pore_background_box = [300,550,100,350]
#pore_dt = 2.11
#pore_pixel_arc = 0.0968063872255489 # lol
#pore__full_extent = [0,pore_pixel_arc*pore_data.shape[1],0,pore_pixel_arc*pore_data.shape[2]]
#pore_cut_extent = [pore_bounding_box[0]*pore_pixel_arc,pore_bounding_box[1]*pore_pixel_arc,
#                   pore_bounding_box[2]*pore_pixel_arc,pore_bounding_box[3]*pore_pixel_arc]

"""
Overview plot of the waveguides.
"""

#fig, (ax1, ax2) = plt.subplots(1,2)
#
#ax1.imshow(sun_data[0], origin='lower',interpolation='nearest',extent=sun__full_extent, cmap=plt.get_cmap('afmhot'))
#ax1.set_xlabel('Distance (Arcseconds)')
#ax1.set_ylabel('Distance (Arcseconds)')
#ax1.axes.axis((16.827111807533818, 77.654734496225558, 24.186586938884577, 84.019074696063981))
#
#ax2.imshow(pore_data[0], origin='lower',interpolation='nearest',extent=pore__full_extent, cmap=plt.get_cmap('gray'))
#ax2.set_xlabel('Distance (Arcseconds)')
#ax2.set_ylabel('Distance (Arcseconds)')
#ax2.axes.axis((23.335643499752834, 89.270284752885019, 24.592770410202537, 91.649173045713155))
#
#fig.tight_layout()
#plt.savefig('/home/nabobalis/GitRepos/PhDThesis/Chapter2/Figs/overview.pdf',dpi=300,bbox_inches='tight')


#####################################
"""
This is need for the following parts!
"""
#####################################

#sun_bound = sun_data[:,sun_bounding_box[2]:sun_bounding_box[3],sun_bounding_box[0]:sun_bounding_box[1]]
#sun_cut_box = sun_data[:,sun_background_box[2]:sun_background_box[3],sun_background_box[0]:sun_background_box[1]]
#sun_cut = sun_cut_box.reshape(sun_cut_box.shape[0],sun_cut_box.shape[1]*sun_cut_box.shape[2])
#sun_time = np.linspace(0,sun_data.shape[0]*sun_dt,sun_data.shape[0])
#
#pore_bound = pore_data[:,pore_bounding_box[2]:pore_bounding_box[3],pore_bounding_box[0]:pore_bounding_box[1]]
#pore_cut_box = pore_data[:,pore_background_box[2]:pore_background_box[3],pore_background_box[0]:pore_background_box[1]]
#pore_cut = pore_cut_box.reshape(pore_cut_box.shape[0],pore_cut_box.shape[1]*pore_cut_box.shape[2])
#pore_time = np.linspace(0,pore_data.shape[0]*pore_dt,pore_data.shape[0])
#
#
#sun_lim_list = [3,3.5,4,4.5]
#pore_lim_list = [2,2.5,3,3.5]
#color = ['Blue','Green','Purple', 'Orange']
#sun_lim = []
#pore_lim = []
#sun_area = []
#pore_area = []
#sun_inten = []
#pore_inten = []
#
#for slim, plim in zip(sun_lim_list,pore_lim_list):
#    sun_lim_inten =  np.mean(sun_cut, axis = 1) - slim*np.std(sun_cut, axis = 1)
#    pore_lim_inten =  np.mean(pore_cut, axis = 1) - plim*np.std(pore_cut, axis = 1)
#    s_area, s_inten = area_inten(sun_bound,sun_lim_inten)
#    p_area, p_inten = area_inten(pore_bound,pore_lim_inten)
#    sun_lim.append(sun_lim_inten)
#    pore_lim.append(pore_lim_inten)
#    sun_area.append(s_area)
#    sun_inten.append(s_inten)
#    pore_area.append(p_area)
#    pore_inten.append(p_inten)
#
"""
Overview of Method
"""

#idx = 1
#fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(figsize=(12,8),nrows=2, ncols=2)
#
##Sunspot
#ax1.imshow(sun_bound[idx], origin = 'lower', interpolation = 'nearest', cmap=plt.get_cmap('gray'), extent=sun_cut_extent)
#ax1.set_xlabel('Distance (Arcseconds)')
#ax1.set_ylabel('Distance (Arcseconds)')
#
#for i,li in enumerate(sun_lim):
#    ax1.contour(sun_bound[idx] <= li[idx], origin = 'lower', interpolation = 'nearest', colors=color[i], extent=sun_cut_extent, levels=[0,1])
#
#ax2.hist([sun_bound[idx].flatten(),sun_cut_box[idx].flatten()],bins=250, color= ['Red', 'Orange'],
#          label=['Boundary', 'Background'], stacked=True, fill=True, edgecolor = 'none')
#ax2.set_xlim(1000,5000)
#
#for i,li in enumerate(sun_lim):
#    ax2.axvline(li[idx], color=color[i], linestyle='dashed', linewidth=2)
#
#ax2.set_xlabel('Intensity bins (counts)')
#ax2.set_ylabel('Frequency')
##ax2.locator_params(axis='x', nbins=6)
#ax2.legend()
#
##Pore
#ax3.imshow(pore_bound[idx], origin = 'lower', interpolation = 'nearest', cmap=plt.get_cmap('gray'), extent=pore_cut_extent)
#ax3.set_xlabel('Distance (Arcseconds)')
#ax3.set_ylabel('Distance (Arcseconds)')
#
#for i,li in enumerate(pore_lim):
#    ax3.contour(pore_bound[idx] <= li[idx], origin = 'lower', interpolation = 'nearest', colors=color[i], extent=pore_cut_extent, levels=[0,1])
#
#ax4.hist([pore_bound[idx].flatten(),pore_cut_box[idx].flatten()],bins=250, color= ['Red', 'Orange'],
#          label=['Boundary', 'Background'], stacked=True, fill=True, edgecolor = 'none')
#ax4.set_xlim(0.2,2)
#
#for i,li in enumerate(pore_lim):
#    ax4.axvline(li[idx], color=color[i], linestyle='dashed', linewidth=2)
#
#ax4.set_xlabel('Intensity bins (counts)')
#ax4.set_ylabel('Number')
#ax4.legend()
#
#fig.tight_layout()
#plt.savefig('/home/nabobalis/GitRepos/PhDThesis/Chapter2/Figs/method_overview.pdf',dpi=300,bbox_inches='tight')
#
"""
Wavelet of Signals
"""

##Sunspot
#wave,scales,sig95,idx,t1,coi,period,power = wavel(sun_area[0],sun_dt)
#wave2,scales2,sig952,idx2,t12,coi2,period2,power2 = wavel(sun_area[-1],sun_dt)

#fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(figsize=(12,8),nrows=2, ncols=2)
# First Signal
#ax1.plot(t1, sun_area[0], 'k', linewidth=1.5)
#ax1.set_xlabel('Time (minutes)')
#ax1.set_ylabel('Pixels')
#extent = [t1.min(),t1.max(),0,max(period)]
#im = NonUniformImage(ax3, interpolation='nearest', extent=extent)
#im.set_cmap('cubehelix_r')
#im.set_data(t1, period[:idx], power[:idx,:])
#ax3.images.append(im)
#ax3.set_ylabel('Period (minutes)')
#ax3.set_xlabel('Time (minutes)')
#ax3.contour(t1, period[:idx], sig95[:idx,:], [-99,1], colors='k', linewidths=2, extent=extent)
#ax3.fill(np.concatenate([t1, t1[-1:]+sun_dt, t1[-1:]+sun_dt,t1[:1]-sun_dt, t1[:1]-sun_dt]),
#        (np.concatenate([coi,[1e-9], period[-1:], period[-1:], [1e-9]])),
#        'k', alpha=0.3,hatch='/', zorder=100, antialiased=True, rasterized=True)
#ax3.set_xlim(t1.min(),t1.max())
#ax3.set_ylim(([min(period), period[idx]]))
## Second Signal
#ax2.plot(t1, sun_area[-1], 'k', linewidth=1.5)
#ax2.set_xlabel('Time (minutes)')
#ax2.set_ylabel('Pixels')
#extent = [t1.min(),t1.max(),0,max(period2)]
#im = NonUniformImage(ax4, interpolation='nearest', extent=extent)
#im.set_cmap('cubehelix_r')
#im.set_data(t1, period2[:idx], power2[:idx,:])
#ax4.images.append(im)
#ax4.set_ylabel('Period (minutes)')
#ax4.set_xlabel('Time (minutes)')
#ax4.contour(t1, period2[:idx], sig952[:idx,:], [-99,1], colors='k', linewidths=2, extent=extent)
#ax4.fill(np.concatenate([t1, t1[-1:]+sun_dt, t1[-1:]+sun_dt,t1[:1]-sun_dt, t1[:1]-sun_dt]),
#        (np.concatenate([coi,[1e-9], period[-1:], period[-1:], [1e-9]])),
#        'k', alpha=0.3,hatch='/', zorder=100, antialiased=True, rasterized=True)
#ax4.set_xlim(t1.min(),t1.max())
#ax4.set_ylim(([min(period), period[idx]]))
#
#fig.tight_layout()
#plt.savefig('/home/nabobalis/GitRepos/PhDThesis/Chapter2/Figs/sunspot_wavelet.pdf',dpi=300,bbox_inches='tight')

##Pore
#wave,scales,sig95,idx,t1,coi,period,power = wavel(pore_area[0],pore_dt)
#wave2,scales2,sig952,idx2,t12,coi2,period2,power2 = wavel(pore_area[-1],pore_dt)
#
#fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(figsize=(12,8),nrows=2, ncols=2)
## First Signal
#ax1.plot(t1, pore_area[0], 'k', linewidth=1.5)
#ax1.set_xlabel('Time (minutes)')
#ax1.set_ylabel('Pixels')
#extent = [t1.min(),t1.max(),0,max(period)]
#im = NonUniformImage(ax3, interpolation='nearest', extent=extent)
#im.set_cmap('cubehelix_r')
#im.set_data(t1, period[:idx], power[:idx,:])
#ax3.images.append(im)
#ax3.set_ylabel('Period (minutes)')
#ax3.set_xlabel('Time (minutes)')
#ax3.contour(t1, period[:idx], sig95[:idx,:], [-99,1], colors='k', linewidths=2, extent=extent)
#ax3.fill(np.concatenate([t1, t1[-1:]+pore_dt, t1[-1:]+pore_dt,t1[:1]-pore_dt, t1[:1]-pore_dt]),
#        (np.concatenate([coi,[1e-9], period[-1:], period[-1:], [1e-9]])),
#        'k', alpha=0.3,hatch='/', zorder=100, antialiased=True, rasterized=True)
#ax3.set_xlim(t1.min(),t1.max())
#ax3.set_ylim(([min(period), period[idx]]))
## Second Signal
#ax2.plot(t1, pore_area[-1], 'k', linewidth=1.5)
#ax2.set_xlabel('Time (minutes)')
#ax2.set_ylabel('Pixels')
#extent = [t1.min(),t1.max(),0,max(period2)]
#im = NonUniformImage(ax4, interpolation='nearest', extent=extent)
#im.set_cmap('cubehelix_r')
#im.set_data(t1, period2[:idx], power2[:idx,:])
#ax4.images.append(im)
#ax4.set_ylabel('Period (minutes)')
#ax4.set_xlabel('Time (minutes)')
#ax4.contour(t1, period2[:idx], sig952[:idx,:], [-99,1], colors='k', linewidths=2, extent=extent)
#ax4.fill(np.concatenate([t1, t1[-1:]+pore_dt, t1[-1:]+pore_dt,t1[:1]-pore_dt, t1[:1]-pore_dt]),
#        (np.concatenate([coi,[1e-9], period[-1:], period[-1:], [1e-9]])),
#        'k', alpha=0.3,hatch='/', zorder=100, antialiased=True, rasterized=True)
#ax4.set_xlim(t1.min(),t1.max())
#ax4.set_ylim(([min(period), period[idx]]))
#
#fig.tight_layout()
#plt.savefig('/home/nabobalis/GitRepos/PhDThesis/Chapter2/Figs/pore_wavelet.pdf',dpi=300,bbox_inches='tight')

#for sig1, sig2 in zip(sun_area,sun_inten):
#    cross_wavelet(sig1, sig2, sun_dt, mother='morlet', plot=True)
#
#for sig1, sig2 in zip(pore_area,pore_inten):
#    cross_wavelet(sig1, sig2, pore_dt, mother='morlet', plot=True)
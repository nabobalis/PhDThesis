# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from pyhht import EMD
import numpy as np
import matplotlib.pyplot as plt

time = np.linspace(0,1000,1000)
data = np.sin(np.pi*2*time/100) + np.sin(np.pi*2*time/200)

exterma = EMD.exterma(data)
top, bot, mean = EMD.envelope(exterma[0],exterma[1],1000)
IMFs = EMD.emd(data)

# PLOT 1 SIGNAL
plt.figure(figsize=(20,8))
plt.plot(time,data,'b')
plt.gca().axes.get_xaxis().set_ticks([])
plt.gca().axes.get_yaxis().set_ticks([])
plt.savefig("C:\\Users\\Nabil\\Downloads\\signal.svg", dpi=300, bbox_inches='tight')

# PLOT 2 EXTERAM POINTS
plt.figure(figsize=(20,8))
plt.plot(time,data,'black')
plt.scatter(exterma[0][0],exterma[0][1],c='red', s=100)
plt.scatter(exterma[1][0],exterma[1][1],c='green', s=100)
plt.xlim(0,1000)
plt.gca().axes.get_xaxis().set_ticks([])
plt.gca().axes.get_yaxis().set_ticks([])
plt.savefig("C:\\Users\\Nabil\\Downloads\\exterma.svg", dpi=300, bbox_inches='tight')

# PLOT 3 ENVELOPE
plt.figure(figsize=(20,8))
plt.plot(time,data,'black')
plt.scatter(exterma[0][0],exterma[0][1],c='red', s=100)
plt.scatter(exterma[1][0],exterma[1][1],c='green', s=100)
plt.plot(time,bot,'green')
plt.plot(time,top,'red')
plt.xlim(0,1000)
plt.gca().axes.get_xaxis().set_ticks([])
plt.gca().axes.get_yaxis().set_ticks([])
plt.savefig("C:\\Users\\Nabil\\Downloads\\envelope.svg", dpi=300, bbox_inches='tight')

# PLOt 4 MEAN ENVELOPE
plt.figure(figsize=(20,8))
plt.plot(time,data,'black')
plt.scatter(exterma[0][0],exterma[0][1],c='red', s=100)
plt.scatter(exterma[1][0],exterma[1][1],c='green', s=100)
plt.plot(time,bot,'green')
plt.plot(time,top,'red')
plt.plot(time,mean,'blue',linestyle='--')
plt.xlim(0,1000)
plt.gca().axes.get_xaxis().set_ticks([])
plt.gca().axes.get_yaxis().set_ticks([])
plt.savefig("C:\\Users\\Nabil\\Downloads\\mean.svg", dpi=300, bbox_inches='tight')

# PLOT 5 IMF
plt.figure(figsize=(20,8))
plt.plot(time,IMFs[:,0],'black')
plt.gca().axes.get_xaxis().set_ticks([])
plt.gca().axes.get_yaxis().set_ticks([])
plt.savefig("C:\\Users\\Nabil\\Downloads\\IMFs.svg", dpi=300, bbox_inches='tight')

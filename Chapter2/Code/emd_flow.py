# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from pyhht import EMD
import numpy as np
import matplotlib.pyplot as plt

time = np.linspace(0,1000,1000)
a = np.sin(np.pi*2*time/100)
b = np.sin(np.pi*2*time/200)
data = a + b

exterma = EMD.exterma(data)
top, bot, mean = EMD.envelope(exterma[0],exterma[1],1000)
#IMFs = EMD.emd(data)

# PLOT 1 SIGNAL
plt.figure(figsize=(20,8))
plt.plot(time,data,'black')
plt.gca().axes.get_xaxis().set_ticks([])
plt.gca().axes.get_yaxis().set_ticks([])
plt.savefig("/home/nabobalis/signal.svg", dpi=300, bbox_inches='tight')
# PLOT 2 EXTERAM POINTS

plt.figure(figsize=(20,8))
plt.plot(time,data,'black')
plt.scatter(exterma[0][0],exterma[0][1],c='red', s=250)
plt.scatter(exterma[1][0],exterma[1][1],c='green', s=250)
plt.xlim(0,1000)
plt.gca().axes.get_xaxis().set_ticks([])
plt.gca().axes.get_yaxis().set_ticks([])
plt.savefig("/home/nabobalis/exterma.svg", dpi=300, bbox_inches='tight')

# PLOT 3 ENVELOPE
plt.figure(figsize=(20,8))
plt.plot(time,data,'black')
plt.scatter(exterma[0][0],exterma[0][1],c='red', s=250)
plt.scatter(exterma[1][0],exterma[1][1],c='green', s=250)
plt.plot(time,bot,'green', linewidth=2)
plt.plot(time,top,'red', linewidth=2)
plt.xlim(0,1000)
plt.gca().axes.get_xaxis().set_ticks([])
plt.gca().axes.get_yaxis().set_ticks([])
plt.savefig("/home/nabobalis/envelope.svg", dpi=300, bbox_inches='tight')

# PLOt 4 MEAN ENVELOPE
plt.figure(figsize=(20,8))
plt.plot(time,data,'black')
plt.scatter(exterma[0][0],exterma[0][1],c='red', s=250)
plt.scatter(exterma[1][0],exterma[1][1],c='green', s=250)
plt.plot(time,bot,'green', linewidth=2)
plt.plot(time,top,'red', linewidth=2)
plt.plot(time,mean,'blue',linestyle='--', linewidth=2)
plt.xlim(0,1000)
plt.gca().axes.get_xaxis().set_ticks([])
plt.gca().axes.get_yaxis().set_ticks([])
plt.savefig("/home/nabobalis/mean.svg", dpi=300, bbox_inches='tight')

# PLOT 5 IMF i.e the signal cuz
plt.figure(figsize=(20,8))
plt.plot(time,a,'black')
plt.gca().axes.get_xaxis().set_ticks([])
plt.gca().axes.get_yaxis().set_ticks([])
plt.savefig("/home/nabobalis/IMF.svg", dpi=300, bbox_inches='tight')

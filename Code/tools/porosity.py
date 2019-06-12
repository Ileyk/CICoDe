import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import matplotlib.cm as cm
from matplotlib.colors import rgb2hex
# from bokeh import mpl
# from bokeh.plotting import output_file, show, ColumnDataSource, figure, vplot
# from bokeh.models import HoverTool
import pandas as pd
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as ticker
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext
import holoviews as hv
# import colorcet as cc
from matplotlib.colors import LinearSegmentedColormap

from matplotlib.patches import Rectangle

import os

import subprocess

import imageio

from parfile import *

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

hv.notebook_extension("matplotlib")

def fmt(x, pos):
	a, b = '{:.1e}'.format(x).split('e')
	b = int(b)
	return r'{} $\cdot $10$^{{{}}}$'.format(a, b)

def fmtCbar(x, pos):
	a, b = '{:.1e}'.format(x).split('e')
	b = int(b)
	if (b!=0 and b!=1 and b!=-1):
		return r'{}E{}'.format(a, b)
	elif (b==0):
		return r'{}'.format(a)
	elif (b==1):
		return '{:2.0f}'.format(x)
	elif (b==-1):
		return '{:3.2f}'.format(x)

font = {'family' : 'normal',
'size'   : fontsize}

plt.rc('font', **font)


filename='../src/output/porosity'
proc = subprocess.Popen(['wc',filename], stdout=subprocess.PIPE)
tmp = proc.communicate()
Nr = int(tmp[0].split()[0])-1
f=open(filename)
tmp = f.readline() # 1st header line
r = np.zeros(Nr)
d = np.zeros(Nr)
Rcl = np.zeros(Nr)
dd = np.zeros(Nr)
fig, axes = plt.subplots(1,1, sharex=True,figsize=(1.5*figsize,figsize))
for j in range(Nr):
    tmp = f.readline()
    r[j]   = float(tmp.split()[0])
    d[j]   = float(tmp.split()[1])
    Rcl[j] = float(tmp.split()[2])
    dd[j]  = float(tmp.split()[3])

# ax1 = axes[0]
ax2 = axes.twinx()

axes.plot(r,d,color='r',linestyle='solid',linewidth=3)
axes.plot(r,Rcl,color='g',linestyle='solid',linewidth=3)
ax2.plot(r,dd,color='b',linestyle='solid',linewidth=3)

# axes[0].plot(r,d,color='r',linestyle='solid',linewidth=3)
# axes[1].plot(r,Rcl,color='g',linestyle='solid',linewidth=3)
# axes[2].plot(r,dd,color='b',linestyle='solid',linewidth=3)

# axes[0].set_ylabel(r'mean dist. between centers / stellar radius', fontweight='bold', fontsize=fontsize)
# axes[1].set_ylabel(r'radii / stellar radius', fontweight='bold', fontsize=fontsize)
# axes[2].set_ylabel(r'mean dist. between shells / local clump radius', fontweight='bold', fontsize=fontsize)

# for i, ax in enumerate(axes.flat):
axes.grid(which='major', linestyle='dotted', linewidth=2, alpha=0.5)
axes.grid(which='minor', linestyle='dotted', linewidth=2, alpha=0.5)
axes.set_xlabel(r'x / R${_*}$',fontweight='bold',fontsize=fontsize)
axes.set_ylabel(r'clump radius / R${_*}$',fontweight='bold',fontsize=fontsize,color='g')
plt.text(-0.18, 0.5, r'mean dist. between centers / R${_*}$',fontsize=fontsize,color='r',fontweight='bold',rotation=90,horizontalalignment='center',verticalalignment='center',transform = axes.transAxes)
ax2.set_ylabel(r'mean dist. between shells / R${_*}$',fontweight='bold',fontsize=fontsize,color='b')

ax2.tick_params(axis='y', labelcolor='b')

# axes.set_xlim(0.,np.max(r))
axes.set_xscale('log')

# fig1.set_ylim(-10.,10.)
# plt.show()
# stop

fig.tight_layout()
plt.savefig('porosity.png',format='png',dpi=140,bbox_inches='tight')

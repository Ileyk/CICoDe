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


filename='../src/output/initial_distributions'
proc = subprocess.Popen(['wc',filename], stdout=subprocess.PIPE)
tmp = proc.communicate()
Ncl = int(tmp[0].split()[0])-2 # 2 lines of header
r  = np.zeros(Ncl)
th = np.zeros(Ncl)
ph = np.zeros(Ncl)
Rcl = np.zeros(Ncl)
f=open(filename)
tmp = f.readline() # 2 lines of header
tmp = f.readline() # 2 lines of header
for i in range(Ncl):
	tmp = f.readline()
	r[i]  = float(tmp.split()[0])
	th[i] = float(tmp.split()[1])
	ph[i] = float(tmp.split()[2])
	Rcl[i] = float(tmp.split()[3])
Nbins=50 # Ncl/10000

fig, fig1 = plt.subplots(1,sharex=True,figsize=(1.1*figsize,figsize))

n, bins, patches = fig1.hist(r,Nbins,density=1)



# fig1.plot(r,mach_parker,color='k')
# fig1.plot(r,mach_beta[:,0],color='b',linestyle='solid')
# fig1.plot(r,mach_beta[:,1],color='r',linestyle='solid')
# fig1.plot(1,1,marker='o',color='k',markersize=10)
# plt.text(0.96, 0.24, 'Isothermal assumption',fontsize=2*fontsize/3,color='k',fontweight='bold',rotation=-90,horizontalalignment='center',verticalalignment='center',transform = fig1.transAxes)
# plt.text(0.5, 0.25, 'Isothermal Parker',fontsize=2*fontsize/3,color='k',fontweight='bold',horizontalalignment='center',verticalalignment='center',transform = fig1.transAxes)
# plt.text(0.5, 0.15, r'C-rich, $\beta=$'+str(beta[0])+r' & v$_{inf}$/c$_s=$'+str(eta[0])+' \n => r$_s\sim$'+'{0:.2g}'.format(rs[0])+r'R$_{dust}$',fontsize=2*fontsize/3,color='b',fontweight='bold',horizontalalignment='center',verticalalignment='center',transform = fig1.transAxes)
# plt.text(0.5, 0.05, r'O-rich, $\beta=$'+str(beta[1])+r' & v$_{inf}$/c$_s=$'+str(eta[1])+' \n => r$_s\sim$'+'{0:.2g}'.format(rs[1])+r'R$_{dust}$',fontsize=2*fontsize/3,color='r',fontweight='bold',horizontalalignment='center',verticalalignment='center',transform = fig1.transAxes)
fig1.set_xlabel(r'r / sonic radius',fontweight='bold',fontsize=fontsize)
fig1.set_ylabel(r'$Mach$', fontweight='bold', fontsize=fontsize)
fig1.grid(which='major', linestyle='dotted', linewidth=2, alpha=0.5)
fig1.grid(which='minor', linestyle='dotted', linewidth=2, alpha=0.5)
# fig1.set_xlim(rmin,rmax)
# fig1.set_ylim(0.,3.)
plt.show()
stop
# fig.tight_layout()
# fig.savefig(outputs+'difference_beta_law_Parker_wind.png',bbox_inches='tight')

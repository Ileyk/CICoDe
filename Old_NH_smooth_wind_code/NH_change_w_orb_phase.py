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

# display
# (M2,P)
# (3,4) | (3,8) | (3,12)
# (2,4) | (2,8) | (2,12)
# (1,4) | (1,8) | (1,12)
# => M2[i//3] and P[i%3]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

hv.notebook_extension("matplotlib")

# x = np.arange(0, np.pi, 0.1)
# y = np.arange(0, 2*np.pi, 0.1)
# X, Y = np.meshgrid(x, y)
# Z = np.cos(X) * np.sin(Y) * 10
# colors = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1)]  # R -> G -> B
# n_bins = [4, 6, 10, 100]  # Discretizes the interpolation into bins
# cmap_name = 'my_list'
# # fig, ax = plt.subplots(1, 1, figsize=(6, 6))
# vmax=25
# cm = LinearSegmentedColormap.from_list(cmap_name, [(0,    'blue'),
#                                               (0.2/vmax, 'red'),
#                                               (0.6/vmax, 'green'),
# 											  (2./vmax, 'black'),
# 											  (6./vmax, 'yellow'),
# 											  (20./vmax, 'orange'),
# 											  (vmax/vmax, 'cyan')], N=7)

# im = ax.imshow(Z, interpolation='nearest', origin='lower', cmap=cm)
# ax.set_title("N bins: %s" % 4)
# fig.colorbar(im, ax=ax)
# plt.show()

# fig.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.05)
# for ax in zip(axs.ravel()):
#     # Create the colormap
#     cm = LinearSegmentedColormap.from_list(
#         cmap_name, colors, N=4)
#     # Fewer bins will result in "coarser" colomap interpolation
#     im = ax.imshow(Z, interpolation='nearest', origin='lower', cmap=cm)
#     ax.set_title("N bins: %s" % 4)
#     fig.colorbar(im, ax=ax)
# plt.show()

# def plot_array():

def fmt(x, pos):
	a, b = '{:.1e}'.format(x).split('e')
	b = int(b)
	return r'{} $\cdot $10$^{{{}}}$'.format(a, b)

# def fmtCbar(x, pos):
# 	#a = '{:5.2f}'.format(x) #.split('e')[0]
# 	a = '{:5.1f}'.format(x)#.split('e')[0]
# 	return a
# 	# return r'{}'.format(a)

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

#'weight' : 'bold',

plt.rc('font', **font)

fig, fig1 = plt.subplots(1,figsize=(1.1*figsize,figsize))

fig1.text(0.7, 0.9, r'For a 27 degrees inclination', fontsize=fontsize, color='k', fontweight='bold', horizontalalignment='center',verticalalignment='center',transform = fig1.transAxes)
fig1.text(0.7, 0.8, r'q=4,beta=2,f=95%', fontsize=fontsize, color='b', fontweight='bold', horizontalalignment='center',verticalalignment='center',transform = fig1.transAxes)
fig1.text(0.7, 0.7, r'q=1,beta=0.5,f=95%', fontsize=fontsize, color='r', fontweight='bold', horizontalalignment='center',verticalalignment='center',transform = fig1.transAxes)

# DOF for Cygnus X-1
al_max=27.*(math.pi/180.) # inclination of LOS w/ respect to orb. ang. mom.
f=0.95 # minor influence within realistic range of variation
# Also, q and beta defined below

maxNH=1.E99

# Numerical D.O.F.
Nstep=1000 # to integrate NH. half of steps before min distance, linear ; half after, log
Nph=32 # # of phases sampled over half orb. period

for p in range(2):
	if (p==0):
		q=1. # increases the peak-to-peak (from 1 to 4)
		beta=0.5 # increases the peak-to-peak (as expected) (from 0.5 to 2)
		color_curve='r'
	else:
		q=4. # increases the peak-to-peak (from 1 to 4)
		beta=2. # increases the peak-to-peak (as expected) (from 0.5 to 2)
		color_curve='b'

	Egg=(0.6*q**(2./3.))/(0.49*q**(2./3.)+np.log(1+q**(1./3.)))

	print 'Max inclination before eclipse : ', (math.pi/2.-np.arctan(f*Egg))*(180./math.pi)

	ph_all=np.zeros(Nph)
	NH_all=np.zeros(Nph) # to store results for plotting
	for j in range(Nph):
		ph=(j+0.5)*(0.5/Nph) # ph goes from 0 to 0.5 and al taken negative close from inferior conjuntion
		# al=al_max-(2.*al_max/0.5)*ph # instantaneous inclination corresponding to phase ph
		al=np.arcsin(np.sin(al_max)*np.cos(2*math.pi*ph)) # instantaneous inclination corresponding to phase ph
		# integration along LOS, from accreting BH
		NH=0.
		dl_lin=np.abs(np.sin(al))/(Nstep/2) # l step for 1st half, linear
		l=dl_lin/2. # centered integration (@ least on the linear 1st half)
		q=(10./np.abs(np.sin(al)))**(2./Nstep) # reason for log serie from l=sin(al) to l=10
		for i in range(Nstep):
			if (l<np.abs(np.sin(al))):
				l=l+dl_lin
				r2=1.-2.*l*np.sin(al)+l**2.
				r =np.sqrt(r2)
				if (r>f*Egg):
					NH=NH+dl_lin/(r2*(1.-f*Egg/r)**beta)
				else:
					NH=maxNH
			else:
				dl_log=l*(q-1.)
				l=l+dl_log
				r2=1.-2.*l*np.sin(al)+l**2.
				r =np.sqrt(r2)
				if (r>f*Egg):
					NH=NH+dl_log/(r2*(1.-f*Egg/r)**beta)
				else:
					NH=maxNH
		ph_all[j]=ph
		NH_all[j]=NH

	bb=np.where(NH_all<maxNH/2.)
	# print np.max(NH_all[bb])
	# print np.max(NH_all), np.max(NH_all(np.where(NH_all<maxNH/2.)))
	# NH_all[np.where(NH_all==maxNH)]=np.max(NH_all[bb])
	# print np.max(NH_all)
	# stop
	# cooling=[r'$\beta=$2',r'$\beta=$3']
	# leg1=2*['']
	# leg1[0] = plt.Line2D((0,0),(0,1),linestyle='solid',linewidth=3.0,color='k')
	# leg1[1] = plt.Line2D((0,0),(0,1),linestyle='dashed',linewidth=3.0,color='k')

	fig1.plot(ph_all,NH_all,linestyle='solid',marker='o',linewidth=3,color=color_curve)

fig1.set_xlim(0.,0.5)
fig1.set_ylim(0.,1.1*np.max(NH_all[bb]))

# extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
# legend_handle = [extra, leg1[0], extra, leg1[1]]
# leg_lab=np.concatenate([[cooling[0],''],[cooling[1],'']])
# fig1.legend(legend_handle,leg_lab,
# 	loc='upper right',fontsize=16,ncol=2,shadow=True,handletextpad=-2,numpoints=1)
fig1.set_xlabel(r'Orbital phase',fontweight='bold',fontsize=fontsize)
fig1.set_ylabel(r'NH (normalized)', fontweight='bold', fontsize=fontsize)
fig1.grid(which='major', linestyle='dotted', linewidth=2, alpha=0.9)
fig1.grid(which='minor', linestyle='dotted', linewidth=2, alpha=0.9)
fig.savefig(outputs+'NH_change_w_orb_phase.png',bbox_inches='tight')

# plt.show()
# stop

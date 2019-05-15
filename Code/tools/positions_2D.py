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

def circles(x, y, s, c='b', vmin=None, vmax=None, **kwargs):
    """
    Make a scatter of circles plot of x vs y, where x and y are sequence
    like objects of the same lengths. The size of circles are in data scale.

    Parameters
    ----------
    x,y : scalar or array_like, shape (n, )
        Input data
    s : scalar or array_like, shape (n, )
        Radius of circle in data unit.
    c : color or sequence of color, optional, default : 'b'
        `c` can be a single color format string, or a sequence of color
        specifications of length `N`, or a sequence of `N` numbers to be
        mapped to colors using the `cmap` and `norm` specified via kwargs.
        Note that `c` should not be a single numeric RGB or RGBA sequence
        because that is indistinguishable from an array of values
        to be colormapped. (If you insist, use `color` instead.)
        `c` can be a 2-D array in which the rows are RGB or RGBA, however.
    vmin, vmax : scalar, optional, default: None
        `vmin` and `vmax` are used in conjunction with `norm` to normalize
        luminance data.  If either are `None`, the min and max of the
        color array is used.
    kwargs : `~matplotlib.collections.Collection` properties
        Eg. alpha, edgecolor(ec), facecolor(fc), linewidth(lw), linestyle(ls),
        norm, cmap, transform, etc.

    Returns
    -------
    paths : `~matplotlib.collections.PathCollection`

    Examples
    --------
    a = np.arange(11)
    circles(a, a, a*0.2, c=a, alpha=0.5, edgecolor='none')
    plt.colorbar()

    License
    --------
    This code is under [The BSD 3-Clause License]
    (http://opensource.org/licenses/BSD-3-Clause)
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    from matplotlib.collections import PatchCollection

    if np.isscalar(c):
        kwargs.setdefault('color', c)
        c = None
    if 'fc' in kwargs: kwargs.setdefault('facecolor', kwargs.pop('fc'))
    if 'ec' in kwargs: kwargs.setdefault('edgecolor', kwargs.pop('ec'))
    if 'ls' in kwargs: kwargs.setdefault('linestyle', kwargs.pop('ls'))
    if 'lw' in kwargs: kwargs.setdefault('linewidth', kwargs.pop('lw'))

    patches = [Circle((x_, y_), s_) for x_, y_, s_ in np.broadcast(x, y, s)]
    collection = PatchCollection(patches, **kwargs)
    if c is not None:
        collection.set_array(np.asarray(c))
        collection.set_clim(vmin, vmax)

    ax = plt.gca()
    ax.add_collection(collection)
    ax.autoscale_view()
    if c is not None:
        plt.sci(collection)
    return collection

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


filename='../src/output/positions'
proc = subprocess.Popen(['wc',filename], stdout=subprocess.PIPE)
tmp = proc.communicate()
Nlines = int(tmp[0].split()[0])-2
f=open(filename)
tmp = f.readline() # 1st header line
tmp = f.readline() # 2nd header line
Nphases = int(tmp.split()[2])
# Ncl     = int(tmp.split()[3])
# x = np.zeros(Ncl)
# y = np.zeros(Ncl)
# Rcl = np.zeros(Ncl)
# for j in range(Nphases):
# 	print j
# 	for i in range(Ncl):
# where the snapshot are stored before being bound together into a GIF
images = []
for j in range(Nphases):
    print j, '/', Nphases
    fig, fig1 = plt.subplots(1,sharex=True,figsize=(1.1*figsize,figsize))
    fig1.set_xlim(-10.,10.)
    fig1.set_ylim(-10.,10.)
    circle3 = plt.Circle((0.,0.), 1., color='b', clip_on=False)
    fig1.add_artist(circle3)
    for i in range(Nlines):
        tmp = f.readline()
        if (tmp.split()[0]=='xxx'):
            break
    	# x[i]  = float(tmp.split()[0])
    	# y[i]  = float(tmp.split()[1])
    	# Rcl[i]  = float(tmp.split()[2])
        xx   =  float(tmp.split()[0])
        yy   =  float(tmp.split()[1])
        Rcll =  float(tmp.split()[2])

        circles(xx,yy,Rcll,c='k',alpha=0.5, edgecolor='none')
	# circles(x,y,Rcl,c='k',alpha=0.5, edgecolor='none')

    fig1.set_xlabel(r'x / stellar radius',fontweight='bold',fontsize=fontsize)
    fig1.set_ylabel(r'y / stellar radius', fontweight='bold', fontsize=fontsize)
    fig.tight_layout()
    plt.savefig('pif.png',format='png',dpi=70,bbox_inches='tight') # overwrite same file each time
    images.append(imageio.imread('pif.png')) # but before, add this file to images

imageio.mimsave('pif.gif',images,duration=15./float(Nphases)) # duration per snapshot such as the whole GIF lasts 10s

	# plt.colorbar()
	# Rcl_in_points = np.diff(fig1.transData.transform(zip([0]*len(Rcl), Rcl)))
	# fig1.scatter(x,y,c='k',s=Rcl_in_points**2/30.,marker='o',edgecolors='none')

	# fig1.plot(x,y,color='k',linewidth=0,markersize=1,marker='o')


# fig1.plot(r,mach_parker,color='k')
# fig1.plot(r,mach_beta[:,0],color='b',linestyle='solid')
# fig1.plot(r,mach_beta[:,1],color='r',linestyle='solid')
# fig1.plot(1,1,marker='o',color='k',markersize=10)
# plt.text(0.96, 0.24, 'Isothermal assumption',fontsize=2*fontsize/3,color='k',fontweight='bold',rotation=-90,horizontalalignment='center',verticalalignment='center',transform = fig1.transAxes)
# plt.text(0.5, 0.25, 'Isothermal Parker',fontsize=2*fontsize/3,color='k',fontweight='bold',horizontalalignment='center',verticalalignment='center',transform = fig1.transAxes)
# plt.text(0.5, 0.15, r'C-rich, $\beta=$'+str(beta[0])+r' & v$_{inf}$/c$_s=$'+str(eta[0])+' \n => r$_s\sim$'+'{0:.2g}'.format(rs[0])+r'R$_{dust}$',fontsize=2*fontsize/3,color='b',fontweight='bold',horizontalalignment='center',verticalalignment='center',transform = fig1.transAxes)
# plt.text(0.5, 0.05, r'O-rich, $\beta=$'+str(beta[1])+r' & v$_{inf}$/c$_s=$'+str(eta[1])+' \n => r$_s\sim$'+'{0:.2g}'.format(rs[1])+r'R$_{dust}$',fontsize=2*fontsize/3,color='r',fontweight='bold',horizontalalignment='center',verticalalignment='center',transform = fig1.transAxes)
# fig1.grid(which='major', linestyle='dotted', linewidth=2, alpha=0.5)
# fig1.grid(which='minor', linestyle='dotted', linewidth=2, alpha=0.5)
# fig1.set_xlim(rmin,rmax)
# fig1.set_ylim(0.,3.)
# plt.show()
# stop
# fig.tight_layout()
# fig.savefig(outputs+'difference_beta_law_Parker_wind.png',bbox_inches='tight')

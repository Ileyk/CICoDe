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

import glob

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

NH_shells_fl='../src/output/NH_shells'
proc = subprocess.Popen(['wc',NH_shells_fl], stdout=subprocess.PIPE)
tmp = proc.communicate()
Nlines = int(tmp[0].split()[0])-1 # 1 line of header

rR=np.zeros(Nlines)
NH_bins=np.zeros(Nlines)
NH_cmltd=np.zeros(Nlines)

f=open(NH_shells_fl)
tmp = f.readline() # 1st header line
for j in range(Nlines):
    tmp = f.readline() # 1st header line
    rR[j]       = float(tmp.split()[2])
    NH_bins[j]  = float(tmp.split()[3])
    NH_cmltd[j] = float(tmp.split()[4])

fig, axes = plt.subplots(1,1, sharex=True,figsize=(1.5*figsize,figsize))
ax2 = axes.twinx()
axes.step(rR,NH_bins,color='g',linestyle='solid',linewidth=2,where='post')
# Normalize cumulated NH to its final value (%) ...
ax2.step(rR,100.*NH_cmltd/NH_cmltd[np.size(NH_cmltd)-1],color='b',linestyle='solid',linewidth=2,where='post')
# ... or not
# ax2.step(rR,NH_cmltd,color='b',linestyle='solid',linewidth=2,where='post')
axes.grid(which='major', linestyle='dotted', linewidth=2, alpha=0.5)
axes.grid(which='minor', linestyle='dotted', linewidth=2, alpha=0.5)
axes.set_xlabel(r'x / R${_*}$',fontweight='bold',fontsize=fontsize)
axes.set_ylabel(r'NH per radial bin',fontweight='bold',fontsize=fontsize,color='g')
axes.tick_params(axis='y',labelcolor='g')
ax2.set_ylabel(r'Cumulated NH (%)',fontweight='bold',fontsize=fontsize,color='b')
# ax2.set_ylabel(r'Cumulated NH',fontweight='bold',fontsize=fontsize,color='b')
ax2.tick_params(axis='y',labelcolor='b')
# axes.set_xscale('log')
fig.tight_layout()
plt.savefig('NH_shells.png',format='png',dpi=140,bbox_inches='tight')

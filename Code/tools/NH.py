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

# 1. Read all the NH bins available up to positions_Nph.dat
# (w/ Nph * dNsave < Nphases)

Nph = 100
NH_fl='../src/output/NH'
log_fl='../src/output/log'
proc = subprocess.Popen(['wc',NH_fl], stdout=subprocess.PIPE)
tmp = proc.communicate()

# 1. Read Nphases from log file (2600)
f=open(log_fl)
tmp = f.read()
tmp_split = tmp.split('\n')
index = [i for i, s in enumerate(tmp_split) if 'The # of phase bins required is' in s]
Nphases = int(tmp_split[index[0]].split()[len(tmp_split[index[0]].split())-1])
f.close()

# 2. Read Nphases from 1st line of NH file (should be the same) (2600)
f=open(NH_fl)
tmp = f.readline() # 1st header line
if (Nphases!=int(tmp.split()[len(tmp.split())-1])):
    print Nphases, int(tmp.split()[len(tmp.split())-1])
    print "Unconsistent # of phase bins between log and NH files"
    stop
f.close()
# Nphases = int(tmp.split()[len(tmp.split())-1])

# 3. Read Nsave from log file (100)
f=open(log_fl)
tmp = f.read()
tmp_split = tmp.split('\n')
index = [i for i, s in enumerate(tmp_split) if 'Among them, we save only' in s]
Nsave = int(tmp_split[index[0]].split()[len(tmp_split[index[0]].split())-1])
f.close()

# 4. Check how many .dat files have already been produced.
# Which is the index of the last saved file?
# Look for how many data files produced by checking # of files
# starting w/ positions_ and ending w/ .dat (includes 0000)
list = glob.glob("../src/output/positions_*.dat")
Nfiles=np.shape(list)[0]-1
if (Nph>Nfiles):
    print 'The data files you want to plot have not been produced yet, lower Nph'
    stop

dNsave=Nphases/Nsave # 26
if (Nph > Nsave):
    print 'The data files you want to plot will not been produced, Nph larger than Nsave!'
    stop

# Read NH chunks, dNsave by dNsave up to dNsave*Nph
NH = np.zeros((dNsave*Nph,Nphases))
f=open(NH_fl)
for i in range(Nph):
    print i+1, '/', Nph
    for j in range(dNsave):
        tmp = f.readline() # 1st header line
        for k in range(Nphases):
            tmp = f.readline() # 1st header line
            NH[i*dNsave+j][k] = float(tmp.split()[0])
        tmp = f.readline() # blank line
f.close()

dph = 1./Nphases
# orbital phase edge
phases = np.linspace(-dph/2.,1.-dph/2.,num=Nphases+1,endpoint=True)
# orbital phase centered
phasesCC = np.linspace(0.,1.-dph,num=Nphases,endpoint=True)

images = []

for j in range(3,Nph):

    print j+1, '/', Nph

    NH_median = np.zeros(Nphases)
    NH_10 = np.zeros(Nphases)
    NH_30 = np.zeros(Nphases)
    NH_70 = np.zeros(Nphases)
    NH_90 = np.zeros(Nphases)

    for i in range(Nphases):
        NH_median[i]=np.median(NH[0:dNsave*(j+1),[i]])
        NH_10[i]=np.percentile(NH[0:dNsave*(j+1),[i]],10.)
        NH_30[i]=np.percentile(NH[0:dNsave*(j+1),[i]],30.)
        NH_70[i]=np.percentile(NH[0:dNsave*(j+1),[i]],70.)
        NH_90[i]=np.percentile(NH[0:dNsave*(j+1),[i]],90.)

    fig, fig1 = plt.subplots(1,sharex=True,figsize=(2*figsize,figsize))
    # fig1.set_xlim(np.min(phases),np.max(phases))
    fig1.set_xlim(0,2.)
    # fig1.set_ylim(0.5*np.min(NH_median[np.where(NH_median>0.)]),2.*np.max(NH_median[np.where(NH_median<1E10)]))
    fig1.set_ylim(0.0001,0.1)
    fig1.set_xlabel(r'Orbital phase',fontweight='bold',fontsize=fontsize)
    fig1.set_ylabel(r'NH', fontweight='bold', fontsize=fontsize)
    # Plot instantaneous NH for X-ray source motion w/ ph=0 @ time = 0
    fig1.plot(j*dNsave*dph,NH[dNsave*j,[dNsave*j]],marker='P',color='r',markersize=10)

    fig1.step(phasesCC,NH_median,where='mid',label='mid',color='k')
    fig1.step(phasesCC,NH_10,where='mid',label='mid',color='grey',alpha=0.3)
    fig1.step(phasesCC,NH_30,where='mid',label='mid',color='grey')
    fig1.step(phasesCC,NH_70,where='mid',label='mid',color='grey')
    fig1.step(phasesCC,NH_90,where='mid',label='mid',color='grey',alpha=0.3)
    # Repeat plot from 1 to 2, for better visualization
    fig1.step(phasesCC+1.,NH_median,where='mid',label='mid',color='k')
    fig1.step(phasesCC+1.,NH_10,where='mid',label='mid',color='grey',alpha=0.3)
    fig1.step(phasesCC+1.,NH_30,where='mid',label='mid',color='grey')
    fig1.step(phasesCC+1.,NH_70,where='mid',label='mid',color='grey')
    fig1.step(phasesCC+1.,NH_90,where='mid',label='mid',color='grey',alpha=0.3)
    # fig1.step(phasesCC+1.,NH_median,where='mid',label='mid',color='b')
    fig1.set_yscale('log')
    fig1.grid(which='major', linestyle='dotted', linewidth=2, alpha=0.5)
    fig1.grid(which='minor', linestyle='dotted', linewidth=2, alpha=0.5)

    # plt.show()
    # fig.tight_layout()
    # plt.savefig('NH.png',format='png',dpi=140,bbox_inches='tight') # overwrite same file each time

    plt.savefig('NH.png',format='png',dpi=90,bbox_inches='tight') # overwrite same file each time
    images.append(imageio.imread('NH.png')) # but before, add this file to images

    plt.close()

imageio.mimsave('NH.gif',images,duration=15./float(Nph)) # duration per snapshot such as the whole GIF lasts 10s

stop

# hist, bins = np.histogram(NH[0],bins)
# print NH[0]
# print hist
# print bins
# print 'ok'
# stop


Nlines = int(tmp[0].split()[0]) /(+2) # +2 for 1 header line and 1 blank line


# 1. Plot all the orbit - - - - - -

orb_fl='../src/output/orbit'
proc = subprocess.Popen(['wc',orb_fl], stdout=subprocess.PIPE)
tmp = proc.communicate()
Nphases = int(tmp[0].split()[0])-1
phases = np.zeros(Nphases)
pos_Xx = np.zeros(Nphases)
pos_Xy = np.zeros(Nphases)
pos_Xz = np.zeros(Nphases)
f=open(orb_fl)
tmp = f.readline() # 1st header line
for i in range(Nphases):
    tmp = f.readline()
    phases[i]   =  float(tmp.split()[0])
    pos_Xx[i]   =  float(tmp.split()[1])
    pos_Xy[i]   =  float(tmp.split()[2])
    pos_Xz[i]   =  float(tmp.split()[3])



# 2. For each snapshot, plot... - - - - - -

images = []

# Max radius up to which we plot
rmax=2.

# Max number of clumps to be plotted
Ncl_max=10000

# 2.1 ... the instantaneous position of the X-ray source

for i in range(Nfiles):

    print i, '/', Nfiles

    fig, fig1 = plt.subplots(1,sharex=True,figsize=(1.1*figsize,figsize))
    fig1.set_xlim(-rmax,rmax)
    fig1.set_ylim(-rmax,rmax)
    fig1.set_xlabel(r'x / stellar radius',fontweight='bold',fontsize=fontsize)
    fig1.set_ylabel(r'y / stellar radius', fontweight='bold', fontsize=fontsize)

    # Plot the orbit - - -
    # inferior part of the orbit
    fig1.plot(pos_Xx[np.where(pos_Xz>0.)],pos_Xy[np.where(pos_Xz>0.)],color='g',linestyle='solid')
    # superior part of the orbit (eclipse not plotted)...
    # ... on the right of the star
    fig1.plot(pos_Xx[np.where((pos_Xz<0.) & (np.sqrt(pos_Xx**2.+pos_Xy**2.)>1.) & (pos_Xx>0.)) ],pos_Xy[np.where((pos_Xz<0.) & (np.sqrt(pos_Xx**2.+pos_Xy**2.)>1.) & (pos_Xx>0.))],color='g',linestyle='dashed')
    # ... on the left of the star
    fig1.plot(pos_Xx[np.where((pos_Xz<0.) & (np.sqrt(pos_Xx**2.+pos_Xy**2.)>1.) & (pos_Xx<0.)) ],pos_Xy[np.where((pos_Xz<0.) & (np.sqrt(pos_Xx**2.+pos_Xy**2.)>1.) & (pos_Xx<0.))],color='g',linestyle='dashed')

    pos_fl='../src/output/positions_'+str(i).zfill(4)+'.dat'
    proc = subprocess.Popen(['wc',pos_fl], stdout=subprocess.PIPE)
    tmp = proc.communicate()
    Ncl = int(tmp[0].split()[0])-4
    f=open(pos_fl)
    tmp = f.readline() # 1st header line
    tmp = f.readline() # 2nd header line
    tmp = f.readline()
    pos_Xx_inst=float(tmp.split()[0])
    pos_Xy_inst=float(tmp.split()[1])
    pos_Xz_inst=float(tmp.split()[2])
    if (pos_Xz_inst>0.):
        circl_CO = plt.Circle((pos_Xx_inst,pos_Xy_inst), 0.2, color='g', clip_on=False, alpha=0.8)
    if (pos_Xz_inst<0.):
        circl_CO = plt.Circle((pos_Xx_inst,pos_Xy_inst), 0.2, color='g', fill=False, clip_on=False, alpha=0.8)

    circl_str = plt.Circle((0.,0.), 1., color='b', clip_on=False, alpha=0.5)
    fig1.add_artist(circl_str)
    # fig1.add_artist(circl_CO)
    fig1.plot(pos_Xx_inst,pos_Xy_inst,color='k',linewidth=0,markersize=10,marker='P')

    tmp = f.readline()
    # Read the 10000 first clump positions
    index=0
    while (index<min(Ncl_max,Ncl)):
        tmp = f.readline()
        # If EOF, quit
        if tmp == '':
            break
        xx   =  float(tmp.split()[0])
        yy   =  float(tmp.split()[1])
        # if outside plotting window, skip
        if (np.abs(xx)>rmax or np.abs(yy)>rmax):
            continue
        Rcll =  float(tmp.split()[3])
        circles(xx,yy,Rcll,c='k',alpha=0.5, edgecolor='none')
        index = index + 1

    fig1.set_xlabel(r'x / stellar radius',fontweight='bold',fontsize=fontsize)
    fig1.set_ylabel(r'y / stellar radius', fontweight='bold', fontsize=fontsize)
    plt.text(0.5, 1.04, str(index)+' in '+str(Ncl)+' clumps represented',fontsize=fontsize,color='k',fontweight='bold',horizontalalignment='center',verticalalignment='center',transform = fig1.transAxes)
    plt.savefig('pif.png',format='png',dpi=140,bbox_inches='tight') # overwrite same file each time
    images.append(imageio.imread('pif.png')) # but before, add this file to images

    plt.close()

imageio.mimsave('pif.gif',images,duration=15./float(Nfiles)) # duration per snapshot such as the whole GIF lasts 10s

plt.show()
stop

# pos_fl='../src/output/positions'
# proc = subprocess.Popen(['wc',pos_fl], stdout=subprocess.PIPE)
# tmp = proc.communicate()
# Nlines = int(tmp[0].split()[0])-2
# f=open(pos_fl)
# tmp = f.readline() # 1st header line
# tmp = f.readline() # 2nd header line
# Nphases = int(tmp.split()[1])
# # Ncl     = int(tmp.split()[3])
# # x = np.zeros(Ncl)
# # y = np.zeros(Ncl)
# # Rcl = np.zeros(Ncl)
# # for j in range(Nphases):
# # 	print j
# # 	for i in range(Ncl):
# # where the snapshot are stored before being bound together into a GIF
# images = []
# for j in range(Nphases):
#     print j+1, '/', Nphases
#     fig, fig1 = plt.subplots(1,sharex=True,figsize=(1.1*figsize,figsize))
#     fig1.set_xlim(-10.,10.)
#     fig1.set_ylim(-10.,10.)
#     for i in range(Nlines):
#         tmp = f.readline()
#         if (tmp.split()[0]=='xxx'):
#             break
#     	# x[i]  = float(tmp.split()[0])
#     	# y[i]  = float(tmp.split()[1])
#     	# Rcl[i]  = float(tmp.split()[2])
#         xx   =  float(tmp.split()[0])
#         yy   =  float(tmp.split()[1])
#         Rcll =  float(tmp.split()[3])
#
#         circles(xx,yy,Rcll,c='k',alpha=0.5, edgecolor='none')
# 	# circles(x,y,Rcl,c='k',alpha=0.5, edgecolor='none')
#     circle3 = plt.Circle((0.,0.), 1., color='b', clip_on=False, alpha=0.5)
#     fig1.add_artist(circle3)
#     fig.tight_layout()
#     plt.savefig('pif.png',format='png',dpi=140,bbox_inches='tight') # overwrite same file each time
#     # plt.show()
#     images.append(imageio.imread('pif.png')) # but before, add this file to images
#
# imageio.mimsave('pif.gif',images,duration=15./float(Nphases)) # duration per snapshot such as the whole GIF lasts 10s

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

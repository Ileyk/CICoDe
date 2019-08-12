import math
import numpy as np
import matplotlib as mpl

outputs='output/'

figsize=7
fontsize=20
colors=['#8B0000','#006400','#00008B','#F08080','#00FF00','#87CEFA']
linestyles=['solid','dashed']
mark=["s","o","^"]
mpl.rc('font',size=fontsize,weight='bold')

typeMean='arithmetic' # ** user-specified ** type of mean (arithmetic, geometric or harmonic)

secinyear=3600.*24.*365.25
msolar=2.E33
mH=1.67E-24 # mass of a proton

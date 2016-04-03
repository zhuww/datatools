from triplot import getsigmalevels, makesubplot2d
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import os, sys
import libstempo 

#data = np.load('1950.before.npz')
data = np.load('1950.after.npz')

fitpars = list(data['fitpars'])
res = data['res']

#psr = libstempo.tempopulsar(parfile='before_DDGR.par', timfile='J1950+2414_T2.tim') 
psr = libstempo.tempopulsar(parfile='after_DDGR.par', timfile='1950.all.tim')
psr.fit()
pars = psr.pars() 
vals = psr.vals()
errs = psr.errs()
fitidx = [i for i,p in enumerate(pars) if p in fitpars ]
#print fitidx
vals0 = vals[fitidx]
errs0 = errs[fitidx]
parsize = vals0.size
#print parsize
#print 'fit for: ', fitpars

plist = fitpars
MarkovChain = res[:,1:]

MCMCSize = len(MarkovChain)

def extract(par):
    ipar = plist.index(par)
    return MarkovChain[:,ipar]*errs0[ipar] + vals0[ipar]

f, ax = plt.subplots()

M2 = extract('M2').astype(float)
MTOT = extract('MTOT').astype(float)
M1 = MTOT - M2

makesubplot2d(ax, M2, M1) 
ax.set_xlim(0., 0.5)
ax.set_ylim(0., 3.0)
show()

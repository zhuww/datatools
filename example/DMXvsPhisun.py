from datatools.tempo import *
from pylab import *
from numpy import genfromtxt
import numpy as np
from scipy.stats import linregress,std
import sys,os

m = model('J1909.zww.par')
t = TOAfile('J1909.zww.tim')

os.system('tempo -f %s %s -a' % (m.parfile, t.toafile)) #run tempo to generate the chisun.tmp file

m.tempofit(t)
m.average()
phisun = genfromtxt('phisun.tmp')
#print np.unique(phisun).size, len(m.dmxlist[0])
print len(phisun), len(m.toa) #should match
#errorbar(phisun, m.res, yerr=m.err)
#print phisun.size, m.res.size,m.err.size
#print phisun, m.res
#plot(phisun, m.res)
#m.plot('date', 'res', LegendOn=True)
#m.plot('date', 'averes', LegendOn=True)
#m.plot('freq', 'res')
nm = model(m.newpar.parfile)
DMX, DMXErr, DMXR1, DMXR2 = nm.dmxlist
DMXgrp = {}
for j in DMXR1.keys():
    DMXgrp[j] = []
for i in range(len(m.toa)):
    for j in DMXR1.keys():
        if float(m.toa[i]) > float(DMXR1[j]) and float(m.toa[i]) < float(DMXR2[j]):
            DMXgrp[j].append(i)
dmxphi = [np.mean([phisun[i] for i in DMXgrp[j]]) for j in DMXgrp.keys()]
#print dmxphi
dmx = np.array([float(DMX[j]) for j in DMX.keys()])
dmxR = np.array([float(DMXR1[j] + DMXR2[j])/2. for j in DMXR1.keys()])
a,b,r,p,s = linregress(dmxR, dmx)
newdmx = dmx - a*dmxR - b
dmx = newdmx#-mean(newdmx)
dmxerr = [float(DMXErr[j]) for j in DMXErr.keys()]
#dmxbar = np.array([float(DMXR2[j] - DMXR1[j])/2. for j in DMXR1.keys()])
#nm.plot('date', 'DMX')
errorbar(dmxR, dmx, yerr= dmxerr, fmt='.')
show()
errorbar(dmxphi, dmx, yerr=dmxerr, fmt='.')
show()
#print m.dof, m.chisq
#print m.wrms, m.avewrms



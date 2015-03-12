#!/usr/bin/python
from datatools.tempo import *
import sys,os

import argparse
if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='comparetempo.py')
    parser.add_argument('par1', help='parfile 1')
    parser.add_argument('par2', help='parfile 2')
    parser.add_argument('--withdmx', help='withdmx', type=bool, default=False)
    parser.add_argument('--plot', help='makeplot', type=bool, default=False)
    parser.add_argument('--LegendLoc', help='LegendLoc', type=int, default=2)
    values = parser.parse_args(sys.argv[1:])

    f1, f2 = values.par1, values.par2
    withdmx = values.withdmx
    makeplot = values.plot
    LegLoc = values.LegendLoc
    #print withdmx, type(withdmx)

    m1 = model(f1)
    m2 = model(f2)

    par1 = [p for p in m1.manifest if not p.startswith('DMX_') and p in m1.parameters]
    dmx1 = [p for p in m1.manifest if p.startswith('DMX_') and p in m1.parameters]
    par2 = [p for p in m2.manifest if not p.startswith('DMX_') and p in m2.parameters]
    dmx2 = [p for p in m2.manifest if p.startswith('DMX_') and p in m2.parameters]

    par = [p for p in par1 if p in par2]
    dmx = [p for p in dmx1 if p in dmx2]
    if withdmx:
        par = par + dmx

    print 'parname:| par1 \t|\t par2 \t|\t abs(par1-par2)/max(sigma1, sigma2)'
    value1 = []
    value2 = []
    error1 = []
    error2 = []
    deviat = []
    redlist = []
    i = 0
    for p in par:
        if p in ['RAJ', 'DECJ']:
            val1 = float(m1.__dict__[p][0].split(':')[-1])
            val2 = float(m2.__dict__[p][0].split(':')[-1])
            err1 = float(m1.__dict__[p][1])
            err2 = float(m2.__dict__[p][1])
        else:
            try:
                val1 = float(m1.__dict__[p][0])
                val2 = float(m2.__dict__[p][0])
                err1 = float(m1.__dict__[p][1])
                err2 = float(m2.__dict__[p][1])
            except:
                pass
        try:
            dev = abs(val1-val2)/max((err1, err2))
        except:
            dev = np.nan
        value1.append(val1)
        value2.append(val2)
        error1.append(err1)
        error2.append(err2)
        deviat.append(dev)
        if dev > 3:
            print p, ':', m1.__dict__[p][0], '\t|\t',  m2.__dict__[p][0], '\t|\t', dev
            redlist.append(i)
        i+=1

    if makeplot:
        #from pylab import *
        import matplotlib.pyplot as plt
        fig, ax1 = plt.subplots()#figsize=(10,6))
        fig.canvas.set_window_title('compare %s %s' % ( f1, f2))
        ax1.yaxis.grid(True, linestyle='-', which='major', color='lightblue')
        ax1.xaxis.grid(True, linestyle='-', which='major', color='lightblue')
        xticks = plt.setp(ax1, xticks=range(len(par)))
        xtickNames = plt.setp(ax1, xticklabels=par)
        #xtickNames = plt.setp(ax1, xticklabels=[('%s' % s) for s in range(len(par))])
        plt.setp(xtickNames, rotation=45, fontsize=8)
        xvals = np.array(range(len(par)))
        numpars = len(par)
        ax1.set_xlim(0.5, numpars+0.5)
        l1 = ax1.errorbar(xvals, (np.array(value2) - np.array(value1))/np.array(error1), yerr=np.array(error2)/np.array(error1), fmt='r.')
        l2 = ax1.errorbar(xvals, (np.array(value1) - np.array(value1))/np.array(error1), yerr=np.ones(len(error1)), fmt='b.')
        ax1.set_xlim(0, len(par))
        ax1.legend([l2,l1], [f1,f2], loc=LegLoc, numpoints=1)
        #print ax1.get_ticks(xvals)
        #print ax1.get_ticklabels(xvals)
        plt.show()

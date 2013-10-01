from datatools.tempo import  *
from pylab import *
import sys, cPickle
from scipy.optimize import fmin_powell as minimize
#from scipy.optimize import leastsq as minimize

ParamMap = ['M31410A:EQUAD', 'M31410B:EQUAD', 'ABPP1410:EQUAD', 'ABPP2380:EQUAD', 'M41410:EQUAD', 'M42380:EQUAD', 'GASP-L:EMIN', 'GUPPI-L:EMIN', 'ASP-L:EMIN', 'PUPPI-L:EMIN', 'GASP-8:EMIN', 'GUPPI-8:EMIN', 'ASP-S:EMIN', 'PUPPI-S:EMIN']
#ParamMap = ['800:EMIN', 'L:EMIN']
#ParamMap = ['GASP-800:EMIN', 'GASP-L:EMIN','GUPPI-800:EMIN', 'GUPPI-L:EMIN']
parfile = '1713.Sep.Paul.par'
timfile = '1713.Sep.Paul.tim'
init_param = []
m = model(parfile)
t = TOAfile(timfile)
for i,PM in enumerate(ParamMap):
    grp,key = PM.split(':')
    toaex = t.toalist[t.groups[grp][0]]
    init_param.append(float(toaex.flags[key]))
print init_param
#init_param = cPickle.load(open('result.pkl', 'rb'))

def calChisq(params):
    for i,PM in enumerate(ParamMap):
        grp,key = PM.split(':')
        value = params[i]
        for j in t.groups[grp]:
            t.toalist[j].flags[key] = Decimal(value)
    #ts = t.subgroup(groups = t.groups.keys())
    with open('1713.adjust_hyperpar.toa', 'w') as fout:
        fout.write(t.tempo2fmt(tempo1use=True))
    ts = TOAfile('1713.adjust_hyperpar.toa')
    m.tempofit(ts)
    print 'overall:', m.dof, m.chisq

    SumFac = 0.
    Output = []
    for grp in t.groups:
        idx = t.groups[grp]
        dof = len(idx) 
        chisq_grp = (m.res[idx]**2/m.err[idx]**2).sum()/dof
        print grp, ':', chisq_grp
        Output.append(1.-chisq_grp)
        SumFac += (1.-chisq_grp)**2
        #toaex = t.toalist[t.groups[grp][0]]
        #toapars = {}
        #for k in ['EMIN', 'EQUAD']:
            #if toaex.flags.has_key(k):
                #toapars[k] = float(toaex.flags[k])
            #else:
                #toapars[k] = 0.
        #if toaex.flags.has_key('EFAC'):
            #toapars['EFAC'] = float(toaex.flags['EFAC'])
        #else:toapars['EFAC'] = 1.
        #SumFac += (toapars['EMIN']**2 + toapars['EQUAD']**2)*toapars['EFAC']**2

    #print m.chisq, (m.res**2/m.err**2).sum()

    #return Output
    return SumFac


print calChisq(init_param)
result = minimize(calChisq,init_param,xtol=0.05, ftol=0.05, maxiter=100, maxfun=500 )
#result,i = leastsq(calChisq,init_param)
#xopt, fopt, iter = result[:3]
#print xopt, fopt, iter
print result
cPickle.dump(result, open('result.pkl', 'wb'), protocol=2)

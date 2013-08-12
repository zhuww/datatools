from datatools.tempo import  *
from pylab import *
import sys, cPickle
from scipy.optimize import fmin_powell as minimize
from scipy.optimize import leastsq

ParamMap = ['M31410A:EQUAD', 'M31410B:EQUAD', 'ABPP1410:EQUAD', 'ABPP2380:EQUAD', 'M41410:EQUAD', 'M42380:EQUAD', 'ASP:EMIN', 'GASP:EMIN', 'GUPPI:EMIN', 'PUPPI:EMIN']
parfile = '1713.ext.moreDMX.par'
timfile = '1713.May.Weighted.tim'
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
    ts = t.subgroup(groups = t.groups.keys())
    m.tempofit(ts)

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
result = minimize(calChisq,init_param)
#result,i = leastsq(calChisq,init_param)
#xopt, fopt, iter = result[:3]
#print xopt, fopt, iter
print result
cPickle.dump(result, open('result.pkl', 'wb'), protocol=2)

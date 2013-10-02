"""
An code for running MCMC simulation to determine the confidence range of tempo parameters. Code in use:
    runmcmc.py : the main driving program
    plotmc.py : the plotting program
    ProgressBar.py : for plotting the progress bar
"""
from datatools.tempo import tempofit, tempo2fit, touchparfile, uniquename, PARfile #, model, TOAfile
from math import *
from decimal import *
import os
from copy import *
from numpy.random import normal , uniform ,seed
import numpy as np


#parfile = '1713.sns.par'
#toafile = '1713.sns.tim'
#chisq, dof = tempofit(parfile, toafile = toafile)
#smallestchisq = chisq

def randomnew(pf, stepsize): #special for 1713
    twopi = 6.283185307179586
    fac = 1.536e-16 * 1.e12
    x = float(str(pf.A1[0]))
    sini = float(str(pf.SINI[0]))
    cosi = np.sqrt(1 - sini**2)
    Omega = float(str(pf.PAASCNODE))
    m2 = float(str(pf.M2[0])) + normal(0,0.03*stepsize)
    cosi = cosi + normal(0, 0.02*stepsize)
    Omega = Omega + normal(0, 0.2*stepsize)
    mu = float(np.sqrt(float(str(pf.PMRA[0]**2+pf.PMDEC[0]**2))))
    if m2 <= 0 or Omega > 360 or Omega < -360 or cosi > 1.:
        return 0
    #sini = sqrt(1 - cosi**2)
    thetamu = 180. + np.arctan(float(str(pf.PMRA[0]/pf.PMDEC[0])))/np.pi*180
    xdot = -1.* fac * x * mu * (cosi/sini) * sin((thetamu-Omega)*twopi/360.)
    sini = np.sqrt(1 - cosi**2)
    pf.SINI[0] = Decimal(str(sini))
    pf.XDOT[0] = Decimal(str(xdot))
    pf.PAASCNODE = Decimal(str(Omega))
    pf.M2[0] = Decimal(str(m2))
    return pf

def probcal(pf):
    pf.write()
    #if m2 <= 0 or Omega > 360 or Omega < -360:
        #return 0
    #sini = sqrt(1 - cosi**2)
    #xdot = -1.* fac * x * mu * (cosi/sini) * sin((thetamu-Omega)*twopi/360.)
    #touchparfile(parfile, NITS=1, PAASCNODE=Omega, SINI = sini, M2 = m2, XDOT = xdot)
    chisq, dof = tempofit(parfile, toafile = toafile, pulsefile = pulsefile)
    pf.chisq = chisq
    #print dof, chisq
    #print parfile, toafile
    #print chisq, smallestchisq
    #if chisq >= 999999999.:return 0
    #smallestchisq = 308.28
    try:
        return exp((smallestchisq - chisq)/2.) #Ingrid/Paul?
    except OverflowError:
        print chisq, smallestchisq
        print pf.parfile
        raise OverflowError 

#print probcal(90, 0.25, 0.28)


from itertools import count
import cPickle as pickle
import os,sys
from tempfile import mkdtemp

#print probcal(iOmega, icosi, im2)

class MChain(object):
    def __enter__(self):
        #try:
            #Chain = pickle.load(open('MChain.p', 'r'))['Chain']
        #except:
            #Chain = []
        Chain = []
        self.Chain = Chain
        self.cwd = os.getcwd()
        return self.Chain
    def __exit__(self, exc_type, exc_value, exc_tb):
        os.chdir(self.cwd)
        try:
            MarkovChain = pickle.load(open('MChain.p', 'rb'))['Chain']
            MarkovChain.extend(self.Chain)
        except:
            MarkovChain = self.Chain
        if len(MarkovChain)>2: 
            if len(MarkovChain[-1]) < len(MarkovChain[-2]):
                MarkovChain = MarkovChain[:-1]
        dict = {'Chain':MarkovChain}
        pickle.dump(dict, open('MChain.p', 'wb'), protocol=2)

        if exc_type is KeyboardInterrupt:
            print '\nManually Stopped\n'
            return True
        else:
            return exc_type is None
        print '\nFinish running\n' 

def motifile(file, cwd, tmpdir):
    os.system('cp %s/%s %s/%s' % (cwd, file, tmpdir, file))
    text = ''
    f = open(file, 'rw')
    for l in f.readlines():
        if not l.find('INCLUDE') == -1:
            a = l.split()
            if a[0] == 'C' or a[0] =='#':
                continue
            if not open(cwd+'/'+a[1],'r').read().find('INCLUDE') == -1: 
                motifile(a[1], '..', '.')
                l = a[0] +' '+a[1]
            else:
                l = a[0] + ' '+cwd+'/'+a[1]
            if not l[-1] == '\n':
                l += '\n'
            text += l
        else:
            if not l[-1] == '\n':
                l += '\n'
            text += l
    f.close()
    f = open(file, 'w')
    f.write(text)
    f.close() #motify the tim file to make sure INCLUDE follow the right files.
    
from ProgressBar import progressBar
def mcmc(Chain, runtime, mixingtime=1000, stepsize=1):
    #mixingtime = 1000
    #runtime = 50000
    pb = progressBar(maxValue = runtime + mixingtime)
    cwd=os.getcwd()
    tmpdir = cwd+'/.'+uniquename()
    if not tmpdir == None:
        if os.path.exists(tmpdir):
            os.chdir(tmpdir)
        else:
            os.mkdir(tmpdir)
            os.chdir(tmpdir)
    os.system('cp %s/%s %s/%s' % (cwd, parfile, tmpdir, parfile))
    os.system('cp %s/%s %s/%s' % (cwd, pulsefile, tmpdir, pulsefile))
    motifile(toafile, cwd, tmpdir)
    MarkovChain = Chain
    pf = PARfile(parfile)
    #pf = model(parfile)
    #pf.thawall()
    pf.freezeall('DMX_0')
    #pf.parameters['SINI'] = '0'
    #pf.parameters['M2'] = '0'
    #pf.parameters['XDOT'] = '0'

    pf.write()
    chisq, dof = tempofit(parfile, toafile = toafile, pulsefile = pulsefile)
    #pf.tempofit(toafile, pulsefile = pulsefile)
    chisq, dof = chisq, dof
    pf.matrix(toafile)
    #pf.freezeall()
    #pf.thawall('JUMP_')
    pf.write()

    #if 'PAASCNODE'in pf.__dict__:
        #plist = [x for x in pf.manifest if x in pf.parameters.keys()] + ['PAASCNODE']
        #dict = {'BEST':[pf.__dict__[p][0] for p in plist[:-1]] + [pf.__dict__[p] for p in plist[-1:]], 'parfile':pf.parfile, 'parameters':plist + ['chisq']}
    #else:
    plist = [x for x in pf.manifest if x in pf.parameters.keys()]

    dict = {'BEST':[pf.__dict__[p][0] for p in plist] + [pf.PAASCNODE, chisq], 'parfile':pf.parfile, 'parameters':plist + ['PAASCNODE', 'chisq']}
    pickle.dump(dict, open('%s/bestpar.p' % cwd, 'w'), protocol=2)
    p0 = probcal(pf)
    p = p0
    #print 'P0', p0
    #try:
        #MChain = pickle.load(open('MChain.p', 'r'))
    #except:
        #MChain = {'Chain':[]}
    #MChain['parameters'] = plist
    #pickle.dump(MChain, open('MChain.p', 'w'))
    n = count()
    m = count()
    while n.next() <= mixingtime + runtime:
        npf = pf.randomnew(stepsize=stepsize)
        randomnew(npf, stepsize)
        p1 = probcal(npf)
        c = m.next()
        if c % 30 == 0:pb(c)
        if p1 > p0:
            if c > mixingtime:
                MarkovChain.append([npf.__dict__[p][0] for p in plist] + [npf.PAASCNODE, npf.chisq])
            pf = npf
            p0 = p1
            if p1 > p:
                p = p1
                #if 'PAASCNODE' in plist:
                    #dict['BEST'] = [pf.__dict__[p][0] for p in plist[:-1]] + [pf.__dict__[p] for p in plist[-1:]] + [npf.chisq]
                #else:
                dict['BEST'] = [npf.__dict__[p][0] for p in plist] + [npf.PAASCNODE,  npf.chisq]
                pickle.dump(dict, open('%s/bestpar.p' % cwd, 'w'), protocol=2)
        else:
            t = uniform(0,1,1)[0]
            if t < p1/p0:
                if c > mixingtime:
                    MarkovChain.append([npf.__dict__[p][0] for p in plist] + [npf.PAASCNODE, npf.chisq])
                pf = npf
                p0 = p1
            else:
                if c > mixingtime:
                    MarkovChain.append([pf.__dict__[p][0] for p in plist] + [npf.PAASCNODE, npf.chisq])
    #print  MarkovChain
    #print best
    print '%d points added to the Chain.' % len(MarkovChain)
    #OldChain = pickle.load(open('MChain.p','r'))['Chain']
    #dict['Chain'] = OldChain + MarkovChain
    #dict['Chain'] = MarkovChain
    os.chdir(cwd)
    #os.rmdir(tmpdir)

        
from optparse import OptionParser
if __name__ == '__main__':
    #main()
    usage = "usage: %prog [options] arg"
    parser = OptionParser()
    parser.add_option("-f", '--parfile', dest="parfile", help="par file")
    parser.add_option("-t", '--timfile', dest="toafile", help="toa file")
    parser.add_option("-n", '--pulsefile', dest="pulsefile", help="pulse number file", default=None)
    parser.add_option("-i", '--iter', type='int', nargs=1, dest='steps', help="number of steps")
    parser.add_option("-m", '--mixing', type='int', nargs=1, dest='mixing', help="number of mixing steps", default=1000)
    parser.add_option("-p", '--parallel', type='int', nargs=1, dest='paral', help="number of parallel processes")
    parser.add_option("-s", '--seed', type='int', nargs=1, dest='seed', default=int(os.getpid()), help="random number seed")
    parser.add_option("-z", '--stepsize', type='float', nargs=1, dest='stepsize', default=1., help="step size")
    (options, args) = parser.parse_args(args=sys.argv[1:])
    print options

    parfile = options.parfile
    toafile = options.toafile
    pulsefile = options.pulsefile
    steps = options.steps
    mixing = options.mixing
    rseed = options.seed
    stepsize = options.stepsize
    px = options.paral
    pf =PARfile(parfile)
    #pf = model(parfile)
    pf.freezeall()
    pf.thawall('JUMP_')
    pf.write('mcmc.par')
    touchparfile('mcmc.par', NITS=1)
    chisq, dof = tempofit('mcmc.par', toafile = toafile, pulsefile = pulsefile)
    #pf.tempofit(TOAfile(toafile), pulsefile = pulsefile)
    smallestchisq = chisq
    #print 'smallestchisq', smallestchisq, dof

    def run(s):
        seed(s) # assigning different initial seed for the random number generator in different threads.
        #steps = 50000
        with MChain() as Chain:
            #print steps, tmpdir
            mcmc(Chain, steps, mixingtime=mixing, stepsize=stepsize)
        return Chain

    if px == None:
        run(rseed)
    else:
        px = options.paral
        from multiprocessing import Pool
        p = Pool(px)
        p.map(run, range(rseed, rseed+px))





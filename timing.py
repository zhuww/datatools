'''A module that handles timing analysis. Start with a list of photon arrival times.'''
import numpy, os
from numpy import mean, std, var, array, double, angle
from math import *
from fileio import *
from SaveLoadable import MetaSaveLoader
from decimal import *
#from scipy.stats import chisqprob as chi2prob

#secperday = Decimal('86400')
secperday = 86400

def phase(toa, F0, F1, PEPOCH, CXOREF, F2=None):
    if isinstance(toa, (float, double, Decimal)): 
        toa = toa - (PEPOCH-CXOREF)*secperday
        if F2 == None:
            return (F0*toa + F1*toa**2/2) % 1
        else:
            return (F0*toa + F1*toa**2/2 + F2*toa**3/6) % 1
    elif isinstance(toa, (list, numpy.ndarray)):
        result = []
        for event in toa:
            if isinstance(event, Decimal) or isinstance((PEPOCH-CXOREF), Decimal):
                event = Decimal(str(event)) - (PEPOCH-CXOREF)*secperday
            else:
                event = event - (PEPOCH-CXOREF)*86400
            if F2 == None:
                phs = (F0*event + F1*event**2/2) % 1
            else:
                phs = (F0*event + F1*event**2/2 + F2*event**3/6) % 1
            if phs > 0 or phs == 0:
                result.append(phs)
            else:
                result.append(1+phs)
        return result
    else:
        print toa.__class__
        raise TypeError
            

def fold(phaselist, Nbin, file=None):
    lc = [0] * Nbin
    for phase in phaselist:
        lc[int(phase*Nbin)] += 1
    if not file == None:
        phase = (array(range(1,Nbin+1))-0.5)/Nbin
        lcerr = []
        for x in lc:
            lcerr.append(sqrt(float(x)))
        tofile(file, phase, lc, lcerr)
    return lc

def constfit(lc):
    N = len(lc)
    try:
        c = N/sum([1./x for x in lc])
    except ZeroDivisionError:
        print 'Warning: two few counts in some bins, unable to calculate chi^2 statistics.'
        c =  N/sum([1./x for x in lc if x > 0])
    chisq = lambda c: sum([(x-c)**2/x for x in lc if x > 0])
    return chisq(c)

#def chisqprob(chisq, dof):
    #return chi2prob(chisq, dof)

def Z2test(m, phaselist):
    PI = 3.14159265
    ''' Perform the Z^2_m test for the photons of phase list, where m is the number of harmonics. '''
    N = len(phaselist)
    Z2 = (2./N)*sum([sum([ cos(2.*PI*k*phs)for phs in phaselist])**2 + sum([ sin(2.*PI*k*phs) for phs in phaselist])**2 for k in range(1,m+1)])
    return Z2

def Htest(phaselist):
    '''Calculate the H test value Max(Zm^2-4m+4 for m < 23). '''
    return max([Z2test(m, phaselist)-4*m+4 for m in range(1,24)])

def Hprob(HorPhs):
    '''Calculate the null-hypothesis  probability based on the H value. '''
    if isinstance(HorPhs, (float,int)):
        return exp(-0.39802*HorPhs)
    elif isinstance(HorPhs, (list, tuple)):
        return exp(-0.39802*Htest(HorPhs))
    else:
        raise TypeError

def fold_ph(toa, F0, F1, PEPOCH, CXOREF, F2=None, Nbin = 16, file=None):
    return fold(phase(toa, F0, F1, PEPOCH, CXOREF, F2), Nbin, file)

from round import *
def shorten(value, error):
    #fig = figure((float(str(value)), float(str(error))))
    #fig = figure((Decimal(str(value)), Decimal(str(error))))
    #valuestr = str(fig)
    #first, last = valuestr.split(')')
    #head, toe = first.split('(')
    #result = head+last
    #return Decimal(result)
    #return str(fig)
    return Decimal(str(value)).quantize(Decimal(str(error)))

class PulseProfile(object):
    def __init__(self, lc):
        self.profile = lc
        Nbin = len(lc)
        self.Nbin = Nbin
        self.chisq = constfit(lc)
        self.phase = (array(range(1,Nbin+1))-0.5)/Nbin

    def sinfit(self):
        lc = self.profile
        powspec = numpy.fft.rfft(lc)
        self.powspec = powspec
        self.PPhase = angle(powspec[1], deg=True)/-360.
        for i in range(2, len(powspec)):
            powspec[i] = 0.+0.j
        self.sinmodel = numpy.fft.irfft(powspec)
        lcerr = [sqrt(x) for x in lc]
        errspec = numpy.fft.rfft(lcerr)
        sigma = abs(errspec[1])
        pow1 = abs(powspec[1])
        self.PPhase_err = asin(sigma/pow1)/3.14159265/2.

    def plot(self, fig=None, color='k'):
        def pl(fig, lc):
            Nbin = len(lc)
            lc+=lc
            lcerr = [sqrt(x) for x in lc]
            phs = (array(range(1,Nbin*2+1))-0.5)/Nbin
            Stepphs = (array(range(0,Nbin*2+1)))*1./Nbin
            fig.errorbar(phs, lc, lcerr, fmt=color+'.')
            fig.step(Stepphs, [0] + lc, c=color)
            if self.__dict__.has_key('sinmodel'):
                Nbin = len(self.sinmodel)
                phs = (array(range(1,Nbin*2+1))-0.5)/Nbin
                model = list(self.sinmodel)+list(self.sinmodel)
                fig.plot(phs, model, color+'-')
            try:
                fig.set_ylabel('Counts/bin')
            except:pass
        lc = self.profile
        if fig == None:
            print '!!!Use subplot to creat a fig object first, then use PulseProfile.plot(fig, color=...)'
            return
            #pl(fig,lc)
            #try:
                #xlabel('Phase')
            #except:pass
            #show()
        else:
            pl(fig,lc)
            try:
                fig.set_xlabel('Phase')
            except:pass


class ephemeris(object):
    __metaclass__ = MetaSaveLoader
    def __init__(self, file):
        """A class object to contain the APIs for extracting ephemeris information from a .par file."""
        if os.access(file, os.R_OK):
            self.parfile = file
        else: 
            raise FileError(file) 
        def floatify(val):
            try: return Decimal(val)
            except: 
                import re
                p = re.compile('-*\d+\.*\d*D[-\+]\d+', re.VERBOSE)
                if p.match(val):
                    return Decimal(val.replace('D','E'))
                else:
                    return str(val)
        self.manifest = []
        self.file = open(file, 'r')
        for lines in self.file.readlines():
            items = lines.split()
            if len(items) == 2 or len(items) == 3:
                self.__dict__[items[0]] = floatify(items[1])
                self.manifest.append(items[0])
            elif len(items) == 4:
                value = floatify(items[1])
                error = floatify(items[3])
                if all([isinstance(value, Decimal), isinstance(error, Decimal)]):
                    #self.__dict__[items[0]] = value
                    #self.__dict__[items[0]] = shorten(items[1],items[3])
                    self.__dict__[items[0]] = shorten(value,error)
                else:
                    self.__dict__[items[0]] = (value, error)
                self.manifest.append(items[0])
            elif len(items) == 0:pass
            elif items[0] == '#':pass
            else:
                print items
                raise IndexError
        if self.__dict__.has_key('PSR'):
            self.psrname = self.PSR
        elif self.__dict__.has_key('PSRJ'): 
            self.psrname = 'J'+self.PSRJ
        elif self.__dict__.has_key('PSRB'): 
            self.psrname = 'B'+self.PSRB
        else:
            raise 'Cant Find PSR name'
        self.file.close()

    def shift(self, TZMJD):
        TZMJD = Decimal(str(TZMJD))
        t = (TZMJD - self.PEPOCH)*secperday
        F0 = self.F0 + self.F1*t + self.F2*t**2/2
        F1 = self.F1 + self.F2*t
        self.PEPOCH = TZMJD
        self.F0 = F0
        self.F1 = F1

    def phase(self, toa, CXOREF, opt=None):
        if opt == None:
            if self.__dict__.has_key('F2'):
                F2 = self.F2
            else:
                F2 = None
            return phase(toa, self.F0, self.F1, self.PEPOCH, CXOREF, F2=F2)
        elif opt == 'Tempo2':pass
        else:
            raise NotImplimented
        
    def fold(self, toa, CXOREF, file=None, Nbin=16):
        if self.__dict__.has_key('F2'):
            F2 = self.F2
        else:
            F2 = None
        return fold_ph(toa, self.F0, self.F1, self.PEPOCH, CXOREF, F2=F2, Nbin=Nbin, file=file)

    def XrayTOA(self, toa, CXOREF, Nbin=16, fmt=None):
        def unify(value):
            if isinstance(CXOREF, (Decimal, int)) and not isinstance(value, Decimal):
                return Decimal(str(value))
            elif isinstance(CXOREF, (float, int)) and not isinstance(value, float):
                return float(str(value))
            else:
                return value
        PEPOCH = unify(self.PEPOCH)
        #absphase = []
        phaselist = []
        toffset = []
        for t in toa:
            t = unify(t) - (PEPOCH-CXOREF)*secperday
            toffset.append(t)
            F0 = unify(self.F0)
            F1 = unify(self.F1)
            if 'F2' in self.__dict__:
                F2 = unify(self.F2)
                p = (F0*t + F1*t**2/2 + F2*t**3/6)
            else:
                p = (F0*t + F1*t**2/2)
            #absphase.append(p)
            phs = p % 1
            if phs > 0 or phs == 0:
                phaselist.append(phs)
            else:
                phaselist.append(1+phs)
        lc = fold(phaselist, Nbin=Nbin)
        PP = PulseProfile(lc)
        PP.sinfit()
        PPhase = unify(PP.PPhase)
        PPhase_err = unify(PP.PPhase_err)
        XTOA = []
        for i in range(len(toa)):
            Dt = (PPhase - phaselist[i])/(F0+F1*toffset[i]+F2*toffset[i]**2/2)
            terr = PPhase_err/(F0+F1*toffset[i]+F2*toffset[i]**2/2)
            if fmt == 'Tempo':
                XTOA.append(' xxxxxxxxxxxxxxxxx   0  0 0000.000  %s    0.00 %7.2f        @' % (str(CXOREF + (unify(toa[i]) + Dt)/secperday)[:20].rjust(20, ' '), terr*1000 ))
            elif fmt == 'Tempo2':
                XTOA.append('@                  0.000 %s %.2f' % (str(CXOREF + (unify(toa[i]) + Dt)/secperday)[:20].rjust(20, ' '), terr*1000000 ))
            else:
                XTOA.append((CXOREF + (unify(toa[i]) + Dt)/secperday, terr/secperday))
        return XTOA

class TOAfile(object):
    def __init__(self, fn):
        with open(fn, 'rt') as f:
            l = 0
            text = f.read()
            for line in text.split('\n'):
                if not line == '' and not line[0] == '#' and not line[0] == 'C':
                    if len(line) > l:
                        l = len(line)
                        lline = line
            self.longestline = lline
            if  lline[1] == ' ':
                self.format = 'Princeton'
            elif lline[0] == ' ':
                self.format = 'Parkes'
            elif not (lline[0]==' ' or lline[1]==' '):
                self.format = 'ITOA'
            else:
                self.format = 'Unknown'

            if self.format == 'Princeton':
                self.observotary = int(self.longestline[0])
                self.freq = float(self.longestline[15:24])
                self.TOAstart = 1000000.
                self.TOAend = 0.
                for line in text.split('\n')[:-1]:
                    if not line == '' and not line[0] == '#' and not line[0] == 'C' and line[1]==' ':
                        toa = float(line[24:44])
                        if self.TOAstart > toa:
                            self.TOAstart = toa
                        if self.TOAend < toa:
                            self.TOAend = toa

            if self.format == 'ITOA':
                self.observotary = self.longestline[57:59]
                self.freq = float(self.longestline[34:45])
                self.TOAstart = 1000000.
                self.TOAend = 0.
                for line in text.split('\n')[:-1]:
                    if not line[0:2] == '  ':
                        toa = float(line[9:28])
                        if self.TOAstart > toa:
                            self.TOAstart = toa
                        if self.TOAend < toa:
                            self.TOAend = toa
                        
            if self.format == 'Parkes':
                self.observotary = int(self.longestline[79])
                self.freq = float(self.longestline[25:34])
                self.TOAstart = 1000000.
                self.TOAend = 0.
                for line in text.split('\n')[:-1]:
                    if line[0] == ' ':
                        toa = float(line[34:55])
                        if self.TOAstart > toa:
                            self.TOAstart = toa
                        if self.TOAend < toa:
                            self.TOAend = toa



            

    #def plot(self):
        #from PyXGraph.pyxgraph import *
        #from datetime import *
        #from MJD import *
        #from numpy import *
        #def convert_to_path(g, x, y, close_path=True):
            #"""Convert the coordinates from the arrays x and y to a (closed) path."""
            #p = pyx.path.path()
            #p.append(pyx.path.moveto(*g.pos(x[0], y[0])))
            #for i in xrange(1, len(x)):
                #p.append(pyx.path.lineto(*g.pos(x[i], y[i])))
            #if close_path:
                #p.append(pyx.path.closepath())
            #return p
        #MJD = 



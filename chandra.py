from Data import *
from fitstools import *
import os, sys, copy, commands, re
#from subprocess import call, Popen
from Exceptions import *
from math import cos, sqrt
import tempfile
from fileio import *
from tools.PyATNF import *
from tools.Coordinate import RA, Dec

CIAO_CMD='. /usr/local/ciao/bin/ciao.sh'

def lookforCIAO():
    PATH = os.getenv('PATH')
    index = PATH.find('ciao')
    if index == -1: 
        sys.exit('***Error: You need to source ciao before using it in python.')
        #p = call(['/bin/bash', CIAO_CMD], shell=True)
    else: print 'Using CIAO:', PATH[index:index+8]

class MetaCiaoEnv(type):
    pass

class ciao(chandraData):
    '''A class that provide the chandra data with the proper analysing tools.'''
    def __init__(self, Data):
        if Data == 'dummy': return None
        if isinstance(Data, basestring):
            chandraData.__init__(self, Data)
        elif isinstance(Data, chandraData):
            for (attr, value) in sorted(Data.__dict__.iteritems()):
                self.__setattr__(attr, value)
        else:
            raise TypeError(Data)
        self.gainfile = fitsfile(self.evtfile).hdread('GAINFILE', block=1)
        self.object = fitsfile(self.evtfile).hdread('Object', block=1).replace('-','m').replace('+','p')

    def verify(self, file):
        if os.access(file, os.R_OK):
            return file
        elif os.access(self.path+'/primary/'+file, os.R_OK):
            return self.path+'/primary/'+file
        else:
            raise FileError(file)

    def psextract(self, src='src.reg', bg='bg.reg', root=None, gtype="NUM_CTS", gspec=15, clobber='no', psr=None):
        if root == None:root = self.object
        try:
            src = self.verify(src)
            bg = self.verify(bg)
        except:
            self.psdetect(psr=psr)
            region = copy.deepcopy(self.region)
            region.regions = region.regions[0:1]
            region.regions[0].shape='circle'
            region.regions[0].R = region.regions[0].Rout
            region.save(self.path+'/'+src)
            region.regions[0].shape='annulus'
            region.regions[0].Rout *= 5
            region.save(self.path+'/'+bg)

        #check_ctiapp.sh takes care the CTI_APP keyword. See http://cxc.harvard.edu/ciao/why/cti.html for detail.
        cmdline = '''
        chmod +w %(evtfile)s
        check_ctiapp.sh %(evtfile)s 
        psextract events="%(evtfile)s[sky=region(%(src)s)]" bgevents="%(evtfile)s[sky=region(%(bg)s)]" root=%(root)s asol="%(asolfile)s" bgasol="" pbkfile="%(pbkfile)s" gtype=%(gtype)s gspec=%(gspec)r clobber=%(clobber)s''' % {
                'evtfile':self.evtfile,
                'src':src, 
                'bg':bg,
                'root':root,
                'asolfile':self.asolfile,
                'pbkfile':self.pbkfile,
                'gtype':gtype, 
                'gspec':gspec,
                'clobber':clobber,}
        os.system(cmdline)
        self.mkacisrmf(root=root)
        path = os.getcwd()
        grpspec = '''grppha infile="%(root)s.pi" outfile="%(root)s_grp.fits" comm="chkey BACKFILE %(path)s/%(root)s_bg.pi&chkey RESPFILE %(path)s/%(root)s.rmf&chkey ANCRFILE %(path)s/%(root)s.arf&group min %(gspec)s&exit !%(root)s_grp.fits"''' % {'root':root, 'path':path, 'gspec':gspec}
        os.system(grpspec)


    def pssharp(self, evtfile=None, outfile='sharperimage.fits', root=None):
        if root == None:root = self.object
        if evtfile == None:evtfile=self.evtfile
        cmdline = """acis_process_events infile=%(infile)s outfile=%(outfile)s \
                    acaofffile=%(asolfile)s rand_pha=no doevtgrade=no gradefile=none \
                    calculate_pi=no rand_pix_size=0.0 \
                    eventdef='{d:time,s:ccd_id,s:node_id,i:expno,s:chip,s:tdet,f:det,f:sky,l:pha,f:energy,l:pi,s:fltgrade,s:grade,x:status}' \
                    clobber=yes""" % { 'infile':evtfile, 'outfile':self.path+'/primary/'+outfile, 'asolfile':self.asolfile}
        os.system(cmdline)
        return outfile

    def psdetect(self, evtfile=None, Energy=None, psr=None):
        if evtfile == None:evtfile=self.evtfile
        if not Energy == None:
           newevt = tempfile.NamedTemporaryFile(mode='w+b',suffix='.fits').name 
           filterevents(evtfile, outfile=newevt, Emin=Energy[0], Emax=Energy[1])
           evtfile=newevt
        cmdline = '''celldetect infile=%(infile)s outfile=%(outfile)s regfile=%(regfile)s clobber=yes''' % {
                'infile':evtfile,
                'outfile': 'srclist.fits',
                'regfile': 'srclist.reg'
                }
        print cmdline
        os.system(cmdline)
        file = fitsfile('srclist.fits')
        table = list(file.gettable())
        cols = file.colinfo()
        for i in range(len(table)):
            table[i] = (table[i], i)

        if not psr == None:
            ra, dec = QatnfPos(psr, in_unit_degree=True)
            def Dist(src):
                srcra = src[0][cols.names.index('RA')]
                srcdec = src[0][cols.names.index('DEC')]
                return sqrt((ra-srcra)**2+(dec-srcdec)**2)*3600 # Distance in arcsecond.

            table = sorted(table, key=Dist, reverse=True)
        else:
            table = sorted(table, key=lambda table: table[0][cols.names.index('SNR')])
        srctable = table[-1][0]
        srcidx = table[-1][1]
        psrname = lambda psr:psr.replace('+','p').replace('-','m')
        regions = open('srclist.reg', 'r').readlines()
        self.region = psregion(text=regions[srcidx])
        if not psr == None:
            self.region.save(self.path+'/'+psrname(psr)+'.reg')
        else:
            self.region.save(self.path+'/'+'MSsrc'+'.reg')
        cols = file.colinfo()
        RA = srctable[cols.names.index('RA')]
        RAERR = srctable[cols.names.index('RA_ERR')]
        DEC = srctable[cols.names.index('DEC')]
        DECERR = srctable[cols.names.index('DEC_ERR')]
        X = srctable[cols.names.index('X')]
        XERR = srctable[cols.names.index('X_ERR')]
        Y = srctable[cols.names.index('Y')]
        YERR = srctable[cols.names.index('Y_ERR')]
        if not self.__dict__.has_key('position'):
            self.position = {}
        self.position['RA'] = RA
        self.position['RA_ERR'] = RAERR
        self.position['DEC'] = DEC
        self.position['DEC_ERR'] = DECERR
        self.position['X'] = X
        self.position['X_ERR'] = XERR
        self.position['Y'] = Y
        self.position['Y_ERR'] = YERR
        if not self.position.has_key('ra'): 
            self.position['ra'] = RA
        if not self.position.has_key('dec'): 
            self.position['dec'] = DEC
        print self.position
        return table


    def coordinates(self, src=None):
        if src:
            output = commands.getoutput("cat %s | grep -v '#'" % src)
            shape=re.split(r'\(',output)
            between=re.split(r'\)',shape[1])
            number=re.split(r',',between[0])
            X = float(number[0])
            Y = float(number[1])
            output = commands.getoutput("dmcoords %(evt)s asolfile=%(asl)s  option=sky x=%(X)f y=%(Y)f celfmt=deg; pget dmcoords chip_id chipx chipy theta phi x y ra dec" % {'evt':self.evtfile, 'asl':self.asolfile, 'X':X, 'Y':Y})
            print output
            coords = output.split()
            if not self.__dict__.has_key('position'):
                self.position = {}
            self.position['chip_id'] = int(coords[0])
            self.position['chipx'] = float(coords[1])
            self.position['chipy'] = float(coords[2])
            self.position['theta'] = float(coords[3])
            self.position['phi'] = float(coords[4])
            self.position['X'] = float(coords[5])
            self.position['Y'] = float(coords[6])
            self.position['ra'] = coords[5]
            self.position['dec'] = coords[6]
        else:
            output = commands.getoutput("dmcoords %(evt)s asolfile=%(asl)s  option=cel ra=%(ra)f dec=%(dec)f celfmt=deg; pget dmcoords chip_id chipx chipy theta phi" % {'evt':self.evtfile, 'asl':self.asolfile, 'ra':self.position['ra'], 'dec':self.position['dec']})
            print output
            coords = output.split()
            if not self.__dict__.has_key('position'):
                self.position = {}
            self.position['chip_id'] = int(coords[0])
            self.position['chipx'] = float(coords[1])
            self.position['chipy'] = float(coords[2])
            self.position['theta'] = float(coords[3])
            self.position['phi'] = float(coords[4])
            #self.position['X'] = float(coords[5])
            #self.position['Y'] = float(coords[6])
            #self.position['ra'] = coords[5]
            #self.position['dec'] = coords[6]


    def mkacisrmf(self, root=None, src='src.reg', bg='bg.reg', clobber='yes'):
        srcreg = self.verify(src)
        bgreg = self.verify(bg)
        if root == None:root = self.object
        CALresp = self.gainfile.replace('gain_cti','p2_resp')
        self.coordinates(src=srcreg)
        chipx = self.position['chipx']
        chipy = self.position['chipy']
        chip_id = self.position['chip_id']
        cmdline = '''mkacisrmf \
      infile="$CALDB/data/chandra/acis/p2_resp/%(resp)s" \
      outfile=%(root)s.rmf \
      energy=0.3:10.0:0.005 \
      channel=1:1024:1 \
      chantype=PI \
      wmap=none \
      ccd_id=%(id)r chipx=%(chipx)r chipy=%(chipy)r \
      gain=$CALDB/data/chandra/acis/det_gain/%(gain)s \
      clobber=%(clobber)s''' % { 'id':chip_id,
              'chipx':chipx,
              'chipy':chipy,
              'clobber':clobber,
              'resp':CALresp,
              'gain':self.gainfile,
              'root':root
              }
        os.system(cmdline)
        cmdmkarf = commands.getoutput("dmhistory %s.arf tool=mkarf" % (root)).split()
        for i in range(len(cmdmkarf)):
            if cmdmkarf[i].startswith('engrid'):
                cmdmkarf[i] = 'engrid="grid(%s.rmf[MATRIX][cols ENERG_LO,ENERG_HI])"' % (root)
            elif cmdmkarf[i].startswith('ENERG_LO'):
                cmdmkarf[i] = ''
            elif cmdmkarf[i].startswith('clobber'):
                cmdmkarf[i] = 'clobber=yes'
        os.system(' '.join(cmdmkarf))

    def marx(self, psffile=None):
        #if psffile == None: sys.exit('***Error: To run marx, you must specify a psf file.')
        if psffile == None:
            psffile = os.getcwd()+'/'+commands.getoutput('ls *.fits')
        (psfdir, file) = os.path.split(psffile)
        if psfdir == '':psfdir=os.getcwd()
        fileroot = ''.join(file.split('.')[:-1])
        if not os.access(psfdir+'/marx.par', os.R_OK): os.system('cp /usr/local/share/marx/pfiles/marx.par %s' % psfdir)
        #get detector_nominal_position
        #ra_nom=fitsfile(self.evtfile).hdread('RA_NOM',block=1)
        #dec_nom=fitsfile(self.evtfile).hdread('DEC_NOM',block=1)
        #roll_nom=fitsfile(self.evtfile).hdread('ROLL_NOM',block=1)
        (ra_nom, dec_nom, roll_nom) =fitsfile(self.evtfile).hdread(('RA_NOM','DEC_NOM','ROLL_NOM'),block=1)
        instrument=fitsfile(self.evtfile).hdread('INSTRUME')
        if instrument == 'ACIS' and self.position['chip_id'] > 4:
            detector = 'ACIS-S'
        #print ra_nom, dec_nom, roll_nom, detector
        if os.access(psfdir+'/psf.fits', os.R_OK): os.system('mv %(path)s/psf.fits %(path)s/psf.fits.old' % {'path':psfdir} )
        marxcmdline = '''
pset %(path)s/marx SAOSACFile=%(psffile)s
pset %(path)s/marx OutputDir=%(outputdir)s
pset %(path)s/marx DitherModel=INTERNAL 
pset %(path)s/marx DitherBlur=0.2
pset %(path)s/marx SourceType=SAOSAC

pset %(path)s/marx RA_Nom=%(ra_nom)f
pset %(path)s/marx Dec_Nom=%(dec_nom)f
pset %(path)s/marx Roll_Nom=%(roll_nom)f
pset %(path)s/marx SourceRA=%(ra)f
pset %(path)s/marx SourceDEC=%(dec)f

pset %(path)s/marx DetectorType=%(detector)s
pset %(path)s/marx GratingType=NONE
pset %(path)s/marx ExposureTime=0.0
marx @@%(path)s/marx.par
marx2fits %(path)s/%(outputdir)s %(path)s/psf.fits
''' % {'path':psfdir,
       'psffile':psffile,
       'outputdir':fileroot+'.dir',
       'ra':self.position['ra'],
       'dec':self.position['dec'],
       'ra_nom':ra_nom,
       'dec_nom':dec_nom,
       'roll_nom':roll_nom,
       'detector':detector,
       }
        os.system(marxcmdline)
        (obsSIM_X, obsSIM_Y, obsSIM_Z) = fitsfile(self.evtfile).hdread(('SIM_X','SIM_Y','SIM_Z'),block=1)
        (psfSIM_X, psfSIM_Y, psfSIM_Z) = fitsfile(psfdir+'/psf.fits').hdread(('SIM_X','SIM_Y','SIM_Z'),block=1)
        DSIM_X = obsSIM_X - psfSIM_X
        DSIM_Y = obsSIM_Y - psfSIM_Y
        DSIM_Z = obsSIM_Z - psfSIM_Z
        newmarxcmd = '''
        pset %(path)s/marx DetOffsetX = %(DSIM_X)s
        pset %(path)s/marx DetOffsetY = %(DSIM_Y)s
        pset %(path)s/marx DetOffsetZ = %(DSIM_Z)s
        pset %(path)s/marx OutputDir=%(outputdir)s
        marx @@%(path)s/marx.par
        marx2fits %(path)s/%(outputdir)s %(path)s/offset_psf.fits
        ''' % {'path': psfdir,
                'outputdir':'offset_'+fileroot+'.dir',
                'DSIM_X':DSIM_X,
                'DSIM_Y':DSIM_Y,
                'DSIM_Z':DSIM_Z,
                }
        os.system(newmarxcmd)

        
    def reproject(self, outfile, infile=None, aspect=None, match=None):
        '''By default, this method reproject the image to match its own best source position, one could use this to project other data image to match the best source position of this data or project this data to match the best source position of other data set.'''
        if infile == None: infile=self.evtfile
        if aspect == None: aspect = self.asolfile
        if match == None: match='%f %f' % (self.position['ra'], self.position['dec'])
        cmdline = 'reproject_events infile=%(infile)s outfile=%(outfile)s aspect=%(aspect)s match="%(match)s" clobber=yes' % {
                'infile':infile,
                'outfile':outfile,
                'aspect':aspect,
                'match':match,
                }
        os.system(cmdline)

    def align_point_source(self, other, outfile=None, position=None):
        '''Adjust the aspect solution file of the other dataset to align its best source position to that of this one.'''
        if not self.__dict__.has_key('imgref'):
            imgref = self.path+'/primary/centered_%s.fits' % self.ObsID
            self.reproject(outfile=imgref)
            self.imgref = imgref
        if position == None:
            (RA, DEC) = self.position['ra'], self.position['dec']
        else:
            (RA, DEC) = position
        (ra, dec) = other.position['ra'], other.position['dec']
        Delta_ra = RA - ra
        Delta_dec = DEC - dec
        delta_ra = Delta_ra * cos(DEC * 2 *3.14159625/360)
        delta_x = delta_ra / -1.3667e-4 #1 pixel = 0.492" = 1.3667e-4 deg
        delta_y = Delta_dec / 1.3667e-4
        newasolfile = other.path+'/primary/newaspectsol.fits'
        cmdline = 'punlearn wcs_update; wcs_update transformfile="" infile=%(infile)s outfile=%(outfile)s wcsfile=%(wcsfile)s deltax=%(delta_x)f deltay=%(delta_y)f clobber=yes' %        {'infile':other.asolfile,
             'outfile':newasolfile,
             'wcsfile':self.imgref,
             'delta_x':delta_x,
             'delta_y':delta_y}
        print cmdline
        os.system(cmdline)
        self.reproject(infile=other.evtfile, outfile=other.path+'/primary/newimg.fits', aspect=other.asolfile) # apply aspect solution to the other data set and shift its tangent point to the source position of this data set.
        if outfile == None:outfile = '%s/primary/aligned_with_%d.fits' % (other.path,self.ObsID)
        self.reproject(infile=other.path+'/primary/newimg.fits', outfile=outfile, match='none', aspect=newasolfile) # apply the corrected aspect solution to the reprojected image of the other data.

    @staticmethod
    def merge(infile, outfile):
        #cmdline = 'dmmerge infile="%s" outfile="%s" lookupTab="/usr/local/ciao-4.2/data/dmmerge_header_lookup.txt" clobber=yes' % (infile, outfile)
        cmdline = 'dmmerge infile="%s" outfile="%s" lookupTab="/usr/local/ciao/data/dmmerge_header_lookup.txt" clobber=yes' % (infile, outfile)
        print cmdline
        os.system(cmdline)

    def barycenter(self, outfile=None):
        if outfile == None:
            outfile = str(self.ObsID)+'barycorr.fits'
        orbitfiles = commands.getoutput('ls %s/primary/orbit*eph1.fits' % self.path)
        if not orbitfiles.find('\n') == -1:
            tstart = fitsfile(self.evtfile).hdread('TSTART')
            ephstart = 0.
            for orbit in orbitfiles.split('\n'):
                nexteph = fitsfile(orbit).hdread('TSTART')
                if nexteph > ephstart and nexteph < tstart:
                    ephstart = nexteph
                    orbeph = orbit
        else:
            orbeph = orbitfiles
        cmdline = 'axbary infile=%(infile)s orbitfile=%(orbitfile)s outfile=%(outfile)s ra=%(ra)f dec=%(dec)f clobber=yes' % {
                'infile':self.evtfile,
                'orbitfile':orbeph,
                'outfile':outfile,
                'ra':self.position['ra'],
                'dec':self.position['dec']}
        print cmdline
        os.system(cmdline)

    def UpdateRegion(self, region, evtfile=None):
        print 'Before update, the region file reads:\n%s' % (region)
        self.psdetect(evtfile=evtfile)
        region.X = self.position['X']
        region.Y = self.position['Y']
        print 'After update, the region file reads:\n%s' % (region)

    def RadioProfile(self, N, R, outfile=None, bg=None):
        try:
            Table = filterevents(self.evtfile, cols=['x','y'])
        except:
            Table = filterevents(self.evtfile, cols=['X','Y'])
        Cnt =  [0.] * N
        Area = [0.]
        PI = 3.14159265
        for i in range(1,N+1):
            r = float(R)*i/N
            Area.append(PI*r**2-Area[-1])
        Area = Area[1:]
        #lastRingArea = 0
        X = self.position['X']
        Y = self.position['Y']
        for line in Table:
            r = sqrt((line[0]-X)**2+(line[1]-Y)**2)
            if r >= R: 
                continue
            Cnt[int(r/R*N)]+=1
        CntPerPix = []
        CntPerPixErr = []
        for i in range(N):
            CntPerPix.append(Cnt[i]/Area[i])
            CntPerPixErr.append(sqrt(Cnt[i])/Area[i])
        radius = [float(R)*x/N for x in range(N+1)]
        #print 'Src Cnt, area, Cnt/pixel', Cnt, Area, CntPerPix
        r = []
        for i in range(1,N+1):
            r.append((radius[i]+radius[i-1])/2.)
        if not bg == None:
            bkg = psregion(file=bg)
            BkgTable = filterevents(self.evtfile, cols=['x','y'], region=bg)
            bgCntRate = len(BkgTable)/bkg.area()
            #print 'Bkg Cnt, area, Cnt/pixel',len(BkgTable), bkg.area(), bgCntRate
            CntPerPix = []
            for i in range(N):
                Cnt[i] += -1.*bgCntRate*Area[i]
                #if Cnt[i] < 0.:
                    #Cnt[i] = 0.
                CntPerPix.append(Cnt[i]/Area[i])
        if not outfile == None:
            tofile(outfile,r, CntPerPix, CntPerPixErr)
        return (r, CntPerPix, CntPerPixErr)



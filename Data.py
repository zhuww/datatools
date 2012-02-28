from SaveLoadable import MetaSaveLoader
import os, fnmatch, re
from fitstools import *
import tarfile 
from Exceptions import *
from numpy import mean
#class RawData(SaveLoadable):
class RawData(object):
    """A class that collect basic data infomation, and provide methods for unpacking and inspecting the zipped/tarred data files and saving the collected informations."""
    __metaclass__ = MetaSaveLoader
    TarFile=''
    ObsID=''
    def __init__(self, Location):
        """Way to loading the data using the path to the zipped/tarred file."""
        if os.access(Location, os.R_OK):
            #print 'File %s readable;' % (Location)
            self.TarFile = Location
        else: 
            raise FileError(Location) 
        self.rootpath = os.path.split(Location)[0]
        self.inspect()
    def unzip(self):
        """Function that unzip/untar the zipped/tarred data file."""
        tar = tarfile.open(self.TarFile)
        tar.extractall()
        tar.close()
    def inspect(self):
        """Function that inspect the data file and collect relative informations like target, telescope, instrument, mode, PI, etc."""
        tar = tarfile.open(self.TarFile, 'r')
        ObsID=[]
        for tarinfo in tar:
            FirstBlock = os.path.split(tarinfo.name)[0].split('/')[0]
            try:
                FirstBlock = int(FirstBlock)
                if any([FirstBlock == ID for ID in ObsID]):pass
                else: ObsID.append(FirstBlock)
            except:
                raise PatternError(FirstBlock, 'ObsID')
        tar.close()
        print 'Find Possible Observation IDs: ', ObsID
        self.ObsID = ObsID
 

class chandraData(object):
    __metaclass__ = MetaSaveLoader
    """A class that point to an unpacked chandra Data directory, and automaticly look up the important files such as the event list file."""
    def __init__(self, Path):
        if os.access(Path, os.R_OK):
            #print 'File %s readable;' % (Location)
            self.path = Path
            if Path.split('/')[-1] == '':
                LastPart = Path.split('/')[-2]
            else:
                LastPart = Path.split('/')[-1]
            try:
                LastPart = int(LastPart)
                self.ObsID = LastPart
                if os.path.exists(Path+'/primary/'):pass
                else: FileError('/primary/')
            except:
                raise PatternError('Base name of the path %s' % (LastPart), 'ObsID')
        else: 
            raise FileError(Path)
        def findfile(file, dir, pattern, name):
            if fnmatch.fnmatch(file, pattern): 
                if os.path.splitext(file)[1] == '.gz':
                    os.system('gunzip '+Path+dir+file)
                    self.__dict__[name] = Path+dir+file[:-3]
                else:
                    self.__dict__[name] = Path+dir+file
                return True
            else: return False
        for file in os.listdir(Path+'/primary/'):
            if findfile(file, '/primary/', '*_evt2.fits*', 'evtfile'):print 'event list found.'
            elif findfile(file, '/primary/', '*_asol*', 'asolfile'):print 'asolfile found'
            else:pass
        for file in os.listdir(Path+'/secondary/'):
            if findfile(file, '/secondary/', '*_pbk*', 'pbkfile'):print 'pbk file found.'
            else:pass


class region(object):
    def __init__(self, regionstr):
        p = re.compile('(?P<shape>circle|annulus|ellipse|box)\(\s*(?P<X>\d+\.*\d+)\s*,\s*(?P<Y>\d+\.*\d+)\s*,\s*(?P<R>\d+\.*\d+)(\s*,\s*(?P<Rout>\d+\.*\d+)(\s*,\s*(?P<Angle>\d+\.*\d+)){0,1}){0,1}\)', re.U)
        for m in p.finditer(regionstr):
            self.shape = m.group('shape')
            self.X = float(m.group('X'))
            self.Y = float(m.group('Y'))
            self.R = float(m.group('R'))
            try:
                self.Rout = float(m.group('Rout'))
                self.Angle = float(m.group('Angle'))
            except:pass
    def area(self):
        PI = 3.14159265
        if self.shape == 'circle':
            return PI*self.R**2
        if self.shape == 'annulus':
            return PI*(self.Rout**2 - self.R**2)
        if self.shape == 'box':
            return self.R*self.Rout
        if self.shape == 'ellipse':
            return PI*self.R*self.Rout
    def __str__(self):
        if self.shape == 'ellipse' or self.shape == 'box':
            return '%s(%f,%f,%f,%f,%f)' % (self.shape, self.X, self.Y, self.R, self.Rout, self.Angle)
        elif self.shape == 'annulus':
            return '%s(%f,%f,%f,%f)' % (self.shape, self.X, self.Y, self.R, self.Rout)
        else:
            return '%s(%f,%f,%f)' % (self.shape, self.X, self.Y, self.R)
    def isinside(self, x, y):
        if self.shape == 'circle':
            if (x-self.X)**2 + (y-self.Y)**2 < self.R**2:
                return True 
            else:
                return False
        elif self.shape == 'annulus':
            rsquared = (x-self.X)**2 + (y-self.Y)**2
            if rsquared < self.Rout**2 and rsquared > self.R**2:
                return True 
            else:
                return False

class psregion(object):
    '''A class that provides APIs for accessing the region files of point sources.'''
    def __init__(self, file=None, text=None):
        #if not file[0] == '/':
            #self.file = os.getcwd()+'/'+file
        #else:
            #self.file = file
        if text == None:
            text = open(file,'r').read()
            #text = file
        self.header = ''
        self.regions = []
        regionstr = ''
        for lines in text.split('\n'):
            if lines.find('(') == -1 or lines.find(')') == -1 or not lines.find('#') == -1:
                self.header +=lines
            else: 
                self.regions.append(region(lines))
                #regionstr +=lines
    def area(self):
        return sum([x.area() for x in self.regions])
    def __str__(self):
        return self.header + '\n' + '\n'.join([str(x) for x in self.regions])
    def save(self, file):
        if file == None:file = self.file
        f = open(file, 'w')
        f.write(str(self))
    def filter(self, table, X_index, Y_index):
        newtable = []
        for row in table:
            x = row[X_index]
            y = row[Y_index]
            if any([reg.isinside(x, y) for reg in self.regions]):
                newtable.append(row)
        return newtable

def _tableandindices_(evtfile):
    file = fitsfile(evtfile)
    table = file.gettable()
    cols = file.colinfo()
    return (table, cols.names)


def improve_position(evtfile, region):
    '''Improve the region file set around a point source by the mean position of all the events in the region.'''
    (table, names) = _tableandindices_(evtfile)
    (X_index, Y_index) = names.index('x'), names.index('y')
    print 'Before improving position, the region reads:\n', str(region)
    old_Radius = region.R
    newtable = table
    for region_enlarge_factor in [1.5, 1.]:
        region.R = region_enlarge_factor * old_Radius
        newtable = region.filter(newtable, X_index, Y_index)
        region.X = mean([row[X_index] for row in newtable])
        region.Y = mean([row[Y_index] for row in newtable])
    print 'After improving position, the region reads:\n', str(region)
    print 'But you still have to save this region using the save() method of region object.'

def filterevents(evtfile, Emin=0.3, Emax=10.0, region=None, outfile=None, cols=None, id='energy'):
    '''Filter for event table or creat filtered file.'''
    if not outfile == None: #outfile='src%s-%s.fits' % (Emin, Emax)
        if region == None: regionfmt = ''
        else:
            regionfmt = '[sky=region(%s)]' % (region)
        if cols == None: columnfmt = ''
        else:
            columnfmt = '[col %s]' % ','.join(cols)
        if  os.access(outfile,os.R_OK) and outfile.find('*') == -1: os.system('rm %s' % outfile)
        cmdline = 'dmcopy "%(evtfile)s[%(id)s=%(Emin).0f:%(Emax).0f]%(region)s%(column)s" %(outfile)s clobber=yes ' %  {
                'evtfile':evtfile, 
                'id':id,
                'Emin':Emin*1000, 
                'Emax':Emax*1000, 
                'region':regionfmt,
                'column':columnfmt,
                'outfile':outfile}
        print cmdline
        os.system(cmdline)
        return fitsfile(outfile).gettable()
    else:
        (table, names) = _tableandindices_(evtfile)
        try:
            (X_index, Y_index, E_index) = [names.index(item) for item in ['x','y','energy']]
        except:
            (X_index, Y_index, E_index) = [names.index(item) for item in ['X','Y','ENERGY']]
        if not region == None: 
            text = open(region, 'r').read()
            psreg = psregion(text=text)
            table = psreg.filter(table, X_index, Y_index)
        newtable = []
        for row in table:
            if row[E_index] > Emin*1000 and row[E_index] < Emax*1000:
                newtable.append(row)
        del(table)
        table = newtable
        if cols == None: return table
        else:
            indics = [names.index(item) for item in cols]
            newtable= []
            for row in table:
                newtable.append([row[i] for i in indics])
        del(table)
        table = newtable
        return table

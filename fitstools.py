import pyfits, os, re
import numpy as np
#from numpy import mean
#from numpy import mean
#import commands
from datatools import Exceptions

class fitsfile(object):
    '''A group of APIs for accessing the FITS files'''
    def __init__(self, file=None):
        '''Open up a FITS file for editing and inspection'''
        if not isinstance(file, basestring):
            raise TypeError, 'Must specify a name for the fits file.'
        elif os.access(file,os.R_OK):
            self.file = pyfits.open(file)
        else:raise  Exceptions.FileError(file)

    def hdread(self, key, block=0):
        '''Read from the header of some block of the FITS file'''
        result = []
        if block == 0:
            blocklist = range(len(self.file))
        else:
            blocklist = [block]
        for block in blocklist:
            try :
                if isinstance(key, basestring):
                    if key == '' or key =='all' or key =='*':
                        result.append(self.file[block].header)
                    else:
                        result.append(self.file[block].header[key])

                elif isinstance(key, (list, tuple)):
                    result.extend([ self.file[block].header[everykey] for everykey in key ])
            except:continue
        return result[0]

    def hdwrite(self, key, value, block=0):
        '''Update the header of some block of the FITS file'''
        self.file[block].header[key] = value
        return self

    def colinfo(self, block=1):
        '''Display informations about the columns of the FITS file'''
        cols = self.file[block].columns
        #cols.info()
        return cols

    def gettable(self,block=1,cols=None):
        '''Return the table of some block of the FITS file'''
        table = self.file[block].data
        if cols == None: 
            return table
        else:
            #names = self.colinfo(block).names
            #indics = [names.index(item) for item in cols]
            #newtable= []
            #for row in table:
                #newtable.append([row[i] for i in indics])
        #del(table)
        #table = newtable
            newtable = np.vstack([table.field(col) for col in cols]).T
        return newtable

    def save(self):
        '''Save the updated FITS file'''
        self.file.close()



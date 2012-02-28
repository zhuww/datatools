class TypeError(Exception):
    def __init__(self, Data):
        print 'Unrecongnizable Data %s' % Data
class FileError(Exception): 
    def __init__(self, file):
        print '%s could not be found.' % file
class PatternError(Exception):
    def __init__(self, target, pattern):
        print '%s does not match the pattern %s.' % (target, pattern)

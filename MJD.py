import datetime

MJDzero = datetime.datetime.utcfromtimestamp(0)+datetime.timedelta(seconds=-float(0x007c95674beb4000)/10000000)

def datetime_to_MJD(d):
    w = d-MJDzero
    return w.days + (w.seconds + w.microseconds/1e6)/86400.

def MJD_to_datetime(m):
    return MJDzero + datetime.timedelta(days=m)

def MJD_to_year(m):
    dayinyear = 365.242199
    t = MJD_to_datetime(m)
    Y= datetime.datetime(t.year, 1, 1)
    f = (t - Y).days/dayinyear
    return t.year + f

if __name__=='__main__':
    print datetime_to_MJD(datetime.datetime.utcnow())

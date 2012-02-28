from decimal import *
from round import *
def shorten(value, error):
    fig = figure((float(str(value)), float(str(error))))
    valuestr = str(fig)
    first, last = valuestr.split(')')
    head, toe = first.split('(')
    result = head+last
    return Decimal(result)

print shorten(Decimal('1.23456e-12'), Decimal('2.34567e-24'))
#q = Decimal(10)** -2
#print Decimal('1.23456e-10').quantize(q).as_tuple()

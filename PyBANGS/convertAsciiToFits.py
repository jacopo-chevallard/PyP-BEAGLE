import numpy as np

import pyfits

fileName = "/Users/efcl/BANGS/data/ECL13/ECL13.cat"
data = np.genfromtxt(fileName, dtype=None, names=True)
outputFile = "/Users/efcl/BANGS/data/ECL13/ECL13.fits" 

print data.dtype

columnNames = data.dtype.names

if columnNames[0] != 'ID' and columnNames[0] != 'id':
    format = 'E'
else:
    format = 'A'
    
print columnNames[0], format

col = pyfits.Column(name=columnNames[0], format=format, array=data[columnNames[0]])
allColumns=pyfits.ColDefs([col])

for i in range(1,len(columnNames)):
    column = columnNames[i]
    print column
    if column != 'ID' and column != 'id':
        format = 'E'
    else:
        format = 'A'
    col = pyfits.Column(name=column, format=format, array=data[column])
    allColumns.add_col(col)

tbhdu = pyfits.BinTableHDU.from_columns(allColumns)
tbhdu.data = data
tbhdu.writeto(outputFile)

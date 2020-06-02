from __future__ import absolute_import
from datetime import datetime
import numpy as np
from numpy.lib.recfunctions import drop_fields, append_fields

class NLDNdataFile(object):
    
    stroke_DC3 = {'columns':[ ('date','S10'), ('time','S20'), 
                              ('lat','f4'), ('lon','f4'), 
                              ('peak_current','f4'), ('ellipse','f4'),
                              ], 
                  'date_dtype':[('year','i2'),('month','i1'),('day','i1')],
                  'time_dtype':[('hour','i1'),('minute','i1'),('second','float64')]
                 }

    stroke_ICCG = {'columns':[ ('date','S10'), ('time','S20'),
                              ('lat','f4'), ('lon','f4'),
                              ('peak_current','f4'), ('ICCG','S1'),
                              ],
                  'date_dtype':[('year','i2'),('month','i1'),('day','i1')],
                  'time_dtype':[('hour','i1'),('minute','i1'),('second','float64')]
                 }
    
    def __init__(self, filename, date_sep='-', time_sep=':', format='stroke_DC3'):
        """ Load NLDN data from a file, into a numpy named array stored in the
            *data* attribute. *data*['time'] is relative to the *basedate* datetime
            attribute
            """
        self.format=format
        
        dtype_specs = getattr(self, format)
        
        
        nldn_initial = np.genfromtxt(filename, dtype=dtype_specs['columns'])
        date_part = np.genfromtxt(nldn_initial['date'],
                        delimiter=date_sep, dtype=dtype_specs['date_dtype'])
        time_part = np.genfromtxt(nldn_initial['time'],
                        delimiter=time_sep, dtype=dtype_specs['time_dtype'])
        dates = [datetime(a['year'], a['month'], a['day'], b['hour'], b['minute']) 
                    for a, b in zip(date_part, time_part)]
        min_date = min(dates)
        min_date = datetime(min_date.year, min_date.month, min_date.day)
        t = np.fromiter( ((d-min_date).total_seconds() for d in dates), dtype='float64')
        t += time_part['second']
        
        self.basedate = min_date
        data = drop_fields(nldn_initial, ('date', 'time'))
        data = append_fields(data, 'time', t)
        
        self.data = data
        
from __future__ import absolute_import
from __future__ import print_function
import numpy as np
from datetime import datetime, timedelta

def to_seconds(dt):
    'convert datetime.timedelta object to float seconds'
    return dt.days * 86400 + dt.seconds + dt.microseconds/1.0e6

def read_flashes(h5LMAfiles, target, base_date=None, min_points=10):
    """ This routine is the data pipeline source, responsible for pushing out 
        events and flashes. Using a different flash data source format is a matter of
        replacing this routine to read in the necessary event and flash data."""
    import tables

    for filename in h5LMAfiles:
        h5 = tables.openFile(filename)
        table_name = list(h5.root.events._v_children.keys())[0]

        # event_dtype = getattr(h5.root.events, table_name)[0].dtype
        # events_all = getattr(h5.root.events, table_name)[:]
        flashes = getattr(h5.root.flashes, table_name)
        events = getattr(h5.root.events, table_name)[:]
        
        if flashes.shape[0] == 0:
            # There are no flashes
            pass
        
        fl_dtype = flashes.dtype
        flashes = np.fromiter((fl[:] for fl in flashes if fl['n_points'] >= min_points), dtype=fl_dtype)

        # get the start date of the file, which is the seconds-of-day reference
        h5_start =  datetime(*getattr(h5.root.events, table_name).attrs['start_time'][0:3])
        if base_date==None:
            base_date = h5_start
        extra_dt = to_seconds(h5_start - base_date)
        # adjust time in seconds of day to be in correct reference
        events['time'] += extra_dt
        flashes['start'] += extra_dt

        n_flashes = len(flashes)
        print('{0} -- {1} flashes > {2} pts; dt+={3}  '.format(filename, n_flashes, min_points, extra_dt))
        
        push_out = (events, flashes)

        target.send(push_out)
        
        # del events
        # del flashes
        del push_out

        h5.close()
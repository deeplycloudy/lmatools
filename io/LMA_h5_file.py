from __future__ import absolute_import
from __future__ import print_function
import numpy as np
from datetime import datetime, timedelta
import tables

def to_seconds(dt):
    'convert datetime.timedelta object to float seconds'
    return dt.days * 86400 + dt.seconds + dt.microseconds/1.0e6
    


class LMAh5Collection(object):
    def __init__(self, filenames, base_date=None, min_points=1):
        """ Load multiple HDF5 files given by filenames. 
            
            h5s = LMAh5Collection(filenames)
        
            >>> for events, flashes in h5s:
            >>>     print(events.shape, flashes.shape)

            Or, if you know a time accurately, you can do:
            >>> from datetime import datetime
            >>> t = datetime(2013,6,6,3,0,0)
            >>> events, flashes = h5s.data_for_time(t)
            Another way of obtaining a valid t is to loop over h5s.times
        """
        self._filenames = filenames
        self.min_points = min_points
        # populate _time_lookup by calling _all_times
        self._time_lookup = {} 
        self.times = [t for t in self._all_times()]
        print("collection times = ", self.times)
        self.times.sort()
        
        if base_date is None:
            t_min = min(self.times)
            y, m, d = t_min.year, t_min.month, t_min.day
            self.base_date = datetime(y,m,d)
        else:
            self.base_date = base_date
    
    def _table_times_for_file(self, fname):
        """ Called once by init to set up frame lookup tables and yield 
            the frame start times. _time_lookup goes from 
            datetime->(h5 filename, table_name).
        """
        h5 = LMAh5File(fname, min_points=self.min_points)
        for t_start, table_name in zip(h5.start_times, h5.table_names):
            self._time_lookup[t_start] = (fname, table_name)
        yield t_start

            
    def _all_times(self):
        for f in self._filenames:
            for t in self._table_times_for_file(f):
                yield t 
                
    def data_for_time(self, t0):
        """ Return events, flashes whose start time matches datetime t0.
        
        events['time'] and flashes['start'] are corrected if necessary
        to be referenced with respect to self.base_date.
        """
        fname, table_name = self._time_lookup[t0]
        h5 = LMAh5File(fname, min_points=self.min_points)
        events, flashes = h5.data_for_table(table_name)
        
        extra_dt = to_seconds(h5.base_date - self.base_date)
        print("data for time extra dt = {0}", extra_dt)
        events['time'] += extra_dt
        flashes['start'] += extra_dt
        
        return events, flashes
        
    def data_nearest_time(self, t0):
        """ Return events, flashes whose start time is just before t0
        
                data start time
                t1     t2    t3            t5
                |      |     |             |
                                  |
                                  t0 
        tn - t0 = -    -     -             +
        
        Want to return data corresponding to time t3
        """
        pass
        
    def __iter__(self):
        for t in self.times:
            events, flashes = self.data_for_time(t)
            yield events, flashes
        
        

class LMAh5File(object):
    def __init__(self, filename, min_points=1):
        """ Open an  HDF5 LMA file so as to read LMA event/flash tables.
            The events and flashes groups may contain more than one table,
            and the list of self.table_names corresponds to self.start_times 
            and are sorted in increasing time order.
        
            self.base_date provides the base_date against which 
            time in seconds is referenced by all other methods.
        
            CAUTION
            This class is thought to support more than one table, but has not
            been tested with HDF5 LMA data files with more than one table.
        """
        self.h5 = tables.open_file(filename)
        self.min_points = min_points
        table_names = list(self.h5.root.events._v_children.keys())
        start_times = [datetime(*getattr(self.h5.root.events, 
                                    the_table).attrs['start_time']) for
                                    the_table in table_names] 
        
        # Sort by start time.                
        st, tn  = zip(*sorted(zip(start_times, table_names)))
        self.start_times, self.table_names = list(st), list(tn)
                       
        # Date of each table             
        self._base_dates = [datetime(start_time.year,
                                     start_time.month,
                                     start_time.day) for start_time
                                     in self.start_times]
        self.base_date = min(self._base_dates)
        self.table_metadata = {}
        for name, t0, date in zip(self.table_names,
                                  self.start_times,
                                  self._base_dates):
            self.table_metadata[name] = {'start_time': t0,
                                         'start_date': date}
                                         
        print(self.start_times, self.table_names, self.table_metadata)
        
    def data_for_table(self, table_name):
        """ Return events, flashes as numpy arrays with a named dtype. 
            Returns all events, but only those flashes with n_points 
            >= min_points. 
        
            events['time'] and flashes['start'] are corrected if necessary
            to be referenced with respect to self.base_date.
        """
        table_date = self.table_metadata[table_name]['start_date']
        extra_dt = to_seconds(table_date - self.base_date)
        print("data for table dt = {0}", extra_dt)
        
        flashes = getattr(self.h5.root.flashes, table_name)
        fl_dtype = flashes.dtype
        flashes = np.fromiter(
                  (fl[:] for fl in flashes if fl['n_points'] >= self.min_points), 
                  dtype=fl_dtype
                  )
        events = getattr(self.h5.root.events, table_name)[:]
        
        events['time'] += extra_dt
        flashes['start'] += extra_dt
        
        return events, flashes

    def gen_events_flashes(self):
        """ Generate events, flashes as numpy arrays with a named dtype. 
            Returns all events, but only those flashes with n_points 
            >= min_points. 
        
            events['time'] and flashes['start'] are corrected if necessary
            to be referenced with respect to self.base_date.
        """
        for table_name in self.table_names:
            events, flashes = data_for_table(self, table_name)            
            yield events, flashes
        

def read_flashes(h5LMAfiles, target, base_date=None, min_points=10):
    """ This routine is the data pipeline source, responsible for pushing out 
        events and flashes. Using a different flash data source format is a matter of
        replacing this routine to read in the necessary event and flash data."""

    # for filename in h5LMAfiles:
        # h5 = LMAh5File(filename, min_points=min_points)
        # for events, flashes in h5.gen_events_flashes():
        #
        #     # get the start date of the file, which is the seconds-of-day reference
        #     h5_start =  h5.base_date
        #     if base_date==None:
        #         base_date = h5_start
        #
        #     # adjust time in seconds of day to be in correct reference
        #     extra_dt = to_seconds(h5_start - base_date)
        #     events['time'] += extra_dt
        #     flashes['start'] += extra_dt
        #
        #     n_flashes = len(flashes)
        #     print('{0} -- {1} flashes > {2} pts; dt+={3}  '.format(filename, n_flashes, min_points, extra_dt))
    
    h5s = LMAh5Collection(h5LMAfiles, base_date=base_date, min_points=min_points)    
    for push_out in h5s: # get a stream of events, flashes 
        target.send(push_out)
        del push_out
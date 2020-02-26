import datetime
from collections import namedtuple

from lmatools.io.LMA_h5_write import write_h5
from lmatools.io.LMAarrayFile import LMAdataFile

class LMADataset(object):    
    def __init__(self, filename=None,  
                 data=None, basedate=None, startdate=None,
                 sec_analyzed=None, header=''):
        """ Create a new LMA Dataset which can be used in a standardized way 
            by the clustering routines and written to a standardized HDF5 
            output format.
                 
            The data attribute may be overwritten by flash sorting, as
                 necessary to maintain a self-consistent dataset.
        
            Data may be loaded in two ways.
        
            1. From an LMA ASCII data file. Use the filename kwarg.
               If the file has a hexadecimal station mask length
               that is different than six characters, use the file_mask_length
               keyword argument.
        
            2. Provide data directly using the data, basedate, startdate
               and sec_analyzed kwargs. data is an array with a named dtype of 
               (minimally, see below)
               the_dtype = [('time', '<f8'), ('latitude', '<f8'), 
                           ('longitude', '<f8'), ('altitude', '<f8'), 
                           ('chi2', '<f8'), ('stations', '<f8')] 
               This may also work with a pandas dataframe with the column names
               as above, but this has not been tested.
               basedate is as defined below.
               sec_analyzed is as defined below.
        
            It is very easy to construct the data array from individual 
                 arrays containing columns of data.
            data = np.empty((N_points,) dtype=the_dtype)
            data['latitiude'] = my_lats
            data['longitude'] = my_lons 
            # etc.
            
            Attributes: 
            basedate: datetime object against which the time coordinate in data
                 is measured. UTC year, month, day.
            startdate: the actual start date of the LMA data, including hour
                 minute, and second. For a collection of LMAData objects, a list
                 of startdates should correspond to the starts of the chunks of 
                 data contained in the LMAdata objects.
            sec_analyzed: the total number of seconds covered by data.
                 May differ from t.max() - t.min(). For example, the standard
                 LMA processing produces ten minute files, which may only contain
                 a flash or two.
            header: header of the LMA data file if so initialized, otherwise
                 the value of the kwarg.
            metadata: dynamically generated from the above attributes. See the
                 docstring for the metadata property for further details.
            data: data table which must minimally contain the following fields
                 in order for flash sorting to work.
                 The data table should use the following names and units:
                     lon (degrees)
                     lat (degrees)
                     alt (meters; above sea level or ellipsoid)
                     time (UTC seconds since start of day given by basedate)
                     chi2 (chi squared quality parameter)
                     stations (number of contributing stations)
                 Optional names and units: 
                     power (received power, dBW)
                     mask (station bitmask, in hexadecimal)
                     charge_id (signed integer, typically <= 1 byte)
                     flash_id (integer, may be added by flash sorting)
                     point_id (integer, unique ID of each LMA source)
             flashes: list of flash objects created during flash sorting, if any

             
        """
        
        if filename is not None:
            # self.load_data_from_LMA_file(filename, mask_length=file_mask_length)
            self.load_data_from_LMA_file(filename)
        elif data is not None:
            if basedate is None:
                e = "Please provide basedate when manually loading data"
                raise AttributeError(e)
            if startdate is None:
                e = "Please provide startdate when manually loading data"
                raise AttributeError(e)
            if sec_analyzed is None:
                e = "Please provide sec_analyzed when manually loading data"
                raise AttributeError(e)
                
            self.data = data
            self.basedate = basedate
            self.sec_analyzed = sec_analyzed
            self.header = header
        self.flashes = []
        print ('Station mask Length =', self.mask_length)

    @property
    def startyear(self):
        return self.startdate.year
    @property
    def startmonth(self):
        return self.startdate.month
    @property
    def startday(self):
        return self.startdate.day
    @property
    def starthour(self):
        return self.startdate.hour
    @property
    def startminute(self):
        return self.startdate.minute
    @property
    def startsecond(self):
        return self.startdate.second
        
    @property
    def metadata(self):
        """ Return file-level LMA metadata (or its equivalent for manually 
            supplied data) in a form representable by basic numeric and string
            types in, for example, HDF5 or NetCDF formats.
        
            Returns: a named tuple of metadata, with the following values
            accesible as attributes by this method:
                header', 'sec_analyzed',
                'startyear', 'startmonth', 'startday', 
                'starthour', 'startminute', 'startsecond'
        """
        metadata_fields = ['header', 'sec_analyzed',
            'startyear', 'startmonth', 'startday', 
            'starthour', 'startminute', 'startsecond']
        Metadata = namedtuple('Metadata', metadata_fields)

        d = {}
        for k in metadata_fields:
            d[k] = getattr(self, k)
        
        m = Metadata(**d)
        return m

            
    def load_data_from_LMA_file(self, filename):
        """ Load data from a single LMA file. Filter data using
            minimum number of stations and maximum chi2 values provided
            in the flashsort params dictionary.
        
            The data are stored in an 
        """
            
        lma=LMAdataFile(filename)
        # lma=LMAdataFile(filename, mask_length = mask_length)
        # make sure the number of stations are calculated from the mask
        # see the implementation of LMAdataFile.__get_attr__
        stns = lma.stations
        self.header = ''.join(lma.header)
        self.basedate = datetime.datetime(lma.startyear, 
                            lma.startmonth, lma.startday)
        self.startdate = datetime.datetime(lma.startyear, 
                            lma.startmonth, lma.startday, 
                            lma.starthour, lma.startminute, 
                            lma.startsecond)
        self.sec_analyzed = lma.sec_analyzed
        self.data = lma.data
        self.mask_length = lma.mask_length
        # return (lma, data)
    
    def filter_data(self, params):
        
        """ params is a dictionary with tuple values that give the lower
            and upper bounds. Assumes stations and chi2 are always provided.
            Also checks the altitude range if provided. 
        
            This prevents someone from working with a ten station LMA setting a
            range of (5,10), and then switching to a sixteen station LMA and
            silently losing solutions retrieved by more than ten stations. A
            similar logic applies to chi2.
        
            TODO: use params to automatically apply a filter
                  if the params key matches an available data type 
        
            Exceptions:
            Only the upper bound is used on chi2
            Only the lower bound is used on stations
        """
        # Filter out noisy or otherwise mislocated points
        good = ((self.data['stations'] >= params['stations'][0]) & 
                (self.data['chi2'] <= params['chi2'][1]) 
                )
        if 'alt' in params:
            good = good & ((self.data['altitude'] < params['altitude'][1]) & 
                           (self.data['altitude'] > params['altitude'][0])
                    )
        
        return self.data[good]
        
    def write_h5_output(self, outfile, orig_LMA_file):
        """ Thin wrapper around common h5 writer; makes use of metadata
            and flashes already stored as attributes.
        """
        write_h5(outfile, self.flashes, self.metadata, orig_LMA_file, self.mask_length)

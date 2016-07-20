from __future__ import absolute_import
from __future__ import print_function
import os, sys, re
import glob
import tempfile
import shutil
import subprocess
import logging, logging.handlers
import datetime
from collections import namedtuple

from .autosort.write_flashes import write_h5
from .autosort.LMAarrayFile import LMAdataFile
# from .autosort.flash_stats import FlashMetadata

LOG_BASEFILENAME = datetime.datetime.now().strftime('Flash-autosort.log')

def logger_setup(logpath):
    logger = logging.getLogger('FlashAutorunLogger')
    logfile = os.path.join(logpath, LOG_BASEFILENAME)
    loghandler = logging.handlers.RotatingFileHandler(logfile, maxBytes=1024*1024, backupCount=3)
    logger.addHandler(loghandler)
    logger.setLevel(logging.DEBUG)


def write_output(outfile, flashes, orig_LMA_file, metadata=None):
    if metadata is None:
        # use metadata from the first flash as the canonical metadata,
        #   since all flashes were sorted fromt the same LYLOUT file
        # This breaks in the case of empty flashes; in that case the calling function should pass in metadata.
        metadata = flashes[0].metadata
    write_h5(outfile, flashes, metadata, orig_LMA_file)

    
def sort_files(files, output_path, clusterer=None):
    
    logger = logging.getLogger('FlashAutorunLogger')
    now = datetime.datetime.now().strftime('Flash autosort started %Y%m%d-%H%M%S')
    logger.info(now)
    
    h5_outfiles = []
    for a_file in files:
        try:
            file_base_name = os.path.split(a_file)[-1].replace('.gz', '')
            outfile = os.path.join(output_path, file_base_name+'.flash')
            
            # clusterer should use the name outfile as the base for any, e.g., ASCII data it would like to save
            
            lmadata = LMADataset(a_file)
            clusterer(lmadata)
            
            # header = ''.join(lmadata.header)
            # fl_metadata = FlashMetadata(header)
            outfile_with_extension = outfile + '.h5'
            h5_outfiles.append(outfile_with_extension)
            
            # make write_output a method of the dataset, and take dataset instead of a_file.
            write_output(outfile_with_extension, lmadata.flashes, 
                a_file, metadata=lmadata.metadata)
        
        except:
            logger.error("Did not successfully sort %s \n Error was: %s" % (a_file, sys.exc_info()[1]))
            raise
    # loghandler.doRollover()
    return h5_outfiles


class LMADataset(object):
    """ Things to accomodate:
    
        Events data table
        Flashes data table
    
        Read LMA ASCII data (delegate)
        Write H5 file (delegate)
    """
    
    def __init__(self, filename=None, file_mask_length=6, 
                 data=None, basedate=None, sec_analyzed=None, header=''):
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
        
            2. Provide data directly using the data, basedate, and sec_analyzed
               kwargs. data is an array with a named dtype of 
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
            self.load_data_from_LMA_file(filename, mask_length=file_mask_length)
        elif data is not None:
            if basedate is None:
                e = "Please provide basedate when manually loading data"
                raise AttributeError(e)
            if sec_analyzed is None:
                e = "Please provide sec_analyzed when manually loading data"
                raise AttributeError(e)
                
            self.data = data
            self.basedate = basedate
            self.sec_analyzed = sec_analyzed
            self.header = header
        self.flashes = []
        
    @property
    def startyear(self):
        return self.basedate.year
    @property
    def startmonth(self):
        return self.basedate.month
    @property
    def startday(self):
        return self.basedate.day
    @property
    def starthour(self):
        return self.basedate.hour
    @property
    def startminute(self):
        return self.basedate.minute
    @property
    def startsecond(self):
        return self.basedate.second
        
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

            
    def load_data_from_LMA_file(self, filename, mask_length=6):
        """ Load data from a single LMA file. Filter data using
            minimum number of stations and maximum chi2 values provided
            in the flashsort params dictionary.
        
            The data are stored in an 
        """
            
        lma=LMAdataFile(filename, mask_length = mask_length)
        # make sure the number of stations are calculated from the mask
        # see the implementation of LMAdataFile.__get_attr__
        stns = lma.stations
        self.header = ''.join(lma.header)
        self.basedate = datetime.datetime(lma.startyear, 
                            lma.startmonth, lma.startday)
        self.sec_analyzed = lma.sec_analyzed
        self.data = lma.data
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

""" Next things to untangle
        
    Clustering code takes a file. Design API for providing data in some standardized way.
    That way the clustering code can be used for data from files or from some other source.
    
    Need to decouple loading (and maybe filtering, but QC might belong with clustering) 
    data, maybe to a separate object that just deals with LMA file IO 
    A separate one for writing out results.
    
    The API for providing data also needs to account for the odd use of the lma object
    to store flash results. Suggest starting there and figuring out a good way to
    report results that isn't so opaque, which will also clarify how closely coupled it is
    to the LMADataFile object.
    
    Might be the key to a standardized hierarchical data model.
    
    Really cluster.create_flash_objs should be part of the Dataset object
        since by OO principles the function that modifies a dataset should
        be grouped with that dataset?
        
    Write instructions about how to convert an old run_files_with_params to
        the new gen_flashes script.
    Update lmaworkshop materials by following the above instructions 
    
    Also TODO:
    1. figure out how to use this method with old run_files_with_parmas 
        without rewriting a bunch of driver scripts
    2. fully decouple from things containing autosort in the input lines.
        Maybe just pull everything out of that directory.
    3. Get rid of old use of FlashMetadata class in mflash and autorun_sklearn, 
        in favor of new style in gen_sklearn. It was never flash metadata anyway,
        and belonged in an LMA dataset metadata class.
    """
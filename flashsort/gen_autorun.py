from __future__ import absolute_import
from __future__ import print_function
import os, sys, re
import glob
import tempfile
import shutil
import subprocess
import logging, logging.handlers
import datetime

from .autosort.write_flashes import write_h5
from .autosort.flash_stats import FlashMetadata

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
            
            lmadata, flashes = clusterer(a_file)
            
            header = ''.join(lmadata.header)
            fl_metadata = FlashMetadata(header)
            outfile_with_extension = outfile + '.h5'
            h5_outfiles.append(outfile_with_extension)
            write_output(outfile_with_extension, flashes, a_file, metadata=fl_metadata)
        
        except:
            logger.error("Did not successfully sort %s \n Error was: %s" % (a_file, sys.exc_info()[1]))
            raise
    # loghandler.doRollover()
    return h5_outfiles


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
    
    
    Also TODO:
    1. figure out how to use this method with old run_files_with_parmas 
        without rewriting a bunch of driver scripts
    2. fully decouple from things containing autosort in the input lines.
        Maybe just pull everything out of that directory.
    """
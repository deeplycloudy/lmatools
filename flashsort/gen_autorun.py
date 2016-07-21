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

from lmatools.io.LMA import LMADataset

LOG_BASEFILENAME = datetime.datetime.now().strftime('Flash-autosort.log')

def logger_setup(logpath):
    logger = logging.getLogger('FlashAutorunLogger')
    logfile = os.path.join(logpath, LOG_BASEFILENAME)
    loghandler = logging.handlers.RotatingFileHandler(logfile, maxBytes=1024*1024, backupCount=3)
    logger.addHandler(loghandler)
    logger.setLevel(logging.DEBUG)

    
def sort_files(files, output_path, clusterer):
    
    logger = logging.getLogger('FlashAutorunLogger')
    now = datetime.datetime.now().strftime('Flash autosort started %Y%m%d-%H%M%S')
    logger.info(now)
    
    h5_outfiles = []
    for a_file in files:
        try:
            file_base_name = os.path.split(a_file)[-1].replace('.gz', '')
            outfile = os.path.join(output_path, file_base_name+'.flash')
            
            lmadata = LMADataset(a_file)
            clusterer(lmadata)
            
            outfile_with_extension = outfile + '.h5'
            h5_outfiles.append(outfile_with_extension)
            
            lmadata.write_h5_output(outfile_with_extension, a_file)
        
        except:
            logger.error("Did not successfully sort %s \n Error was: %s" % (a_file, sys.exc_info()[1]))
            raise
    # loghandler.doRollover()
    return h5_outfiles


"""         
    Write instructions about how to convert an old run_files_with_params to
        the new gen_flashes script.
    Update lmaworkshop materials by following the above instructions 
    
    
    Would like to better reflect the three steps of the normal workflow
        in the structure of the project
        - sort flashes
        - calculate flash statistics and
          produce derived products including gridded data
        - make plots
    Reorganize IO stuff outside flashsort directory, maybe including
        LMADataset model, etc. Also get rid of unnecessary 
        hierarchy in directory structure overall.
        
    Try moving mflash into new infrastructure
        
    Redo make_grids to add 3D stuff
        
    Add flash stats and time series code 
    
    Also TODO:
    1. figure out how to use this method with old run_files_with_parmas 
        without rewriting a bunch of driver scripts
    2. fully decouple from things containing autosort in the input lines.
        Maybe just pull everything out of that directory.
    3. Get rid of old use of FlashMetadata class in mflash and autorun_sklearn, 
        in favor of new style in gen_sklearn. It was never flash metadata anyway,
        and belonged in an LMA dataset metadata class.
    """
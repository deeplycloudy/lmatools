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



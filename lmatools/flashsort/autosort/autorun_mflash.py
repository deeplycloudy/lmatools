from __future__ import absolute_import
from __future__ import print_function
import os, sys, re
import glob
import tempfile
import shutil
import subprocess
import logging, logging.handlers
import datetime

from lmatools.flashsort.flash_stats import calculate_flash_stats, Flash, FlashMetadata
from lmatools.io.LMAarrayFile import cat_LMA, LMAdataFile

from .mflash import write_header

LOG_BASEFILENAME = datetime.datetime.now().strftime('Flash-autosort.log')


def logger_setup(logpath):
    logger = logging.getLogger('FlashAutorunLogger')
    logfile = os.path.join(logpath, LOG_BASEFILENAME)
    loghandler = logging.handlers.RotatingFileHandler(logfile, maxBytes=1024*1024, backupCount=3)
    logger.addHandler(loghandler)
    logger.setLevel(logging.DEBUG)


# assume that source code for the flash program is in the same directory as this module
try:
    src_dir = os.path.dirname(__file__)
except:
    src_dir = os.getcwd()
    
build_files = ( #'mflashcustom.f',
                'makefile',
                )
                
flash_output_dir = 'output'     # created within the build directory
flash_prg_name   = 'mflashcustom'   # created within the build directory

tmp_flashsort_prepend = 'flashsort_'

def build(directory=None, mflash_params=None):
    """ Build the flash program in directory; a temporary directory is created by default. """
    

    if directory == None:
        d = tempfile.mkdtemp('', tmp_flashsort_prepend)
    else:
        d = directory
        
    write_header(os.path.join(d,'mflashcustom.f'), **mflash_params)

    for filename in build_files:
        shutil.copyfile(os.path.join(src_dir, filename),
                        os.path.join(d, filename))
    
    os.chdir(d)
    os.mkdir(flash_output_dir)
    
    make_cmd = ['make','-f','makefile']
    rc = subprocess.call(make_cmd)
    
    return d
    

def cleanup_build(directory):
    logger = logging.getLogger('FlashAutorunLogger')
    
    if tmp_flashsort_prepend in directory:
        shutil.rmtree(directory)
    else:
        logger.warning( 
            "Temporary build directory %s doesn't appear to be a flashsort directory." % (directory,)
            )
    

    

def sort_file(filename, directory):
    """ Sort one LMA data file into flashes. dir is the directory with the flash program"""
    logger = logging.getLogger('FlashAutorunLogger')
    
    f, command, the_input = cat_LMA(filename)
    
    run_cmd = [os.path.join(directory, flash_prg_name)]
    logger.info( 'Running %s' % (run_cmd,)) #, 'with stdin from ', command 
    
    # comment out stdout=subprocess.PIPE to print stdout to the terminal. when uncommented,
    #   stdout is captured to python, which leads to less noise in the terminal
    p = subprocess.Popen(run_cmd, stdin=f.stdout, stdout=subprocess.PIPE)#, preexec_fn=f.stdin.close)
    
    # The communication step is key to not blocking at completion.
    out, err = p.communicate()#input=the_input) #out, err not connected to pipes, so nothing to capture or print
    
    # print out
    # print 'Errors: ', err
    return out, err
    
    

def collect_output(datafile, min_points=1, mask_length=4):
    """ collect all output from the flash program output directory created by
        the flash program in flash_dir and calculate some additional stats on each flash        
        """
    import numpy as np
    
    logger = logging.getLogger('FlashAutorunLogger')
    
    # outdir = os.path.join(flash_dir,flash_output_dir)
    # os.chdir(outdir)
    
    lma = LMAdataFile(datafile, mask_length=mask_length)
    
    # get the mapping from flash_ids to the points
    order = np.argsort(lma.flash_id)
    
    # In the case of no data in the file, lma.data.shape will have length zero, i.e., a 0-d array
    if len(lma.data.shape) == 0:
        # No data
        flashes = []
    else:
        flid = lma.flash_id[order]
        boundaries, = np.where(flid[1:]-flid[:-1])    # where indices are nonzero
        boundaries = np.hstack(([0], boundaries+1))
    
        all_data = lma.data[order]
    
        max_idx = len(flid) #- 1
        slice_lower_edges = tuple(boundaries)
        slice_upper_edges = slice_lower_edges[1:] + (max_idx,)
        slices = list(zip(slice_lower_edges, slice_upper_edges))
        
        flashes = [ Flash(all_data[slice(*sl)]) for sl in slices ]
    
        # calculate extra flash metadata, e.g., initation, centroid
        logtext = "Calculating flash initation, centroid, area, etc. for %d flashes" % (len(flashes), )
        logger.info(logtext)
        # print flashes[0].points.dtype
        for fl in flashes:
            header = ''.join(lma.header)
            fl.metadata = FlashMetadata(header)
            calculate_flash_stats(fl)
                    
    return lma, flashes



def cluster(a_file, output_path, outfile, params, logger, min_points=1, retain_ascii_output=True, cleanup_tmp=True):
    # logger = logging.getLogger('FlashAutorunLogger')
    
    now = datetime.datetime.now().strftime('mflash run started %Y%m%d-%H%M%S')
    logger.info(now)
        
    
    d = build(mflash_params=params) # returns temporary directory name
    logger.debug('Built flash program in %s' % d)

    logger.info('Sorting flashes for %s' % a_file)
    out, err = sort_file(a_file, d) # out contains sorted flash data (via stdout)
    # print out


    # file_base_name = os.path.split(a_file)[-1].replace('.gz', '')
    # print 'filebase =  ', file_base_name
    # outfile = os.path.join(output_path, file_base_name+'.flash.h5')
    # write_output(outfile, flashes, file_base_name)
    # print 'Wrote original data and flashes to ', outfile

    # outfile = os.path.join(output_path, file_base_name+'.flash')
    shutil.copyfile(os.path.join(d, 'flashes_out.dat'),
                    outfile)            
    
    if 'mask_length' in params:
        mask_length = params['mask_length']
    else:
        mask_length = 4
    lmadata, flashes = collect_output(outfile, mask_length=mask_length)#, min_points=min_points)

    if retain_ascii_output==False:
        os.remove(outfile)
    else:
        logger.info('Wrote ascii event and flash data to %s' % outfile)
    if cleanup_tmp == True:
        result = cleanup_build(d)
    else:
        logger.warning('did not delete flashsort directory in /tmp')
        
    return lmadata, flashes


def test_output():
    import tables as T
    h5 = T.openFile('/Users/ebruning/out/LYLOUT_040526_213000_0600.dat.gz.flash.h5')
    flashes = h5.root.flashes.LMA_040526_213000_600
    events  = h5.root.events.LMA_040526_213000_600
    # flashes.cols.n_points[0:100]
    big = [fl['flash_id'] for fl in flashes if fl['n_points'] > 100]
    a_flash = big[0]
    points = [ev['lat'] for ev in events if ev['flash_id'] == a_flash]
    print(flashes.cols.init_lon[0:10])



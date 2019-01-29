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
from lmatools.io.LMAarrayFile import cat_LMA



LOG_BASEFILENAME = datetime.datetime.now().strftime('Flash-autosort.log')



def logger_setup(logpath):
    logger = logging.getLogger('FlashAutorunLogger')
    logfile = os.path.join(logpath, LOG_BASEFILENAME)
    loghandler = logging.handlers.RotatingFileHandler(logfile, maxBytes=1024*1024, backupCount=3)
    logger.addHandler(loghandler)
    logger.setLevel(logging.DEBUG)


def write_output(outfile, flashes, orig_LMA_file, metadata=None):
    from lmatools.io.LMA_h5_write import write_h5
    if metadata is None:
        # use metadata from the first flash as the canonical metadata, 
        #   since all flashes were sorted fromt the same LYLOUT file
        # This breaks in the case of empty flashes; in that case the calling function should pass in metadata.
        metadata = flashes[0].metadata
    write_h5(outfile, flashes, metadata, orig_LMA_file)


def run_files_with_params(files, output_path, params, clusterer=None, min_points=1, retain_ascii_output=True, cleanup_tmp=True):
    if clusterer is None:
        from .autorun_mflash import build, cleanup_build, collect_output
        from .autorun_mflash import cluster
        clusterer = cluster

    logger = logging.getLogger('FlashAutorunLogger')
    
    now = datetime.datetime.now().strftime('Flash autosort started %Y%m%d-%H%M%S')
    logger.info(now)
    
    # Calculate the number of header lines based on the largest data file.
    
    f_sizes = [os.path.getsize(f) for f in files]
    largest_file_index = f_sizes.index(max(f_sizes))
    lma_pipe, command, any_input = cat_LMA(files[largest_file_index])
    lma_text, err = lma_pipe.communicate(input=any_input)
    lma_text = lma_text.decode()
    isDataLine = r"^.*\*+.*data.*\*+.*" #Search for asterisks, data, and asterisks
    matchDataLine = re.compile(isDataLine, re.IGNORECASE)
    split_lma_text = lma_text.split('\n')
    for lineIdx, line in enumerate(split_lma_text):
        if matchDataLine.search(line):
            params['nhead'] = lineIdx+1
            logger.info("Header is %d lines. This length will be used for all files this run." % (params['nhead'],))
            break
            
    # We could parse for the number of sources from the header, but that is wrong sometimes. 
    # Instead, set the number of points to the total number of lines in the file minus the header.
    # Add 10% to the total. The largest file might actually not be the largest when uncompressed due to variable packing efficiency
    params['n_sources'] = int(1.10*(len(split_lma_text) - params['nhead']))
    logger.info('Calculated max source count for this run: {0}'.format(params['n_sources']))
    del lma_text, split_lma_text
    # lma_pipe.close()
    
    logger.info('%s', params)

    h5_outfiles = []
    for a_file in files:
        try:
            file_base_name = os.path.split(a_file)[-1].replace('.gz', '')
            outfile = os.path.join(output_path, file_base_name+'.flash')
            
            # clusterer should use the name outfile as the base for any, e.g., ASCII data it would like to save
            lmadata, flashes = clusterer(a_file, output_path, outfile, params, logger,
                       min_points=min_points, retain_ascii_output=retain_ascii_output, cleanup_tmp=cleanup_tmp )
                        
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


if __name__ == '__main__':
    
    logger_setup('.')
    
    DClat = 38.8888500 # falls church / western tip of arlington, rough centroid of stations
    DClon = -77.1685800
    ctrLat = DClat
    ctrLon = DClon
    # kounLat = 35.23833
    # kounLon = -97.46028
    
    
    params = {'stations':(10,13),
              'chi2':(0,1.0),
              'ascii_flashes_out':'flashes_out.dat',
              'ctr_lat':ctrLat, 'ctr_lon':ctrLon,
              }

    outpath = sys.argv[1]
    files = sys.argv[2:]
    
    # outpath = '/Users/ebruning/out/'
    # files = [# '/data/20040526/LMA/LYLOUT_040526_224000_0600.dat.gz',
    #          # '/data/20040526/LMA/LYLOUT_040526_211000_0600.dat.gz',
    #          # '/data/20040526/LMA/LYLOUT_040526_212000_0600.dat.gz',
    #          '/data/20040526/LMA/LYLOUT_040526_213000_0600.dat.gz',
    #         ]
    # python autorun.py "`pwd`/tmpout" /data/20090329/LYLOUT_090329_200000_3600.dat.gz
    run_files_with_params(files, outpath, params, retain_ascii_output=False, cleanup_tmp=True)
    
    # from initialization import write_header
    # write_header('initialization.h', stations=(0,13))
    # d = build()
    # print 'Built flash program in', d
    # test_file = '/data/20090329/LYLOUT_090329_200000_3600.dat.gz'
    # test_file = '/data/20040526/LMA/LYLOUT_040526_213000_0600.dat.gz'
    # sort_file(test_file, d)
    # print 'Ran file ', test_file
    # flf, flashes = collect_output(d, min_points=3)
    # write_output('test.h5', flashes, test_file)
    
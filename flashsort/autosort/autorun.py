import os, sys, re
import glob
import tempfile
import shutil
import subprocess
import logging, logging.handlers
import datetime

from flash_stats import calculate_flash_stats
from mflash import write_header

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
    
build_files = ( 'mflashcustom.f',
                'makefile',
                )
                
flash_output_dir = 'output'     # created within the build directory
flash_prg_name   = 'mflashcustom'   # created within the build directory

tmp_flashsort_prepend = 'flashsort_'

def build(directory=None):
    """ Build the flash program in directory; a temporary directory is created by default. """


    if directory == None:
        d = tempfile.mkdtemp('', tmp_flashsort_prepend)
    else:
        d = directory
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
    

def fetch_gzipdata_from_url(the_url):
    """ Returns unzipped data inside a gzipped file residing at the_url.

        From http://diveintopython.org/http_web_services/gzip_compression.html
    """
    import urllib2, httplib
    import StringIO
    import gzip

    httplib.HTTPConnection.debuglevel = 1
    request = urllib2.Request(the_url)
    request.add_header('Accept-encoding', 'gzip')
    opener = urllib2.build_opener()
    f = opener.open(request)
    compresseddata = f.read()
    compressedstream = StringIO.StringIO(compresseddata)
    gzipper = gzip.GzipFile(fileobj=compressedstream)      
    data = gzipper.read()
    gzipper.close()                            
    return data

    
def cat_LMA(filename):
    """ Returns subprocess pipe with LMA data file on stdout """
    if filename.find('http://') >= 0:
        # data = fetch_gzipdata_from_url(filename)
        url_copy = subprocess.Popen(['curl', '-s', filename], stdout=subprocess.PIPE)#, stdin=subprocess.PIPE)
        command = ['gunzip', '-c']
        f = subprocess.Popen(command, stdin=url_copy.stdout, stdout=subprocess.PIPE)
        any_input=None
    else:
        if filename.find('.gz') >= 0:
            command = ['gunzip', '-c', filename]
        else:
            command  = ['cat', filename]
        f = subprocess.Popen(command, stdout=subprocess.PIPE)
        any_input = None
    return f, command, any_input

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
    
class Flash(object):
    def __init__(self, points):
        self.points = points

class FlashMetadata(object):

    def __init__(self, headerText):
        #Search the header for info on how the data is written

        self.header = headerText

        isColumnHeaderLine = r"^Data:(.*)"
        matchDataFormatLine = re.compile(isColumnHeaderLine, re.IGNORECASE | re.MULTILINE)

        isDataStartTime = r"^Data start time:(.*)"
        matchDataStartTimeLine = re.compile(isDataStartTime, re.IGNORECASE | re.MULTILINE)

        secAnalyzedLine = r"^Number of seconds analyzed:(.*)"
        matchSecAnalyzedLine = re.compile(secAnalyzedLine, re.IGNORECASE | re.MULTILINE)


        startTimeMatch = matchDataStartTimeLine.search(headerText)
        if startTimeMatch:
            #Looking to deal with something like: " 06/28/04 23:50:00"
            dateAndTime = startTimeMatch.group(1).split()
            self.startmonth, self.startday, self.startyear = [ int(datePart) for datePart in dateAndTime[0].split('/') ]
            self.starthour, self.startminute, self.startsecond = [ int(timePart) for timePart in dateAndTime[1].split(':') ]
            if self.startyear < 1000 and self.startyear > 70:
                self.startyear += 1900
            else:
                self.startyear += 2000

        secAnalyzedMatch=matchSecAnalyzedLine.search(headerText)
        if secAnalyzedMatch:
            self.sec_analyzed = int(secAnalyzedMatch.group(1))


        formatMatch=matchDataFormatLine.search(headerText)
        if formatMatch:
            columns = formatMatch.group(1).split(',')
            self.columns = [columnName.strip() for columnName in columns]
    

def collect_output(datafile, min_points=1):
    """ collect all output from the flash program output directory created by
        the flash program in flash_dir and calculate some additional stats on each flash        
        """
    from LMAarrayFile import LMAdataFile
    import numpy as np
    
    logger = logging.getLogger('FlashAutorunLogger')
    
    # outdir = os.path.join(flash_dir,flash_output_dir)
    # os.chdir(outdir)
    
    lma = LMAdataFile(datafile)
    
    # get the mapping from flash_ids to the points
    order = np.argsort(lma.flash_id)
    flid = lma.flash_id[order]
    boundaries, = np.where(flid[1:]-flid[:-1])    # where indices are nonzero
    boundaries = np.hstack(([0], boundaries+1))
    
    all_data = lma.data[order]
    
    max_idx = len(flid) #- 1
    slice_lower_edges = tuple(boundaries)
    slice_upper_edges = slice_lower_edges[1:] + (max_idx,)
    slices = zip(slice_lower_edges, slice_upper_edges)
    
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


def write_output(outfile, flashes, orig_LMA_file):
    from write_flashes import write_h5
    # use metadata from the first flash as the canonical metadata, 
    #   since all flashes were sorted fromt the same LYLOUT file
    write_h5(outfile, flashes, flashes[0].metadata, orig_LMA_file)


def run_files_with_params(files, output_path, params, min_points=1, retain_ascii_output=True, cleanup_tmp=True):
    logger = logging.getLogger('FlashAutorunLogger')
    
    now = datetime.datetime.now().strftime('Flash autosort started %Y%m%d-%H%M%S')
    logger.info(now)
    
    # Calculate the number of header lines based on the first data file.
    lma_pipe, command, any_input = cat_LMA(files[0])
    lma_text, err = lma_pipe.communicate(input=any_input)
    isDataLine = r"^.*\*+.*data.*\*+.*" #Search for asterisks, data, and asterisks
    matchDataLine = re.compile(isDataLine, re.IGNORECASE)
    for lineIdx, line in enumerate(lma_text.split('\n')):
        if matchDataLine.search(line):
            params['nhead'] = lineIdx+1
            logger.info("Header is %d lines. This length will be used for all files this run." % (params['nhead'],))
            break
    del lma_text
    # lma_pipe.close()
    
    logger.info('%s', params)

    h5_outfiles = []
    for a_file in files:
        try:
            write_header(os.path.join(src_dir,'mflashcustom.f'), **params)

            d = build() # returns temporary directory name
            logger.debug('Built flash program in %s' % d)
        
            logger.info('Sorting flashes for %s' % a_file)
            out, err = sort_file(a_file, d) # out contains sorted flash data (via stdout)
            # print out
        
        
            file_base_name = os.path.split(a_file)[-1].replace('.gz', '')
            # print 'filebase =  ', file_base_name
            # outfile = os.path.join(output_path, file_base_name+'.flash.h5')
            # write_output(outfile, flashes, file_base_name)
            # print 'Wrote original data and flashes to ', outfile
        
            outfile = os.path.join(output_path, file_base_name+'.flash')
            shutil.copyfile(os.path.join(d, 'flashes_out.dat'),
                            outfile)            
            
            lmadata, flashes = collect_output(outfile)#, min_points=min_points)
            outfile_with_extension = outfile + '.h5'
            h5_outfiles.append(outfile_with_extension)
            write_output(outfile_with_extension, flashes, a_file)

            if retain_ascii_output==False:
                os.remove(outfile)
            else:
                logger.info('Wrote ascii event and flash data to %s' % outfile)
            if cleanup_tmp == True:
                result = cleanup_build(d)
            else:
                logger.warning('did not delete flashsort directory in /tmp')
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
    print flashes.cols.init_lon[0:10]


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
    
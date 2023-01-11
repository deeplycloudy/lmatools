from __future__ import absolute_import
from __future__ import print_function
import os, re #, gzip
import numpy as np
import logging
import subprocess

logger = logging.getLogger('FlashAutorunLogger')

from numpy.lib.recfunctions import append_fields

def cat_LMA(filename):
    """ Returns subprocess pipe with LMA data file on stdout """
    if filename.find('http://') >= 0:
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


def dec2bin(dec_number):
    if dec_number == 0: return '0'
    return (dec2bin(dec_number >> 1) + str(dec_number % 2))


def countBits(values):
    # bit shifting routines are in numpy 1.4
    from numpy import array, left_shift, right_shift

    v = array(values).astype('uint32')


    # Bit counting for a 32 bit unsigned integer.
    # there is a fully generic version of this method at
    # http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
    # Binary magic numbers method, also found in  in pages 187-188 of Software Optimization Guide for AMD Athlon 64 and Opteron Processors.

    # The C version is:
    # v = v - ((v >> 1) & 0x55555555);                    # reuse input as temporary
    # v = (v & 0x33333333) + ((v >> 2) & 0x33333333);     # temp
    # c = ((v + (v >> 4) & 0xF0F0F0F) * 0x1010101) >> 24; # count

    fives  = int('0x55555555', base=16)
    threes = int('0x33333333', base=16)
    effs   = int('0xF0F0F0F',  base=16)
    ones   = int('0x1010101',  base=16)

    v = v - ( (right_shift(v, 1)) & fives);                        # reuse input as temporary
    v = (v & threes) + ( (right_shift(v,2)) & threes);             # temp
    c = right_shift(((v + (right_shift(v,4)) & effs) * ones), 24); # count

    return c

def mask_to_int(mask):
    """ Convert object array of mask strings to integers"""
    mask = np.atleast_1d(mask)
    if len(mask.shape) == 0:
        mask_int = np.asarray([], dtype=int)
    else:
        try:
            # mask is a plain integer
            mask_int = np.fromiter((int(v) for v in mask), int)
        except ValueError:
            # mask is a string representing a base-16 (hex) number
            mask_int = np.fromiter((int(v,16) for v in mask), int)
    return mask_int

class LMAdataFile(object):

    def __init__(self, filename, iterator=False):
        """ iterator=True makes self.data an iterator yielding individual event records"""
        self.filename=filename
        self.iterator=iterator
        self.field_names = None
        self.field_formats = None
        self.field_converters = None
        self.header=[]
        self.startmonth=1
        self.startday=1
        self.startyear=0
        self.mask_length = 0
        self.starthour, self.startminute, self.startsecond = 0, 0, 0

        # Could read the mask length here, but that would mean opening/closing the
        #   file twice. Use this if there are issues
        # if self.filename.find('.gz') >= 0:
        #     with gzip.open(filename) as f:
        #         for line_no, line in enumerate(f):
        #             if line.startswith(b'Data format:'):
        #                 self.mask_length = int(line.decode().split(' ')[-1][:-2])
        #                 # print (file_mask_length)
        #                 break
        # f.close()

        # This dictionary is used to match getattr calls for lat, lon, etc.
        #   with column names as indicated in the LMA data file header
        #   To a user of this object, the dictionary keys become attribute
        #   names that she can access.
        #   Add a new entry to the list if the header column names change.
        self.columnNames = {    'time':     ['time (UT sec of day)', 'time'],
                                'lat':      ['lat'],
                                'lon':      ['lon'],
                                'chi2':     ['reduced chi^2', 'chi2'],
                                'power':    ['P(dBW)'],
                                'mask':     ['mask'],
                                'alt':      ['alt(m)', 'alt'],
                                'stations': ['stations', '# of stations contributed'],
                                'charge':   ['charge'],
                                'range':    ['range (km)'],
                                'flash_id': ['flash_id'],
                                'point_id': ['point_id'],
                            }

        self.columnTypes = {
                                'time':      float,
                                'lat':       float,
                                'lon':       float,
                                'chi2':      float,
                                'power':     float,
                                'mask':      'S{0:d}'.format(2),
                                # 'mask':      'S{0:d}'.format(mask_length),
                                # 'mask':      object,
                                # 'mask':      str,
                                'alt':       float,
                                'stations':  int,
                                'charge':    int  ,
                                'range':     float,
                                'flash_id':  int,
                                'point_id':  int,

                            }

        self.columnConverters = {
                                # 'mask':      mask_to_int,
                            }


        if filename: self.read()


    def __getattr__(self, attrname):

        # See __init__ for column names
        try:
            return self.data[attrname]
        except:
            pass

        # If we got here, stations column wasn't in file.
        #   Try getting it from station mask.
        if attrname=='stations':
            stations = self.hexMaskToStationCount()
            # placing self.data in a list due to this bug
            # http://stackoverflow.com/questions/36440557/typeerror-when-appending-fields-to-a-structured-array-of-size-one
            self.data = append_fields([self.data], ('stations',), (stations,))
            return stations

        return None
        #raise AttributeError, attrname

    def hexMaskToStationCount(self, mask=None):
        if mask is None:
            mask=self.mask
        if mask is None: return None

        mask = mask_to_int(mask)
        try:
            stationCount = countBits(mask)
        except:
            logger.warning("Using slow method of getting station counts")
            stationCount = [ dec2bin(onemask).count('1') for onemask in mask]
            #stationCount = map( (lambda onemask: (dec2bin(int(onemask, 16))).count('1') ), mask)

        self.stations = stationCount
        return stationCount

    def get_file_obj(self, notify=True):
        if not (self.filename.find('.dat') >= 0):
            raise IOError("Can only read .dat files")

        if notify:
            logger.info("Loading LMA data from " + self.filename)

        if self.filename.find('.gz') >= 0:
            command = 'gunzip -c ' + self.filename
            thefile = os.popen(command)
            #gzip readlines is slooow ... at least 5x worse.
            #thefile=gzip.GzipFile(self.filename)
        else:
            thefile=open(self.filename, 'r')

        return thefile

    def read(self):


        # if not (self.filename.find('.dat') >= 0):
        #     raise "Can only read .dat files"
        #
        # print "Loading LMA data from " + self.filename
        #
        # if self.filename.find('.gz') >= 0:
        #     command = 'gzcat ' + self.filename
        #     thefile = os.popen(command)
        #     #gzip readlines is slooow ... at least 5x worse.
        #     #thefile=gzip.GzipFile(self.filename)
        # else:
        #     thefile=file(self.filename, 'r')
        thefile = self.get_file_obj()

        #foundDataLine=False
        isDataLine = r"^.*\*+.*data.*\*+.*" #Search for asterisks, data, and asterisks
        matchDataLine = re.compile(isDataLine, re.IGNORECASE)

        # fileLines = thefile.readlines()
        #First read in the header, then the data
        #dataTemp=[]
        # lineIdx = 0
        for lineIdx, line in enumerate(thefile):
            #if not foundDataLine:
            self.header.append(line)

            #Once the ***data*** line is found,
            #   need to treat all lines as data
            if matchDataLine.search(line):
                #print 'Found data'
                #foundDataLine=True
                break
            #else:
            #   dataTemp.append(line.split())

        #Search the header for info on how the data is written
        isColumnHeaderLine = r"^Data:(.*)"
        matchDataFormatLine = re.compile(isColumnHeaderLine, re.IGNORECASE)

        DataFormatLine = r"^Data format:(.*)"
        matchDataFormatLines = re.compile(DataFormatLine, re.IGNORECASE | re.MULTILINE)

        isDataStartTime = r"^Data start time:(.*)"
        matchDataStartTimeLine = re.compile(isDataStartTime, re.IGNORECASE)

        secAnalyzedLine = r"^Number of seconds analyzed:(.*)"
        matchSecAnalyzedLine = re.compile(secAnalyzedLine, re.IGNORECASE | re.MULTILINE)

        for line in self.header:

            startTimeMatch = matchDataStartTimeLine.search(line)
            if startTimeMatch:
                # Looking to deal with something like: " 06/28/04 23:50:00"
                dateAndTime = startTimeMatch.group(1).split()
                self.startmonth, self.startday, self.startyear = [ int(datePart) for datePart in dateAndTime[0].split('/') ]
                self.starthour, self.startminute, self.startsecond = [ int(timePart) for timePart in dateAndTime[1].split(':') ]
                if self.startyear < 1000 and self.startyear > 70:
                    self.startyear += 1900
                else:
                    self.startyear += 2000
            formatMatch=matchDataFormatLine.search(line)
            if formatMatch:
                columns = formatMatch.group(1).split(',')
                columns = [columnName.strip() for columnName in columns]

            secAnalyzedMatch=matchSecAnalyzedLine.search(line)
            if secAnalyzedMatch:
                self.sec_analyzed = int(secAnalyzedMatch.group(1))
            DataFormatMatchs=matchDataFormatLines.search(line)
            if DataFormatMatchs:
                self.mask_length = int(DataFormatMatchs.group(1).split(' ')[-1][:-1])
                self.columnTypes['mask'] = 'S{0:d}'.format(self.mask_length)


        n_columns = len(columns)
        field_names, types = list(range(n_columns)), list(range(n_columns))
        converters = {}
        for column_idx, column in enumerate(columns):
            for field in self.columnNames.keys():
                if column in self.columnNames[field]:
                    field_names[column_idx] = field
                    typecode = self.columnTypes[field]
                    types[column_idx]  = typecode
                    try:
                        converters[column_idx] = self.columnConverters[field]
                    except:
                        pass

        self.field_names   = field_names
        self.field_formats = types
        self.field_converters = converters


        thefile.close()
        thefile = self.get_file_obj(notify=False)


        if not self.iterator:
            # creates a rec-array
            self.data = np.loadtxt(thefile, dtype={'names':field_names, 'formats':types},
                                            converters=converters, skiprows=len(self.header))
            if ('mask' in self.field_names) & ('station' not in self.field_names):
                # Read the stations property, which triggers checking for the
                # adding the column to the data structure
                stations = self.stations
            # Done with file IO
            thefile.close()
        else:
            self._file_object = thefile
            self.data = self._data_record_iterator(skiprows=len(self.header))

        # print "LMA data loaded."

    def _data_record_iterator(self, skiprows=0):
        if ('mask' in self.field_names) & ('station' not in self.field_names):
            calc_stations = True
        else:
            calc_stations = False
        line_count = 0
        for line in self._file_object:
            line_count+=1
            if line_count <= skiprows:
                logger.debug('skipping ' + line[0:20])
                continue
            items = line.split()
            vals = list(map(apply_format, items, self.field_formats))
            record = dict(list(zip(self.field_names, vals)))
            if calc_stations:
                stations = self.hexMaskToStationCount(mask=[record['mask']])
                record['stations'] = stations[0]
            yield record

def apply_format(item, format):
    if format is object: return item
    if (format == 'S4') or (format =='S6'):
        return str(item)
    return format(item)

def runtest():
    lma = LMAdataFile('tmpout/LYLOUT_090329_200000_3600.dat.flash')
    # lma = LMAdataFile('/data/2004/LMA/040526/LYLOUT_040526_230000_0600.dat.gz')
    #lma = LMAdataFile('/Users/ericbruning/othercode/LYLOUT_040526_220000_0600.dat')

    print(lma.startmonth, lma.startday, lma.startyear)
    print(lma.starthour, lma.startminute, lma.startsecond)
    # keys = lma.data.keys()
    # for key in keys:
    #     print key, len(lma.data[key])
    print(lma.lat)
    print(lma.stations)
    print(lma.flash_id)


if __name__ == '__main__':
    runtest()
    # import hotshot
    # from hotshot import stats
    # prof = hotshot.Profile("lmafile_test_profile_array")
    # prof.runcall(runtest)
    # prof.close()
    # s=stats.load("lmafile_test_profile_array")
    # s.sort_stats("time").print_stats()

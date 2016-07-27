from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import os
import datetime

PRELOADED_DATA = 'LYLOUT_111121_145001_0600.dat.npy'

def preload_some_data():
    from lmatools.io.LMAarrayFile import LMAdataFile
    lma = LMAdataFile('/data/WTLMA/FirstLightning/processed7stn/LYLOUT_111121_145001_0600.dat.gz')
    np.save(PRELOADED_DATA, lma.data)

def recycle_LMA_file(duration, seconds_offset):
    import cStringIO
    
    try:
        data = np.load(PRELOADED_DATA)
    except:
        print("Regenerating preloadable LMA data")
        preload_some_data()
        data = np.load(PRELOADED_DATA)
        
    data['time'] -= data['time'][0]
    avail_duration = data['time'][-1] - data['time'][0]
    while (avail_duration < duration):
        # Keep extending the array with itself until we have enough data
        start_len = data.shape[0]
        t_max = data['time'][-1]
        data = np.hstack((data,data))
        data['time'][start_len:] += t_max
        avail_duration = data['time'][-1] - data['time'][0]
    data['time'] += seconds_offset
    in_duration = (data['time'] > seconds_offset) & (data['time'] < (seconds_offset+duration))
    data = data[in_duration]
    n_points = data.shape[0]
    
    # Data: time (UT sec of day), lat, lon, alt(m), reduced chi^2, P(dBW), mask
    # output_format = "{time:15.9f} {lat:10.6f} {lon:11.6f} {alt:7.1f} {chi2:5.2f} {power:5.1f} {mask:s}"
    output_format = "%15.9f %10.6f %11.6f %7.1f %5.2f %5.1f %4s"
    str_out = cStringIO.StringIO()
    np.savetxt(str_out, data, fmt=output_format)
    output_data_string = str_out.getvalue()
    str_out.close()

    return n_points, output_data_string

def fake_LMA_file(year=2012, month=2, day=28, hour=0, minute=0, second=0, 
                  duration=60, outpath=None, outfile=None, gzip=False, header_template="", event_generator=None):
    """ event_rate: events per second
        event_generator: called with (duration, seconds_offset). Is expected to return (n_points, event_multiline_string)
        """

    start_day   = datetime.datetime(year,month,day,0,0,0)
    start_time = datetime.datetime(year,month,day,hour,minute,second)
    
    seconds_offset = (start_time - start_day).total_seconds()
    
    now = datetime.datetime.now()
    
    header_data = { 'analysis_start':now.strftime('%c'),
                    'analysis_end':now.strftime('%c'),
                    'data_start':start_time.strftime('%m/%d/%y %H:%M:%S'),
                    'duration':duration,
                    'location':'WTLMA {0}'.format(start_time.strftime('%Y')),
                    
                    }
                    
    if event_generator is None:
        n_points, output_data_string = recycle_LMA_file(duration, seconds_offset)
    else:
        n_points, output_data_string = event_generator(duration, seconds_offset)

    if outfile == None:
        outfile = "LYLOUT_{0:s}_{1:04d}.dat".format(start_time.strftime("%y%m%d_%H%M%S"), duration)
        
    if outpath is not None:
        outfile = os.path.join(outpath, outfile)

    header_data['n_points'] = n_points
    header_data['analysis_end'] = datetime.datetime.now().strftime('%c')
    
    f_out = open(outfile, 'w')
    f_out.write(header_template.format(**header_data))
    f_out.write(output_data_string)
    f_out.close()
    
    # print header_template.format(**header_data)
    return outfile


real_analysis_program_line =  "Analysis program: /data/rtlma/bin/lma_analysis-10.8.0/lma_analysis -d 20111121 -t 145001 -s 600 -l wt -a -n 5 -o /data/rtlma/processed_data/20111121 -q -v -x 5.00 -y 500.00 [datafiles]"


late2011_header = """New Mexico Tech Lightning Mapping Array - analyzed data
Analysis program: /path/to/fakeLMA.py
Analysis program version: 10.8.0
Analysis started: {analysis_start}
Analysis finished: {analysis_end}
Data start time: {data_start}
Number of seconds analyzed: {duration:d}
Location: {location}
Coordinate center (lat,lon,alt): 33.6000000 -101.8000000 984.00
Number of stations: 8
Number of active stations: 7
Active stations: E G W N O R L
Minimum number of stations per solution: 5
Maximum reduced chi-squared: 5.00
Maximum number of chi-squared iterations: 20
Station information: id, name, lat(d), lon(d), alt(m), delay(ns), board_rev, rec_ch
Sta_info: E  E                  33.6000000  -101.8000000   984.00   26 3  3
Sta_info: G  I                  33.7000000  -101.6000000   992.00   26 3  3
Sta_info: W  L                  33.4000000  -101.7000000   956.85   26 3  3
Sta_info: B  B                  33.7000000  -102.0000000  1007.59   26 3  3
Sta_info: N  N                  33.7000000  -101.8000000   998.45   26 3  3
Sta_info: O  R                  33.5000000  -102.0000000  1018.00   26 3  3
Sta_info: R  R                  33.5000000  -101.6000000   960.00   26 3  3
Sta_info: L  L                  33.6000000  -101.5000000   956.00   26 3  3
Station data: id, name, win(us), dec_win(us), data_ver, rms_error(ns), sources, %, <P/P_m>, active
Sta_data: E  E                   80    80   8  70   510143  87.4  1.00   A
Sta_data: G  I                   80    80   8  70   495303  84.9  0.65   A
Sta_data: W  L                   80    80   8  70   493903  84.6  1.70   A
Sta_data: B  B                    0     0   0  70        0   0.0   0.0  NA
Sta_data: N  N                   80    80   8  70   497490  85.2  1.57   A
Sta_data: O  R                   80    80   8  70   523372  89.7  0.93   A
Sta_data: R  R                   80    80   8  70   385601  66.1  0.76   A
Sta_data: L  L                   80    80   8  70   162450  27.8  0.33   A
Metric file version: 2
Initial metric: lowest of (c*delta_t * (1 + rchi2^2/4.0))
Stop-early combination search: 5 metric repeats
Final metric: lowest of initial_metric * fuzzy_factor * max (1, rchi2/1.0)
Minimum forego ratio (unique/total combinations): 0.699
Fuzzy table factors: log base 2
Fuzzy table: repeats(rows) 1 to 7, num_stations(cols) 4 to 12
Fuzzy_data:   15.0    7.5    0.0   -7.5  -15.0  -22.5  -30.0  -37.5  -45.0
Fuzzy_data:   14.0    6.5   -1.0   -8.5  -16.0  -23.5  -31.0  -38.5  -46.0
Fuzzy_data:   13.0    5.5   -2.0   -9.5  -17.0  -24.5  -32.0  -39.5  -47.0
Minimum forego ratio (unique/total combinations): 0.699
Fuzzy table factors: log base 2
Fuzzy table: repeats(rows) 1 to 7, num_stations(cols) 4 to 12
Fuzzy_data:   15.0    7.5    0.0   -7.5  -15.0  -22.5  -30.0  -37.5  -45.0
Fuzzy_data:   14.0    6.5   -1.0   -8.5  -16.0  -23.5  -31.0  -38.5  -46.0
Fuzzy_data:   13.0    5.5   -2.0   -9.5  -17.0  -24.5  -32.0  -39.5  -47.0
Fuzzy_data:   12.0    4.5   -3.0  -10.5  -18.0  -25.5  -33.0  -40.5  -48.0
Fuzzy_data:   11.0    3.5   -4.0  -11.5  -19.0  -26.5  -34.0  -41.5  -49.0
Fuzzy_data:   10.0    2.5   -5.0  -12.5  -20.0  -27.5  -35.0  -42.5  -50.0
Fuzzy_data:    9.0    1.5   -6.0  -13.5  -21.0  -28.5  -36.0  -43.5  -51.0
Discarded solutions (stations,repeats):
Station mask order: LRONBWGE
Data: time (UT sec of day), lat, lon, alt(m), reduced chi^2, P(dBW), mask
Data format: 15.9f 10.6f 11.6f 7.1f 5.2f 5.1f 4x
Number of events: {n_points:d}
*** data ***
"""



if __name__ == '__main__':
    import argparse

    now = datetime.datetime.now()
    
    parser = argparse.ArgumentParser(description="Generate fake LMA data for arbitrary durations and days.\nExample Usage: python fakeLMA.py -Y 2012 -m5 -d26 -H23 -M59")
    parser.add_argument('-Y', type=int, dest='year', default=now.year,
                        help='Year, 4-digit (e.g., 2012)')
    parser.add_argument('-m', type=int, dest='month', default=now.month,
                        help='Month, digit  (e.g., 2)')
    parser.add_argument('-d', type=int, dest='day', default=now.day,
                        help='Day, digit (e.g., 29)')
    parser.add_argument('-H', type=int, dest='hour', default=now.hour,
                        help='Hour, digit, 24-hour  (e.g., 13)')
    parser.add_argument('-M', type=int, dest='minute', default=now.minute,
                        help='Minute, digit (e.g., 50)')
    parser.add_argument('-S', type=int, default=0, dest='second',
                        help='Second, digit  (e.g., 0)')
    parser.add_argument('--duration', type=int, default=60, dest='duration',
                        help='Duration of data in seconds')


    args = parser.parse_args()
    fake_LMA_file(year=args.year, month=args.month, day=args.day, hour=args.hour, minute=args.minute, second=args.second, 
                    duration=args.duration, header_template=late2011_header)


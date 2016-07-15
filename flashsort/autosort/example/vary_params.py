from __future__ import absolute_import
import numpy as np
import os

from autosort.autorun import run_files_with_params

    
# DClat = 38.8888500 # falls church / western tip of arlington, rough centroid of stations
# DClon = -77.1685800
kounLat = 35.23833
kounLon = -97.46028


params = {'stations':(7,13),
          'chi2':(0,2.0),  
          'lat':(-90., 90.),
          'lon':(-180., 180.),
          'alt':(0.,20000.),     
          'distance':3000.0, 'thresh_duration':3.0, 
          'thresh_critical_time':0.15, 'merge_critical_time':0.5, 
          'ascii_flashes_out':'flashes_out.dat',
          'ctr_lat':kounLat, 'ctr_lon':kounLon,
          }

files = [# '/data/20040526/LMA/LYLOUT_040526_224000_0600.dat.gz',
       # '/data/20040526/LMA/LYLOUT_040526_211000_0600.dat.gz',
       # '/data/20040526/LMA/LYLOUT_040526_212000_0600.dat.gz',
       # '/data/20040526/LMA/LYLOUT_040526_213000_0600.dat.gz',
       # '/data/20040526/LMA/LYLOUT_040526_230000_0600.dat.gz',
       # '/data/20040526/LMA/LYLOUT_040526_231000_0600.dat.gz',
       # '/data/20040526/LMA/LYLOUT_040526_232000_0600.dat.gz',
       # '/data/20040526/LMA/LYLOUT_040526_233000_0600.dat.gz',
       # '/data/20040526/LMA/LYLOUT_040526_234000_0600.dat.gz',
       # '/data/20040526/LMA/LYLOUT_040526_235000_0600.dat.gz',
       # '/data/20040526/LMA/LYLOUT_040527_000000_0600.dat.gz',
       # '/data/20040526/LMA/LYLOUT_040527_001000_0600.dat.gz',
       # '/data/20040526/LMA/LYLOUT_040527_002000_0600.dat.gz',
       '/data/20040526/LMA/LYLOUT_040527_003000_0600.dat.gz',
       # '/data/20040526/LMA/LYLOUT_040527_004000_0600.dat.gz',
       # '/data/20040526/LMA/LYLOUT_040527_005000_0600.dat.gz',
      ]


base_out_dir = '/Users/ebruning/out/flash_sort/'

params['lat'] = (35., 36.)

for thresh_critical_time in np.arange(0.15, 0.35, 0.05):
    for distance in np.arange(500, 10000, 500):
        params['thresh_critical_time'] = thresh_critical_time
        params['distance'] = distance
        
        tag = 'thresh-%s_dist-%s' % (thresh_critical_time, distance)
        outdir = os.path.join(base_out_dir, tag)
        os.mkdir(outdir)
        
        info = open(os.path.join(outdir, 'input_params.py'), 'w')
        info.write(str(params))
        info.close()
        
        run_files_with_params(files, outdir, params)

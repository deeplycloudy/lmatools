# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>



from lmatools.flashsort.autosort.autorun import run_files_with_params, logger_setup


outpath = '/Users/ebruning/out/scratch/'
logger_setup(outpath)


#files = ['/data/20040526/LMA/LYLOUT_040526_224000_0600.dat.gz', '/data/20040526/LMA/LYLOUT_040526_225000_0600.dat.gz' ]
files = ['/data/20040526/LMA/LYLOUT_040526_214000_0600.dat.gz', ]


ctrLat = 35.23833
ctrLon = -97.46028

params = {'stations':(6,13),
          'chi2':(0,1.0),
          'ascii_flashes_out':'flashes_out.dat',
          'ctr_lat':ctrLat, 'ctr_lon':ctrLon,
          'distance':3000.0, 'thresh_critical_time':0.15
          }


from lmatools.flashsort.autosort.autorun_sklearn import cluster
run_files_with_params(files, outpath, params, cluster, retain_ascii_output=False, cleanup_tmp=True)











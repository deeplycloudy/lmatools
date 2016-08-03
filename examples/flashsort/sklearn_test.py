# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
import pkg_resources, os
from lmatools.flashsort.autosort.autorun import run_files_with_params, logger_setup


outpath = '/data/tmp/lmatools_test/'
logger_setup(outpath)


# load 
def sample_data_paths():
    filenames = 'LYLOUT_140526_094000_0600.dat.gz', 'LYLOUT_140526_095000_0600.dat.gz'
    resource_package = 'lmatools'
    for filename in filenames:
        resource_path = os.path.join('sampledata','ASCII_solutions',filename)
        template = pkg_resources.resource_filename(resource_package, resource_path)
        yield template

files = [f for f in sample_data_paths()]
print(files)

ctrLat = 35.23833
ctrLon = -97.46028

params = {'stations':(6,13),
          'chi2':(0,1.0),
          'ascii_flashes_out':'flashes_out.dat',
          'ctr_lat':ctrLat, 'ctr_lon':ctrLon,
          'thresh_duration':3.0,
          'distance':3000.0, 'thresh_critical_time':0.15
          }


from lmatools.flashsort.autosort.autorun_sklearn import cluster
run_files_with_params(files, outpath, params, cluster, retain_ascii_output=False, cleanup_tmp=True)











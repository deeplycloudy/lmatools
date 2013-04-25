# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from lmatools.flashsort.autosort.autorun import run_files_with_params, test_output, logger_setup

# <codecell>

outpath = '/Users/ebruning/out/scratch/'
logger_setup(outpath)

# <codecell>

outpath = '.'

# <codecell>

#files = ['/data/20040526/LMA/LYLOUT_040526_224000_0600.dat.gz', '/data/20040526/LMA/LYLOUT_040526_225000_0600.dat.gz' ]
files = ['/data/20040526/LMA/LYLOUT_040526_214000_0600.dat.gz', ]

# <codecell>

ctrLat = 35.23833
ctrLon = -97.46028

params = {'stations':(6,13),
          'chi2':(0,1.0),
          'ascii_flashes_out':'flashes_out.dat',
          'ctr_lat':ctrLat, 'ctr_lon':ctrLon,
          }

# <codecell>

from lmatools.flashsort.autosort.autorun_sklearn import cluster
run_files_with_params(files, outpath, params, cluster, retain_ascii_output=False, cleanup_tmp=True)

# <codecell>


# <codecell>


# <codecell>



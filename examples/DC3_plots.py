from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import os
import sys
import glob
from datetime import datetime, timedelta
import subprocess

from lmatools.grid.make_grids import grid_h5flashfiles, dlonlat_at_grid_center, write_cf_netcdf_latlon
from lmatools.flashsort.autosort import autorun
from lmatools.vis.multiples_nc import make_plot

# Meant to be run on a just one day's worth of processed data

base_sort_dir = '/data/DC3/flash_sort_prelim/'
h5_dir = os.path.join(base_sort_dir, 'h5_files')
grid_dir = os.path.join(base_sort_dir, 'grid_files')
plot_dir = os.path.join(base_sort_dir, 'plots')
center_ID = 'COLMA'

centers = {'WTLMA':{'ctr_lat':33.606968, 'ctr_lon':-101.822625},
           'COLMA':{'ctr_lat':40.446398, 'ctr_lon':-104.636813},
           }

def tfromfile(name):
    parts = name.split('_')
    y, m, d = list(map(int, (parts[-3][0:2], parts[-3][2:4], parts[-3][4:6])))
    H, M, S = list(map(int, (parts[-2][0:2], parts[-2][2:4], parts[-2][4:6])))
    return y+2000,m,d,H,M,S

autorun.logger_setup(base_sort_dir)
    
params = {'stations':(6,13),
          'chi2':(0,5.0),  
          'lat':(-90., 90.),
          'lon':(-180., 180.),
          'alt':(0.,20000.),     
          'distance':3000.0, 'thresh_duration':3.0, 
          'thresh_critical_time':0.15, 'merge_critical_time':0.5, 
          'ascii_flashes_out':'flashes_out.dat',
          'ctr_lat':centers[center_ID]['ctr_lat'], 'ctr_lon':centers[center_ID]['ctr_lon'],
          'mask_length':6,
          }

files = sys.argv[1:]  # Pass in a list of LMA .dat.gz data files, including the full path name
y,m,d,H,M,S = tfromfile(files[0])
date = datetime(y,m,d, H,M,S)
# date = datetime.datetime(2012, 3, 19, 0) #for use for files not in naming nomenclature
base_out_dir = (h5_dir+"/20%s" %(date.strftime('%y/%b/%d')))

if os.path.exists(base_out_dir) == False:
    os.makedirs(base_out_dir)
    subprocess.call(['chmod', 'a+w', base_out_dir, h5_dir+'/20%s' %(date.strftime('%y/%b')), h5_dir+'/20%s' %(date.strftime('%y'))])

tag = ''
outdir = os.path.join(base_out_dir, tag) 
info = open(os.path.join(outdir, 'input_params.py'), 'w')
info.write(str(params))
info.close()
autorun.run_files_with_params(files, outdir, params, cleanup_tmp=False)

#for thresh_critical_time in np.arange(0.15, 0.20, 0.05):
#    for distance in np.arange(3000, 3500, 500):
#        params['thresh_critical_time'] = thresh_critical_time
#        params['distance'] = distance
#        tag = 'thresh-%s_dist-%s' % (thresh_critical_time, distance) 
#        #tag = ''
#        outdir = os.path.join(base_out_dir, tag)
#        os.mkdir(outdir)
#        info = open(os.path.join(outdir, 'input_params.py'), 'w')
#        info.write(str(params))
#        info.close()
#        autorun.run_files_with_params(files, outdir, params, cleanup_tmp=False)

""" The data used to make the grids are LMA data that have been sorted into flashes. The LMA source data and flash metadata are saved into an HDF5 file."""
#h5_filenames = glob.glob('/home/ebruning/Mar18-19/out/thresh-0.15_dist-3000/LYLOUT*.dat.flash.h5')
h5_filenames = glob.glob(h5_dir+'/20%s/LYLOUT*.dat.flash.h5' %(date.strftime('%y/%b/%d')))
h5_filenames.sort()




frame_interval=60.0*5
ctr_lat =  centers[center_ID]['ctr_lat']
ctr_lon =  centers[center_ID]['ctr_lon']
dx_km=3.0e3
dy_km=3.0e3
x_bnd_km = (-200e3, 200e3)
y_bnd_km = (-200e3, 200e3)

dx, dy, x_bnd, y_bnd = dlonlat_at_grid_center(ctr_lat, ctr_lon, 
                            dx=dx_km, dy=dy_km,
                            x_bnd = x_bnd_km, y_bnd = y_bnd_km )
print(("dx, dy = {0}, {1} deg".format(dx,dy)))
print(("lon_range = {0} deg".format(x_bnd)))
print(("lat_range = {0} deg".format(y_bnd)))

for f in h5_filenames:
    y,m,d,H,M,S = tfromfile(f)
    # print y,m,d,H,M,S
    start_time = datetime(y,m,d, H,M,S)
    end_time   = start_time + timedelta(0,600)
    print(start_time, end_time)
    

    outpath = grid_dir+'/20%s' %(date.strftime('%y/%b/%d'))
    if os.path.exists(outpath) == False:
        os.makedirs(outpath)
        subprocess.call(['chmod', 'a+w', outpath, grid_dir+'/20%s' %(date.strftime('%y/%b')), grid_dir+'/20%s' %(date.strftime('%y'))])
    grid_h5flashfiles(h5_filenames, start_time, end_time, frame_interval=frame_interval, proj_name='latlong',
                dx=dx, dy=dy, x_bnd=x_bnd, y_bnd=y_bnd, ctr_lon=ctr_lon, ctr_lat=ctr_lat, outpath = outpath,
                output_writer = write_cf_netcdf_latlon, output_filename_prefix=center_ID, spatial_scale_factor=1.0
                )

# Make plots
if False:
    n_cols=2
    mapping = { 'source':'lma_source',
                'flash_extent':'flash_extent',
                'flash_init':'flash_initiation',
                'footprint':'flash_footprint'}

    #nc_names = glob.glob('/home/ebruning/Mar18-19/grids/*.nc')
    nc_names = glob.glob(grid_dir+'/20%s/*.nc' %(date.strftime('%y/%b/%d')))
    nc_names.sort()
    outpath = plot_dir+'/20%s' %(date.strftime('%y/%b/%d'))
    if os.path.exists(outpath) == False:
        os.makedirs(outpath)
        subprocess.call(['chmod', 'a+w', outpath, plot_dir+'/20%s' %(date.strftime('%y/%b')), plot_dir+'/20%s' %(date.strftime('%y'))])

    for f in nc_names:
        gridtype = f.split('dx_')[-1].replace('.nc', '')
        var = mapping[gridtype]
        print(f)
        make_plot(f, var, n_cols=n_cols, outpath = outpath)
    
# make_plot('/data/rtlma/flash_sort/LMA_20120319_010000_600_10src_4000.0m-dx_flash_extent.nc', 'flash_extent', n_cols=n_cols)
# make_plot('/data/rtlma/flash_sort/LMA_20120319_010000_600_10src_4000.0m-dx_flash_init.nc', 'flash_initiation', n_cols=n_cols)
# make_plot('/data/rtlma/flash_sort/LMA_20120319_010000_600_10src_4000.0m-dx_source.nc', 'lma_source', n_cols=n_cols)
# make_plot('/data/rtlma/flash_sort/LMA_20120319_010000_600_10src_4000.0m-dx_footprint.nc', 'flash_footprint', n_cols=n_cols)

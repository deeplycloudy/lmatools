# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

from __future__ import absolute_import
import sys, os, glob
from datetime import datetime, timedelta
import subprocess

from lmatools.flashsort.gen_autorun import logger_setup, sort_files
from lmatools.flashsort.gen_sklearn import DBSCANFlashSorter

from lmatools.grid.make_grids import grid_h5flashfiles, dlonlat_at_grid_center, write_cf_netcdf_latlon, write_cf_netcdf_3d_latlon, write_cf_netcdf, write_cf_netcdf_3d
from lmatools.vis.multiples_nc import make_plot, make_plot_3d, read_file_3d


def get_sample_data_list():
    import lmatools.sampledata
    data_path = os.path.abspath(lmatools.sampledata.__path__[0])
    ascii_folder = 'ASCII_solutions' #inside data_path
    ascii_path = os.path.join(data_path, ascii_folder)
    filenames  = glob.glob(os.path.join(ascii_path, 'LYLOUT*.dat.gz'))
    return filenames

def tfromfile(name):
    parts = name.split('_')
    y, m, d = list(map(int, (parts[-3][0:2], parts[-3][2:4], parts[-3][4:6])))
    H, M, S = list(map(int, (parts[-2][0:2], parts[-2][2:4], parts[-2][4:6])))
    return y+2000,m,d,H,M,S    

def test_sort_and_grid_and_plot(outpath):
    """ Given an output path, run sample data included in lmatools through flash sorting and gridding"""
    base_sort_dir = outpath
    logger_setup(outpath)
    files = get_sample_data_list()
    center_ID = 'WTLMA'
    ctr_lat, ctr_lon = 33.5, -101.5
    params = {'stations':(6,13),
              'chi2':(0,1.0),
              'ctr_lat':ctr_lat, 'ctr_lon':ctr_lon,
              'distance':3000.0, 'thresh_critical_time':0.15,
              'thresh_duration':3.0,
              'mask_length':6,
              }
    h5_dir = os.path.join(base_sort_dir, 'h5_files')
    grid_dir = os.path.join(base_sort_dir, 'grid_files')
    plot_dir = os.path.join(base_sort_dir, 'plots')
    
    y,m,d,H,M,S = tfromfile(files[0])
    date = datetime(y,m,d, H,M,S)

    # Create HDF5 flash files
    base_out_dir = (h5_dir+"/20%s" %(date.strftime('%y/%b/%d')))
    if os.path.exists(base_out_dir) == False:
        os.makedirs(base_out_dir)
        subprocess.call(['chmod', 'a+w', base_out_dir, h5_dir+'/20%s' %(date.strftime('%y/%b')), h5_dir+'/20%s' %(date.strftime('%y'))])

    tag = ''
    outdir = os.path.join(base_out_dir, tag) 
    info = open(os.path.join(outdir, 'input_params.py'), 'w')
    info.write(str(params))
    info.close()
    
    if True:
        cluster = DBSCANFlashSorter(params).cluster
        sort_files(files, outdir, cluster)
    # Figure out which HDF5 files were created
    h5_filenames = glob.glob(h5_dir+'/20%s/LYLOUT*.dat.flash.h5' %(date.strftime('%y/%b/%d')))
    h5_filenames.sort()
    
    # Create NetCDF gridded data
    frame_interval=60.0*2 # seconds
    dx_km=3.0e3 # meters
    dy_km=3.0e3
    x_bnd_km = (-200e3, 200e3)
    y_bnd_km = (-200e3, 200e3)
    z_bnd_km = (0.0e3, 15.0e3)

    # There are similar functions in lmatools to grid on a regular x,y grid in some map projection.
    # dx, dy, x_bnd, y_bnd = dlonlat_at_grid_center(ctr_lat, ctr_lon,
    #                             dx=dx_km, dy=dy_km,
    #                             x_bnd = x_bnd_km, y_bnd = y_bnd_km )
    dx, dy, x_bnd, y_bnd = dx_km, dy_km, x_bnd_km, y_bnd_km
    # print("dx, dy = {0}, {1} deg".format(dx,dy))
    # print("lon_range = {0} deg".format(x_bnd))
    # print("lat_range = {0} deg".format(y_bnd))
    
    for f in h5_filenames:
        y,m,d,H,M,S = tfromfile(f)
        # print y,m,d,H,M,S
        start_time = datetime(y,m,d, H,M,S)
        end_time   = start_time + timedelta(0,600)
        # print start_time, end_time
    
        outpath = grid_dir+'/20%s' %(date.strftime('%y/%b/%d'))
        if os.path.exists(outpath) == False:
            os.makedirs(outpath)
            subprocess.call(['chmod', 'a+w', outpath, grid_dir+'/20%s' %(date.strftime('%y/%b')), grid_dir+'/20%s' %(date.strftime('%y'))])
        if True:
            grid_h5flashfiles(h5_filenames, start_time, end_time, frame_interval=frame_interval,
                    dx=dx, dy=dy, x_bnd=x_bnd, y_bnd=y_bnd, z_bnd=z_bnd_km,
                    ctr_lon=ctr_lon, ctr_lat=ctr_lat, outpath = outpath,
                    proj_name='aeqd',
                    output_writer = write_cf_netcdf,
                    output_writer_3d = write_cf_netcdf_3d,
                    output_filename_prefix=center_ID, spatial_scale_factor=1.0e-3,
                    # proj_name='latlong',
                    # output_writer = write_cf_netcdf_latlon,
                    # output_writer_3d = write_cf_netcdf_3d_latlon,
                    # output_filename_prefix=center_ID, spatial_scale_factor=1.0,
                    energy_grids = True
                    )
        
    # Create plots
    n_cols=2
    mapping = { 'source':'lma_source',
                'flash_extent':'flash_extent',
                'flash_init':'flash_initiation',
                'footprint':'flash_footprint',
                'specific_energy':'specific_energy',
                'flashsize_std':'flashsize_std',
                'total_energy': 'total_energy'
               }
    
    nc_names = glob.glob(grid_dir+'/20%s/*.nc' %(date.strftime('%y/%b/%d')))
    nc_names_3d = glob.glob(grid_dir+'/20%s/*_3d.nc' %(date.strftime('%y/%b/%d')))
    nc_names_2d = list(set(nc_names) - set(nc_names_3d))
    nc_names_2d.sort()
    nc_names_3d.sort()
    outpath = plot_dir+'/20%s' %(date.strftime('%y/%b/%d'))
    if os.path.exists(outpath) == False:
        os.makedirs(outpath)
        subprocess.call(['chmod', 'a+w', outpath, plot_dir+'/20%s' %(date.strftime('%y/%b')), plot_dir+'/20%s' %(date.strftime('%y'))])

    for f in nc_names_2d:
        gridtype = f.split('dx_')[-1].replace('.nc', '')
        var = mapping[gridtype]
        # make_plot(f, var, n_cols=n_cols, x_name='longitude', y_name='latitude', outpath = outpath)
        make_plot(f, var, n_cols=n_cols, x_name='x', y_name='y', outpath = outpath)
    
    for f in nc_names_3d:
        gridtype = f.split('dx_')[-1].replace('.nc', '').replace('_3d', '')
        var = mapping[gridtype]
        # grid_range = range_mapping[gridtype]
   
        ###Read grid files, then plot in either 2d or 3d space###
        # grid, grid_name, x, y, z, t, grid_t_idx, grid_x_idx, grid_z_idx = read_file_3d(f, var, x_name='longitude', y_name='latitude', z_name='altitude')
        grid, grid_name, x, y, z, t, grid_t_idx, grid_x_idx, grid_z_idx = read_file_3d(f, var, x_name='x', y_name='y', z_name='z')
        make_plot_3d(grid, grid_name, x, y, z, t, 
                     grid_t_idx, grid_x_idx, grid_z_idx, 
                     n_cols = n_cols, outpath = outpath)
        #, grid_range=grid_range)
    


if __name__ == '__main__':    
    outpath = sys.argv[1]
    outpath = os.path.abspath(outpath)
    test_sort_and_grid_and_plot(outpath)

#!/usr/bin/python

"""
Example script showing how to use lmatools to sort LMA ASCII data into flashes
and how to create gridded imagery from those flashes.

The block of code at the bottom shows how to call sorting and gridding functions.
The params dictionary controls the flash sorting parameters.

USAGE: 
python flash_sort_and_grid.py /path/to/output/ lmatools/sampledata/ASCII_solutions/LYLOUT_140526_*.dat.gz

Directories for the HDF5 files, grids, and plots are created within the 
directory indicated by /path/to/output/ 
"""

import sys, os, glob
from datetime import datetime, timedelta
import subprocess

from lmatools.flashsort.gen_autorun import logger_setup, sort_files
from lmatools.flashsort.gen_sklearn import DBSCANFlashSorter

from lmatools.grid.make_grids import grid_h5flashfiles, dlonlat_at_grid_center, write_cf_netcdf_latlon, write_cf_netcdf_3d_latlon
from lmatools.vis.multiples_nc import make_plot, make_plot_3d, read_file_3d
from six.moves import map


def tfromfile(name):
    parts = name.split('_')
    y, m, d = list(map(int, (parts[-3][0:2], parts[-3][2:4], parts[-3][4:6])))
    H, M, S = list(map(int, (parts[-2][0:2], parts[-2][2:4], parts[-2][4:6])))
    return y+2000,m,d,H,M,S    

def sort_flashes(files, base_sort_dir, params):
    """ Given a list of LMA ASCII data files, created HDF5 flash-sorted data 
        files in base_sort_dir/h5_files.
    
        params is a dictionary with the following format
        params = {'stations':(6,99), # range of allowable numbers of contributing stations
                  'chi2':(0,1.0),    # range of allowable chi-sq values
                  'distance':3000.0, 'thresh_critical_time':0.15, # space and time grouping thresholds
                  'thresh_duration':3.0, # maximum expected flash duration
                  'mask_length':6, # length of the hexadecimal station mask column in HDF5 
                  }
        
    """
    
    # base_sort_dir = outpath
    logger_setup(base_sort_dir)
    # files = lmatools.testing.test_gen_autorun_DBSCAN.get_sample_data_list()
    h5_dir = os.path.join(base_sort_dir, 'h5_files') #Changed for sensitivity analysis from 'h5_files'
    
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
    return h5_filenames

#####################################################    
#Change dx,dy,dz == 1.0e3 the original resolution. ##
#####################################################
def grid_and_plot(h5_filenames, base_sort_dir, dx=1.0e3, dy=1.0e3, dz=1.0e3, frame_interval=60.0,
                  x_bnd=(-200.0e3, 200.0e3), y_bnd=(-200.0e3, 200.0e3), z_bnd=(0.0e3, 20.0e3),
                  ctr_lat=33.5, ctr_lon=-101.5, center_ID='WTLMA',
                  n_cols=2, base_date=None
                  ):
    """ Given a list of HDF5 filenames (sorted by time order) in h5_filenames, 
        create 2D and 3D NetCDF grids with spacing dx, dy, dz in meters, 
        frame_interval in seconds, and tuples of grid edges 
        x_bnd, y_bnd, and z_bnd in meters
                  
        The actual grids are in regular lat,lon coordinates, with spacing at the 
        grid center matched to the dx, dy values given.
                                    
        n_cols controls how many columns are plotted on each page.
                  
        Grids and plots are written to base_sort_dir/grid_files/ and  base_sort_dir/plots/
				  
		base_date is used to optionally set a common reference time for each of the NetCDF grids.
    """              
    # not really in km, just a different name to distinguish from similar variables below.
    dx_km=dx
    dy_km=dy
    x_bnd_km = x_bnd
    y_bnd_km = y_bnd
    z_bnd_km = z_bnd

    grid_dir = os.path.join(base_sort_dir, 'grid_files')
    plot_dir = os.path.join(base_sort_dir, 'plots')

    # There are similar functions in lmatools to grid on a regular x,y grid in some map projection.
    dx, dy, x_bnd, y_bnd = dlonlat_at_grid_center(ctr_lat, ctr_lon, 
                                dx=dx_km, dy=dy_km,
                                x_bnd = x_bnd_km, y_bnd = y_bnd_km )
    # print("dx, dy = {0}, {1} deg".format(dx,dy))
    # print("lon_range = {0} deg".format(x_bnd))
    # print("lat_range = {0} deg".format(y_bnd))
    
    for f in h5_filenames:
        y,m,d,H,M,S = tfromfile(f)
        # print y,m,d,H,M,S
        start_time = datetime(y,m,d, H,M,S)
        end_time   = start_time + timedelta(0,600)
        date = start_time
        # print start_time, end_time
    
        outpath = grid_dir+'/20%s' %(date.strftime('%y/%b/%d'))
        if os.path.exists(outpath) == False:
            os.makedirs(outpath)
            subprocess.call(['chmod', 'a+w', outpath, grid_dir+'/20%s' %(date.strftime('%y/%b')), grid_dir+'/20%s' %(date.strftime('%y'))])
        if True:
            grid_h5flashfiles(h5_filenames, start_time, end_time, frame_interval=frame_interval, proj_name='latlong',
                    base_date = base_date,
					dx=dx, dy=dy, x_bnd=x_bnd, y_bnd=y_bnd, z_bnd=z_bnd_km,
                    ctr_lon=ctr_lon, ctr_lat=ctr_lat, outpath = outpath,
                    output_writer = write_cf_netcdf_latlon, output_writer_3d = write_cf_netcdf_3d_latlon,
                    output_filename_prefix=center_ID, spatial_scale_factor=1.0
                    )
        
    # Create plots
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
        make_plot(f, var, n_cols=n_cols, x_name='longitude', y_name='latitude', outpath = outpath)
    
    # for f in nc_names_3d:
    #     gridtype = f.split('dx_')[-1].replace('.nc', '').replace('_3d', '')
    #     var = mapping[gridtype]
    #     # grid_range = range_mapping[gridtype]
    #
    #     ###Read grid files, then plot in either 2d or 3d space###
    #     grid, grid_name, x, y, z, t, grid_t_idx, grid_x_idx, grid_z_idx = read_file_3d(f, var, x_name='longitude', y_name='latitude', z_name='altitude')
    #     make_plot_3d(grid, grid_name, x, y, z, t,
    #                  grid_t_idx, grid_x_idx, grid_z_idx,
    #                  n_cols = n_cols, outpath = outpath)
        #, grid_range=grid_range)
        
    return nc_names_2d, nc_names_3d
        
        
if __name__ == '__main__':
    params = {'stations':(6,99), # range of allowable numbers of contributing stations
              'chi2':(0,1.0),    # range of allowable chi-sq values
              'distance':3000.0, 'thresh_critical_time':0.15, # space and time grouping thresholds
              'thresh_duration':3.0, # maximum expected flash duration
              'ctr_lat':33.56, 'ctr_lon':-101.88, #center lat/lon to use for flash sorting, gridding
              'mask_length':6, # length of the hexadecimal station mask column in the LMA ASCII files
              }
    center_ID='WTLMA'
    
    data_out = sys.argv[1]
    filenames = sys.argv[2:]
    # h5_filenames = sort_flashes(filenames, data_out, params) #Turn off when .h5 files are made to make grid files!!
    h5_filenames=filenames
    # other keyword arguments control the grid spacing ... see the function definition

    nc_names_2d, nc_names_3d = grid_and_plot(h5_filenames, data_out, base_date=datetime(2012, 1, 1)
        ctr_lat=params['ctr_lat'], ctr_lon=params['ctr_lon'], center_ID=center_ID)


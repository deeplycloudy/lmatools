from __future__ import absolute_import
from __future__ import print_function
import glob 
import os, sys
from datetime import datetime, timedelta

from math import ceil

import numpy as np
import tables

from . import density_to_files
from lmatools.stream.subset import broadcast
from lmatools.io.LMA_h5_file import read_flashes, to_seconds

from lmatools.coordinateSystems import MapProjection, GeographicSystem
from six.moves import range
    

def dlonlat_at_grid_center(ctr_lat, ctr_lon, dx=4.0e3, dy=4.0e3,
    x_bnd = (-100e3, 100e3), y_bnd = (-100e3, 100e3),
    proj_datum = 'WGS84', proj_ellipse = 'WGS84'):
    """ 
    
    Utility function useful for producing a regular grid of lat/lon data, 
    where an approximate spacing (dx, dy) and total span of the grid (x_bnd, y_bnd) 
    is desired. Units are in meters.
    
    There is guaranteed to be distortion away from the grid center, i.e.,
    only the grid cells adjacent to the center location will have area dx * dy. 
    
    Likewise, the lat, lon range is calculated naively using dlat, dlon multiplied
    by the number of grid cells implied by x_bnd/dx, y_bnd/dy. This is the naive approach,
    but probably what's expected when specifying distances in kilometers for
    an inherently distorted lat/lon grid.
    
    Returns: 
    (dlon, dlat, lon_bnd, lat_bnd)
    corresponding to 
    (dx, dy, x_range, y_range)
    
    """

    # Use the Azimuthal equidistant projection as the method for converting to kilometers.
    proj_name = 'aeqd'
    
    mapProj = MapProjection(projection=proj_name, ctrLat=ctr_lat, ctrLon=ctr_lon, lat_ts=ctr_lat, 
                            lon_0=ctr_lon, lat_0=ctr_lat, lat_1=ctr_lat, ellipse=proj_ellipse, datum=proj_datum)
    geoProj = GeographicSystem()
    
    # Get dlat
    lon_n, lat_n, z_n = geoProj.fromECEF(*mapProj.toECEF(0,dy,0))
    dlat = lat_n - ctr_lat
    
    # Get dlon
    lon_e, lat_e, z_e = geoProj.fromECEF(*mapProj.toECEF(dx,0,0)) 
    dlon = lon_e - ctr_lon
    
    lon_min = ctr_lon + dlon * (x_bnd[0]/dx)
    lon_max = ctr_lon + dlon * (x_bnd[1]/dx)
    lat_min = ctr_lat + dlat * (y_bnd[0]/dy)
    lat_max = ctr_lat + dlat * (y_bnd[1]/dy)
    
    # Alternate method: lat lon for the actual distance to the NSEW in the projection
    #lon_range_n, lat_range_n, z_range_n = geoProj.fromECEF(*mapProj.toECEF(0,y_bnd,0))
    #lon_range_e, lat_range_e, z_range_e = geoProj.fromECEF(*mapProj.toECEF(x_bnd,0,0))
    
    return dlon, dlat, (lon_min, lon_max), (lat_min, lat_max)
    

def write_cf_netcdf_latlon(outfile, t_start, t, xloc, yloc, lon_for_x, lat_for_y, ctr_lat, ctr_lon, grid, grid_var_name, grid_description, format='i', **kwargs):
    """ Write a Climate and Forecast Metadata-compliant NetCDF file. 
        
        Grid is regular in lon, lat and so no map projection information is necessary.
    
        Should display natively in conformant packages like McIDAS-V.
        
    """

    import scipy.io.netcdf as nc

    missing_value = -9999
    
    nc_out = nc.NetCDFFile(outfile, 'w')
    nc_out.createDimension('lon', xloc.shape[0])
    nc_out.createDimension('lat', yloc.shape[0])
    nc_out.createDimension('ntimes', t.shape[0])  #unlimited==None

    # declare the coordinate reference system, WGS84 values
    proj = nc_out.createVariable('crs', 'i', ())
    proj.grid_mapping_name = 'latitude_longitude'
    proj.longitude_of_prime_meridian = 0.0 
    proj.semi_major_axis = 6378137.0 
    proj.inverse_flattening = 298.257223563 
    
    y_coord = nc_out.createVariable('latitude', 'f', ('lat',))
    y_coord.units = "degrees_north"
    y_coord.long_name = "latitude"
    y_coord.standard_name = 'latitude'

    x_coord = nc_out.createVariable('longitude', 'f', ('lon',))
    x_coord.units = "degrees_east"
    x_coord.long_name = "longitude"
    x_coord.standard_name = 'longitude'

    times = nc_out.createVariable('time', 'f', ('ntimes',) )#, filters=no_compress)
    times.long_name="time"
    times.units = "seconds since %s" % t_start.strftime('%Y-%m-%d %H:%M:%S')
    
    # lons = nc_out.createVariable('lons', 'd', ('nx','ny') )#, filters=no_compress)
    # lons.long_name="longitude"
    # lons.standard_name="longitude"
    # lons.units = "degrees_east"
    # 
    # lats = nc_out.createVariable('lats', 'd', ('nx','ny') )#, filters=no_compress)
    # lats.long_name="latitude"
    # lats.standard_name="latitude"
    # lats.units = "degrees_north"

    lightning2d = nc_out.createVariable(grid_var_name, format, ('ntimes','lon','lat') )#, filters=no_compress)
    lightning2d.long_name=grid_description #'LMA VHF event counts (vertically integrated)'
    lightning2d.units=kwargs['grid_units']
    # lightning2d.coordinates='time lons lats'
    lightning2d.grid_mapping = "crs"
    lightning2d.missing_value = missing_value

    x_coord[:] = xloc[:]
    y_coord[:] = yloc[:]
    times[:] = t[:]
    # lons[:] = lon_for_x[:]
    # lats[:] = lat_for_y[:]

    for i in range(grid.shape[2]):
        lightning2d[i,:,:] = grid[:,:,i]
    nc_out.close()    

def write_cf_netcdf(outfile, t_start, t, xloc, yloc, lon_for_x, lat_for_y, ctr_lat, ctr_lon, grid, grid_var_name, grid_description, format='i', **kwargs):
    """ Write a Climate and Forecast Metadata-compliant NetCDF file. 
    
        Should display natively in conformant packages like McIDAS-V.
        
    """

    import scipy.io.netcdf as nc

    missing_value = -9999
    
    nc_out = nc.NetCDFFile(outfile, 'w')
    nc_out.createDimension('nx', xloc.shape[0])
    nc_out.createDimension('ny', yloc.shape[0])
    nc_out.createDimension('ntimes', t.shape[0])  #unlimited==None

    proj = nc_out.createVariable('Lambert_Azimuthal_Equal_Area', 'i', ())
    proj.grid_mapping_name = 'lambert_azimuthal_equal_area'
    proj.longitude_of_projection_origin = ctr_lon
    proj.latitude_of_projection_origin = ctr_lat
    proj.false_easting = 0.0
    proj.false_northing = 0.0
    
    # x_coord = nc_out.createVariable('longitude', 'f', ('nx',))
    # x_coord.long_name="longitude"
    # x_coord.standard_name="longitude"
    # x_coord.units = "degrees_east"
    x_coord = nc_out.createVariable('x', 'f', ('nx',))
    x_coord.units = "km"
    x_coord.long_name = "x coordinate of projection"
    x_coord.standard_name = 'projection_x_coordinate'

    # y_coord = nc_out.createVariable('latitude', 'f', ('nx',))
    # y_coord.long_name="latitude"
    # y_coord.standard_name="latitude"
    # y_coord.units = "degrees_north"
    y_coord = nc_out.createVariable('y', 'f', ('ny',))
    y_coord.units = "km"
    y_coord.long_name = "y coordinate of projection"
    y_coord.standard_name = 'projection_y_coordinate'

    times = nc_out.createVariable('time', 'f', ('ntimes',) )#, filters=no_compress)
    times.long_name="time"
    times.units = "seconds since %s" % t_start.strftime('%Y-%m-%d %H:%M:%S')
    
    #Dtype change from 'd' to 'f':
    lons = nc_out.createVariable('lons', 'f', ('nx','ny') )#, filters=no_compress)
    lons.long_name="longitude"
    lons.standard_name="longitude"
    lons.units = "degrees_east"
    
    lats = nc_out.createVariable('lats', 'f', ('nx','ny') )#, filters=no_compress)
    lats.long_name="latitude"
    lats.standard_name="latitude"
    lats.units = "degrees_north"

    lightning2d = nc_out.createVariable(grid_var_name, format, ('ntimes','nx','ny') )#, filters=no_compress)
    lightning2d.long_name=grid_description #'LMA VHF event counts (vertically integrated)'
    lightning2d.units='dimensionless'
    lightning2d.coordinates='time lons lats'
    lightning2d.grid_mapping = "Lambert_Azimuthal_Equal_Area"
    lightning2d.missing_value = missing_value

    x_coord[:] = xloc[:]
    y_coord[:] = yloc[:]
    times[:] = t[:]
    lons[:] = lon_for_x[:]
    lats[:] = lat_for_y[:]

    for i in range(grid.shape[2]):
        lightning2d[i,:,:] = grid[:,:,i]
    nc_out.close()

def write_cf_netcdf_3d(outfile, t_start, t, xloc, yloc, zloc, lon_for_x, lat_for_y, alt_for_z, ctr_lat, ctr_lon, ctr_alt, grid, grid_var_name, grid_description, format='i', **kwargs):
    """ Write a Climate and Forecast Metadata-compliant NetCDF file. 
    
        Should display natively in conformant packages like McIDAS-V.
        
    """

    import scipy.io.netcdf as nc

    missing_value = -9999
    
    nc_out = nc.NetCDFFile(outfile, 'w')
    nc_out.createDimension('nx', xloc.shape[0])
    nc_out.createDimension('ny', yloc.shape[0])
    nc_out.createDimension('nz', zloc.shape[0])
    nc_out.createDimension('ntimes', t.shape[0])  #unlimited==None

    proj = nc_out.createVariable('Lambert_Azimuthal_Equal_Area', 'i', ())
    proj.grid_mapping_name = 'lambert_azimuthal_equal_area'
    proj.longitude_of_projection_origin = ctr_lon
    proj.latitude_of_projection_origin = ctr_lat
    proj.altitude_of_projection_origin = ctr_alt
    proj.false_easting = 0.0
    proj.false_northing = 0.0
    
    # x_coord = nc_out.createVariable('longitude', 'f', ('nx',))
    # x_coord.long_name="longitude"
    # x_coord.standard_name="longitude"
    # x_coord.units = "degrees_east"
    x_coord = nc_out.createVariable('x', 'f', ('nx',))
    x_coord.units = "km"
    x_coord.long_name = "x coordinate of projection"
    x_coord.standard_name = 'projection_x_coordinate'

    # y_coord = nc_out.createVariable('latitude', 'f', ('nx',))
    # y_coord.long_name="latitude"
    # y_coord.standard_name="latitude"
    # y_coord.units = "degrees_north"
    y_coord = nc_out.createVariable('y', 'f', ('ny',))
    y_coord.units = "km"
    y_coord.long_name = "y coordinate of projection"
    y_coord.standard_name = 'projection_y_coordinate'

    z_coord = nc_out.createVariable('z', 'f', ('nz',))
    z_coord.units = 'km'
    z_coord.long_name = "z coordinate of projection"
    z_coord.standard_name = 'projection_z_coordinate'

    times = nc_out.createVariable('time', 'f', ('ntimes',) )#, filters=no_compress)
    times.long_name="time"
    times.units = "seconds since %s" % t_start.strftime('%Y-%m-%d %H:%M:%S')
    
    #Dtype changed due to previous issue with pupynere.
    #---------------------------------------------------
    lons = nc_out.createVariable('lons', 'f', ('nx','ny','nz') )#, filters=no_compress)
    lons.long_name="longitude"
    lons.standard_name="longitude"
    lons.units = "degrees_east"
    
    lats = nc_out.createVariable('lats', 'f', ('nx','ny', 'nz') )#, filters=no_compress)
    lats.long_name="latitude"
    lats.standard_name="latitude"
    lats.units = "degrees_north"
    
    alts = nc_out.createVariable('alts', 'f', ('nx','ny', 'nz') )#, filters=no_compress)
    alts.long_name="altitude"
    alts.standard_name="altitude"
    alts.units = "meters"

    lightning3d = nc_out.createVariable(grid_var_name, format, ('ntimes','nx','ny', 'nz') )#, filters=no_compress)
    lightning3d.long_name=grid_description #'LMA VHF event counts (vertically integrated)'
    lightning3d.units='dimensionless'
    lightning3d.coordinates='time lons lats alts'
    lightning3d.grid_mapping = "Lambert_Azimuthal_Equal_Area"
    lightning3d.missing_value = missing_value

    x_coord[:] = xloc[:]
    y_coord[:] = yloc[:]
    z_coord[:] = zloc[:]
    times[:] = t[:]
    lons[:] = lon_for_x[:]
    lats[:] = lat_for_y[:]
    alts[:] = alt_for_z[:]

    for i in range(grid.shape[3]):
        lightning3d[i,:,:,:] = grid[:,:,:,i]
    nc_out.close()

def write_cf_netcdf_3d_latlon(outfile, t_start, t, xloc, yloc, zloc, lon_for_x, lat_for_y, alt_for_z, ctr_lat, ctr_lon, ctr_alt, grid, grid_var_name, grid_description, format='i', **kwargs):
    """ Write a Climate and Forecast Metadata-compliant NetCDF file. 
    
        Should display natively in conformant packages like McIDAS-V.
        
    """

    import scipy.io.netcdf as nc

    missing_value = -9999
    
    nc_out = nc.NetCDFFile(outfile, 'w')
    nc_out.createDimension('lon', xloc.shape[0])
    nc_out.createDimension('lat', yloc.shape[0])
    nc_out.createDimension('alt', zloc.shape[0])
    nc_out.createDimension('ntimes', t.shape[0])  #unlimited==None

    # declare the coordinate reference system, WGS84 values
    proj = nc_out.createVariable('crs', 'i', ())
    proj.grid_mapping_name = 'latitude_longitude'
    proj.longitude_of_prime_meridian = 0.0 
    proj.semi_major_axis = 6378137.0 
    proj.inverse_flattening = 298.257223563 
    
    y_coord = nc_out.createVariable('latitude', 'f', ('lat',))
    y_coord.units = "degrees_north"
    y_coord.long_name = "latitude"
    y_coord.standard_name = 'latitude'

    x_coord = nc_out.createVariable('longitude', 'f', ('lon',))
    x_coord.units = "degrees_east"
    x_coord.long_name = "longitude"
    x_coord.standard_name = 'longitude'

    z_coord = nc_out.createVariable('altitude', 'f', ('alt',))
    z_coord.units = 'km'
    z_coord.long_name = "height above mean sea level"
    z_coord.standard_name = 'altitude'
    z_coord.positive = 'up'

    times = nc_out.createVariable('time', 'f', ('ntimes',) )#, filters=no_compress)
    times.long_name="time"
    times.units = "seconds since %s" % t_start.strftime('%Y-%m-%d %H:%M:%S')
    
    lightning3d = nc_out.createVariable(grid_var_name, format, ('ntimes','lon','lat','alt') )#, filters=no_compress)
    lightning3d.long_name=grid_description #'LMA VHF event counts (vertically integrated)'
    lightning3d.units=kwargs['grid_units']
    # lightning3d.coordinates='time lons lats alts'
    lightning3d.grid_mapping = "crs"
    lightning3d.missing_value = missing_value

    x_coord[:] = xloc[:]
    y_coord[:] = yloc[:]
    z_coord[:] = zloc[:]
    times[:] = t[:]
    # lons[:] = lon_for_x[:]
    # lats[:] = lat_for_y[:]
    # alts[:] = alt_for_z[:]

    for i in range(grid.shape[3]):
        lightning3d[i,:,:,:] = grid[:,:,:,i]
    nc_out.close()


def time_edges(start_time, end_time, frame_interval):
    """ Return lists cooresponding the start and end times of frames lasting frame_interval
        between start_time and end_time. The last interval may extend beyond end_time, but 
        by no more than frame_interval. This makes each frame the same length.
        
        returns t_edges, duration, where t_edges is a list of datetime objects, and
        duration is the total duration between the start and end times (and not the duration
        of all frames)
    """
    frame_dt = timedelta(0, frame_interval, 0)
    duration = end_time - start_time
    n_frames = int(np.ceil(to_seconds(duration) / to_seconds(frame_dt)))
    t_edges = [start_time + i*frame_dt for i in range(n_frames+1)]
    return t_edges, duration
    
def seconds_since_start_of_day(start_time, t):
    """ For each datetime object t, return the number of seconds elapsed since the 
        start of the date given by start_time. Only the date part of start_time is used. 
    """
    ref_date = start_time.date()
    t_ref = datetime(ref_date.year, ref_date.month, ref_date.day)
    t_edges_seconds = [to_seconds(edge - t_ref) for edge in t]
    return t_ref, t_edges_seconds


class FlashGridder(object):
    def __init__(self, start_time, end_time, 
                    frame_interval=120.0, dx=4.0e3, dy=4.0e3, dz=1.0e3,
					base_date = None,
                    x_bnd = (-100e3, 100e3),
                    y_bnd = (-100e3, 100e3),
                    z_bnd = (0e3, 20e3),
                    ctr_lat = 35.23833, ctr_lon = -97.46028, ctr_alt=0.0,
                    proj_name = 'aeqd',
                    proj_datum = 'WGS84',
                    proj_ellipse = 'WGS84',
                    flash_count_logfile = None,):
        """ Class to support gridding of flash and event data. 
                    
            On init, specify the grid """
            
        # args, kwargs that are saved for the future
        self.start_time = start_time
        self.dx, self.dy, self.dz = dx, dy, dz
        self.end_time = end_time
        self.frame_interval = frame_interval
        self.base_date = base_date
        self.min_points_per_flash = 1
        self.proj_name = 'aeqd'
        self.ctr_lat, self.ctr_lon, self.ctr_alt = ctr_lat, ctr_lon, ctr_alt
        
        if flash_count_logfile is None:
            flash_count_logfile = sys.stdout
    
        t_edges, duration = time_edges(self.start_time, self.end_time, self.frame_interval)
        # reference time is the date part of the start_time, unless the user provides a different date.
        if self.base_date is None:
            t_ref, t_edges_seconds = seconds_since_start_of_day(self.start_time, t_edges)
        else:
            t_ref, t_edges_seconds = seconds_since_start_of_day(self.base_date, t_edges)
        n_frames = len(t_edges)-1
    
        xedge=np.arange(x_bnd[0], x_bnd[1]+dx, dx)
        yedge=np.arange(y_bnd[0], y_bnd[1]+dy, dy)
        zedge=np.arange(z_bnd[0], z_bnd[1]+dz, dz) 
    
        x0 = xedge[0]
        y0 = yedge[0]
        z0 = zedge[0]
    
        if proj_name == 'latlong':
            dx_units = '{0:6.4f}deg'.format(dx)
            mapProj = GeographicSystem()
        else:
            dx_units = '{0:5.1f}m'.format(dx)
            mapProj = MapProjection(projection=proj_name, ctrLat=ctr_lat, ctrLon=ctr_lon, lat_ts=ctr_lat, 
                                lon_0=ctr_lon, lat_0=ctr_lat, lat_1=ctr_lat, ellipse=proj_ellipse, datum=proj_datum)
        geoProj = GeographicSystem()
        
        self.t_ref = t_ref
        self.t_edges_seconds = t_edges_seconds
        self.duration = duration
        self.xedge = xedge
        self.yedge = yedge
        self.zedge = zedge
        self.mapProj = mapProj
        self.geoProj = geoProj
        self.dx_units = dx_units
        
        event_density_grid  = np.zeros((xedge.shape[0]-1, yedge.shape[0]-1, n_frames), dtype='int32')
        init_density_grid   = np.zeros((xedge.shape[0]-1, yedge.shape[0]-1, n_frames), dtype='int32')
        extent_density_grid = np.zeros((xedge.shape[0]-1, yedge.shape[0]-1, n_frames), dtype='int32')
        footprint_grid      = np.zeros((xedge.shape[0]-1, yedge.shape[0]-1, n_frames), dtype='float32')

        specific_energy_grid = np.zeros((xedge.shape[0]-1, yedge.shape[0]-1, n_frames), dtype='float32')
        total_energy_grid = np.zeros((xedge.shape[0]-1, yedge.shape[0]-1, n_frames), dtype='float32')
        flashsize_std_grid  = np.zeros((xedge.shape[0]-1, yedge.shape[0]-1, n_frames), dtype='float32')

        event_density_grid_3d  = np.zeros((xedge.shape[0]-1, yedge.shape[0]-1, zedge.shape[0]-1, n_frames), dtype='int32')
        init_density_grid_3d   = np.zeros((xedge.shape[0]-1, yedge.shape[0]-1, zedge.shape[0]-1, n_frames), dtype='int32')
        extent_density_grid_3d = np.zeros((xedge.shape[0]-1, yedge.shape[0]-1, zedge.shape[0]-1, n_frames), dtype='int32')
        footprint_grid_3d      = np.zeros((xedge.shape[0]-1, yedge.shape[0]-1, zedge.shape[0]-1, n_frames), dtype='float32')
    
        specific_energy_grid_3d = np.zeros((xedge.shape[0]-1, yedge.shape[0]-1, zedge.shape[0]-1, n_frames), dtype='float32')
        total_energy_grid_3d    = np.zeros((xedge.shape[0]-1, yedge.shape[0]-1, zedge.shape[0]-1, n_frames), dtype='float32')
        flashsize_std_grid_3d   = np.zeros((xedge.shape[0]-1, yedge.shape[0]-1, zedge.shape[0]-1, n_frames), dtype='float32')
        
        self.outgrids = (extent_density_grid, 
                    init_density_grid,   
                    event_density_grid,  
                    footprint_grid,
                    specific_energy_grid,
                    flashsize_std_grid,
                    total_energy_grid,
                    )
                
        self.outgrids_3d = (extent_density_grid_3d,
                    init_density_grid_3d,
                    event_density_grid_3d,
                    footprint_grid_3d,
                    specific_energy_grid_3d,
                    flashsize_std_grid_3d,
                    total_energy_grid_3d
                    )
        
        
        all_frames = []
        for i in range(n_frames):
            extent_out = {'name':'extent'}
            init_out   = {'name':'init'}
            event_out  = {'name':'event'}
            std_out    = {'name':'std'}
        
            extent_out_3d = {'name':'extent_3d'}
            init_out_3d   = {'name':'init_3d'}
            event_out_3d  = {'name':'event_3d'}
            std_out_3d    = {'name':'std_3d'}
        
            accum_event_density  = density_to_files.accumulate_points_on_grid(event_density_grid[:,:,i], xedge, yedge,  out=event_out, label='event')
            accum_init_density   = density_to_files.accumulate_points_on_grid(init_density_grid[:,:,i], xedge, yedge,   out=init_out,  label='init')
            accum_extent_density = density_to_files.accumulate_points_on_grid(extent_density_grid[:,:,i], xedge, yedge, out=extent_out,label='extent')
            accum_footprint      = density_to_files.accumulate_points_on_grid(footprint_grid[:,:,i], xedge, yedge, label='footprint')

            accum_specific_energy = density_to_files.accumulate_energy_on_grid(specific_energy_grid[:,:,i], xedge, yedge, out=extent_out, label='specific_energy')
            accum_flashstd       = density_to_files.accumulate_points_on_grid_sdev(flashsize_std_grid[:,:,i], footprint_grid[:,:,i], xedge, yedge, out=extent_out, label='flashsize_std')
            accum_total_energy   = density_to_files.accumulate_energy_on_grid(total_energy_grid[:,:,i], xedge, yedge, out=extent_out, label='total_energy')

            accum_event_density_3d  = density_to_files.accumulate_points_on_grid_3d(event_density_grid_3d[:,:,:,i], xedge, yedge, zedge,  out=event_out_3d, label='event_3d')
            accum_init_density_3d   = density_to_files.accumulate_points_on_grid_3d(init_density_grid_3d[:,:,:,i], xedge, yedge, zedge,   out=init_out_3d,  label='init_3d')
            accum_extent_density_3d = density_to_files.accumulate_points_on_grid_3d(extent_density_grid_3d[:,:,:,i], xedge, yedge, zedge, out=extent_out_3d,label='extent_3d')
            accum_footprint_3d      = density_to_files.accumulate_points_on_grid_3d(footprint_grid_3d[:,:,:,i], xedge, yedge, zedge, label='footprint_3d')

            accum_specific_energy_3d = density_to_files.accumulate_energy_on_grid_3d(specific_energy_grid_3d[:,:,:,i], xedge, yedge, zedge, out=extent_out_3d,label='specific_energy_3d')
            accum_flashstd_3d       = density_to_files.accumulate_points_on_grid_sdev_3d(flashsize_std_grid_3d[:,:,:,i], footprint_grid_3d[:,:,:,i], xedge, yedge, zedge, out=extent_out_3d,label='flashsize_std_3d')
            accum_total_energy_3d   = density_to_files.accumulate_energy_on_grid_3d(total_energy_grid_3d[:,:,:,i], xedge, yedge, zedge, out=extent_out_3d,label='total_energy_3d')

            extent_out['func'] = accum_extent_density
            init_out['func'] = accum_init_density
            event_out['func'] = accum_event_density

            extent_out_3d['func'] = accum_extent_density_3d
            init_out_3d['func'] = accum_init_density_3d
            event_out_3d['func'] = accum_event_density_3d        
        
            event_density_target  = density_to_files.point_density(accum_event_density)
            init_density_target   = density_to_files.point_density(accum_init_density)
            extent_density_target = density_to_files.extent_density(x0, y0, dx, dy, accum_extent_density)
            mean_footprint_target = density_to_files.extent_density(x0, y0, dx, dy, accum_footprint, weight_key='area')
            mean_energy_target    = density_to_files.extent_density(x0, y0, dx, dy, accum_specific_energy, weight_key='specific_energy') #tot_energy
            mean_total_energy_target = density_to_files.extent_density(x0, y0, dx, dy, accum_total_energy, weight_key='total_energy')    #Energy
            std_flashsize_target  = density_to_files.extent_density(x0, y0, dx, dy, accum_flashstd, weight_key='area')

            event_density_target_3d  = density_to_files.point_density_3d(accum_event_density_3d)
            init_density_target_3d   = density_to_files.point_density_3d(accum_init_density_3d)
            extent_density_target_3d = density_to_files.extent_density_3d(x0, y0, z0, dx, dy, dz, accum_extent_density_3d)
            mean_footprint_target_3d = density_to_files.extent_density_3d(x0, y0, z0, dx, dy, dz, accum_footprint_3d, weight_key='area')
            
            mean_energy_target_3d    = density_to_files.extent_density_3d(x0, y0, z0, dx, dy, dz, accum_specific_energy_3d, weight_key='specific_energy') #tot_energy
            mean_total_energy_target_3d = density_to_files.extent_density_3d(x0, y0, z0, dx, dy, dz, accum_total_energy_3d, weight_key='total_energy')     #Energy
            std_flashsize_target_3d  = density_to_files.extent_density_3d(x0, y0, z0, dx, dy, dz, accum_flashstd_3d, weight_key='area')

            spew_to_density_types = broadcast( ( 
                        density_to_files.project('lon', 'lat', 'alt', mapProj, geoProj, event_density_target, use_flashes=False),
                        density_to_files.project('init_lon', 'init_lat', 'init_alt', mapProj, geoProj, init_density_target, use_flashes=True),
                        density_to_files.project('lon', 'lat', 'alt', mapProj, geoProj, extent_density_target, use_flashes=False),
                        density_to_files.project('lon', 'lat', 'alt', mapProj, geoProj, mean_footprint_target, use_flashes=False),
                        density_to_files.project('lon', 'lat', 'alt', mapProj, geoProj, mean_energy_target, use_flashes=False),
                        density_to_files.project('lon', 'lat', 'alt', mapProj, geoProj, std_flashsize_target, use_flashes=False),
                        density_to_files.project('lon', 'lat', 'alt', mapProj, geoProj, mean_total_energy_target, use_flashes=False),
                    
                        density_to_files.project('lon', 'lat', 'alt', mapProj, geoProj, event_density_target_3d, use_flashes=False),
                        density_to_files.project('init_lon', 'init_lat', 'init_alt', mapProj, geoProj, init_density_target_3d, use_flashes=True),
                        density_to_files.project('lon', 'lat', 'alt', mapProj, geoProj, extent_density_target_3d, use_flashes=False),
                        density_to_files.project('lon', 'lat', 'alt', mapProj, geoProj, mean_footprint_target_3d, use_flashes=False),
                        density_to_files.project('lon', 'lat', 'alt', mapProj, geoProj, mean_energy_target_3d, use_flashes=False),
                        density_to_files.project('lon', 'lat', 'alt', mapProj, geoProj, std_flashsize_target_3d, use_flashes=False),
                        density_to_files.project('lon', 'lat', 'alt', mapProj, geoProj, mean_total_energy_target_3d, use_flashes=False),
                        ) )

            all_frames.append( density_to_files.extract_events_for_flashes( spew_to_density_types ) )

        frame_count_log = density_to_files.flash_count_log(flash_count_logfile)
        
        framer = density_to_files.flashes_to_frames(t_edges_seconds, all_frames, time_key='start', time_edges_datetime=t_edges, flash_counter=frame_count_log)
    
        self.framer=framer
        
    def process_flashes(self, h5_filenames, min_points_per_flash=10):
        self.min_points_per_flash = min_points_per_flash
        read_flashes( h5_filenames, self.framer, base_date=self.t_ref, min_points=min_points_per_flash)

    def write_grids(self, outpath = '',
                    output_writer = write_cf_netcdf, 
                    output_writer_3d = write_cf_netcdf_3d,
                    output_filename_prefix="LMA",
                    output_kwargs = {},
                    spatial_scale_factor = 1.0/1000.0,):
        xedge = self.xedge    
        yedge = self.yedge
        zedge = self.zedge
        t_ref = self.t_ref
        t_edges_seconds = self.t_edges_seconds
        mapProj = self.mapProj
        geoProj = self.geoProj
        ctr_lat, ctr_lon, ctr_alt = self.ctr_lat, self.ctr_lon, self.ctr_alt
        dx, dy, dz = self.dx, self.dy, self.dz
        outgrids = self.outgrids
        outgrids_3d = self.outgrids_3d
                
        x_coord = (xedge[:-1] + xedge[1:])/2.0
        y_coord = (yedge[:-1] + yedge[1:])/2.0
        z_coord = (zedge[:-1] + zedge[1:])/2.0
        nx = x_coord.shape[0]
        ny = y_coord.shape[0]
        nz = z_coord.shape[0]
    
        x_all, y_all = (a.T for a in np.meshgrid(x_coord, y_coord))
        assert x_all.shape == y_all.shape
        assert x_all.shape[0] == nx
        assert x_all.shape[1] == ny
        z_all = np.zeros_like(x_all)
    
    
        grid_shape_3d = (nx,ny,nz) 
        x_ones_3d = np.ones(grid_shape_3d, dtype='f4')
        y_ones_3d = np.ones(grid_shape_3d, dtype='f4')
        z_ones_3d = np.ones(grid_shape_3d, dtype='f4')
    
        x_all_3d = x_coord[:, None, None]*x_ones_3d
        y_all_3d = y_coord[None,:,None]*y_ones_3d
        z_all_3d = z_coord[None, None, :]*z_ones_3d
    
            
        lons, lats, alts = x,y,z = geoProj.fromECEF( *mapProj.toECEF(x_all, y_all, z_all) )
        lons.shape=x_all.shape
        lats.shape=y_all.shape
    
        lons_3d, lats_3d, alts_3d = x_3d,y_3d,z_3d = geoProj.fromECEF( *mapProj.toECEF(x_all_3d, y_all_3d, z_all_3d) )
        lons_3d.shape=x_all_3d.shape
        lats_3d.shape=y_all_3d.shape
        alts_3d.shape=z_all_3d.shape
    
    
        outflile_basename = os.path.join(outpath,'%s_%s_%d_%dsrc_%s-dx_' % (
                                                        output_filename_prefix, 
                                                        self.start_time.strftime('%Y%m%d_%H%M%S'), 
                                                        to_seconds(self.duration), 
                                                        self.min_points_per_flash, 
                                                        self.dx_units )
                                                        )
    
        outfiles = (outflile_basename+'flash_extent.nc',
                    outflile_basename+'flash_init.nc',
                    outflile_basename+'source.nc',
                    outflile_basename+'footprint.nc',
                    outflile_basename+'specific_energy.nc',
                    outflile_basename+'flashsize_std.nc',
                    outflile_basename+'total_energy.nc'
                    )
        outfiles_3d = (outflile_basename+'flash_extent_3d.nc',
                    outflile_basename+'flash_init_3d.nc',
                    outflile_basename+'source_3d.nc',
                    outflile_basename+'footprint_3d.nc',
                    outflile_basename+'specific_energy_3d.nc',
                    outflile_basename+'flashsize_std_3d.nc',
                    outflile_basename+'total_energy_3d.nc'
                    )
                                
                            
        field_names = ('flash_extent', 'flash_initiation', 'lma_source', 'flash_footprint', 'specific_energy', 'flashsize_std','total_energy')
    
        field_descriptions = ('LMA flash extent density',
                            'LMA flash initiation density',
                            'LMA source density',
                            'LMA local mean flash area',
                            'LMA flash specific energy (approx)',
                            'LMA local standard deviation of flash size',
                            'LMA flash total energy (approx)')
    
        if self.proj_name=='latlong':
            density_units = "grid"
            density_units_3d = "grid"
        else:
            density_units = "{0:5.1f} km^2".format(dx*spatial_scale_factor * dy*spatial_scale_factor).lstrip()
            density_units_3d = "{0:5.1f} km^3".format(dx*spatial_scale_factor * dy*spatial_scale_factor * dz*spatial_scale_factor).lstrip()
        time_units = "{0:5.1f} min".format(self.frame_interval/60.0).lstrip()
        density_label = 'Count per ' + density_units + " pixel per "+ time_units
        density_label_3d = 'Count per ' + density_units_3d + " pixel per "+ time_units
    
        field_units = ( density_label,
                        density_label,
                        density_label,
                        "km^2 per flash",
                        "J/kg per flash",
                        density_label,
                        "J per flash",
                         )
        field_units_3d = ( density_label_3d,
                        density_label_3d,
                        density_label_3d,
                        "km^2 per flash",
                        "J/kg per flash",
                        density_label_3d,
                        "J per flash",
                         )
    
        output_writer(outfiles[0], t_ref, np.asarray(t_edges_seconds[:-1]),
                        x_coord*spatial_scale_factor, y_coord*spatial_scale_factor, 
                        lons, lats, ctr_lat, ctr_lon, 
                        outgrids[0], field_names[0], field_descriptions[0], 
                        grid_units=field_units[0],
                        **output_kwargs)
        output_writer(outfiles[1], t_ref, np.asarray(t_edges_seconds[:-1]),
                        x_coord*spatial_scale_factor, y_coord*spatial_scale_factor, 
                        lons, lats, ctr_lat, ctr_lon, 
                        outgrids[1], field_names[1], field_descriptions[1], 
                        grid_units=field_units[1],
                        **output_kwargs)
        output_writer(outfiles[2], t_ref, np.asarray(t_edges_seconds[:-1]),
                        x_coord*spatial_scale_factor, y_coord*spatial_scale_factor, 
                        lons, lats, ctr_lat, ctr_lon, 
                        outgrids[2], field_names[2], field_descriptions[2], 
                        grid_units=field_units[2],
                        **output_kwargs)
        output_writer(outfiles[3], t_ref, np.asarray(t_edges_seconds[:-1]),
                        x_coord*spatial_scale_factor, y_coord*spatial_scale_factor, 
                        lons, lats, ctr_lat, ctr_lon, 
                        outgrids[3], field_names[3], field_descriptions[3], format='f', 
                        grid_units=field_units[3],
                        **output_kwargs)
        output_writer(outfiles[4], t_ref, np.asarray(t_edges_seconds[:-1]),
                        x_coord*spatial_scale_factor, y_coord*spatial_scale_factor, 
                        lons, lats, ctr_lat, ctr_lon, 
                        outgrids[4], field_names[4], field_descriptions[4], format='f', 
                        grid_units=field_units[4],
                        **output_kwargs)
        output_writer(outfiles[5], t_ref, np.asarray(t_edges_seconds[:-1]),
                        x_coord*spatial_scale_factor, y_coord*spatial_scale_factor, 
                        lons, lats, ctr_lat, ctr_lon, 
                        outgrids[5], field_names[5], field_descriptions[5], format='f', 
                        grid_units=field_units[5],
                        **output_kwargs)
        output_writer(outfiles[6], t_ref, np.asarray(t_edges_seconds[:-1]),
                        x_coord*spatial_scale_factor, y_coord*spatial_scale_factor, 
                        lons, lats, ctr_lat, ctr_lon, 
                        outgrids[6], field_names[6], field_descriptions[6], format='f', 
                        grid_units=field_units[6],
                        **output_kwargs)
    
        ########3D:
        output_writer_3d(outfiles_3d[0], t_ref, np.asarray(t_edges_seconds[:-1]),
                        x_coord*spatial_scale_factor, y_coord*spatial_scale_factor,
                        z_coord*spatial_scale_factor, 
                        lons_3d, lats_3d, alts_3d, ctr_lat, ctr_lon, ctr_alt,
                        outgrids_3d[0], field_names[0], field_descriptions[0], 
                        grid_units=field_units_3d[0],
                        **output_kwargs)
        output_writer_3d(outfiles_3d[1], t_ref, np.asarray(t_edges_seconds[:-1]),
                        x_coord*spatial_scale_factor, y_coord*spatial_scale_factor,
                        z_coord*spatial_scale_factor, 
                        lons_3d, lats_3d, alts_3d, ctr_lat, ctr_lon, ctr_alt,
                        outgrids_3d[1], field_names[1], field_descriptions[1], 
                        grid_units=field_units_3d[1],
                        **output_kwargs)
        output_writer_3d(outfiles_3d[2], t_ref, np.asarray(t_edges_seconds[:-1]),
                        x_coord*spatial_scale_factor, y_coord*spatial_scale_factor,
                        z_coord*spatial_scale_factor, 
                        lons_3d, lats_3d, alts_3d, ctr_lat, ctr_lon, ctr_alt,
                        outgrids_3d[2], field_names[2], field_descriptions[2], 
                        grid_units=field_units_3d[2],
                        **output_kwargs)
        output_writer_3d(outfiles_3d[3], t_ref, np.asarray(t_edges_seconds[:-1]),
                        x_coord*spatial_scale_factor, y_coord*spatial_scale_factor,
                        z_coord*spatial_scale_factor, 
                        lons_3d, lats_3d, alts_3d, ctr_lat, ctr_lon, ctr_alt,
                        outgrids_3d[3], field_names[3], field_descriptions[3], format='f', 
                        grid_units=field_units_3d[3],
                        **output_kwargs)
        output_writer_3d(outfiles_3d[4], t_ref, np.asarray(t_edges_seconds[:-1]),
                        x_coord*spatial_scale_factor, y_coord*spatial_scale_factor,
                        z_coord*spatial_scale_factor, 
                        lons_3d, lats_3d, alts_3d, ctr_lat, ctr_lon, ctr_alt,
                        outgrids_3d[4], field_names[4], field_descriptions[4], format='f', 
                        grid_units=field_units_3d[4],
                        **output_kwargs)    
        output_writer_3d(outfiles_3d[5], t_ref, np.asarray(t_edges_seconds[:-1]),
                        x_coord*spatial_scale_factor, y_coord*spatial_scale_factor,
                        z_coord*spatial_scale_factor, 
                        lons_3d, lats_3d, alts_3d, ctr_lat, ctr_lon, ctr_alt,
                        outgrids_3d[5], field_names[5], field_descriptions[5], format='f', 
                        grid_units=field_units_3d[5],
                        **output_kwargs)                
        output_writer_3d(outfiles_3d[6], t_ref, np.asarray(t_edges_seconds[:-1]),
                        x_coord*spatial_scale_factor, y_coord*spatial_scale_factor,
                        z_coord*spatial_scale_factor, 
                        lons_3d, lats_3d, alts_3d, ctr_lat, ctr_lon, ctr_alt,
                        outgrids_3d[6], field_names[6], field_descriptions[6], format='f', 
                        grid_units=field_units_3d[6],
                        **output_kwargs)                
                        
        # return x_coord, y_coord, z_coord, lons, lats, alts, extent_density_grid, outfiles, field_names

def grid_h5flashfiles(h5_filenames, start_time, end_time, **kwargs):
    """ Grid LMA data contianed in HDF5-format LMA files. Keyword arguments to this function
        are those to the FlashGridder class and its functions.
    
        This function is provided as a convenience for compatibility with earlier 
        implementations of flash gridding.
    """
    
    process_flash_kwargs = {}
    for prock in ('min_points_per_flash'):
        if prock in kwargs:
            process_flash_kwargs[prock] = kwargs.pop(proc)

    out_kwargs = {}
    for outk in ('outpath', 'output_writer', 'output_writer_3d', 'output_kwargs',
                 'output_filename_prefix', 'spatial_scale_factor'):
        if outk in kwargs:
            out_kwargs[outk] = kwargs.pop(outk)
    
    gridder = FlashGridder(start_time, end_time, **kwargs)
    gridder.process_flashes(h5_filenames, **process_flash_kwargs)
    output = gridder.write_grids(**out_kwargs)
    return output    


if __name__ == '__main__':
    h5_filenames = glob.glob('data/LYL*090610_20*.h5')
    start_time = datetime(2009,6,10, 20,0,0)
    end_time   = datetime(2009,6,10, 21,0,0)
    
    frame_interval=120.0
    dx=4.0e3
    dy=4.0e3
    x_bnd = (-100e3, 100e3)
    y_bnd = (-100e3, 100e3)
    # # KOUN
    ctr_lat = 35.23833
    ctr_lon = -97.46028

    # DC
    # ctr_lat =  38.889444
    # ctr_lon =  -77.035278
    
    grid_h5flashfiles(h5_filenames, start_time, end_time, frame_interval=frame_interval, 
                dx=dx, dy=dy, x_bnd=x_bnd, y_bnd=y_bnd, ctr_lon=ctr_lon, ctr_lat=ctr_lat)

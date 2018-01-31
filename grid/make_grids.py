from __future__ import absolute_import
from __future__ import print_function
import glob 
import os, sys
from datetime import datetime, timedelta

import numpy as np

from . import density_to_files
from lmatools.stream.subset import broadcast
from lmatools.io.LMA_h5_file import read_flashes, to_seconds

from lmatools.coordinateSystems import MapProjection, GeographicSystem
from six.moves import range

from .cf_netcdf import write_cf_netcdf, write_cf_netcdf_3d, write_cf_netcdf_latlon, write_cf_netcdf_3d_latlon, write_cf_netcdf_noproj, write_cf_netcdf_fixedgrid
    

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
    def __init__(self, start_time, end_time, do_3d = True,
                    frame_interval=120.0, dx=4.0e3, dy=4.0e3, dz=1.0e3,
					base_date = None,
                    x_bnd = (-100e3, 100e3),
                    y_bnd = (-100e3, 100e3),
                    z_bnd = (0e3, 20e3),
                    ctr_lat = 35.23833, ctr_lon = -97.46028, ctr_alt=0.0,
                    proj_name = 'aeqd',
                    proj_datum = 'WGS84',
                    proj_ellipse = 'WGS84',
                    pixel_coords = None,
                    flash_count_logfile = None, energy_grids = None,
                    event_grid_area_fraction_key = None):
        """ Class to support gridding of flash and event data. 
                    
            On init, specify the grid 
                    
            If proj_name = 'pixel_grid' then pixel_coords must be an
                instance of lmatools.coordinateSystems.PixelGrid
            If proj_name = 'geos' then pixel_coords must be an instance of 
                lmatools.coordinateSystems.GeostationaryFixedGridSystem
                    
            energy_grids controls which energy grids are saved, default None.
                energy_grids may be True, which will calculate all energy grid 
                types, or it may be one of 'specific_energy', 'total_energy', 
                or a list of one or more of these.
                    
            event_grid_area_fraction_key specifies the name of the variable
                in the events array that gives the fraction of each grid cell
                covered by each event. Used only for pixel-based event 
                detectors.
        """                    
        if energy_grids == True:
            energy_grids = ('specific_energy', 'total_energy')
        self.energy_grids = energy_grids
        
        self.event_grid_area_fraction_key = event_grid_area_fraction_key
        
        # args, kwargs that are saved for the future
        self.do_3d = do_3d
        self.start_time = start_time
        self.dx, self.dy, self.dz = dx, dy, dz
        self.end_time = end_time
        self.frame_interval = frame_interval
        self.base_date = base_date
        self.min_points_per_flash = 1
        self.proj_name = proj_name
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
    
        if self.proj_name == 'latlong':
            dx_units = '{0:6.4f}deg'.format(dx)
            mapProj = GeographicSystem()
        elif self.proj_name == 'pixel_grid':
            dx_units = 'pixel'
            mapProj = pixel_coords
        elif self.proj_name == 'geos':
            dx_units = 'radians'
            mapProj = pixel_coords
        else:
            dx_units = '{0:5.1f}m'.format(dx)
            mapProj = MapProjection(projection=self.proj_name, ctrLat=ctr_lat, ctrLon=ctr_lon, lat_ts=ctr_lat, 
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
        if self.event_grid_area_fraction_key is not None:
            extent_density_grid = np.zeros((xedge.shape[0]-1, yedge.shape[0]-1, n_frames), dtype='float32')
        else:
            extent_density_grid = np.zeros((xedge.shape[0]-1, yedge.shape[0]-1, n_frames), dtype='int32')
        footprint_grid      = np.zeros((xedge.shape[0]-1, yedge.shape[0]-1, n_frames), dtype='float32')

        specific_energy_grid = np.zeros((xedge.shape[0]-1, yedge.shape[0]-1, n_frames), dtype='float32')
        total_energy_grid = np.zeros((xedge.shape[0]-1, yedge.shape[0]-1, n_frames), dtype='float32')
        flashsize_std_grid  = np.zeros((xedge.shape[0]-1, yedge.shape[0]-1, n_frames), dtype='float32')

        if do_3d == True:
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
                
        if do_3d == True:
            self.outgrids_3d = (extent_density_grid_3d,
                    init_density_grid_3d,
                    event_density_grid_3d,
                    footprint_grid_3d,
                    specific_energy_grid_3d,
                    flashsize_std_grid_3d,
                    total_energy_grid_3d
                    )
        else:
            self.outgrids_3d = None
        
        
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
            accum_extent_density = density_to_files.accumulate_points_on_grid(extent_density_grid[:,:,i], xedge, yedge, out=extent_out,label='extent',  grid_frac_weights=True)
            accum_footprint      = density_to_files.accumulate_points_on_grid(footprint_grid[:,:,i], xedge, yedge, label='footprint',  grid_frac_weights=True)

            accum_specific_energy = density_to_files.accumulate_energy_on_grid(specific_energy_grid[:,:,i], xedge, yedge, out=extent_out, label='specific_energy',  grid_frac_weights=True)
            accum_flashstd       = density_to_files.accumulate_points_on_grid_sdev(flashsize_std_grid[:,:,i], footprint_grid[:,:,i], xedge, yedge, out=extent_out, label='flashsize_std',  grid_frac_weights=True)
            accum_total_energy   = density_to_files.accumulate_energy_on_grid(total_energy_grid[:,:,i], xedge, yedge, out=extent_out, label='total_energy',  grid_frac_weights=True)

            if do_3d == True:
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

            if do_3d == True:
                extent_out_3d['func'] = accum_extent_density_3d
                init_out_3d['func'] = accum_init_density_3d
                event_out_3d['func'] = accum_event_density_3d        
        
            event_density_target  = density_to_files.point_density(accum_event_density)
            init_density_target   = density_to_files.point_density(accum_init_density)
            extent_density_target = density_to_files.extent_density(x0, y0, dx, dy, accum_extent_density,
                event_grid_area_fraction_key=event_grid_area_fraction_key)
            mean_footprint_target = density_to_files.extent_density(x0, y0, dx, dy, accum_footprint, weight_key='area',
                event_grid_area_fraction_key=event_grid_area_fraction_key)
            mean_energy_target    = density_to_files.extent_density(x0, y0, dx, dy, accum_specific_energy, weight_key='specific_energy',
                event_grid_area_fraction_key=event_grid_area_fraction_key) #tot_energy
            mean_total_energy_target = density_to_files.extent_density(x0, y0, dx, dy, accum_total_energy, weight_key='total_energy',
                event_grid_area_fraction_key=event_grid_area_fraction_key)  #Energy
            std_flashsize_target  = density_to_files.extent_density(x0, y0, dx, dy, accum_flashstd, weight_key='area',
                event_grid_area_fraction_key=event_grid_area_fraction_key)

            if do_3d == True:
                event_density_target_3d  = density_to_files.point_density_3d(accum_event_density_3d)
                init_density_target_3d   = density_to_files.point_density_3d(accum_init_density_3d)
                extent_density_target_3d = density_to_files.extent_density_3d(x0, y0, z0, dx, dy, dz, accum_extent_density_3d)
                mean_footprint_target_3d = density_to_files.extent_density_3d(x0, y0, z0, dx, dy, dz, accum_footprint_3d, weight_key='area')
            
                mean_energy_target_3d    = density_to_files.extent_density_3d(x0, y0, z0, dx, dy, dz, accum_specific_energy_3d, weight_key='specific_energy') #tot_energy
                mean_total_energy_target_3d = density_to_files.extent_density_3d(x0, y0, z0, dx, dy, dz, accum_total_energy_3d, weight_key='total_energy')     #Energy
                std_flashsize_target_3d  = density_to_files.extent_density_3d(x0, y0, z0, dx, dy, dz, accum_flashstd_3d, weight_key='area')

            broadcast_targets = ( 
                density_to_files.project('lon', 'lat', 'alt', mapProj, geoProj, event_density_target, use_flashes=False),
                density_to_files.project('init_lon', 'init_lat', 'init_alt', mapProj, geoProj, init_density_target, use_flashes=True),
                density_to_files.project('lon', 'lat', 'alt', mapProj, geoProj, extent_density_target, use_flashes=False),
                density_to_files.project('lon', 'lat', 'alt', mapProj, geoProj, mean_footprint_target, use_flashes=False),
                density_to_files.project('lon', 'lat', 'alt', mapProj, geoProj, std_flashsize_target, use_flashes=False),
                )
            if energy_grids is not None:
                if ('specific_energy' == energy_grids) | ('specific_energy' in energy_grids):
                    broadcast_targets += (
                        density_to_files.project('lon', 'lat', 'alt', mapProj, geoProj, mean_total_energy_target, use_flashes=False),
                        )
                if ('total_energy' == energy_grids) | ('total_energy' in energy_grids):
                    broadcast_targets += (
                        density_to_files.project('lon', 'lat', 'alt', mapProj, geoProj, mean_total_energy_target, use_flashes=False),
                        )
            if do_3d == True:
                broadcast_targets += (
                    density_to_files.project('lon', 'lat', 'alt', mapProj, geoProj, event_density_target_3d, use_flashes=False),
                    density_to_files.project('init_lon', 'init_lat', 'init_alt', mapProj, geoProj, init_density_target_3d, use_flashes=True),
                    density_to_files.project('lon', 'lat', 'alt', mapProj, geoProj, extent_density_target_3d, use_flashes=False),
                    density_to_files.project('lon', 'lat', 'alt', mapProj, geoProj, mean_footprint_target_3d, use_flashes=False),
                    density_to_files.project('lon', 'lat', 'alt', mapProj, geoProj, std_flashsize_target_3d, use_flashes=False),
                    )
                if energy_grids is not None:
                    if ('specific_energy' == energy_grids) | ('specific_energy' in energy_grids):
                        broadcast_targets += (
                            density_to_files.project('lon', 'lat', 'alt', mapProj, geoProj, mean_energy_target_3d, use_flashes=False),
                        )
                    if ('total_energy' == energy_grids) | ('total_energy' in energy_grids):
                        broadcast_targets += (
                            density_to_files.project('lon', 'lat', 'alt', mapProj, geoProj, mean_total_energy_target_3d, use_flashes=False),
                        )


            spew_to_density_types = broadcast( broadcast_targets )

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

        energy_grids = self.energy_grids
                    
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
        
        if self.event_grid_area_fraction_key is not None:
            extent_format='f'
        else:
            extent_format='i'
        
        output_writer(outfiles[0], t_ref, np.asarray(t_edges_seconds[:-1]),
                        x_coord*spatial_scale_factor, y_coord*spatial_scale_factor, 
                        lons, lats, ctr_lat, ctr_lon, 
                        outgrids[0], field_names[0], field_descriptions[0],
                        format=extent_format,
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
        if energy_grids is not None:
            if ('specific_energy' == energy_grids) | ('specific_energy' in energy_grids):
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
        if energy_grids is not None:
            if ('total_energy' == energy_grids) | ('total_energy' in energy_grids):
                output_writer(outfiles[6], t_ref, np.asarray(t_edges_seconds[:-1]),
                                x_coord*spatial_scale_factor, y_coord*spatial_scale_factor, 
                                lons, lats, ctr_lat, ctr_lon, 
                                outgrids[6], field_names[6], field_descriptions[6], format='f', 
                                grid_units=field_units[6],
                                **output_kwargs)
    
        ########3D:
        if self.outgrids_3d is not None:
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
            if energy_grids is not None:
                if ('specific_energy' == energy_grids) | ('specific_energy' in energy_grids):
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
            if energy_grids is not None:
                if ('total_energy' == energy_grids) | ('total_energy' in energy_grids):
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
    for prock in ('min_points_per_flash',):
        if prock in kwargs:
            process_flash_kwargs[prock] = kwargs.pop(prock)

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

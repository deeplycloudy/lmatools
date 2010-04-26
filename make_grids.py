import glob 
import os
from datetime import datetime, timedelta

import numpy as np
import tables

import density_to_files

from coordinateSystems import MapProjection, GeographicSystem


def to_seconds(dt):
    'convert datetime.timedelta object to float seconds'
    return dt.days * 86400 + dt.seconds + dt.microseconds/1.0e6


def write_cf_netcdf(outfile, t_start, t, xloc, yloc, lon_for_x, lat_for_y, ctr_lat, ctr_lon, grid, grid_var_name, grid_description, format='i'):
    """ Write a Climate and Forecast Metadata-compliant NetCDF file. 
    
        Should display natively in conformant packages like McIDAS-V.
        
    """

    import pupynere as nc

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
    
    lons = nc_out.createVariable('lons', 'd', ('nx','ny') )#, filters=no_compress)
    lons.long_name="longitude"
    lons.standard_name="longitude"
    lons.units = "degrees_east"
    
    lats = nc_out.createVariable('lats', 'd', ('nx','ny') )#, filters=no_compress)
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



def read_flashes(h5LMAfiles, target, base_date=None, min_points=10):
    """ This routine is the head of the data pipeline, responsible for pushing out 
        events and flashes. Using a different flash data source format is a matter of
        replacing this routine to read in the necessary event and flash data."""

    for filename in h5LMAfiles:
        h5 = tables.openFile(filename)
        table_name = h5.root.events._v_children.keys()[0]

        print filename

        # event_dtype = getattr(h5.root.events, table_name)[0].dtype
        # events_all = getattr(h5.root.events, table_name)[:]
        flashes = getattr(h5.root.flashes, table_name)
        events = getattr(h5.root.events, table_name)[:]

        fl_dtype = flashes[0].dtype
        flashes = np.fromiter((fl[:] for fl in flashes if fl['n_points'] >= min_points), dtype=fl_dtype)

        # get the start date of the file, which is the seconds-of-day reference
        h5_start =  datetime(*getattr(h5.root.events, table_name).attrs['start_time'][0:3])
        if base_date==None:
            base_date = h5_start
        extra_dt = to_seconds(h5_start - base_date)
        # adjust time in seconds of day to be in correct reference
        events['time'] += extra_dt
        flashes['start'] += extra_dt

        n_flashes = len(flashes)
        print '    ', n_flashes, ' total flashes, with extra dt of', extra_dt

        target.send((events,flashes))

        h5.close()
        

def grid_h5flashfiles(h5_filenames, start_time, end_time, 
                        frame_interval=120.0, dx=4.0e3, dy=4.0e3,
                        x_bnd = (-100e3, 100e3),
                        y_bnd = (-100e3, 100e3),
                        z_bnd = (-20e3, 20e3),
                        ctr_lat = 35.23833, ctr_lon = -97.46028,
                        min_points_per_flash=10,
                        outpath = ''
                        ):
    from math import ceil
    """
    
    Create 2D plan-view density grids for events, flash origins, flash extents, and mean flash footprint
    
    frame_interval: Frame time-step in seconds
    dx, dy: horizontal grid size in m
    {x,y,z}_bnd: horizontal grid edges in m
    ctr_lat, ctr_lon: coordinate center
    
    Uses an azimuthal equidistant map projection on the WGS84 ellipsoid.
    
    
    read_flashes
    
    filter_flash
    extract_events
    flash_to_frame
    
    frame0_broadcast, frame1_broadcast, ...
    
    each broadcaster above sends events and flashes to:
    projection( event_location), projection(flash_init_location), projection(event_location)
    which map respectively to:
    point_density->accum_on_grid(event density), point_density->accum_on_grid(flash init density), extent_density->accum_on_grid(flash_extent_density)

    grids are in an HDF5 file. how to handle flushing?
    """
    
    # reference time is the date part of the start_time
    ref_date = start_time.date()
    t_ref = datetime(ref_date.year, ref_date.month, ref_date.day)
    
    frame_dt = timedelta(0, frame_interval, 0)
    
    duration = end_time - start_time
    n_frames = int(ceil(to_seconds(duration) / to_seconds(frame_dt)))

    t_edges = [start_time + i*frame_dt for i in range(n_frames+1)]
    
    t_edges_seconds = [to_seconds(edge - t_ref) for edge in t_edges]
    
    xedge=np.arange(x_bnd[0], x_bnd[1]+dx, dx)
    yedge=np.arange(y_bnd[0], y_bnd[1]+dy, dy)
    
    x0 = xedge[0]
    y0 = yedge[0]
    
    mapProj = MapProjection(projection='aeqd', ctrLat=ctr_lat, ctrLon=ctr_lon, lat_ts=ctr_lat, lon_0=ctr_lon, lat_0=ctr_lat)
    geoProj = GeographicSystem()
    
    
    event_density_grid  = np.zeros((xedge.shape[0]-1, yedge.shape[0]-1, n_frames), dtype='int32')
    init_density_grid   = np.zeros((xedge.shape[0]-1, yedge.shape[0]-1, n_frames), dtype='int32')
    extent_density_grid = np.zeros((xedge.shape[0]-1, yedge.shape[0]-1, n_frames), dtype='int32')
    footprint_grid      = np.zeros((xedge.shape[0]-1, yedge.shape[0]-1, n_frames), dtype='float32')
        
    all_frames = []
    extent_frames = []
    init_frames = []
    event_frames = []
    for i in range(n_frames):
        extent_out = {'name':'extent'}
        init_out   = {'name':'init'}
        event_out  = {'name':'event'}
        accum_event_density  = density_to_files.accumulate_points_on_grid(event_density_grid[:,:,i], xedge, yedge,  out=event_out, label='event')
        accum_init_density   = density_to_files.accumulate_points_on_grid(init_density_grid[:,:,i], xedge, yedge,   out=init_out,  label='init')
        accum_extent_density = density_to_files.accumulate_points_on_grid(extent_density_grid[:,:,i], xedge, yedge, out=extent_out,label='extent')
        accum_footprint      = density_to_files.accumulate_points_on_grid(footprint_grid[:,:,i], xedge, yedge, label='footprint')
        extent_out['func'] = accum_extent_density
        init_out['func'] = accum_init_density
        event_out['func'] = accum_event_density
        extent_frames.append(extent_out)
        init_frames.append(init_out)
        event_frames.append(event_out)
        
        event_density_target  = density_to_files.point_density(accum_event_density)
        init_density_target   = density_to_files.point_density(accum_init_density)
        extent_density_target = density_to_files.extent_density(x0, y0, dx, dy, accum_extent_density)
        mean_footprint_target = density_to_files.extent_density(x0, y0, dx, dy, accum_footprint, weight_key='area')

        spew_to_density_types = density_to_files.broadcast( ( 
                    density_to_files.project('lon', 'lat', 'alt', mapProj, geoProj, event_density_target, use_flashes=False),
                    density_to_files.project('init_lon', 'init_lat', 'init_alt', mapProj, geoProj, init_density_target, use_flashes=True),
                    density_to_files.project('lon', 'lat', 'alt', mapProj, geoProj, extent_density_target, use_flashes=False),
                    density_to_files.project('lon', 'lat', 'alt', mapProj, geoProj, mean_footprint_target, use_flashes=False),
                    ) )

        all_frames.append( density_to_files.extract_events_for_flashes( spew_to_density_types ) )
        
    framer = density_to_files.flashes_to_frames(t_edges_seconds, all_frames, time_key='start')
    
    read_flashes( h5_filenames, framer, base_date=t_ref, min_points=min_points_per_flash)
    
    # print 'event_density_grid ', id(event_density_grid[:,:,-1])
    # print 'extent_density_grid', id(extent_density_grid[:,:,-1])
    # print 'init_density_grid  ', id(init_density_grid[:,:,-1])
    
    
    x_coord = (xedge[:-1] + xedge[1:])/2.0
    y_coord = (yedge[:-1] + yedge[1:])/2.0
    nx = x_coord.shape[0]
    ny = y_coord.shape[0]
    
    x_all, y_all = (a.T for a in np.meshgrid(x_coord, y_coord))
    assert x_all.shape == y_all.shape
    z_all = np.zeros_like(x_all)
            
    lons, lats, alts = x,y,z = geoProj.fromECEF( *mapProj.toECEF(x_all, y_all, z_all) )
    lons.shape=x_all.shape
    lats.shape=y_all.shape
    
    outflile_basename = os.path.join(outpath,'LMA_%s_%d_' % (start_time.strftime('%Y%m%d_%H%M%S'), to_seconds(duration)))
    
    outfiles = (outflile_basename+'flash_extent.nc',
                outflile_basename+'flash_init.nc',
                outflile_basename+'source.nc',
                outflile_basename+'footprint.nc',
                )
                
    outgrids = (extent_density_grid,
                init_density_grid,
                event_density_grid,
                footprint_grid
                )
                
    field_names = ('flash_extent', 'flash_initiation', 'lma_source', 'flash_footprint')
    
    field_descriptions = ('LMA flash extent density',
                        'LMA flash initiation density',
                        'LMA source density',
                        'LMA per-flash footprint')
                        
    write_cf_netcdf(outfiles[0], t_ref, np.asarray(t_edges_seconds[:-1]),
                    x_coord/1.0e3, y_coord/1.0e3, lons, lats, ctr_lat, ctr_lon, 
                    outgrids[0], fieldnames[0], field_descriptions[0])
    write_cf_netcdf(outfiles[1], t_ref, np.asarray(t_edges_seconds[:-1]),
                    x_coord/1.0e3, y_coord/1.0e3, lons, lats, ctr_lat, ctr_lon, 
                    outgrids[1], fieldnames[1], field_descriptions[1])
    write_cf_netcdf(outfiles[2], t_ref, np.asarray(t_edges_seconds[:-1]),
                    x_coord/1.0e3, y_coord/1.0e3, lons, lats, ctr_lat, ctr_lon, 
                    outgrids[2], fieldnames[2], field_descriptions[2])
    write_cf_netcdf(outfiles[3], t_ref, np.asarray(t_edges_seconds[:-1]),
                    x_coord/1.0e3, y_coord/1.0e3, lons, lats, ctr_lat, ctr_lon, 
                    outgrids[3], fieldnames[3], field_descriptions[3], format='f')
                    
    print 'max extent is', extent_density_grid.max()

    return x_coord, y_coord, lons, lats, extent_density_grid, outfiles, field_names


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
from __future__ import absolute_import
from __future__ import print_function

from functools import partial

try:
    from netCDF4 import Dataset as NetCDFFile
except ImportError:
    from scipy.io.netcdf import NetCDFFile

def write_cf_netcdf_latlon(outfile, t_start, t, xloc, yloc, 
    lon_for_x=None, lat_for_y=None, ctr_lat=None, ctr_lon=None, 
    grid=None, grid_var_name='lightning',
    grid_description='some kind of lightning data', 
    format='i', missing_value=-9999, is_latlon=True, is_fixedgrid=False,
    **kwargs):
    """ 
    Write a Climate and Forecast Metadata-compliant NetCDF file. 
        
    Grid is assumed to be regular in lon, lat and so no map projection
    information is necessary. If is_latlon=False, lat_for_y, and 
    lon_for_x, ctr_lat, and ctr_lon will be used as coordinates for a 
    Lambert Azimuthal Equal Area projection with the specified center.
    
    If is_projected=False, and is_latlon=False, then x and y are treated
    as unitless pixel coordinates with no physical meaning. This is useful
    for arbitrarily stretched lat/lon grids, or explicitly mappings between 
    a sensor pixel array and geolocated coordinates.
    
    If is_projected=True and is_fixedgrid=True, the optional keyword argument
    nadir_lon should be passed with the longitude of the satellite sub-point.
    The fixed grid is assumed to be the GOES-R ABI fixed grid, so is not
    readily adapted to other satellite platforms.
    
    Should display natively in conformant packages like McIDAS-V.
    """
    if grid is None:
        raise ValueError("An array must be given for grid")

    grid_units = kwargs.pop('grid_units', 'dimensionless')
    is_projected = kwargs.pop('is_projected', True)
        
    east_dim = 'lon'
    north_dim = 'lat'
    grid_mapping = "crs"
    
    if not is_latlon:
        east_dim = 'nx'
        north_dim = 'ny'
        if is_projected:
            if is_fixedgrid:
                grid_mapping  = "goes_imager_projection"
            else:
                grid_mapping = "Lambert_Azimuthal_Equal_Area"

    nc_out = NetCDFFile(outfile, 'w')
    nc_out.createDimension(east_dim, xloc.shape[0])
    nc_out.createDimension(north_dim, yloc.shape[0])
    nc_out.createDimension('ntimes', t.shape[0])  #unlimited==None

    # declare the coordinate reference system, WGS84 values
    if is_projected:
        proj = nc_out.createVariable(grid_mapping, 'i', ())
        if is_latlon:
            proj.grid_mapping_name = 'latitude_longitude'
            proj.longitude_of_prime_meridian = 0.0 
            proj.semi_major_axis = 6378137.0 
            proj.inverse_flattening = 298.257223563 
        elif is_fixedgrid:
            proj.long_name = 'GOES-R ABI fixed grid projection'
            proj.grid_mapping_name = 'geostationary'
            proj.perspective_point_height = 35786023.0
            proj.semi_major_axis = 6378137.0
            proj.semi_minor_axis = 6356752.31414
            proj.inverse_flattening = 298.2572221
            proj.latitude_of_projection_origin = 0.0
            proj.longitude_of_projection_origin = kwargs['nadir_lon']
            proj.sweep_angle_axis = 'x'
        else:
            proj.grid_mapping_name = 'lambert_azimuthal_equal_area'
            proj.longitude_of_projection_origin = ctr_lon
            proj.latitude_of_projection_origin = ctr_lat
            proj.false_easting = 0.0
            proj.false_northing = 0.0
    
    if is_latlon:
        y_coord = nc_out.createVariable('latitude', 'f', (north_dim,))
        y_coord.units = "degrees_north"
        y_coord.long_name = "latitude"
        y_coord.standard_name = 'latitude'
        
        x_coord = nc_out.createVariable('longitude', 'f', (east_dim,))
        x_coord.units = "degrees_east"
        x_coord.long_name = "longitude"
        x_coord.standard_name = 'longitude'
    else:
        if is_projected:
            var_type = 'f'
        else:
            var_type = 'i' # integer pixel coordinates
        x_coord = nc_out.createVariable('x', var_type, (east_dim,))
        y_coord = nc_out.createVariable('y', var_type, (north_dim,))
        if is_projected:
            if is_fixedgrid:
                x_coord.units = "rad"
                x_coord.long_name = "GOES fixed grid projection x-coordinate"
                x_coord.standard_name = 'projection_x_coordinate'
                y_coord.units = "rad"
                y_coord.long_name = "GOES fixed grid projection y-coordinate"
                y_coord.standard_name = 'projection_y_coordinate'
            else:
                x_coord.units = "km"
                x_coord.long_name = "x coordinate of projection"
                x_coord.standard_name = 'projection_x_coordinate'        
                y_coord.units = "km"
                y_coord.long_name = "y coordinate of projection"
                y_coord.standard_name = 'projection_y_coordinate'

    times = nc_out.createVariable('time', 'f', ('ntimes',) )#, filters=no_compress)
    times.long_name="time"
    times.units = "seconds since %s" % t_start.strftime('%Y-%m-%d %H:%M:%S')
    
    #Dtype change from 'd' to 'f':
    if (not is_latlon) & (not is_fixedgrid):
        lons = nc_out.createVariable('lons', 'f', (east_dim, north_dim))
        lons.long_name="longitude"
        lons.standard_name="longitude"
        lons.units = "degrees_east"
        
        lats = nc_out.createVariable('lats', 'f', (east_dim, north_dim))
        lats.long_name="latitude"
        lats.standard_name="latitude"
        lats.units = "degrees_north"

        lons[:] = lon_for_x[:]
        lats[:] = lat_for_y[:]

    lightning2d = nc_out.createVariable(grid_var_name, format,
        ('ntimes',east_dim,north_dim), zlib=True)
    lightning2d.long_name=grid_description
    lightning2d.units=grid_units
    if (not is_latlon) & (not is_fixedgrid):
        lightning2d.coordinates='time lons lats'
    lightning2d.grid_mapping = grid_mapping
    lightning2d.missing_value = missing_value

    x_coord[:] = xloc[:]
    y_coord[:] = yloc[:]
    times[:] = t[:]

    for i in range(grid.shape[2]):
        lightning2d[i,:,:] = grid[:,:,i]
    nc_out.sync()
    nc_out.close()

write_cf_netcdf = partial(write_cf_netcdf_latlon, is_latlon=False)
write_cf_netcdf_noproj = partial(write_cf_netcdf, is_projected=False)
write_cf_netcdf_fixedgrid = partial(write_cf_netcdf, is_fixedgrid=True)

def write_cf_netcdf_3d(outfile, t_start, t, xloc, yloc, zloc, lon_for_x, lat_for_y, alt_for_z, ctr_lat, ctr_lon, ctr_alt, grid, grid_var_name, grid_description, format='i', **kwargs):
    """ Write a Climate and Forecast Metadata-compliant NetCDF file. 
    
        Should display natively in conformant packages like McIDAS-V.
        
    """


    missing_value = -9999
    
    nc_out = NetCDFFile(outfile, 'w')
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

    lightning3d = nc_out.createVariable(grid_var_name, format, ('ntimes','nx','ny', 'nz'), zlib=True)#, filters=no_compress)
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


    missing_value = -9999
    
    nc_out = NetCDFFile(outfile, 'w')
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
    
    lightning3d = nc_out.createVariable(grid_var_name, format, ('ntimes','lon','lat','alt'), zlib=True )#, filters=no_compress)
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

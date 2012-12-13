
def write_cf_netcdf_latlon(outfile, t_start, t, xloc, yloc, lon_for_x, lat_for_y, ctr_lat, ctr_lon, grid, grid_var_name, grid_description, format='i', **kwargs):
    """ Write a Climate and Forecast Metadata-compliant NetCDF file. 
        
        Grid is regular in lon, lat and so no map projection information is necessary.
    
        Should display natively in conformant packages like McIDAS-V.
        
    """

    import pupynere as nc

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
from __future__ import absolute_import
from __future__ import print_function
import datetime
import scipy.io.netcdf as nc
import glob
from .make_grids import grid_h5flashfiles

def write_AWIPS_netcdf_grid(outfile, t_start, t, xloc, yloc, lon_for_x, lat_for_y, 
                ctr_lat, ctr_lon, grid, grid_var_name, grid_description, format='i',
                refresh_minutes = 1, grid_units = 'dimensionless'):
    """ Write an AWIPS-formatted NetCDF grid for lightning data. 
       
       TYPEMAP = { NC_BYTE:   ('b', 1),
                   NC_CHAR:   ('c', 1),
                   NC_SHORT:  ('h', 2),
                   NC_INT:    ('i', 4), # this is also long
                   NC_FLOAT:  ('f', 4),
                   NC_DOUBLE: ('d', 8) } 
    """

    

    missing_value = -9999
    # all_levels = (1,)
    level_name = "SFC"
    
    # ----------    
    # Dimensions
    # ----------
    nc_out = nc.NetCDFFile(outfile, 'w')
    nc_out.createDimension('record', None)  #unlimited==None
    nc_out.createDimension('n_valtimes', t.shape[0])
    nc_out.createDimension('data_variables', 1)
    nc_out.createDimension('namelen', 132)
    nc_out.createDimension('charsPerLevel', 20)
    nc_out.createDimension('levels', 1)
    # for level in all_levels:
        # nc_out.createDimension('levels_{0}'.format(level), level)
    nc_out.createDimension('x', xloc.shape[0])
    nc_out.createDimension('y', yloc.shape[0])
    
    # -----------------
    # Global attributes 
    # -----------------
    nc_out.Data_record_time = t_start.strftime('%Y%m%d.%H%M')
    nc_out.depictorName = "WTLMA_ltng_grid"
    nc_out.projIndex = 8
    nc_out.projName  = "CYLINDRICAL_EQUIDISTANT" 
    nc_out.centralLon = ctr_lon
    nc_out.centralLat = ctr_lat
    nc_out.xMin    = lon_for_x.min()
    nc_out.xMax    = lon_for_x.max()
    nc_out.yMin    = lat_for_y.min()
    nc_out.yMax    = lat_for_y.max()
    nc_out.lat00   = nc_out.yMin
    nc_out.lon00   = nc_out.xMin
    nc_out.latNxNy = nc_out.yMax
    nc_out.lonNxNy = nc_out.xMax
    nc_out.dxKm    = nc_out.xMin
    nc_out.dyKm    = nc_out.yMax
    nc_out.latDxDy = (nc_out.yMin + nc_out.yMax)/2.0
    nc_out.lonDxDy = (nc_out.xMin + nc_out.xMax)/2.0
    nc_out.UnitSquareSideSize = xloc[1]-xloc[0]
    nc_out.RefreshPeriod = 60*int(refresh_minutes)
    # nc_out.NbRefreshPeriods = int(10*60/nc_out.RefreshPeriod)
    
    
    charsPerLevel = nc_out.createVariable(grid_var_name+'Levels','c',('levels','charsPerLevel'))
    inventory = nc_out.createVariable(grid_var_name+'Inventory','c',('n_valtimes','charsPerLevel'))
    charsPerLevel[0, 0:len(level_name)] = level_name
    for i in range(nc_out.dimensions['n_valtimes']):
        inventory[i,:] = '1'*nc_out.dimensions['charsPerLevel'] # What is this? Why AWIPS, Why?

    the_epoch = datetime.datetime(1970,1,1,0,0,0)
    reftime_since_the_epoch = int((t_start-the_epoch).total_seconds())
    valtimeref_delta = nc_out.createVariable('valtimeMINUSreftime','i',('n_valtimes',))
    valtime = nc_out.createVariable('valtimeMINUSreftime','i',('n_valtimes',))
    reftime = nc_out.createVariable('valtimeMINUSreftime','i',())
    reftime.long_name="reference time"
    reftime.units="seconds since (1970-1-1 00:00:00.0)"
    reftime.assignValue(reftime_since_the_epoch)
    valtime.long_name="valid time"
    valtime.units="seconds since (1970-1-1 00:00:00.0)"
    valtime[:] = t + reftime_since_the_epoch
    valtimeref_delta[:] = t
    
    lightning2d = nc_out.createVariable(grid_var_name, format, ('record','levels','y','x') )#, filters=no_compress)
    lightning2d.long_name=grid_description #'LMA VHF event counts (vertically integrated)'
    lightning2d.units=grid_units
    lightning2d.uiname=grid_var_name
    lightning2d._FillValue = missing_value
    lightning2d._n3D = 1
    lightning2d.levels = level_name
    
    # Loop over times
    for i in range(grid.shape[2]):
        lightning2d[i,0,:,:] = grid[:,:,i].T
    nc_out.close()



def make_AWIPS_grid(h5_filenames=None):
    """ Write LMA data to a Lambert Conformal grid in AWIPS NetCDF Grid format"""
    
    if h5_filenames is None:
        h5_filenames = glob.glob('/data/WTLMA/FirstLightning/processed7stn/thresh-0.15_dist-3000.0/LYL*.flash.h5')
    
    print(h5_filenames)

    start_time = datetime.datetime(2011,11,21, 14,50,00)
    end_time   = datetime.datetime(2011,11,21, 14,51,00)
    frame_interval=60.0
        
    # ctr_lat = 25.0
    # ctr_lon = -95.0
    ctr_lat = 33.7
    ctr_lon = -101.8
    target_datum = 'NAD83'
    target_ellipse = 'GRS80'
    
    
    # gs = GeographicSystem()
    # mp = MapProjection('lcc', 25, -95, 'GRS80', 'NAD83', lat_1=25, lon_0=-95)
    # X,Y,Z = gs.toECEF(-101.8, 33.7, 1000)
    # x,y,z = mp.fromECEF(X,Y,Z)
    # #      (-620594.38006647746, 976685.16351688653, 1000.0000000018626)
    # x0, y0 = -620.0e3, 975.0e3 # Lubbock in system centered on 25, -95
    x0, y0 = 0.0, 0.0
    
    dx=2.0e3
    dy=2.0e3
    x_bnd = (x0-200e3, x0+200e3)
    y_bnd = (y0-200e3, y0+200e3)
    
    # Lambert Conformal Conic, one standard parallel
    # See http://remotesensing.org/geotiff/proj_list/lambert_conic_conformal_1sp.html
    # target_proj = 'lcc'
    target_proj = 'eqc'
        
    x_coord, y_coord, lons, lats, extent_density_grid, outfiles, field_names = grid_h5flashfiles(
                h5_filenames, start_time, end_time, frame_interval=frame_interval, 
                dx=dx, dy=dy, x_bnd=x_bnd, y_bnd=y_bnd, ctr_lon=ctr_lon, ctr_lat=ctr_lat,
                proj_name = target_proj, proj_datum = target_datum, proj_ellipse = target_ellipse,
                output_writer = write_AWIPS_netcdf_grid, output_kwargs = {'refresh_minutes':int(frame_interval/60.0)}
                )
if __name__ == '__main__':
    make_AWIPS_grid()
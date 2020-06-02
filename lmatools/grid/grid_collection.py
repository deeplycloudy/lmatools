from __future__ import absolute_import
from __future__ import print_function
import glob
from datetime import datetime, timedelta
import itertools

import numpy as np

try:
    from netCDF4 import Dataset as NetCDFFile
except ImportError:
    from scipy.io.netcdf import NetCDFFile

from lmatools.vis.multiples_nc import centers_to_edges

class LMAgridFileCollection(object):
    def __init__(self, filenames, grid_name, 
                 x_name='x', y_name='y', t_name='time'):
        """ Load NetCDF files given by filenames and a known grid type.
            The NetCDF files must have been created by lmatools,
            though you might have some luck if the grids are cf-compliant
            more generally.
        
            NCs = LMAgridFileCollection(filenames, 'lma_source')
            
            The data are retrievable by iteration over the 
            collection of files:
            
            >>> for t, xedge, yedge, data in NCs:
            >>>     print(t)

            Or, if you know a time accurately, you can do:
            >>> from datetime import datetime
            >>> t = datetime(2013,6,6,3,0,0)
            >>> xedge, yedge, data = NCs.data_for_time(t)

            The data have been oriented such that a call to matplotlib's
            pcolormesh(xedge,yedge,data)
            will do the expected thing.
        
        """

        self.x_name = x_name
        self.y_name = y_name
        self.t_name = t_name
        self.grid_name = grid_name
        self._filenames = filenames
#         self._ncfiles = [NetCDFFile(f) for f in filenames]
        self._time_lookup = {} # populated by self._frame_times_for_file
        self.times = [t for t in self._all_times()]
        self.times.sort()
        
    def data_for_time(self, t0, return_nc=False):
        """ Read data from the file corresponding to datetime t.
            Returns xedge, yedge, and density for the file 
        
            If return_nc is True, also return the NetCDFFile corresponding
            to the data as a fourth return value
        """
        fname, i = self._time_lookup[t0]  #i is the frame id for this time in NetCDFFile f
        
        f = NetCDFFile(fname)
        
        data = f.variables  # dictionary of variable names to nc_var objects
        dims = f.dimensions # dictionary of dimension names to sizes
        x = data[self.x_name]
        y = data[self.y_name]
        t = data[self.t_name]
        grid = data[self.grid_name]
        indexer = [slice(None),]*len(grid.shape)
        
        grid_dims = grid.dimensions # tuple of dimension names
        name_to_idx = dict((k, i) for i, k in enumerate(grid_dims))
        
        grid_t_idx = name_to_idx[t.dimensions[0]]
        grid_x_idx = name_to_idx[x.dimensions[0]]
        grid_y_idx = name_to_idx[y.dimensions[0]]
                
        xedge = centers_to_edges(x)
        yedge = centers_to_edges(y)
        indexer[grid_t_idx] = i
        density = grid[indexer].transpose()
        out = xedge, yedge, density
        if return_nc:
            out += (f,)
        else:
            f.close()
        return out
        
    def _all_times(self):
        for f in self._filenames:
            for t in self._frame_times_for_file(f):
                yield t    
        
    def _frame_times_for_file(self, fname):
        """ Called once by init to set up frame lookup tables and yield 
            the frame start times. _frame_lookup goes from 
            datetime->(nc file, frame index)"""
        f = NetCDFFile(fname)

        data = f.variables  # dictionary of variable names to nc_var objects
        dims = f.dimensions # dictionary of dimension names to sizes
        t = data[self.t_name]
        
        try:
            base_date = datetime.strptime(t.units, "seconds since %Y-%m-%d %H:%M:%S")
        except ValueError:
            base_date = datetime.strptime(t.units, "seconds since %Y-%m-%d")
        for i in range(np.atleast_1d(t).shape[0]):
            frame_start = base_date + timedelta(0,float(t[i]),0)
            self._time_lookup[frame_start] = (fname, i)
            yield frame_start
        f.close()
        
    def get_projection(self):
        """ Returns GeographicSystem and MapProjection instances from
            lmatools.coordinateSystems corresponding
            to the coordinate system specified by the metadata of the
            first NetCDF file in self._filenames.
        """
        from lmatools.coordinateSystems import GeographicSystem, MapProjection
        geosys = GeographicSystem()
        
        f = NetCDFFile(self._filenames[0])
        
        # Surely someone has written an automated library to parse coordinate 
        # reference data from CF-compliant files.
        if 'Lambert_Azimuthal_Equal_Area' in list(f.variables.keys()):
            nc_proj = f.variables['Lambert_Azimuthal_Equal_Area']
            proj_name = 'laea'
            ctrlon, ctrlat  = (nc_proj.longitude_of_projection_origin,
                                      nc_proj.latitude_of_projection_origin,)
            try:
                ctralt = nc_proj.altitude_of_projection_origin
            except AttributeError:
                print("No altitude attribute in NetCDF projection data, setting to 0.0")
                ctralt = 0.0
            mapproj = MapProjection(proj_name, ctrLat = ctrlat, ctrLon=ctrlon, 
                                     lat_0=ctrlat, lon_0=ctrlon)
            # print geosys.fromECEF(*mapproj.toECEF((0,0), (0,0), (0,0)))
            return geosys, mapproj
        else:
            print("projection not found, assuming lat, lon grid")    
            return geosys, geosys
        f.close()

        
                    
    def __iter__(self):
        for t in self.times:
            xedge, yedge, density = self.data_for_time(t)
            yield t, xedge, yedge, density

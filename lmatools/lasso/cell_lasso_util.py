from __future__ import absolute_import
import os, glob, itertools
from datetime import datetime, timedelta

import logging, json, operator

import numpy as np

from lmatools.io.LMA_h5_file import parse_lma_h5_filename

def gen_polys(filename, time_keys=None):
    """ time_keys is a dictionary mapping from variable names 
        that store time to datetime format strings used to 
        parse those times.

        There are two JSON file formats used. The polygon lasso
        log GUI tool produces a file with one compact-encoded JOSN lasso
        object per line.
        { "created": "2017-08-02T17:19:45.0", ...}
        { "created": "2017-08-02T17:19:50.0", ...}

        An additional format, more suitable for manual editing, has a single
        JSON object "lassos", which has a list of the individual lasso objects.
        {
            "lassos": [
                { "created": "2017-08-02T17:19:45.0", ...}
                { "created": "2017-08-02T17:19:50.0", ...}
            ]
        }
        """

    def parse_one_poly(v):
        # We'll return just the polygon data
        p = v['poly'] 
        
        # Pull the metadata about time of creation into the polygon dict itself
        created = datetime.strptime(v['created'], "%Y-%m-%dT%H:%M:%S.%f")
        p['created']= created

        if time_keys is not None:
            for k in time_keys:
                # Parse into a datetime object
                if k in p:
                    k_date = datetime.strptime(p[k], time_keys[k])
                    p[k] = k_date
        return p
    
    f = open(filename)
    
    s = f.read(100)
    f.seek(0)
    first_json = "".join(s.split()) # strip all whitespace
    if first_json.startswith('{"lassos":[{'):
        polys = json.load(f)
        for v in polys['lassos']:
            yield parse_one_poly(v)
    else:
        for line in f:
            v = json.loads(line)
            yield parse_one_poly(v)
    f.close()
    
def read_polys(filename, sort_key=None, time_keys=None):
    """ Read in a list of polygons, optionally sorting by sort_key.
        time_keys is a dictionary giving datetime.strptime strings to be used
        on the polygon's time data columns named by the keys of the dictionary,
        for example time_keys={'frame_time':'%Y-%m-%dT%H:%M:%S'}
    """
    polys = [p for p in gen_polys(filename, time_keys=time_keys)]
    if sort_key is not None:
        polys.sort(key=operator.itemgetter(sort_key))
    return polys


def apply_filter_poly(filter_poly, polys_log, flash_stat_polys):

    from shapely.geometry import Polygon

    filter_geom = Polygon(filter_poly)
    trimmed_flash_stat_polys = [Polygon(poly).intersection(filter_geom) 
                                for poly in flash_stat_polys]

    # sample_trimmed = trimmed_flash_stat_polys[0]
    updated_coords = [{'x_verts':list(p.exterior.xy[0]), 
                      'y_verts':list(p.exterior.xy[1])} 
                     for p in trimmed_flash_stat_polys ]
    assert len(polys_log) == len(updated_coords)

    # The line below will update the original list of poly dictionaries with new coordinates
    # old_polys = copy.deepcopy(polys)
    for poly, new_coords in zip(polys_log, updated_coords):
        poly.update(new_coords)
        
    # grab the revised coordinates
    flash_stat_polys = tuple(list(zip(p['x_verts'], p['y_verts'])) for p in polys_log)
    return flash_stat_polys


def read_poly_log_file(filename, lat_lon_filter=None):
    """ Import polygons from logger tool.
        
        Returns (polys, t_edges)
        polys: length N sequence of (lon,lat) tuples giving the polygon vertices 
        t_edges: length N+1 sequence of datetime instances giving the 
            start and end times to which each of the N polygons applies. Each
            polygon is valid beginning at the frame_time given in the log file
            for each polygon, and ends when next polygon begins. If the last 
            polygon has frame_end, then that time is used as the end time of the
            last polygon. Otherwise, the last time in t_edges is found by taking
            the frame_time of the last lasso and extending it by the frame_time 
            difference between the second to last and last lassos.
    """

    time_keys = {'frame_time':'%Y-%m-%dT%H:%M:%S',
                 'frame_end':'%Y-%m-%dT%H:%M:%S'}
    polys_log = read_polys(filename, 
                       sort_key = 'frame_time', time_keys=time_keys)

    if 'lon_verts' in polys_log[0]:
        lon_name, lat_name = 'lon_verts', 'lat_verts'
    else:
        "Assuming lon, lat data are in the x_verts, y_verts polygon log entries"
        lon_name, lat_name = 'x_verts', 'y_verts'
    
    flash_stat_polys = tuple(list(zip(p[lon_name], p[lat_name])) for p in polys_log)
    flash_stat_left_time_edges = tuple(p['frame_time'] for p in polys_log)
    if 'frame_end' not in polys_log[-1]:
        dt = flash_stat_left_time_edges[-1] - flash_stat_left_time_edges[-2]
        t_edges = flash_stat_left_time_edges + (flash_stat_left_time_edges[-1] + dt, )
    else:
        flash_stat_right_time_edges = tuple(p['frame_end'] for p in polys_log)
        t_edges = flash_stat_left_time_edges + (flash_stat_right_time_edges[-1],)
        
    if lat_lon_filter is not None:
        flash_stat_polys = apply_filter_poly(lat_lon_filter, polys_log, flash_stat_polys)
    
    return (flash_stat_polys, t_edges)
    
def polys_to_bounding_box(polys):
    """ Given a sequence of (x,y) tuples in polys, return xspan, yspan
        where xspan is a tuple of left, right bounds
        and   yspan is a tuple of bottom, top bounds 
        that encompass all polygons
    """
    all_verts = np.asarray([p for p in itertools.chain(*polys)])
    sw_corner = all_verts.min(axis=0)
    ne_corner = all_verts.max(axis=0)
    lons = (sw_corner[0], ne_corner[0])
    lats = (sw_corner[1], ne_corner[1])
    return lons, lats
    
    
def h5_files_from_standard_path(path_to_sort_results, date_min, date_max):
    """ Path contains a set of directories with h5 files as follows
        path/h5_files/yyyy/Mon/dd/LYLOUT*.h5
        where Mon is the three-letter month name (e.g., Dec).
        This is the directory structure produced by many of the lmatools examples.
        date_min and date_max are datetime objects given the range of dates to examine.
        Return a list of unique HDF5 filenames (sorted by filename) on that path.
    
        Right now, date_min and date_max can only be used to span a single day 
        boundary.
    """
    globber_start = path_to_sort_results+'/h5_files/{0}/*.h5'.format(date_min.strftime('%Y/%b/%d'))
    globber_end = path_to_sort_results+'/h5_files/{0}/*.h5'.format(date_max.strftime('%Y/%b/%d'))
    # print globber_start
    # print globber_end
    h5_filenames = glob.glob(globber_start)
    h5_filenames += glob.glob(globber_end)
    # Keep only the unique filenames
    h5_filenames = set(h5_filenames)

    reduced_h5 = []
    for h5i in h5_filenames:
        h5time, h5endtime = parse_lma_h5_filename(h5i)
        if ( (h5time - date_min).total_seconds() > -60.0 ) & ( (h5time - date_max).total_seconds() < 60.0 ):
            reduced_h5.append(h5i)
    h5_filenames = reduced_h5
    h5_filenames.sort()
    assert len(h5_filenames) > 0
    return h5_filenames

def nc_files_from_standard_path(path_to_sort_results, field_file_code, date_min, date_max):
    filenames = glob.glob(path_to_sort_results+'/grid_files/{1}/*{2}.nc'.format('unused',date_min.strftime('%Y/%b/%d'), field_file_code))
    filenames += glob.glob(path_to_sort_results+'/grid_files/{1}/*{2}.nc'.format('unused',date_max.strftime('%Y/%b/%d'), field_file_code))
    # Keep only the unique filenames
    filenames = set(f for f in filenames if not (f.find('_3d.nc') >=0))

    # nc_names = glob.glob(grid_dir+'/20%s/*.nc' %(date.strftime('%y/%b/%d')))
    # nc_names_3d = glob.glob(grid_dir+'/20%s/*_3d.nc' %(date.strftime('%y/%b/%d')))
    # nc_names_2d = list(set(nc_names) - set(nc_names_3d))

    reduced_grids = []
    for gi in filenames:
        g_parts = os.path.split(gi)[-1].split('_')[1:3]
        g_duration_part = float(os.path.split(gi)[-1].split('_')[3])
        gtime = datetime.strptime('_'.join(g_parts), '%Y%m%d_%H%M%S')
        if ( ( (gtime - date_min).total_seconds() > -g_duration_part ) &
             ( (gtime - date_max).total_seconds() <  g_duration_part ) ):
            reduced_grids.append(gi)
    filenames = reduced_grids
    filenames.sort()
    # print filenames
    assert len(filenames) > 0
    return filenames
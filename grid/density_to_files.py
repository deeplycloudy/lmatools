from __future__ import absolute_import
from __future__ import print_function
import glob
import gc
import numpy as np
from lmatools.stream.subset import coroutine
from lmatools.density_tools import unique_vectors
from six.moves import zip

# -------------------------------------------------------------------------- 
# ----- This section could be replaced with stormdrain.pipeline imports ----
# -------------------------------------------------------------------------- 



            
# class map_projector(object):
#     def __init__(self, ctr_lat, ctr_lon, proj_name='eqc'):
#         self.mapProj = MapProjection(projection=proj_name, ctrLat=ctr_lat, ctrLon=ctr_lon, lat_ts=ctr_lat, lon_0=ctr_lon)
#         self.geoProj = GeographicSystem()
# 
#     def __call__(self, lon, lat, alt):
#         x,y,z = self.mapProj.fromECEF( 
#                 *self.geoProj.toECEF(lon, lat, alt)
#                 )
#         return x, y, z
#         
# @coroutine
# def map_projector(ctr_lat, ctr_lon, target, proj_name='eqc'):
#     mapProj = MapProjection(projection=proj_name, ctrLat=ctr_lat, ctrLon=ctr_lon, lat_ts=ctr_lat, lon_0=ctr_lon)
#     geoProj = GeographicSystem()
#     while True:
#         lon, lat, alt = (yield)
#         x,y,z = self.mapProj.fromECEF( 
#                 *self.geoProj.toECEF(lon, lat, alt)
#                 )
#         target.send((x,y,z))

# -------------------------------------------------------------------------- 
# -------------------------------------------------------------------------- 
# -------------------------------------------------------------------------- 


@coroutine
def flash_count_log(logfile, format_string="%s flashes in frame starting at %s"):
    """ Write flash count for some frame to a file-like object. File open/close should be handled
        by the calling routine."""
        
    # Track flash count for each frame
    frame_times = {}
    try:
        while True:
            # Receive list of flashes, frame start time
            flashes, frame_start_time = (yield)
            n_flashes = len(flashes)
        
            try:
                frame_times[frame_start_time] += n_flashes
            except KeyError:
                # Key doesn't exist, so can't increment flash count 
                frame_times[frame_start_time]  = n_flashes
    except GeneratorExit:
        all_times = list(frame_times.keys())
        all_times.sort()
        for frame_start_time in all_times:
            flash_count_status = format_string % (frame_times[frame_start_time], frame_start_time)
            logfile.write(flash_count_status+'\n')


@coroutine
def filter_flash(target, min_points=10):
    """ Filters flash by minimum number of points.
    """
    while True:
        evs, flash = (yield)           # Receive a flash
        if (flash['n_points'] >= 10):
            target.send((evs, flash))
        del evs, flash

def stack_chopped_arrays(chop_sequence):
    """ Given a sequence of lists of arrays, return an equal length sequence 
        where the arrays have been combined by position in the original sequence.
        The lists of arrays must each be of the same length. This is useful when
        there is a list of arrays corresponding to data subdivided into time
        series chunks.
        
        In the example below, each row is data from a different file (letters)
        and each column is a different time window in a time series. By stacking
        the columns, a combined time series is generated.
        
        ([a0, a1, a2, a3], 
         [b0, b1, b2, b3],
         [c0, c1, c2, c3],)
        becomes
        [a0+b0+c0, a1+b1+c1, a2+b2+c2, a3+b3+b3]
        where plus indicates concatenation        
    """
         
    combined = [np.hstack(a) for a in zip(*chop_sequence)]
    return combined


class ArrayChopper(object):
    """ Initialized with an array of N_+1 edges corresponding to N
        windows. The edges are assumed to be sorted.
    
        Methods
        window_masks(data, edge_key=None): given an array of data with a named dtype,
            return a list of boolean masks that can be used to index data,
            giving the subset of data which corresponds to each window.
            If an edge_key is provided, it is assumed to reference a named array 
            and masking is performed on data[edge_key] 
    
        chop(data, edge_key=None): Returns a list of arrays where the
            masks described above have been applied to chop the data
    
        Generator functions for each of the above are also available
        gen_window_masks, gen_chop
    
    """
    def __init__(self, edges):
        self.edges = edges
        
    def _apply_edge_key(self, data, edge_key):
        if edge_key is not None: 
            d = data[edge_key]
        else:
            d = data
        return d
        
    def gen_edge_pairs(self):
        for l, r in zip(self.edges[:-1], self.edges[1:]):
            yield l, r
        
    def window_masks(self, data, edge_key=None):
        masks = [w for w in self.gen_window_masks(self, data, edge_key)]
        return masks
        
    def gen_window_masks(self, data, edge_key=None):
        d = self._apply_edge_key(data, edge_key)
        for l, r in self.gen_edge_pairs():
            # make sure this is only one-side inclusive to eliminate double-counting
            within = (d >= l) & (d < r)
            yield within
    
    def chop(self, data, edge_key=None):
        chopped = [d for d in self.gen_chop(data, edge_key)]
        return chopped
        
    def gen_chop(self, data, edge_key=None):
        # d = self._apply_edge_key(data, edge_key)
        for mask in self.gen_window_masks(data, edge_key):
            yield data[mask]

@coroutine
def flashes_to_frames(time_edges, targets, time_key='start', time_edges_datetime=None, flash_counter=None):
    """ time_edges_datetime is same len as time_edges but with datetime objects instead of floats.
        allows 
    """
    if time_edges_datetime is None:
        # print "Datetime-style time edges not found, using time edges in seconds for flash count label"
        time_edges_datetime = time_edges
    
    flash_count_messages = []
    
    assert len(time_edges) == (len(time_edges_datetime))
    assert len(time_edges) == (len(targets)+1)
    while True:
        events, flashes = (yield)
        start_times = flashes[time_key]
        sort_idx = np.argsort(start_times) #, order=[time_key])
        idx = np.searchsorted(start_times[sort_idx], time_edges)
        slices = [slice(*i) for i in zip(idx[0:-1], idx[1:])]
        for target, s, frame_start_time in zip(targets, slices, time_edges_datetime[:-1]):
            these_flashes = flashes[sort_idx][s]
            if flash_counter is not None:
                flash_counter.send((these_flashes, frame_start_time))
            # flash_count_status = "Sending %s flashes to frame starting at %s" % (len(these_flashes), frame_start_time)
            # flash_count_messages += flash_count_status
            # print flash_count_status
            target.send((events, these_flashes))
        del events, flashes, start_times, sort_idx, idx, slices
    print(flash_count_messages)

def event_yielder(evs, fls):
    for fl in fls:
        these_events = evs[evs['flash_id'] == fl['flash_id']]
        # if len(these_events) <> fl['n_points']:
            # print 'not giving all ', fl['n_points'], ' events? ', these_events.shape
        for an_ev in these_events:
            yield an_ev


@coroutine
def extract_events_for_flashes(target, flashID_key='flash_id'):
    """ Takes a large table of events and grabs only the events belonging to the flashes.
    """
    
    while True:
        evs, fls = (yield)
        # print 'extracting events'
        # event_dtype = evs[0].dtype
        event_dtype = evs.dtype
        events = np.fromiter( (event_yielder(evs, fls)) , dtype=event_dtype)
        
        # The line below (maybe maybe maybe)
        # events = np.fromiter((evs[evs['flash_id'] == fl['flash_id']] for fl in fls), dtype=event_dtype)
        # does the same thing as the two following lines, but is 10x slower.
        # The 'mapper' could actually be optimized further by calculating it globally, once per events table,
        # but this is fast enough and saves having to pass through another variable.
        # mapper = dict(zip(evs['flash_id'],evs))
        # events = np.fromiter( (mapper[fl['flash_id']] for fl in fls), dtype=event_dtype)

        target.send((events, fls))
        del events, evs, fls

# @coroutine
# def extract_events(target, flashID_key='flash_id'):
#     """ Takes a large table of events and grabs only the events belonging to the flash.
#         This is useful if you filter out a bunch of flashes before going to the trouble of
#         reading the flashes in.
#     """
#     while True:
#         evs, flash = (yield)
#         flash_id = flash[flashID_key]
#         event_dtype = evs[0].dtype
#         # events = [ev[:] for ev in evs if ev[flashID_key] == flash_id]
#         # events = np.asarray(events, dtype=event_dtype)
#         # events = evs[:]
#         events = evs[evs[flashID_key]==flash_id]
#         # events = np.fromiter((ev[:] for ev in evs if ev[flashID_key] == flash_id), dtype=event_dtype)
#         target.send((events, flash))

@coroutine
def no_projection(x_coord, y_coord, z_coord, target, use_flashes=False):
    while True:
        events, flashes = (yield)
        if use_flashes==True:
            points = flashes
        else:
            points = events
        x,y,z = points[x_coord], points[y_coord], points[z_coord]
        target.send((events, flashes, x,y,z))
        del events, flashes, x,y,z, points
    

@coroutine
def project(x_coord, y_coord, z_coord, mapProj, geoProj, target, use_flashes=False):
    """ Adds projected coordinates to the flash and events stream"""
    while True:
        events, flashes = (yield)
        if use_flashes==True:
            points = flashes
        else:
            points = events
        x,y,z = mapProj.fromECEF( 
                *geoProj.toECEF(points[x_coord], points[y_coord], points[z_coord])
                )
        target.send((events, flashes, x,y,z))
        del events, flashes, x,y,z, points

@coroutine
def footprint_mean(flash_id_key='flash_id', area_key='area'):
    """ Takes x, y, z flash locations and gets 
        Extent density unique pixels, average all flashes 
    """
    while True:
        events, flash, x,y,z = (yield)
        # print 'Doing extent density',
        x_i = np.floor( (x-x0)/dx ).astype('int32')
        y_i = np.floor( (y-y0)/dy ).astype('int32')

        if len(x_i) > 0:
            footprints = dict(list(zip(flash[flash_id_key], flash[area_key])))
            # print 'with points numbering', len(x_i)
            unq_idx = unique_vectors(x_i, y_i, events['flash_id'])
            # if x[unq_idx].shape[0] > 1:
            fl_id = events['flash_id'][unq_idx]
            areas = [footprints[fi] for fi in fl_id] #puts areas in same order as x[unq_idx], y[unq_idx]
            # counts normalized by areas 
            target.send((x[unq_idx],y[unq_idx],areas))
            del footprints, unq_idx, fl_id, areas
        # else:
            # print ''
        del events, flash, x, y, z, x_i, y_i

@coroutine
def footprint_mean_3d(flash_id_key='flash_id', area_key='area'):
    """ Takes x, y, z flash locations and gets 
        Extent density unique pixels, average all flashes 
    """
    while True:
        events, flash, x,y,z = (yield)
        # print 'Doing extent density',
        x_i = np.floor( (x-x0)/dx ).astype('int32')
        y_i = np.floor( (y-y0)/dy ).astype('int32')
        z_i = np.floor( (z-z0)/dz ).astype('int32')
        if len(x_i) > 0:
            footprints = dict(list(zip(flash[flash_id_key], flash[area_key])))
            # print 'with points numbering', len(x_i)
            unq_idx = unique_vectors(x_i, y_i, z_i, events['flash_id'])
            # if x[unq_idx].shape[0] > 1:
            fl_id = events['flash_id'][unq_idx]
            areas = [footprints[fi] for fi in fl_id] #puts areas in same order as x[unq_idx], y[unq_idx]
            # counts normalized by areas 
            target.send((x[unq_idx],y[unq_idx],z[unq_idx],areas))
            del footprints, unq_idx, fl_id, areas
        # else:
            # print ''
        del events, flash, x, y, z, x_i, y_i, z_i


@coroutine
def point_density(target):
    """ Sends event x, y, z location directly
    """
    while True:
        events, flash, x, y, z = (yield)
        # print 'Doing point density',
        if len(x) > 0:
            print('with points numbering', len(x))
            target.send((x, y, None))
        del events, flash ,x,y,z

@coroutine
def point_density_3d(target):
    """ Sends event x, y, z location directly
    """
    while True:
        events, flash, x, y, z = (yield)
        # print 'Doing point density',
        if len(x) > 0:
            print('with points numbering', len(x))
            target.send((x, y, z, None))
        del events, flash ,x,y,z

@coroutine
def flash_std(x0, y0, dx, dy, target, flash_id_key='flash_id', weight_key=None):
    """ This function assumes a regular grid in x and y with spacing dx, dy
        
        x0, y0 is the x coordinate of the lower left corner of the lower-left grid cell, 
        i.e., the lower left node of the grid mesh in cartesian space
        
        Eliminates duplicate points in gridded space and sends the reduced
        set of points to the target.
    
        NOTE: Use of this function is to only find the standard deviation of flash size.
    """
    while True:
        # assumes x,y,z are in same order as events
        events, flash, x,y,z = (yield)
        # print 'Doing extent density',
        x_i = np.floor( (x-x0)/dx ).astype('int32')
        y_i = np.floor( (y-y0)/dy ).astype('int32')

        if len(x_i) > 0:
            print('extent with points numbering', len(x_i), ' with weights', weight_key)
            unq_idx = unique_vectors(x_i, y_i, events[flash_id_key])
            # if x[unq_idx].shape[0] > 1:
            if weight_key != None:
                weight_lookup = dict(list(zip(flash[flash_id_key], flash[weight_key]**2.)))
                weights = [weight_lookup[fi] for fi in events[unq_idx]['flash_id']] #puts weights in same order as x[unq_idx], y[unq_idx]
                del weight_lookup
            else:
                weights = None
            target.send((x[unq_idx], y[unq_idx], weights))
            del weights, unq_idx
        # else:
            # print ''
        del events, flash, x, y, z, x_i, y_i
        
@coroutine
def flash_std_3d(x0, y0, z0, dx, dy, dz, target, flash_id_key='flash_id', weight_key=None):
    """ This function assumes a regular grid in x and y with spacing dx, dy
        
        x0, y0 is the x coordinate of the lower left corner of the lower-left grid cell, 
        i.e., the lower left node of the grid mesh in cartesian space
        
        Eliminates duplicate points in gridded space and sends the reduced
        set of points to the target.
    """
    while True:
        # assumes x,y,z are in same order as events
        events, flash, x,y,z = (yield)
        # print('Doing extent density',)
        x_i = np.floor( (x-x0)/dx ).astype('int32')
        y_i = np.floor( (y-y0)/dy ).astype('int32')
        z_i = np.floor( (z-z0)/dz ).astype('int32')
        print(len(x_i))
        if len(x_i) > 0:
            print('extent with points numbering', len(x_i), ' with weights', weight_key)
            unq_idx = unique_vectors(x_i, y_i, z_i, events[flash_id_key])
            # if x[unq_idx].shape[0] > 1:
            if weight_key != None:
                weight_lookup = dict(list(zip(flash[flash_id_key], flash[weight_key]**2.)))
                weights = [weight_lookup[fi] for fi in events[unq_idx]['flash_id']] #puts weights in same order as x[unq_idx], y[unq_idx]
                del weight_lookup
            else:
                weights = None
            target.send((x[unq_idx], y[unq_idx], z[unq_idx], weights)) 
            del weights, unq_idx
        # else:
            # print ''
        del events, flash, x, y, z, x_i, y_i, z_i


@coroutine
def extent_density(x0, y0, dx, dy, target, flash_id_key='flash_id', weight_key=None):
    """ This function assumes a regular grid in x and y with spacing dx, dy
        
        x0, y0 is the x coordinate of the lower left corner of the lower-left grid cell, 
        i.e., the lower left node of the grid mesh in cartesian space
        
        Eliminates duplicate points in gridded space and sends the reduced
        set of points to the target.
    """
    while True:
        # assumes x,y,z are in same order as events
        events, flash, x,y,z = (yield)
        # print 'Doing extent density',
        x_i = np.floor( (x-x0)/dx ).astype('int32')
        y_i = np.floor( (y-y0)/dy ).astype('int32')

        if len(x_i) > 0:
            print('extent with points numbering', len(x_i), ' with weights', weight_key)
            unq_idx = unique_vectors(x_i, y_i, events[flash_id_key])
            # if x[unq_idx].shape[0] > 1:
            if weight_key != None:
                weight_lookup = dict(list(zip(flash[flash_id_key], flash[weight_key])))
                weights = [weight_lookup[fi] for fi in events[unq_idx]['flash_id']] #puts weights in same order as x[unq_idx], y[unq_idx]
                del weight_lookup
            else:
                weights = None
                

            target.send((x[unq_idx], y[unq_idx], weights))
            del weights, unq_idx
        # else:
            # print ''
        del events, flash, x, y, z, x_i, y_i


@coroutine
def extent_density_3d(x0, y0, z0, dx, dy, dz, target, flash_id_key='flash_id', weight_key=None):
    """ This function assumes a regular grid in x and y with spacing dx, dy
        
        x0, y0 is the x coordinate of the lower left corner of the lower-left grid cell, 
        i.e., the lower left node of the grid mesh in cartesian space
        
        Eliminates duplicate points in gridded space and sends the reduced
        set of points to the target.
    """
    while True:
        # assumes x,y,z are in same order as events
        events, flash, x,y,z = (yield)
        # print('Doing extent density',)
        x_i = np.floor( (x-x0)/dx ).astype('int32')
        y_i = np.floor( (y-y0)/dy ).astype('int32')
        z_i = np.floor( (z-z0)/dz ).astype('int32')
        print(len(x_i))
        if len(x_i) > 0:
            print('extent with points numbering', len(x_i), ' with weights', weight_key)
            unq_idx = unique_vectors(x_i, y_i, z_i, events[flash_id_key])
            # if x[unq_idx].shape[0] > 1:
            if weight_key != None:
                weight_lookup = dict(list(zip(flash[flash_id_key], flash[weight_key])))
                weights = [weight_lookup[fi] for fi in events[unq_idx]['flash_id']] #puts weights in same order as x[unq_idx], y[unq_idx]
                del weight_lookup
            else:
                weights = None
            target.send((x[unq_idx], y[unq_idx], z[unq_idx], weights)) 
            del weights, unq_idx
        # else:
            # print ''
        del events, flash, x, y, z, x_i, y_i, z_i


@coroutine
def accumulate_points_on_grid(grid, xedge, yedge, out=None, label=''):
    assert xedge.shape[0] == grid.shape[0]+1
    assert yedge.shape[0] == grid.shape[1]+1
    if out == None:
        out = {}
    # grid = None
    try:
        while True:
            x, y, weights = (yield)
            if len(x) > 0:
                x = np.atleast_1d(x)
                y = np.atleast_1d(y)
                print('accumulating ', len(x), 'points for ', label)
                count, edges = np.histogramdd((x,y), bins=(xedge, yedge), weights=None, normed=False)    
                
                if weights != None:
                    # histogramdd sums up weights in each bin for normed=False
                    total, edges = np.histogramdd((x,y), bins=(xedge, yedge), weights=weights, normed=False)
                    # return the mean of the weights in each bin
                    bad = (count <= 0)
                    count = np.asarray(total, dtype='float32')/count
                    count[bad] = 0.0 
                    
                    del total, edges, bad
                    
                # try:
                    # count, edges = np.histogramdd((x,y), bins=(xedge, yedge), weights=weights)
                # except AttributeError:
                #     # if x,y are each scalars, need to make 1D arrays
                #     x = np.asarray((x,))
                #     y = np.asarray((y,))
                #     count, edges = np.histogramdd((x,y), bins=(xedge, yedge), weights=weights)
                # using += (as opposed to grid = grid + count) is essential
                # so that the user can retain a reference to the grid object
                # outside this routine.
                if grid is None:
                    grid = count
                    out['out'] = grid
                else:
                    grid += count.astype(grid.dtype)
                del count
            del x, y, weights
            gc.collect()
    except GeneratorExit:
        out['out'] = grid

@coroutine
def accumulate_points_on_grid_3d(grid, xedge, yedge, zedge, out=None, label=''):
    assert xedge.shape[0] == grid.shape[0]+1
    assert yedge.shape[0] == grid.shape[1]+1
    assert zedge.shape[0] == grid.shape[2]+1
    if out == None:
        out = {}
    # grid = None
    try:
        while True:
            x, y, z, weights = (yield)
            if len(x) > 0:
                x = np.atleast_1d(x)
                y = np.atleast_1d(y)
                z = np.atleast_1d(z)
                print('accumulating ', len(x), 'points for ', label)
                count, edges = np.histogramdd((x,y,z), bins=(xedge, yedge, zedge), weights=None, normed=False)    
                
                if weights != None:
                    # histogramdd sums up weights in each bin for normed=False
                    total, edges = np.histogramdd((x,y,z), bins=(xedge, yedge, zedge), weights=weights, normed=False)
                    # return the mean of the weights in each bin
                    bad = (count <= 0)
                    count = np.asarray(total, dtype='float32')/count
                    count[bad] = 0.0 
                    
                    del total, edges, bad
        
                # using += (as opposed to grid = grid + count) is essential
                # so that the user can retain a reference to the grid object
                # outside this routine.
                if grid is None:
                    grid = count
                    out['out'] = grid
                else:
                    grid += count.astype(grid.dtype)
                del count
            del x, y, z, weights
            gc.collect()
    except GeneratorExit:
        out['out'] = grid
        
##Repetition of functions below can probably be reduced to the original functions, however
##remain as this was the only way to get new gridded fields that were no the mean.
        
####FOR STANDARD DEVIATION OF A SINGLE FIELD:
@coroutine
def accumulate_points_on_grid_sdev(grid, grid2, xedge, yedge, out=None, label=''):
    assert xedge.shape[0] == grid.shape[0]+1
    assert yedge.shape[0] == grid.shape[1]+1
    if out == None:
        out = {}
    # grid = None
    try:
        while True:
            x, y, weights = (yield)
            if len(x) > 0:
                x = np.atleast_1d(x)
                y = np.atleast_1d(y)
                print('accumulating ', len(x), 'points for ', label)
                count, edges = np.histogramdd((x,y), bins=(xedge, yedge), weights=None, normed=False)    
                
                if weights != None:
                    # histogramdd sums up weights in each bin for normed=False
                    total, edges = np.histogramdd((x,y), bins=(xedge, yedge), weights=np.asarray(weights)**2., normed=False)
                    
                    # return the mean of the weights in each bin
                    bad = (count <= 0)
                    count = np.asarray(total, dtype='float32')/count 
                    count[bad] = 0.0 
                    
                    del total, edges, bad
                    
                # using += (as opposed to grid = grid + count) is essential
                # so that the user can retain a reference to the grid object
                # outside this routine.
                if grid is None:
                    grid = count
                    out['out'] = grid
                else:
                    grid += count.astype(grid.dtype) 
                    grid = np.sqrt(grid - (grid2)**2.)
                del count
            del x, y, weights
            gc.collect()
    except GeneratorExit:
        out['out'] = grid

@coroutine
def accumulate_points_on_grid_sdev_3d(grid, grid2, xedge, yedge, zedge, out=None, label=''):
    assert xedge.shape[0] == grid.shape[0]+1
    assert yedge.shape[0] == grid.shape[1]+1
    assert zedge.shape[0] == grid.shape[2]+1
    if out == None:
        out = {}
    # grid = None
    try:
        while True:
            x, y, z, weights = (yield)
            if len(x) > 0:
                x = np.atleast_1d(x)
                y = np.atleast_1d(y)
                z = np.atleast_1d(z)
                print('accumulating ', len(x), 'points for ', label)
                count, edges = np.histogramdd((x,y,z), bins=(xedge, yedge, zedge), weights=None, normed=False)    
                
                if weights != None:
                    # histogramdd sums up weights in each bin for normed=False
                    total, edges = np.histogramdd((x,y,z), bins=(xedge, yedge, zedge), weights=np.asarray(weights)**2., normed=False)
                    
                    # return the mean of the weights in each bin
                    bad = (count <= 0)
                    count = np.asarray(total, dtype='float32')/count 
                    count[bad] = 0.0 
                    
                    del total, edges, bad
                    
                # using += (as opposed to grid = grid + count) is essential
                # so that the user can retain a reference to the grid object
                # outside this routine.
                if grid is None:
                    grid = count
                    out['out'] = grid
                else:
                    grid += count.astype(grid.dtype)
                    grid = np.sqrt(grid - (grid2)**2.)
                del count
            del x, y, z, weights
            gc.collect()
    except GeneratorExit:
        out['out'] = grid
        
#####FOR TOTAL ENERGY:
@coroutine
def accumulate_energy_on_grid(grid, xedge, yedge, out=None, label=''):
    assert xedge.shape[0] == grid.shape[0]+1
    assert yedge.shape[0] == grid.shape[1]+1
    if out == None:
        out = {}
    # grid = None
    try:
        while True:
            x, y, weights = (yield)
            if len(x) > 0:
                x = np.atleast_1d(x)
                y = np.atleast_1d(y)
                print('accumulating ', len(x), 'points for ', label)
                count, edges = np.histogramdd((x,y), bins=(xedge, yedge), weights=None, normed=False)    
                
                if weights != None:
                    # histogramdd sums up weights in each bin for normed=False
                    total, edges = np.histogramdd((x,y), bins=(xedge, yedge), weights=np.abs(weights), normed=False)
                    # return the mean of the weights in each bin
                    bad = (count <= 0)
                  
                    count = np.asarray(total, dtype='float32')#/count 
                    count[bad] = 0.0 
                    
                    del total, edges, bad
    
                # using += (as opposed to grid = grid + count) is essential
                # so that the user can retain a reference to the grid object
                # outside this routine.
                if grid is None:
                    grid = count
                    out['out'] = grid
                else:
                    grid += count.astype(grid.dtype)
                del count
            del x, y, weights
            gc.collect()
    except GeneratorExit:
        out['out'] = grid

@coroutine
def accumulate_energy_on_grid_3d(grid, xedge, yedge, zedge, out=None, label=''):
    assert xedge.shape[0] == grid.shape[0]+1
    assert yedge.shape[0] == grid.shape[1]+1
    assert zedge.shape[0] == grid.shape[2]+1
    if out == None:
        out = {}
    # grid = None
    try:
        while True:
            x, y, z, weights = (yield)
            if len(x) > 0:
                x = np.atleast_1d(x)
                y = np.atleast_1d(y)
                z = np.atleast_1d(z)
                print('accumulating ', len(x), 'points for ', label)
                count, edges = np.histogramdd((x,y,z), bins=(xedge, yedge, zedge), weights=None, normed=False)    
                
                if weights != None:
                    # histogramdd sums up weights in each bin for normed=False
                    total, edges = np.histogramdd((x,y,z), bins=(xedge, yedge, zedge), weights=np.abs(weights), normed=False)
                    
                    # return the mean of the weights in each bin
                    bad = (count <= 0)
                    count = np.asarray(total, dtype='float32')#/count 
                    count[bad] = 0.0 
                    
                    del total, edges, bad
                    
                # try:
                    # count, edges = np.histogramdd((x,y), bins=(xedge, yedge), weights=weights)
                # except AttributeError:
                #     # if x,y are each scalars, need to make 1D arrays
                #     x = np.asarray((x,))
                #     y = np.asarray((y,))
                #     count, edges = np.histogramdd((x,y), bins=(xedge, yedge), weights=weights)
                # using += (as opposed to grid = grid + count) is essential
                # so that the user can retain a reference to the grid object
                # outside this routine.
                if grid is None:
                    grid = count
                    out['out'] = grid
                else:
                    grid += count.astype(grid.dtype)
                del count
            del x, y, z, weights
            gc.collect()
    except GeneratorExit:
        out['out'] = grid
        
# if __name__ == '__main__':
#     do_profile=False
#     if do_profile:
#         import hotshot
#         from hotshot import stats
#         prof = hotshot.Profile("density_test_profile")
#         prof.runcall(example)
#         prof.close()
#         s=stats.load("density_test_profile")
#         s.sort_stats("time").print_stats()
#     else:
#         x_coord, y_coord, lons, lats, test_grid = example()

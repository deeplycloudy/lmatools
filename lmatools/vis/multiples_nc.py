from __future__ import absolute_import
from __future__ import print_function
from .small_multiples import small_multiples_plot

import os
import gc

from datetime import datetime, timedelta
import numpy as np

try:
    from netCDF4 import Dataset as NetCDFFile
except ImportError:
    from scipy.io.netcdf import NetCDFFile

# from pylab import figure, get_cmap, colorbar
import matplotlib
from matplotlib.figure import figaspect, Figure
from matplotlib.colorbar import ColorbarBase
from matplotlib.cm import get_cmap
from matplotlib.ticker import FuncFormatter
from matplotlib.dates import date2num, DateFormatter
from matplotlib.backends.backend_agg import FigureCanvasAgg  

from math import ceil

#import pytz
#tz=pytz.timezone('US/Eastern') # Why, oh, why, is it using my local time zone?
#time_series_x_fmt = DateFormatter('%H%M:%S', tz=tz)


def kilo(x, pos):
    'The two args are the value and tick position'
    return '%s' % (x/1000.0)

kilo_formatter = FuncFormatter(kilo)

# # /* DC LMA */
# DClat = 38.8888500 # falls church / western tip of arlington, rough centroid of stations
# DClon = -77.1685800
# kounLat = 35.23833
# kounLon = -97.46028
# kounAlt = 377.0
# radarLat=kounLat
# radarLon=kounLon
# radarAlt=0.0

# # ARMOR
# radarLat = 34.6461
# radarLon = -86.7714
# radarAlt = 180.0

# mapProj = MapProjection(projection='eqc', ctrLat=radarLat, ctrLon=radarLon, lat_ts=radarLat, lon_0=radarLon)
# geoProj = GeographicSystem()

def centers_to_edges(x):
    xedge=np.zeros(x.shape[0]+1)
    xedge[1:-1] = (x[:-1] + x[1:])/2.0
    dx = np.mean(np.abs(xedge[2:-1] - xedge[1:-2]))
    xedge[0] = xedge[1] - dx
    xedge[-1] = xedge[-2] + dx
    return xedge        


def multiples_figaspect(n_rows, n_cols, x_range, y_range, fig_width=None, max_height=None):
    """ 
        DIMENSIONS, assuming no margins anywhere.
        Plottable area: w, h (inches)
        Data units: x_range, y_range (km, array shape, etc.)
        Subaxis count: n_cols, n_rows (laid out in a uniformly spaced grid)
        Subaxis size: x_box, y_box (inches)
        
        y_box = h / n_row
        x_box = w / n_col
        
        Assume we're working with a plottable area of fig_width (no margins).
        Also assume we want squared-up units (e.g, square 
        grid pixels if x_range, y_range are nx and ny for an array, or physical
        units, so 1 km north = 1 km east). In that case, we want:
            
            y_range / y_box  =  x_range / x_box  (e.g., km/in)
            
        So, given w = fig_width (taken from rcParams by default):
            
            (y_range * n_rows) / h  =  (x_range * n_col) / w
            
            or 
            
            h = (y_range * n_rows * w) / (x_range * n_col)
            
        If we want to calculate the number of rows that can fit on a page with max_height,
            
            y_box = h / n_rows
            n_rows_perpage = floor(max_height / y_box)
            
        Returns: w, h, n_rows_perpage, n_pages
        
    """
    if fig_width is None:
        fig_width, fig_height = matplotlib.rcParams['figure.figsize'][0:1]
    w = float(fig_width)
    
    n_rows, n_cols, x_range, y_range = (float(i) for i in (n_rows, n_cols, x_range, y_range))
    
    h = (y_range * n_rows * w) / (x_range * n_cols)
    
    if max_height is not None:
        max_height = float(max_height)
        y_box = h / n_row
        n_rows_perpage = floor(max_height / y_box)
        n_pages = ceil(n_rows/n_rows_perpage)
    else:
        n_rows_perpage = n_rows
        n_pages = 1
        
    return w, h, n_rows_perpage, n_pages

def read_file_3d(filename, grid_name, x_name='x', y_name='y', z_name = 'z', t_name='time'):
    """ colormap: a string giving the name of a matplotlib built-in colormap, 
            or a matplotlib.colors.Colormap instance
        """
    
    f = NetCDFFile(filename)
    data = f.variables  # dictionary of variable names to nc_var objects
    dims = f.dimensions # dictionary of dimension names to sizes
    x_idx = data[x_name]
    y_idx = data[y_name]
    z_idx = data[z_name]
    x = x_idx[:]
    y = y_idx[:]
    z = z_idx[:]
    t = data[t_name]
    all_z = z[:]
    grid = data[grid_name]
    
    assert len(x.shape) == 1
    assert len(y.shape) == 1
    assert len(z.shape) == 1
    assert len(t.shape) == 1
      
    grid_dims = grid.dimensions # tuple of dimension names
    name_to_idx = dict((k, i) for i, k in enumerate(grid_dims))
    grid_t_idx = name_to_idx[t.dimensions[0]]
    grid_x_idx = name_to_idx[x_idx.dimensions[0]]
    grid_y_idx = name_to_idx[y_idx.dimensions[0]]
    grid_z_idx = name_to_idx[z_idx.dimensions[0]]
        
    return grid, grid_name, x, y, all_z, t, grid_t_idx, grid_x_idx, grid_z_idx


def make_plot_3d(grid, grid_name, x, y, all_z, t, grid_t_idx, grid_x_idx, grid_z_idx, n_cols = 6,
                 outpath='', filename_prefix='LMA',do_save=True, 
                 image_type='pdf', colormap='cubehelix' , grid_range=None): 

    n_frames = t.shape[0]
    # n_cols = 6
    n_rows = int(ceil( float(n_frames) / n_cols ))
  
    if type(colormap) == type(''):
        colormap = get_cmap(colormap)
    grey_color = (0.5,)*3
    frame_color = (0.2,)*3
    
    density_maxes = []
    total_counts = []
    all_t = []
   
    xedge = centers_to_edges(x)
    x_range = xedge.max() - xedge.min()
    yedge = centers_to_edges(y)
    y_range = yedge.max() - yedge.min()
    dx = (xedge[1]-xedge[0])

    # count_scale_factor = dx # / 1000.0
    # max_count_baseline = 450 * count_scale_factor #/ 10.0
    min_count, max_count = 1, grid[:].max() #max_count_baseline*(t[1]-t[0])
    if (max_count == 0) | (max_count == 1 ):
        max_count = min_count+1

    default_vmin = -0.2
    if np.log10(max_count) <= default_vmin:
        vmin_count = np.log10(max_count) + default_vmin
    else:
        vmin_count = default_vmin

    # If the range of values for the colorbar is manually specified, 
    # overwrite what we did above
    if grid_range is not None:
        vmin_count = grid_range[0]
        max_count = grid_range[1]
            
    # w, h = figaspect(float(n_rows)/n_cols) # breaks for large numbers of frames - has a hard-coded max figure size
    w, h, n_rows_perpage, n_pages = multiples_figaspect(n_rows, n_cols, x_range, y_range, fig_width=8.5, max_height=None)
    fig = Figure(figsize=(w,h))
    canvas = FigureCanvasAgg(fig)
    fig.set_canvas(canvas)
    p = small_multiples_plot(fig=fig, rows=n_rows, columns=n_cols)
    p.label_edges(True)
    pad = 0.0 # for time labels in each frame

    for ax in p.multiples.flat:
        if int(matplotlib.__version__[0]) <= 1:
            ax.set_axis_bgcolor('white')
        else:
            ax.set_facecolor('white')
        ax.spines['top'].set_edgecolor(frame_color)
        ax.spines['bottom'].set_edgecolor(frame_color)
        ax.spines['left'].set_edgecolor(frame_color)
        ax.spines['right'].set_edgecolor(frame_color)
    #   ax.yaxis.set_major_formatter(kilo_formatter)
    #   ax.xaxis.set_major_formatter(kilo_formatter)
    base_date = datetime.strptime(t.units, "seconds since %Y-%m-%d %H:%M:%S")
    time_delta = timedelta(0,float(t[0]),0)
    start_time = base_date + time_delta
    for zi in range(len(all_z)):
        indexer = [slice(None),]*len(grid.shape)
                                
        frame_start_times = []
          
        altitude = all_z[zi]
        for i in range(n_frames):
            p.multiples.flat[i].clear()   # reset (clear) the axes
            frame_start = base_date + timedelta(0,float(t[i]),0)
            frame_start_times.append(frame_start)
            indexer[grid_t_idx] = i
            indexer[grid_z_idx] = zi
            density = grid[indexer]
            density = np.ma.masked_where(density<=0.0, density) # mask grids 0 grids to reveal background color
             
            # density,edges = np.histogramdd((x,y), bins=(xedge,yedge))
            density_plot  = p.multiples.flat[i].pcolormesh(xedge,yedge,
                                       np.log10(density.transpose()),
                                       vmin=vmin_count,
                                       vmax=np.log10(max_count),
                                       cmap=colormap)
            label_string = frame_start.strftime('%H%M:%S')
            x_lab = xedge[0]-pad+x_range*.015
            y_lab = yedge[0]-pad+y_range*.015
            text_label = p.multiples.flat[i].text(x_lab, y_lab, label_string, color=grey_color, size=6)
            density_plot.set_rasterized(True)
            density_maxes.append(density.max())
            total_counts.append(density.sum())
            all_t.append(frame_start)
            print(label_string, x_lab, y_lab, grid_name, density.max(), density.sum())

        color_scale = ColorbarBase(p.colorbar_ax, cmap=density_plot.cmap,
                                           norm=density_plot.norm,
                                           orientation='horizontal')

        # color_scale.set_label('count per pixel')
        color_scale.set_label('log10(count per pixel)')

        view_x = (xedge.min(), xedge.max())
        view_y = (yedge.min(), yedge.max())

        print('making multiples')
        p.multiples.flat[0].axis(view_x+view_y)
        filename = '%s-%s_%s_%05.2fkm_%04.1fkm_%05.1fs.%s' % (filename_prefix, grid_name, start_time.strftime('%Y%m%d_%H%M%S'), dx, altitude, time_delta.seconds, image_type)
        filename = os.path.join(outpath, filename)
        if do_save:
            fig.savefig(filename, dpi=150)
        print('done ', zi)
    return fig, p, frame_start_times, filename
    print(' ... done with loop')

    
def make_plot(filename, grid_name, x_name='x', y_name='y', t_name='time',
                n_cols=6, outpath='', filename_prefix='LMA', 
                do_save=True, image_type='pdf', colormap='gist_earth'):
    """ colormap: a string giving the name of a matplotlib built-in colormap, 
            or a matplotlib.colors.Colormap instance
        """

    f = NetCDFFile(filename)
    data = f.variables  # dictionary of variable names to nc_var objects
    dims = f.dimensions # dictionary of dimension names to sizes
    x = data[x_name]
    y = data[y_name]
    t = data[t_name]
    grid = data[grid_name]
    
    assert len(x.shape) == 1
    assert len(y.shape) == 1
    assert len(t.shape) == 1
        
    grid_dims = grid.dimensions # tuple of dimension names
    name_to_idx = dict((k, i) for i, k in enumerate(grid_dims))
    
    grid_t_idx = name_to_idx[t.dimensions[0]]
    grid_x_idx = name_to_idx[x.dimensions[0]]
    grid_y_idx = name_to_idx[y.dimensions[0]]
        
    n_frames = t.shape[0]
    # n_cols = 6
    n_rows = int(ceil( float(n_frames) / n_cols ))
    
    if type(colormap) == type(''):
        colormap = get_cmap(colormap)
    grey_color = (0.5,)*3
    frame_color = (0.2,)*3
    
    density_maxes = []
    total_counts = []
    all_t = []
    
    xedge = centers_to_edges(x)
    x_range = xedge.max() - xedge.min()
    yedge = centers_to_edges(y)
    y_range = yedge.max() - yedge.min()
    dx = (xedge[1]-xedge[0])
    
    # w, h = figaspect(float(n_rows)/n_cols) # breaks for large numbers of frames - has a hard-coded max figure size
    w, h, n_rows_perpage, n_pages = multiples_figaspect(n_rows, n_cols, x_range, y_range, fig_width=8.5, max_height=None)
    
    # count_scale_factor = dx # / 1000.0
    # max_count_baseline = 450 * count_scale_factor #/ 10.0
    min_count, max_count = 1, grid[:].max() #max_count_baseline*(t[1]-t[0])
    if (max_count == 0) | (max_count == 1 ):
        max_count = min_count+1

    default_vmin = -1.0#0.2
    if np.log10(max_count) <= default_vmin:
        vmin_count = np.log10(max_count) + default_vmin
    else:
        vmin_count = default_vmin
    
    fig = Figure(figsize=(w,h))
    canvas = FigureCanvasAgg(fig)
    fig.set_canvas(canvas)
    p = small_multiples_plot(fig=fig, rows=n_rows, columns=n_cols)
    p.label_edges(True)
    pad = 0.0 # for time labels in each frame
    
    for ax in p.multiples.flat:
        if int(matplotlib.__version__[0]) <= 1:
            ax.set_axis_bgcolor('black')
        else:
            ax.set_facecolor('black')
        ax.spines['top'].set_edgecolor(frame_color)
        ax.spines['bottom'].set_edgecolor(frame_color)
        ax.spines['left'].set_edgecolor(frame_color)
        ax.spines['right'].set_edgecolor(frame_color)
    #     ax.yaxis.set_major_formatter(kilo_formatter)
    #     ax.xaxis.set_major_formatter(kilo_formatter)
    base_date = datetime.strptime(t.units, "seconds since %Y-%m-%d %H:%M:%S")
    time_delta = timedelta(0,float(t[0]),0)
    start_time = base_date + time_delta
        
    indexer = [slice(None),]*len(grid.shape)
    
    
    frame_start_times = []
    for i in range(n_frames):
        frame_start = base_date + timedelta(0,float(t[i]),0)
        frame_start_times.append(frame_start)
        indexer[grid_t_idx] = i
        
        density = grid[indexer]
        
        # density,edges = np.histogramdd((x,y), bins=(xedge,yedge))
        density_plot  = p.multiples.flat[i].pcolormesh(xedge,yedge,
                                   np.log10(density.transpose()), 
                                   vmin=vmin_count,
                                   vmax=np.log10(max_count),
                                   cmap=colormap)
        label_string = frame_start.strftime('%H%M:%S')
        text_label = p.multiples.flat[i].text(xedge[0]-pad+x_range*.015, yedge[0]-pad+y_range*.015, label_string, color=grey_color, size=6)
        density_plot.set_rasterized(True)
        density_maxes.append(density.max())
        total_counts.append(density.sum())
        all_t.append(frame_start)
        print(label_string, x.shape, density.max(), density.sum())
        
    color_scale = ColorbarBase(p.colorbar_ax, cmap=density_plot.cmap,
                                       norm=density_plot.norm,
                                       orientation='horizontal')
    
    # color_scale.set_label('count per pixel')
    color_scale.set_label('log10(count per pixel)')
    
    view_x = (xedge.min(), xedge.max())
    view_y = (yedge.min(), yedge.max())
    
    print('making multiples', end=' ')
    p.multiples.flat[0].axis(view_x+view_y)
    filename = '%s-%s_%s_%05.2fkm_%05.1fs.%s' % (filename_prefix, grid_name, start_time.strftime('%Y%m%d_%H%M%S'), dx, time_delta.seconds, image_type)
    filename = os.path.join(outpath, filename)
    if do_save:
        fig.savefig(filename, dpi=150)
    
    return fig, p, frame_start_times, filename
    
    print(' ... done')
        

        
if __name__ == '__main__':
    do_profile=False
    if do_profile:
        import hotshot
        from hotshot import stats
        prof = hotshot.Profile("multiples_test_profile")
        prof.runcall(runtest)
        prof.close()
        s=stats.load("multiples_test_profile")
        s.sort_stats("time").print_stats()
    else:
        import sys
        res = runtest(sys.argv[1], sys.argv[2])

from small_multiples import small_multiples_plot

from datetime import datetime, timedelta
import numpy as np
import pupynere as nc

from pylab import figure, get_cmap, colorbar
from matplotlib.figure import figaspect
from matplotlib.colorbar import ColorbarBase
from matplotlib.ticker import FuncFormatter
from matplotlib.dates import mx2num, date2num, DateFormatter

from math import ceil

import pytz
tz=pytz.timezone('US/Eastern') # Why, oh, why, is it using my local time zone?
time_series_x_fmt = DateFormatter('%H%M:%S', tz=tz)


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



def make_plot(filename, grid_name, x_name='x', y_name='y', t_name='time', n_cols=6):
    
    f = nc.NetCDFFile(filename)
    data = f.variables  # dictionary of variable names to nc_var objects
    dims = f.dimensions # dictionary of dimension names to sizes
    x = data[y_name]
    y = data[x_name]
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
    
    w, h = figaspect(float(n_rows)/n_cols)

    colormap = get_cmap('gist_earth')
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

    f = figure(figsize=(w,h))
    p = small_multiples_plot(fig=f, rows=n_rows, columns=n_cols)
    p.label_edges(True)
    pad = 0.0 # for time labels in each frame
            
    for ax in p.multiples.flat:
        ax.set_axis_bgcolor('black')
        ax.spines['top'].set_edgecolor(frame_color)
        ax.spines['bottom'].set_edgecolor(frame_color)
        ax.spines['left'].set_edgecolor(frame_color)
        ax.spines['right'].set_edgecolor(frame_color)
    #     ax.yaxis.set_major_formatter(kilo_formatter)
    #     ax.xaxis.set_major_formatter(kilo_formatter)
    base_date = datetime.strptime(t.units, "seconds since %Y-%m-%d %H:%M:%S")
    time_delta = timedelta(0,float(t[1]-t[0]),0)
    start_time = base_date + time_delta
        
    indexer = [None,]*len(grid.shape)

    for i in range(n_frames):
        frame_start = base_date + timedelta(0,float(t[i]),0)
        indexer[grid_t_idx] = i
        
        density = grid[indexer]
        
        # density,edges = np.histogramdd((x,y), bins=(xedge,yedge))
        density_plot  = p.multiples.flat[i].pcolormesh(xedge,yedge,
                                   np.log10(density.transpose()), 
                                   vmin=-0.2,
                                   vmax=np.log10(max_count),
                                   cmap=colormap)
        label_string = frame_start.strftime('%H%M:%S')
        text_label = p.multiples.flat[i].text(xedge[0]-pad+x_range*.015, yedge[0]-pad+y_range*.015, label_string, color=grey_color, size=6)
        density_plot.set_rasterized(True)
        density_maxes.append(density.max())
        total_counts.append(density.sum())
        all_t.append(frame_start)
        print label_string, x.shape, density.max(), density.sum()

    color_scale = ColorbarBase(p.colorbar_ax, cmap=density_plot.cmap,
                                       norm=density_plot.norm,
                                       orientation='horizontal')

    # color_scale.set_label('count per pixel')
    color_scale.set_label('log10(count per pixel)')

    view_x = (xedge.min(), xedge.max())
    view_y = (yedge.min(), yedge.max())
    
    print 'making multiples',
    p.multiples.flat[0].axis(view_x+view_y)
    filename = 'LMA-%s_%s_%5.2fkm_%5.1fs.pdf' % (grid_name, start_time.strftime('%Y%m%d_%H%M%S'), dx, time_delta.seconds)
    f.savefig(filename, dpi=150)
    print ' ... done'
    f.clf()
        

        
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
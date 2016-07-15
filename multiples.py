from __future__ import absolute_import
from __future__ import print_function
from .small_multiples import small_multiples_plot
from acuity.data import Bounds, Data, DataTransform, DataSelection
from acuity.coordinateSystems import MapProjection, GeographicSystem
from acuity.LMA.LMAdataHDF import LMAdataManagerHDF
# from acuity.LMA.LMAarrayFile import LMAdataFile
from acuity.LMA.LMAdata import LMAdataManager
from acuity.views import AcuityView

from .density_tools import extent_density

from mx.DateTime import DateTime, DateTimeDelta
import numpy as np
from pylab import figure, get_cmap, colorbar
from matplotlib.figure import figaspect
from matplotlib.colorbar import ColorbarBase
from matplotlib.ticker import FuncFormatter
from matplotlib.dates import mx2num, date2num, DateFormatter

from math import ceil

import pytz
from six.moves import range
from six.moves import zip
tz=pytz.timezone('US/Eastern') # Why, oh, why, is it using my local time zone?
time_series_x_fmt = DateFormatter('%H%M:%S', tz=tz)


def kilo(x, pos):
    'The two args are the value and tick position'
    return '%s' % (x/1000.0)

kilo_formatter = FuncFormatter(kilo)


loadLMA = True

# /* DC LMA */
DClat = 38.8888500 # falls church / western tip of arlington, rough centroid of stations
DClon = -77.1685800
kounLat = 35.23833
kounLon = -97.46028
kounAlt = 377.0
radarLat=kounLat
radarLon=kounLon
radarAlt=0.0

# # ARMOR
# radarLat = 34.6461
# radarLon = -86.7714
# radarAlt = 180.0

mapProj = MapProjection(projection='eqc', ctrLat=radarLat, ctrLon=radarLon, lat_ts=radarLat, lon_0=radarLon)
geoProj = GeographicSystem()

dx, dy = (8.0e3,)*2
# dx, dy = 500.0, 500.0
minute_intervals = [2.0]
count_scale_factor = dx / 1000.0

b = Bounds()
#Crude limits to the domain of the storm
b.x = (-60e3, 140e3)
b.y = (-150e3, 50e3)
b.z = (-20e3, 20e3)
b.chi2 = (0,1.0)
b.stations = (7,99)
start_time = DateTime(2009,6,10,20,50,0)
end_time   = DateTime(2009,6,10,21,00,0)
max_count_baseline = 450 * count_scale_factor #/ 10.0

pad = 0 #-25e3

# approximate velocity of reference frame. used to adjust the viewport.
# is average speed of LMA density center between 830 and 945 UTC
u = 0 #17.8    # m/s
v = 0 #15.6 
view_dx = b.x[1]-b.x[0] #200.0e3
view_dy = b.y[1]-b.y[0] #200.0e3
# Position at some initial time
x0 = b.x[1]-view_dx/2.0#-150.0e3
y0 = b.y[1]-view_dy/2.0#-150.0e3
t0 = DateTime(2009,6,10,22,40,0)

source_density=False

import glob

if source_density==True:
    LMAfiles = glob.glob("data/LYL*090610_20*.dat.gz") #+ glob.glob("data/LYL*090610_21*.dat.gz")
    lmaManager = LMAdataManager(LMAfiles)
    lma_view = AcuityView(DataSelection(lmaManager.data, b), mapProj, bounds=b)    
else:
    import tables
    LMAfilesHDF = glob.glob('data/LYL*090610_20*.h5')
    LMAtables = []
    for hdffile in LMAfilesHDF:
        h5 = tables.openFile(hdffile)
        table_name = list(h5.root.events._v_children.keys())[0]
        LMAtables.append('/events/'+table_name)
    
        # events = getattr(h5.root.events, table_name)[:]
        # flashes = getattr(h5.root.flashes, table_name)[:]
        # mapping = dict( ( fl, events[events['flash_id'] == fl['flash_id']] ) 
        #                 for fl in flashes if (fl['n_points']>9)
        #               )
        h5.close()
    HDFmanagers = [LMAdataManagerHDF(*args) for args in zip(LMAfilesHDF, LMAtables)]
    # LMAtables   = [ '/events/LMA_080206_080000_3600', '/events/LMA_080206_090000_3600', '/events/LMA_080206_100000_3600' ]

    # h5 = tables.openFile('data/LYLOUT_090610_180000_0600.dat.gz.flash.h5')
    # table_name = h5.root.events._v_children.keys()[0]
    # events = getattr(h5.root.events, table_name)[:]
    # flashes = getattr(h5.root.flashes, table_name)[:]
    # mapping = dict((fl, events[events['flash_id'] == fl['flash_id']]) for fl in flashes)
    # h5.close()



def runtest(lmaManager=None, lma_view=None, HDFmanagers=None):
    # colormap = get_cmap('gist_yarg_r')
    colormap = get_cmap('gist_earth')
    
    density_maxes = []
    total_counts = []
    all_t = []
    
    for delta_minutes in minute_intervals:
        time_delta = DateTimeDelta(0, 0, delta_minutes, 0)
        
        n_frames   = int(ceil((end_time - start_time) / time_delta))
        n_cols = 6
        n_rows = int(ceil( float(n_frames) / n_cols ))
        w, h = figaspect(float(n_rows)/n_cols)

        xedge=np.arange(b.x[0], b.x[1]+dx, dx)
        yedge=np.arange(b.y[0], b.y[1]+dy, dy)
        x_range = b.x[1] - b.x[0]
        y_range = b.y[1] - b.y[0]

        min_count, max_count = 1, max_count_baseline*delta_minutes

        f = figure(figsize=(w,h))
        p = small_multiples_plot(fig=f, rows=n_rows, columns=n_cols)
        p.label_edges(True)
        
        for ax in p.multiples.flat:
            ax.yaxis.set_major_formatter(kilo_formatter)
            ax.xaxis.set_major_formatter(kilo_formatter)

        for i in range(n_frames):
            frame_start = start_time + i*time_delta
            frame_end   = frame_start + time_delta
            b.sec_of_day = (frame_start.abstime, frame_end.abstime)
            b.t = (frame_start, frame_end)
            
            do_plot = False
            flash_extent_density = True
            density = None
            
            if source_density==True:
                lmaManager.refresh(b)
                lma_view.transformed.cache_is_old()
                x,y,t=lma_view.transformed['x','y','t']
                density,edges = np.histogramdd((x,y), bins=(xedge,yedge))
                do_plot=True
            else:
                for lmaManager in HDFmanagers:
                    # yes, loop through every file every time and reselect data.
                    # so wrong, yet so convenient.
                    h5 = lmaManager.h5file
                    if flash_extent_density == False:
                        lmaManager.refresh(b)
                        lma_view = AcuityView(DataSelection(lmaManager.data, b), mapProj, bounds=b)
                        # lma_view.transformed.cache_is_old()
                        x,y,t=lma_view.transformed['x','y','t']
                        if x.shape[0] > 1: do_plot = True
                        break
                    else:
                        # assume here that the bounds sec_of_day day is the same as
                        # the dataset day
                        t0, t1 = b.sec_of_day
                        # events = getattr(h5.root.events, lmaManager.table.name)[:]
                        # flashes = getattr(h5.root.flashes, lmaManager.table.name)[:]
                        
                        event_dtype = getattr(h5.root.events, lmaManager.table.name)[0].dtype
                        events_all = getattr(h5.root.events, lmaManager.table.name)[:]
                        flashes = getattr(h5.root.flashes, lmaManager.table.name)
                        
                        def event_yielder(evs, fls):
                            these_events = []
                            for fl in fls:
                                if (    (fl['n_points']>9) & 
                                        (t0 < fl['start']) & 
                                        (fl['start'] <= t1) 
                                    ):
                                    these_events = evs[evs['flash_id'] == fl['flash_id']]
                                    if len(these_events) != fl['n_points']:
                                        print('not giving all ', fl['n_points'], ' events? ', these_events.shape)
                                    for an_ev in these_events:
                                        yield an_ev

                        
                        # events = np.fromiter((an_ev for an_ev in ( events_all[events_all['flash_id'] == fl['flash_id']] 
                        #                 for fl in flashes if (
                        #                   (fl['n_points']>9) & (t0 < fl['start']) & (fl['start'] <= t1)
                        #                 )
                        #               ) ), dtype=event_dtype)
                        events = np.fromiter(event_yielder(events_all, flashes), dtype=event_dtype)
                        
                        # print events['flash_id'].shape

                        ### Flash extent density ###                        
                        x,y,z = mapProj.fromECEF( 
                                *geoProj.toECEF(events['lon'], events['lat'], events['alt'])
                                )
                                
                        # Convert to integer grid coordinate bins
                        #      0    1    2    3
                        #   |    |    |    |    |
                        # -1.5  0.0  1.5  3.0  4.5
                    
                        if x.shape[0] > 1:
                            density, edges = extent_density(x,y,events['flash_id'].astype('int32'),
                                                            b.x[0], b.y[0], dx, dy, xedge, yedge)
                            do_plot = True                        
                            break
                # print 'density values: ', density.min(), density.max()
                    
            
            if do_plot == True:  # need some data
                # density,edges = np.histogramdd((x,y), bins=(xedge,yedge))
                density_plot  = p.multiples.flat[i].pcolormesh(xedge,yedge,
                                           np.log10(density.transpose()), 
                                           vmin=-0.2,
                                           vmax=np.log10(max_count),
                                           cmap=colormap)
                label_string = frame_start.strftime('%H%M:%S')
                text_label = p.multiples.flat[i].text(b.x[0]-pad+x_range*.01, b.y[0]-pad+y_range*.01, label_string, color=(0.5,)*3, size=6)
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
        
        # moving reference frame correction. all panels will have same limits, based on time of last frame
        view_dt = 0.0 # (frame_start - t0).seconds
        x_ctr = x0 + view_dt*u
        y_ctr = y0 + view_dt*v
        view_x = (x_ctr - view_dx/2.0 - pad, x_ctr + view_dx/2.0 + pad)
        view_y = (y_ctr - view_dy/2.0 - pad, y_ctr + view_dy/2.0 + pad)
        # view_x  = (b.x[0]+view_dt*u, b.x[1]+view_dt*u)
        # view_y  = (b.y[0]+view_dt*v, b.y[1]+view_dt*v)
        
        # print 'making timeseries',
        # time_series = figure(figsize=(16,9))
        # ts_ax = time_series.add_subplot(111)
        # ts_ax.plot_date(mx2num(all_t),total_counts,'-', label='total sources', tz=tz)
        # ts_ax.plot_date(mx2num(all_t),density_maxes,'-', label='max pixel', tz=tz)
        # ts_ax.xaxis.set_major_formatter(time_series_x_fmt)
        # ts_ax.legend()
        # time_filename = 'out/LMA-timeseries_%s_%5.2fkm_%5.1fs.pdf' % (start_time.strftime('%Y%m%d_%H%M%S'), dx/1000.0, time_delta.seconds)
        # time_series.savefig(time_filename)
        # print ' ... done'
        
        print('making multiples', end=' ')
        p.multiples.flat[0].axis(view_x+view_y)
        filename = 'out/LMA-density_%s_%5.2fkm_%5.1fs.pdf' % (start_time.strftime('%Y%m%d_%H%M%S'), dx/1000.0, time_delta.seconds)
        f.savefig(filename, dpi=150)
        print(' ... done')
        f.clf()
        return events
        

        
if __name__ == '__main__':
    do_profile=True
    if do_profile:
        import hotshot
        from hotshot import stats
        prof = hotshot.Profile("multiples_test_profile")
        prof.runcall(runtest, HDFmanagers=HDFmanagers)
        prof.close()
        s=stats.load("multiples_test_profile")
        s.sort_stats("time").print_stats()
    else:
        if source_density:
            res = runtest(lmaManager=lmaManager, lma_view=lma_view)
        else:
            res = runtest(HDFmanagers=HDFmanagers)
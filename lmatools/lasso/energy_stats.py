
# coding: utf-8

# In[1]:

from __future__ import absolute_import
from __future__ import print_function
import os
import numpy as np

# import matplotlib
# fontspec = {'family':'Helvetica', 'weight':'bold', 'size':10}
# matplotlib.rc('font', **fontspec)
# matplotlib.use('Agg')

from stormdrain.pipeline import coroutine
from stormdrain.bounds import Bounds, BoundsFilter
from stormdrain.support.matplotlib.formatters import SecDayFormatter
from matplotlib.ticker import MultipleLocator
from lmatools.flash_stats import raw_moments, raw_moments_for_parameter, central_moments_from_raw, events_flashes_receiver
from lmatools.grid.make_grids import time_edges, seconds_since_start_of_day
import matplotlib.pyplot as plt

def summary_stat_line(start_t, end_t, moments):
    """ Moments should be a tuple of number, mean, variance, skewness, kurtosis """
    line_template = "{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}"
    number, mean, variance, skewness, kurtosis = moments
    energy = number * (mean*mean + variance)
    line = line_template.format(start_t.isoformat(), end_t.isoformat(), 
                                number, mean, variance, skewness, kurtosis, energy, energy/number)
    return line, energy, energy/number


def plot_flash_stat_time_series(basedate, t_edges, stats, major_tick_every=1800):
    t_start, t_end = t_edges[:-1], t_edges[1:]
    fig = plt.figure(figsize=(11,8.5))
    starts = np.fromiter( ((s - basedate).total_seconds() for s in t_start), dtype=float )
    ends = np.fromiter( ((e - basedate).total_seconds() for e in t_end), dtype=float )
    t = (starts+ends) / 2.0
    window_size = (ends-starts)/60.0 # in minutes

    ax_energy = fig.add_subplot(2,2,2)
    ax_energy.plot(t, stats['energy']/window_size, label='Total energy')
    ax_energy.legend()

    ax_count = fig.add_subplot(2,2,1, sharex=ax_energy)
    # ax_count.plot(t, stats['energy_per_flash']/window_size, label='$E_T$ flash$^{-1}$ min$^{-1}$ ')
    ax_count.plot(t, stats['number']/window_size, label='Flash rate (min$^{-1}$)')
    ax_count.legend()#location='upper left')

    ax_mean = fig.add_subplot(2,2,3, sharex=ax_energy)
    mean, std = stats['mean'], stats['variance'], 
    skew, kurt = stats['skewness'], stats['kurtosis']
    sigma = np.sqrt(std)
    ax_mean.plot(t, mean, label = 'Mean flash size (km)')
    ax_mean.plot(t, sigma, label = 'Standard deviation (km)')
    ax_mean.set_ylim(0, 20)
    ax_mean.legend()
#         ax_mean.fill_between(t, mean+sigma, mean-sigma, facecolor='blue', alpha=0.5)

    ax_higher = fig.add_subplot(2,2,4, sharex=ax_energy)
    ax_higher.plot(t, skew, label='Skewness')
    ax_higher.plot(t, kurt, label='Kurtosis')
    ax_higher.set_ylim(-2,10)
    ax_higher.legend()

    for ax in fig.get_axes():
        ax.xaxis.set_major_formatter(SecDayFormatter(basedate, ax.xaxis))  
        ax.set_xlabel('Time (UTC)')
        ax.xaxis.set_major_locator(MultipleLocator(major_tick_every))
        ax.xaxis.set_minor_locator(MultipleLocator(major_tick_every/2.0))
    return fig

class StatResults(object):
    def __init__(self, basedate=None, outdir=None):
        self.basedate = basedate
        self.outdir = outdir
        self.t_start = []
        self.t_end = []
        self.moments = []
        self.energy = []
        self.energy_per = []
        
    def plot(self):
        fig = plt.figure(figsize=(11,8.5))
        starts = np.fromiter( ((s - self.basedate).total_seconds() for s in self.t_start), dtype=float )
        ends = np.fromiter( ((e - self.basedate).total_seconds() for e in self.t_end), dtype=float )
        t = (starts+ends) / 2.0
        window_size = (ends-starts)/60.0 # in minutes
        moments = np.asarray(self.moments) # N by n_moments
        
        ax_energy = fig.add_subplot(2,2,2)
        ax_energy.plot(t, self.energy/window_size, label='Total energy')
        ax_energy.legend()
        
        ax_count = fig.add_subplot(2,2,1, sharex=ax_energy)
        #ax_count.plot(t, self.energy_per/window_size, label='$E_T$ flash$^{-1}$ min$^{-1}$ ')
        ax_count.plot(t, moments[:,0]/window_size, label='Flash rate (min$^{-1}$)')
        ax_count.legend()#location='upper left')
        
        ax_mean = fig.add_subplot(2,2,3, sharex=ax_energy)
        mean, std, skew, kurt = moments[:,1], moments[:,2], moments[:,3], moments[:,4]
        sigma = np.sqrt(std)
        ax_mean.plot(t, mean, label = 'Mean flash size (km)')
        ax_mean.plot(t, sigma, label = 'Standard deviation (km)')
        ax_mean.set_ylim(0, 20)
        ax_mean.legend()
#         ax_mean.fill_between(t, mean+sigma, mean-sigma, facecolor='blue', alpha=0.5)
        
        ax_higher = fig.add_subplot(2,2,4, sharex=ax_energy)
        ax_higher.plot(t, skew, label='Skewness')
        ax_higher.plot(t, kurt, label='Kurtosis')
        ax_higher.set_ylim(-2,10)
        ax_higher.legend()
        
        for ax in fig.get_axes():
            ax.xaxis.set_major_formatter(SecDayFormatter(self.basedate, ax.xaxis))  
            ax.set_xlabel('Time (UTC)')
            ax.xaxis.set_major_locator(MultipleLocator(3600))
            ax.xaxis.set_minor_locator(MultipleLocator(1800))
        return fig
        
    @coroutine
    def stats_for_frame(self, t_start,t_end, target=None):
        while True:
            raw_moments = (yield)
            central, standard = central_moments_from_raw(raw_moments)        
            summary, energy, energy_per = summary_stat_line(t_start, t_end, standard)
            self.t_start.append(t_start)
            self.t_end.append(t_end)
            self.moments.append(standard)
            self.energy.append(energy)
            self.energy_per.append(energy_per)
            if target is not None:
                target.send(summary)
        
    @coroutine
    def stats_printer(self):
        lines = []
        try:
            while True:
                line = (yield)
                lines.append(line)
        except GeneratorExit:
            lines.sort()
            for l in lines:
                print(l)
            fig = self.plot()
            fig.savefig(self.outdir+'/moment-energy-timeseries.pdf')
        

def flash_size_stats_for_intervals(t_start, t_end, dt, lat_bounds=None, lon_bounds=None, outdir=None):
    """ Create a pipeline to process an arbitrary number of flashes which are 
        accumulated in dt-length intervals between t_start and t_end.
        
        t_start, t_end: datetime.datetime instances
        dt: datetime.timedelta
        
        Returns:
        printer_of_stats: 
        all_frame_targets: a list of targets that receive (events, flashes) from the flash file reader

        all_frame_targets can be passed as other_analysis_targets to lmatools.flash_stats.plot_spectra_for_files.
        
        When all data have been read (i.e., plot_spectra_for_files has returned), the accumulated statistics are
        retrieved by closing the targets returned by this function, i.e.
            for target in timeframe_targets:
                target.close()
            printer_of_stats.close()
        
        plot_spectra_for_files itself could be refactored to make use of these same time windows natively, but
        that's a project for later. The dependency here on stormdrain is another issue; it may make sense to have
        lmatools depend on stormdrain.
    """
    t_edges, duration = time_edges(t_start, t_end, dt.total_seconds())
    t_ref, t_edges_seconds = seconds_since_start_of_day(t_start, t_edges)
    n_frames = len(t_edges)-1
    
    results = StatResults(basedate=t_ref, outdir=outdir)
    printer_of_stats = results.stats_printer()
    
    all_frame_targets = []
    for n in range(n_frames):
        # New BoundsFilter object for filtering flashes and flash_stats.stats_for_parameter that sits behind the bounds filter
        t0, t1 = t_edges_seconds[n:n+2]
        b = Bounds()
        
        if lat_bounds is not None:
            b.ctr_lat = lat_bounds
        if lon_bounds is not None:
            b.ctr_lon = lon_bounds

        b.start = t0, t1 
        statframer = results.stats_for_frame(t_edges[n], t_edges[n+1], target=printer_of_stats)
        momentizer = raw_moments_for_parameter('area', preprocess=np.sqrt, output_target=statframer)
        bound_filt = BoundsFilter(bounds=b, target=momentizer)
        ev_fl_rcvr = events_flashes_receiver(target=bound_filt.filter())
        all_frame_targets.append(ev_fl_rcvr)
    return printer_of_stats, all_frame_targets


def flash_size_stats(flashes):
    """ Using the area data in the flashes table, calculate the flash size moments
        return moments, energy, and energy_per_flash, where energy is the dimensional
        energy defined by Bruning and Thomas (2015, JGR).
    
        Returned is a single row array with column names as follows:
        result_dtype = [
            ('number','f4'), 
            ('mean', 'f4')
            ('variance', 'f4')
            ('skewness', 'f4')
            ('kurtosis', 'f4')
            ('energy', 'f4')
            ('energy_per_flash', 'f4')
           ]
    """
    area = flashes['area']
    length_scale = np.sqrt(area) # length scale per Bruning and Thomas (2015, JGR)
    
    result_dtype = [
        ('number','f4'), 
        ('mean', 'f4'),
        ('variance', 'f4'),
        ('skewness', 'f4'),
        ('kurtosis', 'f4'),
        ('energy', 'f4'),
        ('energy_per_flash', 'f4'),
        ('specific_energy', 'f4'),
        ('total_energy', 'f4'),
       ]
    result = np.zeros((1,), dtype=result_dtype)
     
    # Moment calculation reduces all flash lengths into a few representative moments.
    central, standard = central_moments_from_raw(raw_moments(length_scale, n_moments=5))
    number, mean, variance, skewness, kurtosis = standard
    result['number'] = number
    result['mean'] = mean
    result['variance'] = variance
    result['skewness'] = skewness
    result['kurtosis'] = kurtosis
    
    energy = number * (mean*mean + variance) # Bruning and Thomas (2015, JGR)
    result['energy'] = energy
    result['energy_per_flash'] = energy/number
    
    #Moments for energy:
    try:
        result['specific_energy'] = flashes['specific_energy'].sum()
        result['total_energy']  = flashes['total_energy'].sum()
    except ValueError:
        "Print skipping energy quantities"
    return result

####ADDED FUNCTION:   
def plot_energy_stats(size_stats, basedate, t_edges, outdir):
    t_start, t_end = t_edges[:-1], t_edges[1:]
    starts = np.fromiter( ((s - basedate).total_seconds() for s in t_start), dtype=float )
    ends = np.fromiter( ((e - basedate).total_seconds() for e in t_end), dtype=float )
    t = (starts+ends) / 2.0
    
    specific_energy = size_stats
    
    figure = plt.figure(figsize=(15,10))
    ax     = figure.add_subplot(111)
    ax.plot(t,specific_energy,'k-',label='Specific Energy',alpha=0.6)
    plt.legend()
    # ax.set_xlabel('Time UTC')
    ax.set_ylabel('Specific Energy (J/kg)')
    
    for axs in figure.get_axes():
        axs.xaxis.set_major_formatter(SecDayFormatter(basedate, axs.xaxis))  
        axs.set_xlabel('Time (UTC)')
        axs.xaxis.set_major_locator(MultipleLocator(1800))
        axs.xaxis.set_minor_locator(MultipleLocator(1800/2))
    
    return figure

def plot_tot_energy_stats(size_stats, basedate, t_edges, outdir):
    t_start, t_end = t_edges[:-1], t_edges[1:]
    starts = np.fromiter( ((s - basedate).total_seconds() for s in t_start), dtype=float )
    ends = np.fromiter( ((e - basedate).total_seconds() for e in t_end), dtype=float )
    t = (starts+ends) / 2.0
    
    specific_energy = np.abs(size_stats)
    
    figure = plt.figure(figsize=(15,10))
    ax     = figure.add_subplot(111)
    ax.plot(t,specific_energy,'k-',label='Total Energy',alpha=0.6)
    plt.legend()
    # ax.set_xlabel('Time UTC')
    ax.set_ylabel('Total Energy (J)')
    
    for axs in figure.get_axes():
        axs.xaxis.set_major_formatter(SecDayFormatter(basedate, axs.xaxis))  
        axs.set_xlabel('Time (UTC)')
        axs.xaxis.set_major_locator(MultipleLocator(1800))
        axs.xaxis.set_minor_locator(MultipleLocator(1800/2))
    
    return figure
    
# In[6]:

from stormdrain.support.matplotlib.poly_lasso import LassoFilter

class TimeSeriesPolygonLassoFilter(object):
    """ Given a series of N polygons which are valid between N+1 time edges,
        set up a routing structure that filters data with the correct
        polygon given the time of the data.
        
        time_edges is a sequence of datetime objects, which
        are converted to float seconds since basedate.
        
        time_key is the name used to index to retrieve time in
        float seconds from the named array to be filtered. Seconds should be
        relative to the same basedate as used above.
        
        polys is a sequence of N tuples, where each tuple in another length-M tuple 
        of (x,y) polygon vertices. Each of the N polygons may have a different number
        M of vertices.
        
        polys = (
          ( (xa1,ya1), (xa2, ya2), (xa3, ya3) ),
          ( (xb1,yb1), (xb2, yb2), (xb3, yb3), (xb4, yb4) ),
        )
        
        """
        
    
    def __init__(self, *args, **kwargs):
        self.coord_names = kwargs.pop('coord_names', [])
        self.time_edges = kwargs.pop('time_edges', [])
        self.time_key = kwargs.pop('time_key',[])
        self.basedate = kwargs.pop('basedate', [])
        polys = kwargs.pop('polys', [])
        
        self.poly_lookup = {}
        
        for iv, verts in enumerate(polys):
            t0 = (self.time_edges[iv] - self.basedate).total_seconds()
            t1 = (self.time_edges[iv+1] - self.basedate).total_seconds()
            lassofilter = LassoFilter(coord_names=self.coord_names)
            lassofilter.verts = np.asarray(verts)
            self.poly_lookup[(t0, t1)] = lassofilter
                                
    def filter_mask(self, a):
        """ Filter out all points in a that belong to any of the polys.
            
            Returns an array mask that indexes a to filter the data.
        """
        
        # assume nothing makes it through
        good = np.zeros(a.shape, dtype=bool)
        
        t = a[self.time_key]
        for (t0, t1), poly_target in self.poly_lookup.items():
            # make sure this is only one-side inclusive to eliminate double-counting
            in_time = (t >= t0) & (t < t1)
            reduced = a[in_time]
            if reduced.size > 0:
                inpoly = poly_target.filter_mask(reduced)
                # accept anything in this polygon at this time
                good[in_time] |= inpoly
        return good

    




"""
read_flashes 
filterer.filt_ev_fl
brancher(events_flashes_receiever, other_analysis_targets=timeframe_targets)

ev_fl_rx ->
histogrammer
energy calculation
plotter

timeframe_targets is
ev_fl_rx ->
bound_filt (on time)


read_flashes -> filterer.filt_ev_fl -> length_for_these_flashes -> brancher -> timeframe_targets

"""
from __future__ import absolute_import
from __future__ import print_function
import os

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
from numpy.lib.recfunctions import stack_arrays #, append_fields

from lmatools import coordinateSystems
from lmatools.flash_stats import length_from_area, volumetric_length_from_points, vertical_length_distribution, gen_flash_events
from lmatools.grid.make_grids import time_edges, seconds_since_start_of_day

from stormdrain.pipeline import coroutine, Branchpoint
from stormdrain.support.matplotlib.formatters import SecDayFormatter

def gen_fractal_length_for_flashes(events, flashes, D, b_s, alt_bins):
    """ Given events and flashes, calculate the fractal length
        as described in Bruning and MacGorman (2015). 
        D is the fractal dimension, and b_s is the minimum box size.

        This function returns a dictionary with the 2D and 3D length as well as the
        convex hull data 
    """
    GeoSys = coordinateSystems.GeographicSystem()
    for fl_srcs, fl in gen_flash_events(events, flashes):
        if fl_srcs.shape[0] < 5:
            print("Skipping flash because source count less than 5 prevents calculation of volume")
            continue
        x,y,z = GeoSys.toECEF(fl_srcs['lon'], fl_srcs['lat'], fl_srcs['alt'])
        x,y,z = x/1000.0, y/1000.0, z/1000.0
        # These calls will fail if the chi2 and stations criteria are too stringent and reduce the number of points below the minimum needed to construct the geometry (5, according to QHULL's error)
        simplex_centroids, simplex_volumes,  volume, L_fractal, length_weighted = volumetric_length_from_points(x,y,z,D,b_s)
        simplex_lon,simplex_lat,simplex_alt = GeoSys.fromECEF(simplex_centroids[:,0]*1000.0, simplex_centroids[:,1]*1000.0, simplex_centroids[:,2]*1000.0)
        alt_bins, bin_total_src, bin_total_length = vertical_length_distribution(fl_srcs['alt']/1000.0, simplex_alt/1000.0, length_weighted, alt_bins, norm=True)

        L_capacitor = length_from_area(fl['area'],D,b_s)

        results = {'2D':{'length':L_capacitor},
                   '3D':{'length':L_fractal, 'volume':volume,
                         'subvolumes':{'centroids':simplex_centroids, 'volumes':simplex_volumes, 'lengths':length_weighted},
                         'vertical_profile':{'alt_bins':alt_bins, 'src_density':bin_total_src, 'length_density':bin_total_length}
                        },
                   'flash':fl}
        yield results
        
class FractalLengthProfileCalculator(object):
    def __init__(self, D, b_s, alt_bins):
        self.alt_bins = alt_bins
        self.D = D
        self.b_s = b_s
        
    def process_events_flashes(self, events, flashes):
        """ Given a table of flashes and their corresponding events table,
            calculate the Bruning and Thomas (2015, JGR) fractal length
            for each flash and produce a vertical profile of flash lengths
            for all flashes.
        
            Returns:
            per_flash_data: List of dictionaries containing the flash length 
                and triangulation geometry data, including the flash volume.
            IC_totals: array of summary data for each profile; see dtype below.
            CG_totals: as above but for CGs only.
            These data are not normalized by time.
        
            totals_dtype = [
                ('fractal_length_2D_hull', 'f8'),
                ('fractal_length_3D_hull', 'f8'),
                ('source_density_profile', 'f8', alt_bin_shape), 
                ('length_density_profile', 'f8', alt_bin_shape),
                ]   
            
        """
        results = [result for result in gen_fractal_length_for_flashes(
                       events, flashes, self.D, self.b_s, self.alt_bins)
                  ]
                       
        totals_CG = self._sum_profiles(results, is_CG_flag=True)
        totals_IC = self._sum_profiles(results, is_CG_flag=False)

        # IC_totals = self._process_totals(totals_IC, duration)
        # CG_totals = self._process_totals(totals_CG, duration)
        
        return results, totals_IC, totals_CG

    def normalize_profiles(self, totals_seq, durations):
        """ Given a sequence of results from repeated calls to self.process_events_flashes,
            return profiles normalized by durations (an array of durations for each interval
            corresponding to the items in total_seq). also return a source density 
            
            returns:    
            lengths_2D, lengths_3D, scaled_sources, L_profile_rate, max_L_bin
            total 2D and 3D lengths, source and 3D length profiles, and the maximum 
            length per unit time per km altitude across all bins.
        """
        totals = np.asarray(list(totals_seq))
        L_profile_rate = np.squeeze(totals['length_density_profile'].T/durations)
        max_L_bin = L_profile_rate.max()
        src_profile_rate = totals['source_density_profile'].T/durations
        scaled_sources = np.squeeze(max_L_bin * src_profile_rate / src_profile_rate.max())
    
        lengths_2D = np.squeeze(totals['fractal_length_2D_hull'])/durations
        lengths_3D = np.squeeze(totals['fractal_length_3D_hull'])/durations
        return lengths_2D, lengths_3D, scaled_sources, L_profile_rate, max_L_bin

                           
    def _sum_profiles(self, results, is_CG_flag=False):
        
        alt_bins = self.alt_bins
        # one fewer density value (bin interval) than bin edges.
        alt_bin_shape = (alt_bins.shape[0]-1,)
        res_dtype = [ # this the dtype of the return value, i.e., the sums
            ('fractal_length_2D_hull', 'f8'),
            ('fractal_length_3D_hull', 'f8'),
            ('source_density_profile', 'f8', alt_bin_shape), 
            ('length_density_profile', 'f8', alt_bin_shape), #, (alt_bins.shape[0]-1,) ,)
            ]   
        
        if len(results) == 0:
            return np.zeros((1,), dtype=res_dtype)
        
        # Check to see if there are CG flags in the flash data table
        # just look at the first flash and assume the rest have it if one does.
        if 'CG' in results[0]['flash'].dtype.names:
            results_iter = (ri for ri in results if
                                    (ri['flash']['CG'] == is_CG_flag) )
        else:
            results_iter = results
            
        # def get_res_iter():
        res_iter = ( (
             r['2D']['length'], 
             r['3D']['length'],
             tuple(r['3D']['vertical_profile']['src_density']),
             tuple(r['3D']['vertical_profile']['length_density'])
            ) for r in results_iter
            )

        # http://stackoverflow.com/questions/19201868/how-to-set-dtype-for-nested-numpy-ndarray
        each_result = np.fromiter(res_iter, dtype=res_dtype)
        total = np.empty((1,), dtype=res_dtype)
        for colname in total.dtype.names:
            total[colname] = each_result[colname].sum(axis=0)
        return total
        
    def make_time_series_plot(self, basedate, t_edges, 
        lengths_2D, lengths_3D, scaled_sources, L_profile_rate, max_L_bin,
        label_t_every=3600., figsize=(7.5,10)):      
        """ Arguments: 
            t_edges: N+1 bin boundaries corresponding to the N time series values in the
                    other arguments. Units: seconds since start of day.
            lengths_2D, lengths_3D, scaled_sources, L_profile_rate, max_L_bin: 
                    The values returned by normalize_profiles
            label_t_every: interval in seconds at which to label the time axis
            figsize: (width, height) of figure in inches (passed to matplotlib.figure)
            
        """
        import matplotlib.pyplot as plt
        cmap = 'cubehelix_r'        
    
        fig = plt.figure(figsize=figsize)
        ax_L = fig.add_subplot(311)
        ax_prof_L = fig.add_subplot(312)
        ax_prof_src = fig.add_subplot(313)
        
        min_L_bin = 0*max_L_bin
        
        starts, ends = t_edges[:-1], t_edges[1:]
        t_centers = (starts+ends)/2.0
    
        src_pm = ax_prof_src.pcolormesh(t_edges, self.alt_bins, scaled_sources, 
            cmap=cmap, vmin=min_L_bin, vmax=max_L_bin)
        src_pm.set_rasterized(True)
        cbar_src = fig.colorbar(src_pm, ax=ax_prof_src, orientation='horizontal')
        cbar_src.set_label('Source count per height interval per time\n(scaled to max length, km/km/min)')
    
        L_pm = ax_prof_L.pcolormesh(t_edges, self.alt_bins, L_profile_rate, 
            cmap=cmap, vmin=min_L_bin, vmax=max_L_bin)
        L_pm.set_rasterized(True)
        cbar_L = fig.colorbar(L_pm, ax=ax_prof_L, orientation='horizontal')
        cbar_L.set_label('Length per height interval per time\n(km/km/min)')

        for ax in (ax_prof_src, ax_prof_L):
            ax.set_ylabel('Altitude (km)')
                    
        ax_L.plot(t_centers, lengths_2D, 
            label='2D Hull Area')
        ax_L.plot(t_centers, lengths_3D,
            label='3D Hull Volume')
        ax_L.legend()
        ax_L.set_ylabel('Fractal Length (km/min)')
    
        for ax in (ax_L, ax_prof_L, ax_prof_src):
            ax.xaxis.set_major_formatter(SecDayFormatter(basedate, ax.xaxis))  
            ax.set_xlabel('Time (UTC)')
            ax.xaxis.set_major_locator(MultipleLocator(label_t_every))
        return fig
        
    def write_profile_data(self, basedate, t_edge, outfilebase,
            lengths_2D, lengths_3D, scaled_sources, L_profile_rate, max_L_bin,
            partition_kind='total'):
        """ Arguments:
            outfilebase: file base name, including directory.
            basedate: date against which t_edge is referenced.
            t_edge: N+1 bin boundaries corresponding to the N time series values in the
                    other arguments. Units: seconds since start of day given by basedate.
            lengths_2D, lengths_3D, scaled_sources, L_profile_rate, max_L_bin: 
                    The values returned by normalize_profiles
            partition_kind: one of 'total', 'IC', 'CG', or another unique tag
                    classifying this kind of profile. Appended to filename just before.

            The filename is calculated as [outfile_base]_[parition_kind].txt

        """
        starts = t_edge[:-1]
        ends = t_edge[1:]
        header = ""
        header += "# LMA channel length distribution\n"
        header += "# Base date = " + basedate.isoformat() +"\n"
        header += "# Altitude_bins = " + str(self.alt_bins.tolist()) +"\n"
        header += "# starts, ends, lengths_2D, lengths_3D, [scaled_sources x Nbins], [L_profile_rate x Nbins]\n"
        text_dump = open(outfilebase+'_{0}.txt'.format(partition_kind), 'w')
        text_dump.write(header)
        # shape of scaled_sources and L_profile_rate are (N_alt_bins, N_times), and loop is over the first dimension, so take transpose so loop is over N_times
        for s0, e1, l2, l3, Sprof, Lprof in zip(starts, ends, lengths_2D, lengths_3D, scaled_sources.T, L_profile_rate.T):
            text_dump.write('{0}, {1}, {2}, {3}, {4}, {5}\n'.format(
                                                    s0, e1, l2, l3,
                                                    str(Sprof.tolist()), str(Lprof.tolist())
                                                    ))

        


@coroutine
def length_for_these_flashes(D, b_s, alt_bins, chi2=5.0, stations=5, target=None):
    """ Receive events, flashes arrays. Calculate flash length using fractal dimension D and step length b_s
        For each flash, will send out a dictionary with:
            capacitor_length: total length determined from flash area
            volumetric_length: total length determined from flash volume
    """
    GeoSys = coordinateSystems.GeographicSystem()

    while True:
        pts, flashes = (yield) # pts == events
        areas = flashes['area']
    
        capacitor_length = length_from_area(areas,D,b_s)

        for L_capacitor, fl in zip(capacitor_length, flashes):
            fl_id=fl['flash_id']
            this_flash = (pts['flash_id']==fl_id)
            good = (pts['stations'] >= stations) & (pts['chi2']<=chi2)
            fl_srcs = pts[this_flash & good] 
            
            if fl_srcs.shape[0] < 5:
                print(("Skipping flash because original source count reduced from {0}={4} to {1} for flash {2} at {3}".format(
                    pts[this_flash].shape[0], fl_srcs.shape[0], fl_id, fl['start'], fl['n_points']
                )))
                continue

            x,y,z = GeoSys.toECEF(fl_srcs['lon'], fl_srcs['lat'], fl_srcs['alt'])
            x,y,z = x/1000.0, y/1000.0, z/1000.0
            
            # These calls will fail if the chi2 and stations criteria are too stringent and reduce the number of points below the minimum needed to construct the geometry (5, according to QHULL's error)
            simplex_centroids, simplex_volumes,  volume, L_fractal, length_weighted = volumetric_length_from_points(x,y,z,D,b_s)
            simplex_lon,simplex_lat,simplex_alt = GeoSys.fromECEF(simplex_centroids[:,0]*1000.0, simplex_centroids[:,1]*1000.0, simplex_centroids[:,2]*1000.0)
            alt_bins, bin_total_src, bin_total_length = vertical_length_distribution(fl_srcs['alt']/1000.0, simplex_alt/1000.0, length_weighted, alt_bins, norm=True)

            results = {'2D':{'length':L_capacitor},
                       '3D':{'length':L_fractal, 
                             'subvolumes':{'centroids':simplex_centroids, 'volumes':simplex_volumes, 'lengths':length_weighted},
                             'vertical_profile':{'alt_bins':alt_bins, 'src_density':bin_total_src, 'length_density':bin_total_length}
                            },
                       'flash':fl
                       }
                   

            if target is not None:
                 target.send(results)

@coroutine
def in_time_range(t_min, t_max, target):
    results_in_range = []
    try:
        while True:
            results = (yield)
            start = results['flash']['start']
            if (start >= t_min) & (start < t_max): 
                results_in_range.append(results)
    except GeneratorExit:
        target.send((t_min, t_max, results_in_range))

            
class StatResults(object):
    def __init__(self, alt_bins, basedate=None):
        self.basedate = basedate
        # t_start, t_end, flashes_in_range should be equal length
        # and mutually indexable because of how stats_rcvr builds them
        self.t_start = []
        self.t_end = []
        self.results_in_range = []
        
        # it's imperative this be passed in so that _sum_profiles can 
        # return an 0-filled array of the right shape if no flashes are in
        # any time interal.
        self.alt_bins = alt_bins
        
    def _sum_profiles(self, results, is_CG_flag=False):
        
        alt_bins = self.alt_bins
        # one fewer density value (bin interval) than bin edges.
        alt_bin_shape = (alt_bins.shape[0]-1,)
        res_dtype = [ # this the dtype of the return value, i.e., the sums
            ('fractal_length_2D_hull', 'f8'),
            ('fractal_length_3D_hull', 'f8'),
            ('source_density_profile', 'f8', alt_bin_shape), 
            ('length_density_profile', 'f8', alt_bin_shape), #, (alt_bins.shape[0]-1,) ,)
            ]   
        
        if len(results) == 0:
            return np.zeros((1,), dtype=res_dtype)
            
        # def get_res_iter():
        res_iter = ( (
             r['2D']['length'], 
             r['3D']['length'],
             tuple(r['3D']['vertical_profile']['src_density']),
             tuple(r['3D']['vertical_profile']['length_density'])
            ) for r in results if (r['flash']['CG'] == is_CG_flag)
            )
            # return res_iter
        
        # print res_dtype
        # for ri in get_res_iter():
        #     print ri
        # http://stackoverflow.com/questions/19201868/how-to-set-dtype-for-nested-numpy-ndarray
        each_result = np.fromiter(res_iter, dtype=res_dtype)
        # print each_result.shape, each_result.dtype
        total = np.empty((1,), dtype=res_dtype)
        for colname in total.dtype.names:
            print((colname, each_result[colname].shape))
            total[colname] = each_result[colname].sum(axis=0)
        return total

    def total_all_profiles(self):
        """ List of profiles, summed over all flashes for all times. """
        total_CG = np.asarray([self._sum_profiles(r, is_CG_flag=True) for r in self.results_in_range])
        total_IC = np.asarray([self._sum_profiles(r, is_CG_flag=False) for r in self.results_in_range])
        # total_all  = total_IC + total_CG
        
        return total_IC, total_CG

    def process_totals(self, totals, durations):
        L_profile_rate = np.squeeze(totals['length_density_profile'].T/durations)
        max_L_bin = L_profile_rate.max()
        src_profile_rate = totals['source_density_profile'].T/durations
        scaled_sources = np.squeeze(max_L_bin * src_profile_rate / src_profile_rate.max())
    
        lengths_2D = np.squeeze(totals['fractal_length_2D_hull'])/durations
        lengths_3D = np.squeeze(totals['fractal_length_3D_hull'])/durations
        return lengths_2D, lengths_3D, scaled_sources, L_profile_rate, max_L_bin

            
    def plot_stats(self, outfile):
        outfile_base, outfile_ext = os.path.splitext(outfile)
        
        n_frames = len(self.results_in_range)
        starts = np.asarray(self.t_start)
        ends = np.asarray(self.t_end)
        
        t_edges = np.asarray(self.t_start + self.t_end[-1:])
        t_centers = (starts+ends)/2.0
        alt_bins = self.alt_bins
        durations = (t_edges[1:] - t_edges[:-1])/60.0 #minutes
                
        totals_IC, totals_CG = self.total_all_profiles()
        # fig_total = make_plot(totals_total)
        
        IC_totals = process_totals(totals_IC, durations)
        CG_totals = process_totals(totals_CG, durations)
        
        def make_plot(lengths_2D, lengths_3D, scaled_sources, L_profile_rate, max_L_bin):      
            import matplotlib.pyplot as plt
            cmap = 'cubehelix_r'        
        
            fig = plt.figure(figsize=(7.5,10))
            ax_L = fig.add_subplot(311)
            ax_prof_L = fig.add_subplot(312)
            ax_prof_src = fig.add_subplot(313)
            
            min_L_bin = 0*max_L_bin
        
            src_pm = ax_prof_src.pcolormesh(t_edges, alt_bins, scaled_sources, 
                cmap=cmap, vmin=min_L_bin, vmax=max_L_bin)
            src_pm.set_rasterized(True)
            cbar_src = fig.colorbar(src_pm, ax=ax_prof_src, orientation='horizontal')
            cbar_src.set_label('Source count per height interval per time\n(scaled to max length, km/km/min)')
        
            L_pm = ax_prof_L.pcolormesh(t_edges, alt_bins, L_profile_rate, 
                cmap=cmap, vmin=min_L_bin, vmax=max_L_bin)
            L_pm.set_rasterized(True)
            cbar_L = fig.colorbar(L_pm, ax=ax_prof_L, orientation='horizontal')
            cbar_L.set_label('Length per height interval per time\n(km/km/min)')

            for ax in (ax_prof_src, ax_prof_L):
                ax.set_ylabel('Altitude (km)')
                        
            ax_L.plot(t_centers, lengths_2D, 
                label='2D Hull Area')
            ax_L.plot(t_centers, lengths_3D,
                label='3D Hull Volume')
            ax_L.legend()
            ax_L.set_ylabel('Fractal Length (km/min)')
        
            for ax in (ax_L, ax_prof_L, ax_prof_src):
                ax.xaxis.set_major_formatter(SecDayFormatter(self.basedate, ax.xaxis))  
                ax.set_xlabel('Time (UTC)')
                ax.xaxis.set_major_locator(MultipleLocator(3600))
            return fig
        
        
        fig_IC = make_plot(*IC_totals)
        fig_CG = make_plot(*CG_totals)

        # fig_total.savefig(outfile_base+'_total'+outfile_ext)
        # fig_total.clf()
        fig_IC.savefig(outfile_base+'_IC'+outfile_ext)
        fig_IC.clf()
        fig_CG.savefig(outfile_base+'_CG'+outfile_ext)
        fig_CG.clf()
        
        for diagnose in CG_totals:
            print((diagnose.dtype, diagnose.shape))
        
        for partition_kind, kind_totals in zip (('IC', 'CG'), (IC_totals, CG_totals)):
            lengths_2D, lengths_3D, scaled_sources, L_profile_rate, max_L_bin = kind_totals
            header = ""
            header += "# LMA channel length distribution\n"
            header += "# Base date = " + self.basedate.isoformat() +"\n"
            header += "# Altitude_bins = " + str(alt_bins.tolist()) +"\n"
            header += "# starts, ends, lengths_2D, lengths_3D, [scaled_sources x Nbins], [L_profile_rate x Nbins]\n"
            text_dump = open(outfile_base+'_{0}.txt'.format(partition_kind), 'w')
            text_dump.write(header)
            # shape of scaled_sources and L_profile_rate are (N_alt_bins, N_times), and loop is over the first dimension, so take transpose so loop is over N_times
            for s0, e1, l2, l3, Sprof, Lprof in zip(starts, ends, lengths_2D, lengths_3D, scaled_sources.T, L_profile_rate.T):
                text_dump.write('{0}, {1}, {2}, {3}, {4}, {5}\n'.format(
                                                        s0, e1, l2, l3,
                                                        str(Sprof.tolist()), str(Lprof.tolist())
                                                        ))
            text_dump.close()


        
        
        
    @coroutine
    def timeframe_results_rcvr(self):
        t_min, t_max, results_in_range = (yield)
        self.t_start.append(t_min)
        self.t_end.append(t_max)
        self.results_in_range.append(results_in_range)
        
        
                 
def length_stats_for_intervals(t_start, t_end, dt, D, b_s, chi2=5.0, stations=5):
    """ 
    """
    t_edges, duration = time_edges(t_start, t_end, dt.total_seconds())
    t_ref, t_edges_seconds = seconds_since_start_of_day(t_start, t_edges)
    n_frames = len(t_edges)-1
    
    max_alt, d_alt = 20.0, 0.5
    alt_bins = np.arange(0.0,max_alt+d_alt, d_alt)
    
    results_aggregator = StatResults(alt_bins, basedate=t_ref)
    
    all_frame_targets = []
    for n in range(n_frames):
        t0, t1 = t_edges_seconds[n:n+2]
        statframer = results_aggregator.timeframe_results_rcvr()
        this_frame = in_time_range(t0, t1, statframer)
        all_frame_targets.append(this_frame)
        
    brancher = Branchpoint(all_frame_targets)
    ev_fl_rcvr = length_for_these_flashes(D, b_s, alt_bins,
                        chi2=chi2, stations=stations, target=brancher.broadcast())
    
    return ev_fl_rcvr, all_frame_targets, results_aggregator

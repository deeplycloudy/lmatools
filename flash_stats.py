import glob
import os

import numpy as np

from lmatools.lma_io import read_flashes
from lmatools.density_to_files import coroutine, Branchpoint

from lmatools.flashsort.autosort.flash_stats import hull_volume

def length_from_area(A,D,b_s):
    return ( (np.sqrt(A))**D ) / (b_s**(D-1))

def volumetric_length_from_points(x,y,z,D, b_s):
    xyz = np.vstack((x,y,z)).T
    
    volume, vertices, simplex_volumes = hull_volume(xyz)
    simplex_centroids = np.average(vertices[:,:], axis=1)
    # The simplex volumes have negative values since they are oriented (think surface normal direction for a triangle)

    vol_len_scale = volume**(1.0/3.0)
    S_i = np.abs(simplex_volumes)**(1.0/3.0) #S_i
    P_i = (S_i**D)  / (b_s**(D-1.0)) #P_i
    L_3 = (vol_len_scale**D)  / (b_s**(D-1.0)) #L_3
    
    length_weighted = (S_i / P_i) * L_3
    sum_weights =  (S_i/P_i).sum()
    print "The sum of the ratio S_i/P_i is not equal to one, but is {0}".format(sum_weights)
    
    # Therefore, divide by the sum of the weights
    length_weighted =  (S_i / P_i) * L_3 / sum_weights
    
    return simplex_centroids, np.abs(simplex_volumes), volume, L_3, length_weighted

def vertical_length_distribution(src_alt, simplex_alt, simplex_lengths, 
        alt_bins, norm=True):
    """ given input altitudes and lengths in km, create vertical
        profiles of source counts and total length.
        
        Returns alt_bins, bin_total_src, bin_total_length
        
        If norm==True, divide the counts by the bin width, returning
        km, counts/km and km/km. Otherwise just return km, counts and km.
        """
        
    # Not sure why we're not using histogram here, so that's a TODO
    # d_alt = 0.5
    d_alt = alt_bins[1:]-alt_bins[:-1]
    # alt_bins = np.arange(0.0,max_alt+d_alt, d_alt)
    bin_total_length = np.zeros(alt_bins.shape[0]-1, dtype=float)
    bin_total_src = np.zeros(alt_bins.shape[0]-1, dtype=float)
    # bin_total_length_sq = np.zeros(alt_bins.shape[0]-1, dtype=float)
    tri_bin_idx = np.digitize(simplex_alt, alt_bins)
    src_bin_idx = np.digitize(src_alt,alt_bins)
    src_bin_idx[src_bin_idx>(bin_total_src.shape[0]-1)]=bin_total_src.shape[0]-1

    for idx in src_bin_idx:
        bin_total_src[idx] += 1

    for lw,idx in zip(simplex_lengths,tri_bin_idx):
        bin_total_length[idx]+=lw
        # bin_total_length_sq[idx] += lw*lw
    # bin_total_length[tri_bin_idx] += length_weighted
    if norm==True:
        return alt_bins, bin_total_src/d_alt, bin_total_length/d_alt
    else:
        return alt_bins, bin_total_src, bin_total_length


def bin_center(bins):
    return (bins[:-1] + bins[1:]) / 2.0
    

def energy_plot_setup(fig=None, subplot=111, bin_unit='km'):
    """ Create an energy spectrum plot with a 5/3 slope line and an spare line
        to be used for plotting the spectrum. The spare line is intially located
        on top of the 5/3 line 
        If fig is specified, the spectrum axes will be created on that figure
        in the subplot position given by subplot.
        
        Returns 
    """    
    if fig is None:
        import matplotlib.pyplot as plt
        fig = plt.figure()
    
    spectrum_ax = fig.add_subplot(subplot)
    spectrum_ax.set_xlabel('Flash width ($\sqrt{A_h}$, %s)' % (bin_unit,))
    spectrum_ax.set_ylabel('$E(l) \mathrm{(m^2 s^{-2} km^{-1})}$')
    spectrum_ax.set_xlim(10**-1, 10**2)
    spectrum_ax.set_ylim(10**0, 10**8)
    spectrum_ax.set_yscale('log')
    spectrum_ax.set_xscale('log')
    
    #1e-2 to 1e4
    min_pwr = -2
    max_pwr = 4
    delta_pwr = 0.1
    powers = np.arange(min_pwr, max_pwr+delta_pwr, delta_pwr)
    flash_1d_extent = 10**powers
    wavenumber = (2*np.pi)/flash_1d_extent
    inertialsubrange = 10**6 * (wavenumber)**(-5.0/3.0)
    
    spectrum_line_artist = spectrum_ax.loglog(flash_1d_extent, inertialsubrange, 'r')[0]
    fivethirds_line_artist = spectrum_ax.loglog(flash_1d_extent, inertialsubrange, 'k')[0]
    
    return fig, spectrum_ax, fivethirds_line_artist, spectrum_line_artist


def calculate_energy_from_area_histogram(histo, bin_edges, duration, scaling=1.0):
    """ Given a histogram and bin edges for flash area, calculate the specific energy density
        in units of (m^2/s^2) / km, as in Bruning and MacGorman 2013, J. Atmos. Sci. 
        duration is the total number of seconds over which the flashes were counted.
        
        bin_edges are assumed to be in km^2. histo is a count corresponding to the intervals
        specified by bin_edges.
        
        Before return, the spectrum is multipled by scaling (default = 1.0)
    """
    duration=float(duration)
    flash_1d_extent = bin_center(np.sqrt(bin_edges))    
    bin_widths = np.sqrt(bin_edges[1:] - bin_edges[:-1])    
    # This should give   s^-2                 m^2                      km^-1   =  m s^-2 km^-1
    specific_energy = (histo/duration * flash_1d_extent*1000.0)**2.0 / (bin_widths) # flash_1d_extent #bin_widths
    specific_energy *= scaling
    return flash_1d_extent, specific_energy

def plot_energy_from_area_histogram(histo, bin_edges, bin_unit='km', save=False, fig=None, color_cycle_length=1, color_map='gist_earth', duration=600.0):
    """ Histogram for flash width vs. count """
    
    fig, spectrum_ax, fivethirds_line_artist, spectrum_artist = energy_plot_setup()
    spectrum_ax.set_title(save.split('/')[-1].split('.')[0])
    
    flash_1d_extent, specific_energy = calculate_energy_from_area_histogram(histo, bin_edges, duration)
    spectrum_artist.set_data(flash_1d_extent, specific_energy)
    
    if save==False:
        plt.show()
    else:
        # ax.set_title(save)
        fig.savefig(save)
        fig.clf()
    


@coroutine
def histogram_for_parameter(parameter, bin_edges, target=None):
    """ General coroutine that accepts a named numpy array with field parameter
        and calcualtes a histogram using bin_edges. Target is sent histogram, edges.
    """
    while True:
        a = (yield)
        histo, edges = np.histogram(a[parameter], bins=bin_edges)
        if target is not None:
            target.send((histo, edges))

@coroutine
def events_flashes_receiver(target=None):
    """ Passes along only flashes """
    while True:
        events, flashes = (yield)
        if target is not None:
            target.send(flashes)


@coroutine
def histogram_accumulate_plot(plotter, histo_array=None, save=False, fig=None):
    bin_edges=None
    try:
        while True:        
            histo, edges = (yield) 
            
            if bin_edges is None:
                bin_edges = edges
            else:
                assert (bin_edges == edges).all()
            
            if histo_array is None:
                histo_array  = histo
            else:
                histo_array += histo
    except GeneratorExit:
        plotter(histo_array, bin_edges, save=save, fig=fig)

@coroutine        
def raw_moments_for_parameter(parameter, preprocess=None, n_moments=5, output_target=None):
    """ This coroutine builds up raw moments by streaming samples into the coroutine.
        When the samples are exhausted by a GeneratorExit, An array of size n_moments 
        will be sent to the output_target, and will contain the zeroth and n_moments-1 
        higher moments. The higher moments are already divided by the zeroth moment.
        
        Receives a data array, and calculate basic statistics from the distribution
        of the named parameter in that array. Optionally call the preprocess function 
        on the parameter before calculating the moment.
        
    """
    sample_raw_moments = np.zeros(n_moments, dtype='f8')
    try:
        while True:
            data = (yield)
            a = data[parameter]
            if preprocess is not None:
                a = preprocess(a)
            # calculate the sample moments
            for i in range(n_moments):
                sample_raw_moments[i] += (a**i).sum()
            
    except GeneratorExit:
        sample_raw_moments[1:] /= sample_raw_moments[0]
        if output_target is not None:
            output_target.send(sample_raw_moments)
        

def central_moments_from_raw(raw):
    """ Returns ctr, std where
            ctr = the zero through fourth central moments
            std = (number, mean, variance, skewness, kurtosis)
            
            the zeroth and first central moments are set to zero.
    """
    ctr = np.zeros_like(raw)
    # ctr[0] = 0
    # ctr[1] = raw[1] - raw[1]
    ctr[2] = raw[2] - raw[1]*raw[1]
    ctr[3] = 2*(raw[1]**3) - 3*raw[1]*raw[2] + raw[3]
    ctr[4] = -3*(raw[1]**4) + 6*raw[1]*raw[1]*raw[2] - 4*raw[3]*raw[1] + raw[4]
    # number, mean, variance, skewness, kurtosis
    std = raw[0], raw[1], ctr[2], ctr[3]/(ctr[2]**1.5), (ctr[4]/(ctr[2]*ctr[2]) - 3)
    return ctr, std
    
def footprint_stats(h5_filenames, save=False, fig=None, min_points=10, 
                    base_date=None, other_analysis_targets=None, filterer=None):
    """ filter should be a non-running coroutine that receives the (events, flashes)
        arrays emitted by io.read_flashes and sends filtered (events, flashes) arrays to 
        a target defined by this function and passed to filter by keyword argument.
        
    """

    #1e-2 to 1e4
    min_pwr = -2
    max_pwr = 4
    delta_pwr = 0.1
    powers = np.arange(min_pwr, max_pwr+delta_pwr, delta_pwr)
    footprint_bin_edges = 10**powers

    plotter = plot_energy_from_area_histogram
    
    histogram_plot = histogram_accumulate_plot(plotter, save=save, fig=fig)
    histogrammer=histogram_for_parameter('area', footprint_bin_edges, target=histogram_plot)
    ev_fl_rx = events_flashes_receiver(target=histogrammer)
    brancher = Branchpoint([ev_fl_rx])
    if other_analysis_targets is not None:
        for t in other_analysis_targets:
            brancher.targets.add(t)
    if filterer is not None:
        rcvr = filterer(target=brancher.broadcast())
    else:
        rcvr = brancher.broadcast()
    read_flashes(h5_filenames, rcvr, min_points=min_points, base_date=base_date)

            

def plot_spectra_for_files(h5_filenames, min_points, time_criteria, distance_criteria,
                           outdir_template='./figures/thresh-{0}_dist-{1}_pts-{2}/',
                           other_analysis_targets=None, base_date=None, filterer=None):
    """ Make standard plots of the flash energy spectrum. There will be one spectrum created
        for each file in h5_filenames.

        min_pts, time_criteria, distance_critera are tuples of point and time-space thresholds, all the same length.
        They will be looped over, used to generate outdir_template which 
        needs {0},{1},{2}, which will be filled in with time, distance, and min_pts criteria

        The path in outdir_template will be created if it does not exist.
    """
    import matplotlib.pyplot as plt
    fig = plt.figure()
    for dpt in min_points:
        for dt in time_criteria:
            for dx in distance_criteria:
                outdir = outdir_template.format(dt,dx,dpt)
                if not os.path.exists(outdir):
                    os.makedirs(outdir)

                for h5_file in h5_filenames:
                    file_basename = os.path.split(h5_file)[-1].split('.')[:-3][0]
                    figure_file = os.path.join(outdir, '{0}-energy.pdf'.format(file_basename))
                    footprint_stats([h5_file], save=figure_file, fig=fig, min_points=dpt, 
                                    base_date=base_date, other_analysis_targets=other_analysis_targets,
                                    filterer=filterer)


    
if __name__ == '__main__':
    # '/data/20090610/data'
    # min_points = 10

    # --- All times (6+ hours), .15 s and 3 km ---
    # h5_filenames = glob.glob('29may-thresh-0.15_dist-3000.0/LYL*.flash.h5')
    # h5_filenames += glob.glob('30may-thresh-0.15_dist-3000.0/LYL*.flash.h5')
    # h5_filenames = glob.glob('30may-thresh-0.15_dist-3000.0/LYL*0130*.flash.h5')

    # h5_filenames = glob.glob('/Users/ebruning/code/McCaul Flash/test20040529/fixed-area-run/thresh-0.15_dist-3000.0/LYL*.flash.h5')
        
    # footprint_stats(h5_filenames, min_points=min_points)

    # --- 0130 - 0140 UTC, range of different space/time criteria ---
    if True:
        import matplotlib.pyplot as plt
        
        fig = plt.figure()
        
        min_points = (10,)
        time_critera = (0.15,)
        distance_critera = (3000.0,)

        filename_template = '/Users/ebruning/code/McCaul Flash/test20040529/fixed-area-run-expandedtime/thresh-{0}_dist-{1}/LYL*.flash.h5'
        for dpt in min_points:
            for dt in time_critera:
                for dx in distance_critera:
                    pattern = filename_template.format(dt,dx)
                    print pattern
                    h5_filenames = glob.glob(pattern)
                    for h5_file in h5_filenames:
                        file_basename = os.path.split(h5_file)[-1].split('.')[:-3][0]
                        figure_file = '/Users/ebruning/code/McCaul Flash/test20040529/fixed-area-run-expandedtime/thresh-{0}_dist-{1}/histos/{2}-footprint_histogram-{3}pts.pdf'.format(dt,dx,file_basename,dpt)
                        # print figure_file
                        footprint_stats([h5_file], save=figure_file, fig=fig, min_points=dpt)
                # break
            # break
    
    # To open all figures in Preview, you can use a pattern like so from sh/bash
    # open thresh-0.{05,1,15,2,25}_dist*/footprint_histogram-10pts.pdf
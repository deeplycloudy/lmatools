import glob
import os

import numpy as np

from lmatools.io import read_flashes
from lmatools.density_to_files import coroutine


def bin_center(bins):
    return (bins[:-1] + bins[1:]) / 2.0
    


def plot_histogram(histo, bin_edges, bin_unit='km', save=False, fig=None, color_cycle_length=1, color_map='gist_earth'):
    """ Histogram for flash width vs. count """
    import matplotlib.pyplot as plt
    
    flash_1d_extent = bin_center(np.sqrt(bin_edges))
    
    if fig is None:
        fig = plt.figure()
        
    # ax.loglog(flash_1d_extent, histo)
    ax  = fig.add_subplot(111)
    
    bin_widths = np.sqrt(bin_edges[1:] - bin_edges[:-1])
    
    # This should give   s^-2                 m^2                      km^-1   =  m s^-2 km^-1
    specific_energy = (histo/600.0 * flash_1d_extent*1000.0)**2.0 / (bin_widths) # flash_1d_extent #bin_widths
    # print '{0} J total energy'.format((specific_energy*bin_widths).sum())
    ax.loglog(flash_1d_extent, specific_energy, 'r')
        
    # wavelenth ^-5/3
    wavenumber = (2*np.pi)/flash_1d_extent
    inertialsubrange = 10**6 * (wavenumber)**(-5.0/3.0)
    plt.plot(flash_1d_extent, inertialsubrange)
    plt.text(1,1,save.split('/')[-1].split('.')[0], color='g')
    
    plt.xlabel('Flash width ($\sqrt{A_h}$, %s)' % (bin_unit,))
    plt.ylabel('$E(l) \mathrm{(m^2 s^{-2} km^{-1})}$')
    ax.set_ylim(10**0, 10**8)
    
    if save==False:
        plt.show()
    else:
        # ax.set_title(save)
        fig.savefig(save)
        fig.clf()
    


@coroutine
def histogram_for_flash_parameter(paramname, bin_edges, plotter, flash_id_key='flash_id', histo_array=None, save=False, fig=None):
    try:
        while True:
            events, flashes = (yield)
        
            histo, edges = np.histogram(flashes[paramname], bins=bin_edges)
            
            if histo_array is None:
                histo_array  = histo
            else:
                histo_array += histo
    except GeneratorExit:
        plotter(histo_array, bin_edges, save=save, fig=fig)
        

    
def footprint_stats(h5_filenames, save=False, fig=None, min_points=10):
    
    # start_time = datetime(2009,6,10, 20,0,0)
    # end_time   = datetime(2009,6,10, 21,0,0)

    #1e-2 to 1e4
    min_pwr = -2
    max_pwr = 4
    delta_pwr = 0.1
    powers = np.arange(min_pwr, max_pwr+delta_pwr, delta_pwr)
    footprint_bin_edges = 10**powers

    plotter = plot_histogram

    histogrammer = histogram_for_flash_parameter('area', footprint_bin_edges, plotter, save=save, fig=fig)

    read_flashes(h5_filenames, histogrammer, min_points=min_points)
    
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
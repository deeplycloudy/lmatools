from __future__ import absolute_import
from __future__ import print_function
import os
import glob
import tables as T

base_out_dir = '/Users/ebruning/out/flash_sort/'
params_filename = 'input_params.py'

analysis_params = ('thresh_critical_time', 'distance')
counted = 'total_flashes'
data_lists = [[] for i in analysis_params]
plot_data = dict(list(zip(analysis_params, data_lists)))
plot_data[counted] = []
    
for path, dirnames, filenames in os.walk(base_out_dir):
    
    h5_search = os.path.join(path, '*.h5')
    h5files = glob.glob(h5_search)
    
    print("Checking directory ", path, end=' ')
    
    try:
        params = open(os.path.join(path, params_filename), 'r')
        param_data = params.read()
        params.close()
        params = eval(param_data)
        print("... found param data")
        # print type(params), params
    except IOError:
        # Since there is no run-parameters data,
        # this isn't a directory with results of a flash code run
        print("... nothing")
        continue
    
    total_flashes = 0
    for h5filename in h5files:
        h5 = T.openFile(h5filename)
        
        for flash_table in h5.walkNodes('/flashes', 'Leaf'):
            flash_pts = [fl['n_points'] for fl in flash_table if fl['n_points'] > 9]
            n_flashes = len(flash_pts)
            total_flashes +=n_flashes
            # n_flashes = flash_table.shape[0]
            print('%s -- %d flashes' % (path, n_flashes,))
        h5.close()
    
    for param in analysis_params:
        key = param
        value = params[param]
        plot_data[key].append(value)

    plot_data[counted].append(total_flashes)

import matplotlib.pyplot as plt    
plt.scatter(plot_data[analysis_params[1]], plot_data[analysis_params[0]], 
                c=plot_data[counted], s=100 )#, vmin=50, vmax=100)
plt.xlabel(analysis_params[1])
plt.ylabel(analysis_params[0])
plt.title(counted)
plt.colorbar()
plt.show()
        
things_to_consider = """
Decreasing singletons are a good indicator of whether or not flash algorithm is properly grouping at long range - not leaving too many orphans

looking at > 10 points and all flashes makes a very big difference in terms of singletons that get orphaned away from flashes for sensitve changes in distance/time critera

Increasing flash counts at large space/time criteria have to do with many singletons getting
grouped together to flashes with more than 10 (or whatever minimum number) of points. Singleton 
filtering is therefore important prior to flash sorting.

6 stations on 26 May has too much env noise.
"""

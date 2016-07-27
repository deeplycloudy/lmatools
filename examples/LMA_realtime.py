from __future__ import absolute_import
from __future__ import print_function
import datetime
import os
import glob

from subprocess import call
gzip_command = "gzip"

REALTIME_DIR = "/data/realtime-test/"
LMA_SUBDIR = "LMA/"
AWIPS_SUBDIR = "AWIPS/"

LMA_DATA_DIR = os.path.join(REALTIME_DIR, LMA_SUBDIR)
AWIPS_GRID_DIR = os.path.join(REALTIME_DIR, AWIPS_SUBDIR)

# ----- Grid and plot parameters -----
now = datetime.datetime.now()
start_time = datetime.datetime(now.year, now.month, now.day, now.hour, now.minute, 0, 0)
end_time   = start_time + datetime.timedelta(0, 60)
frame_interval=60.0

target_datum = 'NAD83'
target_ellipse = 'GRS80'
target_proj = 'eqc'
ctr_lat = 33.7
ctr_lon = -101.8
dx=1.0e3
dy=1.0e3
x0, y0 = 0.0, 0.0
x_bnd = (x0-200e3, x0+200e3)
y_bnd = (y0-200e3, y0+200e3)


now = datetime.datetime.now()

# ----- Process raw LMA data every minute into a 1-min-long file ----- 
from lmatools.io.fakeLMA import fake_LMA_file, late2011_header
LMA_ASCII_outfile = fake_LMA_file(year=now.year, month=now.month, day=now.day, hour=now.hour, minute=now.minute, second=0, 
                duration=60, header_template=late2011_header, outpath=LMA_DATA_DIR)
print("Wrote ", LMA_ASCII_outfile)

# ----- Sort flashes, creating HDF5 LMA data file ----- 
from lmatools.flashsort.autosort.autorun import run_files_with_params
thresh_critical_time = 0.15
distance = 3000.0

params = {'stations':(5,13),
          'chi2':(0,2.0),  
          'lat':(-90., 90.),
          'lon':(-180., 180.),
          'alt':(0.,20000.),     
          'distance':distance, 'thresh_duration':3.0, 
          'thresh_critical_time':thresh_critical_time, 'merge_critical_time':0.5, 
          'ascii_flashes_out':'flashes_out.dat',
          'ctr_lat':ctr_lat, 'ctr_lon':ctr_lon,
          }

params['thresh_critical_time'] = thresh_critical_time
params['distance'] = distance

# tag = 'thresh-%s_dist-%s' % (thresh_critical_time, distance)
# outdir = os.path.join(base_out_dir, tag)
# os.mkdir(outdir)

# info = open(os.path.join(outdir, 'input_params.py'), 'w')
# info.write(str(params))
# info.close()

h5_filenames = run_files_with_params([LMA_ASCII_outfile], LMA_DATA_DIR, params, retain_ascii_output=False)
# h5_path_searcher = os.path.join(LMA_DATA_DIR,'LYL*.flash.h5')
# h5_filenames = glob.glob(h5_path_searcher)

# ----- Create NetCDF grids ----- 
from lmatools.grid.make_grids import grid_h5flashfiles
from lmatools.grid.AWIPS_tools import write_AWIPS_netcdf_grid

x_coord, y_coord, lons, lats, extent_density_grid, AWIPS_outfiles, field_names = grid_h5flashfiles(
            h5_filenames, start_time, end_time, frame_interval=frame_interval, 
            dx=dx, dy=dy, x_bnd=x_bnd, y_bnd=y_bnd, ctr_lon=ctr_lon, ctr_lat=ctr_lat,
            proj_name = target_proj, proj_datum = target_datum, proj_ellipse = target_ellipse,
            output_writer = write_AWIPS_netcdf_grid, output_kwargs = {'refresh_minutes':int(frame_interval/60.0)},
            outpath = AWIPS_GRID_DIR,
            )
            


# ----- Put ASCII data and AWIPS grids on the LDM queue ----
call([gzip_command, LMA_ASCII_outfile])
call([gzip_command] + list(AWIPS_outfiles))
LDM_files = [f + '.gz' for f in [LMA_ASCII_outfile]+list(AWIPS_outfiles)]
for f in LDM_files:
    ldm_command = "sudo -u ldm /home/ldm/bin/pqinsert -v -f EXP -q /home/ldm/var/queues/ldm.pq {0:s}".format(f)
    print(ldm_command)
    # call(ldm_command.split())


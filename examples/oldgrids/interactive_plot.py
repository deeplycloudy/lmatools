""" This script makes a multi-panel plot from NetCDF data, just like a simple example,
    except that the figure, panels, starttime, and calculated filename are returned for
    further processing.
    
    This is useful if one wishes to add additional data to the plot, e.g., an GPS track.
    Be careful to use the correct map projection information from the NetCDF file when 
    plotting new data.
    
    It also demonstrates how to change the backend from the default (non-windowed) backend so 
    that the figure can be worked with interactively.
     
    It is necessary to edit the "backend" import line to use a backend that can
    run on your computer.

"""

from __future__ import absolute_import
from lmatools.vis.multiples_nc import make_plot
f,p,start,fname=make_plot('LMA_20040622_052000_600_10src_flash_extent.nc', 'flash_extent', n_cols=1, do_save=False)

# Plot some other stuff here, using 
# ax = p.multiples.flat[frame_sequence_id]
# ax.plot(...) 
# etc.

# ---EDIT ME --- Import a backend that can run on your computer
from matplotlib.backends.backend_macosx import FigureManagerMac, FigureCanvasMac

# Replace the existing canvas with one that matches the backend
canvas = FigureCanvasMac(f)
f.set_canvas(canvas)

# Tell the OS to make a window
fm = FigureManagerMac(canvas, 1)
f.show()
# Close figure
# fm.destroy() 

# You can still save the figure after any additions.
# fig.savefig(filename, dpi=150)

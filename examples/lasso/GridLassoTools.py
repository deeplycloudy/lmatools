import netCDF4
import numpy as np
import pandas as pd
import scipy as sci

# from scipy.spatial import *
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.widgets import Widget
from matplotlib.colors import LogNorm, Normalize
from lmatools.coordinateSystems import GeographicSystem, MapProjection
geosys = GeographicSystem()

from ipywidgets import widgets
# from IPython.display import HTML
# from IPython.display import Javascript
# from IPython.display import display

from lmatools.vis import ctables as ct

import logging, json

# from http://stackoverflow.com/questions/3488934/simplejson-and-numpy-array/24375113#24375113
class NumpyAwareJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray): # and obj.ndim == 1:
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

class PolyLasso(Widget):
    """
    A lasso widget that allows the user to define a lasso region by
    clicking to form a polygon.
    """

    def __init__(self, figure, callback=None,
                 line_to_next=True,
                 useblit=True,
                 color='blue'):
        """
        Create a new lasso.

        *ax* is the axis on which to draw

        *callback* is a function that will be called with arguments (*ax*, *line*, *verts*).
        *verts* is a list of (x,y) pairs that define the polygon, with the first and
        last entries equal.

        *line_to_next* controls whether or not a line is drawn to the current
        cursor position as the polygon is created

        *useblit* = *True* is the only thorougly tested option.

        """
        self.axes = None
        self.figure = figure
        self.canvas = self.figure.canvas
        self.useblit = useblit
        self.line_to_next = line_to_next
        self.background = None
        self.color = color
        # if useblit:
        #     self.background = self.canvas.copy_from_bbox(self.axes.bbox)


        self.verts = []
        self.line = None
        self.callback = callback
        self.cids = []
        self.cids.append(self.canvas.mpl_connect('button_release_event', self.onrelease))
        self.cids.append(self.canvas.mpl_connect('motion_notify_event', self.onmove))

        # moved to after axes are chosen
        # self.cids.append(self.canvas.mpl_connect('draw_event', self.ondraw))

    def ondraw(self, event):
        """ draw_event callback, to take care of some blit artifacts """
        self.background = self.canvas.copy_from_bbox(self.axes.bbox)
        if self.line:
            self.axes.draw_artist(self.line)
        self.canvas.blit(self.axes.bbox)

    def do_callback(self, event):
        """ idle_event callback after polygon is finalized. """
        # Clear out callbacks.
        for cid in self.cids:
            self.canvas.mpl_disconnect(cid)
        self.callback(self.axes, self.line, self.verts)
        self.cleanup()

    def cleanup(self):
        """ Remove the lasso line. """
        # Clear out callbacks
        for cid in self.cids:
            self.canvas.mpl_disconnect(cid)
        self.axes.lines.remove(self.line)
        self.canvas.draw()

    def finalize(self):
        """ Called when user makes the final right click """
        # all done with callbacks pertaining to adding verts to the polygon
        for cid in self.cids:
            self.canvas.mpl_disconnect(cid)
        self.cids = []

            # Matplotlib will not draw the closed polygon until we step completely
            # out of this routine. It is desirable to see the closed polygon on screen
            # in the event that *callback* takes a long time. So, delay callback until
            # we get an idle event, which will signal that the closed polygon has also
            # been drawn.
        self.cids.append(self.canvas.mpl_connect('idle_event', self.do_callback))

    def draw_update(self):
        """ Adjust the polygon line, and blit it to the screen """
        self.line.set_data(zip(*self.verts))
        if self.useblit:
            self.canvas.restore_region(self.background)
            if self.axes is not None:
                self.axes.draw_artist(self.line)
                self.canvas.blit(self.axes.bbox)
        else:
            self.canvas.draw_idle()

    def onmove(self, event):
        """ Update the next vertex position """
        if self.line == None: return
        self.verts[-1] = ((event.xdata, event.ydata))
        if self.line_to_next:
            self.draw_update()

    def onrelease(self, event):
        """ User clicked the mouse. Add a vertex or finalize depending on mouse button. """
        if self.verts is None: return

        # detect the axes to draw in automatically on first click
        if (self.axes == None) & (event.inaxes is not None):
            self.axes = event.inaxes
            if self.useblit:
                self.background = self.canvas.copy_from_bbox(self.axes.bbox)
                self.cids.append(self.canvas.mpl_connect('draw_event', self.ondraw))
        else:
            if event.inaxes != self.axes: return

        if event.button == 3:
            # Right click - close the polygon
            # Set the dummy point to the first point
            self.verts[-1] = self.verts[0]
            self.draw_update()
            # Greater than three verts, since a triangle
            # has a duplicate point to close the poly.
            if len(self.verts)>3:
                self.do_callback(event)
                # self.finalize()
            else:
                print('Need at least three vertices to make a polygon')
                self.cleanup()
            return

        # The rest pertains to left click only
        if event.button != 1: return

        if (self.axes is not None):
            if (self.line == None):
                # Deal with the first click
                self.line=Line2D([event.xdata], [event.ydata],
                            linestyle='-', marker='H', color=self.color, lw=1, ms=4, animated=True)
                self.axes.add_line(self.line)
                self.verts.append((event.xdata, event.ydata))

            # finalize vertex at this click, set up a new one that changes as mouse moves
            self.verts[-1] = (event.xdata, event.ydata)
            self.verts.append((event.xdata, event.ydata))

            self.draw_update()
            
class NCgridLasso(object):
    def __init__(self, NCDFcollection, log_path):
        self.NCs = NCDFcollection
        
        self.do_log=False
        self.log_path = log_path

        self.logger = logging.getLogger('.'.join((__name__,'Poly')))
        hdlr = logging.FileHandler(self.log_path) #<--Change to User's Path.
        # JSON-compliant so that each log entry can be read by the json parser
        # json_fmt = '{"created":"%(asctime)s", "level":"%(levelname)s", "poly":%(message)s}'

        # replace the stupid comma in the asctime default format with a period, since it's milliseconds
        # http://stackoverflow.com/questions/6290739/python-logging-use-milliseconds-in-time-format/7517430#7517430
        json_fmt = '{"created":"%(asctime)s.%(msecs)03d000", "level":"%(levelname)s", "poly":%(message)s}'

        formatter = logging.Formatter(json_fmt, datefmt='%Y-%m-%dT%H:%M:%S')
        hdlr.setFormatter(formatter)
        self.logger.addHandler(hdlr)
        self.logger.setLevel(logging.INFO)
        
        self.cmap = ct.NWSRef
        self.vmin = -30#vmin
        self.vmax = 30#vmax
        
        # Defaults. These change depending on input from the notebook widget.
        self.field = self.NCs.grid_name
        self.frame_id = 0
        
        # Make a copy of the data that reflects the current field and slice index.
        self._sync_data()
        
        
        self.f = plt.figure(figsize=(10,7))
        self.ax = self.f.add_subplot(111)
        self.cax = self.f.add_axes([0.93, .1, 0.02, .8])
        
        
        # cc = self.ax.pcolormesh(self.x, self.y, self.field_in,
        #                vmin=self.vmin, vmax=self.vmax, cmap=self.cmap)
        # self.ax.set_ylim(0,20e3)
        # self.ax.set_xlim(0,25e3)
        # self.title = "{1} {2} | Frame {3}".format(self.field.capitalize(),
        #                                                         self.t.isoformat()+" UTC", self.frame_id)
        # self.f.colorbar(cc, cax=self.cax)
        # self.ax.set_title(self.title)
        
        self.ax.set_xlabel('East distance (km)')
        self.ax.set_ylabel('North distance (km)')
        
        
    def _sync_data(self):
        
        t = self.NCs.times[self.frame_id]
        xedge, yedge, data = self.NCs.data_for_time(t)

        # Check for a 3D grid, and assume vertical coordinate is the first dimension
        if len(data.shape) == 3:
            data_3d = data
            if self.field == 'flash_extent':
                data = data_3d.max(axis=0) # composite
            elif self.field == 'flash_footprint':
                # have to exclude zeros before taking min.
                nonzero = np.ma.masked_less(data_3d, 1.0, copy=False)
                data = nonzero.min(axis=0) # composite, but for minimum
            elif self.field == 'flash_initiation':
                data = data_3d.sum(axis=0) # total flash rate

        """ Use chosen field and sweep to update arrays of data"""
        self.t = t
        self.x = xedge
        self.y = yedge
        self.field_in = data
        
        # self.time_text = t.isoformat() + 'Z ' # this should be reconfigured to change with the sweep.
        
    def do_lasso(self, log=False):
        """ Attach to B4D_panel_lasso_drawn exchange to get the panels, axes,
            MPL line artist, and verts for each lasso drawn
        """
        self.do_log = log
        
        lock = self.f.canvas.widgetlock
        if lock.locked()==False:
            self.active_lasso = PolyLasso(self.f, self.callback)
            lock(self.active_lasso)
        elif lock.locked()==True:
            print("Please deselect other tools to use the lasso.")
    
    def callback(self, ax, lasso_line, verts):
        from matplotlib.path import Path
        self.v = verts

        self.lasso_line = lasso_line
        print("releasing lock")
        self.f.canvas.widgetlock.release(self.active_lasso)

#         self.active_lasso=None
        
        x_verts = np.ravel(self.v)[0::2]
        y_verts = np.ravel(self.v)[1::2]

        fname, frame_i = self.NCs._time_lookup[self.t] 
        geosys, mapproj = self.NCs.get_projection()
        nc = netCDF4.Dataset(fname)
        
        # if 'Lambert_Azimuthal_Equal_Area' in nc.variables.keys():
        #     nc_proj = nc.variables['Lambert_Azimuthal_Equal_Area']
        #     print nc_proj
        #     proj_name = 'laea'
        #     ctrlon, ctrlat, ctralt = (nc_proj.longitude_of_projection_origin,
        #                               nc_proj.latitude_of_projection_origin,
        #                               nc_proj.altitude_of_projection_origin,)
        #     mapproj = MapProjection(proj_name, ctrLat = ctrlat, ctrLon=ctrlon,
        #                              lat_0=ctrlat, lon_0=ctrlon)
        lon, lat, alt = geosys.fromECEF(*mapproj.toECEF(x_verts, y_verts, np.zeros_like(x_verts)))
        # else:
        #     print "projection not found"
        #     lon = None
        #     lat = None
        #     alt = None
            
        if self.do_log == True:
            # verts_dict = { 'Verts':self.v,
            #                       'Title':'Lasso Verts'
            #                      }
            lass_verts = {'title':'Lasso Verts','field':self.field,
                          'x_verts':x_verts,
                          'lon_verts':lon,
                          'lat_verts':lat,
                          'y_verts':y_verts,
                          'frame_id':self.frame_id,
                          'frame_time':self.t.isoformat(),
                          'filename':fname}
            self.logger.info(json.dumps(lass_verts, cls=NumpyAwareJSONEncoder))
            # self.logger.info(lass_verts)

        # return self.y, self.z, self.v, self.vel
    
    def data_table(self):
        '''If desired, this function takes the PolyLasso vertices drawn, and returns
        a Pandas Dataframe containing data for each radar field constrained within
        the Polygon.'''

        # mask_shape, xvert, yvert = self.mask()

        mask_shape, xvert, yvert = self.masked_data()

        poly_velocities = (self.vel[mask_shape])
        poly_ref = (self.ref[mask_shape])
        poly_sw = (self.sw[mask_shape])


        titles = ['Velocities', 'Reflectivity', 'Spectrum Width']

        table = pd.DataFrame([poly_velocities, poly_ref, poly_sw]).T
        table.columns = titles

        return table
    

def gridLassoOverplotDefault(lasso_widget, ax):
    pass

class NCgridLassoWidgets(object):
    def __init__(self, ncgridlasso, overplot_func=gridLassoOverplotDefault):
        """ The user may pass a function to overplot_func in order to 
            add additional features to the plot. That function is passed 
            the lasso widget instance which called it, as well as
            the matplotlib axes on which the grid's data was just plotted.
            The default function does nothing.
        
        """
    
    
        self.ncgridlasso = ncgridlasso
        self.overplot_func = overplot_func
        
        self.value_range_for_field = {
            'flash_extent':(1, 1000.0),
            'flash_initiation':(1,100.0),
            'flash_footprint':(1, 100000.0),
            'default':(1,1000.0),
            }
        self.value_colors_for_field = {
            'flash_extent': ct.NWSRef,
            'flash_initiation': ct.NWSRef,
            'flash_footprint': 'cubehelix',
            }
        
        field_type = widgets.Dropdown(description='Field:', options=[self.ncgridlasso.field])
        field_type.observe(self.change_field, 'value')
        
        initial_vmin, initial_vmax = self.value_range_for_field[ncgridlasso.field]
        
        toggle = widgets.Checkbox(description='Log Lasso?', value=False)
        toggle_log_scale = widgets.Checkbox(description='Logarithmic', value=True)
        toggle_log_scale.observe(self.redraw, 'value')

        draw = widgets.Button(description='Draw Lasso')
        draw.on_click(self.lasso_draw)
        
        N_frames = len(self.ncgridlasso.NCs.times)
        frame_selection = widgets.IntSlider(description="Frame", min=0, max=N_frames-1)
        frame_selection.value = 0
        frame_selection.observe(self.change_frame, 'value')

        v_min = widgets.FloatText(description='Min')
        v_min.value = initial_vmin
        v_min.observe(self.v_min_set, 'value')

        v_max = widgets.FloatText(description='Max')
        v_max.value= initial_vmax
        v_max.observe(self.v_max_set, 'value')

        # view_mask = widgets.Button(description='View Masked Region', value=False)
        # view_mask.on_click(self.view_mask)
        #
        revert_plot = widgets.Button(description='Redraw', value=False)
        revert_plot.on_click(self.redraw)
        #
        # filtered = widgets.Button(description='View Filtered Data', value=False)
        # filtered.on_click(self.view_filter)
        
        contain = widgets.VBox(description='NetCDF grid tools')
        contain.children = [draw,toggle,v_min,v_max,frame_selection,
                            # filtered, view_mask,
                             revert_plot, toggle_log_scale]
        
        self.draw = draw
        self.toggle = toggle
        self.field_type = field_type
        self.frame_selection = frame_selection
        self.v_min = v_min
        self.v_max = v_max
        self.revert_plot = revert_plot
        self.toggle_log_scale = toggle_log_scale
        # self.view_mask = view_mask
        # self.filtered = filtered
        self.container = contain
        
        self.plot_update()
        
    def lasso_draw(self, onclick):
        self.lasso = self.ncgridlasso.do_lasso(log=self.toggle.value)
        
    def redraw(self, value):
        self.plot_update()
        
    def plot_update(self):
        '''Draws the updated plot for the user selected values of:
        self.radar_type = value in dropdown widget for all radar scan types.
        self.v_max = max value for the selected radar type from the dictionary object value_range_for_field above.
        self.v_min = min value for the selected radar type from the dictionary object value_range_for_field above.
        self.cmap = color map with respect to the radar_type selected.
        '''

        if self.toggle_log_scale.value == True:
            norm = LogNorm()
        else:
            norm = Normalize()    

        self.initial_cmap = self.value_colors_for_field[self.field_type.value]

        self.ncgridlasso.ax.clear()
        self.ncgridlasso.cax.clear()
        cc = self.ncgridlasso.ax.pcolormesh(self.ncgridlasso.x, self.ncgridlasso.y, self.ncgridlasso.field_in,
                        vmin=np.float(self.v_min.value), vmax=np.float(self.v_max.value),
                        cmap=self.initial_cmap, norm=norm)
        self.title = "Composite {1} at {0} UTC".format(self.ncgridlasso.t.isoformat(), self.field_type.value)

        self.ncgridlasso.ax.set_title(self.title)
        self.ncgridlasso.f.colorbar(cc, self.ncgridlasso.cax)
        self.overplot_func(self, self.ncgridlasso.ax)
        self.ncgridlasso.f.canvas.draw()

    def change_field(self, change):
        the_field = change['new']
        self.ncgridlasso.field = the_field
        if the_field in self.value_range_for_field:
            new_min, new_max = self.value_range_for_field[the_field]
        else:
            new_min, new_max = self.value_range_for_field['default']
            self.value_range_for_field[the_field] = (new_min, new_max)
        self.v_min.value = new_min
        self.v_max.value = new_max
        
        self.ncgridlasso._sync_data()

        self.plot_update()
        
    def change_frame(self, change):
        value = change['new']
        self.ncgridlasso.frame_id = value
        self.ncgridlasso._sync_data()
        self.plot_update()
        
    def v_min_set(self, change):
        vmin = change['new']
        # This lookup will always succeed because the min/max for
        # any field is set from the default (if necessary) in change_radar
        # upon any change to the radar field
        old_min, old_max = self.value_range_for_field[self.field_type.value]
        self.ncgridlasso.vmin = vmin
        self.plot_update()

    def v_max_set(self, change):
        vmax = change['new']
        # This lookup will always succeed because the min/max for
        # any field is set from the default (if necessary) in change_radar
        # upon any change to the radar field
        old_min, old_max = self.value_range_for_field[self.field_type.value]
        self.ncgridlasso.vmax = vmax
        self.plot_update()

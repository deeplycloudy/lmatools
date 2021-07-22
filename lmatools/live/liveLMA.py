from __future__ import absolute_import
from __future__ import print_function
import sys
import numpy as np
from collections import deque

import threading

# 
# # Header data is packed in a structure:
# #
# # typedef struct {
# # time_t header_second;
# # int    num_sources;
# # short  num_passes;
# # short  num_stations;
# # short  min_num_stations;
# # float  max_chisq1;
# # float  max_chisq2;
# # short  decimate;
# # int    dec_window;
# # }
# 
# # This should be 34 bytes, but Python says only 32 bytes
# # are recevied, so I fudged things by making max_chisq2 a short ('H')
# header = struct.Struct('Q I H H H f H H I')
# header_size = 32;
# 
# # Source data is packed in a C structure:
# #
# # typedef struct {
# #   float lat, lon, alt; 
# #   float x, y;            // Distance from center of coordinate system
# #   float t;               // Fractions of seconds after header_second
# #   float p_dbw;
# #   float chisq;
# #   short num_stations;
# #   unsigned short stationmask;
# # }
# 
# source = struct.Struct('f f f f f f f f H H')
# source_size = 36;
try:
    import websocket
except ImportError as ierr:
    print("Please pip install websocket-client or get it manually from from https://github.com/liris/websocket-client")
    raise(ierr)
    
class WebsocketClient(object):
    
    def __init__(self, host=None):
        self.host=host
        self.ws=None
    
    def on_error(self, ws, error):
        print(error)

    def on_close(self, ws):
        print("### closed ###")
        self.ws = None
        
    def connect(self, on_message=None):
    
        self.ws = websocket.WebSocketApp(self.host,
                                on_message = on_message,
                                on_error = self.on_error,
                                on_close = self.on_close)
                                # ws.on_open = on_open
        self.ws.run_forever()
    
    
header_dtype = [ 
    ('header_second', 'u8'), 
    ('num_sources', 'i4'), 
    ('num_passes', 'i2'),
    ('num_stations', 'i2'), ('min_num_stations', 'i2'),
    ('max_chisq1', 'f4'), ('max_chisq2', 'f4'),
    ('decimate', 'i2'), ('dec_window', 'i4'),
    ]
header_size = 32


source_dtype = [ 
    ('lat', 'f4'), ('lon', 'f4'), ('alt', 'f4'),
    ('x', 'f4'), ('y', 'f4'), # Distance from center of coordinate system
    ('t', 'f4'), # Fractions of seconds after header_second
    ('power', 'f4'), ('chi2', 'f4'),
    ('stations', 'i2'), ('mask', 'u2'),
    ]
source_size = 36


class LiveLMAprinter(object):
    def __init__(self):
        pass
        
    def show(self, data):
        """docstring for show"""
        for v in data:
            # print("{0}: {1}, {2}".format(data['t'], data['lat'], data['lon']))
            print(("{1}, {2}".format(data['t'], data['lat'], data['lon'])))

def redraw(canvas):
    canvas.draw()

class LiveLMAplanview(object):
    """ Plot incoming LMA points on a matplotlib axis"""

    def __init__(self, scatter_artist, age=600):
        """ scatter_artist is a matplotlib line artist. age in seconds"""
        self._data = deque([],age)
        self.scatter_artist = scatter_artist
    

    def show(self, data):
        canvas = self.scatter_artist.figure.canvas
        
        # print data, self.scatter_artist.get_array()
        
        self._data.append(data)        
        a=np.hstack([d for d in self._data if (d.shape[0] > 0)])

        xy = np.vstack((a['lon'], a['lat'])).T
        t = a['t']
        # print xy
        # print "is xy"
        # print t
        # print "is t"
        self.scatter_artist.set_offsets(xy)
        self.scatter_artist.set_array(t)
        self.scatter_artist.set_clim(t.min(),t.max())
        # print "about to draw"
        canvas.draw()

class LiveLMAController(object):
    """ Data model for messages recieved from a LiveLMA WebSocket server """
    
    def __init__(self):
        """ self.dtype_mapping is from a key in the source_dtype to a 
            value in self.dtype
        """
        self.views=[]
        
        # Mapping from 
        self.dtype_mapping = {
            'lat':'lat',
            'lon':'lon',
            'alt':'alt',
            # time is handled separately to combine basetime from header
            'power':'power',
            'chi2':'chi2',
            'stations':'stations',
            'mask':'mask',
            }
        self.dtype = [ 
                ('lat', 'f4'), ('lon', 'f4'), ('alt', 'f4'),
                ('t', 'f8'), # Fractions of seconds after header_second
                ('power', 'f4'), ('chi2', 'f4'),
                ('stations', 'i2'), ('mask', 'u2'),
                ]

    def on_message(self, ws, message):
        """ ws is the low-level websocket object, message is the binary string of data"""
        

        header = np.fromstring(message[0:header_size], dtype=header_dtype)
        sources = np.fromstring(message[header_size:], dtype=source_dtype, count=header['num_sources'])
        
        #print("{0} new, {1} stations".format(header['num_sources'][0], header['num_stations'][0]))
                
        data = np.empty_like(sources, dtype=self.dtype) 
        
        # This may not be the optimal precision; perhaps should convert to seconds of day.
        # See time handling in archive_to_LiveLMA for examples
        data['t'] = sources['t'].astype('f8')+header['header_second'].astype('f8')
        
        for k, v in self.dtype_mapping.items():
            data[v] = sources[k]
            
        for view in self.views:
            view.show(header, data)
        
        
    
    
if __name__ == '__main__':
    try:
        server = sys.argv[1]
    except IndexError as e:
        print("Please provide a URL for the LiveLMA server, which typically begins with ws://\nUsage: python liveLMA.py server_url")
        raise(e)

    
    liveDataController = LiveLMAController()

    printer = LiveLMAprinter()
    liveDataController.views.append(printer)
    
    # ----- Start a threaded WebSocket Client -----
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.axis(( -103.5, -99.5, 31.5, 35.5))
    scatter = ax.scatter([33,35], [-102,-101], c=[0,1], s=9, marker='.', linewidth=0)
    plt.draw()
    plotter = LiveLMAplanview(scatter)
    liveDataController.views.append(plotter)
    
    # ----- Start a threaded WebSocket Client -----
    # websocket.enableTrace(True)
    client = WebsocketClient(host=server)
    # client.connect(on_message=liveDataController.on_message)
    sock_thr = threading.Thread(target=client.connect, kwargs={'on_message':liveDataController.on_message})
    sock_thr.daemon=True
    sock_thr.start()
    
    plt.show()
    
    
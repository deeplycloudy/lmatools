from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import logging
import re

from scipy.spatial import Delaunay, ConvexHull
from scipy.special import factorial
from scipy.spatial.qhull import QhullError

from lmatools.lasso import empirical_charge_density as cd

# logger = logging.getLogger('FlashAutorunLogger')
class Flash(object):
    def __init__(self, points):
        self.points = points

class FlashMetadata(object):

    def __init__(self, headerText):
        #Search the header for info on how the data is written

        self.header = headerText

        isColumnHeaderLine = r"^Data:(.*)"
        matchDataFormatLine = re.compile(isColumnHeaderLine, re.IGNORECASE | re.MULTILINE)

        isDataStartTime = r"^Data start time:(.*)"
        matchDataStartTimeLine = re.compile(isDataStartTime, re.IGNORECASE | re.MULTILINE)

        secAnalyzedLine = r"^Number of seconds analyzed:(.*)"
        matchSecAnalyzedLine = re.compile(secAnalyzedLine, re.IGNORECASE | re.MULTILINE)


        startTimeMatch = matchDataStartTimeLine.search(headerText)
        if startTimeMatch:
            #Looking to deal with something like: " 06/28/04 23:50:00"
            dateAndTime = startTimeMatch.group(1).split()
            self.startmonth, self.startday, self.startyear = [ int(datePart) for datePart in dateAndTime[0].split('/') ]
            self.starthour, self.startminute, self.startsecond = [ int(timePart) for timePart in dateAndTime[1].split(':') ]
            if self.startyear < 1000 and self.startyear > 70:
                self.startyear += 1900
            else:
                self.startyear += 2000

        secAnalyzedMatch=matchSecAnalyzedLine.search(headerText)
        if secAnalyzedMatch:
            self.sec_analyzed = int(secAnalyzedMatch.group(1))


        # formatMatch=matchDataFormatLine.search(headerText)
        # if formatMatch:
        #     columns = formatMatch.group(1).split(',')
        #     self.columns = [columnName.strip() for columnName in columns]

 
def barotropic_rho(z):
    rho = 1.225e9 #kg/km^3
    H = 8. #km
    return rho*np.exp(-z/H) 
    
def poly_area(x,y):
    """ Calculate the area of a non-self-intersecting planar polygon.
        x0y1 - x0y0 + x1y2 - x2y1 + x2y3 - x2y2 + ... + xny0 - x0yn 
    """
    det = x[:-1]*y[1:] - x[1:]*y[:-1] # determinant
    area = det.sum()
    area += x[-1]*y[0] - x[0]*y[-1] # wrap-around terms in determinant
    area *= 0.5
    return area
    

def hull_volume(xyz):
    """ Calculate the volume of the convex hull of 3D (X,Y,Z) LMA data.
        xyz is a (N_points, 3) array of point locations in space. """
    assert xyz.shape[1] == 3
        
    tri = Delaunay(xyz[:,0:3])
    vertices = tri.points[tri.vertices]
    
    # This is the volume formula in 
    # https://github.com/scipy/scipy/blob/master/scipy/spatial/tests/test_qhull.py#L106
    # Except the formula needs to be divided by ndim! to get the volume, cf., 
    # http://en.wikipedia.org/wiki/Simplex#Geometric_properties
    # Credit Pauli Virtanen, Oct 14, 2012, scipy-user list
    
    q = vertices[:,:-1,:] - vertices[:,-1,None,:]
    simplex_volumes = (1.0 / factorial(q.shape[-1])) * np.fromiter(
           (np.linalg.det(q[k,:,:]) for k in range(tri.nsimplex)) , dtype=float)
    # print vertices.shape # number of simplices, points per simplex, coords
    # print q.shape
    
    # The simplex volumes have negative values since they are oriented 
    # (think surface normal direction for a triangle
    volume=np.sum(np.abs(simplex_volumes))
    return volume, vertices, simplex_volumes

##############ADDED 01/05/2017 ###############
def energy(area, separation, zinit, constant, eta):
    #Charge separation computed from 27th and 73rd percentiles of 
    #flash altitude source locations - marks where the most sources are typically
    #found from synthetic flashes generated in the NSSL COMMAS.
    #
    #eta = 0.01 is recommended and is a ballpark neutrlization efficiency as found in Salinas et al. [In Progress - 060220]

    distance = separation #np.abs(random)
    density  = cd.rho_retrieve(area, distance, zinit, separation, False, None) #None - No constant charge density specified
    rho,w    = density.calculate()
    return(eta*w)
##############################################   

def calculate_flash_stats(flash, min_pts=2):
    logger = logging.getLogger('FlashAutorunLogger')
    
    Re = 6378.137*1000;           #Earth's radius in m
    pi = np.pi
    
    flash.pointCount = flash.points.shape[0]
    
    fl_id = np.unique(flash.points['flash_id'])
    assert (fl_id.shape[0] == 1)
    flash.id = fl_id[0]
    
    
    lat = np.asarray(flash.points['lat'],dtype=float) 
    lon = np.asarray(flash.points['lon'],dtype=float) 
    alt = np.asarray(flash.points['alt'], dtype=float)
    # 
    # # mean location of all points
    latavg, lonavg, altavg = lat.mean(), lon.mean(), alt.mean()
    x = Re * (np.radians(lonavg) - np.radians(lon)) * np.cos(np.radians(latavg))
    y = Re * (np.radians(latavg) - np.radians(lat))
    z = altavg - alt
    # r_sq = x**2.0 + y**2.0 + z**2.0
    # sigma_sq = r_sq.sum()/r_sq.shape[0]
    # sigma = np.std(r_sq)

    separation     = np.abs(np.percentile(alt,73) - np.percentile(alt,27))
    flash_init_idx = np.argmin(flash.points['time'])
    zinit = alt[flash_init_idx] #in meters
    area  = 0.0
    eta   = 0.01

    if flash.pointCount > 2:
        try:
            # find the convex hull and calculate its area
            cvh = ConvexHull(np.vstack((x,y)).T)
            # NOT cvh.area - it is the perimeter in 2D. 
            # cvh.area is the surface area in 3D.
            area = cvh.volume
        except IndexError:
            # tends to happen when a duplicate point causes the point count to
            # drop to 2, leading to a degenerate polygon with no area
            logger.warning('Setting area to 0 for flash with points %s, %s' % (x, y))
            area=0.0
        except KeyError:
            # hull indexing has problems here
            logger.warning('Setting area to 0 for flash with points %s, %s' % (x, y))
            area=0.0
           
    if area == 0.0:
        energy_estimate = 0.
    else:        
        energy_estimate = energy(area, separation, zinit, False, eta)
            
    volume = 0.0
    if flash.pointCount > 4:
        # Need five points to make at least one tetrahedron.
        try:
            volume, vertices, simplex_volumes = hull_volume(np.vstack((x,y,z)).T)
        except QhullError:
            # this can happen with a degenerate first simplex - all points are
            # coplanar to machine precision. Try again, after adding a tiny amount
            # to the first point.
            print("Perturbing one source to help triangulation for flash with {0} points".format(flash.pointCount))
            # we can tolerate perturbing by no more than 1 m
            machine_eps = 1.0 # np.finfo(x.dtype).eps
            perturb = 2*machine_eps*np.random.random(size=3)-machine_eps
            x[0] += perturb[0]
            y[0] += perturb[1]
            z[0] += perturb[2]

    flash_init_idx = np.argmin(flash.points['time'])
    
    ###ROUGH APPROXIMATION FOR NOW: #######################
    air_density = barotropic_rho(alt[flash_init_idx]*1e-3)
    if volume == 0.:
        specific_energy = 0.
    else:
        specific_energy = energy_estimate / ((volume / 1.0e9) * air_density)
    #######################################################
    
    flash.start   = flash.points[flash_init_idx]['time']
    flash.end     = flash.points['time'].max()
    flash.duration = flash.end - flash.start
    flash.area    = area / 1.0e6  # km^2, 1000x1000
    flash.initLat = lat[flash_init_idx] 
    flash.initLon = lon[flash_init_idx]
    flash.initStd = 0.0
    flash.initAlt = alt[flash_init_idx]
    flash.initPts = (int(flash_init_idx),)
    flash.ctralt  = altavg
    flash.ctrlat  = latavg
    flash.ctrlon  = lonavg
    flash.volume  = volume / 1.0e9 # km^3, 1000x1000x1000 m
    #CHANGED 03-20-17
    flash.total_energy    = energy_estimate    #flash.energy ---> flash.tot_energy
    flash.specific_energy = specific_energy  #flash.tot_energy ---> flash.specific_energy

import numpy as np
import logging

# logger = logging.getLogger('FlashAutorunLogger')

    
def poly_area(x,y):
    """ Calculate the area of a non-self-intersecting planar polygon.
        x0y1 - x0y0 + x1y2 - x2y1 + x2y3 - x2y2 + ... + xny0 - x0yn 
    """
    det = x[:-1]*y[1:] - x[1:]*y[:-1] # determinant
    area = det.sum()
    area += x[-1]*y[0] - x[-1]*y[0] # wrap-around terms in determinant
    area *= 0.5
    return area
    
def calculate_flash_stats(flash, min_pts=2, max_pts=10):
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

    area = 0.0
    if flash.pointCount > 2:
        try:
            # find the convex hull and calculate its area
            # There's a patch to deal with duplicate points - mpl has it as of 0.98.6svn
            # scipy scikits.delaunay has it after r745
            from matplotlib import delaunay
            tri = delaunay.Triangulation(x,y)
            hull = tri.hull
            area = poly_area(tri.x[hull], tri.y[hull])
        except ImportError:
            logger.error("*** Flash area not calculated - requires delaunay from or matplotlib or scipy")
        except IndexError:
            # tends to happen when a duplicate point causes the point count to
            # drop to 2, leading to a degenerate polygon with no area
            logger.warning('Setting area to 0 for flash with points %s, %s' % (x, y))
            area=0.0
        except KeyError:
            # hull indexing has problems here
            logger.warning('Setting area to 0 for flash with points %s, %s' % (x, y))
            area=0.0
            
    flash.start   = flash.points['time'].min()
    flash.end     = flash.points['time'].max()
    flash.duration = flash.end - flash.start
    flash.area    = area / 1e6  # km^2, 1000x1000
    flash.initLat = lat[0] 
    flash.initLon = lon[0]
    flash.initStd = 0.0
    flash.initAlt = alt[0]
    flash.initPts = (0,)
    flash.ctralt  = altavg
    flash.ctrlat  = latavg
    flash.ctrlon  = lonavg
from __future__ import absolute_import
import pyproj as proj4
from numpy import *
from numpy.linalg import norm

# def radians(degrees):
    # return deg2rad(asarray(degrees))
    # return array(degrees) * pi / 180.0

# def degrees(radians):
    # return rad2deg(asarray(radians))
    # return array(radians) * 180.0 / pi


class CoordinateSystem(object):
    """The abstract coordinate system handling provided here works as follows.

    Each coordinate system must be able to convert data to a common coordinate system, which is chosen to be ECEF cartesian.
    data -> common system
    common system -> dislpay coordinates
    This is implemented by the fromECEF and toECEF methods in each coordinate system object.
    User code is responsible for taking data in its native coord system,
        transforming it using to/fromECEF using the a coord system appropriate to the data, and then
        transforming that data to the final coordinate system using another coord system.

    Subclasses should maintain an attribute ERSxyz that can be used in
        transformations to/from an ECEF cartesian system, e.g.
        >>> self.ERSxyz = proj4.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
        >>> self.ERSlla = proj4.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
        >>> projectedData = proj4.transform(self.ERSlla, self.ERSxyz, lat, lon, alt )
    The ECEF system has its origin at the center of the earth, with the +Z toward the north pole,
        +X toward (lat=0, lon=0), and +Y right-handed orthogonal to +X, +Z

    Depends on pyproj, http://code.google.com/p/pyproj/ to handle the ugly details of
    various map projections, geodetic transforms, etc.

    "You can think of a coordinate system as being something like character encodings,
    but messier, and without an obvious winner like UTF-8." - Django OSCON tutorial, 2007
    http://toys.jacobian.org/presentations/2007/oscon/tutorial/
    """

    # WGS84xyz = proj4.Proj(proj='geocent',  ellps='WGS84', datum='WGS84')

    def coordinates():
        """Return a tuple of standarized coordinate names"""
        raise NotImplemented

    def fromECEF(self, x, y, z):
        """Take ECEF x, y, z values and return x, y, z in the coordinate system defined by the object subclass"""
        raise NotImplemented

    def toECEF(self, x, y, z):
        """Take x, y, z in the coordinate system defined by the object subclass and return ECEF x, y, z"""
        raise NotImplemented


class GeographicSystem(CoordinateSystem):
    """
    Coordinate system defined on the surface of the earth using latitude,
    longitude, and altitude, referenced by default to the WGS84 ellipse.

    Alternately, specify the ellipse shape using an ellipse known
    to pyproj, or [NOT IMPLEMENTED] specify r_equator and r_pole directly.
    """
    def __init__(self, ellipse='WGS84', datum='WGS84',
                 r_equator=None, r_pole=None):
        if (r_equator is not None) | (r_pole is not None):
            if r_pole is None:
                r_pole=r_equator
            self.ERSlla = proj4.Proj(proj='latlong', a=r_equator, b=r_pole)
        else:
            # lat lon alt in some earth reference system
            self.ERSlla = proj4.Proj(proj='latlong', ellps=ellipse, datum=datum)
        self.ERSxyz = proj4.Proj(proj='geocent', ellps=ellipse, datum=datum)
    def toECEF(self, lon, lat, alt):
        lat = atleast_1d(lat) # proj doesn't like scalars
        lon = atleast_1d(lon)
        alt = atleast_1d(alt)
        if (lat.shape[0] == 0): return lon, lat, alt # proj doesn't like empties
        projectedData = array(proj4.transform(self.ERSlla, self.ERSxyz, lon, lat, alt ))
        if len(projectedData.shape) == 1:
            return projectedData[0], projectedData[1], projectedData[2]
        else:
            return projectedData[0,:], projectedData[1,:], projectedData[2,:]

    def fromECEF(self, x, y, z):
        x = atleast_1d(x) # proj doesn't like scalars
        y = atleast_1d(y)
        z = atleast_1d(z)
        if (x.shape[0] == 0): return x, y, z # proj doesn't like empties
        projectedData = array(proj4.transform(self.ERSxyz, self.ERSlla, x, y, z ))
        if len(projectedData.shape) == 1:
            return projectedData[0], projectedData[1], projectedData[2]
        else:
            return projectedData[0,:], projectedData[1,:], projectedData[2,:]


class MapProjection(CoordinateSystem):
    """Map projection coordinate system. Wraps proj4, and uses its projecion names. Defaults to
        equidistant cylindrical projection
    """

    def __init__(self, projection='eqc', ctrLat=None, ctrLon=None, ellipse='WGS84', datum='WGS84', **kwargs):
        self.ERSxyz = proj4.Proj(proj='geocent', ellps=ellipse, datum=datum)
        self.projection = proj4.Proj(proj=projection, ellps=ellipse, datum=datum, **kwargs)
        self.ctrLat=ctrLat
        self.ctrLon=ctrLon
        self.ctrAlt=0.0
        self.geoCS = GeographicSystem()
        self.cx, self.cy, self.cz = 0, 0, 0
        self.cx, self.cy, self.cz = self.ctrPosition()

    def ctrPosition(self):
        if (self.ctrLat != None) & (self.ctrLon != None):
            ex, ey, ez = self.geoCS.toECEF(self.ctrLon, self.ctrLat, self.ctrAlt)
            cx, cy, cz = self.fromECEF(ex, ey, ez)
        else:
            cx, cy, cz = 0, 0, 0
        return cx, cy, cz

    def toECEF(self, x, y, z):
        x += self.cx
        y += self.cy
        z += self.cz
        projectedData = array(proj4.transform(self.projection, self.ERSxyz, x, y, z ))
        if len(projectedData.shape) == 1:
            px, py, pz = projectedData[0], projectedData[1], projectedData[2]
        else:
            px, py, pz = projectedData[0,:], projectedData[1,:], projectedData[2,:]
        return px, py, pz

    def fromECEF(self, x, y, z):
        projectedData = array(proj4.transform(self.ERSxyz, self.projection, x, y, z ))
        if len(projectedData.shape) == 1:
            px, py, pz = projectedData[0], projectedData[1], projectedData[2]
        else:
            px, py, pz = projectedData[0,:], projectedData[1,:], projectedData[2,:]
        return px-self.cx, py-self.cy, pz-self.cz

class PixelGrid(CoordinateSystem):
    def __init__(self, lons, lats, lookup, x, y, alts=None, geosys=None):
        """
        Coordinate system for arbitrary pixel coordinates in a 2D pixel array.
        Arguments:
        lons: 2D array of longitudes of pixel centers
        lats: 2D array of longitudes of pixel centers
        alts: 2D array of longitudes of pixel centers. If None, zeros are assumed.
        Each array is of shape (nx, ny) with pixel coordinate (x=0, y=0)
            corresponding to grid index [0, 0]

        lookup is an object with a method 'query' that accepts a single argument,
        a (N,2) array of lats, lons and returns pixel IDs that can be used to
        index lons and lats, as well as the distances between the pixel centers
        and the queried locations. X and Y flattened arrays of pixel coordinates
        that align with indices of the flattened lon and lat arrays used to
        create the lookup table.
        >>> test_events = np.vstack([(-101.5, 33.5), (-102.8, 32.5), (-102.81,32.5)])
        >>> distances, idx = lookup.query(test_events)
        >>> loni, lati = lons[X[idx], Y[idx]], lats[X[idx], Y[idx]]
        An instance of sklearn.neighbors.KDTree is one such lookup.

        If geosys is provided, it should be an instance of GeographicSystem;
        otherwise a GeographicSystem instance with default arguments is created.

        When converting toECEF, which accepts pixel coordinates,
        the z pixel coordinate is ignored, as it has no meaning.
        When converting fromECEF, zeros in the shape of x are returned as the z
        coordinate.

        """
        if geosys is None:
            self.geosys = GeographicSystem()
        else:
            self.geosys = geosys
        self.lookup = lookup
        self.x = x
        self.y = y
        self.lons = lons
        self.lats = lats
        if alts is None:
            alts = zeros_like(lons)
        self.alts = alts

    def toECEF(self, x, y, z):
        x = x.astype('int64')
        y = y.astype('int64')
        lons = self.lons[x, y]
        lats = self.lats[x, y]
        alts = self.alts[x, y]
        return self.geosys.toECEF(lons, lats, alts)

    def fromECEF(self, x, y, z):
        lons, lats, alts = self.geosys.fromECEF(x, y, z)
        locs = vstack((lons.flatten(), lats.flatten())).T
        if locs.shape[0] > 0:
            distances, idx = self.lookup.query(locs)
        else:
            idx = []
        x = squeeze(self.x[idx])
        y = squeeze(self.y[idx])
        return x, y, zeros_like(x)

class GeostationaryFixedGridSystem(CoordinateSystem):

    def __init__(self, subsat_lon=0.0, subsat_lat=0.0, sweep_axis='y',
                 sat_ecef_height=35785831.0,
                 ellipse='WGS84', datum='WGS84'):
        """
        Satellite height is with respect to the ellipsoid. Fixed grid
        coordinates are in radians.
        """
        # self.ECEFxyz = proj4.Proj(proj='cart', ellps=ellipse)#, datum=datum)
        self.ECEFxyz = proj4.Proj(proj='geocent', ellps=ellipse)#, datum=datum)
        self.fixedgrid = proj4.Proj(proj='geos', lon_0=subsat_lon,
            lat_0=subsat_lat, h=sat_ecef_height, x_0=0.0, y_0=0.0,
            units='m', sweep=sweep_axis, ellps=ellipse)
        self.h=sat_ecef_height

    def toECEF(self, x, y, z):
        X, Y, Z = x*self.h, y*self.h, z*self.h
        return proj4.transform(self.fixedgrid, self.ECEFxyz, X, Y, Z)

    def fromECEF(self, x, y, z):
        X, Y, Z = proj4.transform(self.ECEFxyz, self.fixedgrid, x, y, z)
        return X/self.h, Y/self.h, Z/self.h

# class AltitudePreservingMapProjection(MapProjection):
#     def toECEF(self, x, y, z):
#         px, py, pz = super(AltitudePreservingMapProjection, self).toECEF(x, y, z)
#         return px, py, z
#
#     def fromECEF(self, x, y, z):
#         px, py, pz = super(AltitudePreservingMapProjection, self).fromECEF(x, y, z)
#         return px, py, z

class RadarCoordinateSystem(CoordinateSystem):
    """
        Converts spherical (range, az, el) radar coordinates to lat/lon/alt, and then to ECEF.

        An earth's effective radius of 4/3 is assumed to correct for atmospheric refraction.
    """

    def __init__(self, ctrLat, ctrLon, ctrAlt, datum='WGS84', ellps='WGS84', effectiveRadiusMultiplier=4./3.):
        self.ctrLat = float(ctrLat)
        self.ctrLon = float(ctrLon)
        self.ctrAlt = float(ctrAlt)
        self.datum=datum
        self.ellps=ellps

        self.lla = proj4.Proj(proj='latlong', ellps=self.ellps, datum=self.datum)
        self.xyz = proj4.Proj(proj='geocent', ellps=self.ellps, datum=self.datum)

        self.Requator, foo1, foo2 = proj4.transform(self.lla,self.xyz,0,0,0) # Equatorial radius  - WGS-84 value = 6378137.0
        foo1, foo2, self.Rpolar = proj4.transform(self.lla,self.xyz,0,90,0) # Polar radius  - WGS-84 value = 6356752.314
        self.flattening = (self.Requator-self.Rpolar)/self.Requator

        self.eccen = (2.0-self.flattening)*self.flattening   # First eccentricity squared - WGS-84 value = 0.00669437999013
        self.effectiveRadiusMultiplier = effectiveRadiusMultiplier

    def getGroundRangeHeight(self, r, elevationAngle):
        """Convert slant range (along the beam) and elevation angle into
        ground range (great circle distance) and height above the earth's surface
        Follows Doviak and Zrnic 1993, eq. 2.28."""

        #Double precison arithmetic is crucial to proper operation.
        lat = self.ctrLat * pi / 180.0
        elev = array(elevationAngle * pi / 180.0, dtype='float64')
        slantr = array(r, dtype='float64')

        #figure out earth's radius at radar's lat ... non-spherical earth model
        e2 = self.eccen           # First eccentricity squared - WGS-84 value = 0.00669437999013
        a = self.Requator         # Equatorial radius  - WGS-84 value = 6378137.0
        Rearth = a/sqrt(1-e2*(sin(lat))**2) # radius of curvature

        Rprime = self.effectiveRadiusMultiplier * Rearth

        # Eqns 2.28b,c in Doviak and Zrnic 1993
        # Radar altitude is tacked on at the end, which isn't part of their derivation. At 100 km, it's
        #   worth < 10 m range error total for a radar at 500 m MSL. For 250 m gate spacing (typical at S-band),
        #   this is not too important.
        h = sqrt(slantr**2.0 + Rprime**2.0 + 2*slantr*Rprime*sin(elev)) - Rprime
        s = Rprime * arcsin( (slantr*cos(elev)) / (Rprime + h) )

        h += self.ctrAlt

        return s, h

    def getSlantRangeElevation(self, groundRange, z):
        """Convert ground range (great circle distance) and height above
        the earth's surface to slant range (along the beam) and elevation angle.
        Follows Doviak and Zrnic 1993, eq. 2.28"""

        lat = self.ctrLat * pi / 180.0

        #figure out earth's radius at radar's lat ... non-spherical earth model
        e2 = self.eccen           # First eccentricity squared - WGS-84 value = 0.00669437999013
        a = self.Requator         # Equatorial radius  - WGS-84 value = 6378137.0
        Rearth = a/sqrt(1-e2*(sin(lat))**2) # radius of curvature

        Rprime = self.effectiveRadiusMultiplier * Rearth

        h = array(z - self.ctrAlt, dtype='float64')
        s = array(groundRange, dtype='float64')

        # Use law of cosines (Side-Angle-Side triangle theorem) with
        # R', R'+h as sides and s/R' as the angle to get slant range
        r  = sqrt(Rprime**2.0 + (Rprime+h)**2.0 - 2*(Rprime+h)*Rprime*cos(s/Rprime))
        # Inverse of eq. 2.28c in Doviak and Zrnic 1993
        # Will return NaN for r=0, and only positive angles
        el = atleast_1d(arccos((Rprime+h) * sin(s/Rprime) / r))
        # Below gives all negative angles
        # el = arcsin((Rprime+h) * sin(s/Rprime) / r) - (pi/2.0)

        # If elevation angle is negative, the triangle will be acute
        acute = atleast_1d( (Rprime+h)*(Rprime+h) < (Rprime*Rprime + r*r) )
        el[acute] *= -1
        el *= 180.0 / pi

        return r, el

    def toLonLatAlt(self, r, az, el):
        """Convert slant range r, azimuth az, and elevation el to ECEF system"""
        geoSys = GeographicSystem()
        geodetic = proj4.Geod(ellps=self.ellps)

        try:
            n = max((az.size, r.size))
        except AttributeError:
            n = max((len(az), len(r)))

        dist, z = self.getGroundRangeHeight(r,el)
        lon, lat, backAz = geodetic.fwd([self.ctrLon]*n, [self.ctrLat]*n, az, dist)
        return lon, lat, z

    def toECEF(self, r, az, el):
        geoSys = GeographicSystem()
        lon, lat, z = self.toLonLatAlt(r, az, el)
        return geoSys.toECEF(lon, lat, z.ravel())

    def fromECEF(self, x, y, z):
        """Convert ECEF system to slant range r, azimuth az, and elevation el"""
        # x = np.atleast1d(x)
        geoSys = GeographicSystem()
        geodetic = proj4.Geod(ellps=self.ellps)

        try:
            n = x.size
        except AttributeError:
            n = len(x)

        lon, lat, z = geoSys.fromECEF(x, y, z)
        radarToGateAz, gateToRadarAz, dist = geodetic.inv([self.ctrLon]*n, [self.ctrLat]*n, lon, lat)
        az = array(radarToGateAz)   #radarToGateAz may be a list.
        # change negative azimuths to positive
        az[az < 0.0] += 360.0

        #have height, ground range, azimuth. need to get elev angle and slant range from ground range and height
        r, el = self.getSlantRangeElevation(dist, z)

        return r, az, el

class TangentPlaneCartesianSystem(object):
    """ TODO: This function needs to be updated to inherit from CoordinateSystem

    """

    def __init__(self, ctrLat=0.0, ctrLon=0.0, ctrAlt=0.0):
        self.ctrLat = float(ctrLat)
        self.ctrLon = float(ctrLon)
        self.ctrAlt = float(ctrAlt)

        ERSlla = proj4.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
        ERSxyz = proj4.Proj(proj='geocent',  ellps='WGS84', datum='WGS84')
        self.centerECEF = array(proj4.transform(ERSlla, ERSxyz, ctrLon, ctrLat, ctrAlt))

        #location of point directly above local center
        aboveCenterECEF = array(proj4.transform(ERSlla, ERSxyz, ctrLon, ctrLat, self.ctrAlt+1e3))

        #normal vector to earth's surface at the center is the local z direction
        n = aboveCenterECEF - self.centerECEF
        n = n / norm(n)
        localz = n[:,None] #make a column vector

        # n (dot) x = d defines a plane for normal vector n and position vector x on the plane
        d = dot(n, aboveCenterECEF)

        #north = array((northx, northy, northz))

        #http://www.euclideanspace.com/maths/geometry/elements/plane/index.htm
        #matrix to project point onto a plane defined by the normal vector n.
        P = identity(3,float) - transpose(vstack((n,n,n))) * vstack((n,n,n))

        # Point just to the north of the center on earth's surface, projected onto the tangent plane
        # This calculation seems like it should only be done with latitude/north since the local x
        #   direction curves away along a non-straight line when projected onto the plane
        northCenterECEF = array(proj4.transform(ERSlla, ERSxyz, self.ctrLon, self.ctrLat+1.01, self.ctrAlt))
        localy = dot(P, northCenterECEF[:,None] )
        localy = localy / norm(localy)


        #local x is y (cross) z to get an orthogonal system
        localx = transpose(cross(localy.transpose(), localz.transpose()))
        localx = localx / norm(localx)


        ECEFx = array((1.0, 0.0, 0.0))[:,None]
        ECEFy = array((0.0, 1.0, 0.0))[:,None]
        ECEFz = array((0.0, 0.0, 1.0))[:,None]

        #
        # Calculate the transformation matrix TM to go from
        #   the earth-centered earth-fixed (ECEF) system to the local tangent plane system
        # http://www.spenvis.oma.be/spenvis/help/background/coortran/coortran.html, http://mathworld.wolfram.com/DirectionCosine.html
        # (X1, X2, X3) are the direction cosines of the X-direction of the b-system, expressed in function of X, Y and Z of the a-system
        # b system = local tangent plane system     a system = ECEF system
        # [vb_x]   [[x1, x2, x3]  [va_x
        # [vb_y] =  [y1, y2, y3]   va_y
        # [vb_z]    [z1, z2, z3]]  va_z]
        # va = transpose(M) vb
        x1 = dot(localx.transpose(), ECEFx) # / abs(localx) ... don't need since normalized
        x2 = dot(localx.transpose(), ECEFy)
        x3 = dot(localx.transpose(), ECEFz)
        y1 = dot(localy.transpose(), ECEFx) # / abs(localx) ... don't need since normalized
        y2 = dot(localy.transpose(), ECEFy)
        y3 = dot(localy.transpose(), ECEFz)
        z1 = dot(localz.transpose(), ECEFx) # / abs(localx) ... don't need since normalized
        z2 = dot(localz.transpose(), ECEFy)
        z3 = dot(localz.transpose(), ECEFz)
        self.TransformToLocal = array([[x1, x2, x3],
                                       [y1, y2, y3],
                                       [z1, z2, z3]]).squeeze()

    def fromECEF(self, x, y, z):
        """ Transforms 1D arrays of ECEF x, y, z to the local tangent plane system"""
        data = vstack((x, y, z))
        tpXYZ = self.toLocal(data)
        if len(tpXYZ.shape) == 1:
            tpX, tpY, tpZ = tpXYZ[0], tpXYZ[1], tpXYZ[2]
        else:
            tpX, tpY, tpZ = tpXYZ[0,:], tpXYZ[1,:], tpXYZ[2,:]
        return tpX, tpY, tpZ

    def toECEF(self, x, y, z):
        """ Transforms 1D arrays of x, y, z in the local tangent plane system to ECEF"""
        data = vstack((x, y, z))
        ecXYZ = self.fromLocal(data)
        if len(ecXYZ.shape) == 1:
            ecX, ecY, ecZ = ecXYZ[0], ecXYZ[1], ecXYZ[2]
        else:
            ecX, ecY, ecZ = ecXYZ[0,:], ecXYZ[1,:], ecXYZ[2,:]
        return ecX, ecY, ecZ

    def toLocal(self, data):
        """Transforms 3xN array of data (position vectors) in the ECEF system to
           the local tangent plane cartesian system.

           Returns another 3xN array.
        """
        return array( [ dot(self.TransformToLocal, (v-self.centerECEF)[:,None])
                        for v in data[0:3,:].transpose()]
                    ).squeeze().transpose()

    def fromLocal(self, data):
        """Transforms 3xN array of data (position vectors) in the local tangent
           plane cartesian system to the ECEF system.

           Returns another 3xN array.
        """
        #Transform from local to ECEF uses transpose of the TransformToLocal matrix
        return array( [ (dot(self.TransformToLocal.transpose(), v) + self.centerECEF)
                        for v in data[0:3,:].transpose()]
                    ).squeeze().transpose()

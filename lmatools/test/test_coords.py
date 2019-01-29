from lmatools.coordinateSystems import GeostationaryFixedGridSystem, GeographicSystem

from numpy.testing import assert_allclose

def test_fixed_grid_GOESR():
    """ Tests from the GOES-R PUG Volume 3, L1b data """
    sat_lon_nadir = -75.0
    goes_sweep = 'x' # Meteosat is 'y'
    ellipse = 'GRS80'
    datum = 'WGS84'
    sat_ecef_height=35786023.0
        
    geofixcs = GeostationaryFixedGridSystem(subsat_lon=sat_lon_nadir,
                   ellipse=ellipse, datum=datum, sweep_axis=goes_sweep,
                   sat_ecef_height=sat_ecef_height)
    latloncs = GeographicSystem(ellipse=ellipse, datum=datum)
    

    test_lat = 33.846162
    test_lon = -84.690932
    test_alt = 0.0
    
    test_fixx = -0.024052
    test_fixy = 0.095340
    test_fixz = 0.0
    
    atol = 1e-6
    
    # Test forward from geodetic
    X,Y,Z = latloncs.toECEF(test_lon, test_lat, test_alt)
    x, y, z = geofixcs.fromECEF(X, Y, Z)
    assert_allclose(test_fixx, x, rtol=atol)
    assert_allclose(test_fixy, y, rtol=atol)
    assert_allclose(test_fixz, z, rtol=atol)
    
    # print(test_fixx, test_fixy, test_fixz)
    # print(x,y,z)

    # Test inverse from fixed grid angle
    X,Y,Z = geofixcs.toECEF(test_fixx, test_fixy, test_fixz)
    x, y, z = latloncs.fromECEF(X, Y, Z)
    assert_allclose(test_lon, x, atol=atol)
    assert_allclose(test_lat, y, atol=atol)
    assert_allclose(test_alt, z, atol=atol)
    
    # print(test_lon, test_lat, test_alt)
    # print(x,y,z)
    

if __name__ == '__main__':
    test_fixed_grid_GOESR()
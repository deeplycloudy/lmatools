"""
Specification of the GOES Fixed Grid.
"""
from lmatools.coordinateSystems import GeostationaryFixedGridSystem, GeographicSystem

goesr_nadir_lon = {
    'east':-75.0,
    'west':-137.0,
    'test':-89.5,
}

goesr_full = { # all values in radians on the fixed grid
    # The spans here match the values in L1b PUG table 5.1.2.7-4 but
    # when combined with the 10 km resolution give two (2) fewer
    # pixels than claimed in table 5.1.2.6
    'spanEW': 0.151872*2.0,
    'spanNS': 0.151872*2.0,
    'centerEW': 0.0,
    'centerNS': 0.0,
}

goesr_conus = { # all values in radians on the fixed grid
    'spanEW': 0.14,
    'spanNS': 0.084,
}

goesr_meso = {
    'spanEW': 0.028,
    'spanNS': 0.028,
}

goeseast_full = goesr_full.copy()
goeswest_full = goesr_full.copy()
goestest_full = goesr_full.copy()

goeseast_meso = goesr_meso.copy()
goeswest_meso = goesr_meso.copy()
goestest_meso = goesr_meso.copy()

goeseast_conus = goesr_conus.copy()
goeseast_conus.update({
    'centerEW': -0.031360,
    'centerNS': 0.086240,
})

goeswest_conus = goesr_conus.copy()
goeswest_conus.update({
    'centerEW': 0.000000,
    'centerNS': 0.086240,
})

goestest_conus = goesr_conus.copy()
goestest_conus.update({
    'centerEW': -0.005040,
    'centerNS': 0.084560,
})

goesr_resolutions = {
    '0.5km': 14e-6,
    '1.0km': 28e-6,
    '2.0km': 56e-6,
    '4.0km': 112e-6,
    '8.0km': 224e-6, # not technically in the spec, but SYMMETRY.
    '10.0km': 280e-6,
    '100.0km': 2800e-6,
}


def get_GOESR_coordsys(sat_lon_nadir = -75.0):
    """
    Values from the GOES-R PUG Volume 3, L1b data

    Returns geofixcs, grs80lla: the fixed grid coordinate system and the
    latitude, longitude, altitude coordinate system referenced to the GRS80
    ellipsoid used by GOES-R as its earth reference.
    """
    goes_sweep = 'x' # Meteosat is 'y'
    ellipse = 'GRS80'
    datum = 'WGS84'
    sat_ecef_height=35786023.0
    geofixcs = GeostationaryFixedGridSystem(subsat_lon=sat_lon_nadir,
                   ellipse=ellipse, datum=datum, sweep_axis=goes_sweep,
                   sat_ecef_height=sat_ecef_height)
    grs80lla = GeographicSystem(ellipse='GRS80', datum='WGS84')
    return geofixcs, grs80lla

def get_GOESR_grid(position='east', view='full', resolution='2.0km'):
    """
    This helper function returns specifications of the GOES-R fixed grids.

    position is ['east'|'west'|'test']
    resolution is ['0.5km'|'1.0km'|'2.0km'|'4.0km'|'10.0km']
    view is ['full'|'conus'|'meso']

    returns a dictionary with the keys 'resolution', 'spanEW', 'spanNS',
        'pixelsEW', 'pixelsNS', 'nadir_lon', and (for all non-meso sectors)
        'centerEW', 'centerNS'.
        values are radians, except for the integer
    """
    assert position in ['east', 'west', 'test']
    assert view in ['full', 'conus', 'meso']
    assert resolution in goesr_resolutions.keys()

    sector = view

    # namespace = __import__(__name__)
    view = globals()['goes'+position+"_"+view].copy()

    # according to scott's suggestion, change goes_west conus domain to
    # 'mirror of East conus (from Eric) + west boundary of original WEST conus'
    # by feng Jun 11, 2019
    if position == 'west' and sector == 'conus':
       view['spanEW'] = 0.172732
       view['centerEW'] = 0.015652

    view['resolution'] = goesr_resolutions[resolution]
    view['pixelsEW'] = int(view['spanEW']/view['resolution'])
    view['pixelsNS'] = int(view['spanNS']/view['resolution'])
    view['nadir_lon'] = goesr_nadir_lon[position]

    return view

if __name__ == '__main__':
    for pos in ['east', 'west', 'test']:
        for view in ['full', 'conus', 'meso']:
            for resolution in ['0.5km', '1.0km', '2.0km', '4.0km', '8.0km', '10.0km']:
                print('-----\n', pos, view, resolution)
                for k, v in get_GOESR_grid(pos, view, resolution).items():
                    print(k, v)

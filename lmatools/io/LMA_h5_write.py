from __future__ import absolute_import
import tables as T
import numpy as np
import logging

# Write this many flashes at a time
flash_chunk = 100 # was 10000 rows

logger = logging.getLogger('FlashAutorunLogger')

def countBits(values):
    # bit shifting routines are in numpy 1.4
    from numpy import array, left_shift, right_shift
    
    v = array(values).astype('uint32')
    
    # Bit counting for a 32 bit unsigned integer.
    # there is a fully generic version of this method at
    # http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
    # Binary magic numbers method, also found in pages 187-188 of Software Optimization Guide for AMD Athlon 64 and Opteron Processors.

    # The C version is:
    # v = v - ((v >> 1) & 0x55555555);                    # reuse input as temporary
    # v = (v & 0x33333333) + ((v >> 2) & 0x33333333);     # temp
    # c = ((v + (v >> 4) & 0xF0F0F0F) * 0x1010101) >> 24; # count

    fives  = int('0x55555555', base=16)
    threes = int('0x33333333', base=16)
    effs   = int('0xF0F0F0F',  base=16)
    ones   = int('0x1010101',  base=16)
    
    v = v - ( (right_shift(v, 1)) & fives);                        # reuse input as temporary
    v = (v & threes) + ( (right_shift(v,2)) & threes);             # temp
    c = right_shift(((v + (right_shift(v,4)) & effs) * ones), 24); # count

    return c

def mask_strings_to_stations(masks):
    mask_int = np.fromiter((int(v,16) for v in masks), int)
    stations = countBits(mask_int)
    return stations


class Flash(T.IsDescription):
    flash_id = T.Int32Col()   # Flash ID
    n_points = T.Int16Col()
    start    = T.Float64Col() # flash start
    duration = T.Float32Col() # flash duration    
    ctr_lat  = T.Float32Col() # centroid latitude
    ctr_lon  = T.Float32Col() # centroid longitude
    ctr_alt  = T.Float32Col() # centroid altitude
    init_lat = T.Float32Col() # initation latitude
    init_lon = T.Float32Col() # initation longitude
    init_alt = T.Float32Col() # initation altitude
    init_pts = T.StringCol(256) # Indices of the points (first point in flash is id=0) used to calc init location
    area     = T.Float32Col() # area of convex hull of the points comprising the flash
    volume   = T.Float32Col()
    #Changed variable names: 03-20-17 ---> Check make_grids.py in case of inconsistancies.
    total_energy   =  T.Float32Col()    #Energy
    specific_energy = T.Float32Col()    #tot_energy
        

def write_h5(outfile, flashes, metadata, orig_LMA_file, mask_length):
    # file_parts = lma.filename.split('.')[0].split('_')
    # time_code  = 'LMA_'+'_'.join(file_parts[1:])
    # LMA_090329_180000_3600
    m=metadata # flashes[0].metadata
    time_code = 'LMA_%s%02d%02d_%02d%02d%02d_%d' % (str(m.startyear)[-2:], m.startmonth, m.startday,
                                m.starthour, m.startminute, m.startsecond, m.sec_analyzed)
    # orig_columns_LYLOUT = m.columns
    
    class Event(T.IsDescription):
        # ascii line for all these data are 64 bytes + 4 for station count + 4 for charge = 72
        
        # 64+32*5+8+8+4*8 = 272 bytes per record

        time = T.Float64Col()       # Seconds elapsed since start of day
        lat  = T.Float32Col()       # Decimal latitude
        lon  = T.Float32Col()       # Decimal longitude
        alt  = T.Float32Col()       # Altitude, km MSL, WGS84
        chi2 = T.Float32Col()       # Chi-squared solution quality
        power= T.Float32Col()       # Radiated power
        stations = T.UInt8Col()     # Station count
        charge   = T.Int8Col()      # Inferred storm charge
        flash_id    = T.Int32Col()     # Flash ID
        mask     = T.StringCol(mask_length)   # Station mask

    h5file = T.open_file(outfile, mode='w', title='Flash-sorted New Mexico Tech LMA Data')
    group  = h5file.create_group('/', 'events', 'Analyzed detected events')
    table  = h5file.create_table(group, time_code, Event, time_code)
    fl_group = h5file.create_group('/', 'flashes', 'Sorted LMA flash data')
    fl_table  = h5file.create_table(fl_group, time_code, Flash, time_code)
    
    table.attrs.header = m.header
    table.attrs.filename = orig_LMA_file
    table.attrs.start_time = ( m.startyear,
                               m.startmonth,
                               m.startday,
                               m.starthour,
                               m.startminute,
                               m.startsecond,
                             )    
    event = table.row
    fl_meta = fl_table.row
        
    for i, flash in enumerate(flashes):
        fl_meta['flash_id'] = flash.id
        fl_meta['n_points'] = flash.pointCount
        fl_meta['start']    = flash.start
        fl_meta['duration'] = flash.duration
        fl_meta['area']     = flash.area
        fl_meta['ctr_lat']  = flash.ctrlat
        fl_meta['ctr_lon']  = flash.ctrlon
        fl_meta['ctr_alt']  = flash.ctralt
        fl_meta['init_lat'] = flash.initLat
        fl_meta['init_lon'] = flash.initLon
        fl_meta['init_alt'] = flash.initAlt
        fl_meta['volume']   = flash.volume
        #Changed 03-20-17
        fl_meta['total_energy']   = flash.total_energy  #flash.energy
        fl_meta['specific_energy'] = flash.specific_energy #flash.tot_energy
        
        # init_table = fl_meta['init_points']
        # init_event = init_table.row
        # for idx in flash.initPts:
            # init_event.append(idx)
        fl_meta.append()

        for j in range(flash.pointCount):
        # for j, record in enumerate(flash.points):
            event['time']     = flash.points['time'][j]
            event['lat']      = flash.points['lat'][j]
            event['lon']      = flash.points['lon'][j]
            event['alt']      = flash.points['alt'][j]
            event['chi2']     = flash.points['chi2'][j]
            event['power']    = flash.points['power'][j]
            event['stations'] = flash.points['stations'][j]
            event['mask']     = flash.points['mask'][j]
            # if 'charge' in record.keys():
                # event['charge'] = record['charge']
            event['flash_id'] = flash.id
            event.append()
        
        # write a flash_chunk flashes at a time
        if (i % flash_chunk == 0) & (i!=0):
            logger.debug(i)
            table.flush()
            
    table.flush()    
    
    
    # # read back in all the event station masks and save the number of stations
    # masks = np.fromiter((row['mask'] for row in table), dtype='|S4')
    # stations = mask_strings_to_stations(masks)
    # # n_stations = stations.shape[0]
    # table.cols.stations[:] = stations.astype('uint8')
    
    
    table.flush()
    if len(flashes) == 0:
        i = -1
        
    logger.info('total flashes: %d' % (i+1,))
    
    h5file.close()
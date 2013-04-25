import logging

import numpy as np
from numpy.lib.recfunctions import append_fields
            
from sklearn.cluster import DBSCAN

from lmatools.stream.subset import coroutine, stream, chunk
from lmatools.flashsort.autosort.LMAarrayFile import LMAdataFile
from lmatools.coordinateSystems import GeographicSystem

from flash_stats import calculate_flash_stats, Flash, FlashMetadata



@coroutine
def cluster_chunk_pairs(clustered_output_target, min_points=10):
    db = DBSCAN(eps=1.0, min_samples=min_points, metric='euclidean')
    
    """Receive chunks, and process overlapping pairs"""
    chunk1 = (yield)
    try:
        while True:
            chunk2 = (yield)
            len1 = chunk1.shape[0]
            len2 = chunk2.shape[0]
            print len1+len2
            
            # do stuff with chunk 1 and 2
            clusters = db.fit(np.vstack((chunk1, chunk2)))
            labels = clusters.labels_.astype(int)
            
            clustered_output_target.send((chunk1, labels[:len1]))
            
            # pull data out of chunk2 that was clustered as part of chunk 1
            chunk1_labelset = set(labels[:len1])
            if -1 in chunk1_labelset:
                chunk1_labelset.remove(-1) # remove the singleton cluster ID - we want to retain these from chunk 2.
            clustered_in_chunk2 = np.fromiter( ( True if label in chunk1_labelset else False for i,label in enumerate(labels[len1:])) , dtype=bool)
            clustered_output_target.send((chunk2[clustered_in_chunk2], labels[len1:][clustered_in_chunk2]))  
            residuals = chunk2[clustered_in_chunk2==False]
            
            # optimization TODO: pull clusters out of chunk 2 whose final point is greater 
            # than the distance threshold from the end of the second chunk interval. They're already clustered
            # and don't need to be clustered again.
            
            # prepare for another chunk
            if len(residuals) == 0:
                residuals = chunk1[0:0,:] # empty array that preserves the number of dimensions in the data vector - no obs.
            del chunk1
            chunk1 = np.asarray(residuals)
            del residuals
    except GeneratorExit:
        clusters = db.fit(chunk1)
        labels = clusters.labels_.astype(int)
        clustered_output_target.send((chunk1, labels))
        
# @coroutine
# def chunk_printer():
#     total = 0
#     try:
#         n_last = 0
#         while True:
#             v = (yield)
#             total += v.shape[0]
#             n,x = v[:,-1].min(), v[:,-1].max()
#             print v.shape, n,x,x-n_last
#             n_last=n
#             del v
#     except GeneratorExit:
#         print total
        
@coroutine
def aggregate_ids(target):
    unique_labels = set([-1])
    total = 0
    point_labels = []
    # all_v = []
    try:
        n_last = 0
        while True:
            (v, orig_labels) = (yield)
            labels = np.atleast_1d(orig_labels).copy()
            if len(unique_labels) > 0:
                # Only add those labels that represent valid clusters (nonnegative) to the unique set.
                # Make sure labels increment continuously across all chunks received
                nonsingleton = (labels >= 0)
                labels[nonsingleton] = labels[nonsingleton] + (max(unique_labels) + 1)
                
            for l in set(labels):
                unique_labels.add(l)

            point_labels.append(labels)
            total += v.shape[0]

            del v, orig_labels, labels
    except GeneratorExit:
        print "done with {0} total points".format(total)
        point_labels = np.concatenate(point_labels)
        print "sending {0} total points".format(total)
        target.send((unique_labels, point_labels))
        print "sent {0} total points".format(total)

@coroutine
def create_flash_objs(lma, good_data):
    """ lma is an LMAdataFile object. Its data instance gets overwritten with the sorted, qc'd, flash_id'd data.
    
        very similar to collect_output in autorun_mflash
        
    """
    logger = logging.getLogger('FlashAutorunLogger')
    
    
    try:
        while True:
            (unique_labels, point_labels) = (yield)
            
            # add flash_id column
            data = append_fields(good_data, ('flash_id',), (point_labels,))
            
            
            # In the case of no data in the file, lma.data.shape will have length zero, i.e., a 0-d array
            if len(data.shape) == 0:
                # No data
                flashes = []
            else:
                # work first with non-singleton flashes to have strictly positive flash ids
                print data.shape
                singles = (data['flash_id'] == -1)
                non_singleton = data[ np.logical_not(singles) ]
                print non_singleton['flash_id'].shape
                order = np.argsort(non_singleton['flash_id'])
                
                ordered_data = non_singleton[order]
                flid = ordered_data['flash_id']                
                max_flash_id = flid[-1]
                try:
                    assert max_flash_id == max(unique_labels)
                except AssertionError:
                    print "Max flash ID {0} is not as expected from unique labels {1}".format(max_flash_id, max(unique_labels))
                    
                boundaries, = np.where(flid[1:]-flid[:-1])    # where indices are nonzero
                boundaries = np.hstack(([0], boundaries+1))
                
                max_idx = len(flid) #- 1
                slice_lower_edges = tuple(boundaries)
                slice_upper_edges = slice_lower_edges[1:] + (max_idx,)
                slices = zip(slice_lower_edges, slice_upper_edges)

                flashes = [ Flash(ordered_data[slice(*sl)]) for sl in slices ]
                
                print "finished non-singletons"
                
                # now deal with the nonsingleton points. Each singleton point will have a high flash_id,
                # starting with the previous maximum flash id.
                singleton = data[singles]
                n_singles = singleton.shape[0]

                # this operation works on a view of the original data array, so it modifies the original data array
                singleton['flash_id'] += max_flash_id + 1 + np.arange(n_singles, dtype=int)
                
                singleton_flashes = [ Flash(singleton[i:i+1]) for i in range(n_singles)]
                
                print "finished singletons"
                
                flashes += singleton_flashes
                
                logtext = "Calculating flash initation, centroid, area, etc. for %d flashes" % (len(flashes), )
                logger.info(logtext)
                # print flashes[0].points.dtype
                for fl in flashes:
                    header = ''.join(lma.header)
                    fl.metadata = FlashMetadata(header)
                    calculate_flash_stats(fl)
                    # logger.info(fl.points.shape[0])
                logger.info('finished setting flash metadata')
                
                lma.raw_data = lma.data
                lma.data = data
                assert (lma.data['flash_id'].min() >=0) # this should be true since the singletons were modified in the original data array above
                lma.sort_status = 'got some flashes'
                
    except GeneratorExit:
        lma.flash_objects = flashes


def cluster(a_file, output_path, outfile, params, logger, min_points=1, **kwargs):
    """
    There is no intermediate ASCII output or temporary file for this code, since all data remains as native Python objects.
    
    """
    logger = logging.getLogger('FlashAutorunLogger')
    
    
    if 'mask_length' in params:
        mask_length = params['mask_length']
    else:
        mask_length = 4
    
    lma=LMAdataFile(a_file, mask_length = mask_length)
    # for line in lma.header:
        # print line

    ctr_lat, ctr_lon, ctr_alt =  params['ctr_lat'], params['ctr_lon'], 0.0

    good = (lma.stations >= params['stations'][0]) & (lma.chi2 <= params['chi2'][1]) 
    if 'alt' in params:
        good = good & (lma.alt < params['alt'][1])
    
    
    data = lma.data[good]
    geoCS = GeographicSystem()
    X,Y,Z = geoCS.toECEF(data['lon'], data['lat'], data['alt'])
    Xc, Yc, Zc = geoCS.toECEF( ctr_lon, ctr_lat, ctr_alt)
    X, Y, Z = X-Xc, Y-Yc, Z-Zc
    
    print "sorting {0} total points".format(data.shape[0])

    D_max, t_max = 3.0e3, 0.15 # m, s

    X_vector = np.hstack((X[:,None],Y[:,None],Z[:,None])) / D_max
    T_vector = data['time'][:,None] / t_max
    XYZT = np.hstack((X_vector, T_vector-T_vector.min()))
    
    lma.sort_status = 'in process'
    
    # Maximum 3 s flash length, normalized to the time separation scale
    chunker = chunk(XYZT[:,-1].min(), 3.0/.15, cluster_chunk_pairs(aggregate_ids(create_flash_objs(lma, data)), min_points=min_points) )
    stream(XYZT.astype('float32'),chunker)
    
    print lma.sort_status
    print len(lma.flash_objects)
                    
    return lma, lma.flash_objects
    

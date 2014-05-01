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
    #db = DBSCAN(eps=1.0, min_samples=min_points, metric='euclidean')
    
    
    """Receive chunks, and process overlapping pairs"""
    chunk1, id1 = (yield)
    try:
        while True:
            chunk2, id2 = (yield)
            len1 = chunk1.shape[0]
            len2 = chunk2.shape[0]
            if len2 == 0:
                conc = chunk1
                concID = id1
                chunk2 = chunk1[0:0,:]
                id2 = id1[0:0]
            elif len1 == 0:
                conc = chunk2
                concID = id2
                chunk1 = chunk2[0:0,:]
                id1 = id2[0:0]
            else:
                print id1.shape, id2.shape
                conc = np.vstack((chunk1, chunk2)) 
                concID = np.concatenate((id1, id2))
                        
            # do stuff with chunk 1 and 2
            
            db = DBSCAN(eps=1.0, min_samples=min_points, metric='euclidean')
            clusters = db.fit(conc)
            labels = clusters.labels_.astype(int)
            
            # defer sending these in one bundle ... need to ensure all labels
            # from this run of DBSCAN stay together
            # clustered_output_target.send((chunk1, labels[:len1]))
            
            # pull data out of chunk2 that was clustered as part of chunk 1
            chunk1_labelset = set(labels[:len1])
            if -1 in chunk1_labelset:
                chunk1_labelset.remove(-1) # remove the singleton cluster ID - we want to retain these from chunk 2.
            clustered_in_chunk2 = np.fromiter( ( True if label in chunk1_labelset else False for i,label in enumerate(labels[len1:])) , dtype=bool)
            clustered_in_chunk1 = np.ones(chunk1.shape[0], dtype = bool)
            clustered_mask = np.hstack((clustered_in_chunk1, clustered_in_chunk2))
            bundle_chunks = conc[clustered_mask,:]
            bundle_IDs = concID[clustered_mask]
            
            bundle_labels = np.concatenate((labels[:len1], labels[len1:][clustered_in_chunk2]))
            assert bundle_chunks.shape[0] == bundle_labels.shape[0]
            clustered_output_target.send((bundle_chunks, bundle_labels, bundle_IDs))
            del bundle_chunks, bundle_labels
            # clustered_output_target.send((chunk2[clustered_in_chunk2], labels[len1:][clustered_in_chunk2]))  
            
            residuals = conc[clustered_mask==False,:]
            # Because we pull some points from chunk2 and combine them with
            # flashes that started in chunk1, the data are now out of their
            # original order. Therefore, send along the data IDs that go with the
            # pulled points so that the original order is still tracked.
            residualIDs = concID[clustered_mask==False]

            # optimization TODO: pull clusters out of chunk 2 whose final point is greater 
            # than the distance threshold from the end of the second chunk interval. They're already clustered
            # and don't need to be clustered again.
            
            # prepare for another chunk
            if len(residuals) == 0:
                residuals = chunk1[0:0,:] # empty array that preserves the number of dimensions in the data vector - no obs.
                residualIDs = id1[0:0]
            del chunk1, id1
            chunk1 = np.asarray(residuals)
            id1 = np.asarray(residualIDs)
            del residuals, residualIDs
    except GeneratorExit:
        if chunk1.shape[0] != 0:
            db = DBSCAN(eps=1.0, min_samples=min_points, metric='euclidean')
            clusters = db.fit(chunk1)
            labels = clusters.labels_.astype(int)
            clustered_output_target.send((chunk1, labels, id1))
        clustered_output_target.close()
        
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
    all_IDs = []
    # all_v = []
    try:
        # n_last = 0
        while True:
            (v, orig_labels, IDs) = (yield)
            labels = np.atleast_1d(orig_labels).copy()
            if len(unique_labels) > 0:
                # Only add those labels that represent valid clusters (nonnegative) to the unique set.
                # Make sure labels increment continuously across all chunks received
                nonsingleton = (labels >= 0)
                labels[nonsingleton] = labels[nonsingleton] + (max(unique_labels) + 1)
                
            for l in set(labels):
                unique_labels.add(l)

            all_IDs.append(np.asarray(IDs))
            point_labels.append(labels)
            total += v.shape[0]

            del v, orig_labels, labels, IDs
    except GeneratorExit:
        print "done with {0} total points".format(total)
        point_labels = np.concatenate(point_labels)
        all_IDs = np.concatenate(all_IDs)
        print "sending {0} total points".format(total)
        target.send((unique_labels, point_labels, all_IDs))
        print "sent {0} total points".format(total)
        target.close()

@coroutine
def create_flash_objs(lma, good_data):
    """ lma is an LMAdataFile object. Its data instance gets overwritten with the sorted, qc'd, flash_id'd data.
    
        very similar to collect_output in autorun_mflash
        
    """
    logger = logging.getLogger('FlashAutorunLogger')
    
    
    try:
        while True:
            (unique_labels, point_labels, all_IDs) = (yield)
            
            # add flash_id column
            empty_labels = np.empty_like(point_labels)
            data = append_fields(good_data, ('flash_id',), (empty_labels,))

            # all_IDs gives the index in the original data array to
            # which each point_label corresponds
            data['flash_id'][all_IDs] = point_labels
            
            # In the case of no data in the file, lma.data.shape will have
            # length zero, i.e., a 0-d array
            if len(data.shape) == 0:
                # No data
                flashes = []
            else:
                # work first with non-singleton flashes 
                # to have strictly positive flash ids
                print data.shape
                singles = (data['flash_id'] == -1)
                non_singleton = data[ np.logical_not(singles) ]
                print non_singleton['flash_id'].shape
                order = np.argsort(non_singleton['flash_id'])
                
                ordered_data = non_singleton[order]
                flid = ordered_data['flash_id']                
                if (flid.shape[0]>0):
                    max_flash_id = flid[-1]
                else: 
                    max_flash_id = 0
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
                
                # Now deal with the nonsingleton points. 
                # Each singleton point will have a high flash_id,
                # starting with the previous maximum flash id.
                singleton = data[singles]
                n_singles = singleton.shape[0]

                # this operation works on a view of the original data array, 
                # so it modifies the original data array
                singleton['flash_id'] += max_flash_id + 1 + np.arange(n_singles, dtype=int)
                
                singleton_flashes = [ Flash(singleton[i:i+1]) for i in range(n_singles)]
                
                data[singles] = singleton
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

    IDs = np.arange(X.shape[0])
    X_vector = np.hstack((X[:,None],Y[:,None],Z[:,None])) / D_max
    T_vector = data['time'][:,None] / t_max
    XYZT = np.hstack((X_vector, T_vector))
    
    lma.sort_status = 'in process'
    
    # Maximum 3 s flash length, normalized to the time separation scale

    flash_object_maker = create_flash_objs(lma, data)
    label_aggregator = aggregate_ids(flash_object_maker)
    clusterer = cluster_chunk_pairs(label_aggregator, min_points=min_points)
    chunker = chunk(XYZT[:,-1].min(), 3.0/.15,  clusterer)
    stream(XYZT.astype('float64'), IDs,chunker)
    flash_object_maker.close()
    
    # These are handled by target.close in each coroutine's GeneratorExit handler
    # clusterer.close()
    # label_aggregator.close()
    # flash_object_maker.close()
    
    print lma.sort_status
    print len(lma.flash_objects)
                    
    return lma, lma.flash_objects
    

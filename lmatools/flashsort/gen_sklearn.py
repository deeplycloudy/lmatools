from __future__ import absolute_import
from __future__ import print_function
import logging

import numpy as np
from numpy.lib.recfunctions import append_fields
            
from sklearn.cluster import DBSCAN

from lmatools.coordinateSystems import GeographicSystem

from lmatools.flashsort.flash_stats import calculate_flash_stats, Flash


def gen_stream(vec, IDs): #<1>
    for v, vi in zip(vec, IDs):
        yield (v, vi)
        
def reset_buffer():
    buf = []
    return buf, buf.append
    
def gen_chunks(stream, start_time, max_duration, t_idx=-1):
    """ Generator function that consumes a stream of points, one at a 
        time, and their unique index. These points are bundled together 
        into a chunks of length max_duration along the time coordinate.
    
        For each point vector v, the time coordinate is given by v[t_idx]
    """
    next_time = start_time + max_duration
    v_buffer, append = reset_buffer() # slight optimization since attr lookup is avoided
    i_buffer, append_idx = reset_buffer()
    for v, vi in stream:
        append(v)
        append_idx(vi)
        t = v[t_idx]
        if t >= next_time:
            yield (np.asarray(v_buffer), np.asarray(i_buffer))
            v_buffer, append = reset_buffer()
            i_buffer, append_idx = reset_buffer()
            next_time = t+max_duration
    yield (np.asarray(v_buffer), np.asarray(i_buffer))
    

    
class ChunkedFlashSorter(object):
    """ 
    Sort LMA data from points to flashes using many small chunks
    of points. Allows for algorithms that do not scale efficiently with
    large numbers of points. 
    
    The __init__ and geo_to_cartesian
    methods are more generally useful, and could be factored out into a
    generic flash sorting class.
    
    The actual clustering algorithm must be implemented in identify_clusters.
    A prototype method is provided below which indicates the necessary call
    signature.
    """
    def __init__(self, params, min_points=1, **kwargs):
        """ 
            params: dictionary of parameters used to perform data QC and clustering
            min_points: the minimum number of points allowed in a cluster
            
        """
        self.logger = logging.getLogger('FlashAutorunLogger')
        self.logger.info('%s', params)
        
        self.params = params
        self.min_points = min_points
        self.ctr_lat, self.ctr_lon, self.ctr_alt =  (
                            params['ctr_lat'], params['ctr_lon'], 0.0)
        

    def geo_to_cartesisan(self, lon, lat, alt):
        """ Convert lat, lon in degrees and altitude in meters to 
            Earth-centered, Earth-fixed cartesian coordinates. Translate
            to coordinate center location. Returns X,Y,Z in meters.
        """
        geoCS = GeographicSystem()
        X,Y,Z = geoCS.toECEF(lon, lat, alt)
        Xc, Yc, Zc = geoCS.toECEF(self.ctr_lon, self.ctr_lat, self.ctr_alt)
        X, Y, Z = X-Xc, Y-Yc, Z-Zc
        return (X, Y, Z)
    
    def identify_clusters(self, data):
        """ For data with shape (N, D) in D dimensions, return 
            a vector of labels of length N. 
    
            min_points is the minimum number of points required to form a 
            a cluster. For the DBSCAN algorithm, this is min_samples for
            a core cluster.

            This function adopts the convention that clusters labeled
            with an ID of -1 are singleton points not belonging to a 
            cluster, consistent with the convention of sklearn.cluster.DBSCAN
        """
        err = "Please create a new subclass and implement this method"
        raise NotImplementedError(err)
    
    def gen_cluster_chunk_pairs(self, stream):
        """ Generator function that consumes a stream of chunks of data, 
            and processes overlapping pairs. The stream is to consist of 
            tuples of (chunk, pt_id), where pt_id is a unique index for 
            each vector in chunk. 
            Chunk is of shape (N, D) for N point vectors in D dimensions
            pt_id has shape (N,)
    
            Calls self.identify_clusters, which returns a vector N labels.
            The labels are presumed to adopt the convention that clusters labeled
            with an ID of -1 are singleton points not belonging to a 
            cluster, consistent with the convention of sklearn.cluster.DBSCAN
        """
        chunk1, id1 = next(stream)
        for chunk2, id2 in stream:
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
                print(id1.shape, id2.shape)
                conc = np.vstack((chunk1, chunk2)) 
                concID = np.concatenate((id1, id2))

            # do stuff with chunk 1 and 2
            labels = self.identify_clusters(conc)
        
            # defer sending these in one bundle ... need to ensure all labels
            # from this run of clustering stay together
            # clustered_output_target.send((chunk1, labels[:len1])) IS BAD
        
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
            yield (bundle_chunks, bundle_labels, bundle_IDs)
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
        if chunk1.shape[0] != 0:
            labels = self.identify_clusters(chunk1)
            yield (chunk1, labels, id1)
        
    def aggregate_ids(self, stream):
        """ Final step in streamed clustering: consume clustered output from
            one or more chunks of data, ensuring that the IDs increment 
            across chunk boundaries.
        """
        # TODO: remove v from loop below; not needed. 
        unique_labels = set([-1])
        total = 0
        point_labels = []
        all_IDs = []
        # all_v = []
            # n_last = 0
        for (v, orig_labels, IDs) in stream:
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
        print("done with {0} total points".format(total))
        if total == 0:
            point_labels = np.asarray(point_labels, dtype=int)
            point_labels = np.asarray(all_IDs, dtype=int)
        else:
            point_labels = np.concatenate(point_labels)
            all_IDs = np.concatenate(all_IDs)
        print("returning {0} total points".format(total))
        return (unique_labels, point_labels, all_IDs)

    def create_flash_objs(self, lma, good_data, unique_labels, point_labels, all_IDs):
        """ lma is an LMADataset object. Its data instance gets overwritten 
            with the qc'd, flash_id'd data, and it gains a flashes attribute
            with a list of flash objects resulting from flash sorting.
        
            all_IDs gives the index in the original data array to
            which each point_label corresponds. 
            unique_labels is the set of all labels produced by previous stages 
            in the flash sorting algorithm, including a -1 ID for all singleton flashes.
        """
        logger = self.logger
                
        # add flash_id column
        empty_labels = np.empty_like(point_labels)
        # placing good_data in a list due to this bug when good_data has length 1
        # http://stackoverflow.com/questions/36440557/typeerror-when-appending-fields-to-a-structured-array-of-size-one
        if 'flash_id' not in good_data.dtype.names:
            data = append_fields([good_data], ('flash_id',), (empty_labels,))
        else:
            data = good_data.copy()

        # all_IDs gives the index in the original data array to
        # which each point_label corresponds
        data['flash_id'][all_IDs] = point_labels
    
        # In the case of no data, lma.data.shape will have
        # length zero, i.e., a 0-d array
    
        if (len(data.shape) == 0) | (data.shape[0] == 0):
            # No data
            flashes = []
        else:
            # work first with non-singleton flashes 
            # to have strictly positive flash ids
            print(data.shape)
            singles = (data['flash_id'] == -1)
            non_singleton = data[ np.logical_not(singles) ]
            print(non_singleton['flash_id'].shape)
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
                print("Max flash ID {0} is not as expected from unique labels {1}".format(max_flash_id, max(unique_labels)))
            
            boundaries, = np.where(flid[1:]-flid[:-1])    # where indices are nonzero
            boundaries = np.hstack(([0], boundaries+1))
        
            max_idx = len(flid) #- 1
            slice_lower_edges = tuple(boundaries)
            slice_upper_edges = slice_lower_edges[1:] + (max_idx,)
            slices = list(zip(slice_lower_edges, slice_upper_edges))

            flashes = [ Flash(ordered_data[slice(*sl)]) for sl in slices ]
        
            print("finished non-singletons")
        
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
            print("finished singletons")
        
            flashes += singleton_flashes
        
            logtext = "Calculating flash initation, centroid, area, etc. for %d flashes" % (len(flashes), )
            logger.info(logtext)
            # print flashes[0].points.dtype
            for fl in flashes:
                # header = ''.join(lma.header)
                fl.metadata = lma.metadata #FlashMetadata(header)
                calculate_flash_stats(fl)
                # logger.info(fl.points.shape[0])
            logger.info('finished setting flash metadata')
        
            lma.raw_data = lma.data
            lma.data = data
            assert (lma.data['flash_id'].min() >=0) # this should be true since the singletons were modified in the original data array above
        
        lma.flashes = flashes
        
    def perform_chunked_clustering(self, XYZT, ptIDs, chunk_duration):
        """ Perform clustering of a 4D data vector in overlapping chunks of
            data, 
    
        XYZT: (N,4) array of N 4D point vectors
        ptIDs: array of N unique identifiers of each point vector.
        chunk_duration: duration of chunk in the units along the T coordinate
    
        """
        if XYZT.shape[0] < 1:
            # no data, so minimum time is zero. Assume nothing is done with the 
            # data, so that time doesn't matter. No flashes can result.
            t_min = 0
        else:
            t_min = XYZT[:,-1].min()

        point_stream = gen_stream(XYZT.astype('float64'), ptIDs)
        chunk_stream = gen_chunks(point_stream, t_min, chunk_duration) 
        cluster_stream = self.gen_cluster_chunk_pairs(chunk_stream)
        unique_labels, point_labels, all_IDs = self.aggregate_ids(cluster_stream)
        return unique_labels, point_labels, all_IDs


    def cluster(self, dataset, **kwargs):
        """ Cluster an lmatools LMADataset provided in the dataset argument.
            Basic filtering of the data is performed by calling the filter_data
            method of the dataset, which returns a filtered data array. The 
            params are provided by the argument used to initialize this class.
            
            This method modifies dataset as a side effect.
        """
    
        data = dataset.filter_data(self.params)
        print("sorting {0} total points".format(data.shape[0]))
        
        # Transform to cartesian coordiantes
        X, Y, Z = self.geo_to_cartesisan(data['lon'], data['lat'], data['alt'])
        
        # Assemble a normalized data vector using flash sorting parameters
        D_max, t_max = (self.params['distance'], # meters
                        self.params['thresh_critical_time']) # seconds
        duration_max = self.params['thresh_duration'] # seconds

        IDs = np.arange(X.shape[0]) # Vector of unique point IDs
        X_vector = np.hstack((X[:,None],Y[:,None],Z[:,None])) / D_max
        T_vector = data['time'][:,None] / t_max
        XYZT = np.hstack((X_vector, T_vector))
    
        # Perform chunked clustering of the data
        normed_chunk_duration = duration_max/t_max
        unique_labels, point_labels, all_IDs = self.perform_chunked_clustering(XYZT, IDs, normed_chunk_duration)

        # Calculate flash metadata and store it in flash objects
        # This should be factored out into somehting that modifies the data table
        # and something that creates the flash objects themselves
        self.create_flash_objs(dataset, data, unique_labels, point_labels, all_IDs)


class DBSCANFlashSorter(ChunkedFlashSorter):
    def identify_clusters(self, data):
        """ For data with shape (N, D) in D dimensions, return 
            a vector of labels of length N. 
        
            min_points is the minimum number of points required to form a 
            a cluster. For the DBSCAN algorithm, this is min_samples for
            a core cluster.
    
            This function adopts the convention that clusters labeled
            with an ID of -1 are singleton points not belonging to a 
            cluster, consistent with the convention of sklearn.cluster.DBSCAN
        """
        db = DBSCAN(eps=1.0, min_samples=self.min_points, metric='euclidean')
        clusters = db.fit(data)
        labels = clusters.labels_.astype(int)
        return labels
        

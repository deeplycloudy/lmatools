""" 

Lightning Mapping Array datasets are often comprised of as many as N=O(10^5)
points per minute. It is often of interest to group points into flashes. Some
of the clustering algorithms in Python's scikit-learn package require a
pairwise distance calculation in the form of a N^2 distance matrix, which
exceeds available memory. At the same time, it is clear that there is an upper
physical bound to flash size, such that the data could be processed in chunks,
with

The purpose of this code is to demonstrate a method for streaming points, from
which are built a distance matrix, in this case internally by the DBSCAN
algorithm. Given an assumed maximum duration of the physical process t_max,
processing the first half of a 2*t_max length buffer is guaranteed to capture
all points. After clustering, all points and belonging to cluster members from the
first half chunk are pruned, retaining residual points within the second half.

Probably should take a careful look at memory leaks; based on past experience
a few more del statements might be needed.

"""


from __future__ import absolute_import
from __future__ import print_function
import numpy as np

from sklearn.cluster import DBSCAN

from lmatools.stream.subset import coroutine, stream, chunk

from lmatools.io.LMAarrayFile import LMAdataFile
from lmatools.coordinateSystems import GeographicSystem, TangentPlaneCartesianSystem, RadarCoordinateSystem


def split_clusters(data,labels):
    no_single = set(labels)
    if -1 in no_single:
        no_single.remove(-1)
    for l in set(no_single):
        idx = (labels==l)
        t = data[idx,-1]
        print(l, t.shape, t.min(), t.max())


@coroutine
def cluster_chunk_pairs(clustered_output_target):
    db = DBSCAN(eps=1.0, min_samples=10, metric='euclidean')
    
    """Receive chunks, and process overlapping pairs"""
    chunk1 = (yield)
    try:
        while True:
            chunk2 = (yield)
            len1 = chunk1.shape[0]
            len2 = chunk2.shape[0]
            print(len1+len2)
            
            # do stuff with chunk 1 and 2
            clusters = db.fit(np.vstack((chunk1, chunk2)))
            labels = clusters.labels_
            
            clustered_output_target.send((chunk1, labels[:len1]))
            
            # pull data out of chunk2 that was clustered as part of chunk 1
            chunk1_labelset = set(labels[:len1])
            if -1 in chunk1_labelset:
                chunk1_labelset.remove(-1) # remove the singleton cluster ID - we want to retain these from chunk 2.
            clustered_in_chunk2 = np.fromiter( ( True if label in chunk1_labelset else False for i,label in enumerate(labels[len1:])) , dtype=bool)
            clustered_output_target.send((chunk2[clustered_in_chunk2], labels[len1:][clustered_in_chunk2]))  
            residuals = chunk2[clustered_in_chunk2==False]
            
            # prepare for another chunk
            if len(residuals) == 0:
                residuals = chunk1[0:0,:] # empty array that preserves the number of dimensions in the data vector - no obs.
            del chunk1
            chunk1 = np.asarray(residuals)
            del residuals
    except GeneratorExit:
        clusters = db.fit(chunk1)
        labels = clusters.labels_
        clustered_output_target.send((chunk1, labels))
        
@coroutine
def chunk_printer():
    total = 0
    try:
        n_last = 0
        while True:
            v = (yield)
            total += v.shape[0]
            n,x = v[:,-1].min(), v[:,-1].max()
            print(v.shape, n,x,x-n_last)
            n_last=n
            del v
    except GeneratorExit:
        print(total)
        
@coroutine
def cluster_printer():
    unique_labels = set([-1.0])
    total = 0
    all_labels = []
    all_v = []
    try:
        n_last = 0
        while True:
            (v, orig_labels) = (yield)
            labels = np.atleast_1d(orig_labels).copy()
            if len(unique_labels) > 0:
                # Only add those labels that represent valid clusters (nonnegative) to the unique set.
                nonsingleton = (labels >= 0)
                labels[nonsingleton] = labels[nonsingleton] + (max(unique_labels) + 1)
                
            for l in set(labels):
                unique_labels.add(l)
            total += v.shape[0]
            #print "---\n{0:s}\n{1:s}\n{2:s}\n{3:s}\n".format(v.shape, set(orig_labels), set(labels), unique_labels)
            #t = v[:,-1]
            #if t.shape[0] > 0:
            #    print t.min(), t.max()
            #for vec, lbl in zip(all_v, all_labels):
            split_clusters(v, labels)

            
            #all_labels.append(labels)
            #all_v.append(v)
            del v, orig_labels, labels
    except GeneratorExit:
        print(total)




#lma=LMAdataFile('/Users/ebruning/Documents/Lightning\ interpretation/Flash-length/Thomas/LYLOUT_120412_01817.exported.dat.gz')
#ctr_lat, ctr_lon, ctr_alt =  40.4463980, -104.6368130, 1000.00

lma=LMAdataFile('/data/20040526/LMA/LYLOUT_040526_224000_0600.dat.gz')
# for line in lma.header:
    # print line

ctr_lat, ctr_lon, ctr_alt =  35.2791257, -97.9178678, 417.90 # OKLMA
#ctr_lat, ctr_lon, ctr_alt =  40.4463980, -104.6368130, 1000.00 # COLMA
good = (lma.stations >= 6) & (lma.chi2 <= 2.0) & (lma.alt < 20e3)
data = lma.data[good]
geoCS = GeographicSystem()
X,Y,Z = geoCS.toECEF(data['lon'], data['lat'], data['alt'])
Xc, Yc, Zc = geoCS.toECEF( ctr_lon, ctr_lat, ctr_alt)
X, Y, Z = X-Xc, Y-Yc, Z-Zc


D_max, t_max = 3.0e3, 0.15 # m, s

X_vector = np.hstack((X[:,None],Y[:,None],Z[:,None])) / D_max
T_vector = data['time'][:,None] / t_max
XYZT = np.hstack((X_vector, T_vector-T_vector.min()))

# Maximum 3 s flash length, normalized to the time separation scale
chunker = chunk(XYZT[:,-1].min(), 3.0/.15, cluster_chunk_pairs(cluster_printer()))
stream(XYZT.astype('float32'),chunker)

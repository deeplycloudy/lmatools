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


import numpy as np


def coroutine(func):
    def start(*args,**kwargs):
        cr = func(*args,**kwargs)
        cr.next()
        return cr
    return start



def split_clusters(data,labels):
    no_single = set(labels)
    if -1 in no_single:
        no_single.remove(-1)
    for l in set(no_single):
        idx = (labels==l)
        t = data[idx,-1]
        print l, t.shape, t.min(), t.max()


def stream(vec, IDs, target):
    for v, vi in zip(vec, IDs):
        target.send((v, vi))
    target.close()
    
def reset_buffer():
    buf = []
    return buf, buf.append
        
@coroutine
def chunk(start_time, max_duration, target, t_idx=-1):
    """ Receive points, assumed in order. Send out chunks"""
    next_time = start_time + max_duration
    v_buffer, append = reset_buffer() # slight optimization since attr lookup is avoided
    i_buffer, append_idx = reset_buffer()
    try:    
        while True:
            v, vi = (yield)
            append(v)
            append_idx(vi)
            t = v[t_idx]
            if t >= next_time:
                target.send((np.asarray(v_buffer), np.asarray(i_buffer)))
                v_buffer, append = reset_buffer()
                i_buffer, append_idx = reset_buffer()
                next_time = t+max_duration
            del v, vi
    except GeneratorExit:
        target.send((np.asarray(v_buffer), np.asarray(i_buffer)))
        target.close()


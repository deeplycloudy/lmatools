""" Fast multiscale clustering

Kushnir et al 2006 (K06)


"""
from __future__ import absolute_import
from __future__ import print_function
import itertools

import numpy as np
from scipy.spatial import cKDTree as KDTree
from scipy import sparse

from .InterpolationMatrix import interpolation_matrix


def initial_weights(data, k, min_dist=None):
    """ Calculate a matrix of weights for the k nearest points.
    
    For N data points in D dimensions, data has shape (N, D).
    
    min_dist (c_dist in K06) is 1.0/(some known minimum distance to a nearest neighbor)
    
    k is the number of nearest neighbors for which weights (K06 eq. 1) should be calculated
    
    """
    
    knn = KDTree(data)
    # D, i contains the point itself and the associated zero distance, of course.
    # D are distances (0, d1, d2, d3, ... , dk-1)
    # i are indices (i0, i1, i2, i3, ..., ik-1) where i[:,0] increment by one
    # Both have shape (N points, k neighbors)

    D, conn = knn.query(knn.data, k)
    if min_dist is None:
        min_dist = D[D>0].min()
    W = np.exp(-min_dist * D)
    return W, conn

def dilute_weights(W, conn, dilution_gamma=0.1):
    """ Using weights and connectivity from initial_weights, return diluted weights
    
        This assumes that we're suppose to dilute (set the weight to zero) only when both
        
        weight of some pair i,j
        -----------divided by ------------
        sum of all weights of all k neighbors to i
        
                     AND
                     
        weight of some pair i,j
        -----------divided by ------------
        sum of all weights of all k neighbors to j

        are less than dilution_gamma 
        
    """
    diluted_weights = W.copy()
    
    # sum weights for the k neareast neighbors of each point
    ik_totals = W.sum(axis=1) 
    
    print("Diluting ...")
    dilution_counter = 0
    for w, k in zip(W, conn):
        # loop over each point, looking at all d distances for k connections
        i = k[0]
        
        for second_index, (ij_weight, j) in enumerate(zip(w[1:], k[1:])):
            ik_dilution_test = ( (ij_weight/ik_totals[i]) < dilution_gamma )
            jk_dilution_test = ( (ij_weight/ik_totals[j]) < dilution_gamma )
            if ik_dilution_test and jk_dilution_test:
                dilution_counter += 1
                diluted_weights[i, second_index+1] = 0.0
    
    print("           ... {0} points".format(dilution_counter))
    return diluted_weights        
    
def sparsify(W,conn):
    """ Construct a sparse matrix representing the k nearest connections 
        between N points given values in W and connectivity conn with shape (N, k). 
        One of the k connections is the trivial case of zero distance to self. """
        
    N = W.shape[0]
    W_sparse = sparse.lil_matrix((N,N)) 
    for row, weights in itertools.izip(conn, W):
        W_sparse[row[0],row[1:]] = weights[1:]
    return W_sparse
    
def fmsc(data, neighbors_to_find=10):
    print("Finding {0} nearest neighbors".format(neighbors_to_find))

    W_knn, conn = initial_weights(data, neighbors_to_find)
    W_dilute = dilute_weights(W_knn, conn)
    print("Making sparse array")
    W = sparsify(W_dilute, conn)
    V = np.arange(W.shape[0], dtype=int) # every point
    G = [(V, W)]
    s = 1
    print("Calculating interpolation matrix")
    P, C_idx, F_idx = interpolation_matrix(W, 0.2)
    
    N = P.shape[0]
    W_s1 = sparse.lil_matrix((N,N))
    W_csc = W.tocsc()
    P_csc = P.tocsc()
    for p_i,q_i in itertools.permutations(list(range(len(C_idx))),2):
        p, q = C_idx[p_i], C_idx[q_i]
        w_pq = 0
        for k in range(N):
            this_P_csc = P_csc[k,p_i]
            if this_P_csc > 0: # subscripting the other two is expensive
                l = list(range(0,k)) + list(range(k+1,N))
                w_pq += (this_P_csc * W_csc[k,l] * P_csc[l,q_i]).sum()
        # print p,q, w_pq
        W_s1[p,q] = w_pq
        
    W_s1_other = P.T * W * P
    
    print(W_s1)
    print(W_s1_other)
    
    
    
        
    return P, W, W_s1, C_idx, F_idx
    
    
    
    
if __name__ == '__main__':
    
    data = np.random.randn(100,2) # (N pts, D-dimensional)
    P, W, W_s1, C_idx, F_idx = fmsc(data, neighbors_to_find=10)
    
    # Coarse aggregate points p and q are connected by
    # W_pq_coarse = sum(P_kp * W_kl * P_lq)
    # for k not equal to l
    
        
        
    
    # All points should be either in the coarse or fine array
    assert len(np.intersect1d(F_idx, C_idx)) == 0
    
    print("Done.")
    
    s = 1
    
    
    
    import matplotlib.pyplot as plt
    ax1=plt.subplot(111)
    ax1.scatter(data[F_idx,0], data[F_idx,1], c='k', marker='s', edgecolors='none')    
    ax1.scatter(data[C_idx,0], data[C_idx,1], c='r', s=36, marker='s', edgecolors='none')
    
    # loop over every point in the fine array, getting corresponding link weights to the coarse points:
    for f_i in F_idx:
        weights = P[f_i,:].data[0]
        coarse_ids = C_idx[P[f_i,:].rows[0]]
        
        for w, c_i in zip(weights, coarse_ids):
            endpoints = [f_i, c_i]
            if w == 1:
                color = (0.0,0.0,1.0)
            else:
                color = (0.8 - 0.8*w, )*3
            ax1.plot(data[endpoints , 0], data[endpoints, 1], '-', color=color)
    
    # loop over new weights, showing the links between the remaining coarse points
    for p,q in itertools.permutations(C_idx,2):
        w = W_s1[p,q]
        coarse_ids = C_idx[P[f_i,:].rows[0]]
        endpoints = [p, q]
        if w > 0:
            color = (0, 1.0, 0)
            foo=ax1.plot(data[endpoints , 0], data[endpoints, 1], '-', color=color)




    plt.show()
    
    # print W
    # print W_dilute
    
    # print "Getting interpolation matrix..."
    # P, C_idx = interpolation_matrix(A, 0.2)
    # print "...done"
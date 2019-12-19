from __future__ import absolute_import
from __future__ import print_function
import itertools

import numpy as np
from scipy import sparse


def interpolation_matrix(A_in, alpha, ignore_negative=True, should_dilute=True, dilutThresh=0.1, dilutLimit=1.0, is_laplac=True):
    """ Create an interpolation matrix according to the given matrix A and a threshold alpha.
        
        For choosing coarse variables, ignore_negative (positive entries only) should be True.
        
        Don't want to use the diagonal in the computation.
        
        Returns
        -------
        Matrix P[number of points, number of coarse points], the likelihood that a data point (first index)
            belongs to the coarse array
        
        C_idx, indices of points on the coarse array
        F_idx, indices of fine points, i.e., those not covered by C_idx
    """

    idx_coarse = None
    P = None


    num_pts = A_in.shape[0]
    # A = A_in.copy()
    
    # throw away the diagonal of A
    # A[np.diag_indices_from(A)] = 0.0
    
    # COO matrices have A.data, A.row, A.col all of the same shape
    A = A_in #.tocoo()
    
    if A.sum() == 0:
        # if A is a diagonal matrix, A - diag(A) == 0
        print("All distances are zero")
        return P, idx_coarse
    
    # Choosing the coarse set variables
    C = np.zeros((1, num_pts))
    F = np.zeros((1, num_pts))
    
    # Initialize with the first node by default
    C[0,0] = 1.0
    
    for i in range(1, num_pts): # we already have manually chosen C[0,0] = 1 as the first seed
        coarse_entries = C[0,:] > 0
        idx_coarse, = np.where(coarse_entries)
        
        # print A[i, idx_coarse].sum(axis=1) / A[i,:].sum(axis=1) # is the result of this line supposed to be stored somewhere?
        
        if is_laplac==True:
            a_row = A.data[i]
            j_indices = A.rows[i]
            if ignore_negative==True:
                idx = [j for j, a_ij in itertools.izip(j_indices, a_row) if a_ij > 0]                    
                # idx = A[i,:] > 0
                # sum_ngbs = (np.abs(A[i, idx])).sum(axis=1)
                # sign_coarse_idx = A[i, idx_coarse] > 0 # find idx of coarse nodes that have pos weight to i
            else:
                idx = [j for j, a_ij in itertools.izip(j_indices, a_row) if a_ij < 0]
                # idx = A[i,:] < 0
                # sum_ngbs = A[i, idx].sum(axis=1)
                # sign_coarse_idx = A[i, idx_coarse] < 0
            # sum_coarse_ngbs = A[i,(idx_coarse & sign_coarse_idx) ].sum(axis=1)   # compute sum of pos weights to coarse nodes    

            coarse_idx = np.intersect1d(idx, idx_coarse) # could also do this with coarse_idx=set(idx) & set(idx_coarse)
            try:
                sum_ngbs = A[i, idx].sum()
            except ValueError:
                # No non-zero weight indices
                sum_ngbs = 0
            try:
                sum_coarse_ngbs = A[i, coarse_idx].sum()
            except ValueError:
                # No non-zero weight indices that are also tagged as coarse
                sum_coarse_ngbs = 0
            
        else:
            # Special test code for Dan's LLE (local linear embedding) purposes
            raise NotImplemented
            # sum_coarse_ngbs = A[i,idx_coarse]   # compute sum of pos weights to coarse nodes    
            # sum_ngbs = np.abs(A[i,:]).sum(axis=1)
            
        if (abs(sum_ngbs) == 0):
            C[0,i] = i
        elif (abs(sum_coarse_ngbs) / abs(sum_ngbs) < alpha):
            C[0,i] = i
        else:
            F[0,i] = i
    
    idx_coarse, = np.where(C[0,:] > 0)
    idx_fine,   = np.where(F[0,:] > 0)
    # P = np.zeros( A[:, :len(idx_coarse)].shape )
    P = sparse.lil_matrix((num_pts,len(idx_coarse)))
    
    for cur in idx_fine:
        A_cur_coarse = A[cur, idx_coarse]
        a_row = A_cur_coarse.data[0]
        j_indices = A_cur_coarse.rows[0]
        
        if is_laplac==True:
            if ignore_negative==True:
                # Find coarse indices that are nonzero positive. Will be intersection of coarse and nonzero cur row entries
                idx = [j for j, a_ij in itertools.izip(j_indices, a_row) if a_ij > 0]
                # idx_coarse_cur = A[cur, idx_coarse] > 0
            else:
                # Find coarse indices that are nonzero negative
                idx = [j for j, a_ij in itertools.izip(j_indices, a_row) if a_ij < 0]
                # idx_coarse_cur = A[cur, idx_coarse] < 0 
                
            # This is unnecessary - only looking at coarse indices of A.
            # sign_filtered_idx = np.intersect1d(idx, idx_coarse)
            
            # We assume here that all indices are either positive or negative, allowing abs at the end.
            P[cur, idx] = (A_cur_coarse[0, idx] / float(abs(A_cur_coarse[0, idx].sum(axis=1)))) #.data[0,:] #use this for P as regular array
        else:
            # Special test code for Dan's LLE purposes
            raise NotImplemented
            # P[cur, :] = A[cur, idx_coarse] / A[cur,idx_coarse].sum(axis=1)
        
        if (dilutThresh > 0) and (dilutThresh < 1):
            a_row = A.data[cur]
            j_indices = A.rows[cur]
            
            fraction = np.floor(len(a_row) * dilutLimit)
            if fraction > 0:
                idx_dilute = [j for j, a_ij in itertools.izip(j_indices, a_row) if a_ij < dilutThresh]
                min_num_to_dilute = min(fraction, len(idx_dilute))
                P[cur, idx_dilute] = 0
                P[cur, :] = P[cur, :] / P[cur,:].sum()
                
            # idx = P[cur,:] > 0
            # fraction = np.floor( idx.shape[0] * dilutLimit )
            # if fraction > 0:
            #     idx_dilute = P[cur, idx] < dilutThresh
            #     min_num_to_dilute = min(fraction, idx_dilute.shape[0])
            #     P[cur, idx[idx_dilute[:min_num_to_dilute]]] = 0
            #     P[cur, :] = P[cur, :] / P[cur,:].sum()
                
    for i in range(len(idx_coarse)):
        P[idx_coarse[i], i] = 1
            
    return P, idx_coarse, idx_fine
    
    
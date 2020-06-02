from __future__ import absolute_import
from __future__ import print_function
import numpy as np

def unique_vectors(*args, **kwargs):
    """ Given D, N-element arrays of vector components return the
        unique vectors as a (D, N_reduced) array. If return_indices_only=True
        is true (the default) then return the N_reduced indices of the original
        arrays corresponding to the unique vectors.
        
        args: x0_i, x1_i, x2_i, ...
            where each x0, x1, etc. are discretized (bin index) values of point
            locations.
        """
    try:
        return_indices_only=kwargs['return_indices_only']
    except KeyError:
        return_indices_only=True
        
    vector_len = len(args)
    itemsize = max((x_i.itemsize for x_i in args))
    inttype = 'i{0:d}'.format(itemsize)
    
    vec_cast = tuple((x_i.astype(inttype) for x_i in args))
        
    # Because we do casting below, ascontiguousarray is important
    locs = np.ascontiguousarray(np.vstack(vec_cast).T)

	# Find unique rows (unique grid cells occupied per flash) by casting to a set
	# of strings comprised of the bytes that make up the rows. Strings are not serialized, just ints cast to byte strings.
	# Based on the example at
	# http://www.mail-archive.com/numpy-discussion@scipy.org/msg04176.html
    # which doesn't quite get it right: it returns unique elements.
    vectorbytes = 'S{0:d}'.format(itemsize*vector_len)
    unq, unq_index = np.unique(locs.view(vectorbytes), return_index=True)

    if return_indices_only==True:
        return unq_index
    else:
        # this is unnecessary if we only need the row indices
        unq_restored = unq.view(inttype).reshape(unq.shape[0],vector_len)
        return unq_restored, unq_index

def extent_density(x, y, ids, x0, y0, dx, dy, xedge, yedge):
    x_i = np.floor( (x-x0)/dx ).astype('int32')
    y_i = np.floor( (y-y0)/dy ).astype('int32')
    
    unq_idx = unique_vectors(x_i, y_i, ids)
    # if x[unq_idx].shape[0] > 1:
    density,edges = np.histogramdd((x[unq_idx],y[unq_idx]), bins=(xedge,yedge))
    return density, edges
    


def test_extent_density():
    # create a set of x,y,id, that have one point in each grid location defined by x0,x1,dx
    x = np.asarray((0.5, 1.5, 2.5, 3.5, 4.5))
    x0, x1 = 0.0, 5.0
    dx = 1.0
    y = x - 2
    y0, y1 = x0 - 2, x1 - 2
    dy = dx
    x_cover, y_cover = np.meshgrid(x,y)
    xedge = np.arange(x0, x1+dx, dx)
    yedge = np.arange(y0, y1+dy, dy)
    ids = np.ones(y_cover.flatten().shape[0], dtype=int)
    
    # ------------------------------
    # replicate the x and y points to have two points in each grid location
    x_doubled = np.hstack((x_cover.flatten(), x_cover.flatten()))
    y_doubled = np.hstack((y_cover.flatten(), y_cover.flatten()))

    # replicate the ids such that the doubled points in each grid belong to the same entity
    ids_doubled = np.hstack((ids, ids))
    density, edges = extent_density(x_doubled, y_doubled, ids_doubled, 
                                    x0, y0, dx, dy, xedge, yedge)
    assert (density == 1).all()

    # replicate the ids such that the doubled points in each grid belong to different entities
    ids_doubled = np.hstack((ids, ids+1))
    density, edges = extent_density(x_doubled, y_doubled, ids_doubled, 
                                    x0, y0, dx, dy, xedge, yedge)
    assert (density == 2).all()
    
    # ------------------------------
    # replicate the x and y points to have two points in each grid location, but leave out
    # one of the points (0.5, -1.5); lower-left in space, upper left in the density array printout
    x_doubled = np.hstack((x_cover.flatten()[1:], x_cover.flatten()))
    y_doubled = np.hstack((y_cover.flatten()[1:], y_cover.flatten()))

    # replicate the ids such that the doubled points in each grid belong to the different entities
    ids_doubled = np.hstack((ids[1:], ids+1))
    density, edges = extent_density(x_doubled, y_doubled, ids_doubled, 
                                    x0, y0, dx, dy, xedge, yedge)
    assert (density == np.array([[ 1.,  2.,  2.,  2.,  2.],
                                [ 2.,  2.,  2.,  2.,  2.],
                                [ 2.,  2.,  2.,  2.,  2.],
                                [ 2.,  2.,  2.,  2.,  2.],
                                [ 2.,  2.,  2.,  2.,  2.]] ) ).all()
    
    # replicate the ids such that the doubled points in each grid belong to the same entity
    ids_doubled = np.hstack((ids[1:], ids))
    density, edges = extent_density(x_doubled, y_doubled, ids_doubled, 
                                    x0, y0, dx, dy, xedge, yedge)
    assert (density == 1).all()


def test_unq():
    locs = np.array([[ 0,  1,  2],
                     [ 0,  0,  0],
                     [ 3,  4,  5],
                     [ 8, 10, 11],
                     [ 3,  5,  5],
                     [ 3,  4,  5],
                     [ 3,  4,  6],
                     [ 9, 10, 11]])


    x_i = locs[:,0]
    y_i = locs[:,1]
    g_id = locs[:,2]


    # vector_len = 3
    # itemsize = max((x_i.itemsize,y_i.itemsize,g_id.itemsize))
    # inttype = 'i{0:d}'.format(itemsize)
    # 
    # 
    # locs = np.ascontiguousarray(np.vstack( (
    #                        x_i.astype(inttype), 
    #                        y_i.astype(inttype), 
    #                        g_id.astype(inttype)
    #                        ) ).T )
    # 
    # vectorbytes = 'S{0:d}'.format(itemsize*vector_len)
    # unq, unq_index = np.unique(locs.view(vectorbytes), return_index=True)
    # 
    # unq_restored = unq.view(inttype).reshape(unq.shape[0],vector_len)
    
    # print unq_index
    # print unq_restored
    # print locs[unq_index,:]
    
    unq_restored, unq_index = unique_vectors(x_i,y_i,g_id, return_indices_only=False)
    assert (unq_index == [1, 0, 2, 6, 4, 3, 7,]).all()
    assert (locs[unq_index,:] == unq_restored).all()
    assert locs[unq_index].shape == (7,3)


def test_unq_func():
    pass
        
if __name__ == '__main__':
    test_unq()
    test_extent_density()
    print("Tests complete.")

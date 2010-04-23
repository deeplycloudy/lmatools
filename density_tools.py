import numpy as np

def unique_vectors(x_i, y_i, g_id):
    """ x_i, y_i: Discretized (bin index) values of point locations.
        id:       Entity ids that group points
        
        Returns a (3, n_points) array of unique points and ids.
        """
        
    # Because we do casting below, ascontiguousarray is important
    locs = np.ascontiguousarray(np.vstack( (
                           x_i.astype('int32'), 
                           y_i.astype('int32'), 
                           g_id.astype('int32')
                       ) ).T )

	# Find unique rows (unique grid cells occupied per flash) by casting to a set
	# of strings comprised of the bytes that make up the rows. Strings are not serialized, just ints cast to byte strings.
	# Based on the example at
	# http://www.mail-archive.com/numpy-discussion@scipy.org/msg04176.html
    # which doesn't quite get it right: it returns unique elements.
    unq, unq_idx = np.unique(locs.view(
                    ('S%d'%(locs.itemsize*locs.shape[1]))*locs.shape[0]), 
                    return_index=True)
    # these are unnecessary since we only need the row indices
    # unq = unq.view('int32')
    # unq = unq.reshape((len(unq)/locs.shape[1], locs.shape[1]))
    
    return unq_idx

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
                      [ 3,  4,  5],
                      [ 8, 10, 11],
                      [ 3,  5,  5],
                      [ 3,  4,  5],
                      [ 3,  4,  6],
                      [ 9, 10, 11]])
     
     # this returns the unique elements
     # unq, unq_idx = np.unique(locs.view('S%d'%locs.itemsize*locs.shape[0]), return_index=True)
     
     unq, unq_idx = np.unique(locs.view(('S%d'%(locs.itemsize*locs.shape[1]))*locs.shape[0]), return_index=True)
     unq = unq.view('int32')
     unq = unq.reshape((len(unq)/locs.shape[1], locs.shape[1]))
     # print unq
     # print locs[unq_idx] == unq
     assert (locs[unq_idx] == unq).all()

        
if __name__ == '__main__':
    test_unq()
    test_extent_density()
    print "Tests complete."

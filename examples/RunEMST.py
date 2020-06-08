'''
RunEMST:
===========================================================================
RunEMST uses AstroPy's Euclidean Minimum Spanning Tree Algorithm to connect 
LMA source points in 3D, allowing for estimation of the flash size (channel length).

The FlashConnection class takes arguments:

   -)positive_sources: [n,3] array for the positive xyz coordinates.
   -)negative_sources: [n,3] array for the negative xyz coordinates.
   -)Nearest Neighbor Max: The maximum number of points from which a single point is branched to, depending on the proximity of the    point.
   -)Edge Cut-Off: (From AstroPy) A fraction (from 0 to 1) of edges to keep for a cluster of source points.
   -)Minimum Cluster Size: Sets the minimum number of points for a single cluster. For small flashes, a smaller cluster size may be desired.
=======
   -)Nearest Neighbor Max: The maximum number of points from which a single point is branched to, 
    depending on the proximity of the point.
   -)Edge Cut-Off: (From AstroPy) A fraction (from 0 to 1) of edges to keep for a cluster of source points.
   -)Minimum Cluster Size: Sets the minimum number of points for a single cluster. 
     For small flashes, a smaller cluster size may be desired.

To run, simply initialize FlashConnection. Functions that produce any desired output are:

    -)flash_size: returns array of channel lengths
    -)plot_emst <--- chose either 2D or 3D visualization.

    -)Running MST_ThreeDim will return the connection segments if those are desired.
    
'''

import numpy as np
from scipy import sparse
from mst_clustering_3D import HierarchicalClustering, get_graph_segments

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class FlashConnection(object):
    
    def __init__(self, positive_sources, 
                       negative_sources,
                       neighbors=10,
                       edgecutoff=0.9,
                       min_cluster=10):
        
        self.positive_sources = positive_sources
        self.negative_sources = negative_sources
        self.channel_length = 0
        self.len_array = []
        
        # channel_selection = ['Positive Channel', 'Negative Channel', 'Total']
        # assert desired_channel in channel_selection, 'Channel Selection Invalid!'
        # self.desired_channel = desired_channel
        
        self.neighbors = neighbors
        self.edgeCutoff = edgecutoff
        self.min_cluster = min_cluster
        
        self.pos_x = 0
        self.pos_y = 0
        self.poz_z = 0
        
        self.neg_x = 0
        self.neg_y = 0
        self.neg_z = 0
        
        
        
        
    def channel_lengths(self, branch):
        '''
        Calculate the channel length of positive, negative leaders, or
        the total channel length using the Euclidean Minimum Spanning Tree 
        connectivity points. 
        
        The (branch) variable is given as (x,y,z) for these connection points.
        '''
        x,y,z = branch
        length = np.sqrt((x[1,:] - x[0,:])**2. + 
                         (y[1,:]-y[0,:])**2. + 
                         (z[1,:]-z[0,:])**2.).sum()
                         
        return length
        
    def flash_size(self):
        '''
        Calculate the flash size of the desired channel found from 
        the Euclidean Minimum Spanning Tree, or for the total channel length.
        
        Channel Lengths are returned in a (3,1) array [len_array] for the individual channel lengths, and 
        the total channel length.
        '''
 
        (self.pos_x, self.pos_y, self.pos_z, 
            self.neg_x, self.neg_y, self.neg_z) = self.MST_ThreeDim()
        
         
        #####SAVE TO ARRAY:
        pos_branch = (self.pos_x, self.pos_y, self.pos_z)
        pos_len = self.channel_lengths(pos_branch)
        
        neg_branch = (self.neg_x, self.neg_y, self.neg_z)
        neg_len = self.channel_lengths(neg_branch)
        
        tot_len = pos_len + neg_len
        
        self.len_array = np.zeros([3,1])
        self.len_array[0,:] = pos_len
        self.len_array[1,:] = neg_len
        self.len_array[2,:] = tot_len
        
        print('Length of positive channel: {0} km'.format(self.len_array[0][0]))
        print('Length of negative channel: {0} km'.format(self.len_array[1][0]))
        print('Total estimated channel length: {0} km'.format(self.len_array[2][0]))
        
        return self.len_array
                
    def MST_ThreeDim(self):
        '''
        This function is similar to the 2-Dimensional Minimum Spanning Tree
        function, however, requires 3-Dimensional data. For LMA data, this is indexed as (for example):

            X_3d_pos = stats.xyzt[stats.pos_charge_mask,0:3]
            X_3d_neg = stats.xyzt[stats.neg_charge_mask,0:3]

        Clustering connects the nearest neighbors for both positive and negative
        charge mask analyzed data. This can be changed to connect all X,Y,Z analyzed charge
        by using the following charge mask:

            X_3d_neg = stats.xyzt[stats.charge_mask,0:3]
        '''

        #--------------------------------
        #Fit data into clustering models:

        nearest_neighbors = self.neighbors
        edge_cutoff = self.edgeCutoff
        cluster_cutoff = self.min_cluster
        
        model_3d = HierarchicalClustering(n_neighbors=nearest_neighbors,
                                          edge_cutoff=edge_cutoff,
                                          min_cluster_size=cluster_cutoff)

        model_3d_neg = HierarchicalClustering(n_neighbors=nearest_neighbors,
                                              edge_cutoff=edge_cutoff,
                                              min_cluster_size=cluster_cutoff)

        model_3d.fit(self.positive_sources)
        model_3d_neg.fit(self.negative_sources)

        n_components = model_3d.n_components_
        labels = model_3d.labels_
        
        nn_components = model_3d_neg.n_components_
        labels_n = model_3d_neg.labels_

        #------------------------------------------
        #Set the x,y,z coordinates from the models:

        self.pos_x, self.pos_y, self.pos_z = get_graph_segments(model_3d.X_train_,
                                                                model_3d.full_tree_)

        self.neg_x, self.neg_y, self.neg_z = get_graph_segments(model_3d_neg.X_train_,
                                                                model_3d_neg.full_tree_)
                                                                
        return (self.pos_x, self.pos_y, self.pos_z, self.neg_x, self.neg_y, self.neg_z)

        
    def emst_compute(self):
        self.MST_ThreeDim()
        
    def plot_emst(self, ax_type):
        
        (self.pos_x, self.pos_y, self.pos_z, 
            self.neg_x, self.neg_y, self.neg_z) = self.MST_ThreeDim()
            
        self.len_array = self.flash_size()
        
        import matplotlib.gridspec as gridspec
        fig = plt.figure(figsize=(10,8))
        
        axis_types = ['3d', '2d', '3D', '2D']
        
        assert ax_type in axis_types, 'Selection must be (3d, 3D, 2D, or 2d)'
        
        if ax_type == '3d' or ax_type=='3D': 
        
            fig.subplots_adjust(hspace=0.17, left=0.1, right=0.95, bottom=0.1, top=0.9)

            ax_3d = fig.add_subplot(111,projection='3d')
            for i in range(self.pos_x.shape[1]):
                ax_3d.plot3D(self.pos_x[:,i],self.pos_y[:,i],self.pos_z[:,i],'-',
                                    zdir=self.pos_z,color='r',linewidth=2,alpha=0.6)
                                
                sources = ax_3d.scatter(self.pos_x[:,i], 
                                        self.pos_y[:,i], 
                                        self.pos_z[:,i], color='red', alpha=0.5)
                                    
            for j in range(self.neg_x.shape[1]):
                ax_3d.plot3D(self.neg_x[:,j],self.neg_y[:,j],self.neg_z[:,j],'-',
                                    zdir=self.neg_z,color='b',linewidth=2,alpha=0.6)
                                                
                nsources = ax_3d.scatter(self.neg_x[:,j], 
                                         self.neg_y[:,j], 
                                         self.neg_z[:,j], color='red',alpha=0.5)
        
            ax_3d.set_xlabel('East Distance (km)')
            ax_3d.set_ylabel('North Distance (km)')
            ax_3d.set_zlabel('Altitude (km)')
            ax_3d.set_title('Lightning Channel Connectivity',fontsize=15)
            
        elif ax_type == '2d' or ax_type =='2D':
            gs = gridspec.GridSpec(3, 3)

            ax_2d= plt.subplot(gs[1:, 0:2]) 
            ax_2db = plt.subplot(gs[0,0:2],sharex=ax_2d)             #<---Cross Section Axis
            ax_2dc = plt.subplot(gs[1:,2],sharey=ax_2d)              #<---Cross Section Axis

            for i in range(self.pos_x.shape[1]):
                #Segments:
                ax_2d.plot(self.pos_x[:,i], self.pos_y[:,i], '-', 
                            color='red', alpha=0.5)
                ax_2db.plot(self.pos_x[:,i], self.pos_z[:,i], '-',
                            color='red', alpha=0.5)
                ax_2dc.plot(self.pos_z[:,i], self.pos_y[:,i], '-',
                            color='red', alpha=0.5)
         
                ax_2d.scatter(self.pos_x[:,i], self.pos_y[:,i], 
                            color='darkred', alpha=0.4)
                ax_2db.scatter(self.pos_x[:,i], self.pos_z[:,i],
                            color='darkred', alpha=0.4)
                ax_2dc.scatter(self.pos_z[:,i], self.pos_y[:,i],
                            color='darkred', alpha=0.4)
                            
                            
            for j in range(self.neg_x.shape[1]):
                #Segments
                ax_2d.plot(self.neg_x[:,j], self.neg_y[:,j], '-', 
                            color='blue', alpha=0.5)
                ax_2db.plot(self.neg_x[:,j], self.neg_z[:,j], '-',
                            color='blue', alpha=0.5)
                ax_2dc.plot(self.neg_z[:,j], self.neg_y[:,j], '-',
                            color='blue', alpha=0.5)
                            
           
                ax_2d.scatter(self.neg_x[:,j], self.neg_y[:,j], 
                            color='darkblue', alpha=0.4)
                ax_2db.scatter(self.neg_x[:,j], self.neg_z[:,j],
                            color='darkblue', alpha=0.4)
                ax_2dc.scatter(self.neg_z[:,j], self.neg_y[:,j],
                            color='darkblue', alpha=0.4)
                            
            ax_2d.set_xlabel('East-West Distance (km)')
            ax_2d.set_ylabel('North-South Distance (km)')
    
            ax_2db.set_ylabel('Altitude (km)')
            ax_2dc.set_xlabel('Altitude (km)')
        
        
            fig.suptitle('Flash Connectivity - Total Length: {0} km'.format(self.len_array[2][0]))
        
        plt.show()
        
        

        

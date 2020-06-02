from __future__ import absolute_import
import matplotlib 
matplotlib.rc('xtick', labelsize=6)
matplotlib.rc('ytick', labelsize=6) 


from numpy import arange

class small_multiples_plot(object):
    
    def __init__(self, fig=None, *args, **kwargs):
        if fig is None:
            raise AssertionError("A valid figure must be passed in.")
            # fig = figure()
        self.fig = fig
        self.fig.subplots_adjust(bottom=0.20, left = 0.1, right=0.9, top=0.9)
        self.colorbar_ax = fig.add_axes((0.1, 0.1, 0.8, 0.05))
        self.multiples = small_multiples(self.fig, **kwargs)
        
    def label_edges(self, bool_val):
        m = self.multiples
        leftside = m[:,0]
        for ax in leftside:
            ax.yaxis.tick_left()
            ax.yaxis.set_visible(bool_val)

        #last row
        bottomedge = m[-1,:]
        for ax in bottomedge:
            ax.xaxis.tick_bottom()
            ax.xaxis.set_visible(bool_val)

        
def small_multiples(f, rows=4, columns=5, margin=(0.0,0.0), zoom_together=True):
    """ Given a figure f, create linked subplots with given number of rows and columns.
        Returns an object array of axes instances [rows, columns], with top left being [0,0].
    """
    # rows = 4    #number in y direction
    # columns = 5 #number in x direction
    
    f.subplots_adjust(wspace=margin[0], hspace=margin[1])

    # should use N.empty((rows,columns),dtype=object)
    # and attribute name should perhaps be changed
    multiples = arange(rows*columns, dtype=object)
    multiples.shape=(rows, columns)
    
    # No axis defined to start with
    commonaxis=None

    for row in range(rows):
        for column in range(columns):
            nth_plot = row*columns + column 
            ax = f.add_subplot(rows, columns, nth_plot + 1, sharex=commonaxis, sharey=commonaxis)
            if not commonaxis and zoom_together:
                commonaxis = ax

            # leaves axes frame, but turns off axis labels and ticks
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            multiples[row, column] = ax
            # ax.plot(range(10), range(10))
            # ax.text(1,1,'%i, %i, %i' % (row, column, nth_plot))
            # print row, column, nth_plot
            
    return multiples



if __name__ == '__main__':
    from pylab import figure, show #, subplot, show
    f=figure()
    m=small_multiples(f)

    #first column
    leftside = m[:,0]
    for ax in leftside:
        ax.yaxis.set_visible(True)
    
    #last row
    bottomedge = m[-1,:]
    for ax in bottomedge:
        ax.xaxis.set_visible(True)
    
    show()

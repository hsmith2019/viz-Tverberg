import numpy as np
from cycler import cycler
from scipy.spatial import ConvexHull
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib.artist import Artist
from matplotlib.mlab import dist_point_to_segment
import matplotlib.pyplot as plt
import matplotlib.path as mpath
from matplotlib import collections, colors, transforms
import matplotlib.animation as animation


class PlotPoints(object):
    def __init__(self,figure,ax,pts):
        self.fig=figure
        self.ax=ax
        self.pts=pts
        self.numPartition=None
        self.center=None
        self.patches=[]
        self.color_order= color_order = np.array([
        [     0,    0.4470,    0.7410],
        [0.8500,    0.3250,    0.0980],
        [0.9290,    0.6940,    0.1250],
        [0.4940,    0.1840,    0.5560],
        [0.4660,    0.6740,    0.1880],
        [0.3010,    0.7450,    0.9330],
        [0.6350,    0.0780,    0.1840],
        ])
        self.legendPatch=[]
        self.anim =None
                
    def Draw(self):
        self.RunTverberg()
        self.fig.canvas.draw()
        plt.show()
    def RunTverberg(self):
        results=self.find_tverberg()
        points = np.array(self.pts)
        self.ax.set_title('Visualizing Tverberg')
        result = next(results)
        partition = result['partition']
        pivot_pt_id = result['pivot_pt_id']
        center = result['center']
        residual = result['residual']
        not_done = True
        iteration = 0;
        self.center = center
        while not_done:
            full_partition = np.append(partition,np.zeros(points.shape[0]-len(partition),'int32'))

            print('iter =',iteration)
            print('\tresidual =',residual)
            print('\tpartition =',full_partition)
            print('\tpivot =',pivot_pt_id)
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_title('Click To Create Points')
            ax.set_xlim((-1, 1))
            ax.set_ylim((-1, 1))

            r = partition.max()+1 # number of partitions
            for i in range(r):
                pts = points[full_partition==i,:]
                if pts.shape[0] > 3:
                    hull = ConvexHull(pts)
                    hull_pts =  pts[hull.vertices[:],:]
                else:
                    hull_pts = pts
                pivot_pt = points[pivot_pt_id,:]

                color = self.color_order[i % self.color_order.shape[0],:]
                pgon_plt = plt.Polygon(
                    hull_pts, 
                    facecolor=color, 
                    edgecolor=np.square(color),
                    alpha=0.2,
                    zorder=1
                    )
                ax.add_patch(pgon_plt)
                pt_plt = plt.scatter(
                    pts[:,0], 
                    pts[:,1], 
                    s=25, 
                    c=[np.square(color)], 
                    zorder=2,
                    edgecolor=np.square(color)
                    )
            pt_plt = plt.scatter(
                center[0], 
                center[1], 
                s=100, 
                c='k',
                zorder=3,
                marker="+"
                )
            if pivot_pt_id >= 0:
                pt_plt = plt.scatter(
                    pivot_pt[0], 
                    pivot_pt[1], 
                    s=150, 
                    zorder=3,
                    facecolors='none', 
                    edgecolors='k'
                    )            
            plt.title('iter = %i, residual = %5.4e' % (iteration, residual) )
            if residual < 1e-16:
                not_done = False
            else:
                plt.show()
                result = next(results)
                partition = result['partition']
                pivot_pt_id = result['pivot_pt_id']
                center = result['center']
                residual = result['residual']
                iteration += 1

    def find_tverberg(self):
        import cvxpy as cp
        import numpy as np
        P=self.pts
        P = np.array(P)

        n,d = P.shape # n number of points, d space dimension
        r = np.ceil(n/(d+1.0)).astype('int') # number of partitions (colors)
        m = (r-1)*(d+1)+1 # number of relevant points

        # reduce to relevant points
        P = P[0:m,:]

        V = np.ones((m,d+1)); V[:,0:d] = P
        W = np.zeros((r,r-1)); W[0:(r-1),:] = np.eye(r-1); W[-1,:] = -1
        VxW = V.reshape( (1,m,1,d+1) ) * W.reshape( (r,1,r-1,1) )

        S = np.transpose( VxW.reshape( (r,m,m-1) ), (1,0,2) ) # S.shape = (r,m-1,d)

        # Now we find colorful Carathéodory set s from S using Barany's Algorithm
        # Barany, I. (1982). A generalization of Carathéodory’s theorem. Discrete 
        # Mathematics, 40(23), 141152. https://doi.org/10.1016/0012-365X(82)90115-7

        I = np.zeros((m,), np.int32); # pt ids of colorful set
        min_sqr_norm = 1
        sqr_norm_tol = 1e-16 # close enough to be considered unchanged
        residual_tol = 1e-16 # close enough to be considered zero
        return_data = dict()
        while min_sqr_norm > residual_tol:
            s = np.zeros( (m-1,m) ) # the colorful set
            for i in range(m):
               s[:,i] = S[i,I[i],:]

            # solve a constrained least squares minimization to find the min-norm point:
            x = cp.Variable(m)
            objective = cp.Minimize(cp.sum_squares(s*x))
            constraints = [0 <= x, x <= 1, sum(x) == 1]
            prob = cp.Problem(objective, constraints)
            min_sqr_norm = prob.solve( eps_abs = 1.0e-10, eps_rel = 1.0e-10, verbose=False)
            coeff = x.value
            min_norm_pt = s.dot(coeff)

            # yield the current guess at the partition, center-point, and residual
            return_data['partition'] = np.array(I)
            return_data['center'] = np.array(P).transpose().dot(coeff)
            return_data['residual'] = min_sqr_norm
            return_data['pivot_pt_id'] = -1
            if min_sqr_norm <= residual_tol:
                break

            # then pivot out a point that does not increase the min_sqr_norm
            pivot_order = coeff.argsort()
            x2 = cp.Variable(m-1)

            for pivot in pivot_order:
                not_pivot = np.setdiff1d(range(m), pivot, assume_unique=True)
                s2 = s[:,not_pivot]
                objective2 = cp.Minimize(cp.sum_squares(s2*x2))
                constraints2 = [0 <= x2, x2 <= 1, sum(x2) == 1]
                prob2 = cp.Problem(objective2, constraints2)
                min_sqr_norm2 = prob2.solve( eps_abs = 1.0e-10, eps_rel = 1.0e-10, verbose=False)
                min_norm_pt2 = s2.dot(x2.value)
                if sum(np.square(min_norm_pt2 - min_norm_pt)) < sqr_norm_tol:
                    prev_pt_id = I[pivot]
                    t = S[pivot,:,:] 
                    not_prev_pt_id = np.setdiff1d(range(r), prev_pt_id, assume_unique=True)
                    pt_id = not_prev_pt_id[t.dot(min_norm_pt)[not_prev_pt_id].argmin()]
                    I[pivot] = pt_id
                    break
            
            # yield return_data
            return_data['pivot_pt_id'] = int(pivot)
            yield return_data
        # A valid Tverberg Partition and Tverberg point have been found
        yield return_data
    
def print_header():
    header_txt = '''
           _        _____                _
    __   _(_)____  |_   _|_   _____ _ __| |__   ___ _ __ __ _
    \ \ / / |_  /____| | \ \ / / _ \ '__| '_ \ / _ \ '__/ _` |
     \ V /| |/ /_____| |  \ V /  __/ |  | |_) |  __/ | | (_| |
      \_/ |_/___|    |_|   \_/ \___|_|  |_.__/ \___|_|  \__, |
                                                        |___/ 
    A python script for visualizing Tverberg's Theorem in the 
    plane. 
    A Tverberg partition of a point set is found via Sarkaria's 
    proof of Tverberg's theorem. Sarkaria's proof is essentially 
    a polynomial time reduction of Tverberg's problem to finding 
    a colorful Caratheodory set in a higher dimensional space.
    The colorful Caratheodory set is found via Barany's pivoting 
    algorithm.
    To demonstrate the progress of the algorithm, we pause before 
    each pivot step and plot the current partition of the given 
    point set, and the current guess at a Tverberg point corres-
    ponding to the current colorful set in the higher dimensonal
    space. Once the colorful set contains the origin in its convex
    closure, the algorithm terminates and a Tverberg partition 
    has been found.
    key: 
        o : next pivot point
        + : most centered point
    Close the current figure to proceed to the next iteration.
    '''
    print(header_txt)
    
if __name__ == '__main__':
    import numpy as np
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('Click To Create Points')
    ax.set_xlim((-1, 1))
    ax.set_ylim((-1, 1))
    num = input ("Enter number Points :")
    n= int(num)
    points = plt.ginput(n, show_clicks= True)
    pts=np.asarray(points)
    PPlot = PlotPoints(fig,ax,pts)
    PPlot.Draw()



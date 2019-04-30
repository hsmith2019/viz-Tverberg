import numpy as np
import scipy
from cycler import cycler
from scipy.spatial import ConvexHull
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Wedge, Polygon
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib.artist import Artist
from matplotlib.mlab import dist_point_to_segment
import matplotlib.pyplot as plt
import matplotlib.path as mpath
from matplotlib import collections, colors, transforms
import random
import math

class PlotPoints(object):
    def __init__(self,center,numpoly,figure,ax,poly=None):
        self.fig=figure
        self.ax=ax
        self.numpoly=numpoly
        self.center=center
        self.poly = []
    
    def Draw(self,points,partitions):
        sub = self.IndixMatToPart(points,partitions)
        self.Cover(sub)
        self.ax.scatter(self.center[0],self.center[1],facecolors='none', edgecolors='r')
        self.fig.canvas.draw()
    def IndixMatToPart(self,points,partitions):
        subsets=[]
        part=0
        while(part < self.numpoly):
            subsets.append([])
            part += 1
        i=0
        j=0
        for p in points:
            while(i < self.numpoly and j<len(points)):
                if(partitions[j][i]):
                    subsets[i].append(points[j])
                    j +=1
                    i= 0
                else:
                    i +=1
        return subsets   
    def DrawPoints(self,Points):
        temp = []
        x = []
        y = []
        dt = np.dtype(float,float)
        temp=np.asarray(Points,dt)
        j=0
        for i in temp:
            x.append((temp[j])[0])
            y.append((temp[j])[1])
            j= j+1
        return (x,y)
    
    def CoverPolyPoints(self,PolyPoints):
        polygon = Polygon(PolyPoints, True)
        patches.append(polygon)
        self.ison = True

    def Cover(self,Sublist):
        for h in Sublist:
            if (len(h) > 3):
                HullList=self.ListHullPoints(h)
                self.CoverPolyPoints(HullList)
            else:
                self.CoverPolyPoints(h)
        p = PatchCollection(patches, alpha=0.4)
        PatchColor=[]
        LegendPatch=[]
        i = 1
        for patch in patches:
            c1=np.random.rand(1,1)[0]
            c2=np.random.rand(1,1)[0]
            c3=np.random.rand(1,1)[0]
            c=(c1[0],c2[0],c3[0])
            PatchColor.append(c)
            pol=Polygon(patch.get_xy(),closed=True,ec=c,lw=1,fill=False, animated=True)
            ax.add_patch(pol)
            self.RunPolyEditor(pol)
            PatchPoints=patch.get_xy()
            x=(self.DrawPoints(PatchPoints))[0]
            y=(self.DrawPoints(PatchPoints))[1]
            PointColor=colors.to_rgba_array(c)
            Label = 'Conv '+ str(i)
            LegendPatch.append(mpatches.Patch(color=c,label =Label))
            self.ax.scatter(x,y,s=15,c=PointColor,marker='o',label='Poly')
            i+=1
        p.set_facecolors(PatchColor)
        self.ax.add_collection(p)
        self.ax.legend(handles=LegendPatch)
        
    def SubConvexHull(self,SubPoints):
        return ConvexHull(SubPoints)
    
    def ListHullPoints(self,Points):
        HullPoints=[]
        hull=self.SubConvexHull(Points)
        for v in hull.vertices:
            HullPoints.append(hull.points[v])
        return HullPoints
    def RunPolyEditor(self,pol): 
        p = PolygonInteractor(self.ax,pol)
    


def find_tverberg(P):
    import cvxpy as cp
    import numpy as np

    print_output = False

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
    sqr_norm_tol = 1e-12 # close enough to be considered unchanged
    min_sqr_norm_tol = 1e-16 # close enough to be considered zero
    residual_history = np.zeros((0,))
    while min_sqr_norm > min_sqr_norm_tol:
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
            if  np.abs( min_sqr_norm - min_sqr_norm2 ) < sqr_norm_tol:
                prev_pt_id = I[pivot]
                t = S[pivot,:,:] 
                not_prev_pt_id = np.setdiff1d(range(r), prev_pt_id, assume_unique=True)
                pt_id = not_prev_pt_id[t.dot(min_norm_pt)[not_prev_pt_id].argmin()]
                I[pivot] = pt_id
                break
        if print_output:
            print("residual =",min_sqr_norm)
        residual_history = np.append(residual_history, min_sqr_norm)

    # create incidence matrix from I
    color_bool = np.zeros((m,r), dtype=bool)
    for i in range(r):
        color_bool[I==i,i] = True;


    # assign the irrelevant points to partitions
    if (n > m):
        color_bool = np.concatenate((color_bool, color_bool[0:n-m,:]), axis= 0)

    # compute the center-point
    center_point = np.array(P).transpose().dot(coeff)

    # output
    if print_output:
        print("incidence_mat =\n",color_bool+0)
        print("center_point =",center_point)
        # import matplotlib.pyplot as plt
        # fig, ax = plt.subplots()
        # y = residual_history
        # ax.semilogy(y, 'o-')
        # ax.set(xlabel='iteration', ylabel='residual', title='residual history')
        # plt.show()
    return (color_bool, center_point);

"""
===========
Poly Editor
===========

This is an example to show how to build cross-GUI applications using
Matplotlib event handling to interact with objects on the canvas.
"""
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.artist import Artist
from matplotlib.mlab import dist_point_to_segment


class PolygonInteractor(object):
    """
    A polygon editor.

    Key-bindings

      't' toggle vertex markers on and off.  When vertex markers are on,
          you can move them, delete them

      'd' delete the vertex under point

      'i' insert a vertex at point.  You must be within epsilon of the
          line connecting two existing vertices

    """

    showverts = True
    epsilon = 5  # max pixel distance to count as a vertex hit

    def __init__(self, ax, poly):
        if poly.figure is None:
            raise RuntimeError('You must first add the polygon to a figure '
                               'or canvas before defining the interactor')
        self.ax = ax
        canvas = poly.figure.canvas
        self.poly = poly

        x, y = zip(*self.poly.xy)
        self.line = Line2D(x, y, marker='o', animated=True)
        self.ax.add_line(self.line)

        self.cid = self.poly.add_callback(self.poly_changed)
        self._ind = None  # the active vert

        canvas.mpl_connect('draw_event', self.draw_callback)
        canvas.mpl_connect('button_press_event', self.button_press_callback)
        canvas.mpl_connect('key_press_event', self.key_press_callback)
        canvas.mpl_connect('button_release_event', self.button_release_callback)
        canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        self.canvas = canvas

    def draw_callback(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.ax.draw_artist(self.poly)
        self.ax.draw_artist(self.line)
        # do not need to blit here, this will fire before the screen is
        # updated

    def poly_changed(self, poly):
        'this method is called whenever the polygon object is called'
        # only copy the artist props to the line (except visibility)
        vis = self.line.get_visible()
        Artist.update_from(self.line, poly)
        self.line.set_visible(vis)  # don't use the poly visibility state

    def get_ind_under_point(self, event):
        'get the index of the vertex under point if within epsilon tolerance'

        # display coords
        xy = np.asarray(self.poly.xy)
        xyt = self.poly.get_transform().transform(xy)
        xt, yt = xyt[:, 0], xyt[:, 1]
        d = np.hypot(xt - event.x, yt - event.y)
        indseq, = np.nonzero(d == d.min())
        ind = indseq[0]

        if d[ind] >= self.epsilon:
            ind = None

        return ind

    def button_press_callback(self, event):
        'whenever a mouse button is pressed'
        if not self.showverts:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        self._ind = self.get_ind_under_point(event)

    def button_release_callback(self, event):
        'whenever a mouse button is released'
        if not self.showverts:
            return
        if event.button != 1:
            return
        self._ind = None

    def key_press_callback(self, event):
        'whenever a key is pressed'
        if not event.inaxes:
            return
        if event.key == 't':
            self.showverts = not self.showverts
            self.line.set_visible(self.showverts)
            if not self.showverts:
                self._ind = None
        elif event.key == 'd':
            ind = self.get_ind_under_point(event)
            if ind is not None:
                self.poly.xy = np.delete(self.poly.xy,
                                         ind, axis=0)
                self.line.set_data(zip(*self.poly.xy))
        elif event.key == 'i':
            xys = self.poly.get_transform().transform(self.poly.xy)
            p = event.x, event.y  # display coords
            for i in range(len(xys) - 1):
                s0 = xys[i]
                s1 = xys[i + 1]
                d = dist_point_to_segment(p, s0, s1)
                if d <= self.epsilon:
                    self.poly.xy = np.insert(
                        self.poly.xy, i+1,
                        [event.xdata, event.ydata],
                        axis=0)
                    self.line.set_data(zip(*self.poly.xy))
                    break
        if self.line.stale:
            self.canvas.draw_idle()

    def motion_notify_callback(self, event):
        'on mouse movement'
        if not self.showverts:
            return
        if self._ind is None:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        x, y = event.xdata, event.ydata

        self.poly.xy[self._ind] = x, y
        if self._ind == 0:
            self.poly.xy[-1] = x, y
        elif self._ind == len(self.poly.xy) - 1:
            self.poly.xy[0] = x, y
        self.line.set_data(zip(*self.poly.xy))

        self.canvas.restore_region(self.background)
        self.ax.draw_artist(self.poly)
        self.ax.draw_artist(self.line)
        self.canvas.blit(self.ax.bbox)

if __name__ == '__main__':
    import numpy as np
    import matplotlib.pyplot as plt
    patches = []
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('Click To Create Points')
    ax.set_xlim((-1, 1))
    ax.set_ylim((-1, 1))
    num = input ("Enter number Points :")
    n= int(num)
    points = plt.ginput(n, show_clicks= True)
    (patitions,center) = find_tverberg(points)
    
    NumPartitions=patitions[0].size
    PPlot = PlotPoints(center,NumPartitions,fig,ax)
    PPlot.Draw(points,patitions)
    fig.canvas.draw()

    plt.show()


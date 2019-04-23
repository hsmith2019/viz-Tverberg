import numpy as np
import sys, select, os
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
from matplotlib.artist import Artist
from matplotlib.mlab import dist_point_to_segment
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
import matplotlib.path as mpath
from matplotlib import collections, colors, transforms
import random
import math


class PlotPoints(object):
    ison=False
    def __init__(self,figure,ax):
        self.ison= False
        self.fig=figure
        self.ax=ax
        
    def Draw(self, Points):
        x=(self.DrawPoints(Points))[0]
        y=(self.DrawPoints(Points))[1]
        self.ax.scatter(x,y,s=1)
        sub = self.Polygons(3,Points)
        self.cover(sub)
        self.fig.canvas.draw()
        
    def Polygons(self,NumPolygons,Points):
        sublist=[]
        length = len(Points)
        num =math.floor(length/NumPolygons)
        i=0
        lent=length
        while(i<length):
            if(i+num < length):
                SubPoints = self.PartionPoints(Points,i,num)
                sublist.append(SubPoints)
            else:
                SubPoints = self.PartionPoints(Points,i,num)
                sublist.append(SubPoints)
            i+=num
        return sublist
        
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
        
    def CoverOn(self):
        if(self.ison):
            return True
        else:
            return False
        
    def PartionPoints(self,Points,firstPoint,numPoints):
        SubSetPoints=[]
        i=firstPoint
        while(i<numPoints + firstPoint):
            SubSetPoints.append(Points[i])
            i +=1
        return SubSetPoints
    
    def cover(self,Sublist):
        for s in Sublist:
            self.CoverPolyPoints(s)
        if(self.CoverOn()== True):
            p = PatchCollection(patches, alpha=0.4)
            color=[]
            for patch in patches:
                c1=np.random.rand(1,1)[0]
                c2=np.random.rand(1,1)[0]
                c3=np.random.rand(1,1)[0]
                c=(c1[0],c2[0],c3[0])
                color.append(c)
                ax.add_patch(Polygon(patch.get_xy(),closed=True,ec=c,lw=1,fill=False))
            p.set_facecolors(color)
            self.ax.add_collection(p)
        
if __name__ == '__main__':
    import numpy as np
    import matplotlib.pyplot as plt
    patches = []
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('Click To Create Points')
    ax.set_xlim((-10, 10))
    ax.set_ylim((-10, 10))
    num = input ("Enter number :")
    n= int(num)
    Points = plt.ginput(n, show_clicks= True)
    PPlot = PlotPoints(fig,ax)
    PPlot.Draw(Points)
    fig.canvas.draw()
    
   

"""Plot colloid system and Voronoi"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib import cm
from matplotlib.colors import ListedColormap, Normalize

class Colloid_Visualisation:


    def __init__(self,prefix):
        """Read in data for visualisation"""

        # Read particle coordinates
        self.particle_crds = np.genfromtxt('{}_particle_crds.dat'.format(prefix)).astype(float)

        # Read ring sizes
        self.ring_sizes = np.genfromtxt('{}_ring_sizes.dat'.format(prefix)).astype(int)

        # Read ring coordinates
        self.ring_crds = []
        with open('{}_ring_crds.dat'.format(prefix),'r') as f:
            for i,line in enumerate(f):
                crds = np.array([float(c) for c in line.split()])
                crds = crds.reshape((self.ring_sizes[i],2))
                self.ring_crds.append(crds)
        self.ring_crds = np.array(self.ring_crds)

        # Initialise figure
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)


    def particles(self):
        """Particle positions"""

        self.ax.scatter(self.particle_crds[:,0],self.particle_crds[:,1],c='k',s=0.5,marker='o',zorder=1)


    def rings(self):
        """Voronoi rings"""

        # Set up colour maps
        map_lower = cm.get_cmap('Blues_r', 128)
        map_upper = cm.get_cmap('Reds', 128)
        colour_mean = 'floralwhite'
        map_lower = ListedColormap(map_lower(np.arange(30,100)))
        map_upper = ListedColormap(map_upper(np.arange(30,100)))
        norm_lower = Normalize(vmin=3,vmax=6)
        norm_upper = Normalize(vmin=6,vmax=12)
        size_colours = []
        for i in range(np.max(self.ring_sizes+1)):
            if i<3:
                size_colours.append('white')
            elif i<6:
                size_colours.append(map_lower(norm_lower(i)))
            elif i==6:
                size_colours.append(colour_mean)
            else:
                size_colours.append(map_upper(norm_upper(i)))

        # Make colour list
        colours = []
        for s in self.ring_sizes:
            colours.append(size_colours[s])

        # Make polygons
        polygons = []
        for ring in self.ring_crds:
            polygons.append(Polygon(ring,True))

        # Add polygons to plot
        self.ax.add_collection(PatchCollection(polygons, facecolors=colours, edgecolor='k', linewidths=0.1, zorder=0))


    def display(self):
        """Show visualisation"""

        plt.show()

if __name__ == '__main__':
    vis = Colloid_Visualisation('./voronoi')
    vis.particles()
    vis.rings()
    vis.display()
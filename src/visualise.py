import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.patches import Circle, Polygon
from matplotlib.collections import PatchCollection
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, Normalize

"""
Visualise hard disk Monte Carlo simulation.
"""


class Visualisation:
    """Visualise coordinate snapshot from binary monte carlo simulation"""


    def __init__(self):
        """Read simulation parameters and coordinate file"""

        self.read_input_file()
        self.read_simulation_files()


    def read_input_file(self):
        """Read visualise setup file"""

        # Check if input file exists in current directory, if not kill process
        if not os.path.isfile('./visualise.inpt'):
            print('Cannot find input file "visualise.inpt" in current directory')
            sys.exit()

        # Read input file and analysis options and parameters
        print('Reading input file')
        with open('visualise.inpt','r') as f:
            f.readline()
            self.prefix = f.readline().split()[0]
            f.readline()
            f.readline()
            self.frame = int(f.readline().split()[0])
            f.readline()
            f.readline()
            self.vis_particles = int(f.readline().split()[0])
            self.vis_vortype = int(f.readline().split()[0])
            self.vis_cellcolour = int(f.readline().split()[0])
            self.vis_save = int(f.readline().split()[0])


    def read_simulation_files(self):
        """Read XYZ and visualisation files from simulation"""

        # Check if simulation files exist in current directory, if not kill process
        if not os.path.isfile('{}.xyz'.format(self.prefix)):
            print('Cannot find simulation file "{}.xyz"'.format(self.prefix))
            sys.exit()
        if not os.path.isfile('{}_vis2d.dat'.format(self.prefix)):
            print('Cannot find simulation file "{}_vis2d.dat"'.format(self.prefix))
            sys.exit()
        if not os.path.isfile('{}_dia.dat'.format(self.prefix)):
            print('Cannot find simulation file "{}_dia.dat"'.format(self.prefix))
            sys.exit()

        # Read coordinate file
        print('Reading simulation xyz file')
        with open('{}.xyz'.format(self.prefix),'r') as f:
            self.n = int(f.readline().split()[0])
            self.crds = np.zeros((self.n,2))
            f.seek(0,0)
            for i in range(self.frame):
                for j in range(2+self.n):
                    f.readline()
            f.readline()
            f.readline()
            for j in range(self.n):
                self.crds[j,:] = np.array([float(c) for c in f.readline().split()[1:3]])

        # Read rings file
        print('Reading simulation ring file')
        with open('{}_vis2d.dat'.format(self.prefix),'r') as f:
            self.rings = []
            if self.vis_vortype != 0:
                while True:
                    frame = int(f.readline().split()[0])
                    vor_type = int(f.readline().split()[0])
                    self.param = float(f.readline().split()[0])
                    self.m = int(f.readline().split()[0])
                    if frame==self.frame and vor_type==self.vis_vortype:
                        for i in range(self.m):
                            ring = np.array([float(c) for c in f.readline().split()])
                            self.rings.append(ring.reshape(ring.shape[0]//2,2))
                        break
                    else:
                        for i in range(self.m):
                            f.readline()

        # Read diameter file
        print('Reading simulation radii and weights file')
        data = np.genfromtxt('{}_dia.dat'.format(self.prefix)).astype(float)
        self.radii = data[:self.n]/2
        self.weights = data[self.n:]


    def visualise(self):
        """Visualise particles and polygons"""

        # Initialise figure
        params = {"figure.figsize": (5, 5)}
        pylab.rcParams.update(params)
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)

        # Add particles if selected
        print(self.crds)
        cmap=cm.get_cmap('coolwarm')
        norm=Normalize(0,20)
        print(np.max(self.radii))
        print(np.max(self.weights))
        if self.vis_particles:
            if self.vis_vortype==-2:
                radii=self.weights
                if self.param>10:
                    self.param=(self.param-10)/2+10
                colour=cmap(norm(self.param))
            else:
                radii=self.radii
                radii=self.weights
                colour='orange'
                colour=(0.8,0.687,0.287,1)
                colour='gold'
            patches = []
            patches_pnts = []
            patches_absent = []
            for i,c in enumerate(self.crds):
                patches.append(Circle(c,radius=radii[i]))
                if radii[i]>0:
                    patches_pnts.append(Circle(c,radius=0.1))
                else:
                    patches_absent.append(Circle(c,radius=0.1))
            self.ax.add_collection(PatchCollection(patches, facecolor=colour, edgecolor='k', alpha=0.5))
            self.ax.add_collection(PatchCollection(patches_pnts, facecolor='k', alpha=1,zorder=1))
            if self.vis_vortype==2:
                self.ax.add_collection(PatchCollection(patches_absent, facecolor='k', alpha=0.5,zorder=1))
            else:
                self.ax.add_collection(PatchCollection(patches_absent, facecolor='k', alpha=1,zorder=1))

        # Add voronoi
        if self.vis_vortype!=0:
            patches = []
            colours = []
            if self.vis_cellcolour==1:
                cell_colours = self.init_cell_colours()
            else:
                cell_colours = [(0,0,0,0)]*100
            for i in range(self.m):
                patches.append(Polygon(self.rings[i],True))
                colours.append(cell_colours[self.rings[i][:,0].size])
            self.ax.add_collection(PatchCollection(patches, facecolor=colours, edgecolor='k', linewidth=1, zorder=0))

        # Sandbox
        # print(np.max(self.radii))
        # cmap=cm.get_cmap('coolwarm')
        # norm=Normalize(0,np.max(20))
        sandbox=False
        if sandbox:
        #     z=16
        #     w=np.zeros_like(self.radii)
        #     mask=2*self.radii>z
        #     w[mask]=z**0.5*np.sqrt(2*self.radii[mask]-z)
        #     patches = []
        #     for i,c in enumerate(self.crds):
        #         patches.append(Circle(c,radius=w[i]))
        #     self.ax.add_collection(PatchCollection(patches, facecolor=cmap(norm(z)), edgecolor='k'))
            with open('./phi.dat','w') as f:
                for z in np.arange(0,np.max(self.radii)*2+0.5,0.01):
                    w=np.zeros_like(self.radii)
                    mask=2*self.radii>z
                    w[mask]=z**0.5*np.sqrt(2*self.radii[mask]-z)
                    phi=np.sum(np.pi*w**2)/52359.9
                    # phi=np.sum(np.pi*w**2)/1309
                    f.write('{:.6f} {:.6f}\n'.format(z,phi))



        # Set axes
        buffer = 1.6
        lim = buffer*np.max(np.abs(self.crds))
        self.ax.set_xlim((-lim,lim))
        self.ax.set_ylim((-lim,lim))
        self.ax.set_axis_off()

        # Show figure
        if self.vis_save:
            plt.savefig('{}_{}_{}.png'.format(self.prefix,self.frame,self.vis_vortype),dpi=400)
        plt.show()


    def init_cell_colours(self):

        av_ring_size=6
        map_lower = cm.get_cmap('Blues_r', 128)
        map_upper = cm.get_cmap('Reds', 128)
        map_mean=cm.get_cmap("Greys")
        map_lower=ListedColormap(map_lower(np.arange(20,100)))
        map_upper=ListedColormap(map_upper(np.arange(20,100)))
        map_lower=ListedColormap(map_lower(np.arange(0,60)))
        map_upper=ListedColormap(map_upper(np.arange(20,120)))

        norm_lower=Normalize(vmin=av_ring_size-3,vmax=av_ring_size)
        norm_upper=Normalize(vmin=av_ring_size,vmax=av_ring_size+6)
        colour_mean=map_mean(50)
        ring_colours=[]
        for i in range(30):
            if i < 3:
                ring_colours.append("white")
            elif np.abs(i-av_ring_size)<1e-6:
                ring_colours.append(colour_mean)
            elif i<av_ring_size:
                ring_colours.append(map_lower(norm_lower(i)))
            else:
                ring_colours.append(map_upper(norm_upper(i)))

        print(ring_colours[4])
        return ring_colours


if __name__ == '__main__':
    vis = Visualisation()
    vis.visualise()

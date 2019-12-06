import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
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


    def read_simulation_files(self):
        """Read XYZ and visualisation files from simulation"""

        # Check if simulation files exist in current directory, if not kill process
        if not os.path.isfile('{}.xyz'.format(self.prefix)):
            print('Cannot find simulation file "{}.xyz"'.format(self.prefix))
            sys.exit()
        if not os.path.isfile('{}_vis3d.dat'.format(self.prefix)):
            print('Cannot find simulation file "{}_vis3d.dat"'.format(self.prefix))
            sys.exit()
        if not os.path.isfile('{}_dia.dat'.format(self.prefix)):
            print('Cannot find simulation file "{}_dia.dat"'.format(self.prefix))
            sys.exit()

        # Read coordinate file
        print('Reading simulation xyz file')
        with open('{}.xyz'.format(self.prefix),'r') as f:
            self.n = int(f.readline().split()[0])
            self.crds = np.zeros((self.n,3))
            f.seek(0,0)
            for i in range(self.frame):
                for j in range(2+self.n):
                    f.readline()
            f.readline()
            f.readline()
            for j in range(self.n):
                self.crds[j,:] = np.array([float(c) for c in f.readline().split()[1:]])

        # Read rings file
        print('Reading simulation ring file')
        with open('{}_vis3d.dat'.format(self.prefix),'r') as f:
            self.rings = []
            if self.vis_vortype != 0:
                while True:
                    frame = int(f.readline().split()[0])
                    vor_type = int(f.readline().split()[0])
                    self.z_limits = [float(z) for z in f.readline().split()[:2]]
                    self.m = int(f.readline().split()[0])
                    if frame==self.frame and vor_type==self.vis_vortype:
                        for i in range(self.m):
                            ring = np.array([float(c) for c in f.readline().split()])
                            self.rings.append(ring.reshape(ring.shape[0]//3,3))
                        break
                    else:
                        for i in range(self.m):
                            f.readline()

        # Read diameter file
        print('Reading simulation radii file')
        self.radii = np.genfromtxt('{}_dia.dat'.format(self.prefix)).astype(float)/2.0


    def visualise(self):
        """Visualise particles and polygons"""

        # Initialise figure
        params = {"figure.figsize": (5, 5)}
        pylab.rcParams.update(params)
        fig = plt.figure()
        self.ax = Axes3D(fig)
        self.ax.set_axis_off()

        # Add particles if selected
        # if self.vis_particles:
        #     self.ax.scatter(self.crds[:,0],self.crds[:,1],self.crds[:,2], color="r")

        # Add voronoi
        if self.vis_vortype!=0:
            for i in range(self.m):
                if np.all(np.abs(self.rings[i][:,2]-self.z_limits[0])<1e-6):
                    colour = 'r'
                elif np.all(np.abs(self.rings[i][:,2]-self.z_limits[1])<1e-6):
                    colour = 'b'
                else:
                    colour = 'w'

                # self.ax.add_collection3d(Poly3DCollection([self.rings[i]],facecolor=(0,0,0,0),edgecolor="k"))
                self.ax.add_collection3d(Poly3DCollection([self.rings[i]],facecolor=colour,edgecolor="k"))

    #     # Add polygons if selected
    #     if self.vis_polygons:
    #         if self.voronoi is not None:
    #             patches_p = []
    #             for i in range(self.n):
    #                 p = self.voronoi.power_polygons[i]
    #                 if p.size > 0:
    #                     patches_p.append(Polygon(self.voronoi.power_crds[p], True))
    #             self.ax.add_collection(PatchCollection(patches_p, facecolor=(0, 0, 0, 0), edgecolor='grey'))
    #         if self.power is not None:
    #             patches_p = []
    #             for i in range(self.n):
    #                 p = self.power.power_polygons[i]
    #                 if p.size > 0:
    #                     patches_p.append(Polygon(self.power.power_crds[p], True))
    #             self.ax.add_collection(PatchCollection(patches_p, facecolor=(0, 0, 0, 0), edgecolor='k'))
    #
    #     # Add size labels if selected
    #     if self.vis_sizelabel:
    #         if self.voronoi is not None and self.power is None:
    #             for i in range(self.n):
    #                 p = self.voronoi.power_polygons[i]
    #                 if p.size>0:
    #                     label_x = self.voronoi.crds[i,0]
    #                     label_y = self.voronoi.crds[i,1]
    #                     label=p.size
    #                     self.ax.text(label_x,label_y,label,size=5,ha='center',va='center')
    #         elif self.power is not None and self.voronoi is None:
    #             for i in range(self.n):
    #                 p = self.power.power_polygons[i]
    #                 if p.size>0:
    #                     label_x = self.power.crds[i,0]
    #                     label_y = self.power.crds[i,1]
    #                     label=p.size
    #                     self.ax.text(label_x,label_y,label,size=5,ha='center',va='center')
    #         elif self.voronoi is not None and self.power is not None:
    #             for i in range(self.n):
    #                 p = self.power.power_polygons[i]
    #                 label_x = self.power.crds[i,0]
    #                 label_y = self.power.crds[i,1]
    #                 label=p.size-self.voronoi.power_polygons[i].size
    #                 if label<0:
    #                     self.ax.text(label_x,label_y,label,size=5,ha='center',va='center')
    #                 elif label>0:
    #                     self.ax.text(label_x,label_y,'+{}'.format(label),size=5,ha='center',va='center')
    #
    #
        # Set axes
        buffer = 1.5
        lim = buffer*np.max(np.abs(self.crds))
        self.ax.set_xlim((-lim,lim))
        self.ax.set_ylim((-lim,lim))
        self.ax.set_zlim((-lim,lim))
        self.ax.set_axis_off()

        # Show figure
        plt.savefig('vis.png',dpi=400)
        plt.show()

if __name__ == '__main__':
    vis = Visualisation()
    vis.visualise()

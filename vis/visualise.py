""""Visualise binary non-additive hard disc monte carlo simulation"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.patches import Circle, Polygon
from matplotlib.collections import PatchCollection
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, Normalize
from binary_power import Periodic_Binary_Power


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
            self.vis_type = int(f.readline().split()[0])
            self.vis_particles = int(f.readline().split()[0])
            self.vis_delaunay = int(f.readline().split()[0])
            self.vis_polygons = int(f.readline().split()[0])
            self.vis_sizelabel = int(f.readline().split()[0])
            self.vis_polycolour = int(f.readline().split()[0])


    def read_simulation_files(self):
        """Read auxilary and coordinate files from simulation"""

        # Check if simulation files exist in current directory, if not kill process
        if not os.path.isfile('{}.aux'.format(self.prefix)):
            print('Cannot find simulation file "{}.aux"'.format(self.prefix))
            sys.exit()
        if not os.path.isfile('{}.xyz'.format(self.prefix)):
            print('Cannot find simulation file "{}.xyz"'.format(self.prefix))
            sys.exit()

        # Read aux file and get simulation parameters
        print('Reading simulation auxilary file')
        with open('{}.aux'.format(self.prefix),'r') as f:
            self.additive = int(f.readline().split()[0])
            self.n_a = int(f.readline().split()[0])
            self.n_b = int(f.readline().split()[0])
            self.r_a = float(f.readline().split()[0])
            self.r_b = float(f.readline().split()[0])
            f.readline()
            f.readline()
            self.cell_length = float(f.readline().split()[0])
            # f.readline()
            # f.readline()
            # frames_total = int(float(f.readline().split()[0]))
        self.n = self.n_a+self.n_b
        self.w_a = 2.0*self.r_a*np.sqrt(self.r_a*self.r_b)/(self.r_a+self.r_b)
        self.w_b = 2.0*np.sqrt(self.r_a*self.r_b)*(1.0-self.r_a/(self.r_a+self.r_b))
        # if self.frame>=frames_total:
        #     print('Frame outside simulation range ({})'.format(frames_total))
        #     sys.exit()

        # Read coordinate file
        print('Reading simulation xyz file')
        self.crds_a = np.zeros((self.n_a,2))
        self.crds_b = np.zeros((self.n_b,2))
        with open('{}.xyz'.format(self.prefix),'r') as f:
            for i in range(self.frame-1):
                for j in range(2+self.n):
                    f.readline()
            f.readline()
            f.readline()
            for j in range(self.n_a):
                self.crds_a[j,:] = np.array([float(c) for c in f.readline().split()[1:3]])
            for j in range(self.n_b):
                self.crds_b[j,:] = np.array([float(c) for c in f.readline().split()[1:3]])


    def visualise(self):
        """Visualise particles and polygons"""

        # Initialise figure
        params = {"figure.figsize": (5, 5)}
        pylab.rcParams.update(params)
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)

        # Calculate partitions
        self.voronoi = None
        self.power = None
        if self.vis_type==0 or self.vis_type==2:
            self.voronoi = Periodic_Binary_Power(self.crds_a,self.crds_b,0.0,0.0,self.cell_length)
            self.voronoi.delaunay()
            self.voronoi.power()
        if self.vis_type==1 or self.vis_type==2:
            self.power = Periodic_Binary_Power(self.crds_a,self.crds_b,self.w_a,self.w_b,self.cell_length)
            self.power.delaunay()
            self.power.power()

        # Add particles if selected
        if self.vis_particles:
            patches_a = []
            patches_b = []
            if self.voronoi is not None and self.power is None:
                r_a = 1.0
                r_b = 1.0
                for c in self.voronoi.crds_a:
                    patches_a.append(Circle(c, radius=r_a))
                for c in self.voronoi.crds_b:
                    patches_b.append(Circle(c, radius=r_b))
            if self.power is not None:
                r_a = self.w_a
                r_b = self.w_b
                for c in self.power.crds_a:
                    patches_a.append(Circle(c, radius=r_a))
                for c in self.power.crds_b:
                    patches_b.append(Circle(c, radius=r_b))
            self.ax.add_collection(PatchCollection(patches_a, facecolor='blue', alpha=0.5))
            self.ax.add_collection(PatchCollection(patches_b, facecolor='red', alpha=0.5))
        
        # Add delaunay if selected
        if self.vis_delaunay:
            if self.voronoi is not None and self.power is None:
                self.ax.scatter(self.voronoi.crds[:,0],self.voronoi.crds[:,1],marker='o',s=20,color='k',zorder=2)
                for i,nl in enumerate(self.voronoi.neighbour_list):
                    for j in nl:
                        d = np.sqrt((self.voronoi.crds[i,0]-self.voronoi.crds[j,0])**2+(self.voronoi.crds[i,1]-self.voronoi.crds[j,1])**2)
                        if i<j and d<self.cell_length*0.5:
                            self.ax.plot([self.voronoi.crds[i,0],self.voronoi.crds[j,0]],[self.voronoi.crds[i,1],self.voronoi.crds[j,1]],lw=1,ls='--',color='k',zorder=2)    
            if self.power is not None:
                self.ax.scatter(self.power.crds[:,0],self.power.crds[:,1],marker='o',s=20,color='k',zorder=2)
                for i,nl in enumerate(self.power.neighbour_list):
                    for j in nl:
                        d = np.sqrt((self.power.crds[i,0]-self.power.crds[j,0])**2+(self.power.crds[i,1]-self.power.crds[j,1])**2)
                        if i<j and d<self.cell_length*0.5:
                            self.ax.plot([self.power.crds[i,0],self.power.crds[j,0]],[self.power.crds[i,1],self.power.crds[j,1]],lw=1,ls='--',color='k',zorder=2)

        # Add polygons if selected
        if self.vis_polygons:
            if self.voronoi is not None:
                polygon_colours = self.generate_polygon_colours(self.voronoi.power_polygons)
                patches_p = []
                for i in range(self.n):
                    p = self.voronoi.power_polygons[i]
                    if p.size > 0:
                        patches_p.append(Polygon(self.voronoi.power_crds[p], True))
                self.ax.add_collection(PatchCollection(patches_p, facecolor=polygon_colours, edgecolor='k'))
            if self.power is not None:
                polygon_colours = self.generate_polygon_colours(self.power.power_polygons)
                patches_p = []
                for p in self.power.power_polygons:
                    if p.size > 0:
                        patches_p.append(Polygon(self.power.power_crds[p], True))
                self.ax.add_collection(PatchCollection(patches_p, facecolor=polygon_colours, edgecolor='k'))

        # Add size labels if selected
        if self.vis_sizelabel:
            if self.voronoi is not None and self.power is None:
                for i in range(self.n):
                    p = self.voronoi.power_polygons[i]
                    if p.size>0:
                        label_x = self.voronoi.crds[i,0]
                        label_y = self.voronoi.crds[i,1]
                        label=p.size
                        self.ax.text(label_x,label_y,label,size=5,ha='center',va='center')
            elif self.power is not None and self.voronoi is None:
                for i in range(self.n):
                    p = self.power.power_polygons[i]
                    if p.size>0:
                        label_x = self.power.crds[i,0]
                        label_y = self.power.crds[i,1]
                        label=p.size
                        self.ax.text(label_x,label_y,label,size=5,ha='center',va='center')
            elif self.voronoi is not None and self.power is not None:
                for i in range(self.n):
                    p = self.power.power_polygons[i]
                    label_x = self.power.crds[i,0]
                    label_y = self.power.crds[i,1]
                    label=p.size-self.voronoi.power_polygons[i].size
                    if label<0:
                        self.ax.text(label_x,label_y,label,size=5,ha='center',va='center')
                    elif label>0:
                        self.ax.text(label_x,label_y,'+{}'.format(label),size=5,ha='center',va='center')

        #### Testing Section ####
        # self.convert_to_polydisperse()
        # self.polydisperse_power = Periodic_Binary_Power(self.crds_p,np.zeros((0,2)),self.r_p,0.0,self.cell_length)
        # self.polydisperse_power.delaunay()
        # self.polydisperse_power.power()
        # patches_p = []
        # for i,c in enumerate(self.crds_p):
        #     patches_p.append(Circle(c, radius=self.r_p[i]))
        # self.ax.add_collection(PatchCollection(patches_p, facecolor='blue', alpha=0.5))
        # polygon_colours = self.generate_polygon_colours(self.polydisperse_power.power_polygons)
        # patches_p = []
        # for p in self.polydisperse_power.power_polygons:
        #     if p.size > 0:
        #         patches_p.append(Polygon(self.polydisperse_power.power_crds[p], True))
        # self.ax.add_collection(PatchCollection(patches_p, facecolor=polygon_colours, edgecolor='k'))
        # self.convert_to_monodisperse()
        self.power.reverse()
        ####


        # Set axes
        buffer = 0.6
        lim = buffer*self.cell_length
        self.ax.set_xlim((-lim,lim))
        self.ax.set_ylim((-lim,lim))
        self.ax.set_axis_off()
        
        # Show figure
        plt.savefig('vis.png',dpi=400)
        plt.show()


    def convert_to_monodisperse(self):

        pass




    def convert_to_polydisperse(self):

        self.r_p = np.zeros(self.n)
        self.crds_p = np.zeros((self.n,2))
        for i in range(self.n):
            self.crds_p[i,:] = self.power.crds[i,:]
            poly_crds = self.power.power_crds[self.power.power_polygons[i]]
            dist = 0.0
            v = np.zeros(2)
            while True:
                d,v,e = self.min_dist_to_edge(self.crds_p[i]+0.01*v,poly_crds)
                if d<dist:
                    break
                else:
                    self.crds_p[i] += 0.01*v
                    self.r_p[i] = d
                    dist = d


    def min_dist_to_edge(self,p,poly):

        n = int(poly.size/2)
        a = np.zeros((n,2))
        u = np.zeros((n,2))
        v = np.zeros((n,2))
        d = np.inf
        e = -1

        a[:,:] = poly[:,:]
        u = np.array([(poly[(i+1)%n]-poly[i]) for i in range(n)])
        for i in range(n):
            aa = a[i]
            uu = u[i]/np.sqrt(np.sum(u[i]**2))
            ap = p-aa
            vv = ap - np.sum(ap*uu)*uu
            dd = np.sqrt(np.sum(vv*vv))
            if dd<d:
                d = dd
                v = vv
                e = i

        return d,v/np.sqrt(np.sum(v**2)),e



    def generate_polygon_colours(self,polygons):
        """Generate colours for polygons - either transparent or coloured by number of edges"""

        # Initialise container
        polygon_colours = []

        # Transparent polygons
        if self.vis_polycolour == 0:
            for p in polygons:
                if p.size > 0:
                    polygon_colours.append((0,0,0,0))

        # Coloured by number of edges
        if self.vis_polycolour == 1:
            # Set up colours, blue<6, grey=6, red>6
            map_lower = cm.get_cmap('Blues_r', 128)
            map_upper = cm.get_cmap('Reds', 128)
            map_mean=cm.get_cmap("Greys")
            map_lower=ListedColormap(map_lower(np.arange(30,100)))
            map_upper=ListedColormap(map_upper(np.arange(30,100)))
            norm_lower=Normalize(vmin=3,vmax=6)
            norm_upper=Normalize(vmin=6,vmax=12)
            colour_mean=map_mean(50)
            # Make list of colours based on size
            for p in polygons:
                if p.size < 3:
                    polygon_colours.append((0,0,0,0))
                elif p.size == 6:
                    polygon_colours.append(colour_mean)
                elif p.size < 6:
                    polygon_colours.append(map_lower(norm_lower(p.size)))
                else:
                    polygon_colours.append(map_upper(norm_upper(p.size)))

        return polygon_colours


if __name__ == '__main__':
    vis = Visualisation()
    vis.visualise()

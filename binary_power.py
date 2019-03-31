import os
import numpy as np
import sys
from numpy import linalg
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pylab as pylab
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


class Binary_Power:
    """Base class for power diagram for binary hard-disc system"""


    def __init__(self):
        """Base class should not be instantiated"""

        raise TypeError('Cannot instantiate base Binary_Power class')


    def delaunay(self):
        """Calculates Delaunay triangulation - all information required for network analysis"""

        # Make polyhedron coordinates for convex hull
        polyhedron_crds = self.make_polyhedron_coordinates()

        # Calculate convex hull triangulation
        self.lower_hull_triangulation(polyhedron_crds)

        # Calculate unordered neighbour list
        self.generate_neighbour_list()

        # Calculate number of connections to each node
        self.calculate_coordination()


    def lower_hull_triangulation(self,crds):
        """Calculate convex hull and triangulate lower hull"""

        # Calculate convex hull
        convex_hull = ConvexHull(crds)

        # Triangulate lower hull with triangles oriented anticlockwise
        hull_triangulation = []
        for i, eqn in enumerate(convex_hull.equations):
            if eqn[2] <= 0.0:
                s = convex_hull.simplices[i]
                mat = np.concatenate([np.stack([crds[s[0],:2], crds[s[1],:2], crds[s[2],:2]]), np.ones((3, 1))], axis=1)
                if linalg.det(mat) > 0:
                    hull_triangulation.append(s)
                else:
                    hull_triangulation.append(np.array([s[0],s[2],s[1]]))
        self.hull_triangulation = np.array(hull_triangulation)


    def calculate_coordination(self):
        """Calculate number of connections to each node in the Delaunay"""

        self.cnd = np.zeros(self.n,dtype=int)
        for i,nl in enumerate(self.neighbour_list):
            self.cnd[i] = nl.size
        self.max_cnd = np.max(self.cnd)


    def delaunay_node_distribution(self,k_lim=None):
        """Calculate Delaunay node distribution for each component and combined"""

        # Calculate k-range - use supplied limit if necessary (easier for analysis)
        warning = False
        if k_lim is None:
            k_lim = self.max_cnd
        elif self.max_cnd>k_lim:
            warning = True
        k = np.arange(k_lim+1,dtype=int)

        # Calculate partial and total distributions
        p_ka = np.zeros_like(k,dtype=float)
        p_kb = np.zeros_like(k,dtype=float)
        cnds_a = self.cnd[:self.n_a]
        cnds_b = self.cnd[self.n_a:]
        for i,j in enumerate(k):
            p_ka[i] = (cnds_a==j).sum()
            p_kb[i] = (cnds_b==j).sum()
        p_kt = p_ka+p_kb
        p_ka /= p_ka.sum()
        p_kb /= p_kb.sum()
        p_kt /= p_kt.sum()

        # Pack into dictionary
        node_dist = {'k':k,'p_ka':p_ka,'p_kb':p_kb,'p_kt':p_kt}

        return node_dist,warning


    # def delaunay_node_distributions(self,k_lim=None):
    #     """Calculate Delaunay node distribution, edge distribution and joint-edge distribution"""
    #
    #
    #
    #     # Calculate edge distributions
    #     e_jk = np.zeros((k_lim+1,k_lim+1))
    #     for i,nl in enumerate(self.neighbour_list):
    #         k_i = num_cnxs[i]
    #         for j in nl:
    #             k_j = num_cnxs[j]
    #             e_jk[k_i,k_j] += 1
    #     e_jk /= e_jk.sum()
    #     q_k = np.sum(e_jk,axis=1)
    #
    #     return k,p_k,q_k,e_jk,warning


class Periodic_Binary_Power(Binary_Power):
    """Periodic power diagram for binary system"""


    def __init__(self,crds_a,crds_b,r_a,r_b,cell_length,**kwargs):
        """Initialise with particle coordinates, radii and periodicity information"""

        # Get input variables
        self.crds_a = crds_a
        self.crds_b = crds_b
        self.r_a = r_a
        self.r_b = r_b
        self.cell_length = cell_length
        self.n_a = self.crds_a[:,0].size
        self.n_b = self.crds_b[:,0].size
        self.n = self.n_a + self.n_b

        # Combine coordinates and radii
        self.crds = np.zeros((self.n,2))
        self.crds[:self.n_a] = self.crds_a
        self.crds[self.n_a:] = self.crds_b
        self.r = np.zeros(self.n)
        self.r[:self.n_a] = r_a
        self.r[self.n_a:] = r_b

        # Subtract centre of mass
        centre_of_mass = np.average(self.crds,axis=0)
        self.crds -= centre_of_mass
        self.crds_a -= centre_of_mass
        self.crds_b -= centre_of_mass


    def make_polyhedron_coordinates(self):
        """Make 9 periodic images of coordinates with central image as original, raise coordinates to form polyhedron."""

        # Make periodic images
        polyhedron_crds = np.zeros((9*self.n,3))
        polyhedron_crds[:self.n,:2] = self.crds[:]
        image = 1
        for i in [-1,0,1]:
            translation = np.zeros_like(self.crds)
            translation[:,0] = i*self.cell_length
            for j in [-1,0,1]:
                translation[:,1] = j*self.cell_length
                if i==0 and j==0:
                    pass
                else:
                    polyhedron_crds[image*self.n:(image+1)*self.n,:2] = self.crds+translation
                    image += 1
        # Raise
        for i in range(9):
            n = i*self.n
            m = (i+1)*self.n
            polyhedron_crds[n:m,2] = np.sum(polyhedron_crds[n:m,:2]**2,axis=1)-self.r

        return polyhedron_crds


    def generate_neighbour_list(self):
        """Generate list of adjacent particles"""

        # List will be unordered i.e. neighbours cannot be traced to give polygon
        neighbour_list = [[] for i in range(self.n)]

        # Make non-unique list
        for tri in self.hull_triangulation:
            a,b,c = tri
            aa,bb,cc = tri%self.n
            if a<self.n:
                neighbour_list[a].append(bb)
                neighbour_list[a].append(cc)
            if b<self.n:
                neighbour_list[b].append(aa)
                neighbour_list[b].append(cc)
            if c<self.n:
                neighbour_list[c].append(aa)
                neighbour_list[c].append(bb)

        # Make unique list
        self.neighbour_list = []
        for nl in neighbour_list:
            self.neighbour_list.append(np.array(np.unique(nl),dtype=int))

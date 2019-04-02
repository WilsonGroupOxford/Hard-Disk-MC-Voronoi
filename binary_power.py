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


    def delaunay_node_distributions(self,k_lim=None):
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
        p_kc = p_ka+p_kb
        p_ka /= p_ka.sum()
        p_kb /= p_kb.sum()
        p_kc /= p_kc.sum()

        # Pack into dictionary
        node_dist = {'k':k,'p_ka':p_ka,'p_kb':p_kb,'p_kc':p_kc}

        return node_dist,warning


    def delaunay_edge_distributions(self,k_lim=None):
        """Calculate Delaunay edge distribution and joint-edge distribution"""

        # Calculate k-range - use supplied limit if necessary (easier for analysis)
        warning = False
        if k_lim is None:
            k_lim = self.max_cnd
        elif self.max_cnd>k_lim:
            warning = True
        k = np.arange(k_lim+1,dtype=int)

        # Calculate edge distributions
        e_jk = np.zeros((k_lim+1,k_lim+1))
        for i,nl in enumerate(self.neighbour_list):
            k_i = self.cnd[i]
            for j in nl:
                k_j = self.cnd[j]
                e_jk[k_i,k_j] += 1
        e_jk /= e_jk.sum()
        q_k = np.sum(e_jk,axis=1)

        # Pack into dictionary
        edge_dist = {'k':k,'q_k':q_k,'e_jk':e_jk}

        return edge_dist,warning


    def power(self):
        """Generates power diagram - more expensive and only for visualisation"""

        # Regenerate polyhedron coordinates
        polyhedron_crds = self.make_polyhedron_coordinates()

        # Calculate power vertex coordinates
        self.power_crds = []
        for tri in self.hull_triangulation:
            self.power_crds.append(self.power_circumcentre(polyhedron_crds[tri[0]],polyhedron_crds[tri[1]],polyhedron_crds[tri[2]]))
        self.power_crds = np.array(self.power_crds)

        # Generate map of triangles to power coordinate id
        tri_power_map = {}
        for i,tri in enumerate(self.hull_triangulation):
            sorted_tri = np.sort(tri)
            key = '#{}#{}#{}'.format(*sorted_tri)
            tri_power_map[key] = i

        # Generate particle neighbour list
        n_crds = np.unique(self.hull_triangulation.flatten()).size # Account for (a)periodicity
        unordered_nl = [[] for i in range(n_crds)]
        for tri in self.hull_triangulation:
            unordered_nl[tri[0]].append(tri[1])
            unordered_nl[tri[0]].append(tri[2])
            unordered_nl[tri[1]].append(tri[0])
            unordered_nl[tri[1]].append(tri[2])
            unordered_nl[tri[2]].append(tri[0])
            unordered_nl[tri[2]].append(tri[1])
        for i in range(n_crds):
            unordered_nl[i] = np.array(np.unique(np.array(unordered_nl[i])),dtype=int)

        # Order neighbour list by tracing path through neighbours
        ordered_nl = []
        edges = np.zeros(n_crds,dtype=bool)
        for i in range(n_crds):
            nl = unordered_nl[i]
            pairs = []
            edge = False
            multi_cnx = False
            for j in nl:
                pair = np.intersect1d(nl,unordered_nl[j])
                if pair.size==1:
                    edge = True
                    break
                elif pair.size==2:
                    pairs.append(pair)
                elif pair.size==3:
                    multi_cnx = True
                    pairs.append(pair)
                else:
                    print("Warning: algorithm failure")
            if edge:
                ordered_nl.append(np.array([]))
                edges[i] = 1
            else:
                if multi_cnx:
                    triples = []
                    triple_pos = []
                    for j,pair in enumerate(pairs):
                        if pair.size == 3:
                            triples.append(pair)
                            triple_pos.append(j)
                    for j,triple_a in enumerate(triples):
                        for k,triple_b in enumerate(triples):
                            link = np.intersect1d(triple_a,triple_b)
                            if link.size == 1:
                                triple_a = triple_a[triple_a!=nl[triple_pos[k]]]
                                break
                        if triple_a.size!=2:
                            print("Warning: algorithm failure")
                            sys.exit()
                        pairs[triple_pos[j]] = triple_a[:]
                pairs = np.array(pairs)
                for j,val in enumerate(nl):
                    pairs[pairs==val] = j
                path = np.zeros_like(nl)
                path_nl = np.zeros_like(nl)
                path_nl[0] = nl[0]
                for j in range(1,path.size):
                    a,b = pairs[path[j-1]]
                    if a not in path:
                        path[j]=a
                        path_nl[j]=nl[a]
                    else:
                        path[j]=b
                        path_nl[j]=nl[b]
                ordered_nl.append(np.array(path_nl))

        # Generate power polygons
        power_polygons = []
        for i,nl in enumerate(ordered_nl):
            crd_ids = []
            for j in range(-1,nl.size-1):
                tri = np.sort([i,nl[j],nl[j+1]])
                key = '#{}#{}#{}'.format(*tri)
                crd_ids.append(tri_power_map[key])
            power_polygons.append(np.array(crd_ids,dtype=int))

        # Add to edge particles any whose poylgon extends beyond cell limits
        cell_limits = self.cell_limits
        cell_limits[:,0]-=10
        cell_limits[:,1]+=10
        for i,poly in enumerate(power_polygons):
            crds = self.power_crds[poly]
            if np.any(crds[:,0]<cell_limits[0,0]):
                edges[i] = 1
            elif np.any(crds[:,0]>cell_limits[0,1]):
                edges[i] = 1
            elif np.any(crds[:,1]<cell_limits[1,0]):
                edges[i] = 1
            elif np.any(crds[:,1]>cell_limits[1,1]):
                edges[i] = 1

        # Only include polygons not on edges
        self.power_polygons = []
        for i,poly in enumerate(power_polygons):
            if edges[i]==0:
                self.power_polygons.append(poly)


    def power_circumcentre(self,a,b,c):
        normal = np.cross(a,b)+np.cross(b,c)+np.cross(c,a)
        unit_normal = normal/np.sum(normal**2)
        return (-0.5 / unit_normal[2]) * unit_normal[:2]


    def visualise(self,ax=None,particles=True,polygons=True):
        """Visualise power diagram and particle positions"""

        if ax is None:
            params = {"figure.figsize": (5, 5)}
            pylab.rcParams.update(params)
            fig, ax = plt.subplots()

        if particles:
            patches_a = []
            patches_b = []
            for c in self.crds_a:
                patches_a.append(Circle(c,radius=self.r_a))
            for c in self.crds_b:
                patches_b.append(Circle(c,radius=self.r_b))
            ax.add_collection(PatchCollection(patches_a,facecolor='blue',alpha=0.5))
            ax.add_collection(PatchCollection(patches_b,facecolor='red',alpha=0.5))
            # self.ax.scatter(self.particle_crds[:,0],self.particle_crds[:,1],color='k',s=5)
        if polygons:
            patches_p = []
            for i,p in enumerate(self.power_polygons):
                if p.size>0:
                    patches_p.append(Polygon(self.power_crds[p],True))
            ax.add_collection(PatchCollection(patches_p,facecolor=(0,0,0,0),edgecolor='k'))

        ax.set_xlim(self.cell_limits[0])
        ax.set_ylim(self.cell_limits[1])

        return ax



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
        self.cell_limits = np.array([[0,cell_length],[0,cell_length]])
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
        self.cell_limits -= centre_of_mass


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
            polyhedron_crds[n:m,2] = np.sum(polyhedron_crds[n:m,:2]**2,axis=1)-self.r**2

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

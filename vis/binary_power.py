import os
import numpy as np
import sys
from numpy import linalg
from scipy.spatial import ConvexHull
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pylab as pylab
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import Axes3D


class Binary_Power:
    """Base class for power diagram for binary hard-disc system"""


    def __init__(self):
        """Base class should not be instantiated"""

        raise TypeError('Cannot instantiate base Binary_Power class')


    def generate(self):
        """Calculate Delaunay and Power diagrams"""

        self.delaunay()
        self.power()


    def delaunay(self):
        """Calculates Delaunay triangulation - all information required for network analysis"""

        # Make polyhedron coordinates for convex hull
        self.make_polyhedron_coordinates()

        # Calculate convex hull triangulation
        self.lower_hull_triangulation()

        # Calculate unordered neighbour list
        self.generate_neighbour_list()

        # Calculate number of connections to each node
        self.calculate_coordination()


    def lower_hull_triangulation(self):
        """Calculate convex hull and triangulate lower hull"""

        # Calculate convex hull
        crds = self.polyhedron_crds
        convex_hull = ConvexHull(crds)

        # Triangulate lower hull with triangles oriented anticlockwise
        hull_triangulation = []
        for i, eqn in enumerate(convex_hull.equations):
            if eqn[2] < 0.0:
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


    def analyse_delaunay(self):
        """Network and geometrical analysis of delaunay diagram"""

        # Network analysis - ring statistics
        self.cnd[self.edge_particles[:self.n]==1]=0
        self.max_cnd = np.max(self.cnd)
        k = np.arange(3,self.max_cnd+1,dtype=int)
        k_indices = np.arange(k[-1]+1) - k[0]
        p_ka = np.zeros_like(k,dtype=float)
        p_kb = np.zeros_like(k,dtype=float)
        cnds_a = self.cnd[:self.n_a]
        cnds_b = self.cnd[self.n_a:]
        for i,j in enumerate(k):
            p_ka[i] = (cnds_a==j).sum()
            p_kb[i] = (cnds_b==j).sum()
        p_kc = p_ka+p_kb

        # Network analysis - ring connections
        e_jk = np.zeros((k.size,k.size))
        for i,nl in enumerate(self.neighbour_list):
            k_i = self.cnd[i]
            if self.edge_particles[i]==0:
                for j in nl:
                    k_j = self.cnd[j]
                    if self.edge_particles[j]==0:
                        e_jk[k_indices[k_i],k_indices[k_j]] += 1
        self.delaunay_analysis = {}
        self.delaunay_analysis['k'] = k
        self.delaunay_analysis['pka'] = p_ka
        self.delaunay_analysis['pkb'] = p_kb
        self.delaunay_analysis['pkc'] = p_kc
        self.delaunay_analysis['ejk'] = e_jk


    def analyse_power(self):
        """Geometrical analysis of power diagram"""

        # Loop over polygons and calculate side lengths and angles
        lengths = []
        angles = []
        for poly in self.power_polygons:
            n = poly.size
            poly_crds = self.power_crds[poly]
            for i in range(n):
                j = (i+1)%n
                k = (i+2)%n
                v1 = poly_crds[i,:] - poly_crds[j,:]
                v2 = poly_crds[k,:] - poly_crds[j,:]
                n1 = np.sqrt(np.sum(v1**2))
                n2 = np.sqrt(np.sum(v2**2))
                u1 = v1/n1
                u2 = v2/n2
                theta = np.arccos(np.sum(u1*u2))
                lengths.append(n1)
                angles.append(theta)
        lengths = np.array(lengths)
        angles = np.array(angles)

        self.power_analysis = {}
        self.power_analysis['len_x'] = np.sum(lengths)
        self.power_analysis['len_xx'] = np.sum(lengths*lengths)
        self.power_analysis['len_n'] = np.sum(lengths.size)
        self.power_analysis['ang_x'] = np.sum(angles)
        self.power_analysis['ang_xx'] = np.sum(angles*angles)
        self.power_analysis['ang_n'] = np.sum(angles.size)


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
        """Generates power diagram"""

        # Regenerate polyhedron coordinates
        polyhedron_crds = self.polyhedron_crds

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
            # Generate adjacency matrix
            nl = unordered_nl[i]
            nl_map = {}
            for j,nlj in enumerate(nl): nl_map[nlj] = j
            adjacency_matrix = np.zeros((nl.size,nl.size),dtype=int)
            for j,nlj in enumerate(nl):
                cnxs = np.intersect1d(nl, unordered_nl[nlj])
                for k, nlk in enumerate(cnxs):
                    adjacency_matrix[nl_map[nlj],nl_map[nlk]] = 1
            row_sum = np.sum(adjacency_matrix,axis=1)
            # Order
            if np.any(row_sum==1):
                # Edge case
                ordered_nl.append(np.array([]))
                edges[i] = 1
            elif np.all(row_sum==2):
                # Simple connections
                n = np.arange(nl.size)
                indices = [0]
                j=0
                k=n[adjacency_matrix[j,:]==1][0]
                adjacency_matrix[j,k] = 0
                adjacency_matrix[k,j] = 0
                indices.append(k)
                for l in range(n.size-2):
                    j=k
                    k=n[adjacency_matrix[j,:]==1][0]
                    adjacency_matrix[k,j] = 0
                    indices.append(k)
                ordered_nl.append(np.array([nl[j] for j in indices],dtype=int))
            elif np.any(row_sum==3):
                # Untangle
                n = np.arange(nl.size)
                mask_3 = row_sum==3
                n_3 = n[mask_3]
                sub_mat = adjacency_matrix[mask_3,:]
                for i in range(sub_mat[:,0].size-1):
                    for j in range(i+1,sub_mat[:,0].size):
                        untangle = (sub_mat[i,:]*sub_mat[j,:]).sum()
                        if untangle:
                            adjacency_matrix[n_3[i],n_3[j]] = 0
                            adjacency_matrix[n_3[j],n_3[i]] = 0
                # Now simple connections
                n = np.arange(nl.size)
                indices = [0]
                j=0
                k=n[adjacency_matrix[j,:]==1][0]
                adjacency_matrix[j,k] = 0
                adjacency_matrix[k,j] = 0
                indices.append(k)
                for l in range(n.size-2):
                    j=k
                    k=n[adjacency_matrix[j,:]==1][0]
                    adjacency_matrix[k,j] = 0
                    indices.append(k)
                ordered_nl.append(np.array([nl[j] for j in indices],dtype=int))
            else:
                print('Error: algorithm failure')

        # Generate power polygons
        self.warning_particles=[]
        power_polygons = []
        for i,nl in enumerate(ordered_nl):
            crd_ids = []
            for j in range(-1,nl.size-1):
                tri = np.sort([i,nl[j],nl[j+1]])
                key = '#{}#{}#{}'.format(*tri)
                # Very rarely get edge particle without defined vertex
                try:
                    crd_ids.append(tri_power_map[key])
                except:
                    print("Warning: particle with undefined vertex {}".format(i))
                    self.warning_particles.append(i)
                    pass
            power_polygons.append(np.array(crd_ids,dtype=int))

        # Add to edge particles any whose polygon extends beyond cell limits
        # cell_limits = self.cell_limits
        # for i,poly in enumerate(power_polygons):
        #     crds = self.power_crds[poly]
        #     if np.any(crds[:,0]<cell_limits[0,0]):
        #         edges[i] = 1
        #     elif np.any(crds[:,0]>cell_limits[0,1]):
        #         edges[i] = 1
        #     elif np.any(crds[:,1]<cell_limits[1,0]):
        #         edges[i] = 1
        #     elif np.any(crds[:,1]>cell_limits[1,1]):
        #         edges[i] = 1

        # Only include polygons not on edges
        self.power_polygons = []
        for i,poly in enumerate(power_polygons[:self.n]):
            if edges[i]==0:
                self.power_polygons.append(poly)
        self.edge_particles = edges


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
            patches_warn=[]
            for w in self.warning_particles:
                if w<self.n_a:
                    patches_warn.append(Circle(self.crds_a[w],radius=self.r_b))
                else:
                    patches_warn.append(Circle(self.crds_b[w-self.n_a],radius=self.r_b))
            ax.add_collection(PatchCollection(patches_warn,facecolor='darkgreen',alpha=0.5))
        if polygons:
            patches_p = []
            for i,p in enumerate(self.power_polygons):
                if p.size>0:
                    patches_p.append(Polygon(self.power_crds[p],True))
            ax.add_collection(PatchCollection(patches_p,facecolor=(0,0,0,0),edgecolor='k'))

        ax.set_xlim(self.cell_limits[0])
        ax.set_ylim(self.cell_limits[1])
        plt.show()

        return ax


    def visualise_construction(self):
        """Visualise construction from convex hull"""

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        polyhedron_crds = self.make_polyhedron_coordinates()
        polyhedron_crds[:,2] += 200

        triangles = []
        for tri in self.hull_triangulation:
            tri_crds = np.array([polyhedron_crds[i] for i in tri])
            triangles.append(tri_crds)
        # ax.add_collection3d(Poly3DCollection(triangles, facecolor='ivory', edgecolor="k", alpha=1.0))
        for i in range(len(triangles)):
            triangles[i][:,2]=0.0
        # ax.add_collection3d(Poly3DCollection(triangles, facecolor='ivory', edgecolor="k", alpha=1.0))

        power_poly = []
        power_crds = np.zeros((self.power_crds[:,0].size,3))
        power_crds[:,:2] = self.power_crds[:,:]
        for poly in self.power_polygons:
            poly_crds = np.array([power_crds[i] for i in poly])
            power_poly.append(poly_crds)
        ax.add_collection3d(Poly3DCollection(power_poly, facecolor='ivory', edgecolor="k", alpha=0.0))
        for i in range(len(power_poly)):
            power_poly[i][:,2]=np.sum(power_poly[i][:,:2]**2+200,axis=1)
        # ax.add_collection3d(Poly3DCollection(power_poly, facecolor='ivory', edgecolor="k", alpha=0.0))

        ax.scatter(polyhedron_crds[:,0],polyhedron_crds[:,1],np.zeros(polyhedron_crds[:,0].size),c='red',s=8)
        # ax.scatter(polyhedron_crds[:,0],polyhedron_crds[:,1],polyhedron_crds[:,2],c='red',s=8)

        x_lim = (np.min(polyhedron_crds[:,0]),np.max(polyhedron_crds[:,0]))
        y_lim = (np.min(polyhedron_crds[:,1]),np.max(polyhedron_crds[:,1]))
        z_lim = (np.min(polyhedron_crds[:,2]),np.max(polyhedron_crds[:,2]))
        ax.set_xlim(x_lim)
        ax.set_ylim(y_lim)
        ax.set_zlim(z_lim)

        ax.set_axis_off()
        ax.view_init(32, -160)
        plt.savefig('power.png',dpi=400)
        plt.show()



class Periodic_Binary_Power(Binary_Power):
    """Periodic power diagram for binary system"""


    def __init__(self,crds_a,crds_b,r_a,r_b,cell_length,**kwargs):
        """Initialise with particle coordinates, radii and periodicity information"""

        # Get input variables
        self.crds_a = np.zeros_like(crds_a)
        self.crds_b = np.zeros_like(crds_b)
        self.crds_a[:,:] = crds_a
        self.crds_b[:,:] = crds_b
        self.r_a = r_a
        self.r_b = r_b
        self.cell_length = cell_length
        self.cell_limits = np.array([[-cell_length*0.5,cell_length*0.5],[-cell_length*0.5,cell_length*0.5]])
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
        # centre_of_mass = np.average(self.crds,axis=0)
        # self.crds -= centre_of_mass
        # self.crds_a -= centre_of_mass
        # self.crds_b -= centre_of_mass
        # self.cell_limits -= centre_of_mass


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

        self.polyhedron_crds = polyhedron_crds


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


    def reverse(self):

        weights = np.ones(self.n)
        power_vertices = []
        for poly in self.power_polygons:
            for i in poly:
                power_vertices.append(i%(self.n*2))
        power_vertices = np.array(power_vertices)
        print(power_vertices)
        print(np.unique(power_vertices).size)

        # x = np.zeros((self.n,2))
        # x[:,:] = self.crds[:,:]
        # mask = (self.power_crds[:,0]>self.cell_limits[0,0])*(self.power_crds[:,0]<self.cell_limits[0,1])*(self.power_crds[:,1]>self.cell_limits[1,0])*(self.power_crds[:,1]<self.cell_limits[1,1])
        # power_vertices = np.arange(self.power_crds[:,0].size)[mask]
        # triangles = self.hull_triangulation[power_vertices]%self.n
        # tri_x = triangles*2
        # tri_y = triangles*2+1
        # weights = np.ones(self.n)
        # x = x.flatten()
        # weights[0] = 1.1
        # c0 = reverse_obj_periodic(x,power_vertices,tri_x,tri_y,weights,self)
        # acceptance = 0
        # delta = 1
        # for i in range(int(1e5)):
        #     id = np.random.randint(self.n)
        #     disp = (np.random.rand()-0.5)*delta
        #     x[id] += disp
        #     c1 = reverse_obj_periodic(x,power_vertices,tri_x,tri_y,weights,self)
        #     if c1>c0:
        #         if np.random.rand()>np.exp(-(c1-c0)/1e-3):
        #             x[id] -= disp
        #         else:
        #             c0 = c1
        #             acceptance += 1
        #     else:
        #         c0 = c1
        #         acceptance += 1
        #     if(i%1e3==0):
        #         acceptance /= 1e3
        #         if acceptance<0.4:
        #             delta*=0.5
        #         acceptance = 0
        #
        #     print(c0,delta)
        #
        # bounds = [(self.cell_limits[0,0],self.cell_limits[0,1]) for i in range(self.n*2)]
        # res = minimize(reverse_obj_periodic,weights,args=(power_vertices,triangles,tri_x,tri_y,x,self),options={'maxiter' : 10})
        # print(res)


    def reverse_optimise_particle(self,xy):
        pass



def reverse_obj_periodic(weights,power_vertices,tri,tri_x,tri_y,xy,power):

    f = 0.0

    n = tri_x[:,0].size

    for i in range(n):
        x = xy[tri_x[i]]
        y = xy[tri_y[i]]
        x0,y0 = power.power_crds[power_vertices[i]]
        dx = np.abs(x-x0)
        dy = np.abs(y-y0)
        dx[dx>power.cell_length*0.5] -= power.cell_length
        dy[dy>power.cell_length*0.5] -= power.cell_length
        dSq = dx*dx + dy*dy - weights[tri[i]]**2
        f += np.abs(dSq[0]-dSq[1])
        f += np.abs(dSq[1]-dSq[2])
        f += np.abs(dSq[2]-dSq[0])
    print(f)

    return f


class Aperiodic_Binary_Power(Binary_Power):
    """Aperiodic power diagram for binary system"""


    def __init__(self,crds_a,crds_b,r_a,r_b,cell_dim,cell_cut=0,**kwargs):
        """Initialise with particle coordinates, radii and cell information"""

        # Get input variables
        self.crds_a = np.zeros_like(crds_a)
        self.crds_b = np.zeros_like(crds_b)
        self.crds_a[:,:] = crds_a
        self.crds_b[:,:] = crds_b
        self.r_a = r_a
        self.r_b = r_b
        self.cell_limits = np.array([[0,cell_dim[0]],[0,cell_dim[1]]])
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
        centre_of_mass = cell_dim/2.0
        self.crds -= centre_of_mass
        self.crds_a -= centre_of_mass
        self.crds_b -= centre_of_mass
        self.cell_limits[0] -= centre_of_mass[0]
        self.cell_limits[1] -= centre_of_mass[1]
        self.cell_limits *= (1.0-cell_cut)


    def make_polyhedron_coordinates(self):
        """Raise coordinates to form polyhedron."""

        # Copy xy coordinates and raise
        polyhedron_crds = np.zeros((self.n,3))
        polyhedron_crds[:,:2] = self.crds[:]
        polyhedron_crds[:,2] = np.sum(polyhedron_crds[:,:2]**2,axis=1)-self.r**2

        self.polyhedron_crds = polyhedron_crds


    def generate_neighbour_list(self):
        """Generate list of adjacent particles"""

        # List will be unordered i.e. neighbours cannot be traced to give polygon
        neighbour_list = [[] for i in range(self.n)]

        # Make non-unique list
        for tri in self.hull_triangulation:
            a,b,c = tri
            neighbour_list[a].append(b)
            neighbour_list[a].append(c)
            neighbour_list[b].append(a)
            neighbour_list[b].append(c)
            neighbour_list[c].append(a)
            neighbour_list[c].append(b)

        # Make unique list
        self.neighbour_list = []
        for nl in neighbour_list:
            self.neighbour_list.append(np.array(np.unique(nl),dtype=int))

if __name__ == '__main__':

    #for i in range(100):
        crds_a = np.array((np.random.rand(500,2)-0.5)*1000)
        crds_b = np.array((np.random.rand(500,2)-0.5)*1000)

        voronoi = Periodic_Binary_Power(crds_a,crds_b,1.0,1.0,1000)
        voronoi.delaunay()
        voronoi.power()
        voronoi.analyse_delaunay()
        k = voronoi.delaunay_analysis['k']
        ejk = voronoi.delaunay_analysis['ejk']
        ejk /= ejk.sum()
        qk = np.sum(ejk,axis=1)
        r = 0.0
        for i,ki in enumerate(k):
            for j,kj in enumerate(k):
                r += ki*kj*(ejk[i,j]-qk[i]*qk[j])
        r /= np.sum(k*k*qk) - np.sum(k*qk)**2
        print(r)
    #voronoi.visualise()
    #voronoi.visualise_construction()

"""Voronoi analysis of 2D system"""
import numpy as np
from scipy.spatial import Voronoi
from scipy.stats import linregress


class Colloid_Voronoi:
    """Base class for periodic and aperiodic subclasses"""


    def __init__(self):
        """Base class should not be instantiated as contains virtual methods"""
        raise TypeError('Cannot instantiate base Voronoi class')


    def network_analysis(self,assortative_mixing=True,aboav_weaire=True):
        """Analyse ring statistics and correlations"""

        # Ring statistics
        k_indices = np.arange(self.k[-1]+1) - self.k[0]
        self.p_k = np.zeros_like(self.k,dtype=float)
        for i,k in enumerate(self.k):
            k_count = (self.ring_sizes == k).sum()
            self.p_k[i] = k_count
        self.p_k /= self.p_k.sum()
        self.mean = (self.p_k*self.k).sum()
        self.var = (self.p_k*self.k*self.k).sum() - self.mean*self.mean

        # Break if only one ring size present
        if self.k.size==1:
            self.r = np.nan
            self.aw = np.zeros(3)
            self.aw[:] = np.nan
            return

        # Assortative mixing
        if assortative_mixing or aboav_weaire:
            # Calculate correlation matrix
            self.e = np.zeros((self.k.size,self.k.size),dtype=float)
            for i,cnxs in enumerate(self.ring_connections):
                k_i = self.ring_sizes[i]
                for j in cnxs:
                    k_j = self.ring_sizes[j]
                    self.e[k_indices[k_i],k_indices[k_j]] += 1.0
            self.e /= self.e.sum()
            # Calculate assortativity
            self.r = 0.0
            q = np.sum(self.e,axis=1)
            for i,k_i in enumerate(self.k):
                for j,k_j in enumerate(self.k):
                    self.r += k_i*k_j*(self.e[i,j] - q[i]*q[j])
            self.r /= (q*self.k*self.k).sum() - ((q*self.k).sum())**2

        # Aboav-Weaire
        if aboav_weaire:
            # Calculate mean ring sizes about each ring (noting some intermediate ring stats might be zero)
            m_k = []
            kk = []
            for i,k in enumerate(self.k):
                if q[i] > 0.0:
                    m_k.append((self.e[i,:]*self.k/q[i]).sum())
                    kk.append(k)
            kk = np.array(kk)
            m_k = np.array(m_k)
            # Perform linear fit
            x = 6.0*(kk-6.0)
            y = kk*m_k
            grad, y_int, r_value, p_value, std_err = linregress(x,y)
            self.aw = np.zeros(3,dtype=float)
            self.aw[0] = 1.0-grad
            self.aw[1] = y_int - 36.0
            self.aw[2] = r_value*r_value


    def voronoi_analysis(self):
        """Analyse voronoi ring areas and particle distances"""

        # Polygon areas
        self.areas = np.zeros_like(self.ring_sizes,dtype=float)
        for i,ring in enumerate(self.ring_vertex_crds):
            self.areas[i] = convex_polygon_area(ring)

        # Packing fraction (if sigma supplied)
        if self.sigma is not None:
            self.phi = np.pi*self.sigma*self.sigma/np.average(self.areas)


    def write(self,prefix='./voronoi'):
        """Write Voronoi to files"""

        # Write particle coordinates
        np.savetxt('{}_particle_crds.dat'.format(prefix),self.particle_crds)

        # Write ring coordinates
        with open('{}_ring_crds.dat'.format(prefix),'w') as f:
            for ring in self.ring_vertex_crds:
                for c in ring.flatten():
                    f.write('{:12.6f}  '.format(c))
                f.write('\n')

        # Write ring sizes
        np.savetxt('{}_ring_sizes.dat'.format(prefix),self.ring_sizes)


class Colloid_Periodic_Voronoi(Colloid_Voronoi):
    """Voronoi analysis with periodic boundary conditions"""


    def __init__(self,**kwargs):
        """Initialise with coordinates and periodicity information"""

        # Get input variables
        self.particle_crds = kwargs.get('crds')
        self.num_particle_crds = crds[:,0].size
        self.box_size = kwargs.get('box_size')
        self.sigma = kwargs.get('sigma',None)


    def calculate_voronoi(self):
        """Calculate Voronoi tessellation"""

        # Generate periodic coordinates using 9 images
        periodic_crds = self.particle_crds
        for i in [-1,0,1]:
            translation = np.zeros_like(self.particle_crds)
            translation[:,0] = i*self.box_size[0]
            for j in [-1,0,1]:
                translation[:,1] = j*self.box_size[1]
                if i==0 and j==0:
                    pass
                else:
                    periodic_crds = np.concatenate((periodic_crds,self.particle_crds+translation),axis=0)

        # Calculate voronoi with periodic coordinates
        self.voronoi = Voronoi(periodic_crds)

        # Find rings associated with original aperiodic particle points
        rings = self.voronoi.point_region[:self.num_particle_crds]

        # Get unique ring vertex coordinates and sizes
        self.ring_vertex_crds = []
        self.ring_sizes = np.zeros(self.num_particle_crds,dtype=int)
        for i,ring_id in enumerate(rings):
            crds = []
            for vertex_id in self.voronoi.regions[ring_id]:
                crds.append(self.voronoi.vertices[vertex_id])
            crds = np.array(crds)
            self.ring_vertex_crds.append(crds)
            self.ring_sizes[i] = crds[:,0].size

        # Get ring-ring connectivities
        ring_connections =[[] for i in range(self.num_particle_crds)]
        for pair in self.voronoi.ridge_points:
            id_a = pair[0]
            id_b = pair[1]
            if id_a < self.num_particle_crds:
                ring_connections[id_a].append(id_b)
            if id_b < self.num_particle_crds:
                ring_connections[id_b].append(id_a)
        self.ring_connections = []
        for cnxs in ring_connections:
            self.ring_connections.append(np.array(cnxs))
        self.ring_connections = np.array(self.ring_connections)%self.num_particle_crds

        # Get ring size range
        k_min = np.min(self.ring_sizes)
        k_max = np.max(self.ring_sizes)
        self.k = np.arange(k_min,k_max+1,dtype=int)


class Colloid_Aperiodic_Voronoi(Colloid_Voronoi):
    """Voronoi analysis with aperiodic system"""


    def __init__(self,**kwargs):
        """Initialise with coordinates"""

        # Get input variables
        self.particle_crds = kwargs.get('crds')
        self.num_particle_crds = crds[:,0].size
        self.sigma = kwargs.get('sigma',None)

        # Calculate cell dimensions
        min_x = np.min(crds[:,0])
        max_x = np.max(crds[:,0])
        min_y = np.min(crds[:,1])
        max_y = np.max(crds[:,1])
        self.box_limits = np.array([min_x, max_x, min_y, max_y])


    def calculate_voronoi(self):
        """Calculate Voronoi tessellation"""

        # Calculate voronoi with particle coordinates
        self.voronoi = Voronoi(self.particle_crds)

        # Extract rings
        rings = self.voronoi.point_region

        # Get ring vertex coordinates and sizes, noting rings which have vertices outside the box area
        keep_rings = np.ones(self.num_particle_crds,dtype=bool)
        self.ring_vertex_crds = []
        self.ring_sizes = np.zeros(self.num_particle_crds,dtype=int)
        for i,ring_id in enumerate(rings):
            crds = []
            neglect = False
            for vertex_id in self.voronoi.regions[ring_id]:
                if vertex_id == -1:
                    crds.append(crds[-1])
                    neglect = True
                else:
                    crds.append(self.voronoi.vertices[vertex_id])
            crds = np.array(crds)
            if np.any(crds[:,0]<self.box_limits[0]):
                neglect = True
            elif np.any(crds[:,0]>self.box_limits[1]):
                neglect = True
            elif np.any(crds[:,1]<self.box_limits[2]):
                neglect = True
            elif np.any(crds[:,1]>self.box_limits[3]):
                neglect = True
            if neglect:
                keep_rings[i] = 0
            self.ring_vertex_crds.append(crds)
            self.ring_sizes[i] = crds[:,0].size
        self.ring_vertex_crds = np.array(self.ring_vertex_crds)

        # Remove edge rings and reindex
        self.ring_vertex_crds = self.ring_vertex_crds[keep_rings]
        self.ring_sizes = self.ring_sizes[keep_rings]
        ring_indices = {}
        for i in range(self.num_particle_crds):
            ring_indices[i] = keep_rings[:i+1].sum()-1

        # Get ring-ring connectivities, ignoring edge rings
        ring_connections =[[] for i in range(self.ring_sizes.size)]
        for pair in self.voronoi.ridge_points:
            id_a = pair[0]
            id_b = pair[1]
            if keep_rings[id_a] and keep_rings[id_b]:
                id_c = ring_indices[id_a]
                id_d = ring_indices[id_b]
                ring_connections[id_c].append(id_d)
                ring_connections[id_d].append(id_c)
        self.ring_connections = []
        for cnxs in ring_connections:
            self.ring_connections.append(np.array(cnxs))
        self.ring_connections = np.array(self.ring_connections)

        # Get ring size range
        k_min = np.min(self.ring_sizes)
        k_max = np.max(self.ring_sizes)
        self.k = np.arange(k_min,k_max+1,dtype=int)


def convex_polygon_area(crds):
    """Area of a convex polygon using shoelace method"""

    x = crds[:,0]
    y = crds[:,1]
    area = 0.0
    area += (x[:-1]*y[1:]).sum()
    area -= (x[1:]*y[:-1]).sum()
    area += x[-1]*y[0]
    area -= x[0]*y[-1]

    return 0.5*np.abs(area)


if __name__ == "__main__":
    # crds = np.random.rand(100,2)*100.0
    crds = np.genfromtxt('test_particle_crds.dat')
    # v = Colloid_Periodic_Voronoi(crds=crds,box_size=(808.06421225, 699.8041357),sigma=10.0)
    v = Colloid_Aperiodic_Voronoi(crds=crds,sigma=10.0)
    v.calculate_voronoi()
    v.network_analysis()
    v.voronoi_analysis()
    v.write()




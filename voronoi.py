"""Voronoi analysis of 2D system"""
import numpy as np
from scipy.spatial import Voronoi
from scipy.stats import linregress


class Colloid_Periodic_Voronoi():


    def __init__(self,**kwargs):
        """Initialise with coordinates or coordinate file and periodicity information"""

        # Get coordinates
        self.acrds = kwargs.get('crds',None)
        if self.acrds is None:
            self.read_coordinates(kwargs.get('file','crds.dat'))
        self.num_acrds = self.acrds[:,0].size

        # Get periodicity information and generate periodic coordinates
        self.box_size = kwargs.get('box_size')
        self.pcrds = self.acrds
        for i in [-1,0,1]:
            translation = np.zeros_like(self.acrds)
            translation[:,0] = i*self.box_size[0]
            for j in [-1,0,1]:
                translation[:,1] = j*self.box_size[1]
                if i==0 and j==0:
                    pass
                else:
                    self.pcrds = np.concatenate((self.pcrds,self.acrds+translation),axis=0)


    def read_coordinates(self,file):
        """Read in x-y coordinates from file"""

        self.acrds = np.genfromtxt(file).astype(float)


    def calculate_voronoi(self):
        """Calculate Voronoi tesselation"""

        # Calculate voronoi
        self.voronoi = Voronoi(self.pcrds)

        # Find rings associated with aperiodic points
        arings = self.voronoi.point_region[:self.num_acrds]

        # Get unique ring vertex coordinates
        self.ring_vertex_crds = []
        self.ring_sizes = np.zeros(self.num_acrds,dtype=int)
        for i,ring_id in enumerate(arings):
            crds = []
            for vertex_id in self.voronoi.regions[ring_id]:
                crds.append(self.voronoi.vertices[vertex_id])
            crds = np.array(crds)
            self.ring_vertex_crds.append(crds)
            self.ring_sizes[i] = crds[:,0].size

        # Get ring-ring connectivities
        ring_connections =[[] for i in range(self.num_acrds)]
        for pair in self.voronoi.ridge_points:
            id_a = pair[0]
            id_b = pair[1]
            if id_a < self.num_acrds:
                ring_connections[id_a].append(id_b)
            if id_b < self.num_acrds:
                ring_connections[id_b].append(id_a)
        self.ring_connections = []
        for cnxs in ring_connections:
            self.ring_connections.append(np.array(cnxs))
        self.ring_connections = np.array(self.ring_connections)%self.num_acrds


    def network_analysis(self,assortative_mixing=True,aboav_weaire=True):
        """Analyse ring statistics and correlations"""

        # Ring statistics
        k_min = np.min(self.ring_sizes)
        k_max = np.max(self.ring_sizes)
        k_indices = np.arange(k_max+1)-k_min
        self.k = np.arange(k_min,k_max+1,dtype=int)
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


    def write(self,prefix='./voronoi'):
        """Write Voronoi to files"""

        # Write particle coordinates
        np.savetxt('{}_particle_crds.dat'.format(prefix),self.acrds)

        # Write ring coordinates
        with open('{}_ring_crds.dat'.format(prefix),'w') as f:
            for ring in self.ring_vertex_crds:
                for c in ring.flatten():
                    f.write('{:12.6f}  '.format(c))
                f.write('\n')

        # Write ring sizes
        np.savetxt('{}_ring_sizes.dat'.format(prefix),self.ring_sizes)


if __name__ == "__main__":
    crds = np.random.rand(100,2)*100.0
    v = Colloid_Periodic_Voronoi(crds=crds,box_size=(100,100))
    v.calculate_voronoi()
    v.network_analysis()
    v.write()




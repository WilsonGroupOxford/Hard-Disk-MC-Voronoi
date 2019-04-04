""""Analyse binary non-additive hard disc monte carlo simulation"""
import os
import numpy as np
from logfile import Logfile
from binary_power import Periodic_Binary_Power

class Sample:
    """Analyse coordinate snapshots from binary monte carlo simulation"""


    def __init__(self):
        """Read simulation parameters and coordinate file"""

        self.log = Logfile(name='analysis')
        self.log('Analyser For Binary Non-additive Hard Disc Monte Carlo')
        self.log('David Ormrod Morley')
        self.log('Wilson Group 2019',dash=True)
        self.log('Initialisation')
        self.read_input_file()
        self.read_simulation_files()
        self.log('Initialisation Complete',dash=True)


    def read_input_file(self):
        """Read analysis setup file"""

        # Check if input file exists in current directory, if not kill process
        if not os.path.isfile('./analysis.inpt'):
            self.log.error('Cannot find input file "analysis.inpt" in current directory')

        # Read input file and analysis options and parameters
        self.log('Reading input file')
        with open('analysis.inpt','r') as f:
            f.readline()
            self.prefix = f.readline().split()[0]
            f.readline()
            f.readline()
            self.cut_proportion = float(f.readline().split()[0])
            f.readline()
            f.readline()
            self.analysis_rdf = int(f.readline().split()[0])
            self.analysis_net = int(f.readline().split()[0])
            f.readline()
            f.readline()
            self.rdf_delta = float(f.readline().split()[0])
            self.rdf_extent = float(f.readline().split()[0])
            f.readline()
            f.readline()
            self.net_partition = int(f.readline().split()[0])

    def read_simulation_files(self):
        """Read auxilary and coordinate files from simulation"""

        # Check if simulation files exist in current directory, if not kill process
        if not os.path.isfile('{}.aux'.format(self.prefix)):
            self.log.error('Cannot find simulation file "{}.aux"'.format(self.prefix))
        if not os.path.isfile('{}.xyz'.format(self.prefix)):
            self.log.error('Cannot find simulation file "{}.xyz"'.format(self.prefix))

        # Read aux file and get simulation parameters
        self.log('Reading simulation auxilary file')
        with open('{}.aux'.format(self.prefix),'r') as f:
            l = f.readline().split()
            self.n_a = int(l[0])
            self.n_b = int(l[1])
            l = f.readline().split()
            self.r_a = float(l[0])
            self.r_b = float(l[1])
            l = f.readline().split()
            self.cell_length = float(f.readline().split()[0])
            self.frames_total = int(float(f.readline().split()[0]))
        self.n = self.n_a+self.n_b
        self.frames_cut = int(self.frames_total*self.cut_proportion)
        self.frames = self.frames_total - self.frames_cut
        self.ndensity_a = self.n_a / (self.cell_length**2)
        self.ndensity_b = self.n_b / (self.cell_length**2)
        self.rdf_extent *= self.r_b

        # Read coordinate file
        self.log('Reading simulation xyz file')
        self.crds_eqm_a = np.zeros((self.frames_cut,self.n_a,2))
        self.crds_eqm_b = np.zeros((self.frames_cut,self.n_b,2))
        self.crds_a = np.zeros((self.frames,self.n_a,2))
        self.crds_b = np.zeros((self.frames,self.n_b,2))
        with open('{}.xyz'.format(self.prefix),'r') as f:
            for i in range(self.frames_cut):
                f.readline()
                f.readline()
                for j in range(self.n_a):
                    self.crds_eqm_a[i,j,:] = np.array([float(c) for c in f.readline().split()[1:3]])
                for j in range(self.n_b):
                    self.crds_eqm_b[i,j,:] = np.array([float(c) for c in f.readline().split()[1:3]])
            for i in range(self.frames):
                f.readline()
                f.readline()
                for j in range(self.n_a):
                    self.crds_a[i,j,:] = np.array([float(c) for c in f.readline().split()[1:3]])
                for j in range(self.n_b):
                    self.crds_b[i,j,:] = np.array([float(c) for c in f.readline().split()[1:3]])


    def analyse(self):
        """Perform required analyses"""

        self.log('Analysis')
        if self.analysis_rdf:
            self.analyse_rdf()
        if self.analysis_net:
            self.analyse_network()
        self.log('Analysis complete')


    def analyse_rdf(self):
        """Calculate total and partial rdfs"""

        # RDF parameters
        self.log('Radial Distribution Functions',indent=0)
        if self.rdf_extent < self.cell_length/2:
            cut = np.floor(self.rdf_extent / self.rdf_delta) * self.rdf_delta
        else:
            cut = np.floor(0.5 * self.cell_length / self.rdf_delta) * self.rdf_delta
        cut_sq = cut*cut
        bin_lower = np.arange(0.0, cut, self.rdf_delta)
        bin_upper = bin_lower+self.rdf_delta
        bin_centre = (bin_lower+bin_upper)/2
        mic_length = self.cell_length/2
        self.log('Bins (num,delta): {} {:6.4f}'.format(bin_centre.size,self.rdf_delta),indent=1)
        self.log('Cutoff: {:6.2f}'.format(cut),indent=1)

        # Calculate raw rdf
        rdf_aa = np.zeros((bin_centre.size,2),dtype=float)
        rdf_ab = np.zeros_like(rdf_aa)
        rdf_bb = np.zeros_like(rdf_aa)
        rdf_tot = np.zeros_like(rdf_aa)
        rdf_aa[:,0] = bin_centre
        rdf_ab[:,0] = bin_centre
        rdf_bb[:,0] = bin_centre
        rdf_tot[:,0] = bin_centre
        for frame in range(self.frames):
            print(frame)
            crds_a = self.crds_a[frame]
            crds_b = self.crds_b[frame]
            # A-A
            for i in range(self.n_a-1):
                v = crds_a[i+1:]-crds_a[i]
                v[v>mic_length] -= self.cell_length
                v[v<-mic_length] += self.cell_length
                d = np.sum(v**2,axis=1)
                d = np.sqrt(d[d<cut_sq])
                bins = np.array(np.floor(d / self.rdf_delta),dtype=int)
                for bin in bins:
                    rdf_aa[bin,1] += 2
            # A-B
            for i in range(self.n_a):
                v = crds_b - crds_a[i]
                v[v>mic_length] -= self.cell_length
                v[v<-mic_length] += self.cell_length
                d = np.sum(v**2,axis=1)
                d = np.sqrt(d[d<cut_sq])
                bins = np.array(np.floor(d / self.rdf_delta),dtype=int)
                for bin in bins:
                    rdf_ab[bin,1] += 2
            # B-B
            for i in range(self.n_b-1):
                v = crds_b[i+1:]-crds_b[i]
                v[v>mic_length] -= self.cell_length
                v[v<-mic_length] += self.cell_length
                d = np.sum(v**2,axis=1)
                d = np.sqrt(d[d<cut_sq])
                bins = np.array(np.floor(d / self.rdf_delta),dtype=int)
                for bin in bins:
                    rdf_bb[bin,1] += 2
        rdf_tot[:,1] = rdf_aa[:,1] + rdf_ab[:,1] + rdf_bb[:,1]

        # Normalise
        rdf_aa[:,1] /= np.pi*(bin_upper**2-bin_lower**2)*self.frames*self.n_a*self.ndensity_a
        rdf_ab[:,1] /= np.pi*(bin_upper**2-bin_lower**2)*self.frames*(self.n_a*self.ndensity_b+self.n_b*self.ndensity_a)
        rdf_bb[:,1] /= np.pi*(bin_upper**2-bin_lower**2)*self.frames*self.n_b*self.ndensity_b
        rdf_tot[:,1] /= np.pi*(bin_upper**2-bin_lower**2)*self.frames*((self.n_a+self.n_b)/self.cell_length)**2

        # Write
        np.savetxt('{}_rdf_aa.dat'.format(self.prefix),rdf_aa)
        np.savetxt('{}_rdf_ab.dat'.format(self.prefix),rdf_ab)
        np.savetxt('{}_rdf_bb.dat'.format(self.prefix),rdf_bb)
        np.savetxt('{}_rdf.dat'.format(self.prefix),rdf_tot)
        self.log('RDFs written to files',indent=1)


    def analyse_network(self):
        """Partition space to determine neighbours and perform network analysis"""

        # Parameters based on partition type: 0 Voronoi, 1 Power/Laguerre-Voronoi
        self.log('Network Analysis',indent=0)
        if self.net_partition==0:
            self.log('Voronoi Partition',indent=1)
            f_a = open('{}_vor_a.dat'.format(self.prefix),'w')
            f_b = open('{}_vor_b.dat'.format(self.prefix),'w')
            f_c = open('{}_vor_c.dat'.format(self.prefix),'w')
            f_d = open('{}_vor_d.dat'.format(self.prefix),'w')
            f_e = open('{}_vor_e.dat'.format(self.prefix),'w')
            f_s = open('{}_vor_summary.dat'.format(self.prefix),'w')
            weight_a = 0.0
            weight_b = 0.0
        elif self.net_partition==1:
            self.log('Power Partition',indent=1)
            f_a = open('{}_pow_a.dat'.format(self.prefix),'w')
            f_b = open('{}_pow_b.dat'.format(self.prefix),'w')
            f_c = open('{}_pow_c.dat'.format(self.prefix),'w')
            f_d = open('{}_pow_d.dat'.format(self.prefix),'w')
            f_e = open('{}_pow_e.dat'.format(self.prefix),'w')
            f_s = open('{}_pow_summary.dat'.format(self.prefix),'w')
            weight_a = 2.0*self.r_a*np.sqrt(self.r_a*self.r_b)/(self.r_a+self.r_b)
            weight_b = 2.0*np.sqrt(self.r_a*self.r_b)*(1.0-self.r_a/(self.r_a+self.r_b))
        self.log('Weights (a,b): {:6.4f} {:6.4f}'.format(weight_a,weight_b),indent=2)

        # Set up total distributions i.e. using all frames
        k_max1 = 51
        k = np.arange(0,k_max1) # Limit ring sizes to 21 - increase this if get warning
        total_p_ka = np.zeros(k_max1)
        total_p_kb = np.zeros(k_max1)
        total_p_kc = np.zeros(k_max1)
        total_ejk = np.zeros((k_max1,k_max1))

        # Loop over frames, calculate diagram, extract distributions and calculate metrics
        for frame in range(self.frames):
            # Calculate diagram
            crds_a = self.crds_a[frame]
            crds_b = self.crds_b[frame]
            p_ka, p_kb, p_kc, e_jk = self.laguerre_voronoi(crds_a,crds_b,w_a=weight_a,w_b=weight_b)
            # Update total distributions
            total_p_ka += p_ka
            total_p_kb += p_kb
            total_p_kc += p_kc
            total_ejk += e_jk
            # Write frame distribution summaries
            dist_a = self.summarise_distribution(k,p_ka)
            dist_b = self.summarise_distribution(k,p_kb)
            dist_c = self.summarise_distribution(k,p_kc)
            dist_e = self.summarise_joint_distribution(k,e_jk)
            f_a.write(('{:8.6f}  '*4+'\n').format(*dist_a))
            f_b.write(('{:8.6f}  '*4+'\n').format(*dist_b))
            f_c.write(('{:8.6f}  '*4+'\n').format(*dist_c))
            f_e.write(('{:8.6f}  '*4+'\n').format(*dist_e))

        # Normalise total distributions and make summary files
        total_p_ka /= total_p_ka.sum()
        total_p_kb /= total_p_kb.sum()
        total_p_kc /= total_p_kc.sum()
        total_ejk /= total_ejk.sum()
        f_d.write(('{:10.8f}  '*k_max1+'\n').format(*total_p_ka))
        f_d.write(('{:10.8f}  '*k_max1+'\n').format(*total_p_kb))
        f_d.write(('{:10.8f}  '*k_max1+'\n').format(*total_p_kc))
        for i in range(total_ejk[:,0].size):
            f_d.write(('{:10.8f}  '*k_max1+'\n').format(*total_ejk[i,:]))
        dist_a = self.summarise_distribution(k,total_p_ka)
        dist_b = self.summarise_distribution(k,total_p_kb)
        dist_c = self.summarise_distribution(k,total_p_kc)
        dist_e = self.summarise_joint_distribution(k,total_ejk)
        f_s.write(('{:8.6f}  '*4+'\n').format(*dist_a))
        f_s.write(('{:8.6f}  '*4+'\n').format(*dist_b))
        f_s.write(('{:8.6f}  '*4+'\n').format(*dist_c))
        f_s.write(('{:8.6f}  '*4+'\n').format(*dist_e))

        # Close files
        f_a.close()
        f_b.close()
        f_c.close()
        f_d.close()
        f_e.close()
        f_s.close()
        self.log('Network analysis written to files',indent=1)


    def laguerre_voronoi(self,crds_a,crds_b,w_a=0.0,w_b=0.0):
        """Calculate Laguerre-Voronoi diagram and node distribution"""

        # Calcualte diagram
        diagram = Periodic_Binary_Power(crds_a,crds_b,w_a,w_b,self.cell_length)
        diagram.delaunay()

        # Extract distributions
        node_distributions,warning_n = diagram.delaunay_node_distributions(k_lim=50)
        edge_distributions,warning_e = diagram.delaunay_edge_distributions(k_lim=50)
        if warning_n or warning_e:
            self.log.error('Coordination number exceeded maximum expected - set higher')
        p_ka = node_distributions['p_ka']
        p_kb = node_distributions['p_kb']
        p_kc = node_distributions['p_kc']
        e_jk = edge_distributions['e_jk']

        return p_ka,p_kb,p_kc,e_jk


    def summarise_distribution(self,k,p_k):
        """Calculate mean, variance, proportion of hexagons and entropy"""

        mean = np.sum(k*p_k)
        var = np.sum(k*k*p_k)-mean*mean
        p6 = p_k[k==6]
        p = p_k[p_k>0]
        s = -np.sum(p*np.log(p))

        return np.array([mean,var,p6,s])


    def summarise_joint_distribution(self,k,e_jk):

        # Set up additional quantities
        q_k = np.sum(e_jk,axis=1)
        qq = np.zeros_like(e_jk)
        kk = np.zeros_like(e_jk)
        for i,k_i in enumerate(k):
            for j,k_j in enumerate(k):
                kk[i,j] = k_i*k_j
                qq[i,j] = q_k[i]*q_k[j]

        # Assortativity
        r = np.sum(kk*(e_jk-qq))
        ss = np.sum(k*k*q_k)-np.sum(k*q_k)**2
        r/=ss

        # Entropy of edge distribution, joint distribution and information transfer
        q = q_k[q_k>0]
        e = e_jk[e_jk>0]
        qq = qq[e_jk>0]
        s0 = -np.sum(q*np.log(q))
        s1 = -np.sum(e*np.log(e))
        it = np.sum(e*np.log(e/qq))

        return np.array([r,s0,s1,it])


if __name__ == '__main__':
    sample = Sample()
    sample.analyse()
""""Analyse binary non-additive hard disc monte carlo simulation"""
import os
import numpy as np
from logfile import Logfile


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
            f.readline()
            f.readline()
            self.rdf_delta = float(f.readline().split()[0])
            self.rdf_extent = float(f.readline().split()[0])


    def read_simulation_files(self):
        """Read auxilary and coordinate files from simulation"""

        # Check if simulation files exist in current directory, if not kill process
        if not os.path.isfile('./{}.aux'.format(self.prefix)):
            self.log.error('Cannot find simulation file "{}.aux" in current directory'.format(self.prefix))
        if not os.path.isfile('./{}.xyz'.format(self.prefix)):
            self.log.error('Cannot find simulation file "{}.xyz" in current directory'.format(self.prefix))

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


    def analyse_rdf(self):
        """Calculate total and partial rdfs"""

        # RDF parameters
        self.log('Calculating Radial Distribution Functions',indent=1)
        if self.rdf_extent < self.cell_length/2:
            cut = np.floor(self.r_b * 20 / self.rdf_delta) * self.rdf_delta
        else:
            cut = np.floor(0.5 * self.cell_length / self.rdf_delta) * self.rdf_delta
        cut_sq = cut*cut
        bin_lower = np.arange(0.0, cut, self.rdf_delta)
        bin_upper = bin_lower+self.rdf_delta
        bin_centre = (bin_lower+bin_upper)/2
        mic_length = self.cell_length/2

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
        rdf_aa[:,1] /= np.pi*self.ndensity_a*(bin_upper**2-bin_lower**2)*self.frames*self.n_a
        rdf_ab[:,1] /= np.pi*self.ndensity_b*(bin_upper**2-bin_lower**2)*self.frames*self.n_a
        rdf_bb[:,1] /= np.pi*self.ndensity_b*(bin_upper**2-bin_lower**2)*self.frames*self.n_b
        rdf_tot[:,1] /= np.pi*(bin_upper**2-bin_lower**2)*self.frames*((self.n_a+self.n_b)/self.cell_length)**2

        # Write
        np.savetxt('{}_rdf_aa.dat'.format(self.prefix),rdf_aa)
        np.savetxt('{}_rdf_ab.dat'.format(self.prefix),rdf_ab)
        np.savetxt('{}_rdf_bb.dat'.format(self.prefix),rdf_bb)
        np.savetxt('{}_rdf.dat'.format(self.prefix),rdf_tot)
        self.log('RDFs written to files',indent=2)






if __name__ == '__main__':
    sample = Sample()
    sample.analyse()
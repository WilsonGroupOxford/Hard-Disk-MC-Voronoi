"""Binary non-additive hard disc monte carlo simulation"""
import numpy as np
import os
from logfile import Logfile


# Tolerances
hd_tol = 1e-12 # Allowed hard disc overlap


class Binary_Colloid_Monte_Carlo:


    def __init__(self):
        """Initialise log file, read input parameters and initialise starting configuration"""

        # Initialise log file and write header
        self.log = Logfile(name='binary')
        self.log('Binary Non-additive Hard Disc Monte Carlo')
        self.log('David Ormrod Morley')
        self.log('Wilson Group 2019',dash=True)
        self.log('Initialisation')
        self.read_input_file()
        self.initial_configuration()
        self.log('Initialisation Complete',dash=True)


    def read_input_file(self):
        """Read binary colloid sample properties and monte carlo parameters"""

        # Check if input file exists in current directory, if not kill process
        if not os.path.isfile('./binary.inpt'):
            self.log.error('Cannot find input file "binary.inpt" in current directory')

        # Read input file
        self.log('Reading input file')
        with open('./binary.inpt','r') as f:
            f.readline()
            self.output_prefix = f.readline().split()[0]
            self.output_xyz_freq = int(f.readline().split()[0])
            f.readline()
            f.readline()
            self.n = int(f.readline().split()[0]) # Total number of particles
            self.radius_ratio = float(f.readline().split()[0]) # Sigma ratio b/a
            self.q = float(f.readline().split()[0]) # Composition
            self.phi = float(f.readline().split()[0]) # Total packing fraction
            f.readline()
            f.readline()
            self.random_seed = int(f.readline().split()[0])
            self.mc_eqm_moves = int(f.readline().split()[0])
            self.mc_sample_moves = int(f.readline().split()[0])
            self.mc_delta = float(f.readline().split()[0])
        self.log('Monte Carlo input settings',indent=1)
        self.log('Equilibrium moves: {}'.format(self.mc_eqm_moves),indent=2)
        self.log('Sampling moves: {}'.format(self.mc_sample_moves),indent=2)
        self.log('Maximum displacement: {}'.format(self.mc_delta),indent=2)
        self.log('Mersenne-Twister seed: {}'.format(self.random_seed),indent=2)
        self.log('System input settings',indent=1)
        self.log('Total particles: {}'.format(self.n),indent=2)
        self.log('Size ratio: {}'.format(self.radius_ratio),indent=2)
        self.log('Composition: {}'.format(self.q),indent=2)
        self.log('Total packing fraction: {}'.format(self.phi),indent=2)

        # Initialise additional parameters
        self.log.write('Calculating system properties')
        p_b = self.q/(self.radius_ratio**2-self.q*self.radius_ratio**2+self.q)
        self.n_b = int(np.round(p_b*self.n)) # Number of type b particles
        self.n_a = self.n - self.n_b # Number of type a particles
        self.r_a = 1.0 # Radius of a
        self.r_b = self.radius_ratio*self.r_a # Radius of b
        self.r_ab = np.sqrt(self.r_a*self.r_b) # Non-additive radius
        self.hd_aa = (2.0*self.r_a)**2 # Hard disc a-a squared interaction distance
        self.hd_ab = (2.0*self.r_ab)**2 # Hard disc a-b squared interaction distance
        self.hd_bb = (2.0*self.r_b)**2 # Hard disc b-b squared interaction distance
        self.cell_area = ((self.n_a*np.pi*self.r_a**2)+(self.n_b*np.pi*self.r_b**2))/self.phi # Periodic cell area
        self.cell_length = np.sqrt(self.cell_area) # Periodic cell length
        self.min_image_distance = self.cell_length/2.0 # Minimum image distance
        self.phi_a = (self.n_a*np.pi*self.r_a**2)/self.cell_area # Partial packing fraction of type a
        self.phi_b = (self.n_b*np.pi*self.r_b**2)/self.cell_area # Partial packing fraction of type b
        self.q = self.phi_b / self.phi # Actual composition - may differ from input as limited on number of particles
        self.log('Number of particles (A,B,Total): {} {} {}'.format(self.n_a,self.n_b,self.n),indent=1)
        self.log('Radii of particles (A,B,AB): {} {} {}'.format(self.r_a,self.r_b,self.r_ab),indent=1)
        self.log('Packing fractions (A,B,Total): {:6.4f} {:6.4f} {:6.4f}'.format(self.phi_a,self.phi_b,self.phi),indent=1)
        self.log('Composition: {:6.4f}'.format(self.q),indent=1)

        # Write aux file for analysis
        self.log('Writing auxilary file')
        with open('{}.aux'.format(self.output_prefix),'w') as f:
            f.write('{}  {}  n_a n_b  \n'.format(self.n_a,self.n_b))
            f.write('{}  {}  r_a r_b  \n'.format(self.r_a,self.r_b))
            f.write('{}  {}  phi_a phi_b  \n'.format(self.phi_a,self.phi_b))
            f.write('{}  cell length  \n'.format(self.cell_length))
            f.write('{}  sample frames  \n'.format(self.mc_sample_moves/self.output_xyz_freq+1))


    def initial_configuration(self):
        """Set up initial ordered lattice"""

        # Calculate blocks needed to construct lattice
        self.log('Constructing ordered lattice')
        block_b_capacity = 1
        block_a_capacity = int(np.floor(self.radius_ratio)**2)
        num_block_b = self.n_b//block_b_capacity
        num_block_a = self.n_a//block_a_capacity
        remainder_block_b_capacity = self.n_b%block_b_capacity
        remainder_block_a_capacity = self.n_a%block_a_capacity
        num_rem_block_b = int(remainder_block_b_capacity>0)
        num_rem_block_a = int(remainder_block_a_capacity>0)
        total_blocks = num_block_a + num_block_b + num_rem_block_a + num_rem_block_b
        total_blocks = int(np.ceil(np.sqrt(total_blocks))**2)
        spacing_blocks = total_blocks - num_block_a - num_block_b - num_rem_block_a - num_rem_block_b
        required_area = total_blocks*(2.0*self.r_b)**2
        self.log('Full blocks (A,B): {} {}'.format(num_block_a,num_block_b),indent=1)
        self.log('Partial blocks (A,B): {} {}'.format(num_rem_block_a,num_rem_block_b),indent=1)
        self.log('Empty blocks: {}'.format(spacing_blocks),indent=1)
        self.log('Area (required,cell): {:8.2f} {:8.2f}'.format(required_area,self.cell_area),indent=1)

        # Cannot construct if required area is larger than cell limits, kill process
        if required_area>self.cell_area:
            self.log.error('Algorithmic failure - cannot construct initial lattice at this composition and packing fraction')

        # Arrange blocks in unmixed state
        block_dim = np.sqrt(self.cell_area/total_blocks)
        block_a_crds = np.zeros((block_a_capacity,2))
        block_b_crds = np.zeros((block_b_capacity,2))
        k = 0
        for i in range(int(np.sqrt(block_a_capacity))):
            for j in range(int(np.sqrt(block_a_capacity))):
                block_a_crds[k,0] = j*2*self.r_a+self.r_a
                block_a_crds[k,1] = i*2*self.r_a+self.r_a
                k += 1
        k = 0
        for i in range(int(np.sqrt(block_b_capacity))):
            for j in range(int(np.sqrt(block_b_capacity))):
                block_b_crds[k,0] = j*2*self.r_b+self.r_b
                block_b_crds[k,1] = i*2*self.r_b+self.r_b
                k += 1
        self.crds_a = np.zeros((self.n_a,2))
        self.crds_b = np.zeros((self.n_b,2))
        x_blocks = int(np.sqrt(total_blocks))
        y_blocks = int(np.sqrt(total_blocks))
        # A blocks
        block = 0
        count_a = 0
        for i in range(num_block_a):
            x = block%x_blocks
            y = block//y_blocks
            block_crd = np.array([x*block_dim,y*block_dim])
            block_crds = block_a_crds + block_crd
            for crd in block_crds:
                self.crds_a[count_a,:] = crd
                count_a += 1
            block += 1
        # A remainder blocks
        for i in range(num_rem_block_a):
            x = block%x_blocks
            y = block//y_blocks
            block_crd = np.array([x*block_dim,y*block_dim])
            block_crds = block_a_crds + block_crd
            for crd in block_crds[:remainder_block_a_capacity,:]:
                self.crds_a[count_a,:] = crd
                count_a += 1
            block += 1
        # B blocks
        count_b = 0
        for i in range(num_block_b):
            x = block%x_blocks
            y = block//y_blocks
            block_crd = np.array([x*block_dim,y*block_dim])
            block_crds = block_b_crds + block_crd
            for crd in block_crds:
                self.crds_b[count_b,:] = crd
                count_b += 1
            block += 1
        # B remainder blocks
        for i in range(num_rem_block_b):
            x = block%x_blocks
            y = block//y_blocks
            block_crd = np.array([x*block_dim,y*block_dim])
            block_crds = block_b_crds + block_crd
            for crd in block_crds[:remainder_block_b_capacity,:]:
                self.crds_b[count_b,:] = crd
                count_b += 1
            block += 1
        self.log('Lattice constructed',indent=1)

        # Check no overlaps in starting configuration, if so kill and write configuration
        overlap = False
        for i in range(self.n_a):
            if overlap:
                break
            overlap = self.hard_disc_overlap(i,0)
        for i in range(self.n_b):
            if overlap:
                break
            overlap = self.hard_disc_overlap(i,1)
        if overlap:
            self.f_xyz = open('{}.xyz'.format(self.output_prefix),'w')
            self.write_xyz()
            self.f_xyz.close()
            self.log.error('Lattice contains initial overlaps, starting configuration written to .xyz file')
        else:
            self.log('Lattice checked and contains no overlapping discs',indent=1)


    def hard_disc_overlap(self,ref_id,ref_type):
        """Check if any hard-disc overlap between given particle and all others"""

        # Get coordinate of reference particle
        if ref_type == 0:
            ref_crd = self.crds_a[ref_id,:]
            hd_a = self.hd_aa
            hd_b = self.hd_ab
            a_lim = 1
            b_lim = 0
        else:
            ref_crd = self.crds_b[ref_id,:]
            hd_a = self.hd_ab
            hd_b = self.hd_bb
            a_lim = 0
            b_lim = 1

        # Calculate distances to other particles of type a, applying minimum image convention
        dx_a = self.crds_a[:, 0] - ref_crd[0]
        dy_a = self.crds_a[:, 1] - ref_crd[1]
        dx_a[dx_a < -self.min_image_distance] += self.cell_length
        dx_a[dx_a > self.min_image_distance] -= self.cell_length
        dy_a[dy_a < -self.min_image_distance] += self.cell_length
        dy_a[dy_a > self.min_image_distance] -= self.cell_length
        d_sq_a = dx_a * dx_a + dy_a * dy_a

        # Check overlap and return if found
        overlap = False
        if np.sum(d_sq_a-hd_a<-hd_tol)>a_lim:
            overlap = True
            return overlap

        # Calculate distances to other particles of type b, applying minimum image convention
        dx_b = self.crds_b[:, 0] - ref_crd[0]
        dy_b = self.crds_b[:, 1] - ref_crd[1]
        dx_b[dx_b < -self.min_image_distance] += self.cell_length
        dx_b[dx_b > self.min_image_distance] -= self.cell_length
        dy_b[dy_b < -self.min_image_distance] += self.cell_length
        dy_b[dy_b > self.min_image_distance] -= self.cell_length
        d_sq_b = dx_b * dx_b + dy_b * dy_b

        # Check overlap
        if np.sum(d_sq_b-hd_b<-hd_tol)>b_lim:
            overlap = True

        return overlap


    def monte_carlo(self):
        """Perform Monte Carlo simulation"""

        self.log('Monte Carlo Simulation')

        # Initialise Mersenne-Twister random number generator
        self.log('Initialising random number generator')
        self.random_generator = np.random.RandomState(self.random_seed)

        # Initialise output files
        self.log('Results will be written to: {}.xyz'.format(self.output_prefix))
        self.f_xyz = open('{}.xyz'.format(self.output_prefix),'w')

        # Perform equilibration moves
        self.log('Equilibrating')
        self.mc_acceptance = 0
        for i in range(1,self.mc_eqm_moves+1):
            self.monte_carlo_move()
            if i%self.output_xyz_freq==0:
                self.log('Moves: {}'.format(i),indent=1)
                print('Moves: {}'.format(i))

        # Perform sampling moves
        self.log('Sampling')
        self.mc_acceptance = 0
        self.write_xyz()
        for i in range(1,self.mc_sample_moves+1):
            self.monte_carlo_move()
            if i%self.output_xyz_freq==0:
                self.write_xyz()
                self.log('Moves and acceptance: {} {:4.3f}'.format(i,self.mc_acceptance/i),indent=1)
                print('Moves and acceptance: {} {:4.3f}'.format(i,self.mc_acceptance/i))

        # Close files
        self.log('Simulation complete',dash=True)
        self.f_xyz.close()
        self.log.close()


    def monte_carlo_move(self):
        """Single Monte Carlo displacement move"""

        # Select random particle of type a or b
        particle_id = self.random_generator.randint(0,self.n)
        if particle_id<self.n_a:

            # Get random displacement and trial coordinate
            crd_delta = self.random_generator.uniform(-self.mc_delta,self.mc_delta,size=2)
            crd_prev = np.zeros(2,dtype=float)
            crd_prev[:] = self.crds_a[particle_id,:]
            crd_trial = crd_prev + crd_delta

            # Account for periodic boundary
            crd_trial[crd_trial<0.0] += self.cell_length
            crd_trial[crd_trial>self.cell_length] -= self.cell_length

            # Set trial coordinate
            self.crds_a[particle_id,:] = crd_trial

            # Evaluate condition (hard-disc overlap)
            reject = self.hard_disc_overlap(particle_id,0)

            # Accept/reject move
            if not reject:
                self.mc_acceptance += 1
            else:
                self.crds_a[particle_id,:] = crd_prev
        else:
            particle_id -= self.n_a

            # Get random displacement and trial coordinate
            crd_delta = self.random_generator.uniform(-self.mc_delta,self.mc_delta,size=2)
            crd_prev = np.zeros(2,dtype=float)
            crd_prev[:] = self.crds_b[particle_id,:]
            crd_trial = crd_prev + crd_delta

            # Account for periodic boundary
            crd_trial[crd_trial<0.0] += self.cell_length
            crd_trial[crd_trial>self.cell_length] -= self.cell_length

            # Set trial coordinate
            self.crds_b[particle_id,:] = crd_trial

            # Evaluate condition (hard-disc overlap)
            reject = self.hard_disc_overlap(particle_id,1)

            # Accept/reject move
            if not reject:
                self.mc_acceptance += 1
            else:
                self.crds_b[particle_id,:] = crd_prev


    def write_xyz(self):
        """Write coordinates in xyz file format"""

        # Number of particles, blank line
        self.f_xyz.write('{} \n \n'.format(self.n))

        # Type X Y Z position (set z so base of particles in same plane)
        for i in range(self.n_a):
            self.f_xyz.write('A  {:12.6f}  {:12.6f}  {:12.6f}  \n'.format(self.crds_a[i,0],self.crds_a[i,1],self.r_a))
        for i in range(self.n_b):
            self.f_xyz.write('B  {:12.6f}  {:12.6f}  {:12.6f}  \n'.format(self.crds_b[i,0],self.crds_b[i,1],self.r_b))



if __name__ == "__main__":
    mc = Binary_Colloid_Monte_Carlo()
    mc.monte_carlo()

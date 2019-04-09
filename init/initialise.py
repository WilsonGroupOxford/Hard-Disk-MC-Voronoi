import os
import numpy as np
from logfile import Logfile


# Tolerances
hd_tol = 1e-12 # Allowed hard disc overlap
p_tol = 1e-2 # Acceptance probability tolerance


class Initialiser:
    """
     - Read user options from input file.
     - Calculate parameters for simulation.
     - Generate initial crystal configuration.
     - Write parameters and configuration to files.
    """


    def __init__(self,log_file=None):
        """Initialise log file if not provided and write header"""

        # Initialise log file and write header
        if log_file is None:
            self.log = Logfile(name='initialise')
            self.log("Binary Non-Additive Hard Sphere Monte Carlo Initialisation");
            self.log("Written By: David OM, Wilson Group, 2019",dash=True);
        else:
            self.log = log_file

    def complete(self):
        """Write final log and return log file"""

        self.log('Initialisation complete',dash=True)
        self.log.close()


    def read_input_file(self):
        """Read binary colloid sample properties and Monte Carlo parameters"""

        # Check if input file exists in current directory, if not kill process
        if not os.path.isfile('./initialise.inpt'):
            self.log.error('Cannot find input file "initialise.inpt" in current directory')

        # Read input file
        self.log('Reading input file')
        with open('./initialise.inpt','r') as f:
            f.readline()
            self.output_prefix = f.readline().split()[0]
            f.readline()
            f.readline()
            self.n = int(f.readline().split()[0]) # Total number of particles
            self.radius_ratio = float(f.readline().split()[0]) # Sigma ratio b/a
            self.q = float(f.readline().split()[0]) # Composition
            self.phi = float(f.readline().split()[0]) # Total packing fraction
            f.readline()
            f.readline()
            self.random_seed = int(f.readline().split()[0])
            self.cycle_preeqm = int(f.readline().split()[0])
            self.cycle_eqm = int(f.readline().split()[0])
            self.cycle_prod = int(f.readline().split()[0])
            self.write_freq = int(f.readline().split()[0])
            self.move_disp = int(f.readline().split()[0])*self.n
            self.move_clst = int(f.readline().split()[0])
        self.log('Monte Carlo input settings',indent=1)
        self.log('Mersenne-Twister seed: {}'.format(self.random_seed),indent=2)
        self.log('Pre-equilibrium cycles: {}'.format(self.cycle_preeqm),indent=2)
        self.log('Equilibrium cycles: {}'.format(self.cycle_eqm),indent=2)
        self.log('Production cycles: {}'.format(self.cycle_prod),indent=2)
        self.log('Displacement moves per cycle: {}'.format(self.move_disp),indent=2)
        self.log('Cluster moves per cycle: {}'.format(self.move_clst),indent=2)
        self.log('System input settings',indent=1)
        self.log('Total particles: {}'.format(self.n),indent=2)
        self.log('Size ratio: {}'.format(self.radius_ratio),indent=2)
        self.log('Composition: {}'.format(self.q),indent=2)
        self.log('Total packing fraction: {}'.format(self.phi),indent=2)


    def calculate_parameters(self):
        """Calculate simulation parameters"""

        # Calculatee additional parameters
        self.log.write('Calculating system parameters')
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
        self.log('Radii of particles (A,B,AB): {} {} {:6.4f}'.format(self.r_a,self.r_b,self.r_ab),indent=1)
        self.log('Packing fractions (A,B,Total): {:6.4f} {:6.4f} {:6.4f}'.format(self.phi_a,self.phi_b,self.phi),indent=1)
        self.log('Composition: {:6.4f}'.format(self.q),indent=1)


    def setup_lattice(self):
        """Set up initial lattice"""

        # Make lattice out of blocks of same size
        # Each block contains type a or b in square arrangement
        # Pad with void blocks - limits maximum possible packing fraction
        self.log('Constructing lattice')

        # Calculate occupancy of each block from size ratio
        if self.r_b-int(self.r_b)<1e-6:
            multiplier=1.0
        else:
            multiplier = np.round(1.0/(self.r_b-int(self.r_b)))
        block_a_dim = int(self.r_b*multiplier)
        block_b_dim = int(self.r_a*multiplier)
        block_a_capacity = block_a_dim**2
        block_b_capacity = block_b_dim**2
        num_block_a = self.n_a//block_a_capacity
        num_block_b = self.n_b//block_b_capacity
        excess_a = self.n_a%block_a_capacity
        excess_b = self.n_b%block_b_capacity
        num_part_block_a = int(excess_a>0)
        num_part_block_b = int(excess_b>0)
        total_blocks = num_block_a+num_block_b+num_part_block_a+num_part_block_b
        total_blocks = int(np.ceil(np.sqrt(total_blocks))**2)
        num_void_blocks = total_blocks-num_block_a-num_block_b-num_part_block_a-num_part_block_b
        check_na = num_block_a*block_a_capacity+num_part_block_a*excess_a!=self.n_a
        check_nb = num_block_b*block_b_capacity+num_part_block_b*excess_b!=self.n_b
        area_needed = total_blocks*(2.0*self.r_a*block_a_dim)**2
        # Check no errors in procedure and number of particles add up
        if check_na or check_nb:
            self.log.error('Algorithmic failure 1 - cannot construct initial lattice at this composition and packing fraction')
        # Write summary
        self.log('Full block capacity (A,B): {} {}'.format(block_a_capacity,block_b_capacity),indent=1)
        self.log('Partial block capacity (A,B): {} {}'.format(excess_a,excess_b),indent=1)
        self.log('Full blocks (A,B): {} {}'.format(num_block_a,num_block_b),indent=1)
        self.log('Partial blocks (A,B): {} {}'.format(num_part_block_b,num_part_block_b),indent=1)
        self.log('Void blocks: {}'.format(num_void_blocks),indent=1)
        self.log('Area (required,cell): {:8.2f} {:8.2f}'.format(area_needed,self.cell_area),indent=1)
        # Cannot construct if required area is larger than cell limits, kill process
        if area_needed>self.cell_area:
            self.log.error('Algorithmic failure 2 - cannot construct initial lattice at this composition and packing fraction')

        # Set up blocks
        block_dim = np.sqrt(self.cell_area/total_blocks)
        block_a_crds = np.zeros((block_a_capacity,2))
        block_b_crds = np.zeros((block_b_capacity,2))
        k = 0
        for i in range(block_a_dim):
            for j in range(block_a_dim):
                block_a_crds[k,0] = j*2*self.r_a+self.r_a
                block_a_crds[k,1] = i*2*self.r_a+self.r_a
                k += 1
        k = 0
        for i in range(block_b_dim):
            for j in range(block_b_dim):
                block_b_crds[k,0] = j*2*self.r_b+self.r_b
                block_b_crds[k,1] = i*2*self.r_b+self.r_b
                k += 1
        block_stretch = block_dim / (2*block_b_dim*self.r_b)
        block_a_crds *= block_stretch
        block_b_crds *= block_stretch
        self.crds_a = np.zeros((self.n_a,2))
        self.crds_b = np.zeros((self.n_b,2))
        x_blocks = int(np.sqrt(total_blocks))
        y_blocks = int(np.sqrt(total_blocks))
        # Arrange blocks randomly
        random_generator = np.random.RandomState(self.random_seed)
        block_list = []
        for i in range(num_block_a):
            block_list.append(0)
        for i in range(num_part_block_a):
            block_list.append(1)
        for i in range(num_block_b):
            block_list.append(2)
        for i in range(num_part_block_b):
            block_list.append(3)
        for i in range(num_void_blocks):
            block_list.append(4)
        count_a = 0
        count_b = 0
        for i in range(total_blocks):
            block_type = block_list.pop(random_generator.randint(0,total_blocks-i))
            x = i%x_blocks
            y = i//y_blocks
            block_crd = np.array([x*block_dim,y*block_dim])
            if block_type == 0:
                block_crds = block_a_crds + block_crd
                for crd in block_crds:
                    self.crds_a[count_a,:] = crd
                    count_a += 1
            elif block_type == 1:
                block_crds = block_a_crds + block_crd
                for crd in block_crds[:excess_a,:]:
                    self.crds_a[count_a,:] = crd
                    count_a += 1
            elif block_type == 2:
                block_crds = block_b_crds + block_crd
                for crd in block_crds:
                    self.crds_b[count_b,:] = crd
                    count_b += 1
            elif block_type == 3:
                block_crds = block_b_crds + block_crd
                for crd in block_crds[:excess_b,:]:
                    self.crds_b[count_b,:] = crd
                    count_b += 1
            else:
                pass

        # # Fully ordered arrangement
        # # A blocks
        # block = 0
        # count_a = 0
        # for i in range(num_block_a):
        #     x = block%x_blocks
        #     y = block//y_blocks
        #     block_crd = np.array([x*block_dim,y*block_dim])
        #     block_crds = block_a_crds + block_crd
        #     for crd in block_crds:
        #         self.crds_a[count_a,:] = crd
        #         count_a += 1
        #     block += 1
        # # A partial blocks
        # for i in range(num_part_block_a):
        #     x = block%x_blocks
        #     y = block//y_blocks
        #     block_crd = np.array([x*block_dim,y*block_dim])
        #     block_crds = block_a_crds + block_crd
        #     for crd in block_crds[:excess_a,:]:
        #         self.crds_a[count_a,:] = crd
        #         count_a += 1
        #     block += 1
        # # B blocks
        # count_b = 0
        # for i in range(num_block_b):
        #     x = block%x_blocks
        #     y = block//y_blocks
        #     block_crd = np.array([x*block_dim,y*block_dim])
        #     block_crds = block_b_crds + block_crd
        #     for crd in block_crds:
        #         self.crds_b[count_b,:] = crd
        #         count_b += 1
        #     block += 1
        # # B remainder blocks
        # for i in range(num_part_block_b):
        #     x = block%x_blocks
        #     y = block//y_blocks
        #     block_crd = np.array([x*block_dim,y*block_dim])
        #     block_crds = block_b_crds + block_crd
        #     for crd in block_crds[:excess_b,:]:
        #         self.crds_b[count_b,:] = crd
        #         count_b += 1
        #     block += 1
        # Recentre on origin
        self.crds_a -= self.min_image_distance
        self.crds_b -= self.min_image_distance
        self.log('Lattice constructed',indent=1)

        # Check for overlaps in starting configuration
        overlap = False
        for i in range(self.n_a):
            if overlap:
                break
            overlap = self.hard_disc_overlap_cell(i,0)
        for i in range(self.n_b):
            if overlap:
                break
            overlap = self.hard_disc_overlap_cell(i,1)
        if overlap:
            self.write_xyz(keyword='err')
            self.log.error('Lattice contains initial overlaps, starting configuration written to .xyz file')
        else:
            self.log('Lattice checked and contains no overlapping discs',indent=1)


    def hard_disc_overlap_cell(self,ref_id,ref_type):
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

        # Calculate distances to particles of type a, applying minimum image convention
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

        # Calculate distances to particles of type b, applying minimum image convention
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


    def write_files(self):
        """Write files for cpp Monte Carlo and analysis"""

        self.log('Writing files')

        # Write lattice
        self.write_xyz(keyword='init')
        self.log('Initial xyz configuration written')

        # Write aux file
        self.write_aux()
        self.log('Auxilary file written')


    def write_xyz(self,keyword=None):
        """Write coordinates in xyz file format"""

        # Open file
        if keyword is None:
            f_xyz = open('{}.xyz'.format(self.output_prefix),'w')
        else:
            f_xyz = open('{}_{}.xyz'.format(self.output_prefix,keyword),'w')

        # Number of particles, blank line
        f_xyz.write('{} \n \n'.format(self.n))

        # Type X Y Z position (set z so base of particles in same plane)
        for i in range(self.n_a):
            f_xyz.write('A  {:12.6f}  {:12.6f}  {:12.6f}  \n'.format(self.crds_a[i,0],self.crds_a[i,1],self.r_a))
        for i in range(self.n_b):
            f_xyz.write('B  {:12.6f}  {:12.6f}  {:12.6f}  \n'.format(self.crds_b[i,0],self.crds_b[i,1],self.r_b))

        # Close file
        f_xyz.close()


    def write_aux(self):
        """Write parameters for Monte Carlo process and analysis"""

        with open('{}.aux'.format(self.output_prefix),'w') as f:
            f.write('{}  \n'.format(self.n_a))
            f.write('{}  \n'.format(self.n_b))
            f.write('{}  \n'.format(self.r_a))
            f.write('{}  \n'.format(self.r_b))
            f.write('{}  \n'.format(self.phi_a))
            f.write('{}  \n'.format(self.phi_b))
            f.write('{}  \n'.format(self.cell_length))
            f.write('{}  \n'.format(self.cycle_preeqm))
            f.write('{}  \n'.format(self.cycle_eqm))
            f.write('{}  \n'.format(self.cycle_prod))
            f.write('{}  \n'.format(self.write_freq))
            f.write('{}  \n'.format(self.move_disp))
            f.write('{}  \n'.format(self.move_clst))
            f.write('{}  \n'.format(self.random_seed))


if __name__ == "__main__":
    init = Initialiser()
    init.read_input_file()
    init.calculate_parameters()
    init.setup_lattice()
    init.write_files()
    init.complete()

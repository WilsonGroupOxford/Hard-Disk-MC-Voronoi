"""Hard sphere monte carlo simulation"""
import numpy as np

class Colloid_Monte_Carlo:


    def __init__(self):
        """Read input file and set up regular lattice of colloid particles"""

        self.read_input_file()
        self.regular_lattice()


    def read_input_file(self):
        """Read sample properties and monte carlo parameters"""

        # Read parameters
        with open('./monte_carlo.inpt','r') as f:
            f.readline()
            self.output_prefix = f.readline().split()[0]
            f.readline()
            f.readline()
            self.n = int(f.readline().split()[0]) # Number of particles
            self.sigma = float(f.readline().split()[0]) # Particle radius
            self.phi = f.readline().split()[0] # Packing fraction
            f.readline()
            f.readline()
            self.random_seed = int(f.readline().split()[0])
            self.mc_moves = int(f.readline().split()[0])
            self.mc_max_trial_distance = float(f.readline().split()[0])

        # Initialise additional parameters
        self.diameter_sq = (2.0*self.sigma)**2
        if self.phi == 'max':
            self.phi = np.pi/(2.0*np.sqrt(3.0))
        else:
            self.phi = float(self.phi)
        self.mc_max_trial_distance *= self.sigma


    def regular_lattice(self):
        """Set up colloid particles on regular hexagonal lattice at given packing fraction"""

        # Calculate distances between particles on regular lattice and simulation box size
        dim = int(np.sqrt(self.n))
        sf = np.sqrt(np.pi/(2.0*np.sqrt(3.0)*self.phi))
        dx = self.sigma*2.0*sf
        dy = self.sigma*np.sqrt(3.0)*sf
        self.box_size = np.array([dim*dx,dim*dy])

        # Set up regular hexagonal lattice
        crds = np.zeros((dim,dim,2),dtype=float)
        x_0 = np.arange(0,dim)*dx
        x_1 = x_0+0.5*dx
        y = np.arange(0,dim)*dy
        for i in range(dim):
            if i%2==0:
                crds[i,:,0] = x_0
            else:
                crds[i,:,0] = x_1
            crds[i,:,1] = y[i]
        self.crds = np.reshape(crds,(self.n,2))


    def monte_carlo(self):
        """Perfrom Monte Carlo simulation"""

        # Initialise Mersenne-Twister random number generator
        self.random_generator = np.random.RandomState(self.random_seed)

        # Initialise Monte Carlo progress
        self.mc_acceptance = 0

        # Perform required moves
        for i in range(self.mc_moves):
            self.monte_carlo_move()
            print(i,self.mc_acceptance/(i+1))


    def monte_carlo_move(self):
        """Single Monte Carlo displacement move"""

        # Select random particle
        particle = self.random_generator.randint(0,self.n)

        # Get random displacement and trial coordinate
        delta = self.random_generator.uniform(0.0,self.mc_max_trial_distance,size=2)
        crd = self.crds[particle,:] + delta
        # Account for periodic boundary
        if crd[0] < 0.0:
            crd[0] += self.box_size[0]
        elif crd[0] > self.box_size[0]:
            crd[0] -= self.box_size[0]
        if crd[1] < 0.0:
            crd[1] += self.box_size[1]
        elif crd[1] > self.box_size[1]:
            crd[1] -= self.box_size[1]

        # Calculate distances to other particles
        dx_dy = np.abs(self.crds - crd)
        dx_dy[dx_dy[:,0]>self.box_size[0]*0.5] -= self.box_size[0]
        dx_dy[dx_dy[:,1]>self.box_size[1]*0.5] -= self.box_size[1]
        d_sq = dx_dy[:,0]**2 + dx_dy[:,1]**2

        # Accept if does not violate overlap condition
        if np.sum(d_sq<self.diameter_sq)<=1.0:
            self.crds[particle,:] = crd
            self.mc_acceptance += 1


if __name__ == "__main__":
    mc = Colloid_Monte_Carlo()
    mc.monte_carlo()
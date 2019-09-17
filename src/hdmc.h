#ifndef HDMC_HDMC_H
#define HDMC_HDMC_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <random>
#include "outputfile.h"
#include "vecf.h"
#include "vecr.h"
#include "vec_func.h"

class HDMC {
    //Hard disk Monte Carlo class

public:

    //Particle and system parameters
    int n; //number of particles
    int interaction; //additive or non-additive interactions
    double phi; //packing fraction
    double cellLen,rCellLen,cellLen_2; //cell length, reciprocal and half
    VecF<double> x,y,r; //particle x coords, y coords and radii

    //Random number generation
    mt19937 mtGen; //mersenne twister random generator
    uniform_int_distribution<int> randParticle; //uniform distribution for particle selection
    uniform_real_distribution<double> rand01; //uniform distribution between 0 and 1

    //Monte Carlo parameters
    int eqCycles,prodCycles; //number of monte carlo cycles for equilibrium and production
    double transProb,swapProb; //translation and swap move probability
    double acceptTarget; //move acceptance target
    double transDelta; //shift for translations

    //Constructor and setters
    HDMC();
    int setParticles(int num, double packFrac, int disp, VecF<double> dispParams, int interact); //set particle properties
    int setRandom(int seed); //set random number generation
    int setSimulation(int eq, int prod, double swap, double accTarg); //set simulation parameters

    //Member functions
    int initMono(VecF<double> dispParams); //initialise monodisperse particle system
    void equilibration(Logfile &logfile); //equilibration Monte Carlo
    void production(Logfile &logfile); //production Monte Carlo
    int optimalDelta(double &deltaMin, double &deltaMax, double &accProb); //find optimal translational delta
    int mcCycle(); //set of n-particle Monte Carlo moves
    void mcAdditiveMove(int &counter); //single Monte Carlo move with additive distances
    void writeXYZ(OutputFile &xyzFile); //write configuration to xyz file

};


#endif //HDMC_HDMC_H

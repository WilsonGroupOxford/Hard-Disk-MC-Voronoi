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
    double phi; //packing fraction
    double cellLen,rCellLen,cellLen_2; //cell length, reciprocal and half
    VecF<double> x,y,r; //particle x coords, y coords and radii

    //Random number generation
    mt19937 mtGen; //mersenne twister random generator
    uniform_int_distribution<int> randParticle; //uniform distribution for particle selection
//    uniform_real_distribution<double> rand; //uniform distribution for displacement move
//    uniform_real_distribution<double> randClst; //uniform distribution for cluster move

    //Monte Carlo parameters
    int peqMoves,eqMoves,prodMoves; //number of moves per particle fpr pre-equilibrium,equilibrium and production
    double transProb,swapProb; //translation and swap move probability
    double acceptTarget; //move acceptance target

    //Constructor and setters
    HDMC();
    int setParticles(int num, double packFrac, int disp, VecF<double> dispParams); //set particle properties
    int setRandom(int seed); //set random number generation
    int setSimulation(int preEq, int eq, int prod, double swap, double accTarg); //set simulation parameters

    //Member functions
    int initMono(VecF<double> dispParams); //initialise monodisperse particle system

};


#endif //HDMC_HDMC_H

#ifndef BNAHDMC_MONTECARLO_H
#define BNAHDMC_MONTECARLO_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <random>
#include "outputfile.h"
#include "vecf.h"

class MonteCarlo {

private:

    //System variables
    int nA, nB; //number of type a and b
    double rA, rB, rAB; //radii of type a,b and non-additive ab
    double hdAA,hdBB,hdAB; //hard disc a-a, b-b and a-b squared interaction distance
    double cellLen,minImageDist; //length of periodic cell and minimum image distance
    VecF< VecF<double> > crdsA, crdsB; //coordinates of type a and b

    //Monte carlo
    mt19937 mtGen; //mersenne twister random generator
    uniform_int_distribution<int> randParticle; //uniform distribution for particle selection
    uniform_real_distribution<double> randDisp; //uniform distribution for displacement move
    int nCyclePreEq,nCycleEq,nCycleProd; //number of cycles for (pre)equilibrium and production
    int cycleWriteFreq; //write coordinate frequency
    int cycleDispMoves,cycleClstMoves; //number of displacement and cluster moves per mc cycle
    double dispMoveDelta; //maximum displacement in x,y directions

    //Member functions
    double monteCarloCycle(); //single monte carlo cycle
    inline int displacementMove(); //single displacement move

public:

    //Constructors
    MonteCarlo();
    MonteCarlo(int numA, int numB, double radA, double radB, double cellLength, Logfile& logfile);

    //Member functions
    void loadCoordinates(string prefix, Logfile& logfile); //load coordinates from xyz file
    void setParameters(int seed, int preEq, int eq, int prod, int writeFreq, int disp, int clst, Logfile& logfile); //set mc parameters
    void preEquilibration(Logfile& logfile);
    void equilibration(Logfile& logfile);
    void production(Logfile& logfile);
};


#endif //BNAHDMC_MONTECARLO_H

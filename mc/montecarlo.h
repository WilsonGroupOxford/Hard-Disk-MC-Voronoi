#ifndef BNAHDMC_MONTECARLO_H
#define BNAHDMC_MONTECARLO_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <random>
#include "outputfile.h"
#include "vecf.h"
#include "vecr.h"
#include "vec_func.h"

class MonteCarlo {

private:

    //System variables
    int nA, nB; //number of type a and b
    double rA, rB, rAB; //radii of type a,b and non-additive ab
    double hdAA,hdBB,hdAB; //hard disc a-a, b-b and a-b squared interaction distance
    double cellLen,rCellLen,cellLen_2; //length/reciprocal/half of periodic cell length
    VecF<double> xA,yA,xB,yB; //x and y coordinates of type a and b

    //Monte carlo
    mt19937 mtGen; //mersenne twister random generator
    uniform_int_distribution<int> randParticle; //uniform distribution for particle selection
    uniform_real_distribution<double> randDisp; //uniform distribution for displacement move
    uniform_real_distribution<double> randClst; //uniform distribution for cluster move
    int nCyclePreEq,nCycleEq,nCycleProd; //number of cycles for (pre)equilibrium and production
    int cycleWriteFreq; //write coordinate frequency
    int cycleDispMoves,cycleClstMoves; //number of displacement and cluster moves per mc cycle
    double dispMoveDelta; //maximum displacement in x,y directions
    double hdTol; //tolerance for hard-disc overlap

    //Member functions
    VecF<double> monteCarloCycle(); //single monte carlo cycle
    inline int displacementMove(); //single displacement move
    inline int clusterMove(); //single cluster move
    inline bool nonAdditiveHardDiscOverlapAA(double &x, double &y, int &refId); //overlap of type a with all type a
    inline bool nonAdditiveHardDiscOverlapBB(double &x, double &y, int &refId); //overlap of type b with all type b
    inline bool nonAdditiveHardDiscOverlapAB(double &x, double &y); //overlap of type a with all type b
    inline bool nonAdditiveHardDiscOverlapBA(double &x, double &y); //overlap of type b with all type a
    inline void nonAdditiveHardDiscOverlapA(double &x, double &y, VecF<int> &overlapA, VecF<int> &overlapB); //overlap of type a with a and b, store ids
    inline void nonAdditiveHardDiscOverlapB(double &x, double &y, VecF<int> &overlapA, VecF<int> &overlapB); //overlap of type b with a and b, store ids

public:

    //Constructors
    MonteCarlo();
    MonteCarlo(int numA, int numB, double radA, double radB, double cellLength, int nonAdditive, Logfile& logfile);

    //Member functions
    void loadCoordinates(string prefix, Logfile& logfile); //load coordinates from xyz file
    void setParameters(int seed, int preEq, int eq, int prod, int writeFreq, int disp, int clst, Logfile& logfile); //set mc parameters
    void preEquilibration(Logfile& logfile);
    void equilibration(Logfile& logfile);
    void production(string prefix, Logfile& logfile);
    void checkAllOverlaps(Logfile& logfile);
};


#endif //BNAHDMC_MONTECARLO_H

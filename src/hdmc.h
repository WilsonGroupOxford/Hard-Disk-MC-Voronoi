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
#include "voronoi2d.h"
#include "pot2d.h"
#include "opt.h"

class HDMC {
    //Hard disk Monte Carlo class

public:

    //Particle and system parameters
    int n,nA,nB; //total number of particles, number of type A and B
    int interaction; //additive or non-additive interactions
    int dispersity; //mono/bi/poly disperse
    VecF<double> dispersityParams; //dispersity parameters
    double phi; //packing fraction
    double cellLen,rCellLen,cellLen_2; //cell length, reciprocal and half
    VecF<double> x,y,r,w; //particle x coords, y coords, radii and weights for radical voronoi

    //Random number generation
    mt19937 mtGen; //mersenne twister random generator
    uniform_int_distribution<int> randParticle; //uniform distribution for particle selection
    uniform_real_distribution<double> rand01; //uniform distribution between 0 and 1

    //Monte Carlo parameters
    int eqCycles,prodCycles; //number of monte carlo cycles for equilibrium and production
    double transProb,swapProb; //translation and swap move probability
    double acceptTarget; //move acceptance target
    double transDelta; //shift for translations

    //Analysis and output parameters
    string outputPrefix; //output file path and prefix
    bool xyzWrite,vorWrite,rdfCalc,rdfNorm,vorCalc,radCalc; //flags to write xyz/vor file, calculate and normalise RDF, calculate (radical) Voronoi
    int xyzWriteFreq, vorWriteFreq, analysisFreq; //frequency of xyz/voronoi write and analysis
    int analysisConfigs, xyzConfigs; //number of analysis/xyz configurations
    double rdfDelta; //RDF bin width
    VecF<int> rdfHist,prdfHistAA,prdfHistAB,prdfHistBB; //RDF histogram
    int maxVertices; //set maximum on number of vertices
    VecF<int> vorSizesA,vorSizesB,radSizesA,radSizesB; //voronoi/radical cell sizes
    VecF< VecF<int> > vorAdjs,radAdjs; //voronoi/radical cell size adjacencies
    VecF<double> vorAreasA,vorAreasB,radAreasA,radAreasB; //voronoi/radical cell areas by size
    VecF<int> vorNNCount,radNNCount; //voronoi/radical nearest neighbour type count
    VecF<double> vorNNSep,radNNSep; //voronoi/radical nearest neighbour separations

    //Constructor and setters
    HDMC();
    int setParticles(int num, double packFrac, int disp, VecF<double> dispParams, int interact); //set particle properties
    int setRandom(int seed); //set random number generation
    int setSimulation(int eq, int prod, double swap, double accTarg); //set simulation parameters
    int setAnalysis(string path, int xyzFreq, int vorFreq, int anFreq, int rdf, double rdfDel, int vor); //set analysis parameters

    //Member functions
    int initialiseConfiguration(Logfile &logfile, double maxIt); //generate initial particle positions
    bool rsaPositions(double maxIt); //generate positions using rsa algorithm
    void generateRandomPositions(); //generate random particle positions
    void randomPosition(double &xx, double &yy); //generate random particle position
    bool resolvePositions(); //resolve overlaps using steepest descent minimisation
    int initAnalysis(); //initialise analysis tools
    void equilibration(Logfile &logfile, OutputFile &xyzFile); //equilibration Monte Carlo
    void production(Logfile &logfile, OutputFile &xyzFile, OutputFile &vorFile, OutputFile &radFile, OutputFile &visFile); //production Monte Carlo
    void analyseConfiguration(OutputFile &vorFile, OutputFile &radFile, OutputFile &visFile, bool vorWrite); //analyse current configuration
    void calculateRDF(); //calculate RDF for current configuration
    void calculateVoronoi(OutputFile &vorFile, OutputFile &visFile, bool vis); //calculate Voronoi and analyse
    void calculateRadical(OutputFile &radFile, OutputFile &visFile, bool vis); //calculate Radical Voronoi and analyse
    VecF<double> networkAnalysis(VecF<int> &sizes, VecF< VecF<int> > &adjs); //network analysis of sizes
    int optimalDelta(double &deltaMin, double &deltaMax, double &accProb); //find optimal translational delta
    int mcCycle(); //set of n-particle Monte Carlo moves
    void mcAdditiveMove(int &counter); //single Monte Carlo move with additive distances
    void mcNonAdditiveMove(int &counter); //single Monte Carlo move with non-additive distances
    void writeXYZ(OutputFile &xyzFile); //write configuration to xyz file
    void writeVor(Voronoi2D &vor, OutputFile &visFile, int vorCode); //write voronoi visualisation
    void writeAnalysis(Logfile &logfile, OutputFile &vorFile, OutputFile &radFile, OutputFile &diaFile); //write analysis results to file
};


#endif //HDMC_HDMC_H

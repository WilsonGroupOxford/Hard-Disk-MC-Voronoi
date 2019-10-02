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
#include "voronoi.h"
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

    //Analysis and output parameters
    string outputPrefix; //output file path and prefix
    bool xyzWrite,rdfCalc,rdfNorm,vorCalc,radCalc; //flags to write xyz file, calculate and normalise RDF, calculate (radical) Voronoi
    int xyzWriteFreq, analysisFreq; //frequency of xyz write and analysis
    int analysisConfigs; //number of analysis configurations
    double rdfDelta; //RDF bin width
    VecF<int> rdfHist,prdfHistAA,prdfHistAB,prdfHistBB; //RDF histogram
    int maxVertices; //set maximum on number of vertices
    VecF<int> vorSizes,radSizes; //voronoi/radical cell sizes
    VecF< VecF<int> > vorAdjs,radAdjs; //voronoi/radical cell size adjacencies

    //Constructor and setters
    HDMC();
    int setParticles(int num, double packFrac, int disp, VecF<double> dispParams, int interact); //set particle properties
    int setRandom(int seed); //set random number generation
    int setSimulation(int eq, int prod, double swap, double accTarg); //set simulation parameters
    int setAnalysis(string path, int xyzFreq, int anFreq, int rdf, double rdfDel, int vor); //set analysis parameters

    //Member functions
    int initialiseConfiguration(Logfile &logfile); //generate initial particle positions
    void generateRandomPositions(); //generate random particle positions
    bool resolvePositions(); //resolve overlaps using steepest descent minimisation
    int initAnalysis(); //initialise analysis tools
    void equilibration(Logfile &logfile, OutputFile &xyzFile); //equilibration Monte Carlo
    void production(Logfile &logfile, OutputFile &xyzFile, OutputFile &vorFile, OutputFile &radFile); //production Monte Carlo
    void analyseConfiguration(OutputFile &vorFile, OutputFile &radFile); //analyse current configuration
    void calculateRDF(); //calculate RDF for current configuration
    void calculateVoronoi(OutputFile &vorFile); //calculate Voronoi and analyse
    VecF<double> networkAnalysis(VecF<int> &sizes, VecF< VecF<int> > &adjs); //network analysis of sizes
    int optimalDelta(double &deltaMin, double &deltaMax, double &accProb); //find optimal translational delta
    int mcCycle(); //set of n-particle Monte Carlo moves
    void mcAdditiveMove(int &counter); //single Monte Carlo move with additive distances
    void mcNonAdditiveMove(int &counter); //single Monte Carlo move with non-additive distances
    void writeXYZ(OutputFile &xyzFile); //write configuration to xyz file
    void writeAnalysis(Logfile &logfile, OutputFile &vorFile, OutputFile &radFile); //write analysis results to file
};


#endif //HDMC_HDMC_H

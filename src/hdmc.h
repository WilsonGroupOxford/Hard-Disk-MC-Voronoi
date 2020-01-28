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
#include "voronoi3d.h"
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
    double radCut; //cut for radical voronoi
    VecF<double> x,y,z,r,w; //particle x coords, y coords, z coords, radii and weights for radical voronoi
    VecF<bool> rad2DInclude; //particles to include for radical tessellation

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
    bool rdfCalc,rdfNorm, adfCalc, adfNorm; //RDF/ADF flags
    bool vorCalc2D,radCalc2D,radCalc2DCircle,vorCalc3D,radCalc3D; //Voronoi type flags
    bool visXYZ,visVor2D,visVor3D; //visualisation flags
    int analysisFreq,visFreq; //frequency of analysis/visualisation
    int analysisConfigs, xyzConfigs; //number of analysis/xyz configurations
    double rdfDelta,adfDelta; //RDF/ADF bin width
    VecF<int> rdfHist,prdfHistAA,prdfHistAB,prdfHistBB; //RDF histogram
    VecF<int> adfHistVor2D,adfHistRad2D,adfHistVor3D,adfHistRad3D; //ADF histograms
    int maxVertices; //set maximum on number of vertices
    VecF<int> vor2DSizesA,vor2DSizesB,rad2DSizesA,rad2DSizesB; //voronoi/radical cell sizes
    VecF<int> vor3DSizesA,vor3DSizesB,rad3DSizesA,rad3DSizesB; //voronoi/radical cell sizes
    VecF< VecF<int> > vor2DAdjs,rad2DAdjs; //voronoi/radical cell size adjacencies
    VecF< VecF<int> > vor3DAdjs,rad3DAdjs; //voronoi/radical cell size adjacencies
    VecF<double> vor2DAreasA,vor2DAreasB,rad2DAreasA,rad2DAreasB; //voronoi/radical cell areas by size
    VecF<double> vor3DAreasA,vor3DAreasB,rad3DAreasA,rad3DAreasB; //voronoi/radical cell areas by size
    VecF<int> vor2DNNCount,rad2DNNCount; //voronoi/radical nearest neighbour type count
    VecF<int> vor3DNNCount,rad3DNNCount; //voronoi/radical nearest neighbour type count
    VecF<double> vor2DNNSep,rad2DNNSep; //voronoi/radical nearest neighbour separations
    VecF<double> vor3DNNSep,rad3DNNSep; //voronoi/radical nearest neighbour separations

    //Constructor and setters
    HDMC();
    int setParticles(int num, double packFrac, int disp, VecF<double> dispParams, int interact); //set particle properties
    int setRandom(int seed); //set random number generation
    int setSimulation(int eq, int prod, double swap, double accTarg); //set simulation parameters
    int setAnalysis(string path, int anFreq, int rdf, double rdfDel, int adf, double adfDel, VecF<int> vor, double radZ, int visF, int vis3); //set analysis parameters

    //Member functions
    int initialiseConfiguration(Logfile &logfile, string initType, double maxIt); //generate initial particle positions
    bool rsaPositions(double maxIt); //generate positions using rsa algorithm
    void generateRandomPositions(); //generate random particle positions
    void randomPosition(double &xx, double &yy); //generate random particle position
    bool resolvePositions(); //resolve overlaps using steepest descent minimisation
    int initAnalysis(); //initialise analysis tools
    void equilibration(Logfile &logfile, OutputFile &xyzFile); //equilibration Monte Carlo
    void production(Logfile &logfile, OutputFile &xyzFile, OutputFile &vor2DFile, OutputFile &rad2DFile, OutputFile &vor3DFile, OutputFile &rad3DFile, OutputFile &vis2DFile, OutputFile &vis3DFile); //production Monte Carlo
    void analyseConfiguration(OutputFile &vor2DFile, OutputFile &rad2DFile, OutputFile &vor3DFile, OutputFile &rad3DFile, OutputFile &vis2DFile,  OutputFile &vis3DFile, bool vis); //analyse current configuration
    void calculateRDF(); //calculate RDF for current configuration
    void calculateVoronoi2D(OutputFile &vor2DFile, OutputFile &visFile, bool vis); //calculate Voronoi and analyse
    void calculateRadical2D(OutputFile &rad2DFile, OutputFile &visFile, bool vis); //calculate Radical Voronoi and analyse
    void calculateVoronoi3D(OutputFile &vor3DFile, OutputFile &vis2DFile, OutputFile &vis3DFile, bool vis); //calculate Voronoi and analyse
    void calculateRadical3D(OutputFile &rad3DFile, OutputFile &vis2DFile, OutputFile &vis3DFile, bool vis); //calculate Radical Voronoi and analyse
    VecF<double> networkAnalysis(VecF<int> &sizes, VecF< VecF<int> > &adjs); //network analysis of sizes
    int optimalDelta(double &deltaMin, double &deltaMax, double &accProb); //find optimal translational delta
    int mcCycle(); //set of n-particle Monte Carlo moves
    void mcAdditiveMove(int &counter); //single Monte Carlo move with additive distances
    void mcNonAdditiveMove(int &counter); //single Monte Carlo move with non-additive distances
    void writeXYZ(OutputFile &xyzFile); //write configuration to xyz file
    void writeVor(Voronoi2D &vor, OutputFile &vis2DFile, int vorCode, double param=0.0); //write voronoi visualisation
    void writeVor(Voronoi3D &vor, OutputFile &vis2DFile, OutputFile &vis3DFile, int vorCode, double param=0.0); //write voronoi visualisation
    void writeAnalysis(Logfile &logfile, OutputFile &vor2DFile, OutputFile &rad2DFile, OutputFile &vor3DFile, OutputFile &rad3DFile, OutputFile &diaFile); //write analysis results to file
};


#endif //HDMC_HDMC_H

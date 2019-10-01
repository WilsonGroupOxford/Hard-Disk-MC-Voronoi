#include <iostream>
#include <sstream>
#include "outputfile.h"
#include "vecf.h"
#include "hdmc.h"

using namespace std;

int main(int argc, char **argv) {

    //Set up logfile
    Logfile logfile("./hdmc.log");
    logfile.datetime("Simulation begun at: ");
    logfile.write("Hard Disk Monte Carlo");
    logfile.write("Written By: David OM, Wilson Group, 2019");
    logfile.separator();

    //Read input parameters
    logfile.write("Reading input parameters");
    ifstream inputFile("./hdmc.inpt", ios::in);
    if(!inputFile.good()) logfile.criticalError("Cannot find input file hdmc.inpt");
    string skip,line;
    ++logfile.currIndent;
    //Particle parameters
    int n; //number of particles
    string disp; //particle dispersity (mono,bi,poly-disperse radii)
    VecF<double> dispParams; //dispersity parameters
    string interaction; //particle interactions (additive,non-additive distance)
    double packFrac; //packing fraction
    int dispCode,intCode; //numeric codes for dispersity/interaction types
    logfile.write("Reading particle parameters");
    ++logfile.currIndent;
    for(int i=0; i<3; ++i) getline(inputFile,skip);
    getline(inputFile,line);
    istringstream(line)>>n;
    logfile.write("Number of particles:",n);
    getline(inputFile,line);
    istringstream(line)>>disp;
    logfile.write("Particle dispersity:",disp);
    if(disp.substr(0,2)=="bi") cout<<"XXX"<<endl;
    else if(disp.substr(0,4)=="mono"){
        dispParams=VecF<double>(1);
        getline(inputFile,line);
        istringstream(line)>>dispParams[0];
        dispCode=1;
        logfile.write("Particle radii:",dispParams[0]);
    }
    else if(disp.substr(0,4)=="poly") cout<<"ZZZ"<<endl;
    else logfile.criticalError("Error reading particle dispersity code");
    getline(inputFile,line);
    istringstream(line)>>interaction;
    if(interaction.substr(0,3)=="add"){
        intCode=0;
        logfile.write("Particle interactions: additive");
    }
    else if(interaction.substr(0,3)=="nonadd"){
        intCode=1;
        logfile.write("Particle interactions: non-additive");
    }
    getline(inputFile,line);
    istringstream(line)>>packFrac;
    logfile.write("Packing fraction:",packFrac);
    --logfile.currIndent;
    //Simulation parameters
    logfile.write("Reading simulation parameters");
    ++logfile.currIndent;
    for(int i=0; i<2; ++i) getline(inputFile,skip);
    int randomSeed; //seed for random number generator
    int eqCycles, prodCycles; //number of equilibration and production cycles
    double swapProb,accTarget; //swap probability and acceptance probability target
    getline(inputFile,line);
    istringstream(line)>>randomSeed;
    logfile.write("Random seed:",randomSeed);
    getline(inputFile,line);
    istringstream(line)>>eqCycles;
    logfile.write("Equilibration moves per particle:",eqCycles);
    getline(inputFile,line);
    istringstream(line)>>prodCycles;
    logfile.write("Production moves per particle:",prodCycles);
    getline(inputFile,line);
    istringstream(line)>>swapProb;
    logfile.write("Swap move probability:",swapProb);
    getline(inputFile,line);
    istringstream(line)>>accTarget;
    logfile.write("Target acceptance probability:",accTarget);
    --logfile.currIndent;
    //Analysis and output parameters
    logfile.write("Reading analysis and output parameters");
    ++logfile.currIndent;
    for(int i=0; i<2; ++i) getline(inputFile,skip);
    string outputPrefix;
    int xyzWriteFreq,analysisFreq,rdfAnalysis,vorAnalysis;
    double rdfDelta;
    getline(inputFile,line);
    istringstream(line)>>outputPrefix;
    logfile.write("Output prefix:",outputPrefix);
    getline(inputFile,line);
    istringstream(line)>>xyzWriteFreq;
    logfile.write("XYZ file write frequency (cycles):",xyzWriteFreq);
    getline(inputFile,line);
    istringstream(line)>>analysisFreq;
    logfile.write("Analysis frequency (cycles):",analysisFreq);
    getline(inputFile,line);
    istringstream(line)>>rdfAnalysis;
    logfile.write("Radial distribution function calculation:",rdfAnalysis);
    getline(inputFile,line);
    istringstream(line)>>rdfDelta;
    logfile.write("Radial distribution function bin width:",rdfDelta);
    getline(inputFile,line);
    istringstream(line)>>vorAnalysis;
    logfile.write("Voronoi analysis:",vorAnalysis);
    --logfile.currIndent;
    logfile.separator();

    //Initialise Monte Carlo simulation
    logfile.write("Initialising Monte Carlo simulation");
    ++logfile.currIndent;
    int latticeInit;
    HDMC simulation;
    latticeInit=simulation.setParticles(n,packFrac,dispCode,dispParams,intCode);
    if(latticeInit==1) logfile.criticalError("Packing fraction too high to form initial lattice");
    logfile.write("Starting configuration constructed");
    simulation.setRandom(randomSeed);
    logfile.write("Random number generators initialised");
    simulation.setSimulation(eqCycles,prodCycles,swapProb,accTarget);
    logfile.write("Simulation parameters set");
    simulation.setAnalysis(outputPrefix,xyzWriteFreq,analysisFreq,rdfAnalysis,rdfDelta,vorAnalysis);
    logfile.write("Analysis and write parameters set");
    --logfile.currIndent;
    logfile.separator();

    //Set up output files
    OutputFile xyzFile(outputPrefix+".xyz");
    OutputFile vorFile(outputPrefix+"_vor.dat");
    OutputFile radFile(outputPrefix+"_rad.dat");

    //Run Monte Carlo simulation (xyz written only for production atm)
    simulation.equilibration(logfile,xyzFile);
    simulation.production(logfile,xyzFile,vorFile,radFile);

    //Write analysis to files
    simulation.writeAnalysis(logfile,vorFile,radFile);

    return 0;
}
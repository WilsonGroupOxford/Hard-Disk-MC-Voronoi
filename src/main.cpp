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
    ifstream inputAuxFile("./hdmc.inpt", ios::in);
    if(!inputAuxFile.good()) logfile.criticalError("Cannot find input file hdmc.inpt");
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
    for(int i=0; i<3; ++i) getline(inputAuxFile,skip);
    getline(inputAuxFile,line);
    istringstream(line)>>n;
    logfile.write("Number of particles:",n);
    getline(inputAuxFile,line);
    istringstream(line)>>disp;
    logfile.write("Particle dispersity:",disp);
    if(disp.substr(0,2)=="bi") cout<<"XXX"<<endl;
    else if(disp.substr(0,4)=="mono"){
        dispParams=VecF<double>(1);
        getline(inputAuxFile,line);
        istringstream(line)>>dispParams[0];
        dispCode=1;
        logfile.write("Particle radii:",dispParams[0]);
    }
    else if(disp.substr(0,4)=="poly") cout<<"ZZZ"<<endl;
    else logfile.criticalError("Error reading particle dispersity code");
    getline(inputAuxFile,line);
    istringstream(line)>>interaction;
    if(interaction.substr(0,3)=="add"){
        intCode=0;
        logfile.write("Particle interactions: additive");
    }
    else if(interaction.substr(0,3)=="nonadd"){
        intCode=1;
        logfile.write("Particle interactions: non-additive");
    }
    getline(inputAuxFile,line);
    istringstream(line)>>packFrac;
    logfile.write("Packing fraction:",packFrac);
    logfile.currIndent -= 2;
    logfile.separator();
    //Simulation parameters
    logfile.write("Reading simulation parameters");
    ++logfile.currIndent;
    for(int i=0; i<2; ++i) getline(inputAuxFile,skip);
    int randomSeed; //seed for random number generator
    int preEqMoves, eqMoves, prodMoves; //number of pre-equilibration,equilibration and production moves per particle
    double swapProb,accTarget; //swap probability and acceptance probability target
    getline(inputAuxFile,line);
    istringstream(line)>>randomSeed;
    logfile.write("Random seed:",randomSeed);
    getline(inputAuxFile,line);
    istringstream(line)>>preEqMoves;
    logfile.write("Pre-equilibration moves per particle:",preEqMoves);
    getline(inputAuxFile,line);
    istringstream(line)>>eqMoves;
    logfile.write("Equilibration moves per particle:",eqMoves);
    getline(inputAuxFile,line);
    istringstream(line)>>prodMoves;
    logfile.write("Production moves per particle:",prodMoves);
    getline(inputAuxFile,line);
    istringstream(line)>>swapProb;
    logfile.write("Swap move probability:",swapProb);
    getline(inputAuxFile,line);
    istringstream(line)>>accTarget;
    logfile.write("Target acceptance probability:",accTarget);
    --logfile.currIndent;
    logfile.separator();

    //Initialise Monte Carlo simulation
    logfile.write("Initialising Monte Carlo simulation");
    ++logfile.currIndent;
    int latticeInit;
    HDMC simulation;
    latticeInit=simulation.setParticles(n,packFrac,dispCode,dispParams);
    if(latticeInit==1) logfile.criticalError("Packing fraction too high to form initial lattice");
    logfile.write("Starting configuration constructed");
    simulation.setRandom(randomSeed);
    logfile.write("Random number generators initialised");
    simulation.setSimulation(preEqMoves,eqMoves,prodMoves,swapProb,accTarget);
    logfile.write("Simulation parameters set");
    --logfile.currIndent;
    logfile.separator();



    return 0;
}
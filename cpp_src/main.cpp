#include <iostream>
#include <sstream>
#include "outputfile.h"
#include "montecarlo.h"

using namespace std;

int main(int argc, char **argv) {

    //Set up logfile
    Logfile logfile("./mc.log");
    logfile.datetime("Simulation begun at: ");
    logfile.write("Binary Non-Additive Hard Sphere Monte Carlo");
    logfile.write("Written By: David OM, Wilson Group, 2019");
    logfile.separator();

    //Read auxilary file
    logfile.write("Reading auxilary parameters");
    ifstream inputAuxFile(string(argv[1])+".aux", ios::in);
    if(!inputAuxFile.good()) logfile.criticalError("Cannot find auxilary file");
    string skip,line;
    int nA,nB; //number of type a and b
    double rA,rB; //radius of type a and b
    double cellLength; //periodic cell length
    int nCyclePreEq, nCycleEq, nCycleProd; //number of pre-equilibrium, equilibrium and production cycles
    int cycleWriteFreq; //write frequency of sample
    int cycleDispMoves,cycleClstMoves; //number of displacement and cluster moves per cycle
    int randomSeed; //for random number generator
    getline(inputAuxFile,line);
    istringstream(line)>>nA;
    getline(inputAuxFile,line);
    istringstream(line)>>nB;
    getline(inputAuxFile,line);
    istringstream(line)>>rA;
    getline(inputAuxFile,line);
    istringstream(line)>>rB;
    getline(inputAuxFile,skip);
    getline(inputAuxFile,skip);
    getline(inputAuxFile,line);
    istringstream(line)>>cellLength;
    getline(inputAuxFile,line);
    istringstream(line)>>nCyclePreEq;
    getline(inputAuxFile,line);
    istringstream(line)>>nCycleEq;
    getline(inputAuxFile,line);
    istringstream(line)>>nCycleProd;
    getline(inputAuxFile,line);
    istringstream(line)>>cycleWriteFreq;
    getline(inputAuxFile,line);
    istringstream(line)>>cycleDispMoves;
    getline(inputAuxFile,line);
    istringstream(line)>>cycleClstMoves;
    getline(inputAuxFile,line);
    istringstream(line)>>randomSeed;

    //Initialise Monte Carlo Simulation
    MonteCarlo simulation(nA,nB,rA,rB,cellLength,logfile);
    simulation.loadCoordinates(string(argv[1]),logfile);
    simulation.setParameters(randomSeed,nCyclePreEq,nCycleEq,nCycleProd,cycleWriteFreq,cycleDispMoves,cycleClstMoves,logfile);
    simulation.checkAllOverlaps(logfile);
    logfile.separator();

    //Pre-equilibration to determine ideal move displacement
    simulation.preEquilibration(logfile);
    logfile.separator();

    //Equilibration
    simulation.equilibration(logfile);
    logfile.separator();

    //Production
    simulation.production(string(argv[1]),logfile);
    logfile.separator();

    //Check final configuration
    simulation.checkAllOverlaps(logfile);
    logfile.separator();

    //Close files
    logfile.datetime("Simulation complete at: ");
    return 0;
}
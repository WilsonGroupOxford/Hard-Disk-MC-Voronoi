#include <iostream>
#include <sstream>
#include <cmath>
#include "outputfile.h"

using namespace std;

int main(int argc, char **argv) {

    //Set up logfile
    Logfile logfile("./analysis.log");
    logfile.datetime("Simulation begun at: ");
    logfile.write("Binary Non-Additive Hard Sphere Monte Carlo Analysis");
    logfile.write("Written By: David OM, Wilson Group, 2019");
    logfile.separator();

    //Read input file
    logfile.write("Reading input file");
    ifstream inputFile("analysis.inpt", ios::in);
    if(!inputFile.good()) logfile.criticalError("Cannot find input file");
    string skip,line;
    double cutProportion; //proportion of configurations to cut
    int rdfAnalysis,vorAnalysis,radAnalysis; //analysis types
    getline(inputFile,skip);
    getline(inputFile,line);
    istringstream(line)>>cutProportion;
    getline(inputFile,skip);
    getline(inputFile,skip);
    getline(inputFile,line);
    istringstream(line)>>rdfAnalysis;
    getline(inputFile,line);
    istringstream(line)>>vorAnalysis;
    getline(inputFile,line);
    istringstream(line)>>radAnalysis;

    //Read auxilary file
    logfile.write("Reading auxilary parameters");
    ifstream inputAuxFile(string(argv[1])+".aux", ios::in);
    if(!inputAuxFile.good()) logfile.criticalError("Cannot find auxilary file");
    int nA,nB; //number of type a and b
    double rA,rB; //radius of type a and b
    double cellLength; //periodic cell length
    int nCyclePreEq, nCycleEq, nCycleProd; //number of pre-equilibrium, equilibrium and production cycles
    int cycleWriteFreq; //write frequency of sample
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
    int nConfigs = floor(nCycleProd/cycleWriteFreq)+1;

    //Summarise selection in log file
    logfile.write("Analysis selections");
    ++logfile.currIndent;
    logfile.write("RDF: ",rdfAnalysis);
    logfile.write("Voronoi: ",vorAnalysis);
    logfile.write("Radical/Laguerre: ",radAnalysis);
    --logfile.currIndent;

    logfile.separator();
    return 0;
}

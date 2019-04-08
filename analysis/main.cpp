#include <iostream>
#include <sstream>
#include <cmath>
#include "outputfile.h"
#include "configuration.h"

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
    double rdfDelta,rdfExtent; //bin increment and maximum value
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
    getline(inputFile,skip);
    getline(inputFile,skip);
    getline(inputFile,line);
    istringstream(line)>>rdfDelta;
    getline(inputFile,line);
    istringstream(line)>>rdfExtent;

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

    //Open input file
    logfile.write("Analysis");
    string prefix = string(argv[1]);
    ifstream xyzFile(prefix+"_prod.xyz", ios::in);
    if(!xyzFile.good()) logfile.criticalError("Cannot find input file");
    logfile.write("XYZ file opened");

    //Set up configuration and output files
    ofstream vorFileA,vorFileB,vorFileC;
    Configuration config(nA,nB,rA,rB,cellLength); //set up configuration
    if(rdfAnalysis) config.setRdf(rdfDelta,rdfExtent);
    if(vorAnalysis){
        vorFileA = ofstream(prefix+"_vor_pa.dat", ios::in | ios::trunc);
        vorFileB = ofstream(prefix+"_vor_pb.dat", ios::in | ios::trunc);
        vorFileC = ofstream(prefix+"_vor_pc.dat", ios::in | ios::trunc);
        config.setVoronoi(logfile);
    }

    //Analyse frame by frame
    for(int i=0; i<nConfigs; ++i){
        config.setCoordinates(xyzFile,logfile);
        ++logfile.currIndent;
        if(rdfAnalysis) config.rdf(logfile);
        if(vorAnalysis) config.voronoi(vorFileA,vorFileB,vorFileC,logfile);
        --logfile.currIndent;
    }

    //Finalise analyses
    config.rdfFinalise(prefix,logfile);

    logfile.separator();
    return 0;
}

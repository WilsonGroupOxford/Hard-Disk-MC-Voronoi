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
    ofstream vorFilePKA,vorFilePKB,vorFilePKC,vorFileARA,vorFileARB,vorFileNet;
    Configuration config(nA,nB,rA,rB,cellLength); //set up configuration
    if(rdfAnalysis) config.setRdf(rdfDelta,rdfExtent);
    if(vorAnalysis){
        vorFilePKA = ofstream(prefix+"_vor_pka.dat", ios::in | ios::trunc);
        vorFilePKB = ofstream(prefix+"_vor_pkb.dat", ios::in | ios::trunc);
        vorFilePKC = ofstream(prefix+"_vor_pkc.dat", ios::in | ios::trunc);
        vorFileARA = ofstream(prefix+"_vor_ara.dat", ios::in | ios::trunc);
        vorFileARB = ofstream(prefix+"_vor_arb.dat", ios::in | ios::trunc);
        vorFileNet = ofstream(prefix+"_vor_net.dat", ios::in | ios::trunc);
        vorFilePKA << fixed << showpoint << setprecision(8);
        vorFilePKB << fixed << showpoint << setprecision(8);
        vorFilePKC << fixed << showpoint << setprecision(8);
        vorFileARA << fixed << showpoint << setprecision(8);
        vorFileARB << fixed << showpoint << setprecision(8);
        vorFileNet << fixed << showpoint << setprecision(8);
        config.setVoronoi(logfile);
    }

    //Analyse frame by frame
    for(int i=0; i<nConfigs; ++i){
        config.setCoordinates(xyzFile,logfile);
        ++logfile.currIndent;
        if(rdfAnalysis) config.rdf(logfile);
        if(vorAnalysis) config.voronoi(vorFilePKA,vorFilePKB,vorFilePKC,vorFileARA,vorFileARB,vorFileNet,logfile);
        logfile.write("Configuration: ",i);
        cout<<i<<endl;
        --logfile.currIndent;
    }

    //Finalise analyses
    if(rdfAnalysis) config.rdfFinalise(prefix,logfile);
    if(vorAnalysis) config.voronoiFinalise(vorFilePKA,vorFilePKB,vorFilePKC,vorFileARA,vorFileARB,vorFileNet,logfile);

    logfile.separator();
    return 0;
}

#include <iostream>
#include <sstream>
#include <stdio.h>
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
    int nConfigs; //total number of configurations to analyse
    int rdfAnalysis,vorAnalysis,radAnalysis; //analysis types
    double rdfDelta,rdfExtent; //bin increment and maximum value
    getline(inputFile,skip);
    getline(inputFile,line);
    istringstream(line)>>nConfigs;
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
    int nonAdditive; //non-additivity enabled
    int nA,nB; //number of type a and b
    double rA,rB; //radius of type a and b
    double cellLength; //periodic cell length
    int nCyclePreEq, nCycleEq, nCycleProd; //number of pre-equilibrium, equilibrium and production cycles
    int cycleWriteFreq; //write frequency of sample
    getline(inputAuxFile,line);
    istringstream(line)>>nonAdditive;
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
    //read in number of configurations from input file to allow for merging of results files
//    int nConfigs = floor(nCycleProd/cycleWriteFreq)+1;

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
    ifstream xyzFile(prefix+".xyz", ios::in);
    if(!xyzFile.good()) logfile.criticalError("Cannot find input file");
    logfile.write("XYZ file opened");

    //Set up configuration and output files
    ofstream vorFilePKA(prefix+"_vor_pka.dat", ios::in | ios::trunc);
    ofstream vorFilePKB(prefix+"_vor_pkb.dat", ios::in | ios::trunc);
    ofstream vorFilePKC(prefix+"_vor_pkc.dat", ios::in | ios::trunc);
    ofstream vorFileEJK(prefix+"_vor_ejk.dat", ios::in | ios::trunc);
    ofstream vorFileARA(prefix+"_vor_ara.dat", ios::in | ios::trunc);
    ofstream vorFileARB(prefix+"_vor_arb.dat", ios::in | ios::trunc);
    ofstream vorFileELA(prefix+"_vor_ela.dat", ios::in | ios::trunc);
    ofstream vorFileELB(prefix+"_vor_elb.dat", ios::in | ios::trunc);
    ofstream vorFileLAN(prefix+"_vor_lan.dat", ios::in | ios::trunc);
    ofstream vorFileNet(prefix+"_vor_net.dat", ios::in | ios::trunc);
    ofstream radFilePKA(prefix+"_rad_pka.dat", ios::in | ios::trunc);
    ofstream radFilePKB(prefix+"_rad_pkb.dat", ios::in | ios::trunc);
    ofstream radFilePKC(prefix+"_rad_pkc.dat", ios::in | ios::trunc);
    ofstream radFileEJK(prefix+"_rad_ejk.dat", ios::in | ios::trunc);
    ofstream radFileARA(prefix+"_rad_ara.dat", ios::in | ios::trunc);
    ofstream radFileARB(prefix+"_rad_arb.dat", ios::in | ios::trunc);
    ofstream radFileELA(prefix+"_rad_ela.dat", ios::in | ios::trunc);
    ofstream radFileELB(prefix+"_rad_elb.dat", ios::in | ios::trunc);
    ofstream radFileLAN(prefix+"_rad_lan.dat", ios::in | ios::trunc);
    ofstream radFileNet(prefix+"_rad_net.dat", ios::in | ios::trunc);
    Configuration config(nA,nB,rA,rB,cellLength,nonAdditive); //set up configuration
    if(rdfAnalysis) config.setRdf(rdfDelta,rdfExtent);
    if(vorAnalysis){
        vorFilePKA << fixed << showpoint << setprecision(8);
        vorFilePKB << fixed << showpoint << setprecision(8);
        vorFilePKC << fixed << showpoint << setprecision(8);
        vorFileEJK << fixed << showpoint << setprecision(8);
        vorFileARA << fixed << showpoint << setprecision(8);
        vorFileARB << fixed << showpoint << setprecision(8);
        vorFileELA << fixed << showpoint << setprecision(8);
        vorFileELB << fixed << showpoint << setprecision(8);
        vorFileLAN << fixed << showpoint << setprecision(8);
        vorFileNet << fixed << showpoint << setprecision(8);
        config.setVoronoi(logfile);
    }
    if(radAnalysis){
        radFilePKA << fixed << showpoint << setprecision(8);
        radFilePKB << fixed << showpoint << setprecision(8);
        radFilePKC << fixed << showpoint << setprecision(8);
        radFileEJK << fixed << showpoint << setprecision(8);
        radFileARA << fixed << showpoint << setprecision(8);
        radFileARB << fixed << showpoint << setprecision(8);
        radFileELA << fixed << showpoint << setprecision(8);
        radFileELB << fixed << showpoint << setprecision(8);
        radFileLAN << fixed << showpoint << setprecision(8);
        radFileNet << fixed << showpoint << setprecision(8);
        config.setRadical(logfile);
    }

    //Analyse frame by frame
    for(int i=0; i<nConfigs; ++i){
        config.setCoordinates(xyzFile,logfile);
        ++logfile.currIndent;
        if(rdfAnalysis) config.rdf(logfile);
        if(vorAnalysis) config.voronoi(vorFilePKA,vorFilePKB,vorFilePKC,vorFileEJK,vorFileARA,vorFileARB,vorFileELA,vorFileELB,vorFileNet,logfile);
        if(radAnalysis) config.radical(radFilePKA,radFilePKB,radFilePKC,radFileEJK,radFileARA,radFileARB,radFileELA,radFileELB,radFileNet,logfile);
        logfile.write("Configuration: ",i);
        cout<<i<<endl;
        --logfile.currIndent;
    }

    //Finalise analyses
    if(rdfAnalysis) config.rdfFinalise(prefix,logfile);
    if(vorAnalysis) config.voronoiFinalise(vorFilePKA,vorFilePKB,vorFilePKC,vorFileEJK,vorFileARA,vorFileARB,vorFileELA,vorFileELB,vorFileLAN,vorFileNet,logfile);
    if(radAnalysis) config.radicalFinalise(radFilePKA,radFilePKB,radFilePKC,radFileEJK,radFileARA,radFileARB,radFileELA,radFileELB,radFileLAN,radFileNet,logfile);

    //Close files and remove any unnecessary files
    vorFilePKA.close();
    vorFilePKB.close();
    vorFilePKC.close();
    vorFileEJK.close();
    vorFileARA.close();
    vorFileARB.close();
    vorFileLAN.close();
    vorFileNet.close();
    radFilePKA.close();
    radFilePKB.close();
    radFilePKC.close();
    radFileEJK.close();
    radFileARA.close();
    radFileARB.close();
    radFileLAN.close();
    radFileNet.close();
    if(!vorAnalysis){
        string removeFiles = "rm " + prefix + "_vor*.dat";
        const char *rm = (removeFiles).c_str();
        system(rm);
    }
    if(!radAnalysis){
        string removeFiles = "rm " + prefix + "_rad*.dat";
        const char *rm = (removeFiles).c_str();
        system(rm);
    }

    logfile.separator();
    return 0;
}

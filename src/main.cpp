#include <iostream>
#include <sstream>
#include "outputfile.h"
#include "vecf.h"

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
        istringstream(line)>>disp[0];
        logfile.write("Particle radii:",disp[0]);
    }
    else if(disp.substr(0,4)=="poly") cout<<"ZZZ"<<endl;
    else logfile.criticalError("Error reading particle dispersity code");
    getline(inputAuxFile,line);
    istringstream(line)>>interaction;
    if(interaction.substr(0,3)=="add") logfile.write("Particle interactions: additive");
    else if(interaction.substr(0,3)=="nonadd") logfile.write("Particle interactions: non-additive");
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
    getline(inputAuxFile,line);
    istringstream(line)>>randomSeed;
    logfile.write("Random seed:",randomSeed);
    --logfile.currIndent;
    logfile.separator();



    return 0;
}
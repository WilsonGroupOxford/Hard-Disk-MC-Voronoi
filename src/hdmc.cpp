#include "hdmc.h"

//-------- HARD DISK MONTE CARLO CLASS --------


HDMC::HDMC() {
    //Default constructor

    n=0;
    phi=0;
}


int HDMC::setParticles(int num, double packFrac, int disp, VecF<double> dispParams) {
    //Set particle parameters

    //Set parameters
    n=num;
    phi=packFrac;

    //Allocate vectors
    x=VecF<double>(n);
    y=VecF<double>(n);
    r=VecF<double>(n);

    //Initialise particle lattice
    int lattice;
    if(disp==1) lattice=initMono(dispParams);

    return lattice;
}


int HDMC::initMono(VecF<double> dispParams) {
    //Initialise monodisperse particle configuration

    //Set radii and calculate cell parameters from packing fraction
    r=dispParams[0];
    double area=(n*M_PI*pow(dispParams[0],2))/phi;
    cellLen=sqrt(area);
    rCellLen=1.0/cellLen;
    cellLen_2=cellLen/2.0;

    //Set up regular square lattice
    int rows=floor(sqrt(n));
    if(rows*dispParams[0]*2>cellLen) return 1; //packing fraction too high to form lattice
    double spacing=cellLen/rows;
    for(int i=0; i<rows; ++i){
        y[i]=i*spacing;
        for(int j=0;j<rows; ++j){
            x[i]=j*spacing;
        }
    }

    return 0;
}


int HDMC::setRandom(int seed) {
    //Set random seed and generators

    mtGen.seed(seed);
    randParticle=uniform_int_distribution<int>(0,n-1);

    return 0;
}


int HDMC::setSimulation(int preEq, int eq, int prod, double swap, double accTarg) {
    //Set simulation parameters

    peqMoves=preEq;
    eqMoves=eq;
    prodMoves=prod;
    swapProb=swap;
    transProb=1.0-swapProb;
    acceptTarget=accTarg;

    return 0;
}
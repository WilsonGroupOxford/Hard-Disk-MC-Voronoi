#include "montecarlo.h"

//Default Constructor
MonteCarlo::MonteCarlo() {}

MonteCarlo::MonteCarlo(int numA, int numB, double radA, double radB, double cellLength, Logfile& logfile) {
    //Construct with particle type information

    //Assign parameters
    nA = numA;
    nB = numB;
    rA = radA;
    rB = radB;
    cellLen = cellLength;

    //Calculate additional parameters
    rAB = sqrt(rA*rB);
    hdAA = pow(2.0*rA,2);
    hdBB = pow(2.0*rA,2);
    hdAB = pow(2.0*rAB,2);
    minImageDist = cellLen/2.0;

    //Summarise in log file
    ++logfile.currIndent;
    logfile.write("Number of type a:",nA);
    logfile.write("Number of type b:",nB);
    logfile.write("Radius of type a:",rA);
    logfile.write("Radius of type b:",rB);
    logfile.write("Radius of type ab:",rAB);
    logfile.write("Cell length:",cellLength);
    --logfile.currIndent;
}

void MonteCarlo::loadCoordinates(string prefix, Logfile &logfile) {
    //Load x,y coordinates from xyz file

    //Set up coordinate arrays
    crdsA = VecF< VecF<double> >(nA);
    crdsB = VecF< VecF<double> >(nB);
    VecF<double> crd(2);

    //Read xyz file
    logfile.write("Reading coordinates");
    ifstream inputXYZFile(prefix+".xyz", ios::in);
    if(!inputXYZFile.good()) logfile.criticalError("Cannot find xyz file");
    string skip,line;
    istringstream ss("");
    getline(inputXYZFile,skip);
    getline(inputXYZFile,skip);
    for(int i=0; i<nA; ++i){
        getline(inputXYZFile,line);
        ss.str(line);
        ss>>skip;
        ss>>crd[0];
        ss>>crd[1];
        crdsA[i]=crd;
    }
    ++logfile.currIndent;
    logfile.write("Coordinates of type a read");
    for(int i=0; i<nB; ++i){
        getline(inputXYZFile,line);
        ss.str(line);
        ss>>skip;
        ss>>crd[0];
        ss>>crd[1];
        crdsB[i]=crd;
    }
    logfile.write("Coordinates of type b read");
    --logfile.currIndent;
}

void MonteCarlo::setParameters(int seed, int preEq, int eq, int prod, int writeFreq, int disp, int clst, Logfile& logfile) {
    //Set Monte Carlo variables

    //Initialise random number generators
    mtGen.seed(seed);
    randParticle=uniform_int_distribution<int>(0,nA+nB-1);
    randDisp=uniform_real_distribution<double>(-1.0,1.0);
    logfile.write("Initialising Monte Carlo process");
    ++logfile.currIndent;

    //Assign mc cycle information
    nCyclePreEq=preEq;
    nCycleEq=eq;
    nCycleProd=prod;
    cycleWriteFreq=writeFreq;
    cycleDispMoves=disp;
    cycleClstMoves=clst;

    //Summarise in log file
    logfile.write("Number of pre-equilibration cycles",nCyclePreEq);
    logfile.write("Number of equilibration cycles",nCycleEq);
    logfile.write("Number of production cycles",nCycleProd);
    logfile.write("Cycle write frequency",cycleWriteFreq);
    logfile.write("Displacement moves per cycle",cycleDispMoves);
    logfile.write("Cluster moves per cycle",cycleClstMoves);
    --logfile.currIndent;
}

void MonteCarlo::preEquilibration(Logfile &logfile) {
    //Perform Monte Carlo cycles, adjusting displacement delta to optimal

    monteCarloCycle();

}

double MonteCarlo::monteCarloCycle() {
    //Cycle of displacement and cluster moves

    //Displacement moves
    double pAccept = 0.0;
    for(int i=0; i<cycleDispMoves; ++i){
        pAccept += displacementMove();
    }

    return 0.0;
}

int MonteCarlo::displacementMove() {
    //Single Monte Carlo displacement move

    //Select random particle of type a or b
    int particleId = randParticle(mtGen);
    VecF<double> crdDisp(2),crdPrev(2),crdTrial(2);
    crdDisp[0] = randDisp(mtGen)*dispMoveDelta;
    crdDisp[1] = randDisp(mtGen)*dispMoveDelta;
    if(particleId<nA) {
        cout << "A " << particleId << endl;
        crdPrev = crdsA[particleId];
        crdTrial = crdPrev+crdDisp;

    }
    else{
        particleId -= nA;
        cout<<"B "<<particleId<<endl;
    }

    return 0;
}

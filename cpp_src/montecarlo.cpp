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
    rCellLen = 1.0/cellLen;
    cellLen_2 = cellLen/2.0;

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
    xA = VecF<double>(nA);
    yA = VecF<double>(nA);
    xB = VecF<double>(nB);
    yB = VecF<double>(nB);

    //Read xyz file
    logfile.write("Reading coordinates");
    ifstream inputXYZFile(prefix+"_init.xyz", ios::in);
    if(!inputXYZFile.good()) logfile.criticalError("Cannot find xyz file");
    string skip,line;
    istringstream ss("");
    getline(inputXYZFile,skip);
    getline(inputXYZFile,skip);
    for(int i=0; i<nA; ++i){
        getline(inputXYZFile,line);
        ss.str(line);
        ss>>skip;
        ss>>xA[i];
        ss>>yA[i];
    }
    ++logfile.currIndent;
    logfile.write("Coordinates of type a read");
    for(int i=0; i<nB; ++i){
        getline(inputXYZFile,line);
        ss.str(line);
        ss>>skip;
        ss>>xB[i];
        ss>>yB[i];
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
    dispMoveDelta=1.0;
    hdTol=-1e-12;

    //Summarise in log file
    logfile.write("Number of pre-equilibration cycles",nCyclePreEq);
    logfile.write("Number of equilibration cycles",nCycleEq);
    logfile.write("Number of production cycles",nCycleProd);
    logfile.write("Cycle write frequency",cycleWriteFreq);
    logfile.write("Displacement moves per cycle",cycleDispMoves);
    logfile.write("Cluster moves per cycle",cycleClstMoves);
    --logfile.currIndent;
}

void MonteCarlo::checkAllOverlaps(Logfile &logfile) {
    //Check all particle interactions for overlaps

    logfile.write("Checking all particle interactions");
    ++logfile.currIndent;
    bool overlap;
    for(int i=0; i<nA; ++i){
        overlap = nonAdditativeHardDiscOverlapAA(xA[i],yA[i],i);
        if(overlap) logfile.criticalError("Overlap detected A-A");
        overlap = nonAdditativeHardDiscOverlapAB(xA[i],yA[i]);
        if(overlap) logfile.criticalError("Overlap detected A-B");
    }
    for(int i=0; i<nB; ++i) {
        overlap = nonAdditativeHardDiscOverlapBB(xB[i], yB[i], i);
        if (overlap) logfile.criticalError("Overlap detected B-B");
    }
    logfile.write("All interactions satisfied");
    --logfile.currIndent;
}

void MonteCarlo::preEquilibration(Logfile &logfile) {
    //Perform Monte Carlo cycles, adjusting displacement delta to optimal

    //Use trial and improvement to find ideal displacement delta
    logfile.write("Pre-equilibration");
    ++logfile.currIndent;

    //Begin with smallest displacement to move particles from initial close packed positions
    int cycleCount=0;
    double pLow,dLow=rA/100.0;
    bool stage1=false;
    dispMoveDelta=dLow;
    for(;;){
        pLow=monteCarloCycle();
        ++cycleCount;
        if(pLow>0.9){
            stage1=true;
            break;
        }
        if(cycleCount==nCyclePreEq) break;
    }
    if(!stage1){
        logfile.write("Warning: pre-equilibration too short to determine ideal step size");
        logfile.write("Arbitrarily setting displacement delta to 1");
        dispMoveDelta=1.0;
        return;
    }

    //Next set largest displacement as half cell length
    double pUp,dUp=cellLen_2;
    dispMoveDelta=dUp;
    pUp=monteCarloCycle();
    ++cycleCount;
    bool stage2=cycleCount!=nCyclePreEq;
    if(!stage2){
        logfile.write("Warning: pre-equilibration too short to determine ideal step size");
        logfile.write("Arbitrarily setting displacement delta to 1");
        dispMoveDelta=1.0;
        return;
    }

    //Adjust delta each cycle by trial and improvement
    bool stage3;
    double p,d=0.5*(dLow+dUp);
    for(;;){
        p=monteCarloCycle();
        ++cycleCount;
        if(p>0.4) dLow=d;
        else dUp=d;
        d=0.5*(dLow+dUp);
        dispMoveDelta=d;
        if(abs(p-0.4)<0.01){
            stage3=true;
            break;
        }
        if(cycleCount==nCyclePreEq) break;
    }
    if(!stage3) logfile.write("Displacement delta set to: ",dispMoveDelta);

    //Second round of trial and improvement in pre-equilibrium steps remaining
    dLow=dispMoveDelta/2;
    dUp=dispMoveDelta*2;
    if(dUp>cellLen_2) dUp=cellLen_2;
    bool stage4;
    for(;;){
        p=monteCarloCycle();
        if(p>0.4) dLow=d;
        else dUp=d;
        d=0.5*(dLow+dUp);
        dispMoveDelta=d;
        if(abs(p-0.4)<0.001){
            stage4=true;
            break;
        }
        if(cycleCount==nCyclePreEq) break;
    }
    if(!stage4) logfile.write("Displacement delta set to: ",dispMoveDelta);
    else{
        logfile.write("Pre-equilibration finished early");
        logfile.write("Displacement delta set to: ",dispMoveDelta);
    }
    --logfile.currIndent;
}

void MonteCarlo::equilibration(Logfile &logfile) {
    //Perform Monte Carlo cycles without write out

    logfile.write("Equilibration");
    ++logfile.currIndent;
    double pAcc=0.0;
    for(int i=0; i<nCycleEq; ++i){
        pAcc=monteCarloCycle();
    }
    logfile.write("Equilibration acceptance probability: ",pAcc/nCycleEq);
    --logfile.currIndent;
}

void MonteCarlo::production(string prefix, Logfile &logfile) {
    //Perform Monte Carlo cycles with write out

    //Open xyz results file
    logfile.write("Production");
    ++logfile.currIndent;
    logfile.write("Opening .xyz file");
    ++logfile.currIndent;
    ofstream xyzFile(prefix + "_prod.xyz", ios::in | ios::trunc);
    xyzFile << fixed << showpoint << setprecision(8);
    int nC = nA+nB;

    //Write first configuration
    int nConfigs=0;
    xyzFile<<nC<<endl;
    xyzFile<<""<<endl;
    for(int i=0; i<nA; ++i) {
        xyzFile << setw(10) << left << "A";
        xyzFile << setw(20) << left << xA[i];
        xyzFile << setw(20) << left << yA[i];
        xyzFile << setw(20) << left << rA << endl;
    }
    for(int i=0; i<nB; ++i) {
        xyzFile << setw(10) << left << "B";
        xyzFile << setw(20) << left << xB[i];
        xyzFile << setw(20) << left << yB[i];
        xyzFile << setw(20) << left << rB << endl;
    }
    ++nConfigs;
    logfile.write("Configuration written: ",nConfigs);

    //Monte Carlo moves, writing at given frequency
    double pAcc=0.0;
    for(int i=1; i<=nCycleEq; ++i){
        pAcc=monteCarloCycle();
        if(i%cycleWriteFreq==0){
            xyzFile<<nC<<endl;
            xyzFile<<""<<endl;
            for(int i=0; i<nA; ++i) {
                xyzFile << setw(10) << left << "A";
                xyzFile << setw(20) << left << xA[i];
                xyzFile << setw(20) << left << yA[i];
                xyzFile << setw(20) << left << rA << endl;
            }
            for(int i=0; i<nB; ++i) {
                xyzFile << setw(10) << left << "B";
                xyzFile << setw(20) << left << xB[i];
                xyzFile << setw(20) << left << yB[i];
                xyzFile << setw(20) << left << rB << endl;
            }
            ++nConfigs;
            logfile.write("Configuration written: ",nConfigs);
        }
    }
    --logfile.currIndent;
    logfile.write("Production acceptance probability: ",pAcc/nCycleEq);

    //Close file
    xyzFile.close();
    --logfile.currIndent;
}

double MonteCarlo::monteCarloCycle() {
    //Cycle of displacement and cluster moves

    //Displacement moves
    double pAccept = 0.0;
    for(int i=0; i<cycleDispMoves; ++i){
        pAccept += displacementMove();
    }
    return pAccept/cycleDispMoves;
}

int MonteCarlo::displacementMove() {
    // Single Monte Carlo displacement move

    //Select random particle of type a or b
    int particleId = randParticle(mtGen); //random particle
    double dx,dy; //random displacement
    dx=randDisp(mtGen)*dispMoveDelta;
    dy=randDisp(mtGen)*dispMoveDelta;
    double xTrial,yTrial; //previous and trial coordinates
    bool reject;
    if(particleId<nA) {
        xTrial = xA[particleId]+dx;
        yTrial = yA[particleId]+dy;
        xTrial -= cellLen*floor(xTrial*rCellLen);
        yTrial -= cellLen*floor(yTrial*rCellLen);
        reject = nonAdditativeHardDiscOverlapAA(xTrial,yTrial,particleId);
        if(!reject) reject = nonAdditativeHardDiscOverlapAB(xTrial,yTrial);
        if(!reject){
            xA[particleId] = xTrial;
            yA[particleId] = yTrial;
        }
    }
    else{
        particleId -= nA;
        xTrial = xB[particleId]+dx;
        yTrial = yB[particleId]+dy;
        xTrial -= cellLen*floor(xTrial*rCellLen);
        yTrial -= cellLen*floor(yTrial*rCellLen);
        reject = nonAdditativeHardDiscOverlapBB(xTrial,yTrial,particleId);
        if(!reject) reject = nonAdditativeHardDiscOverlapBA(xTrial,yTrial);
        if(!reject){
            xB[particleId] = xTrial;
            yB[particleId] = yTrial;
        }
    }

    if(reject) return 0;
    else return 1;
}

bool MonteCarlo::nonAdditativeHardDiscOverlapAA(double& x, double& y, int &refId) {
    //Check for overlap of particle of type a with particles of type a

    double dx,dy;
    double dSq;
    for(int i=0; i<nA; ++i){
        dx = x-xA[i];
        dy = y-yA[i];
        dx -= cellLen*nearbyint(dx*rCellLen);
        dy -= cellLen*nearbyint(dy*rCellLen);
        dSq = dx*dx+dy*dy-hdAA;
        if(dSq<hdTol){
            if(i!=refId) return true;
        }
    }
    return false;
}

bool MonteCarlo::nonAdditativeHardDiscOverlapBB(double& x, double& y, int &refId) {
    //Check for overlap of particle of type b with particles of type b

    double dx,dy;
    double dSq;
    for(int i=0; i<nB; ++i){
        dx = x-xB[i];
        dy = y-yB[i];
        dx -= cellLen*nearbyint(dx*rCellLen);
        dy -= cellLen*nearbyint(dy*rCellLen);
        dSq = dx*dx+dy*dy-hdBB;
        if(dSq<hdTol){
            if(i!=refId) return true;
        }
    }
    return false;
}

bool MonteCarlo::nonAdditativeHardDiscOverlapAB(double& x, double& y) {
    //Check for overlap of particle of type a with particles of type b

    double dx,dy;
    double dSq;
    for(int i=0; i<nB; ++i){
        dx = x-xB[i];
        dy = y-yB[i];
        dx -= cellLen*nearbyint(dx*rCellLen);
        dy -= cellLen*nearbyint(dy*rCellLen);
        dSq = dx*dx+dy*dy-hdAB;
        if(dSq<hdTol){
            return true;
        }
    }
    return false;
}

bool MonteCarlo::nonAdditativeHardDiscOverlapBA(double& x, double& y) {
    //Check for overlap of particle of type b with particles of type a

    double dx,dy;
    double dSq;
    for(int i=0; i<nA; ++i){
        dx = x-xA[i];
        dy = y-yA[i];
        dx -= cellLen*nearbyint(dx*rCellLen);
        dy -= cellLen*nearbyint(dy*rCellLen);
        dSq = dx*dx+dy*dy-hdAB;
        if(dSq<hdTol){
            return true;
        }
    }
    return false;
}

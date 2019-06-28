#include "montecarlo.h"

//Default Constructor
MonteCarlo::MonteCarlo() {}

MonteCarlo::MonteCarlo(int numA, int numB, double radA, double radB, double cellLength, int nonAdditive, Logfile& logfile) {
    //Construct with particle type information

    //Assign parameters
    nA = numA;
    nB = numB;
    rA = radA;
    rB = radB;
    cellLen = cellLength;

    //Calculate additional parameters
    if(nonAdditive==1) rAB = sqrt(rA*rB);
    else rAB=0.5*(rA+rB);
    hdAA = pow(2.0*rA,2);
    hdBB = pow(2.0*rB,2);
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
    randClst=uniform_real_distribution<double>(-cellLen_2,cellLen_2);
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
        overlap = nonAdditiveHardDiscOverlapAA(xA[i], yA[i], i);
        if(overlap) logfile.criticalError("Overlap detected A-A");
        overlap = nonAdditiveHardDiscOverlapAB(xA[i], yA[i]);
        if(overlap) logfile.criticalError("Overlap detected A-B");
    }
    for(int i=0; i<nB; ++i) {
        overlap = nonAdditiveHardDiscOverlapBB(xB[i], yB[i], i);
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
        pLow=monteCarloCycle()[0];
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
    pUp=monteCarloCycle()[0];
    ++cycleCount;
    bool stage2=cycleCount!=nCyclePreEq;
    if(!stage2){
        logfile.write("Warning: pre-equilibration too short to determine ideal step size");
        logfile.write("Arbitrarily setting displacement delta to 1");
        dispMoveDelta=1.0;
        return;
    }

    if(pUp>0.4) {
        logfile.write("Too diffuse to determine maximum ideal move displacement");
        logfile.write("Setting displacement delta to L/2");
        dispMoveDelta=cellLen_2;
        return;
    }

    //Adjust delta each cycle by trial and improvement
    for(;;){
        double p;
        double d=0.5*(dLow+dUp);
        bool conv = false;
        dispMoveDelta=d;
        for(int i=0; i<100; ++i){
            p=monteCarloCycle()[0];
            ++cycleCount;
            if(abs(p-0.4)<0.01){
                conv=true;
                break;
            }
            else if(p<0.4) dUp=d;
            else dLow=d;
            d=0.5*(dLow+dUp);
            dispMoveDelta=d;
        }
        if(conv){
            dLow=0.5*d;
            dUp=1.5*d;
        }
        else{
            dLow=0.1*d;
            dUp=10.0*d;
        }
        if(dUp>cellLen_2) dUp=cellLen_2;
        if(cycleCount>nCyclePreEq*0.9) break; //leave cycles to find minimum again
    }
    logfile.write("Displacement delta set to: ",dispMoveDelta);
    --logfile.currIndent;
}

void MonteCarlo::equilibration(Logfile &logfile) {
    //Perform Monte Carlo cycles without write out

    logfile.write("Equilibration");
    ++logfile.currIndent;
    VecF<double> moveAnalysis(2);
    moveAnalysis = 0;
    for(int i=0; i<nCycleEq; ++i){
        moveAnalysis += monteCarloCycle();
    }
    moveAnalysis /= nCycleEq;
    logfile.write("Equilibration acceptance probability: ",moveAnalysis[0]);
    logfile.write("Equilibration average cluster size: ",moveAnalysis[1]);
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
    VecF<double> moveAnalysis(2);
    moveAnalysis = 0;
    for(int i=1; i<=nCycleProd; ++i){
        moveAnalysis += monteCarloCycle();
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
    moveAnalysis /= nCycleProd;
    --logfile.currIndent;
    logfile.write("Production acceptance probability: ",moveAnalysis[0]);
    logfile.write("Production average cluster size: ",moveAnalysis[1]);

    //Close file
    xyzFile.close();
    --logfile.currIndent;
}

VecF<double> MonteCarlo::monteCarloCycle() {
    //Cycle of displacement and cluster moves

    //Move analysis: displacement acceptance probability and cluster size
    VecF<double> moveAnalysis(2);
    moveAnalysis = 0;

    //Displacement moves
    for(int i=0; i<cycleDispMoves; ++i){
        moveAnalysis[0] += displacementMove();
    }

    //Cluster moves
    for(int i=0; i<cycleClstMoves; ++i){
        moveAnalysis[1] += clusterMove();
    }

    //Average move analysis
    moveAnalysis[0] /= cycleDispMoves;
    moveAnalysis[1] /= cycleClstMoves;

    return moveAnalysis;
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
        xTrial -= cellLen*nearbyint(xTrial*rCellLen);
        yTrial -= cellLen*nearbyint(yTrial*rCellLen);
        reject = nonAdditiveHardDiscOverlapAA(xTrial, yTrial, particleId);
        if(!reject) reject = nonAdditiveHardDiscOverlapAB(xTrial, yTrial);
        if(!reject){
            xA[particleId] = xTrial;
            yA[particleId] = yTrial;
        }
    }
    else{
        particleId -= nA;
        xTrial = xB[particleId]+dx;
        yTrial = yB[particleId]+dy;
        xTrial -= cellLen*nearbyint(xTrial*rCellLen);
        yTrial -= cellLen*nearbyint(yTrial*rCellLen);
        reject = nonAdditiveHardDiscOverlapBB(xTrial, yTrial, particleId);
        if(!reject) reject = nonAdditiveHardDiscOverlapBA(xTrial, yTrial);
        if(!reject){
            xB[particleId] = xTrial;
            yB[particleId] = yTrial;
        }
    }

    if(reject) return 0;
    else return 1;
}

int MonteCarlo::clusterMove() {
    //Single Monte Carlo cluster move - invert coordinates in random pivot

    //Vectors containing information on which particles to invert
    VecF<int> invertA(nA),invertB(nB); //0=leave,1=invert,2=inverted
    VecF<int> overlapA(nA),overlapB(nB);
    invertA = 0;
    invertB = 0;
    overlapA = 0;
    overlapB = 0;

    //Select random particle of type a or b and pivot point
    int particleId = randParticle(mtGen);
    if(particleId<nA) invertA[particleId]=1;
    else invertB[particleId-nA]=1;
    double xPivot = randClst(mtGen);
    double yPivot = randClst(mtGen);

    //Perform cluster move: invert particles, find overlapping and loop until no overlaps remain
    int nInverted = 0;
    double dx,dy;
    for(;;){
        //A: invert particles and find overlapping
        for(int i=0; i<nA; ++i){
            if(invertA[i]==1){
                dx = xA[i] - xPivot;
                dy = yA[i] - yPivot;
                xA[i] = xPivot - dx;
                yA[i] = yPivot - dy;
                xA[i] -= cellLen*nearbyint(xA[i]*rCellLen);
                yA[i] -= cellLen*nearbyint(yA[i]*rCellLen);
                nonAdditiveHardDiscOverlapA(xA[i],yA[i],overlapA,overlapB);
                ++nInverted;
                invertA[i]=2;
            }
        }
        //B: invert particles and find overlapping
        for(int i=0; i<nB; ++i){
            if(invertB[i]==1){
                dx = xB[i] - xPivot;
                dy = yB[i] - yPivot;
                xB[i] = xPivot - dx;
                yB[i] = yPivot - dy;
                xB[i] -= cellLen*nearbyint(xB[i]*rCellLen);
                yB[i] -= cellLen*nearbyint(yB[i]*rCellLen);
                nonAdditiveHardDiscOverlapB(xB[i],yB[i],overlapA,overlapB);
                ++nInverted;
                invertB[i]=2;
            }
        }

        //Update inversion list
        bool complete = true;
        for(int i=0; i<nA; ++i){
            if(overlapA[i] && invertA[i]==0){
                invertA[i] = 1;
                complete = false;
            }
        }
        for(int i=0; i<nB; ++i){
            if(overlapB[i] && invertB[i]==0){
                invertB[i] = 1;
                complete = false;
            }
        }
        overlapA=0;
        overlapB=0;

        //Complete when no more overlaps
        if(complete) break;
    }

    return nInverted;
}

bool MonteCarlo::nonAdditiveHardDiscOverlapAA(double &x, double &y, int &refId) {
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

bool MonteCarlo::nonAdditiveHardDiscOverlapBB(double &x, double &y, int &refId) {
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

bool MonteCarlo::nonAdditiveHardDiscOverlapAB(double &x, double &y) {
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

bool MonteCarlo::nonAdditiveHardDiscOverlapBA(double &x, double &y) {
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

void MonteCarlo::nonAdditiveHardDiscOverlapA(double &x, double &y, VecF<int> &overlapA, VecF<int> &overlapB) {
    //Check for overlap of particle of type a with all others, storing overlap ids

    //AA
    double dx,dy;
    double dSq;
    for(int i=0; i<nA; ++i){
        dx = x-xA[i];
        dy = y-yA[i];
        dx -= cellLen*nearbyint(dx*rCellLen);
        dy -= cellLen*nearbyint(dy*rCellLen);
        dSq = dx*dx+dy*dy-hdAA;
        if(dSq<hdTol) overlapA[i]=1;
    }

    //AB
    for(int i=0; i<nB; ++i) {
        dx = x - xB[i];
        dy = y - yB[i];
        dx -= cellLen * nearbyint(dx * rCellLen);
        dy -= cellLen * nearbyint(dy * rCellLen);
        dSq = dx * dx + dy * dy - hdAB;
        if (dSq < hdTol) overlapB[i]=1;
    }
}

void MonteCarlo::nonAdditiveHardDiscOverlapB(double &x, double &y, VecF<int> &overlapA, VecF<int> &overlapB) {
    //Check for overlap of particle of type b with all others, storing overlap ids

    //BB
    double dx,dy;
    double dSq;
    for(int i=0; i<nB; ++i){
        dx = x-xB[i];
        dy = y-yB[i];
        dx -= cellLen*nearbyint(dx*rCellLen);
        dy -= cellLen*nearbyint(dy*rCellLen);
        dSq = dx*dx+dy*dy-hdBB;
        if(dSq<hdTol) overlapB[i]=1;
    }

    //AB
    for(int i=0; i<nA; ++i) {
        dx = x - xA[i];
        dy = y - yA[i];
        dx -= cellLen * nearbyint(dx * rCellLen);
        dy -= cellLen * nearbyint(dy * rCellLen);
        dSq = dx * dx + dy * dy - hdAB;
        if (dSq < hdTol) overlapA[i]=1;
    }
}

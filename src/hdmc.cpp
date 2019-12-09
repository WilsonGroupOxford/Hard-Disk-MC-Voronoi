#include "hdmc.h"

//-------- HARD DISK MONTE CARLO CLASS --------


//-------- CONSTRUCTORS, SETTERS --------


HDMC::HDMC() {
    //Default constructor

    n=0;
    nA=0;
    nB=0;
    phi=0;
}


int HDMC::setParticles(int num, double packFrac, int disp, VecF<double> dispParams, int interact) {
    //Set particle parameters

    //Set parameters
    n=num;
    interaction=interact;
    phi=packFrac;
    dispersity=disp;
    dispersityParams=dispParams;
    if(dispersity==1){//monodisperse neglect particle types
        nA=n;
        nB=0;
    }
    else if(dispersity==2){//bidisperse utilise particle types
        nA=nearbyint(n*dispParams[2]);
        nB=n-nA;
    }
    else if(dispersity==3){//polydisperse neglect particle types
        nA=n;
        nB=0;
    }

    return 0;
}


int HDMC::setRandom(int seed) {
    //Set random seed and generators

    mtGen.seed(seed);
    randParticle=uniform_int_distribution<int>(0,n-1);
    rand01=uniform_real_distribution<double>(0,1);

    return 0;
}


int HDMC::setSimulation(int eq, int prod, double swap, double accTarg) {
    //Set simulation parameters

    eqCycles=eq;
    prodCycles=prod;
    swapProb=swap;
    transProb=1.0-swapProb;
    acceptTarget=accTarg;
    transDelta=1.0;

    return 0;
}


int HDMC::setAnalysis(string path, int anFreq, int rdf, double rdfDel, VecF<int> vor, double radCut, int visF, int vis3) {
    //Set analysis parameters

    outputPrefix=path;
    analysisFreq=anFreq;

    //Set rdf type
    if(rdf==0) rdfCalc=false;
    else if(rdf==1){
        rdfCalc=true;
        rdfNorm=true;
        rdfDelta=rdfDel;
    }
    else if(rdf==2){
        rdfCalc=true;
        rdfNorm= false;
        rdfDelta=rdfDel;
    }

    //Set visualisation output frequency with 0 preventing write
    visFreq=visF;
    if(visFreq==0){
        visXYZ=false;
        visVor2D=false;
        visVor3D=false;
        visFreq=1;
    }
    else{
        visXYZ=true;
        visVor2D=true;
        visVor3D=vis3;
    }

    //Set voronoi type
    vorCalc2D=false;
    radCalc2D=false;
    vorCalc3D=false;
    radCalc3D=false;
    if(vor[0]==1) vorCalc2D=true;
    if(vor[1]==1) radCalc2D=true;
    if(vor[2]==1) vorCalc3D=true;
    if(vor[3]==1) radCalc3D=true;
    rad2DCut=radCut;
    maxVertices=21;

    return 0;
}


int HDMC::initAnalysis() {
    //Initialise analysis tools

    analysisConfigs=0;
    xyzConfigs=0;

    //RDF histogram
    if(rdfCalc){
        int maxBin=floor(cellLen_2/rdfDelta)+1; //max distance is half cell size
        rdfHist=VecF<int>(maxBin);
        if(dispersity==2){//bidisperse calculate partial rdfs
            prdfHistAA=VecF<int>(maxBin);
            prdfHistAB=VecF<int>(maxBin);
            prdfHistBB=VecF<int>(maxBin);
        }
    }

    //Voronoi distributions
    if(vorCalc2D){
        vor2DSizesA=VecF<int>(maxVertices);
        vor2DSizesB=VecF<int>(maxVertices);
        vor2DAdjs=VecF< VecF<int> >(maxVertices);
        for(int i=0; i<maxVertices; ++i) vor2DAdjs[i]=VecF<int>(maxVertices);
        vor2DAreasA=VecF<double>(maxVertices+1);
        vor2DAreasB=VecF<double>(maxVertices+1);
        vor2DNNCount=VecF<int>(3);
        vor2DNNSep=VecF<double>(3);
    }
    if(radCalc2D){
        rad2DSizesA=VecF<int>(maxVertices);
        rad2DSizesB=VecF<int>(maxVertices);
        rad2DAdjs=VecF< VecF<int> >(maxVertices);
        for(int i=0; i<maxVertices; ++i) rad2DAdjs[i]=VecF<int>(maxVertices);
        rad2DAreasA=VecF<double>(maxVertices+1);
        rad2DAreasB=VecF<double>(maxVertices+1);
        rad2DNNCount=VecF<int>(3);
        rad2DNNSep=VecF<double>(3);
    }
    if(vorCalc3D){
        vor3DSizesA=VecF<int>(maxVertices);
        vor3DSizesB=VecF<int>(maxVertices);
        vor3DAdjs=VecF< VecF<int> >(maxVertices);
        for(int i=0; i<maxVertices; ++i) vor3DAdjs[i]=VecF<int>(maxVertices);
        vor3DAreasA=VecF<double>(maxVertices+1);
        vor3DAreasB=VecF<double>(maxVertices+1);
        vor3DNNCount=VecF<int>(3);
        vor3DNNSep=VecF<double>(3);
    }
    if(radCalc3D){
        rad3DSizesA=VecF<int>(maxVertices);
        rad3DSizesB=VecF<int>(maxVertices);
        rad3DAdjs=VecF< VecF<int> >(maxVertices);
        for(int i=0; i<maxVertices; ++i) rad3DAdjs[i]=VecF<int>(maxVertices);
        rad3DAreasA=VecF<double>(maxVertices+1);
        rad3DAreasB=VecF<double>(maxVertices+1);
        rad3DNNCount=VecF<int>(3);
        rad3DNNSep=VecF<double>(3);
    }

    return 0;
}


//---------- INITIAL CONFIGURATION --------


int HDMC::initialiseConfiguration(Logfile &logfile, double maxIt) {
    //Generate initial particle positions

    logfile.write("Generating Initial Configuration");
    cout<<"Generating Initial Configuration"<<endl;
    ++logfile.currIndent;

    //Allocate vectors
    x=VecF<double>(n);
    y=VecF<double>(n);
    z=VecF<double>(n);
    r=VecF<double>(n);
    w=VecF<double>(n);

    //Calculate simulation cell parameters
    double area;
    if(dispersity==1) area=(M_PI*n*pow(dispersityParams[0],2))/phi;
    else if(dispersity==2) area=M_PI*(nA*pow(dispersityParams[0],2)+nB*pow(dispersityParams[1],2))/phi;
    else if(dispersity==3) area=M_PI*(nA*pow(dispersityParams[0],2))/phi;
    cellLen=sqrt(area);
    rCellLen=1.0/cellLen;
    cellLen_2=cellLen/2.0;

    //Generate particle radii
    if(dispersity==1){
        r=dispersityParams[0];
    }
    else if(dispersity==2){
        double rA,rB,wA,wB;
        rA=dispersityParams[0];
        rB=dispersityParams[1];
        for(int i=0; i<nA; ++i){
            r[i]=rA;
        }
        for(int i=nA; i<n; ++i){
            r[i]=rB;
        }
    }
    else if(dispersity==3){
        //Randomly generate radii from gaussian
        normal_distribution<double> normalDist(dispersityParams[0],dispersityParams[1]);
        for(int i=0; i<n; ++i){
            double randR;
            for(;;){
                randR=normalDist(mtGen);
                if(randR>0) break;
            }
            r[i]=randR;
        }
    }

    //Generate weights for 2D Radical
    rad2DInclude=VecF<bool>(n);
    for(int i=0; i<r.n; ++i){
        if(2*r[i]<rad2DCut){
            w[i]=0.0;
            rad2DInclude[i]=false;
        }
        else{
            w[i]=pow(2*rad2DCut*r[i]-rad2DCut*rad2DCut,0.5);
            rad2DInclude[i]=true;
        }
    }
    cout<<"+++ "<<vSum(rad2DInclude)<<endl;


    //Generate initial configuration
    bool success;
    int attempt=1;
    for(;;){
        success=rsaPositions(maxIt);
        logfile.write("Attempt "+to_string(attempt)+" successful:",success);
        cout<<"Attempt "+to_string(attempt)+" successful: "<<success<<endl;
        if(success) break;
        ++attempt;
    }

    //Set z coordinate
    //z only used for analysis not Monte Carlo
    if(interaction==0) z=0.0;
    else if(interaction==1) z=r;

    //Exit if cannot generate
    if(!success) logfile.criticalError("Could not generate starting configuration");

    //Initialise analysis tools here as require cell length
    initAnalysis();

    --logfile.currIndent;
    logfile.separator();

    return !success;
}


bool HDMC::rsaPositions(double maxIt) {
    //Generate random particle positions using Random Sequential Adsorption algorithm

    //Get sort of radii positions - largest to smallest
    bool success=true;
    VecF<int> sort=vArgSort(r,true);

    //Add random particles without overlap
    int maxIterations=pow(n,maxIt);
    randomPosition(x[sort[0]],y[sort[0]]);
    int added=1,ii,jj,iterations=0;
    double xx,yy,rr,dx,dy,dSq,rSq;
    bool accept;
    if(interaction==0){
        while(added<n){
            ii=sort[added];
            randomPosition(xx,yy);
            rr=r[ii];
            accept=true;
            for(int j=0; j<added; ++j){
                jj=sort[j];
                dx=xx-x[jj];
                dy=yy-y[jj];
                dx-=cellLen*nearbyint(dx*rCellLen);
                dy-=cellLen*nearbyint(dy*rCellLen);
                dSq=dx*dx+dy*dy;
                rSq=pow((rr+r[jj]),2);
                if(dSq<rSq){
                    accept=false;
                    break;
                }
            }
            if(accept){
                x[ii]=xx;
                y[ii]=yy;
                ++added;
                cout<<"Particles placed: "<<added<<endl;
            }
            ++iterations;
            if(iterations==maxIterations){
                success=false;
                break;
            }
        }
    }
    else if(interaction==1){
        while(added<n){
            ii=sort[added];
            randomPosition(xx,yy);
            rr=r[ii];
            accept=true;
            for(int j=0; j<added; ++j){
                jj=sort[j];
                dx=xx-x[jj];
                dy=yy-y[jj];
                dx-=cellLen*nearbyint(dx*rCellLen);
                dy-=cellLen*nearbyint(dy*rCellLen);
                dSq=dx*dx+dy*dy;
                rSq=4*rr*r[jj];
                if(dSq<rSq){
                    accept=false;
                    break;
                }
            }
            if(accept){
                x[ii]=xx;
                y[ii]=yy;
                ++added;
                cout<<"Particles placed: "<<added<<endl;
            }
            ++iterations;
            if(iterations==maxIterations){
                success=false;
                break;
            }
        }
    }

    return success;
}


inline void HDMC::randomPosition(double &xx, double &yy) {
    //Generate random particle position inside periodic box

    xx=rand01(mtGen)*cellLen;
    yy=rand01(mtGen)*cellLen;
    xx-=cellLen*nearbyint(xx*rCellLen);
    yy-=cellLen*nearbyint(yy*rCellLen);
}


void HDMC::generateRandomPositions() {
    //Generate random particle positions inside periodic box

    //Generate random overlaps
    for(int i=0; i<n; ++i){
        x[i]=rand01(mtGen)*cellLen;
        y[i]=rand01(mtGen)*cellLen;
        x[i]-=cellLen*nearbyint(x[i]*rCellLen);
        y[i]-=cellLen*nearbyint(y[i]*rCellLen);
    }
}


bool HDMC::resolvePositions() {
    //Resolve particle overlaps using steepest descent minimisation with LJ repulsive particles

    //Set up coordinate vector
    VecF<double> xy(2*n);
    for(int i=0; i<n; ++i){
        xy[2*i]=x[i];
        xy[2*i+1]=y[i];
    }

    //Generate repulsive pairs
    int nReps=n*(n-1)/2;
    VecF<int> repPairs(2*nReps);
    int repCount=0;
    for(int i=0; i<n-1; ++i){
        for(int j=i+1; j<n; ++j){
            repPairs[2*repCount]=i;
            repPairs[2*repCount+1]=j;
            ++repCount;
        }
    }

    //Generate repulsive epsilon parameters
    VecF<double> repParams(2*nReps);
    for(int i=0, j=1; i<2*nReps; i+=2, j+=2) repParams[j]=0.001;

    //Set up potential model and optimiser
    HLJ2DP potModel(cellLen, cellLen);
    SteepestDescentArmijoMultiDim<HLJ2DP> optimiser(1000,0.5,1e-12);

    //Increment radii and minimise iteratively
    int swellSteps=100; //number of steps to swell particles
    double swellFactor=1.0/swellSteps; //amount to swell particles by each step
    double overSwell=0.01; //amount to over-swell by
    OutputFile xyzFile("swell.xyz");
    for(int k=0; k<=swellSteps; ++k){
        cout<<"Swelling particles step: "<<k<<endl;
        for(int i=0,j=1; i<2*nReps; i+=2, j+=2) repParams[i]=pow((k*swellFactor+overSwell)*(r[repPairs[i]]+r[repPairs[j]]),2);
        potModel.setRepulsions(repPairs,repParams);
        optimiser(potModel,xy);
        writeXYZ(xyzFile);
    }

    //Update coordinates
    for(int i=0; i<n; ++i){
        x[i]=xy[2*i];
        y[i]=xy[2*i+1];
    }

    //Check overlaps have been resolved
    double dx,dy,dSq,rSq;
    bool resolved=true;
    for(int i=0; i<n-1; ++i){
        for(int j=i+1; j<n; ++j){
            dx=x[i]-x[j];
            dy=y[i]-y[j];
            dx-=cellLen*nearbyint(dx*rCellLen);
            dy-=cellLen*nearbyint(dy*rCellLen);
            dSq=dx*dx+dy*dy;
            rSq=pow((r[i]+r[j]),2);
            if(dSq<rSq){
                resolved=false;
                break;
            }
        }
    }

    return resolved;
}


//-------- MONTE CARLO MOVES --------


inline int HDMC::mcCycle() {
    //Cycle of n-particle Monte Carlo moves

    int accCount=0;
    if(interaction==0){
        for(int i=0; i<n; ++i) mcAdditiveMove(accCount);
    }
    else if(interaction==1){
        for(int i=0; i<n; ++i) mcNonAdditiveMove(accCount);
    }

    return accCount;
}


inline void HDMC::mcAdditiveMove(int &counter) {
    //Single Monte Carlo move

    //Choose random particle and get position and radius
    int pI=randParticle(mtGen);
    double xI=x[pI];
    double yI=y[pI];
    double rI=r[pI];

    //Perform move
    if(rand01(mtGen)<transProb){
        //Translation move

        //Apply translation
        xI+=transDelta*(2*rand01(mtGen)-1);
        yI+=transDelta*(2*rand01(mtGen)-1);
        xI-=cellLen*nearbyint(xI*rCellLen);
        yI-=cellLen*nearbyint(yI*rCellLen);

        //Check for overlap with other particles
        double dx,dy,dSq,rSq;
        bool accept=true;
        for(int i=0; i<n; ++i){
            dx=xI-x[i];
            dy=yI-y[i];
            dx-=cellLen*nearbyint(dx*rCellLen);
            dy-=cellLen*nearbyint(dy*rCellLen);
            dSq=dx*dx+dy*dy;
            rSq=pow((rI+r[i]),2);
            if(dSq<rSq && i!=pI){
                accept=false;
                break;
            }
        }

        if(accept){
            x[pI]=xI;
            y[pI]=yI;
            ++counter;
        }
    }
    else{
        //Swap move

        //Choose second random particle
        int pJ=pI;
        while(pI==pJ) pJ=randParticle(mtGen);

        //Swap coordinates and radii
        double xJ=xI;
        double yJ=yI;
        double rJ=r[pJ];
        xI=x[pJ];
        yI=y[pJ];

        //Apply translations
        xI+=transDelta*(2*rand01(mtGen)-1);
        yI+=transDelta*(2*rand01(mtGen)-1);
        xI-=cellLen*nearbyint(xI*rCellLen);
        yI-=cellLen*nearbyint(yI*rCellLen);
        xJ+=transDelta*(2*rand01(mtGen)-1);
        yJ+=transDelta*(2*rand01(mtGen)-1);
        xJ-=cellLen*nearbyint(xJ*rCellLen);
        yJ-=cellLen*nearbyint(yJ*rCellLen);

        //Check for overlap with other particles
        double dx,dy,dSq,rSq;
        bool accept=true;
        dx=xI-xJ;
        dy=yI-yJ;
        dx-=cellLen*nearbyint(dx*rCellLen);
        dy-=cellLen*nearbyint(dy*rCellLen);
        dSq=dx*dx+dy*dy;
        rSq=pow((rI+rJ),2);
        if(dSq<rSq) accept=false;
        if(accept){
            for(int i=0; i<n; ++i){
                dx=xI-x[i];
                dy=yI-y[i];
                dx-=cellLen*nearbyint(dx*rCellLen);
                dy-=cellLen*nearbyint(dy*rCellLen);
                dSq=dx*dx+dy*dy;
                rSq=pow((rI+r[i]),2);
                if(dSq<rSq && i!=pI && i!=pJ){
                    accept=false;
                    break;
                }
            }
        }
        if(accept){
            for(int i=0; i<n; ++i){
                dx=xJ-x[i];
                dy=yJ-y[i];
                dx-=cellLen*nearbyint(dx*rCellLen);
                dy-=cellLen*nearbyint(dy*rCellLen);
                dSq=dx*dx+dy*dy;
                rSq=pow((rJ+r[i]),2);
                if(dSq<rSq && i!=pI && i!=pJ){
                    accept=false;
                    break;
                }
            }
        }

        if(accept){
            x[pI]=xI;
            y[pI]=yI;
            x[pJ]=xJ;
            y[pJ]=yJ;
            ++counter;
        }
    }
}


inline void HDMC::mcNonAdditiveMove(int &counter) {
    //Single Monte Carlo move

    //Choose random particle and get position and radius
    int pI=randParticle(mtGen);
    double xI=x[pI];
    double yI=y[pI];
    double rI=r[pI];

    //Perform move
    if(rand01(mtGen)<transProb){
        //Translation move

        //Apply translation
        xI+=transDelta*(2*rand01(mtGen)-1);
        yI+=transDelta*(2*rand01(mtGen)-1);
        xI-=cellLen*nearbyint(xI*rCellLen);
        yI-=cellLen*nearbyint(yI*rCellLen);

        //Check for overlap with other particles
        double dx,dy,dSq,rSq;
        bool accept=true;
        for(int i=0; i<n; ++i){
            dx=xI-x[i];
            dy=yI-y[i];
            dx-=cellLen*nearbyint(dx*rCellLen);
            dy-=cellLen*nearbyint(dy*rCellLen);
            dSq=dx*dx+dy*dy;
            rSq=4*rI*r[i];
            if(dSq<rSq && i!=pI){
                accept=false;
                break;
            }
        }

        if(accept){
            x[pI]=xI;
            y[pI]=yI;
            ++counter;
        }
    }
    else{
        //Swap move

        //Choose second random particle
        int pJ=pI;
        while(pI==pJ) pJ=randParticle(mtGen);

        //Swap coordinates and radii
        double xJ=xI;
        double yJ=yI;
        double rJ=r[pJ];
        xI=x[pJ];
        yI=y[pJ];

        //Apply translations
        xI+=transDelta*(2*rand01(mtGen)-1);
        yI+=transDelta*(2*rand01(mtGen)-1);
        xI-=cellLen*nearbyint(xI*rCellLen);
        yI-=cellLen*nearbyint(yI*rCellLen);
        xJ+=transDelta*(2*rand01(mtGen)-1);
        yJ+=transDelta*(2*rand01(mtGen)-1);
        xJ-=cellLen*nearbyint(xJ*rCellLen);
        yJ-=cellLen*nearbyint(yJ*rCellLen);

        //Check for overlap with other particles
        double dx,dy,dSq,rSq;
        bool accept=true;
        dx=xI-xJ;
        dy=yI-yJ;
        dx-=cellLen*nearbyint(dx*rCellLen);
        dy-=cellLen*nearbyint(dy*rCellLen);
        dSq=dx*dx+dy*dy;
        rSq=4*rI*rJ;
        if(dSq<rSq) accept=false;
        if(accept){
            for(int i=0; i<n; ++i){
                dx=xI-x[i];
                dy=yI-y[i];
                dx-=cellLen*nearbyint(dx*rCellLen);
                dy-=cellLen*nearbyint(dy*rCellLen);
                dSq=dx*dx+dy*dy;
                rSq=4*rI*r[i];
                if(dSq<rSq && i!=pI && i!=pJ){
                    accept=false;
                    break;
                }
            }
        }
        if(accept){
            for(int i=0; i<n; ++i){
                dx=xJ-x[i];
                dy=yJ-y[i];
                dx-=cellLen*nearbyint(dx*rCellLen);
                dy-=cellLen*nearbyint(dy*rCellLen);
                dSq=dx*dx+dy*dy;
                rSq=4*rJ*r[i];
                if(dSq<rSq && i!=pI && i!=pJ){
                    accept=false;
                    break;
                }
            }
        }

        if(accept){
            x[pI]=xI;
            y[pI]=yI;
            x[pJ]=xJ;
            y[pJ]=yJ;
            ++counter;
        }
    }
}


//-------- MONTE CARLO SIMULATION --------


void HDMC::equilibration(Logfile &logfile, OutputFile &xyzFile) {
    //Equilibration Monte Carlo

    //Header
    logfile.write("Equilibration Monte Carlo");
    cout<<"Equilibration"<<endl;
    ++logfile.currIndent;

    //Determine ideal translation delta
    logfile.write("Finding optimal displacement delta for acceptance probability:",acceptTarget);
    ++logfile.currIndent;
    //First some loops to remove any initial ordering
    logfile.write("Disrupting any initial ordering");
    for(int i=0; i<100; ++i){
        double deltaMin=0.01*vMinimum(r);
        double deltaMax=cellLen_2;
        double accProb;
        optimalDelta(deltaMin,deltaMax,accProb);
    }
    //Loop until converged
    bool converged=false;
    int optCode;
    double deltaMin=0.01*vMinimum(r);
    double deltaMax=cellLen_2;
    double accProb;
    int iteration=0;
    while(!converged){
        optCode=optimalDelta(deltaMin,deltaMax,accProb);
        if(optCode==1 && iteration==0){
            logfile.write("System too dense to achieve target");
            break;
        }
        else if(optCode==2 && iteration==0){
            logfile.write("System too dilute to achieve target");
            break;
        }
        else logfile.write("Delta: "+to_string(transDelta)+" acceptance: "+to_string(accProb));
        if(abs(accProb-acceptTarget)<0.005) break;
        if(iteration>100){
            logfile.write("Iteration limit hit");
            break;
        }
        ++iteration;
    }
    --logfile.currIndent;
    logfile.write("Translation delta set to:",transDelta);

    //Equilibration
    logfile.write("Running equilibration");
    ++logfile.currIndent;
    int logMoves=eqCycles/100;
    int accCount=0;
    for (int i = 1; i<=eqCycles; ++i) {
        accCount+=mcCycle();
        if(i%logMoves==0){
            logfile.write("Move cycles and acceptance:",i,double(accCount)/(i*n));
            cout<<"Move cycles and acceptance: "<<i<<" "<<double(accCount)/(i*n)<<endl;
        }
    }
    logfile.currIndent-=2;
    logfile.separator();
}


int HDMC::optimalDelta(double &deltaMin, double &deltaMax, double &accProb) {
    //Find optimal translation delta by trial and improvement

    //Generate trial delta values
    double logDeltaMin=log10(deltaMin);
    double logDeltaMax=log10(deltaMax);
    VecF<double> trialDelta(11),trialProb(11);
    for(int i=0; i<11; ++i) trialDelta[i]=pow(10,logDeltaMin+i*(logDeltaMax-logDeltaMin)/10.0);

    //Calculate acceptance probabilities for trial delta values
    for(int i=0; i<11; ++i){
        transDelta=trialDelta[i];
        int accCount=0;
        for(int j=0; j<10; ++j) accCount+=mcCycle();
        trialProb[i]=double(accCount)/(10*n);
    }

    //Find limiting delta values which surround target acceptance
    int optCode;
    if(trialProb[0]<acceptTarget){
        //lower bound too low
        transDelta=trialDelta[0];
        optCode=1;
    }
    if(trialProb[10]>acceptTarget){
        //upper bound too high
        transDelta=trialDelta[10];
        optCode=2;
    }
    else{
        for(int i=0; i<11; ++i){
            if(trialProb[i]>acceptTarget) deltaMin=trialDelta[i];
            else if(trialProb[i]<acceptTarget){
                deltaMax=trialDelta[i];
                break;
            }
        }
        transDelta=pow(10,0.5*(log10(deltaMin)+log10(deltaMax)));
        optCode=0;
    }

    //Calculate best guess for delta from current iteration
    int accCount=0;
    for(int j=0; j<10; ++j) accCount+=mcCycle();
    accProb=double(accCount)/(10*n);

    return optCode;
}


void HDMC::production(Logfile &logfile, OutputFile &xyzFile, OutputFile &vor2DFile, OutputFile &rad2DFile,
                      OutputFile &vor3DFile, OutputFile &rad3DFile, OutputFile &vis2DFile, OutputFile &vis3DFile) {
    //Production Monte Carlo

    //Production cycles
    logfile.write("Production Monte Carlo");
    cout<<"Production"<<endl;
    ++logfile.currIndent;
    int logMoves=prodCycles/100;
    int accCount=0;
    for (int i = 1; i<=prodCycles; ++i) {
        accCount+=mcCycle();
        if(i%logMoves==0){
            logfile.write("Move cycles and acceptance:",i,double(accCount)/(i*n));
            cout<<"Move cycles and acceptance: "<<i<<" "<<double(accCount)/(i*n)<<endl;
        }
        if(i%analysisFreq==0){
            bool vis=(i%visFreq==0)*visVor2D;
            if(vis && visXYZ) writeXYZ(xyzFile);
            analyseConfiguration(vor2DFile,rad2DFile,vor3DFile,rad3DFile,vis2DFile,vis3DFile,vis);
        }
    }
    logfile.currIndent-=2;
    logfile.separator();
}


//--------- ANALYSIS ----------


void HDMC::analyseConfiguration(OutputFile &vor2DFile, OutputFile &rad2DFile, OutputFile &vor3DFile,
                                OutputFile &rad3DFile, OutputFile &vis2DFile, OutputFile &vis3DFile, bool vis) {
    //Control analysis of current configuration

    if(rdfCalc) calculateRDF();
    if(vorCalc2D) calculateVoronoi2D(vor2DFile,vis2DFile,vis);
    if(radCalc2D) calculateRadical2D(rad2DFile,vis2DFile,vis);
    if(vorCalc3D) calculateVoronoi3D(vor3DFile,vis2DFile,vis3DFile,vis);
    if(radCalc3D) calculateRadical3D(rad3DFile,vis2DFile,vis3DFile,vis);

    ++analysisConfigs;
}


void HDMC::calculateRDF() {
    //Calculate RDF for current configuration

    //Calculate pairwise distances and bin
    double xI,yI,b;
    double dx,dy,dSq,d;
    if(dispersity==1) {//monodisperse calculate total rdf
        for (int i = 0; i < n - 1; ++i) {
            xI = x[i];
            yI = y[i];
            for (int j = i + 1; j < n; ++j) {
                dx = xI - x[j];
                dy = yI - y[j];
                dx -= cellLen * nearbyint(dx * rCellLen);
                dy -= cellLen * nearbyint(dy * rCellLen);
                dSq = dx * dx + dy * dy;
                d = sqrt(dSq);
                if (d < cellLen_2) {
                    b = floor(d / rdfDelta);
                    rdfHist[b] += 2;
                }
            }
        }
    }
    else if(dispersity==2) {//bidisperse calculate partial then total rdfs
        //AA
        for (int i=0; i<nA-1; ++i) {
            xI = x[i];
            yI = y[i];
            for (int j=i+1; j<nA; ++j) {
                dx = xI - x[j];
                dy = yI - y[j];
                dx -= cellLen * nearbyint(dx * rCellLen);
                dy -= cellLen * nearbyint(dy * rCellLen);
                dSq = dx * dx + dy * dy;
                d = sqrt(dSq);
                if (d < cellLen_2) {
                    b = floor(d / rdfDelta);
                    prdfHistAA[b] += 2;
                }
            }
        }
        //AB
        for (int i=0; i<nA; ++i) {
            xI = x[i];
            yI = y[i];
            for (int j=nA; j<n; ++j) {
                dx = xI - x[j];
                dy = yI - y[j];
                dx -= cellLen * nearbyint(dx * rCellLen);
                dy -= cellLen * nearbyint(dy * rCellLen);
                dSq = dx * dx + dy * dy;
                d = sqrt(dSq);
                if (d < cellLen_2) {
                    b = floor(d / rdfDelta);
                    prdfHistAB[b] += 2;
                }
            }
        }
        //BB
        for (int i=nA; i<n-1; ++i) {
            xI = x[i];
            yI = y[i];
            for (int j=i+1; j<n; ++j) {
                dx = xI - x[j];
                dy = yI - y[j];
                dx -= cellLen * nearbyint(dx * rCellLen);
                dy -= cellLen * nearbyint(dy * rCellLen);
                dSq = dx * dx + dy * dy;
                d = sqrt(dSq);
                if (d < cellLen_2) {
                    b = floor(d / rdfDelta);
                    prdfHistBB[b] += 2;
                }
            }
        }
    }
    if(dispersity==3) {//polydisperse calculate total rdf
        for (int i = 0; i < n - 1; ++i) {
            xI = x[i];
            yI = y[i];
            for (int j = i + 1; j < n; ++j) {
                dx = xI - x[j];
                dy = yI - y[j];
                dx -= cellLen * nearbyint(dx * rCellLen);
                dy -= cellLen * nearbyint(dy * rCellLen);
                dSq = dx * dx + dy * dy;
                d = sqrt(dSq);
                if (d < cellLen_2) {
                    b = floor(d / rdfDelta);
                    rdfHist[b] += 2;
                }
            }
        }
    }
}


void HDMC::calculateVoronoi2D(OutputFile &vor2DFile, OutputFile &vis2DFile, bool vis) {
    //Calculate Voronoi and analyse

    //Make voronoi and calculate cell sizes and neighbours
    VecF<int> cellSizeDistA,cellSizeDistB,nnCount;
    VecF<VecF<int> > cellAdjDist;
    VecF<double> cellAreaA,cellAreaB,nnSep;
    Voronoi2D vor(x, y, r, cellLen_2, nA, false, VecF<bool>(0), maxVertices);
    vor.analyse(maxVertices, cellSizeDistA, cellSizeDistB, cellAdjDist, cellAreaA, cellAreaB);
    vor.nnDistances(x,y,cellLen,rCellLen,nnSep,nnCount);

    //Add results to global results
    vor2DSizesA += cellSizeDistA;
    vor2DSizesB += cellSizeDistB;
    vor2DAreasA += cellAreaA;
    vor2DAreasB += cellAreaB;
    vor2DNNCount += nnCount;
    vor2DNNSep += nnSep;
    for (int i = 0; i < cellAdjDist.n; ++i) vor2DAdjs[i] += cellAdjDist[i];

    //Get network analysis and write for type A configuration
    for(int i=0; i<cellSizeDistA.n; ++i) if(cellSizeDistA[i]>0) cellAreaA[i]/=cellSizeDistA[i];
    VecF<double> resA = networkAnalysis(cellSizeDistA, cellAdjDist);
    vor2DFile.writeRowVector(resA);
    vor2DFile.writeRowVector(cellAreaA);

    //Get network analysis and write for type B configuration
    if(dispersity==2) {
        for(int i=0; i<cellSizeDistB.n; ++i) if(cellSizeDistB[i]>0) cellAreaB[i]/=cellSizeDistB[i];
        VecF<double> resB = networkAnalysis(cellSizeDistB, cellAdjDist);
        vor2DFile.writeRowVector(resB);
        vor2DFile.writeRowVector(cellAreaB);
    }

    //Write nearest neighbour distances and frequency
    VecF<double> nn(6);
    for(int i=0; i<3; ++i){
        if(nnCount[i]>0) nn[i]=nnSep[i]/nnCount[i];
        nn[3+i]=nnCount[i];
    }
    vor2DFile.writeRowVector(nn);

    //Write Voronoi visualisation
    if(vis) writeVor(vor,vis2DFile,1);
}


void HDMC::calculateRadical2D(OutputFile &rad2DFile, OutputFile &vis2DFile, bool vis) {
    //Calculate radical and analyse

    //Make voronoi and calculate cell sizes and neighbours
    VecF<int> cellSizeDistA,cellSizeDistB,nnCount;
    VecF<VecF<int> > cellAdjDist;
    VecF<double> cellAreaA,cellAreaB,nnSep;
    Voronoi2D rad(x, y, w, cellLen_2, nA, true, rad2DInclude, maxVertices);
    rad.analyse(maxVertices, cellSizeDistA, cellSizeDistB, cellAdjDist, cellAreaA, cellAreaB);
    rad.nnDistances(x,y,cellLen,rCellLen,nnSep,nnCount);

    //Add results to global results
    rad2DSizesA += cellSizeDistA;
    rad2DSizesB += cellSizeDistB;
    rad2DAreasA += cellAreaA;
    rad2DAreasB += cellAreaB;
    rad2DNNCount += nnCount;
    rad2DNNSep += nnSep;
    for (int i = 0; i < cellAdjDist.n; ++i) rad2DAdjs[i] += cellAdjDist[i];

    //Get network analysis and write for type A configuration
    for(int i=0; i<cellSizeDistA.n; ++i) if(cellSizeDistA[i]>0) cellAreaA[i]/=cellSizeDistA[i];
    VecF<double> resA = networkAnalysis(cellSizeDistA, cellAdjDist);
    rad2DFile.writeRowVector(resA);
    rad2DFile.writeRowVector(cellAreaA);

    //Get network analysis and write for type B configuration
    if(dispersity==2) {
        for(int i=0; i<cellSizeDistB.n; ++i) if(cellSizeDistB[i]>0) cellAreaB[i]/=cellSizeDistB[i];
        VecF<double> resB = networkAnalysis(cellSizeDistB, cellAdjDist);
        rad2DFile.writeRowVector(resB);
        rad2DFile.writeRowVector(cellAreaB);
    }

    //Write nearest neighbour distances and frequency
    VecF<double> nn(6);
    for(int i=0; i<3; ++i){
        if(nnCount[i]>0) nn[i]=nnSep[i]/nnCount[i];
        nn[3+i]=nnCount[i];
    }
    rad2DFile.writeRowVector(nn);

    //Write radical visualisation
    if(vis) writeVor(rad,vis2DFile,2,rad2DCut);
}


void HDMC::calculateVoronoi3D(OutputFile &vor3DFile, OutputFile &vis2DFile, OutputFile &vis3DFile, bool vis) {
    //Calculate Voronoi and analyse

    //Make voronoi and calculate cell sizes and neighbours
    VecF<int> cellSizeDistA,cellSizeDistB,nnCount;
    VecF<VecF<int> > cellAdjDist;
    VecF<double> cellAreaA,cellAreaB,nnSep;
    Voronoi3D vor(x, y, z, r, cellLen_2, nA, false, maxVertices);
    vor.analyse(maxVertices, cellSizeDistA, cellSizeDistB, cellAdjDist, cellAreaA, cellAreaB);
    vor.nnDistances(x,y,cellLen,rCellLen,nnSep,nnCount);

    //Add results to global results
    vor3DSizesA += cellSizeDistA;
    vor3DSizesB += cellSizeDistB;
    vor3DAreasA += cellAreaA;
    vor3DAreasB += cellAreaB;
    vor3DNNCount += nnCount;
    vor3DNNSep += nnSep;
    for (int i = 0; i < cellAdjDist.n; ++i) vor3DAdjs[i] += cellAdjDist[i];

    //Get network analysis and write for type A configuration
    for(int i=0; i<cellSizeDistA.n; ++i) if(cellSizeDistA[i]>0) cellAreaA[i]/=cellSizeDistA[i];
    VecF<double> resA = networkAnalysis(cellSizeDistA, cellAdjDist);
    vor3DFile.writeRowVector(resA);
    vor3DFile.writeRowVector(cellAreaA);

    //Get network analysis and write for type B configuration
    if(dispersity==2) {
        for(int i=0; i<cellSizeDistB.n; ++i) if(cellSizeDistB[i]>0) cellAreaB[i]/=cellSizeDistB[i];
        VecF<double> resB = networkAnalysis(cellSizeDistB, cellAdjDist);
        vor3DFile.writeRowVector(resB);
        vor3DFile.writeRowVector(cellAreaB);
    }

    //Write nearest neighbour distances and frequency
    VecF<double> nn(6);
    for(int i=0; i<3; ++i){
        if(nnCount[i]>0) nn[i]=nnSep[i]/nnCount[i];
        nn[3+i]=nnCount[i];
    }
    vor3DFile.writeRowVector(nn);

    //Write Voronoi visualisation
    if(vis) writeVor(vor,vis2DFile,vis3DFile,3);
}


void HDMC::calculateRadical3D(OutputFile &rad3DFile, OutputFile &vis2DFile, OutputFile &vis3DFile, bool vis) {
    //Calculate Voronoi and analyse

    //Make voronoi and calculate cell sizes and neighbours
    VecF<int> cellSizeDistA,cellSizeDistB,nnCount;
    VecF<VecF<int> > cellAdjDist;
    VecF<double> cellAreaA,cellAreaB,nnSep;
    Voronoi3D rad(x, y, z, r, cellLen_2, nA, true, maxVertices);
    rad.analyse(maxVertices, cellSizeDistA, cellSizeDistB, cellAdjDist, cellAreaA, cellAreaB);
    rad.nnDistances(x,y,cellLen,rCellLen,nnSep,nnCount);

    //Add results to global results
    rad3DSizesA += cellSizeDistA;
    rad3DSizesB += cellSizeDistB;
    rad3DAreasA += cellAreaA;
    rad3DAreasB += cellAreaB;
    rad3DNNCount += nnCount;
    rad3DNNSep += nnSep;
    for (int i = 0; i < cellAdjDist.n; ++i) rad3DAdjs[i] += cellAdjDist[i];

    //Get network analysis and write for type A configuration
    for(int i=0; i<cellSizeDistA.n; ++i) if(cellSizeDistA[i]>0) cellAreaA[i]/=cellSizeDistA[i];
    VecF<double> resA = networkAnalysis(cellSizeDistA, cellAdjDist);
    rad3DFile.writeRowVector(resA);
    rad3DFile.writeRowVector(cellAreaA);

    //Get network analysis and write for type B configuration
    if(dispersity==2) {
        for(int i=0; i<cellSizeDistB.n; ++i) if(cellSizeDistB[i]>0) cellAreaB[i]/=cellSizeDistB[i];
        VecF<double> resB = networkAnalysis(cellSizeDistB, cellAdjDist);
        rad3DFile.writeRowVector(resB);
        rad3DFile.writeRowVector(cellAreaB);
    }

    //Write nearest neighbour distances and frequency
    VecF<double> nn(6);
    for(int i=0; i<3; ++i){
        if(nnCount[i]>0) nn[i]=nnSep[i]/nnCount[i];
        nn[3+i]=nnCount[i];
    }
    rad3DFile.writeRowVector(nn);

    //Write Voronoi visualisation
    if(vis) writeVor(rad,vis2DFile,vis3DFile,4);
}


VecF<double> HDMC::networkAnalysis(VecF<int> &sizes, VecF< VecF<int> > &adjs) {
    //Calculate normalised size distribution and assortativity

    //Normalised size distribution and moments
    VecF<double> res(maxVertices+1);
    double normSize=vSum(sizes);
    for(int i=0; i<sizes.n; ++i) res[i]=sizes[i]/normSize;

    //Assortativity
    double normAdj;
    VecF<double> q(adjs.n);
    for(int i=0; i<adjs.n; ++i) q[i]=vSum(adjs[i]);
    normAdj=vSum(q);
    q /= normAdj;
    double sigA=0.0,sigB=0.0;
    for(int i=0; i<q.n; ++i){
        sigA += i*i*q[i];
        sigB += i*q[i];
    }
    sigB = pow(sigB,2);
    double r=0.0;
    for(int i=0; i<adjs.n; ++i){
        for(int j=0; j<adjs.n; ++j){
            r+=i*j*(adjs[i][j]/normAdj-q[i]*q[j]);
        }
    }
    r=r/(sigA-sigB);
    res[maxVertices]=r;

    return res;
}


void HDMC::writeXYZ(OutputFile &xyzFile) {
    //Write configuration to XYZ file

    xyzFile.write(n);
    xyzFile.write("");
    if(dispersity==1 or dispersity==3){
        for(int i=0; i<n; ++i) xyzFile.write("Ar"+to_string(i)+" "+to_string(x[i])+" "+to_string(y[i])+" "+to_string(z[i]));
    }
    else if(dispersity==2){
        for(int i=0; i<nA; ++i) xyzFile.write("O "+to_string(x[i])+" "+to_string(y[i])+" "+to_string(z[i]));
        for(int i=nA; i<n; ++i) xyzFile.write("S "+to_string(x[i])+" "+to_string(y[i])+" "+to_string(z[i]));
    }
    ++xyzConfigs;
}


void HDMC::writeVor(Voronoi2D &vor, OutputFile &vis2DFile, int vorCode, double param) {
    //Write voronoi visualisation to file

    //Write voronoi frame, type and additional parameter
    vis2DFile.write(xyzConfigs-1);
    vis2DFile.write(vorCode);
    vis2DFile.write(param);

    //Calculate rings
    VecF< VecR<double> > rings;
    vor.getRings(x,y,rings);

    //Write rings
    vis2DFile.write(rings.n);
    for(int i=0; i<rings.n; ++i) vis2DFile.writeRowVector(rings[i]);
}


void HDMC::writeVor(Voronoi3D &vor, OutputFile &vis2DFile, OutputFile &vis3DFile, int vorCode, double param) {
    //Write voronoi visualisation to file

    //Write voronoi frame, type and additional parameter
    vis2DFile.write(xyzConfigs-1);
    vis2DFile.write(vorCode);
    vis2DFile.write(param);

    //Calculate rings
    VecF< VecR<double> > rings;
    rings=vor.getProjectedRings();

    //Write rings
    vis2DFile.write(rings.n);
    for(int i=0; i<rings.n; ++i) vis2DFile.writeRowVector(rings[i]);

    //Write all 3D faces
    if(visVor3D){
        vis3DFile.write(xyzConfigs-1);
        vis3DFile.write(vorCode);
        VecR< VecR<double> > faces;
        VecF<double> zLimits;
        faces=vor.getFaces(zLimits);
        vis3DFile.writeRowVector(zLimits);
        vis3DFile.write(faces.n);
        for(int i=0; i<faces.n; ++i) vis3DFile.writeRowVector(faces[i]);
    }
}


void HDMC::writeAnalysis(Logfile &logfile, OutputFile &vor2DFile, OutputFile &rad2DFile, OutputFile &vor3DFile,
                         OutputFile &rad3DFile, OutputFile &diaFile) {
    //Write analysis results to files

    //RDF
    if(rdfCalc){
        OutputFile rdfFile(outputPrefix+"_rdf.dat");
        VecF<double> bins(rdfHist.n),rdf(rdfHist.n);
        for(int i=0; i<bins.n; ++i) bins[i]=rdfDelta*(i+0.5);
        if(dispersity==1){//monodisperse case total rdf only
            //Copy total rdf
            for(int i=0; i<rdf.n; ++i) rdf[i]=rdfHist[i];
            //Normalise if required
            if(rdfNorm){
                double norm=n*(n/pow(cellLen,2))*M_PI*analysisConfigs; //n*density*pi*configs
                for(int i=0; i<rdf.n; ++i){
                    rdf[i]/=norm*(pow((i+1)*rdfDelta,2)-pow(i*rdfDelta,2));
                }
            }
            //Write
            for(int i=0; i<rdf.n; ++i) rdfFile.write(bins[i],rdf[i]);
        }
        else if(dispersity==2){//bidisperse partials and total rdf
            VecF<double> prdfAA(rdf.n),prdfAB(rdf.n),prdfBB(rdf.n);
            //Copy partial rdfs
            for(int i=0; i<rdf.n; ++i){
                prdfAA[i]=prdfHistAA[i];
                prdfAB[i]=prdfHistAB[i];
                prdfBB[i]=prdfHistBB[i];
            }
            //Sum partials to obtain total rdf
            for(int i=0; i<rdf.n; ++i) rdf[i]=prdfAA[i]+prdfAB[i]+prdfBB[i];
            //Normalise if required
            if(rdfNorm){
                double norm=n*(n/pow(cellLen,2))*M_PI*analysisConfigs; //n*density*pi*configs
                double normAA=nA*(nA/pow(cellLen,2))*M_PI*analysisConfigs;
                double normAB=2*nA*(nB/pow(cellLen,2))*M_PI*analysisConfigs;
                double normBB=nB*(nB/pow(cellLen,2))*M_PI*analysisConfigs;
                for(int i=0; i<bins.n; ++i){
                    rdf[i]/=norm*(pow((i+1)*rdfDelta,2)-pow(i*rdfDelta,2));
                    prdfAA[i]/=normAA*(pow((i+1)*rdfDelta,2)-pow(i*rdfDelta,2));
                    prdfAB[i]/=normAB*(pow((i+1)*rdfDelta,2)-pow(i*rdfDelta,2));
                    prdfBB[i]/=normBB*(pow((i+1)*rdfDelta,2)-pow(i*rdfDelta,2));
                }
            }
            //Write
            VecF<double> row(5);
            for(int i=0; i<rdfHist.n; ++i){
                row[0]=bins[i];
                row[1]=rdf[i];
                row[2]=prdfAA[i];
                row[3]=prdfAB[i];
                row[4]=prdfBB[i];
                rdfFile.writeRowVector(row);
            }
        }
        else if(dispersity==3){//polydisperse case total rdf only
            //Copy total rdf
            for(int i=0; i<rdf.n; ++i) rdf[i]=rdfHist[i];
            //Normalise if required
            if(rdfNorm){
                double norm=n*(n/pow(cellLen,2))*M_PI*analysisConfigs; //n*density*pi*configs
                for(int i=0; i<rdf.n; ++i){
                    rdf[i]/=norm*(pow((i+1)*rdfDelta,2)-pow(i*rdfDelta,2));
                }
            }
            //Write
            for(int i=0; i<rdf.n; ++i) rdfFile.write(bins[i],rdf[i]);
        }
    }

    //Voronoi
    if(vorCalc2D){
        for(int i=0; i<vor2DSizesA.n; ++i) if(vor2DSizesA[i]>0) vor2DAreasA[i]/=vor2DSizesA[i];
        vor2DAreasA[vor2DAreasA.n-1] /= analysisConfigs;
        VecF<double> resA=networkAnalysis(vor2DSizesA,vor2DAdjs);
        vor2DFile.writeRowVector(resA);
        vor2DFile.writeRowVector(vor2DAreasA);
        if(dispersity==2){
            for(int i=0; i<vor2DSizesB.n; ++i) if(vor2DSizesB[i]>0) vor2DAreasB[i]/=vor2DSizesB[i];
            vor2DAreasB[vor2DAreasB.n-1] /= analysisConfigs;
            VecF<double> resB=networkAnalysis(vor2DSizesB,vor2DAdjs);
            vor2DFile.writeRowVector(resB);
            vor2DFile.writeRowVector(vor2DAreasB);
        }
        VecF<double> nn(6);
        for(int i=0; i<3; ++i){
            if(vor2DNNCount[i]>0) nn[i]=vor2DNNSep[i]/vor2DNNCount[i];
            nn[3+i]=double(vor2DNNCount[i])/analysisConfigs;
        }
        vor2DFile.writeRowVector(nn);
    }

    //Radical
    if(radCalc2D){
        for(int i=0; i<rad2DSizesA.n; ++i) if(rad2DSizesA[i]>0) rad2DAreasA[i]/=rad2DSizesA[i];
        rad2DAreasA[rad2DAreasA.n-1] /= analysisConfigs;
        VecF<double> resA=networkAnalysis(rad2DSizesA,rad2DAdjs);
        rad2DFile.writeRowVector(resA);
        rad2DFile.writeRowVector(rad2DAreasA);
        if(dispersity==2){
            for(int i=0; i<rad2DSizesB.n; ++i) if(rad2DSizesB[i]>0) rad2DAreasB[i]/=rad2DSizesB[i];
            rad2DAreasB[rad2DAreasB.n-1] /= analysisConfigs;
            VecF<double> resB=networkAnalysis(rad2DSizesB,rad2DAdjs);
            rad2DFile.writeRowVector(resB);
            rad2DFile.writeRowVector(rad2DAreasB);
        }
        VecF<double> nn(6);
        for(int i=0; i<3; ++i){
            if(rad2DNNCount[i]>0) nn[i]=rad2DNNSep[i]/rad2DNNCount[i];
            nn[3+i]=double(rad2DNNCount[i])/analysisConfigs;
        }
        rad2DFile.writeRowVector(nn);
    }

    //Voronoi
    if(vorCalc3D){
        for(int i=0; i<vor3DSizesA.n; ++i) if(vor3DSizesA[i]>0) vor3DAreasA[i]/=vor3DSizesA[i];
        vor3DAreasA[vor3DAreasA.n-1] /= analysisConfigs;
        VecF<double> resA=networkAnalysis(vor3DSizesA,vor3DAdjs);
        vor3DFile.writeRowVector(resA);
        vor3DFile.writeRowVector(vor3DAreasA);
        if(dispersity==2){
            for(int i=0; i<vor3DSizesB.n; ++i) if(vor3DSizesB[i]>0) vor3DAreasB[i]/=vor3DSizesB[i];
            vor3DAreasB[vor3DAreasB.n-1] /= analysisConfigs;
            VecF<double> resB=networkAnalysis(vor3DSizesB,vor3DAdjs);
            vor3DFile.writeRowVector(resB);
            vor3DFile.writeRowVector(vor3DAreasB);
        }
        VecF<double> nn(6);
        for(int i=0; i<3; ++i){
            if(vor3DNNCount[i]>0) nn[i]=vor3DNNSep[i]/vor3DNNCount[i];
            nn[3+i]=double(vor3DNNCount[i])/analysisConfigs;
        }
        vor3DFile.writeRowVector(nn);
    }

    //Radical
    if(radCalc3D){
        for(int i=0; i<rad3DSizesA.n; ++i) if(rad3DSizesA[i]>0) rad3DAreasA[i]/=rad3DSizesA[i];
        rad3DAreasA[rad3DAreasA.n-1] /= analysisConfigs;
        VecF<double> resA=networkAnalysis(rad3DSizesA,rad3DAdjs);
        rad3DFile.writeRowVector(resA);
        rad3DFile.writeRowVector(rad3DAreasA);
        if(dispersity==2){
            for(int i=0; i<rad3DSizesB.n; ++i) if(rad3DSizesB[i]>0) rad3DAreasB[i]/=rad3DSizesB[i];
            rad3DAreasB[rad3DAreasB.n-1] /= analysisConfigs;
            VecF<double> resB=networkAnalysis(rad3DSizesB,rad3DAdjs);
            rad3DFile.writeRowVector(resB);
            rad3DFile.writeRowVector(rad3DAreasB);
        }
        VecF<double> nn(6);
        for(int i=0; i<3; ++i){
            if(rad3DNNCount[i]>0) nn[i]=rad3DNNSep[i]/rad3DNNCount[i];
            nn[3+i]=double(rad3DNNCount[i])/analysisConfigs;
        }
        rad3DFile.writeRowVector(nn);
    }

    //Diameters and weights
    for(int i=0; i<n; ++i) diaFile.write(2.0*r[i]);
    for(int i=0; i<n; ++i) diaFile.write(w[i]);

}

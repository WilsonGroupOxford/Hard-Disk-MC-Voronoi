#include "hdmc.h"

//-------- HARD DISK MONTE CARLO CLASS --------


//-------- CONSTRUCTORS, SETTERS --------


HDMC::HDMC() {
    //Default constructor

    n=0;
    phi=0;
}


int HDMC::setParticles(int num, double packFrac, int disp, VecF<double> dispParams, int interact) {
    //Set particle parameters

    //Set parameters
    n=num;
    interaction=interact;
    phi=packFrac;

    //Allocate vectors
    x=VecF<double>(n);
    y=VecF<double>(n);
    r=VecF<double>(n);

    //Initialise particle lattice
    int lattice;
    dispersity=disp;
    if(disp==1) lattice=initMono(dispParams);
    if(lattice==0){//apply periodic boundary condition
        for(int i=0; i<n; ++i){
            x[i]-=cellLen*nearbyint(x[i]*rCellLen);
            y[i]-=cellLen*nearbyint(y[i]*rCellLen);
        }
    }

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
        for(int j=0;j<rows; ++j){
            x[i*rows+j]=j*spacing;
            y[i*rows+j]=i*spacing;
        }
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


int HDMC::setAnalysis(string path, int xyzFreq, int anFreq, int rdf, double rdfDel, int vor) {
    //Set analysis parameters

    outputPrefix=path;
    analysisFreq=anFreq;

    //Set xyz output frequency with 0 preventing write
    if(xyzFreq==0){
        xyzWrite=false;
        xyzWriteFreq=1;
    }
    else{
        xyzWrite=true;
        xyzWriteFreq=xyzFreq;
    }

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

    //Set voronoi type
    vorCalc=false;
    radCalc=false;
    if(vor==1) vorCalc=true;
    else if(vor==2) radCalc=true;
    else if(vor==3){
        vorCalc=true;
        radCalc=true;
    }
    maxVertices=21;

    //Initialise analysis tools
    initAnalysis();

    return 0;
}


int HDMC::initAnalysis() {
    //Initialise analysis tools

    analysisConfigs=0;

    //RDF histogram
    if(rdfCalc) rdfHist=VecF<int>(floor(cellLen_2/rdfDelta)+1); //max distance is half cell size

    //Voronoi distributions
    if(vorCalc){
        vorSizes=VecF<int>(maxVertices);
        vorAdjs=VecF< VecF<int> >(maxVertices);
        for(int i=0; i<maxVertices; ++i) vorAdjs[i]=VecF<int>(maxVertices);
    }
    if(radCalc){
        radSizes=VecF<int>(maxVertices);
        radAdjs=VecF< VecF<int> >(maxVertices);
        for(int i=0; i<maxVertices; ++i) radAdjs[i]=VecF<int>(maxVertices);
    }

    return 0;
}


//-------- MONTE CARLO MOVES --------


inline int HDMC::mcCycle() {
    //Cycle of n-particle Monte Carlo moves

    int accCount=0;
    if(interaction==0){
        for(int i=0; i<n; ++i) mcAdditiveMove(accCount);
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
        double rJ=rI;
        xI=x[pJ];
        yI=y[pJ];
        rI=r[pJ];

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
            r[pI]=rI;
            x[pJ]=xJ;
            y[pJ]=yJ;
            r[pJ]=rJ;
            ++counter;
        }
    }
}


//-------- MONTE CARLO SIMULATION --------


void HDMC::equilibration(Logfile &logfile, OutputFile &xyzFile) {
    //Equilibration Monte Carlo

    //Header
    logfile.write("Equilibration Monte Carlo");
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
        if(i%logMoves==0) logfile.write("Moves and acceptance:",i,double(accCount)/(i*n));
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


void HDMC::production(Logfile &logfile, OutputFile &xyzFile, OutputFile &vorFile, OutputFile &radFile) {
    //Production Monte Carlo

    //Production cycles
    logfile.write("Production Monte Carlo");
    ++logfile.currIndent;
    int logMoves=prodCycles/100;
    int accCount=0;
    for (int i = 1; i<=prodCycles; ++i) {
        accCount+=mcCycle();
        if(i%logMoves==0) logfile.write("Moves and acceptance:",i,double(accCount)/(i*n));
        if(xyzWrite && i%xyzWriteFreq==0) writeXYZ(xyzFile);
        if(i%analysisFreq==0) analyseConfiguration(vorFile,radFile);
    }
    logfile.currIndent-=2;
    logfile.separator();
}


//--------- ANALYSIS ----------


void HDMC::analyseConfiguration(OutputFile &vorFile, OutputFile &radFile) {
    //Control analysis of current configuration

    if(rdfCalc) calculateRDF();
    if(vorCalc) calculateVoronoi(vorFile);

    ++analysisConfigs;
}


void HDMC::calculateRDF() {
    //Calculate RDF for current configuration

    double xI,yI,b;
    double dx,dy,dSq,d;
    for(int i=0; i<n-1; ++i){
        xI=x[i];
        yI=y[i];
        for(int j=i+1; j<n; ++j){
            dx=xI-x[j];
            dy=yI-y[j];
            dx-=cellLen*nearbyint(dx*rCellLen);
            dy-=cellLen*nearbyint(dy*rCellLen);
            dSq=dx*dx+dy*dy;
            d=sqrt(dSq);
            if(d<cellLen_2){
                b=floor(d/rdfDelta);
                rdfHist[b]+=2;
            }
        }
    }
}


void HDMC::calculateVoronoi(OutputFile &vorFile) {
    //Calculate Voronoi and analyse

    //Make voronoi and calculate cell sizes and neighbours
    VecF<int> cellSizeDist;
    VecF<VecF<int> > cellAdjDist;
    Voronoi vor(x, y, r, cellLen_2, false);
    vor.analyse(maxVertices, cellSizeDist, cellAdjDist);

    //Add results to global results
    vorSizes += cellSizeDist;
    for (int i = 0; i < cellAdjDist.n; ++i) vorAdjs[i] += cellAdjDist[i];

    //Get network analysis for configuration
    VecF<double> res = networkAnalysis(cellSizeDist, cellAdjDist);

    vorFile.writeRowVector(res);
}


VecF<double> HDMC::networkAnalysis(VecF<int> &sizes, VecF< VecF<int> > &adjs) {
    //Calculate normalised size distribution and assortativity

    //Normalised size distribution and moments
    VecF<double> res(maxVertices+1);
    double normSize=vSum(sizes);
    for(int i=0; i<sizes.n; ++i) res[i]=sizes[i]/normSize;
    double k1=0.0,k2=0.0,k3=0.0;
    for(int i=0; i<sizes.n; ++i){
        k1+=i*res[i];
        k2+=i*i*res[i];
        k3+=i*i*i*res[i];
    }

    //Assortativity
    double normAdj=0.0;
    for(int i=0; i<adjs.n; ++i) normAdj+=vSum(adjs[i]);
    double r=0.0;
    for(int i=0; i<adjs.n; ++i){
        for(int j=0; j<adjs.n; ++j){
            r+=i*j*adjs[i][j];
        }
    }
    r=(k1*k1*r/normAdj-k2*k2)/(k1*k3-k2*k2);
    res[maxVertices]=r;

    return res;
}


void HDMC::writeXYZ(OutputFile &xyzFile) {
    //Write configuration to XYZ file

    xyzFile.write(n);
    xyzFile.write("");
    for(int i=0; i<n; ++i){
        xyzFile.write("Ar "+to_string(x[i])+" "+to_string(y[i])+" 0.0");
    }
}


void HDMC::writeAnalysis(Logfile &logfile, OutputFile &vorFile, OutputFile &radFile) {
    //Write analysis results to files

    //RDF
    if(rdfCalc){
        OutputFile rdfFile(outputPrefix+"_rdf.dat");
        VecF<double> rdfVals(rdfHist.n),rdfBins(rdfHist.n);
        for(int i=0; i<rdfHist.n; ++i){
            rdfBins[i]=rdfDelta*(i+0.5);
            rdfVals[i]=rdfHist[i];
        }
        if(rdfNorm){
            double norm=n*(n/pow(cellLen,2))*M_PI*analysisConfigs; //n*density*pi*configs
            for(int i=0; i<rdfVals.n; ++i){
                rdfVals[i]/=norm*(pow((i+1)*rdfDelta,2)-pow(i*rdfDelta,2));
            }
        }
        for(int i=0; i<rdfBins.n; ++i) rdfFile.write(rdfBins[i],rdfVals[i]);
    }

    //Voronoi
    if(vorCalc){
        VecF<double> res=networkAnalysis(vorSizes,vorAdjs);
        vorFile.writeRowVector(res);
    }
}

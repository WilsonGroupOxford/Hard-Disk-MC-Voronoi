#include "configuration.h"

Configuration::Configuration() {}

Configuration::Configuration(int numA, int numB, double radA, double radB, double cellLength) {
    //Construct with configuration information

    //Assign
    nA = numA;
    nB = numB;
    rA = radA;
    rB = radB;
    cellLen = cellLength;

    //Additional
    nCrdSets = 0;
    nC = nA+nB;
    rCellLen = 1.0/cellLen;
    cellLen_2 = cellLen/2;
    xA = VecF<double>(nA);
    yA = VecF<double>(nA);
    xB = VecF<double>(nB);
    yB = VecF<double>(nB);
}

void Configuration::setCoordinates(ifstream& xyzFile, Logfile& logfile) {
    //Read next configuration from XYZ file

    //Load configuration
    logfile.write("Configuration loaded: ",nCrdSets);
    string skip,line;
    istringstream ss("");
    getline(xyzFile,skip);
    getline(xyzFile,skip);

    //Read type a
    for(int i=0; i<nA; ++i){
        getline(xyzFile,line);
        ss.str(line);
        ss>>skip;
        ss>>xA[i];
        ss>>yA[i];
    }

    //Read type b
    for(int i=0; i<nB; ++i){
        getline(xyzFile,line);
        ss.str(line);
        ss>>skip;
        ss>>xB[i];
        ss>>yB[i];
    }

    //Increment total coordinate sets added
    ++nCrdSets;
}

void Configuration::setRdf(double delta, double extent) {
    //Initialise partial and total rdf histograms

    //Histogram Bin properties
    rdfDelta = delta;
    double rdfExtent = extent*rB;
    if(rdfExtent>cellLen_2) rdfExtent = cellLen_2;
    rdfMaxSq = rdfExtent*rdfExtent;

    //Make histograms
    int n = int(floor(rdfExtent/rdfDelta));
    rdfR = VecF<double>(n); //bin central values
    rdfR = rdfDelta/2;
    for(int i=0; i<n; ++i) rdfR[i] += i*rdfDelta;
    rdfAA = VecF<double>(n);
    rdfBB = VecF<double>(n);
    rdfAB = VecF<double>(n);
    rdfC = VecF<double>(n);
    rdfAA = 0;
    rdfBB = 0;
    rdfAB = 0;
    rdfC = 0;
}

void Configuration::rdf(Logfile &logfile) {
    //Add to unnormalised partial rdfs

    //AA
    int r;
    double x0,y0,x1,y1,dx,dy,rSq;
    for(int i=0; i<nA-1; ++i){
        x0 = xA[i];
        y0 = yA[i];
        for(int j=i+1; j<nA; ++j){
            x1 = xA[j];
            y1 = yA[j];
            dx = x0-x1;
            dy = y0-y1;
            dx -= cellLen*nearbyint(dx*rCellLen);
            dy -= cellLen*nearbyint(dy*rCellLen);
            rSq = dx*dx+dy*dy;
            if(rSq<rdfMaxSq){
                r=int(floor(sqrt(rSq)/rdfDelta));
                rdfAA[r] += 2;
            }
        }
    }

    //BB
    for(int i=0; i<nB-1; ++i){
        x0 = xB[i];
        y0 = yB[i];
        for(int j=i+1; j<nB; ++j){
            x1 = xB[j];
            y1 = yB[j];
            dx = x0-x1;
            dy = y0-y1;
            dx -= cellLen*nearbyint(dx*rCellLen);
            dy -= cellLen*nearbyint(dy*rCellLen);
            rSq = dx*dx+dy*dy;
            if(rSq<rdfMaxSq){
                r=int(floor(sqrt(rSq)/rdfDelta));
                rdfBB[r] += 2;
            }
        }
    }

    //AB
    for(int i=0; i<nA; ++i){
        x0 = xA[i];
        y0 = yA[i];
        for(int j=0; j<nB; ++j){
            x1 = xB[j];
            y1 = yB[j];
            dx = x0-x1;
            dy = y0-y1;
            dx -= cellLen*nearbyint(dx*rCellLen);
            dy -= cellLen*nearbyint(dy*rCellLen);
            rSq = dx*dx+dy*dy;
            if(rSq<rdfMaxSq){
                r=int(floor(sqrt(rSq)/rdfDelta));
                rdfAB[r] += 2;
            }
        }
    }
}

void Configuration::rdfFinalise(string prefix, Logfile &logfile) {
    //Combine partial rdfs, normalise and write to file

    //Combine partial RDFs
    logfile.write("Finalising RDFs");
    ++logfile.currIndent;
    rdfC = rdfAA+rdfBB+rdfAB;
    logfile.write("Partial RDFs combined");

    //Normalise
    double ndensityA,ndensityB,ndensityC;
    ndensityA = nA/(cellLen*cellLen);
    ndensityB = nB/(cellLen*cellLen);
    ndensityC = nC/(cellLen*cellLen);
    VecF<double> norm,normA,normB,normAB,normC;
    norm = ((rdfR+rdfDelta/2)*(rdfR+rdfDelta/2)-(rdfR-rdfDelta/2)*(rdfR-rdfDelta/2))*M_PI*nCrdSets;
    rdfAA /= norm*ndensityA*nA;
    rdfBB /= norm*ndensityB*nB;
    rdfAB /= norm*(ndensityA*nB+ndensityB*nA);
    rdfC /= norm*ndensityC*nC;
    logfile.write("RDFs normalised");

    //Write to file
    ofstream rdfFile(prefix+"_rdf.dat", ios::in | ios::trunc);
    rdfFile << fixed << showpoint << setprecision(8);
    for(int i=0; i<rdfR.n; ++i){
        rdfFile << setw(20) << left << rdfR[i];
        rdfFile << setw(20) << left << rdfAA[i];
        rdfFile << setw(20) << left << rdfBB[i];
        rdfFile << setw(20) << left << rdfAB[i];
        rdfFile << setw(20) << left << rdfC[i]<<endl;
    }
    logfile.write("RDFs written (AA,BB,AB,Total");
    --logfile.currIndent;
}

void Configuration::setVoronoi(Logfile &logfile) {
    //Initialise ring distributions for type a, b and total

    //Maximum size arbitrary - but random Voronoi unlikely to be >16
    int maxK = 20;
    vorPKA = VecF<double>(maxK+1);
    vorPKB = VecF<double>(maxK+1);
    vorPKC = VecF<double>(maxK+1);
    vorARA = VecF<double>(maxK+1);
    vorARB = VecF<double>(maxK+1);
    vorE = VecF< VecF<double> >(maxK+1);
    vorPKA = 0;
    vorPKB = 0;
    vorPKC = 0;
    vorARA = 0;
    vorARB = 0;
    for(int i=0; i<=maxK; ++i){
        vorE[i] = VecF<double>(maxK+1);
        vorE[i] = 0;
    }
}

void Configuration::voronoi(ofstream &vorFilePKA, ofstream &vorFilePKB, ofstream &vorFilePKC, ofstream &vorFileARA,
                            ofstream &vorFileARB, ofstream &vorFileNet, Logfile &logfile) {
    //Voronoi analysis

    //Initialise Voronoi and add configuration coordinates
    voro::container vor3D(-cellLen_2,cellLen_2,-cellLen_2,cellLen_2,-0.5,0.5,3,3,3,true,true,false,nC);
    for(int i=0; i<nA; ++i) vor3D.put(i,xA[i],yA[i],0.0);
    for(int i=0, j=nA; i<nB; ++i,++j) vor3D.put(j,xB[i],yB[i],0.0);

    //Calculate Voronoi and generate temporary files
    vor3D.print_custom("%i","./vorI.tmp");
    vor3D.print_custom("%a","./vorK.tmp");
    vor3D.print_custom("%f","./vorA.tmp");
    vor3D.print_custom("%n","./vorN.tmp");
    vor3D.print_custom("%l","./vorV.tmp");

    //Read into 2D Voronoi and extract unnormalised distributions
    int maxK = 20;
    VecF<double> pA,pB,pC,aA,aB;
    VecF< VecF<double> > e;
    VoronoiBinary2D vor2D(nA,nB,maxK,logfile);
    vor2D.getDistributions(pA,pB,pC,aA,aB,e);

    //Add to running total distributions
    vorPKA += pA;
    vorPKB += pB;
    vorPKC += pC;
    vorARA += aA;
    vorARB += aB;
    for(int i=0; i<=maxK; ++i) vorE[i] += e[i];

    //Normalise distributions
    VecF<double> k(maxK+1);
    for(int i=0; i<=maxK; ++i) k[i]=i;
    for(int i=0; i<=maxK; ++i){
        if(pA[i]>0) aA[i] /= pA[i];
        if(pB[i]>0) aB[i] /= pB[i];
    }
    if(nA>0) pA /= nA;
    if(nB>0) pB /= nB;
    pC /= nC;
    VecF<double> q(maxK+1);
    for(int i=0; i<=maxK; ++i) q[i] = vSum(e[i]);
    double normE = vSum(q);
    for(int i=0; i<=maxK; ++i) e[i] /= normE;
    q /= normE;

    //Calculate distirbution metrics
    double meanA,meanB,meanC,varA,varB,varC;
    meanA = vSum(k*pA);
    meanB = vSum(k*pB);
    meanC = vSum(k*pC);
    varA = vSum(k*k*pA)-meanA*meanA;
    varB = vSum(k*k*pB)-meanB*meanB;
    varC = vSum(k*k*pC)-meanC*meanC;
    double assortativity=0.0;
    for(int i=0; i<=maxK; ++i){
        for(int j=0; j<=maxK; ++j){
            assortativity += i*j*(e[i][j]-q[i]*q[j]);
        }
    }
    assortativity /= vSum(k*k*q)-pow(vSum(k*q),2);
    VecR<double> awX(0,maxK+1),awY(0,maxK+1);
    for(int i=0; i<=maxK; ++i){
        if(pC[i]>0){
            awX.addValue(6.0*(i-6.0));
            awY.addValue(i*vSum(k*e[i])/q[i]);
        }
    }
    VecR<double> aw=vLinearRegression(awX,awY);
    double awA = 1.0-aw[0];
    double awVar = aw[1]-36.0;
    double awRSq = aw[2];

    //Write to files
    for(int i=0; i<=maxK; ++i){
        vorFilePKA<<setw(20)<<left<<pA[i];
        vorFilePKB<<setw(20)<<left<<pB[i];
        vorFilePKC<<setw(20)<<left<<pC[i];
        vorFileARA<<setw(20)<<left<<aA[i];
        vorFileARB<<setw(20)<<left<<aB[i];
    }
    vorFilePKA<<endl;
    vorFilePKB<<endl;
    vorFilePKC<<endl;
    vorFileARA<<endl;
    vorFileARB<<endl;
    vorFileNet<<setw(20)<<left<<meanA;
    vorFileNet<<setw(20)<<left<<meanB;
    vorFileNet<<setw(20)<<left<<meanC;
    vorFileNet<<setw(20)<<left<<varA;
    vorFileNet<<setw(20)<<left<<varB;
    vorFileNet<<setw(20)<<left<<varC;
    vorFileNet<<setw(20)<<left<<assortativity;
    vorFileNet<<setw(20)<<left<<awA;
    vorFileNet<<setw(20)<<left<<awVar;
    vorFileNet<<setw(20)<<left<<awRSq;
    vorFileNet<<endl;
}

void Configuration::voronoiFinalise(ofstream &vorFilePKA, ofstream &vorFilePKB, ofstream &vorFilePKC, ofstream &vorFileARA,
                            ofstream &vorFileARB, ofstream &vorFileNet, Logfile &logfile) {
    //Analyse distributions from all frames

    //Normalise distributions
    int maxK = 20;
    VecF<double> k(maxK+1);
    for(int i=0; i<=maxK; ++i) k[i]=i;
    for(int i=0; i<=maxK; ++i){
        if(vorPKA[i]>0) vorARA[i] /= vorPKA[i];
        if(vorPKB[i]>0) vorARB[i] /= vorPKB[i];
    }
    if(nA>0) vorPKA /= nA*nCrdSets;
    if(nB>0) vorPKB /= nB*nCrdSets;
    vorPKC /= nC*nCrdSets;
    VecF<double> q(maxK+1);
    for(int i=0; i<=maxK; ++i) q[i] = vSum(vorE[i]);
    double normE = vSum(q);
    for(int i=0; i<=maxK; ++i) vorE[i] /= normE;
    q /= normE;

    //Calculate distirbution metrics
    double meanA,meanB,meanC,varA,varB,varC;
    meanA = vSum(k*vorPKA);
    meanB = vSum(k*vorPKB);
    meanC = vSum(k*vorPKC);
    varA = vSum(k*k*vorPKA)-meanA*meanA;
    varB = vSum(k*k*vorPKB)-meanB*meanB;
    varC = vSum(k*k*vorPKC)-meanC*meanC;
    double assortativity=0.0;
    for(int i=0; i<=maxK; ++i){
        for(int j=0; j<=maxK; ++j){
            assortativity += i*j*(vorE[i][j]-q[i]*q[j]);
        }
    }
    assortativity /= vSum(k*k*q)-pow(vSum(k*q),2);
    VecR<double> awX(0,maxK+1),awY(0,maxK+1);
    for(int i=0; i<=maxK; ++i){
        if(vorPKC[i]>0){
            awX.addValue(6.0*(i-6.0));
            awY.addValue(i*vSum(k*vorE[i])/q[i]);
        }
    }
    VecR<double> aw=vLinearRegression(awX,awY);
    double awA = 1.0-aw[0];
    double awVar = aw[1]-36.0;
    double awRSq = aw[2];

    //Write to files
    for(int i=0; i<=maxK; ++i){
        vorFilePKA<<setw(20)<<left<<vorPKA[i];
        vorFilePKB<<setw(20)<<left<<vorPKB[i];
        vorFilePKC<<setw(20)<<left<<vorPKC[i];
        vorFileARA<<setw(20)<<left<<vorARA[i];
        vorFileARB<<setw(20)<<left<<vorARB[i];
    }
    vorFilePKA<<endl;
    vorFilePKB<<endl;
    vorFilePKC<<endl;
    vorFileARA<<endl;
    vorFileARB<<endl;
    vorFileNet<<setw(20)<<left<<meanA;
    vorFileNet<<setw(20)<<left<<meanB;
    vorFileNet<<setw(20)<<left<<meanC;
    vorFileNet<<setw(20)<<left<<varA;
    vorFileNet<<setw(20)<<left<<varB;
    vorFileNet<<setw(20)<<left<<varC;
    vorFileNet<<setw(20)<<left<<assortativity;
    vorFileNet<<setw(20)<<left<<awA;
    vorFileNet<<setw(20)<<left<<awVar;
    vorFileNet<<setw(20)<<left<<awRSq;
    vorFileNet<<endl;

    //Close files
    vorFilePKA.close();
    vorFilePKB.close();
    vorFilePKC.close();
    vorFileARA.close();
    vorFileARB.close();
    vorFileNet.close();

    //Remove any .tmp files
    remove("*.tmp");
}
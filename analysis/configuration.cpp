#include "configuration.h"

Configuration::Configuration() {}

Configuration::Configuration(int numA, int numB, double radiusA, double radiusB, double cellLength, int nonAdditive) {
    //Construct with configuration information

    //Assign
    nA = numA;
    nB = numB;
    rA = radiusA;
    rB = radiusB;
    cellLen = cellLength;

    //Additional
    nCrdSets = 0;
    nC = nA+nB;
    rCellLen = 1.0/cellLen;
    cellLen_2 = cellLen/2;
    if(nonAdditive==1) {
        wA = 2 * rA * sqrt(rA * rB) / (rA + rB);
        wB = 2 * sqrt(rA * rB) * (1 - rA / (rA + rB));
    }
    else{
        wA = rA;
        wB = rB;
    }
    
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
    int n = int(floor(rdfExtent/rdfDelta))+1;
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
    vorELA = VecF<double>(maxK+1);
    vorELB = VecF<double>(maxK+1);
    vorE = VecF< VecF<double> >(maxK+1);
    vorLEA = VecF<double>(3);
    vorLEB = VecF<double>(3);
    vorANA = VecF<double>(3);
    vorANB = VecF<double>(3);
    vorPKA = 0;
    vorPKB = 0;
    vorPKC = 0;
    vorARA = 0;
    vorARB = 0;
    vorELA = 0;
    vorELB = 0;
    vorLEA = 0;
    vorLEB = 0;
    vorANA = 0;
    vorANB = 0;
    for(int i=0; i<=maxK; ++i){
        vorE[i] = VecF<double>(maxK+1);
        vorE[i] = 0;
    }
}

void Configuration::voronoi(ofstream &vorFilePKA, ofstream &vorFilePKB, ofstream &vorFilePKC, ofstream &vorFileEJK,
                            ofstream &vorFileARA, ofstream &vorFileARB, ofstream &vorFileELA, ofstream &vorFileELB,
                            ofstream &vorFileNet, Logfile &logfile) {
    //Voronoi analysis

    //Initialise Voronoi and add configuration coordinates
    voro::container vor3D(-cellLen_2,cellLen_2,-cellLen_2,cellLen_2,-0.5,0.5,10,10,1,true,true,false,nC);
    for(int i=0; i<nA; ++i) vor3D.put(i,xA[i],yA[i],0.0);
    for(int i=0, j=nA; i<nB; ++i,++j) vor3D.put(j,xB[i],yB[i],0.0);

    //Calculate Voronoi and generate temporary files
    vor3D.print_custom("%i","./vorI.tmp");
    vor3D.print_custom("%a","./vorK.tmp");
    vor3D.print_custom("%f","./vorA.tmp");
    vor3D.print_custom("%e","./vorE.tmp");
    vor3D.print_custom("%n","./vorN.tmp");
    vor3D.print_custom("%l","./vorV.tmp");
    vor3D.print_custom("%t","./vorT.tmp");
    vor3D.print_custom("%p","./vorP.tmp");

    //Read into 2D Voronoi and extract unnormalised distributions
    int maxK = 20;
    VecF<double> pA,pB,pC,aA,aB,eA,eB,lA,lB,anA,anB;
    VecF< VecF<double> > e;
    VoronoiBinary2D vor2D(nA,nB,maxK,false,logfile);
    vor2D.getDistributions(pA,pB,pC,aA,aB,eA,eB,e,lA,lB,anA,anB);

    //Add to running total distributions
    vorPKA += pA;
    vorPKB += pB;
    vorPKC += pC;
    vorARA += aA;
    vorARB += aB;
    vorELA += eA;
    vorELB += eB;
    vorLEA += lA;
    vorLEB += lB;
    vorANA += anA;
    vorANB += anB;
    for(int i=0; i<=maxK; ++i) vorE[i] += e[i];

    //Normalise distributions
    VecF<double> k(maxK+1);
    for(int i=0; i<=maxK; ++i) k[i]=i;
    for(int i=0; i<=maxK; ++i){
        if(pA[i]>0){
            aA[i] /= pA[i];
            eA[i] /= pA[i];
        }
        if(pB[i]>0){
            aB[i] /= pB[i];
            eB[i] /= pB[i];
        }
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
    VecF<double> p6I(maxK+1),meanI(maxK+1),varI(maxK+1);
    for(int i=0; i<=maxK; ++i){
        if(pC[i]>0){
            p6I[i]=e[i][6]/q[i];
            meanI[i]=vSum(k*e[i]/q[i]);
            varI[i]=vSum(k*k*e[i]/q[i])-meanI[i]*meanI[i];
        }
    }

    //Write to files
    for(int i=0; i<=maxK; ++i){
        vorFilePKA<<setw(20)<<left<<pA[i];
        vorFilePKB<<setw(20)<<left<<pB[i];
        vorFilePKC<<setw(20)<<left<<pC[i];
        vorFileARA<<setw(20)<<left<<aA[i];
        vorFileARB<<setw(20)<<left<<aB[i];
        vorFileELA<<setw(20)<<left<<eA[i];
        vorFileELB<<setw(20)<<left<<eB[i];
//        vorFileEJK<<setw(20)<<left<<p6I[i];
//        vorFileEJK<<setw(20)<<left<<meanI[i];
//        vorFileEJK<<setw(20)<<left<<varI[i];
//        vorFileEJK<<endl;
    }
    vorFilePKA<<endl;
    vorFilePKB<<endl;
    vorFilePKC<<endl;
    vorFileARA<<endl;
    vorFileARB<<endl;
    vorFileELA<<endl;
    vorFileELB<<endl;
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

void Configuration::voronoiFinalise(ofstream &vorFilePKA, ofstream &vorFilePKB, ofstream &vorFilePKC,
                                    ofstream &vorFileEJK, ofstream &vorFileARA, ofstream &vorFileARB,
                                    ofstream &vorFileELA, ofstream &vorFileELB, ofstream &vorFileLAN,
                                    ofstream &vorFileNet, Logfile &logfile) {
    //Analyse distributions from all frames

    //Normalise distributions
    int maxK = 20;
    VecF<double> k(maxK+1);
    for(int i=0; i<=maxK; ++i) k[i]=i;
    for(int i=0; i<=maxK; ++i){
        if(vorPKA[i]>0){
            vorARA[i] /= vorPKA[i];
            vorELA[i] /= vorPKA[i];
        }
        if(vorPKB[i]>0){
            vorARB[i] /= vorPKB[i];
            vorELB[i] /= vorPKB[i];
        }
    }
    if(nA>0){
        vorPKA /= nA*nCrdSets;
        vorLEA[0] /= vorLEA[2];
        vorLEA[1] /= vorLEA[2];
        vorLEA[1] -= vorLEA[0]*vorLEA[0];
        vorANA[0] /= vorANA[2];
        vorANA[1] /= vorANA[2];
        vorANA[1] -= vorANA[0]*vorANA[0];
    }
    if(nB>0){
        vorPKB /= nB*nCrdSets;
        vorLEB[0] /= vorLEB[2];
        vorLEB[1] /= vorLEB[2];
        vorLEB[1] -= vorLEB[0]*vorLEB[0];
        vorANB[0] /= vorANB[2];
        vorANB[1] /= vorANB[2];
        vorANB[1] -= vorANB[0]*vorANB[0];
    }
    vorPKC /= nC*nCrdSets;
    VecF<double> q(maxK+1);
    for(int i=0; i<=maxK; ++i) q[i] = vSum(vorE[i]);
    double normE = vSum(q);
    for(int i=0; i<=maxK; ++i) vorE[i] /= normE;
    q /= normE;

    //Calculate distribution metrics
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
    VecF<double> p6I(maxK+1),meanI(maxK+1),varI(maxK+1);
    for(int i=0; i<=maxK; ++i){
        if(vorPKC[i]>0){
            p6I[i]=vorE[i][6]/q[i];
            meanI[i]=vSum(k*vorE[i]/q[i]);
            varI[i]=vSum(k*k*vorE[i]/q[i])-meanI[i]*meanI[i];
        }
    }

    //Write to files
    for(int i=0; i<=maxK; ++i){
        vorFilePKA<<setw(20)<<left<<vorPKA[i];
        vorFilePKB<<setw(20)<<left<<vorPKB[i];
        vorFilePKC<<setw(20)<<left<<vorPKC[i];
        vorFileARA<<setw(20)<<left<<vorARA[i];
        vorFileARB<<setw(20)<<left<<vorARB[i];
        vorFileELA<<setw(20)<<left<<vorELA[i];
        vorFileELB<<setw(20)<<left<<vorELB[i];
//        vorFileEJK<<setw(20)<<left<<p6I[i];
//        vorFileEJK<<setw(20)<<left<<meanI[i];
//        vorFileEJK<<setw(20)<<left<<varI[i];
//        vorFileEJK<<endl;
    }
    for(int i=0; i<=maxK; ++i) {
        for (int j = 0; j <= maxK; ++j) vorFileEJK << setw(20) << left << vorE[i][j];
        vorFileEJK<<endl;
    }
    for(int i=0; i<=maxK; ++i) {
        for (int j = 0; j <= maxK; ++j) {
            if (q[i] > 0.0) vorFileEJK << setw(20) << left << vorE[i][j] / q[i];
            else vorFileEJK << setw(20) << left << 0.0;
        }
        vorFileEJK<<endl;
    }
    vorFileEJK<<endl;
    vorFilePKA<<endl;
    vorFilePKB<<endl;
    vorFilePKC<<endl;
    vorFileARA<<endl;
    vorFileARB<<endl;
    vorFileELA<<endl;
    vorFileELB<<endl;
    vorFileLAN<<setw(20)<<left<<vorLEA[0];
    vorFileLAN<<setw(20)<<left<<vorLEA[1]<<endl;
    vorFileLAN<<setw(20)<<left<<vorLEB[0];
    vorFileLAN<<setw(20)<<left<<vorLEB[1]<<endl;
    vorFileLAN<<setw(20)<<left<<vorANA[0];
    vorFileLAN<<setw(20)<<left<<vorANA[1]<<endl;
    vorFileLAN<<setw(20)<<left<<vorANB[0];
    vorFileLAN<<setw(20)<<left<<vorANB[1]<<endl;
    vorFileLAN<<setw(20)<<left<<cellLen;
    vorFileLAN<<setw(20)<<left<<cellLen*cellLen/nC<<endl;
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

    //Remove any .tmp files
    system("rm -f *.tmp");
}

void Configuration::setRadical(Logfile &logfile) {
    //Initialise ring distributions for type a, b and total

    //Maximum size arbitrary
    int maxK = 20;
    radPKA = VecF<double>(maxK+1);
    radPKB = VecF<double>(maxK+1);
    radPKC = VecF<double>(maxK+1);
    radARA = VecF<double>(maxK+1);
    radARB = VecF<double>(maxK+1);
    radELA = VecF<double>(maxK+1);
    radELB = VecF<double>(maxK+1);
    radE = VecF< VecF<double> >(maxK+1);
    radLEA = VecF<double>(3);
    radLEB = VecF<double>(3);
    radANA = VecF<double>(3);
    radANB = VecF<double>(3);
    radPKA = 0;
    radPKB = 0;
    radPKC = 0;
    radARA = 0;
    radARB = 0;
    radELA = 0;
    radELB = 0;
    radLEA = 0;
    radLEB = 0;
    radANA = 0;
    radANB = 0;
    for(int i=0; i<=maxK; ++i){
        radE[i] = VecF<double>(maxK+1);
        radE[i] = 0;
    }
}

void Configuration::radical(ofstream &radFilePKA, ofstream &radFilePKB, ofstream &radFilePKC, ofstream &radFileEJK,
                            ofstream &radFileARA, ofstream &radFileARB, ofstream &radFileELA, ofstream &radFileELB,
                            ofstream &radFileNet, Logfile &logfile) {
    //Radical Voronoi analysis

    //Initialise radical Voronoi and add configuration coordinates
    voro::container_poly rad3D(-cellLen_2,cellLen_2,-cellLen_2,cellLen_2,-0.5,0.5,10,10,1,true,true,false,nC);
    for(int i=0; i<nA; ++i) rad3D.put(i,xA[i],yA[i],0.0,wA);
    for(int i=0, j=nA; i<nB; ++i,++j) rad3D.put(j,xB[i],yB[i],0.0,wB);

    //Calculate Voronoi and generate temporary files
    rad3D.print_custom("%i","./radI.tmp");
    rad3D.print_custom("%a","./radK.tmp");
    rad3D.print_custom("%f","./radA.tmp");
    rad3D.print_custom("%e","./radE.tmp");
    rad3D.print_custom("%n","./radN.tmp");
    rad3D.print_custom("%l","./radV.tmp");
    rad3D.print_custom("%t","./radT.tmp");
    rad3D.print_custom("%p","./radP.tmp");

    //Read into 2D Voronoi and extract unnormalised distributions
    int maxK = 20;
    VecF<double> pA,pB,pC,aA,aB,eA,eB,lA,lB,anA,anB;
    VecF< VecF<double> > e;
    VoronoiBinary2D rad2D(nA,nB,maxK,true,logfile);
    rad2D.getDistributions(pA,pB,pC,aA,aB,eA,eB,e,lA,lB,anA,anB);

    //Add to running total distributions
    radPKA += pA;
    radPKB += pB;
    radPKC += pC;
    radARA += aA;
    radARB += aB;
    radELA += eA;
    radELB += eB;
    radLEA += lA;
    radLEB += lB;
    radANA += anA;
    radANB += anB;
    for(int i=0; i<=maxK; ++i) radE[i] += e[i];

    //Normalise distributions
    VecF<double> k(maxK+1);
    for(int i=0; i<=maxK; ++i) k[i]=i;
    for(int i=0; i<=maxK; ++i){
        if(pA[i]>0){
            aA[i] /= pA[i];
            eA[i] /= pA[i];
        }
        if(pB[i]>0){
            aB[i] /= pB[i];
            eB[i] /= pB[i];
        }
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
        radFilePKA<<setw(20)<<left<<pA[i];
        radFilePKB<<setw(20)<<left<<pB[i];
        radFilePKC<<setw(20)<<left<<pC[i];
        radFileARA<<setw(20)<<left<<aA[i];
        radFileARB<<setw(20)<<left<<aB[i];
        radFileELA<<setw(20)<<left<<eA[i];
        radFileELB<<setw(20)<<left<<eB[i];
    }
    radFilePKA<<endl;
    radFilePKB<<endl;
    radFilePKC<<endl;
    radFileARA<<endl;
    radFileARB<<endl;
    radFileELA<<endl;
    radFileELB<<endl;
    radFileNet<<setw(20)<<left<<meanA;
    radFileNet<<setw(20)<<left<<meanB;
    radFileNet<<setw(20)<<left<<meanC;
    radFileNet<<setw(20)<<left<<varA;
    radFileNet<<setw(20)<<left<<varB;
    radFileNet<<setw(20)<<left<<varC;
    radFileNet<<setw(20)<<left<<assortativity;
    radFileNet<<setw(20)<<left<<awA;
    radFileNet<<setw(20)<<left<<awVar;
    radFileNet<<setw(20)<<left<<awRSq;
    radFileNet<<endl;
}

void Configuration::radicalFinalise(ofstream &radFilePKA, ofstream &radFilePKB, ofstream &radFilePKC,
                                    ofstream &radFileEJK, ofstream &radFileARA, ofstream &radFileARB,
                                    ofstream &radFileELA, ofstream &radFileELB, ofstream &radFileLAN,
                                    ofstream &radFileNet, Logfile &logfile) {
    //Analyse distributions from all frames

    //Normalise distributions
    int maxK = 20;
    VecF<double> k(maxK+1);
    for(int i=0; i<=maxK; ++i) k[i]=i;
    for(int i=0; i<=maxK; ++i){
        if(radPKA[i]>0){
            radARA[i] /= radPKA[i];
            radELA[i] /= radPKA[i];
        }
        if(radPKB[i]>0){
            radARB[i] /= radPKB[i];
            radELB[i] /= radPKB[i];
        }
    }
    if(nA>0){
        radPKA /= nA*nCrdSets;
        radLEA[0] /= radLEA[2];
        radLEA[1] /= radLEA[2];
        radLEA[1] -= radLEA[0]*radLEA[0];
        radANA[0] /= radANA[2];
        radANA[1] /= radANA[2];
        radANA[1] -= radANA[0]*radANA[0];
    }
    if(nB>0){
        radPKB /= nB*nCrdSets;
        radLEB[0] /= radLEB[2];
        radLEB[1] /= radLEB[2];
        radLEB[1] -= radLEB[0]*radLEB[0];
        radANB[0] /= radANB[2];
        radANB[1] /= radANB[2];
        radANB[1] -= radANB[0]*radANB[0];
    }
    radPKC /= nC*nCrdSets;
    VecF<double> q(maxK+1);
    for(int i=0; i<=maxK; ++i) q[i] = vSum(radE[i]);
    double normE = vSum(q);
    for(int i=0; i<=maxK; ++i) radE[i] /= normE;
    q /= normE;

    //Calculate distirbution metrics
    double meanA,meanB,meanC,varA,varB,varC;
    meanA = vSum(k*radPKA);
    meanB = vSum(k*radPKB);
    meanC = vSum(k*radPKC);
    varA = vSum(k*k*radPKA)-meanA*meanA;
    varB = vSum(k*k*radPKB)-meanB*meanB;
    varC = vSum(k*k*radPKC)-meanC*meanC;
    double assortativity=0.0;
    for(int i=0; i<=maxK; ++i){
        for(int j=0; j<=maxK; ++j){
            assortativity += i*j*(radE[i][j]-q[i]*q[j]);
        }
    }
    assortativity /= vSum(k*k*q)-pow(vSum(k*q),2);
    VecR<double> awX(0,maxK+1),awY(0,maxK+1);
    for(int i=0; i<=maxK; ++i){
        if(radPKC[i]>0){
            awX.addValue(6.0*(i-6.0));
            awY.addValue(i*vSum(k*radE[i])/q[i]);
        }
    }
    VecR<double> aw=vLinearRegression(awX,awY);
    double awA = 1.0-aw[0];
    double awVar = aw[1]-36.0;
    double awRSq = aw[2];

    //Write to files
    for(int i=0; i<=maxK; ++i){
        radFilePKA<<setw(20)<<left<<radPKA[i];
        radFilePKB<<setw(20)<<left<<radPKB[i];
        radFilePKC<<setw(20)<<left<<radPKC[i];
        radFileARA<<setw(20)<<left<<radARA[i];
        radFileARB<<setw(20)<<left<<radARB[i];
        radFileELA<<setw(20)<<left<<radELA[i];
        radFileELB<<setw(20)<<left<<radELB[i];
    }
    for(int i=0; i<=maxK; ++i) {
        for (int j = 0; j <= maxK; ++j) radFileEJK << setw(20) << left << radE[i][j];
        radFileEJK<<endl;
    }
    for(int i=0; i<=maxK; ++i) {
        for (int j = 0; j <= maxK; ++j) {
            if (q[i] > 0.0) radFileEJK << setw(20) << left << radE[i][j] / q[i];
            else radFileEJK << setw(20) << left << 0.0;
        }
        radFileEJK<<endl;
    }
    radFileEJK<<endl;
    radFilePKA<<endl;
    radFilePKB<<endl;
    radFilePKC<<endl;
    radFileARA<<endl;
    radFileARB<<endl;
    radFileELA<<endl;
    radFileELB<<endl;
    radFileLAN<<setw(20)<<left<<radLEA[0];
    radFileLAN<<setw(20)<<left<<radLEA[1]<<endl;
    radFileLAN<<setw(20)<<left<<radLEB[0];
    radFileLAN<<setw(20)<<left<<radLEB[1]<<endl;
    radFileLAN<<setw(20)<<left<<radANA[0];
    radFileLAN<<setw(20)<<left<<radANA[1]<<endl;
    radFileLAN<<setw(20)<<left<<radANB[0];
    radFileLAN<<setw(20)<<left<<radANB[1]<<endl;
    radFileLAN<<setw(20)<<left<<cellLen;
    radFileLAN<<setw(20)<<left<<cellLen*cellLen/nC<<endl;
    radFileNet<<setw(20)<<left<<meanA;
    radFileNet<<setw(20)<<left<<meanB;
    radFileNet<<setw(20)<<left<<meanC;
    radFileNet<<setw(20)<<left<<varA;
    radFileNet<<setw(20)<<left<<varB;
    radFileNet<<setw(20)<<left<<varC;
    radFileNet<<setw(20)<<left<<assortativity;
    radFileNet<<setw(20)<<left<<awA;
    radFileNet<<setw(20)<<left<<awVar;
    radFileNet<<setw(20)<<left<<awRSq;
    radFileNet<<endl;

    //Remove any .tmp files
    system("rm -f *.tmp");
}

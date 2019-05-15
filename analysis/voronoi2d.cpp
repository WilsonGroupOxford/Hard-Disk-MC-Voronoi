#include "voronoi2d.h"

VoronoiBinary2D::VoronoiBinary2D() {}

VoronoiBinary2D::VoronoiBinary2D(int numA, int numB, int maxS, bool radical, Logfile &logfile) {
    //Convert output files from voro++ to 2D voronoi/delaunay triangulation

    //Assign
    nA = numA;
    nB = numB;

    //Calculate additional
    nC = nA+nB;
    maxSize = maxS;
    sizeA = VecF<int>(nA);
    sizeB = VecF<int>(nB);
    areaA = VecF<double>(nA);
    areaB = VecF<double>(nB);
    edgeA = VecF<double>(nA);
    edgeB = VecF<double>(nB);
    nbListAB = VecF< VecF<int> >(nC);
    lengthA = VecF<double>(3); //x,x^2 and nx
    lengthB = VecF<double>(3);
    angleA = VecF<double>(3);
    angleB = VecF<double>(3);
    lengthA = 0.0;
    lengthB = 0.0;
    angleA = 0.0;
    angleB = 0.0;

    //Set up filenames
    string filenameI,filenameV,filenameK,filenameA,filenameE,filenameN,filenameT,filenameP;
    if(!radical){
        filenameI = "./vorI.tmp";
        filenameV = "./vorV.tmp";
        filenameK = "./vorK.tmp";
        filenameA = "./vorA.tmp";
        filenameE = "./vorE.tmp";
        filenameN = "./vorN.tmp";
        filenameT = "./vorT.tmp";
        filenameP = "./vorP.tmp";
    }
    else{
        filenameI = "./radI.tmp";
        filenameV = "./radV.tmp";
        filenameK = "./radK.tmp";
        filenameA = "./radA.tmp";
        filenameE = "./radE.tmp";
        filenameN = "./radN.tmp";
        filenameT = "./radT.tmp";
        filenameP = "./radP.tmp";
    }

    //Read id data from .tmp files as faces not necessarily in the correct order
    string line;
    VecF<int> ids(nC);
    ifstream idFile(filenameI, ios::in);
    for(int i=0; i<nC; ++i){
        getline(idFile,line);
        istringstream(line)>>ids[i];
    }

    //To convert from 3D->2D need to retain only faces with normal as z-axis
    double tol=1e-6;
    ifstream vecFile(filenameV, ios::in);
    ifstream sizeFile(filenameK, ios::in);
    ifstream areaFile(filenameA, ios::in);
    ifstream edgeFile(filenameE, ios::in);
    ifstream nlFile(filenameN, ios::in);
    ifstream vertexFile(filenameT, ios::in);
    ifstream vertexCrdFile(filenameP, ios::in);

    for (int i=0; i<nC; ++i){
        //Extract normal vectors as (x,y,z) and convert to x y z
        getline(vecFile,line);
        int j = 0;
        while(j<line.size()){
            if(line[j] == '(' || line[j] == ')') line.erase(j,1);
            else if(line[j] == ','){
                line.replace(j,1," ");
                ++j;
            }
            else ++j;
        }
        istringstream ss(line);
        double normalComponent;
        int xyz=0; //0=x,1=y,2=z
        int pos=0;
        VecR<int> facePos2D(0,20);
        //Save position of faces with nz as -1 or 1
        while (ss >> normalComponent){
            if(xyz==2){
                if(fabs(fabs(normalComponent)-1.0)<tol) {
                    facePos2D.addValue(pos);
                }
                else if(fabs(normalComponent)>tol){
                    logfile.criticalError("Cannot identify rings in 2D Voronoi");
                }
                xyz=0;
                ++pos;
            }
            else ++xyz;
        }

        //Extract ring sizes, double check they agree
        int k,k0,k1;
        pos = 0;
        getline(sizeFile,line);
        istringstream ssk(line);
        while(ssk >> k){
            if(pos==facePos2D[0]) k0=k;
            else if(pos==facePos2D[1]) k1=k;
            ++pos;
        }
        if(k0!=k1) logfile.criticalError("Error in conversion to 2D Voronoi");

        //Extract ring areas, double check they agree
        double a,a0,a1;
        pos = 0;
        getline(areaFile,line);
        istringstream ssa(line);
        while(ssa >> a){
            if(pos==facePos2D[0]) a0=a;
            else if(pos==facePos2D[1]) a1=a;
            ++pos;
        }
        if(fabs(a0-a1)>tol) logfile.criticalError("Error in conversion to 2D Voronoi");

        //Extract ring perimeters, double check they agree
        double e,e0,e1;
        pos = 0;
        getline(edgeFile,line);
        istringstream ssp(line);
        while(ssp >> e){
            if(pos==facePos2D[0]) e0=e;
            else if(pos==facePos2D[1]) e1=e;
            ++pos;
        }
        if(fabs(e0-e1)>tol) logfile.criticalError("Error in conversion to 2D Voronoi");

        //Extract neighbour list
        VecF<int> nl(k0);
        int id;
        pos = 0;
        getline(nlFile,line);
        istringstream ssn(line);
        while(ssn >> id){
            if(id>=0){
                nl[pos] = id;
                ++pos;
            }
        }
        if(pos!=k0) logfile.criticalError("Error in conversion to 2D Voronoi");

        //Extract vertex path
        getline(vertexFile,line);
        pos = -1;
        j = 0;
        string vertexPos;
        while(j<line.size()){
            if(line[j] == '(') ++pos;
            if(pos==facePos2D[0]){
                k = j+1;
                vertexPos="";
                for(;;){
                    if(line[k] == ')') break;
                    else if(line[k] != ',') vertexPos+=line[k];
                    else if(line[k] == ',') vertexPos+=" ";
                    ++k;
                }
                break;
            }
            ++j;
        }
        istringstream vss(vertexPos);
        VecF<int> vPath(k0);
        pos=0;
        while(vss>>k){
            vPath[pos]=k;
            ++pos;
        }
        //Extract vertex path coordinates
        getline(vertexCrdFile,line);
        j = 0;
        while(j<line.size()){
            if(line[j] == '(' || line[j] == ')') line.erase(j,1);
            else if(line[j] == ','){
                line.replace(j,1," ");
                ++j;
            }
            else ++j;
        }
        VecR< VecF<double> > vertexCrds(k0*2);
        VecF<double> crd(2);
        double c;
        istringstream vCrd(line);
        xyz=0;
        pos=0;
        while(vCrd >> c){
            if(xyz%3==0) crd[0]=c;
            else if (xyz%3==1) crd[1]=c;
            else{
                vertexCrds[pos]=crd;
                ++pos;
            }
            ++xyz;
        }
        //Calculate lengths and angles
        VecF<double> v0,v1;
        VecR<double> angles(0,k0);
        VecR<double> lengths(0,k0);
        int l;
        for(int j=0; j<vPath.n; ++j){
            k=(j+1)%vPath.n;
            l=(j+2)%vPath.n;
            v0 = vertexCrds[vPath[j]] - vertexCrds[vPath[k]];
            v1 = vertexCrds[vPath[l]] - vertexCrds[vPath[k]];
            double n0,n1;
            double theta = vAngle(v0,v1,n0,n1);
            if(theta>-1) {//extremely occasional degenerate nodes where n0 or n1=0
                angles.addValue(theta);
                lengths.addValue(n0);
            }
        }

        //Store depending on type
        id = ids[i];
        if(id<nA){
            sizeA[id] = k0;
            areaA[id] = a0;
            edgeA[id] = e0/k0;
            for(int j=0; j<angles.n; ++j){
                lengthA[0] += lengths[j];
                lengthA[1] += lengths[j]*lengths[j];
                ++lengthA[2];
                angleA[0] += angles[j];
                angleA[1] += angles[j]*angles[j];
                ++angleA[2];
            }
        }
        else{
            sizeB[id-nA] = k0;
            areaB[id-nA] = a0;
            edgeB[id-nA] = e0/k0;
            for(int j=0; j<angles.n; ++j){
                lengthB[0] += lengths[j];
                lengthB[1] += lengths[j]*lengths[j];
                ++lengthB[2];
                angleB[0] += angles[j];
                angleB[1] += angles[j]*angles[j];
                ++angleB[2];
            }
        }
        nbListAB[id] = nl;
    }

    //Check for mutual connections in neighbour lists
    bool mutualCnx;
    for(int i=0; i<nC; ++i){
        VecF<int> nl = nbListAB[i];
        for(int j=0; j<nl.n; ++j){
            mutualCnx = vContains(nbListAB[nl[j]],i);
            if(!mutualCnx) logfile.criticalError("Error in conversion to 2D Voronoi");
        }
    }

    //Close files
    vecFile.close();
    sizeFile.close();
    areaFile.close();
    nlFile.close();

    //Calculate size, area and connection distirbutions
    calculateDistributions(logfile);
}

void VoronoiBinary2D::calculateDistributions(Logfile &logfile) {
    //Calculate unnormalised size and area distributions for different types and connection distribution

    //Set up vectors
    double tol=1e-6;
    VecF<double> k(maxSize+1);
    sizeDistA = VecF<double>(maxSize+1);
    sizeDistB = VecF<double>(maxSize+1);
    sizeDistC = VecF<double>(maxSize+1);
    areaDistA = VecF<double>(maxSize+1);
    areaDistB = VecF<double>(maxSize+1);
    edgeDistA = VecF<double>(maxSize+1);
    edgeDistB = VecF<double>(maxSize+1);
    cnxDist = VecF< VecF<double> >(maxSize+1);
    sizeDistA = 0;
    sizeDistB = 0;
    sizeDistC = 0;
    areaDistA = 0;
    areaDistB = 0;
    edgeDistA = 0;
    edgeDistB = 0;
    for(int i=0; i<=maxSize; ++i){
        cnxDist[i] = VecF<double>(maxSize+1);
        cnxDist[i] = 0;
    }

    //A size, area and edge length
    for(int i=0; i<nA; ++i){
        int size = sizeA[i];
        ++sizeDistA[size];
        areaDistA[size] += areaA[i];
        edgeDistA[size] += edgeA[i];
    }

    //B size, area and edge lengt
    for(int i=0; i<nB; ++i){
        int size = sizeB[i];
        ++sizeDistB[size];
        areaDistB[size] += areaB[i];
        edgeDistB[size] += edgeB[i];
    }

    //Combined size
    sizeDistC = sizeDistA+sizeDistB;

    //Connections
    for(int i=0; i<nA; ++i){
        int s0 = sizeA[i];
        VecF<int> nl=nbListAB[i];
        for(int j=0; j<s0; ++j){
            int k = nl[j];
            int s1;
            if(k<nA) s1 = sizeA[k];
            else s1 = sizeB[k-nA];
            ++cnxDist[s0][s1];
        }
    }
    for(int i=0; i<nB; ++i){
        int s0 = sizeB[i];
        VecF<int> nl=nbListAB[nA+i];
        for(int j=0; j<s0; ++j){
            int k = nl[j];
            int s1;
            if(k<nA) s1 = sizeA[k];
            else s1 = sizeB[k-nA];
            ++cnxDist[s0][s1];
        }
    }

    //Check connections symmetric
    bool symmetric = true;
    for(int i=0; i<=maxSize; ++i){
        for(int j=0; j<=maxSize; ++j){
            if(fabs(cnxDist[i][j]-cnxDist[j][i])>tol) symmetric=false;
            if(!symmetric) break;
        }
    }
    if(!symmetric) logfile.criticalError("Error in conversion to 2D Voronoi");
}

void VoronoiBinary2D::getDistributions(VecF<double> &sA, VecF<double> &sB, VecF<double> &sC, VecF<double> &aA,
                                       VecF<double> &aB, VecF<double> &eA, VecF<double> &eB, VecF<VecF<double> > &e,
                                       VecF<double> &lA, VecF<double> &lB, VecF<double> &anA, VecF<double> &anB) {
    //Get unnormalised distributions

    sA = sizeDistA;
    sB = sizeDistB;
    sC = sizeDistC;
    aA = areaDistA;
    aB = areaDistB;
    eA = edgeDistA;
    eB = edgeDistB;
    e = cnxDist;
    lA = lengthA;
    lB = lengthB;
    anA = angleA;
    anB = angleB;
}


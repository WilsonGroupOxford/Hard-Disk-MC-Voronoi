#include "voronoi2d.h"

VoronoiBinary2D::VoronoiBinary2D() {}

VoronoiBinary2D::VoronoiBinary2D(int numA, int numB, int maxS, Logfile &logfile) {
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
    nbListAB = VecF< VecF<int> >(nC);

    //Read id data from .tmp files as faces not necessarily in the correct order
    string line;
    VecF<int> ids(nC);
    ifstream idFile("./vorI.tmp", ios::in);
    for(int i=0; i<nC; ++i){
        getline(idFile,line);
        istringstream(line)>>ids[i];
    }

    //To convert from 3D->2D need to retain only faces with normal as z-axis
    double tol=1e-6;
    ifstream vecFile("./vorV.tmp", ios::in);
    ifstream sizeFile("./vorK.tmp", ios::in);
    ifstream areaFile("./vorA.tmp", ios::in);
    ifstream nlFile("./vorN.tmp", ios::in);
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

        //Store depending on type
        id = ids[i];
        if(id<nA){
            sizeA[id] = k0;
            areaA[id] = a0;
        }
        else{
            sizeB[id-nA] = k0;
            areaB[id-nA] = a0;
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
    cnxDist = VecF< VecF<double> >(maxSize+1);
    sizeDistA = 0;
    sizeDistB = 0;
    sizeDistC = 0;
    areaDistA = 0;
    areaDistB = 0;
    for(int i=0; i<=maxSize; ++i){
        cnxDist[i] = VecF<double>(maxSize+1);
        cnxDist[i] = 0;
    }

    //A size and area
    for(int i=0; i<nA; ++i){
        int size = sizeA[i];
        ++sizeDistA[size];
        areaDistA[size] += areaA[i];
    }

    //B size and area
    for(int i=0; i<nB; ++i){
        int size = sizeB[i];
        ++sizeDistB[size];
        areaDistB[size] += areaB[i];
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
                                       VecF<double> &aB, VecF<VecF<double> > &e) {
    //Get unnormalised distributions

    sA = sizeDistA;
    sB = sizeDistB;
    sC = sizeDistC;
    aA = areaDistA;
    aB = areaDistB;
    e = cnxDist;
}


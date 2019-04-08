#include "voronoi2d.h"

Voronoi2D::Voronoi2D() {}

Voronoi2D::Voronoi2D(int numA, int numB, int maxS, Logfile &logfile) {
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
        istringstream ss{regex_replace(line, regex{R"(\(|\)|,)"}, " ")};
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
        if(pos!=k1) logfile.criticalError("Error in conversion to 2D Voronoi");

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
}

//void Voronoi2D::networkAnalysis(int kMax, Logfile &logfile) {
//    //Network analysis of Delaunay simplex
//
//    //Calculate size and area distributions
//    double tol=1e-6;
//    VecF<double> k(kMax+1);
//    sizeDist = VecF<double>(kMax+1);
//    areaDist = VecF<double>(kMax+1);
//    for(int i=0; i<=kMax; ++i) k[i]=i;
//    sizeDist = 0;
//    areaDist = 0;
//    int size;
//    //Calculate raw distirbutions
//    for(int i=0; i<nF; ++i){
//        size = fSize[i];
//        ++sizeDist[size];
//        areaDist[size] += fArea[i];
//    }
//    //Normalise
//    for(int i=0; i<=kMax; ++i){
//        if(sizeDist[i]>0.0) areaDist[i] /= sizeDist[i];
//    }
//    sizeDist /= nF;
//    meanSize = vSum(k*sizeDist);
//    varSize = vSum(k*k*sizeDist)-meanSize*meanSize;
//
//    //Calculate assortativity
//    cnxDist = VecF< VecF<double> >(kMax+1);
//    for(int i=0; i<=kMax; ++i){
//        cnxDist[i] = VecF<double>(kMax+1);
//        cnxDist[i] = 0.0;
//    }
//    for(int i=0; i<nF; ++i){
//        int size0 = fSize[i];
//        VecF<int> nl = fNL[i];
//        for(int j=0; j<size0; ++j){
//            int size1 = fSize[nl[j]];
//            ++cnxDist[size0][size1];
//        }
//    }
//    //Calculate q, check symmetric and q matches p
//    VecF<double> q(kMax+1);
//    q=0.0;
//    for(int i=0; i<=kMax; ++i){
//        for(int j=0; j<=kMax; ++j){
//            q[i] += cnxDist[i][j];
//            if(fabs(cnxDist[i][j]-cnxDist[j][i])>tol) logfile.criticalError("Error in network analysis");
//        }
//    }
//    double norm = vSum(q);
//    q /= norm;
//    for(int i=0; i<=kMax; ++i) cnxDist[i] /= norm;
//    double qCheck=0.0;
//    for(int i=0; i<=kMax; ++i){
//        qCheck += fabs(q[i]-i*sizeDist[i]/meanSize);
//    }
//    if(qCheck>tol) logfile.criticalError("Error in network analysis");
//    //Calculate assortativity
//    assortativity = 0.0;
//    for(int i=0; i<=kMax; ++i){
//        for(int j=0; j<=kMax; ++j){
//            assortativity += i*j*(cnxDist[i][j]-q[i]*q[j]);
//        }
//    }
//    norm = vSum(k*k*q)-pow(vSum(k*q),2);
//    assortativity /= norm;
//}

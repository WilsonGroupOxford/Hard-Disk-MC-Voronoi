#include "voronoi.h"


Voronoi::Voronoi(VecF<double> &x, VecF<double> &y, VecF<double> &w, double cellLen_2, int numA, bool radical) {
    //Initialise with x,y coordinates and radii

    //Make periodic container in xy
    n=x.n;
    nA=numA;
    nB=n-nA;
    int blocks=sqrt(n);
    dz=cellLen_2*2;
    con=make_shared<voro::container_poly>(-cellLen_2,cellLen_2,-cellLen_2,cellLen_2,-cellLen_2,cellLen_2,
            blocks,blocks,1,true,true,false,n);

    //Add particles and radii if radical
    if(radical){
        for(int i=0; i<n; ++i) con->put(i,x[i],y[i],0.0,w[i]);
    }
    else{
        for(int i=0; i<n; ++i) con->put(i,x[i],y[i],0.0,0.0);
    }

//    voro::container testvor(-cellLen_2,cellLen_2,-cellLen_2,cellLen_2,-0.5,0.5,10,10,1,true,true,false,1);
//    for(int i=0; i<n; ++i) testvor.put(i,x[i],y[i],0.0);
//    testvor.print_custom("%t","./t.tmp");
}


void Voronoi::analyse(int maxSize, VecF<int> &cellSizeDistA, VecF<int> &cellSizeDistB, VecF<VecF<int> > &cellAdjDist,
                      VecF<double> &cellAreaA, VecF<double> &cellAreaB) {
    //Analyse Voronoi cell sizes and adjacencies

    //Compute cell neighbours
    computeCells(maxSize);

    //Calculate cell sizes and size distribution
    VecF<int> cellSizes(n);
    for(int i=0; i<n; ++i){
        cellSizes[i]=cellNbs[i].n;
    }
    cellSizeDistA=VecF<int>(maxSize);
    cellSizeDistB=VecF<int>(maxSize);
    for(int i=0; i<nA; ++i) ++cellSizeDistA[cellSizes[i]];
    for(int i=nA; i<n; ++i) ++cellSizeDistB[cellSizes[i]];

    //Calculate cell adjacencies distribution
    cellAdjDist=VecF< VecF<int> >(maxSize);
    for(int i=0; i<maxSize; ++i) cellAdjDist[i]=VecF<int>(maxSize);
    for(int i=0; i<n; ++i){
        int sizeI=cellSizes[i];
        for(int j=0; j<cellNbs[i].n; ++j){
            int sizeJ=cellSizes[cellNbs[i][j]];
            ++cellAdjDist[sizeI][sizeJ];
            //Account for self interactions
//            if(i==cellNbs[i][j]) ++cellAdjDist[sizeI][sizeJ];
        }
    }

    //Calculate cell size average areas
    cellAreaA=VecF<double>(maxSize+1);
    cellAreaB=VecF<double>(maxSize+1);
    for(int i=0; i<nA; ++i) cellAreaA[cellSizes[i]]+=cellAreas[i];
    for(int i=nA; i<n; ++i) cellAreaB[cellSizes[i]]+=cellAreas[i];
    cellAreaA[maxSize]=vSum(cellAreaA);
    cellAreaB[maxSize]=vSum(cellAreaB);
}


void Voronoi::nnDistances(VecF<double> &x, VecF<double> &y, double cellLen, double rCellLen, VecF<double> &nnSep, VecF<int> &nnCount) {
    //Calculate nearest neighbour distances

    //Calculate nearest neighbour distances between different cell types
    nnSep=VecF<double>(3); //A-A, A-B, B-B
    nnCount=VecF<int>(3);
    VecF<int> nType(n);
    for(int i=0; i<nA; ++i) nType[i]=0;
    for(int i=nA; i<n; ++i) nType[i]=1;
    double dx,dy,d;
    for(int i=0; i<n; ++i){
        for(int j=0; j<cellNbs[i].n; ++j){
            dx=x[i]-x[cellNbs[i][j]];
            dy=y[i]-y[cellNbs[i][j]];
            dx-=cellLen*nearbyint(dx*rCellLen);
            dy-=cellLen*nearbyint(dy*rCellLen);
            d=sqrt(dx*dx+dy*dy);
            ++nnCount[nType[i]+nType[cellNbs[i][j]]];
            nnSep[nType[i]+nType[cellNbs[i][j]]]+=d;
        }
    }
}


void Voronoi::computeCells(int maxSize) {
    //Calculate neighbouring particles for each particle, and area of cells

    //Make looper
    voro::c_loop_all looper(*con);
    looper.start();

    //Resize neighbour vectors
    cellNbs=VecF< VecR<int> >(n);
    for(int i=0; i<n; ++i) cellNbs[i]=VecR<int>(0,maxSize);
    cellAreas=VecF<double>(n);

    //Loop over each cell and extract neighbour ids
    do{
        int id=looper.pid(); //central id
        voro::voronoicell_neighbor cell;
        con->compute_cell(cell,looper);
        vector<int> nbs;
        cell.neighbors(nbs);
        for(int i=0; i<nbs.size(); ++i){
            cellNbs[id].addValue(nbs[i]);
        }
        cellNbs[id].delValue(-5); //remove z cell boundary
        cellNbs[id].delValue(-6); //remove z cell boundary
//        for(int i=0; i<cellNbs[id].n; ++i) if(cellNbs[id][i]<0) cellNbs[id][i]=id; //add self interaction
        cellAreas[id]=cell.volume()/dz;
    } while(looper.inc());
}


void Voronoi::getRings(VecF<double> &x, VecF<double> &y, VecF< VecR<double> > &rings) {
    //Find rings as vertex coordinates

    //Make looper
    voro::c_loop_all looper(*con);
    looper.start();

    //Resize vectors
    rings=VecF< VecR<double> >(cellNbs.n);
    for(int i=0; i<cellNbs.n; ++i) rings[i]=VecR<double>(0,2*(cellNbs[i].n+1));

    //Loop over each cell and extract rings
    do{
        int id=looper.pid(); //central id
        voro::voronoicell_neighbor cell;
        con->compute_cell(cell,looper);
        //Find faces, normals and vertices
        vector<int> faceSizes; //number of vertices in each face
        cell.face_orders(faceSizes);
        int numFaces=faceSizes.size(); //number of faces
        vector<double> normals; //normals for each face (x,y,z)
        cell.normals(normals);
        int keyFace=-1; //id of face of interest
        for(int i=0; i<numFaces; ++i){
            if(fabs(normals[3*i+2]-1)<1e-12){
                keyFace=i;
                break;
            }
        }
        vector<int> vertexIds; //ids of vertices that make up faces
        cell.face_vertices(vertexIds);
        vector<double> vertexCrds; //coordinates of vertices
        cell.vertices(vertexCrds);
        cell.vertices(x[id],y[id],0.0,vertexCrds);
        //Extract vertices for face
        VecR<int> keyVertexIds(0,100);
        int k=0;
        for(int i=0; i<numFaces; ++i){
            if(i==keyFace){
                for(int j=0; j<faceSizes[i]+1; ++j){
                    keyVertexIds.addValue(vertexIds[k]);
                    ++k;
                }
                break;
            }
            else k+=faceSizes[i]+1;
        }
        for(int i=1; i<keyVertexIds.n; ++i){//first index gives face size
            rings[id].addValue(vertexCrds[3*keyVertexIds[i]]);
            rings[id].addValue(vertexCrds[3*keyVertexIds[i]+1]);
        }
    } while(looper.inc());
}
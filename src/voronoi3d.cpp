#include "voronoi3d.h"


Voronoi3D::Voronoi3D(VecF<double> &x, VecF<double> &y, VecF<double> &z, VecF<double> &r, double cellLen_2, int numA, bool radical, int maxV) {
    //Initialise with x,y,z coordinates and radii

    //Make periodic container in xy
    n=x.n;
    nA=numA;
    nB=n-nA;
    maxVertices=maxV;
    int blocks=sqrt(n);
    dz=2*vMaximum(r);
    pbc=cellLen_2*2;
    rpbc=1.0/pbc;
    con=make_shared<voro::container_poly>(-cellLen_2,cellLen_2,-cellLen_2,cellLen_2,0,dz,
            blocks,blocks,1,true,true,false,n);

    //Add particles and radii if radical
    if(radical){
        for(int i=0; i<n; ++i) con->put(i,x[i],y[i],z[i],r[i]);
    }
    else{
        for(int i=0; i<n; ++i) con->put(i,x[i],y[i],z[i],0.0);
    }

    //Calculate cell projections
    computeCellProjections(x,y,z);
}


void Voronoi3D::analyse(int maxSize, VecF<int> &cellSizeDistA, VecF<int> &cellSizeDistB, VecF<VecF<int> > &cellAdjDist,
                      VecF<double> &cellAreaA, VecF<double> &cellAreaB) {
    //Analyse Voronoi cell sizes and adjacencies

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


void Voronoi3D::nnDistances(VecF<double> &x, VecF<double> &y, double cellLen, double rCellLen, VecF<double> &nnSep, VecF<int> &nnCount) {
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


void Voronoi3D::computeCellProjections(VecF<double> &x, VecF<double> &y, VecF<double> &z) {
    //Calculate cell projections onto lower wall
    //Calculate neighbouring particles for each particle, and area of cells

    //Make looper
    voro::c_loop_all looper(*con);
    looper.start();

    //Resize containers
    VecF< VecR<int> > cellNbs3D(n); //all 3D neighbours
    cellNbs3D=VecF< VecR<int> >(n);
    for(int i=0; i<n; ++i) cellNbs3D[i]=VecR<int>(0,12*maxVertices);
    ringCrds=VecF< VecR<double> >(cellNbs3D.n);
    for(int i=0; i<cellNbs3D.n; ++i) ringCrds[i]=VecR<double>(0,12*maxVertices);
    cellAreas=VecF<double>(n);
    faceCrds=VecR< VecR<double> >(0,n*20);

    //Loop over each cell and extract neighbour ids
    do{
        //Find all 3D neighbours
        int id=looper.pid(); //central id
        voro::voronoicell_neighbor cell;
        con->compute_cell(cell,looper);
        vector<int> nbs;
        cell.neighbors(nbs);
        for(int i=0; i<nbs.size(); ++i){
            cellNbs3D[id].addValue(nbs[i]);
        }
        if(vContains(cellNbs3D[id],-5)){
            cellNbs3D[id].delValue(-5); //remove z cell boundary
        }
        else{
            cout<<z[id]<<endl;
        }
        if(vContains(cellNbs3D[id],-6)) cellNbs3D[id].delValue(-6); //remove z cell boundary if present

        //Find faces, normals and vertices
        vector<int> faceSizes; //number of vertices in each face
        cell.face_orders(faceSizes);
        int numFaces=faceSizes.size(); //number of faces
        vector<double> normals; //normals for each face (x,y,z)
        cell.normals(normals);
        int keyFace=-1; //id of face on the lower plane
        for(int i=0; i<numFaces; ++i){
            if(fabs(normals[3*i+2]+1)<1e-12){
                keyFace=i;
                break;
            }
        }
        vector<int> vertexIds; //ids of vertices that make up faces
        cell.face_vertices(vertexIds);
        vector<double> vertexCrds; //coordinates of vertices
        cell.vertices(vertexCrds);
        cell.vertices(x[id], y[id], z[id], vertexCrds);

        //Store coordinates for all faces
        int k=0;
        for(int i=0; i<numFaces; ++i){
            VecR<double> fCrds(0,maxVertices*3);
            for(int j=0; j<faceSizes[i]+1;++j){
                if(j>0){
                    fCrds.addValue(vertexCrds[3*vertexIds[k]]);
                    fCrds.addValue(vertexCrds[3*vertexIds[k]+1]);
                    fCrds.addValue(vertexCrds[3*vertexIds[k]+2]);
                }
                ++k;
            }
            faceCrds.addValue(fCrds);
        }

        //Project coordinates of face on lower wall
        if(keyFace!=-1) {
            //Extract vertices for face
            VecR<int> keyVertexIds(0,2*maxVertices);
            k = 0;
            for (int i = 0; i < numFaces; ++i) {
                if (i == keyFace) {
                    for (int j = 0; j < faceSizes[i] + 1; ++j) {
                        keyVertexIds.addValue(vertexIds[k]);
                        ++k;
                    }
                    break;
                } else k += faceSizes[i] + 1;
            }
            for (int i = 1; i < keyVertexIds.n; ++i) {//first index gives face size
                ringCrds[id].addValue(vertexCrds[3 * keyVertexIds[i]]);
                ringCrds[id].addValue(vertexCrds[3 * keyVertexIds[i] + 1]);
            }
        }
    } while(looper.inc());

    //Find neighbours of projected face only
    cellNbs=VecF< VecR<int> >(n); //only projected neighbours
    double tol=1e-10; //tolerance for floating point comparison
    for(int i=0; i<n; ++i) cellNbs[i]=VecR<int>(0,maxVertices);
    for(int id=0; id<n; ++id){
        for(int i=0; i<cellNbs3D[id].n; ++i){
            bool neighbour=false;
            int nbId=cellNbs3D[id][i];
            for(int j=0; j<ringCrds[nbId].n/2; ++j){
                double cx0=ringCrds[nbId][2*j];
                double cy0=ringCrds[nbId][2*j+1];
                for(int k=0; k<ringCrds[id].n/2; ++k){
                    double cx1=ringCrds[id][2*k];
                    double cy1=ringCrds[id][2*k+1];
                    double dx=cx0-cx1;
                    double dy=cy0-cy1;
                    dx-=pbc*nearbyint(dx*rpbc);
                    dy-=pbc*nearbyint(dy*rpbc);
                    if(fabs(dx)<tol && fabs(dy)<tol){
                        neighbour=true;
                        break;
                    }
                }
                if(neighbour) break;
            }
            if(neighbour) cellNbs[id].addValue(nbId);
        }
        if(cellNbs[id].n!=ringCrds[id].n/2) cout<<"Error in 3D Voronoi neighbours"<<endl;
    }

    //Calculate cell areas using shoelace formula
    for(int id=0; id<n; ++id){
        int nn=cellNbs[id].n;
        for(int i=0,j=1; i<nn-1; ++i,++j){
            cellAreas[id]+=ringCrds[id][2*i]*ringCrds[id][2*j+1];
            cellAreas[id]-=ringCrds[id][2*j]*ringCrds[id][2*i+1];
        }
        cellAreas[id]+=ringCrds[id][2*(nn-1)]*ringCrds[id][1];
        cellAreas[id]-=ringCrds[id][0]*ringCrds[id][2*(nn-1)+1];
        cellAreas[id]=0.5*fabs(cellAreas[id]);
    }
}


VecF< VecR<double> > Voronoi3D::getProjectedRings() {
    //Get projected ring coordinates

    return ringCrds;
}


VecR< VecR<double> > Voronoi3D::getFaces(VecF<double> &zLimits) {
    //Get face coordinates

    zLimits=VecF<double>(2);
    zLimits[1]=dz;
    return faceCrds;
}
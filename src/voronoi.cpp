#include "voronoi.h"


Voronoi::Voronoi(VecF<double> &x, VecF<double> &y, VecF<double> &r, double cellLen_2, bool radical) {
    //Initialise with x,y coordinates and radii

    //Make periodic container in xy
    n=x.n;
    int blocks=sqrt(n);
    con=make_shared<voro::container_poly>(-cellLen_2,cellLen_2,-cellLen_2,cellLen_2,-cellLen_2,cellLen_2,
            blocks,blocks,1,true,true,false,n);

    //Add particles and radii if radical
    if(radical){
        for(int i=0; i<n; ++i) con->put(i,x[i],y[i],0.0,r[i]);
    }
    else{
        for(int i=0; i<n; ++i) con->put(i,x[i],y[i],0.0,0.0);
    }
}


void Voronoi::analyse(int maxSize, VecF<int> &cellSizeDist, VecF< VecF<int> > &cellAdjDist) {
    //Analyse Voronoi cell sizes and adjacencies

    //Compute cell neighbours
    computeNeighbours(maxSize);

    //Calculate cell sizes and size distribution
    VecF<int> cellSizes(n);
    for(int i=0; i<n; ++i){
        cellSizes[i]=cellNbs[i].n;
    }
    cellSizeDist=VecF<int>(maxSize);
    for(int i=0; i<n; ++i) ++cellSizeDist[cellSizes[i]];

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
}


void Voronoi::computeNeighbours(int maxSize) {
    //Calculate neighbouring particles for each particle

    //Make looper
    voro::c_loop_all looper(*con);
    looper.start();

    //Resize neighbour vectors
    cellNbs=VecF< VecR<int> >(n);
    for(int i=0; i<n; ++i) cellNbs[i]=VecR<int>(0,maxSize);

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
    } while(looper.inc());
}
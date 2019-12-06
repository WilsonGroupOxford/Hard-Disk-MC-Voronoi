#ifndef HDMC_VORONOI3D_H
#define HDMC_VORONOI3D_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <memory>
#include "outputfile.h"
#include "vecf.h"
#include "vecr.h"
#include "vec_func.h"
#include "../voro++/src/voro++.hh"

using namespace std;

class Voronoi3D {
    //Voronoi analysis of 3D system using Voro++

private:

    //Data members
    shared_ptr<voro::container_poly> con;
    int n,nA,nB; //total number of particles and of type A, B
    int maxVertices; //maximum number of vertices per cell
    double dz; //height of cells
    double pbc,rpbc; //periodic boundary conditions
    VecF< VecR<int> > cellNbs; //neighbours of each cell
    VecF<double> cellAreas; //areas of each cell
    VecF< VecR<double> > ringCrds; //coordinates of rings
    VecR< VecR<double> > faceCrds; //coordinates of all faces

    //Member functions
    void computeCellProjections(VecF<double> &x, VecF<double> &y, VecF<double> &z); //find neighbours for each cell projection

public:

    //Constructor
    Voronoi3D(VecF<double> &x, VecF<double> &y, VecF<double> &z, VecF<double> &r, double cellLen_2, int numA, bool radical, int maxV); //3D coordinates and radii, cell info

    //Member functions
    void analyse(int maxSize, VecF<int> &cellSizeDistA, VecF<int> &cellSizeDistB, VecF< VecF<int> > &cellAdjDist, VecF<double> &cellAreaA, VecF<double> &cellAreaB);
    void nnDistances(VecF<double> &x, VecF<double> &y, double cellLen, double rCellLen, VecF<double> &nnSep, VecF<int> &nnCount);
    VecF< VecR<double> > getProjectedRings();
    VecR< VecR<double> > getFaces(VecF<double> &zLimits);
};


#endif //HDMC_VORONOI3D_H

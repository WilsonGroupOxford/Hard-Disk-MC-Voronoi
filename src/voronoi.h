#ifndef HDMC_VORONOI_H
#define HDMC_VORONOI_H

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

class Voronoi {
    //Voronoi analysis of 2D system using Voro++

private:

    //Data members
    shared_ptr<voro::container_poly> con;
    int n,nA,nB; //total number of particles and of type A, B
    double dz; //height of cells
    VecF< VecR<int> > cellNbs; //neighbours of each cell
    VecF<double> cellAreas; //areas of each cell

    //Member functions
    void computeCells(int maxSize); //find neighbours for each cell

public:

    //Constructor
    Voronoi(VecF<double> &x, VecF<double> &y, VecF<double> &w, double cellLen_2, int numA, bool radical); //2D coordinates and weights, cell info

    //Member functions
    void analyse(int maxSize, VecF<int> &cellSizeDistA, VecF<int> &cellSizeDistB, VecF< VecF<int> > &cellAdjDist, VecF<double> &cellAreaA, VecF<double> &cellAreaB);
    void nnDistances(VecF<double> &x, VecF<double> &y, double cellLen, double rCellLen, VecF<double> &nnSep, VecF<int> &nnCount);
    void getRings(VecF<double> &x, VecF<double> &y, VecF< VecR<double> > &rings);
};


#endif //HDMC_VORONOI_H

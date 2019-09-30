#ifndef HDMC_VORONOI_H
#define HDMC_VORONOI_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
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
    int n; //number of particles
    VecF< VecR<int> > cellNbs; //neighbours of each cell

    //Member functions
    void computeNeighbours(int maxSize); //find neighbours for each cell

public:

    //Constructor
    Voronoi(VecF<double> &x, VecF<double> &y, VecF<double> &r, double cellLen_2, bool radical); //2D coordinates and radii, cell info

    //Member functions
    void analyse(int maxSize, VecF<int> &cellSizeDist, VecF< VecF<int> > &cellAdjDist);
};


#endif //HDMC_VORONOI_H

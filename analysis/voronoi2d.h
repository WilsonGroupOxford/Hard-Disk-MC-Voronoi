//Read in outputs of Voro++, convert to 2D and analyse

#ifndef BNAHDMC_ANALYSIS_VORONOI2D_H
#define BNAHDMC_ANALYSIS_VORONOI2D_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <regex>
#include <cmath>
#include "outputfile.h"
#include "vecf.h"
#include "vecr.h"
#include "vec_func.h"

using namespace std;

class VoronoiBinary2D {
    //2D binary Voronoi diagram converted from voro++ output files

private:

    //Voronoi variables
    int nA,nB,nC; //number of type a,b and total particles
    int maxSize; //maximum ring sizes to include
    VecF<int> sizeA, sizeB; //sizes of type a and b
    VecF<double> areaA, areaB; //areas of rings of each type
    VecF< VecF<int> > nbListAB; //neighbour list for both type a and b

    //Analysis variables
    VecF<double> sizeDistA, sizeDistB, sizeDistC; //unnormalised ring size distributions
    VecF<double> areaDistA, areaDistB; //unnormalised ring area distributions
    VecF< VecF<double> > cnxDist; //unnormalised distribution of edges in delaunnay

    //Member functions
    void calculateDistributions(Logfile &logfile); //calculate size,area and cnx distributions

public:

    //Constructors
    VoronoiBinary2D();
    VoronoiBinary2D(int numA, int numB, int maxS, Logfile &logfile);

    //Getters
    void getDistributions(VecF<double> &sA, VecF<double> &sB, VecF<double> &sC,
            VecF<double> &aA, VecF<double> &aB, VecF< VecF<double> > &e);
};


#endif //BNAHDMC_ANALYSIS_VORONOI2D_H

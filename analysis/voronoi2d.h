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

class Voronoi2D {
    //2D binary Voronoi diagram converted from voro++ output files

private:

    //Voronoi variables
    int nA,nB,nC; //number of type a,b and total particles
    int maxSize; //maximum ring sizes to include
    VecF<int> sizeA, sizeB; //sizes of type a and b
    VecF<double> areaA, areaB; //areas of rings of each type
    VecF< VecF<int> > nbListAB; //neighbour list for both type a and b





    //Analysis variables
    double meanSize,varSize; //mean and variance of face distribution
    double assortativity; //assortative mixing of Delaunnay
    VecF<double> sizeDist, areaDist; //distribution of sizes and areas
    VecF< VecF<double> > cnxDist; //distribution of edges in Delaunnay

public:

    //Constructors
    Voronoi2D();
    Voronoi2D(int numA, int numB, int maxS, Logfile &logfile);
    void networkAnalysis(int kMax, Logfile &logfile); //calculate network properties


};


#endif //BNAHDMC_ANALYSIS_VORONOI2D_H

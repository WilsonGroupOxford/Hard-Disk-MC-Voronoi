#include "configuration.h"

Configuration::Configuration() {}

Configuration::Configuration(int numA, int numB, double radA, double radB, double cellLength) {
    //Construct with configuration information

    //Assign
    nA = numA;
    nB = numB;
    rA = radA;
    rB = radB;
    cellLen = cellLength;

    //Additional
    nCrdSets = 0;
    nC = nA+nB;
    rCellLen = 1.0/cellLen;
    cellLen_2 = cellLen/2;
    xA = VecF<double>(nA);
    yA = VecF<double>(nA);
    xB = VecF<double>(nB);
    yB = VecF<double>(nB);
}

void Configuration::setCoordinates(ifstream& xyzFile, Logfile& logfile) {
    //Read next configuration from XYZ file

    //Load configuration
    logfile.write("Configuration loaded: ",nCrdSets);
    string skip,line;
    istringstream ss("");
    getline(xyzFile,skip);
    getline(xyzFile,skip);

    //Read type a
    for(int i=0; i<nA; ++i){
        getline(xyzFile,line);
        ss.str(line);
        ss>>skip;
        ss>>xA[i];
        ss>>yA[i];
    }

    //Read type b
    for(int i=0; i<nB; ++i){
        getline(xyzFile,line);
        ss.str(line);
        ss>>skip;
        ss>>xB[i];
        ss>>yB[i];
    }

    //Increment total coordinate sets added
    ++nCrdSets;
}

void Configuration::voronoi(Logfile& logfile) {
    //Voronoi analysis

    //Initialise container
    logfile.write("Voronoi");
    voro::container voronoi(0.0,cellLen,0.0,cellLen,-0.5,0.5,1,1,1,true,true,false,nC);

    //Add coordinates
    for(int i=0; i<nA; ++i) voronoi.put(i,xA[i],yA[i],0.0);
    for(int i=0, j=nA; i<nB; ++i,++j) voronoi.put(j,xB[i],yB[i],0.0);

    //Generate Voronoi
    cout<<voronoi.sum_cell_volumes()<<" "<<cellLen*cellLen<<endl;
}
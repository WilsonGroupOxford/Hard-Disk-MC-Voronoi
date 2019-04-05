#include <iostream>
#include <sstream>
#include "outputfile.h"

using namespace std;

int main() {

    //Set up logfile
    Logfile logfile("./mc.log");
    logfile.datetime("Simulation begun at: ");
    logfile.write("Binary Non-Additive Hard Sphere Monte Carlo");
    logfile.write("Written By: David OM, Wilson Group, 2019");
    logfile.separator();



    //Close files
    logfile.datetime("Simulation complete at: ");
    return 0;
}
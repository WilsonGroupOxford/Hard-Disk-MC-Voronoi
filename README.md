## Hard Disk Monte Carlo

Hard disk simulation code orignally designed to study colloidal monolayers.

It has the following features:

* Simulation of mono- or bi-disperse hard disk systems
* Additive or non-additive interactions
* On-the-fly Voronoi and structural analysis

### Requirements

* Voronoi analysis is performed by the excellent library [Voro++](http://math.lbl.gov/voro++/about.html). 
This should be downloaded, compiled and placed in the directory above src.

### Compilation

The code can be compiled using CMake. 

* Ensure the path to the Voro++ library is correct in ```CMakeLists.txt```
* If ```CMakeCache.txt``` exists remove this file
* Run ```cmake .```
* Edit the ```CMakeCache.txt``` file with the following:
```
CMAKE_BUILD_TYPE:STRING=Release
CMAKE_CXX_FLAGS:STRING=-std=c++11
```
* Run ```make```

This should generate the executable ```hdmc.x``` . 

### Input

Simulation parameters are controlled through the input file ```hdmc.inpt```,
which should be located in the same directory as the executable.

### Running

The code can be run with ```./hdmc.x```. 

As the simulation progresses a log file ```hdmc.log``` will be written 
containing simulation parameters and progress.

### Output

The following outputs will be produced if selected in the input file:

* RDF data is contained in ```rdf.dat```.
For monodisperse systems this is the distance and total RDF.
For bidisperse systems this is the distance, total RDF and partial RDFS in the order 1-1, 1-2, 2-2.
* Voronoi analysis is contained in ```vor.dat```. 
For monodisperse systems each line gives the ring statistics for a given configuration, with the assortativity in the final column.
For bidisperse systems lines alternate between each partial type giving the partial ring statistics and overall assortativity.
The final line(s) in each case give the results averaged over all configurations.
* Radical tessellation analysis is contained in ```rad.dat```. 
It follows the same format as the Voronoi analysis.


 


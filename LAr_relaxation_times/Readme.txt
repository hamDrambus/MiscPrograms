Dependencies:
1) c++ Boost. Headers are enough, no need for compiling it.

Required for the project overall but not used in compilation of C++ code:

=================================================
To compile the code run cmake to create makefile:

mkdir Build
cd Build
cmake ../
make

Or run script to create Eclipse project:

bash GenEclipseProject.sh

Open in eclipse as File->Improt->General->Existing project in folder->browse to generated -build folder. Build via Project Explorer->Build Targets. Debug as C/C++ remote Application, with set binary location and disabled auto build. When any files or libraris are added to/revoved from the project, it must be regenerated with GenEclipseProject.sh.

=================================================
This small program extracts data from 2D tabulated functions and prints 1D functions of interest to file
Parameters (input file, output destination and X or Y values at which 2D f is to be profiled) are set in setting.xml file
Run the program as:
...build/Geant_simulation path/to/settings.xml | tee path/to/log.txt

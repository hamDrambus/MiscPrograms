Dependencies:
1) c++ Boost. Headers are enough, no need for compiling it.

Required for the project overall but not used in compilation of C++ code:
Files with transport parameters (f, w drift) from LAr NBrS geant4 project.

=================================================
To compile the code run cmake to create makefile:

mkdir Build
cd Build
cmake ../
make

Or run .sh script to compile and run:

./CompileAndRun.sh

=================================================
This small program calculates relaxation times and distances for electron distribution function in liquid argon depending on the electric field.
The maximum electric field gradient (dE/dx) to maintain quasi-stationary f(e) is also estimated.
dE/dx real must be << dE/dx from this program.

Run the program as:
...build/LAr_XS_recalculation path/to/settings.xml | tee path/to/log.txt
Or
./Run.sh

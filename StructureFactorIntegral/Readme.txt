Dependencies:
1) c++ Boost. Headers are enough, no need for compiling it.

Required for the project overall but not used in compilation of C++ code:

=================================================
To compile the code run cmake to create makefile:

mkdir Build
cd Build
cmake ../
make

Or run script to compile and run:

bash CompileAndRun.sh

=================================================
This small program calculates S(e) (Atrazhev1985 https://doi.org/10.1088/0022-3719/18/6/015 page 1206) using tabulated structure factor S(k)

Run the program as:
...build/Geant_simulation path/to/settings.xml | tee path/to/log.txt

/*  Global parameters such as filenames, detector dimensions, etc., shared between all classes.
 *  The difference from GlobalData is that this class only stores simple parameters which all
 *  may be easily changed without silent breaking of simulation.
 *  All parameters must be read only during simulation after their initialization.
 */
#ifndef GlobalParameters_h
#define GlobalParameters_h
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <deque>

#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string.hpp>

#include "utilities/PolynomialFit.hh"

//#define TEMP_CODE_

namespace gPars
{
  struct ProgramSetups {
    std::string this_path;
    std::string settings_filename;
    std::string gnuplot_bin;
    std::string input_filename;
    std::string output_filename_template_X;
    std::string output_filename_template_Y;
    std::string output_X_title;
    std::string output_Y_title;
    std::vector<double> Xs_to_print; //Not scaled.
    std::vector<double> Ys_to_print; //Not scaled.
    double X_input_scale; //when reading data, it is multiplied (scaled) by this value
    double Y_input_scale;
    double Z_input_scale;
    std::vector<std::string> output_filenames_X;
    std::vector<std::string> output_filenames_Y;
  };

	bool InitGlobals(std::string filename);
	bool LoadSettings(std::string filename);

	extern ProgramSetups general;
}

#endif //GlobalParameters_h

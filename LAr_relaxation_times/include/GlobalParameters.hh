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

    std::string input_f_distr;
    std::string input_v_drift;
    std::string input_XS_energy_transfer;
    std::string input_XS_momentum_transfer;

    std::string output_folder;
    std::string output_tau_relax;
    std::string output_tau_relax_N;
    std::string output_max_E_gradient;
    std::string output_max_E_gradient_N;
    std::string output_l_relax;
    std::string output_l_relax_N;
  };

	bool InitGlobals(std::string filename);
	bool LoadSettings(std::string filename);

	extern ProgramSetups general;
}

#endif //GlobalParameters_h

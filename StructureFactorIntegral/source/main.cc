
#include <GlobalParameters.hh>
#include <utilities/FunctionTable.hh>

constexpr double joule = 1.0;
constexpr double eV = joule * 1.60217663e-19;
constexpr double m = 1.0;
constexpr double mm = 0.001*m;
constexpr double Angstrem = 1e-10*m;
constexpr double second = 1.0;
constexpr double hbar = 1.054571817e-34 * joule * second;
constexpr double e_mass = 9.1093837e-31;

int main(int argc, char** argv)
{
  std::string settings_fname;
  if (argc != 1) {
    if (argc>2)
      std::cout << "Warning! Only single parameter (settings filename) is used." << std::endl;
    settings_fname = argv[1];
  } else {
    std::cout << "Using default settings file \"settings.xml\"" << std::endl;
    settings_fname = "settings.xml";
  }
	if (!gPars::InitGlobals(settings_fname)) {
	  std::cerr<<"Failed to initialize globals."<<std::endl;
	  return -1;
	}

  DataVector StructureFactor;
  std::ifstream str;
  std::string in_filename = gPars::general.this_path + gPars::general.input_filename;
  str.open(in_filename);
  if (!str.is_open()) {
    std::cout << "Failed to open data file \""<<in_filename<<"\"" << std::endl;
    return -2;
  }
  StructureFactor.read(str);
  str.close();
  StructureFactor.scaleXY(gPars::general.X_input_scale, gPars::general.Y_input_scale); // must be in SI now

  IntegrationRange energy_range = IntegrationInterval(0*eV, 36 * eV, 0.02 * eV);
  IntegrationRange x_range = IntegrationInterval(0.0, 2.0, 1e-4);
  DataVector result;
  for (long int i = 0, i_end_ = energy_range.NumOfIndices(); i!=i_end_; ++i) {
    double energy_SI = energy_range.Value(i);
    double energy_eV = energy_SI / eV;
    double integral = 0;
    double Y_prev = 0;
    double Y;
    double X_prev = x_range.Value(0);
    double X;
    for (long int j = 0, j_end_ = x_range.NumOfIndices(); j!=j_end_; ++j) {
      X = x_range.Value(j);
      Y = 0.5 * X * StructureFactor( 2 * sqrt(e_mass * energy_SI * X) / hbar);
      integral += 0.5*(Y + Y_prev) * (X - X_prev);
      X_prev = X;
      Y_prev = Y;
    }
    result.insert(energy_eV, integral); // integral is dimensionless
  }
  result.write(gPars::general.output_filename, "Energy[eV]\tS(eps)");
	return 0;
}

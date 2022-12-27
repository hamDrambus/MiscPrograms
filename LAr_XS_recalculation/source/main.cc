
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
  std::string in_filename_S = gPars::general.this_path + gPars::general.input_S_filename;
  str.open(in_filename_S);
  if (!str.is_open()) {
    std::cout << "Failed to open data file \""<<in_filename_S<<"\"" << std::endl;
    return -2;
  }
  StructureFactor.read(str);
  str.close();

  DataVector inXS;
  std::string in_filename_XS = gPars::general.this_path + gPars::general.input_XS_filename;
  str.open(in_filename_XS);
  if (!str.is_open()) {
    std::cout << "Failed to open data file \""<<in_filename_XS<<"\"" << std::endl;
    return -2;
  }
  inXS.read(str);
  str.close();

  IntegrationRange x_range = IntegrationInterval(0.0, 1.0, 1e-3);
  x_range += IntegrationInterval(1.0, 10.0, 0.1);
  x_range += IntegrationInterval(10.0, 50.0, 0.5);
  DataVector outXS(inXS.getOrder(), inXS.getNused());
  outXS.use_leftmost(true); outXS.use_rightmost(true);
  for (std::size_t i = 0, i_end_ = x_range.NumOfIndices(); i!=i_end_; ++i) {
    double energy = x_range.Value(i);
    outXS.push_back(energy, inXS(energy)
       * (gPars::general.input_is_effective_XS ? StructureFactor(energy) : 1 / StructureFactor(energy)));
  }
  outXS.write(gPars::general.this_path + gPars::general.output_XS_filename, "Energy[eV]\tXS[cm2]");
  std::cout<<"Finished."<<std::endl;
	return 0;
}

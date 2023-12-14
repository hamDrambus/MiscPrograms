#include "GlobalParameters.hh"

namespace gPars
{
  ProgramSetups general;

  bool LoadSettings(std::string fname)
  {
    std::cout << "Loading settings \"" << fname << "\"..." << std::endl;
    // Create an empty property tree object
    using boost::property_tree::ptree;
    using boost::property_tree::ptree_bad_data;
    using boost::property_tree::ptree_bad_path;
    using boost::property_tree::ptree_error;
    ptree pt;
    try {
      read_xml(fname, pt);
      {
        ptree gen = pt.get_child("Settings.General");
        
        general.input_f_distr = gen.get<std::string>("input_f_distr");
        general.input_v_drift = gen.get<std::string>("input_v_drift");
        general.input_XS_energy_transfer = gen.get<std::string>("input_XS_energy_transfer");
        general.input_XS_momentum_transfer = gen.get<std::string>("input_XS_momentum_transfer");

        general.output_folder = gen.get<std::string>("output_folder", "");
        general.output_tau_relax = gen.get<std::string>("output_tau_relax", "tau_relax.txt");
        general.output_tau_relax_N = gen.get<std::string>("output_tau_relax_N", "tau_relax_N.txt");
        general.output_max_E_gradient = gen.get<std::string>("output_max_E_gradient", "max_E_gradient.txt");
        general.output_max_E_gradient_N = gen.get<std::string>("output_max_E_gradient_N", "max_E_gradient_N.txt");
        general.output_l_relax = gen.get<std::string>("output_l_relax", "output_l_relax.txt");
        general.output_l_relax_N = gen.get<std::string>("output_l_relax_N", "output_l_relax_N.txt");
      }
    } catch (ptree_bad_path& e) {
      std::cout << "LoadSettings: ptree_bad_path exception:" << std::endl;
      std::cout << e.what() << std::endl;
      goto fail_load;
    } catch (ptree_bad_data& e) {
      std::cout << "LoadSettings: ptree_bad_data exception:" << std::endl;
      std::cout << e.what() << std::endl;
      goto fail_load;
    } catch (ptree_error& e) {
      std::cout << "LoadSettings: ptree_error exception:" << std::endl;
      std::cout << e.what() << std::endl;
      goto fail_load;
    } catch (std::exception& e) {
      std::cout << "LoadSettings: std::exception:" << std::endl;
      std::cout << e.what() << std::endl;
      goto fail_load;
    }
    std::cout << "LoadSettings: Successfully loaded \"" << fname << "\"" << std::endl;
    return true;
  fail_load:
    std::cout << "LoadSettings: Failed to load settings \"" << fname << "\"" << std::endl;
    return false;
  }

  bool InitGlobals(std::string filename)
  {
    // Possible defaults are set in LoadSettings
    // They may be overwritten in LoadSettings().
    // If default is not possible and is absent in settings,
    // initialization fails.
    char path[FILENAME_MAX];
#if defined(__WIN32__)
    general.this_path = _getcwd(path, FILENAME_MAX);
#else
    general.this_path = getcwd(path, FILENAME_MAX);
#endif
    if (!general.this_path.empty())
        if (general.this_path.back()!='/')
          general.this_path.push_back('/');

    general.settings_filename = filename;
    if (!LoadSettings(filename))
      return false;

    //=====================================================================================================================
    // Consistency checks below

    std::cout<<"This path: \""<<general.this_path<<"\""<<std::endl;
    return true;
  }
}

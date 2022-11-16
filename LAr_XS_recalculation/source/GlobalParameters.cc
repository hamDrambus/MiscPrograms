#include <GlobalParameters.hh>

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
      general.input_XS_filename = gen.get<std::string>("input_XS_filename");
			general.input_S_filename = gen.get<std::string>("input_S_filename");
			general.output_XS_filename = gen.get<std::string>("output_XS_filename");
	    general.input_is_effective_XS = gen.get<bool>("input_is_effective_XS");
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

#if defined(_WIN32)||defined(_WIN64)
		general.gnuplot_bin = "\"%GNUPLOT%\\gnuplot.exe\"";
#else
		general.gnuplot_bin = "gnuplot";
#endif
		general.settings_filename = filename;
		if (!LoadSettings(filename))
		  return false;

		//=====================================================================================================================
		// Consistency checks below

    std::cout<<"This path: \""<<general.this_path<<"\""<<std::endl;
    std::cout<<"Input XS file: \""<<general.this_path+general.input_XS_filename<<"\""<<std::endl;
		std::cout<<"Input structure factor file: \""<<general.this_path+general.input_S_filename<<"\""<<std::endl;
		std::cout<<"Output XS file: \""<<general.this_path+general.output_XS_filename<<"\""<<std::endl;
    return true;
	}
}

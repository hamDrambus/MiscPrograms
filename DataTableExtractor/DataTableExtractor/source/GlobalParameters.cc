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
      general.input_filename = gen.get<std::string>("input_filename");
			general.output_filename_template_X = gen.get<std::string>("output_filename_template_X");
	    general.output_filename_template_Y = gen.get<std::string>("output_filename_template_Y");
	    general.output_X_title = gen.get<std::string>("output_X_title", "");
	    general.output_Y_title = gen.get<std::string>("output_Y_title", "");
	    general.X_input_scale = gen.get<double>("X_input_scale", 1.0);
	    general.Y_input_scale = gen.get<double>("Y_input_scale", 1.0);
	    general.Z_input_scale = gen.get<double>("Z_input_scale", 1.0);
			BOOST_FOREACH(ptree::value_type &w, gen.get_child("Plots")) {
				if (w.first == "X") {
					double X_val = w.second.get_value<double>();
					std::string x_str = w.second.get_value<std::string>();
					std::map<std::string, std::string> run_envvars;
					run_envvars["($X)"] = x_str;
					run_envvars["($Y)"] = "";
					std::string out_fname = general.output_filename_template_X;
					for (auto mm = run_envvars.begin(), mm_end_ = run_envvars.end(); mm != mm_end_; ++mm) {
						boost::replace_all(out_fname, mm->first, mm->second);
					}
					general.Xs_to_print.push_back(X_val);
					general.output_filenames_X.push_back(out_fname);
					continue;
				}
				if (w.first == "Y") {
					double Y_val = w.second.get_value<double>();
					std::string y_str = w.second.get_value<std::string>();
					std::map<std::string, std::string> run_envvars;
					run_envvars["($Y)"] = y_str;
					run_envvars["($X)"] = "";
					std::string out_fname = general.output_filename_template_Y;
					for (auto mm = run_envvars.begin(), mm_end_ = run_envvars.end(); mm != mm_end_; ++mm) {
						boost::replace_all(out_fname, mm->first, mm->second);
					}
					general.Ys_to_print.push_back(Y_val);
					general.output_filenames_Y.push_back(out_fname);
					continue;
				}
				std::cout << "Settings::Load: Warning: Unknown key \"" << w.first << "\" in \"Settings.General.Plots\" is ignored." << std::endl;
			}
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
    std::cout<<"Input file: \""<<general.this_path+general.input_filename<<"\""<<std::endl;
    return true;
	}
}

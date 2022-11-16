
#include <GlobalParameters.hh>
#include <utilities/FunctionTable.hh>

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

  FunctionTable data;
  std::ifstream str;
  std::string in_filename = gPars::general.this_path + gPars::general.input_filename;
  str.open(in_filename, std::ios_base::binary);
  if (!str.is_open()) {
    std::cout << "Failed to open data file \""<<in_filename<<"\"" << std::endl;
    return -2;
  }
  data.read(str);
  str.close();
  for (std::size_t i = 0, i_end_ = gPars::general.Xs_to_print.size(); i!=i_end_; ++i) {
    // Profiling 2D data table at X. Y values are derived from both DataVectors surrounding X.
    double X = gPars::general.Xs_to_print[i] / gPars::general.X_input_scale;
    boost::optional<std::pair<std::size_t, std::size_t>> ind = data.getX_indices(X);
    if (boost::none == ind)
      continue;
    if (data.getX(ind->first) > X || data.getX(ind->second) < X)
      continue;
    if (ind->first == ind->second) {
      DataVector profile = data[ind->first];
      profile.scaleXY(gPars::general.Y_input_scale, gPars::general.Z_input_scale);
      profile.write(gPars::general.output_filenames_X[i], gPars::general.output_X_title);
      continue;
    }
    DataVector profile1 = data[ind->first], profile2 = data[ind->second];
    DataVector profile;
    double X1 = data.getX(ind->first), X2 = data.getX(ind->second);
    // Deriving Y values at X from two DataVectors at X1 < X < X2.
    double Y, Y1, Y2, Y1_nxt, Y2_nxt, Y_nxt1, Y_nxt2;
    for (std::size_t j1 = 0, j2 = 0, j1_end_ = profile1.size(), j2_end_ = profile2.size(); j1!=j1_end_ && j2!=j2_end_;) {
      Y1 = profile1[j1].first;
      Y2 = profile2[j2].first;
      Y = (X2==X1 ? 0.5*(Y1+Y2) : Y1 + (Y2-Y1)*(X-X1)/(X2-X1));
      profile.insert(Y, data(X, Y));
      Y1_nxt = ((j1+1)!=j1_end_ ? profile1[j1 + 1].first : DBL_MAX);
      Y2_nxt = ((j2+1)!=j2_end_ ? profile1[j2 + 1].first : DBL_MAX);
      if (Y1_nxt == DBL_MAX && Y2_nxt == DBL_MAX)
        break;
      if (Y1_nxt != DBL_MAX)
        Y_nxt1 = (X2==X1 ? 0.5*(Y1_nxt+Y2) : Y1_nxt + (Y2-Y1_nxt)*(X-X1)/(X2-X1));
      else {
        ++j2;
        //profile.insert(Y2, data(X, Y2));
        continue;
      }
      if (Y2_nxt != DBL_MAX)
        Y_nxt2 = (X2==X1 ? 0.5*(Y2_nxt+Y1) : Y1 + (Y2_nxt-Y1)*(X-X1)/(X2-X1));
      else {
        ++j1;
        //profile.insert(Y1, data(X, Y1));
        continue;
      }
      if (Y_nxt2 < Y_nxt1) {
        ++j2;
        //profile.insert(Y2, data(X, Y2));
      } else {
        ++j1;
        //profile.insert(Y1, data(X, Y1));
      }
    }
    profile.scaleXY(gPars::general.Y_input_scale, gPars::general.Z_input_scale);
    profile.write(gPars::general.output_filenames_X[i], gPars::general.output_X_title);
  }

  for (std::size_t i = 0, i_end_ = gPars::general.Xs_to_print.size(); i!=i_end_; ++i) {
    // Profiling 2D data table at Y. Y values are derived from both data points surrounding Y (similar to the case above).
    std::vector<std::pair<double, double>> lower_Y, higher_Y;
    double Y = gPars::general.Ys_to_print[i] / gPars::general.Y_input_scale;
    for (std::size_t j = 0, j_end_=data.size(); j!=j_end_; ++j) {
      DataVector profileX = data[j];
      boost::optional<std::pair<std::size_t, std::size_t>> ind = profileX.getY_indices(Y);
      if (boost::none == ind)
        continue;
      if (profileX.getX(ind->first) > Y) {
        higher_Y.push_back(std::pair<double, double>(data.getX(j), profileX.getY(ind->first)));
        continue;
      }
      if (profileX.getX(ind->second) < Y) {
        lower_Y.push_back(std::pair<double, double>(data.getX(j), profileX.getY(ind->second)));
        continue;
      }
      lower_Y.push_back(std::pair<double, double>(data.getX(j), profileX.getY(ind->first)));
      higher_Y.push_back(std::pair<double, double>(data.getX(j), profileX.getY(ind->second)));
    }
    DataVector profile;
    double Y1, Y2, Y1_nxt, Y2_nxt;
    // Deriving X values at Y from two point arrays.
    double X, X1, X2, X1_nxt, X2_nxt, X_nxt1, X_nxt2;
    for (std::size_t j1 = 0, j2 = 0, j1_end_ = lower_Y.size(), j2_end_ = higher_Y.size(); j1!=j1_end_ && j2!=j2_end_;) {
      X1 = lower_Y[j1].first, X2 = higher_Y[j2].first;
      Y1 = lower_Y[j1].second; Y2 = higher_Y[j2].second;
      X = (Y2==Y1 ? 0.5*(X1+X2) : X1 + (X2-X1)*(Y-Y1)/(Y2-Y1));
      profile.insert(X, data(X, Y));
      X1_nxt = ((j1+1)!=j1_end_ ? lower_Y[j1 + 1].first : DBL_MAX);
      X2_nxt = ((j2+1)!=j2_end_ ? higher_Y[j2 + 1].first : DBL_MAX);
      Y1_nxt = ((j1+1)!=j1_end_ ? lower_Y[j1 + 1].second : DBL_MAX);
      Y2_nxt = ((j2+1)!=j2_end_ ? higher_Y[j2 + 1].second : DBL_MAX);
      if (X1_nxt == DBL_MAX && X2_nxt == DBL_MAX)
        break;
      if (X1_nxt != DBL_MAX)
        X_nxt1 = (Y2==Y1_nxt ? 0.5*(X1_nxt+X2) : X1_nxt + (X2-X1_nxt)*(Y-Y1_nxt)/(Y2-Y1_nxt));
      else {
        ++j2;
        profile.insert(X2, data(X2, Y));
        continue;
      }
      if (X2_nxt != DBL_MAX)
        X_nxt1 = (Y2_nxt==Y1 ? 0.5*(X1+X2_nxt) : X1 + (X2_nxt-X1)*(Y-Y1)/(Y2_nxt-Y1));
      else {
        ++j1;
        profile.insert(X1, data(X1, Y));
        continue;
      }
      if (X_nxt2 < X_nxt1) {
        ++j2;
        profile.insert(X2, data(X2, Y));
      } else {
        ++j1;
        profile.insert(X1, data(X1, Y));
      }
    }
    profile.scaleXY(gPars::general.X_input_scale, gPars::general.Z_input_scale);
    profile.write(gPars::general.output_filenames_Y[i], gPars::general.output_Y_title);
  }
	return 0;
}

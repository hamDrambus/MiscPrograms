#include <exception>
#include <vector>
#include <GlobalParameters.hh>
#include "utilities/FunctionTable.hh"
#include "utilities/GlobalUtilities.hh"

const double m_to_M = 1.373247e-5; // For argon M = 39.948.
const double atomic_density = 2.10e28; // For argon at 87 K, in m^-3
const double e_charge = 1.602176e-19;
const double eV_to_vel = sqrt(2*e_charge/9.109383e-31); // ev_to_vel * sqrt(kinetic_energy_in_eV) === velocity in SI
const double kT_eV = 0.00752293;// kT in eVs

struct equilibrium_data {
  FunctionTable f_distr;
  DataVector v_drift;
  DataVector XS_energy;
  DataVector XS_momentum;
  std::vector<double> fields;
  equilibrium_data(void)
  {
    std::ifstream str;
    std::string in_filename = gPars::general.this_path + gPars::general.input_f_distr;
    str.open(in_filename, std::ios_base::binary);
    if (!str.is_open())
      throw std::runtime_error(std::string("Failed to open data file \"")+in_filename+"\"");
    f_distr.read(str);
    str.close();
    f_distr.scale(1e4, 1, 1); // Convert field from geant4 units to kV/cm
    // Convert f' to simple f. f' = f(e)/sqrt(e). Int_-inf^+inf f(e)de = 1.
    for (std::size_t i = 0, i_end_ = f_distr.size(); i!=i_end_; ++i) {
      DataVector& f = f_distr[i];
      for (std::size_t e = 0, e_end_ = f.size(); e!=e_end_; ++e) {
        f[e].second = f[e].second * sqrt(f[e].first);
      }
    }
    fields.resize(f_distr.size());
    for (std::size_t i = 0, i_end_ = f_distr.size(); i!=i_end_; ++i)
      fields[i] = f_distr.getX(i);

    in_filename = gPars::general.this_path + gPars::general.input_v_drift;
    str.open(in_filename);
    if (!str.is_open())
      throw std::runtime_error(std::string("Failed to open data file \"")+in_filename+"\"");
    v_drift.read(str);
    str.close();
    v_drift.scaleXY(1e4, 1e-3*1e9); // Convert from geant4 units to [m/s] (SI) vs kV/cm
    if (v_drift.size() != fields.size())
      throw std::runtime_error(std::string("Fields mismatch between f and v drift."));
    for (std::size_t i = 0, i_end_ = v_drift.size(); i!=i_end_; ++i)
      if (v_drift.getX(i)!=fields[i])
        throw std::runtime_error(std::string("Fields mismatch between f and v drift, i = ") + std::to_string(i));

    in_filename = gPars::general.this_path + gPars::general.input_XS_energy_transfer;
    str.open(in_filename);
    if (!str.is_open())
      throw std::runtime_error(std::string("Failed to open data file \"")+in_filename+"\"");
    XS_energy.read(str);
    XS_energy.scaleY(1e-4); // Convert to m^2 vs eV.
    str.close();

    in_filename = gPars::general.this_path + gPars::general.input_XS_energy_transfer;
    str.open(in_filename);
    if (!str.is_open())
      throw std::runtime_error(std::string("Failed to open data file \"")+in_filename+"\"");
    XS_momentum.read(str);
    XS_momentum.scaleY(1e-4); // Convert to m^2 vs eV.
    str.close();
  }
};

// Unstable
double get_deriv(const DataVector& data, const std::size_t i)
{
  std::size_t sz = data.size();
  if (sz < 3)
    return 0;
  if (i >= sz)
    return 0;
  if (0 == i)
    return (data.getY(1) - data.getY(0))/(data.getX(1) - data.getX(0));
  if (sz == (i+1))
    return (data.getY(i) - data.getY(i-1))/(data.getX(i) - data.getX(i-1));
  return 0.5 * ((data.getY(i) - data.getY(i-1))/(data.getX(i) - data.getX(i-1))
          + (data.getY(i+1) - data.getY(i))/(data.getX(i+1) - data.getX(i)));
}

// Unstable
double get_deriv2(const DataVector& data, const std::size_t i)
{
  std::size_t sz = data.size();
  if (sz < 3)
    return 0;
  if (i >= sz)
    return 0;
  if (0 == i)
    return get_deriv2(data, i+1);
  if (sz == (i+1))
    return get_deriv2(data, i-1);
  return ((data.getY(i+1) - data.getY(i))/(data.getX(i+1) - data.getX(i)) - 
          (data.getY(i) - data.getY(i-1))/(data.getX(i) - data.getX(i-1))) / (data.getX(i+1) - data.getX(i-1));
}

double get_deriv(const DataVector& data, const double energy)
{
  std::size_t sz = data.size();
  double de = std::max(1e-3, energy * 0.01);
  double e_limit = data.getX(data.size()-1);
  double e_min = energy < de ? 0 : (energy > e_limit - de ? e_limit - de : energy - de);
  double e_max = energy < de ? de : (energy > e_limit - de ? e_limit : energy + de);
  return (data(e_max) - data(e_min))/(e_max - e_min);
}

double get_deriv2(const DataVector& data, const double energy)
{
  std::size_t sz = data.size();
  double de = std::max(2e-3, energy * 0.02);
  double e_limit = data.getX(data.size()-1);
  double e_min = energy < de ? 0 : (energy > e_limit - de ? e_limit - de : energy - de);
  double e_max = energy < de ? de : (energy > e_limit - de ? e_limit : energy + de);
  return (get_deriv(data, e_max) - get_deriv(data, e_min))/(e_max - e_min);
}

// Get derivative of z=data(field, energy) over field (x) at point energy (y).
double get_deriv(const FunctionTable& data, const std::size_t i_field, double energy)
{
  std::size_t sz = data.size();
  if (sz < 3)
    return DBL_MIN;
  if (i_field >= sz)
    return 0;
  if (0 == i_field)
    return get_deriv(data, i_field+1, energy);
  if (sz == (i_field+1))
    return get_deriv(data, i_field-1, energy);
  const DataVector& f0 = data[i_field - 1];
  const DataVector& f1 = data[i_field];
  const DataVector& f2 = data[i_field + 1];
  double E0 = data.getX(i_field - 1);
  double E1 = data.getX(i_field);
  double E2 = data.getX(i_field + 1);
  double out = 0.5*((f1(energy) - f0(energy))/(E1 - E0) + (f2(energy) - f1(energy))/(E2 - E1));
  return std::fabs(out) < 1e-100 ? 1e-100 : out;
}

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

  equilibrium_data eq_data;
  std::vector<FILE*> plot_pipes;
  
  std::size_t fsz = eq_data.fields.size();
  std::vector<double> out_tau0(fsz), out_tau1(fsz), out_tau(fsz),
        out_l(fsz), out_grad(fsz), out_N_mt(fsz), out_N_et(fsz),
        out_grad_N(fsz), out_l_N(fsz), out_tau_N(fsz);
  for (std::size_t i = 0, i_end_ = fsz; i!=i_end_; ++i) {
    double E = eq_data.fields[i];
    double E_SI = E * 1e5;
    std::string title = "E=" + std::to_string(E)+"kV/cm";
    DataVector& f = eq_data.f_distr[i];
    std::size_t sz = 200; // f.size();
    std::vector<double> es(sz), f_vals(sz), taus0(sz), taus0_rhs(sz), taus0_lhs(sz), taus1(sz);
    std::vector<double> max_grad(sz), l_relax(sz);
    std::vector<double> deriv(sz);
    double en_max = f.getX(f.size()-1), en_min = f.getX(0);
    double en_prev = en_min;
    double tau0_avg = 0, tau1_avg = 0, tau_avg = 0, d_log_f_avg = 0, en_avg = 0;
    double N_mt_freq_avg = 0, N_et_freq_avg = 0; // frequency of momentum- and energy-transfer collisions.
    for (std::size_t e = 0, e_end_ = sz; e!=e_end_; ++e) {
      double en = en_min + e * (en_max - en_min) / e_end_; // energy
      es[e] = en;
      f_vals[e] = f(en);
      double xs_mt = eq_data.XS_momentum(en);
      double xs_mt_deriv = get_deriv(eq_data.XS_momentum, en);
      double xs_et = eq_data.XS_energy(en);
      double xs_et_deriv = get_deriv(eq_data.XS_energy, en);

      double f0 = f(en);
      double f1 = get_deriv(f, en);
      double f2 = get_deriv2(f, en);

      // Calculating tau1 from right-hand side of Boltzmann equation
      double inv_tau1 = (eV_to_vel * sqrt(en) * atomic_density * xs_mt);
      taus1[e] = 1.0 / inv_tau1; // In SI

      // Calculating tau0 from left-hand side of Boltzmann equation (term with E)
      double inv_tau0 = std::pow(e_charge * E_SI, 2) * sqrt(en) * eV_to_vel / f0 / (atomic_density) / 3.0;
      inv_tau0 *= (f1 / en / xs_mt
                  - f1*xs_mt_deriv/std::pow(xs_mt, 2)
                  + f2/xs_mt);
      inv_tau0 /= e_charge*e_charge; // (1/eV)s from f to SI.
      taus0_lhs[e] = 1.0 / std::fabs(inv_tau0); // In SI

      // Calculating tau0 from right-hand side of Boltzmann equation
      inv_tau0 = 2 * m_to_M * sqrt(en) * eV_to_vel * atomic_density / f0;
      inv_tau0 *= (2*(f0 + 2 * kT_eV * f1) * xs_et
                  + en*(f1 + 2 * kT_eV *f2) * xs_et
                  + en*(f0 + 2 * kT_eV * f1) * xs_et_deriv);
      inv_tau0 = std::fabs(inv_tau0);
      inv_tau0 = f0 < 1e-10 ? 0 : inv_tau0;
      taus0_rhs[e] = 1.0 / inv_tau0; // In SI
      taus0[e] = taus0_rhs[e];

      double tau_max = std::max(taus0[e], taus1[e]);
      double inv_tau_min = std::min(inv_tau0, inv_tau1);
      l_relax[e] = tau_max * eq_data.v_drift.getY(i); // in m
      double inv_l_relax = inv_tau_min / eq_data.v_drift.getY(i); // in m
      deriv[e] = std::fabs(get_deriv(eq_data.f_distr, i, en));
      double d_log_f = std::fabs(get_deriv(eq_data.f_distr, i, en)) / f0;
      d_log_f = f0 < 1e-10 ? 0 : d_log_f;
      max_grad[e] = d_log_f < 1e-10 ? 0 : 1.0 / (l_relax[e] * 1e2) / d_log_f; // in kV/cm/cm

      // Taking avegare values (integral with f(e)de).
      tau0_avg += inv_tau0 * f0 * (en-en_prev);
      tau1_avg += inv_tau1 * f0 * (en-en_prev);
      tau_avg += inv_tau_min * f0 * (en-en_prev);
      d_log_f_avg += d_log_f * f0 * (en-en_prev);
      en_avg += en * f0 * (en-en_prev);
      N_mt_freq_avg += eV_to_vel * sqrt(en) * atomic_density * xs_mt * f0 * (en-en_prev);
      N_et_freq_avg += eV_to_vel * sqrt(en) * atomic_density * 2 * m_to_M * xs_et * f0 * (en-en_prev);

      en_prev = en;
    }

    out_tau_N[i] = 1.0 / N_et_freq_avg;
    out_l_N[i] = out_tau_N[i] * eq_data.v_drift.getY(i); // in m
    out_tau0[i] = 1.0/tau0_avg;
    out_tau1[i] = 1.0/tau1_avg;
    out_tau[i] = 1.0/tau_avg;
    out_l[i] = out_tau[i] * eq_data.v_drift.getY(i);
    out_grad[i] = 1.0/d_log_f_avg/(out_l[i]*1e2); // in kV/cm/cm
    out_N_mt[i] = N_mt_freq_avg * out_tau[i];
    out_N_et[i] = N_et_freq_avg * out_tau[i];
    out_grad_N[i] = 1.0/d_log_f_avg/(out_l_N[i]*1e2); // in kV/cm/cm
    out_l_N[i] = out_l_N[i];

    if (i == 10000) {
      plot_pipes.push_back(plot(es, f_vals, title, "e (eV)", "f(e)"));
      plot_pipes.push_back(plot(es, taus0, title, "e (eV)", "tau0(e)"));
      plot_pipes.push_back(plot(es, taus0_rhs, title, "e (eV)", "tau0 rhs(e)"));
      plot_pipes.push_back(plot(es, taus0_lhs, title, "e (eV)", "tau0 lhs(e)"));
      plot_pipes.push_back(plot(es, taus1, title, "e (eV)", "tau1(e)"));
      plot_pipes.push_back(plot(es, deriv, title, "e (eV)", "df/dE (cm/kV/eV)"));
      plot_pipes.push_back(plot(es, max_grad, title, "e (eV)", "max dE/dx (kV/cm^2)"));
    }
  }
  std::string x_label = "E (kV/cm)";
  plot_pipes.push_back(plot(eq_data.fields, out_tau0, "tau0(E)", x_label, "tau0(E)"));
  plot_pipes.push_back(plot(eq_data.fields, out_tau1, "tau1(E)", x_label, "tau1(E)"));
  plot_pipes.push_back(plot(eq_data.fields, out_tau, "tau(E)", x_label, "tau(E)"));
  plot_pipes.push_back(plot(eq_data.fields, out_l, "l relax.", x_label, "l relax.(E)"));
  plot_pipes.push_back(plot(eq_data.fields, out_l_N, "l relax. from N collisions", x_label, "l relax.(E)"));
  plot_pipes.push_back(plot(eq_data.fields, out_grad, "Max dE/dx", x_label, "Max dE/dx (kV/cm^2)"));
  plot_pipes.push_back(plot(eq_data.fields, out_grad_N, "Max dE/dx from N collisions", x_label, "Max dE/dx (kV/cm^2)"));
  plot_pipes.push_back(plot(eq_data.fields, out_N_mt, "N momentum-transfer collisions", x_label, "N momentum-transfer collisions"));
  plot_pipes.push_back(plot(eq_data.fields, out_N_et, "N energy-transfer collisions", x_label, "N energy-transfer collisions"));

  std::string fname = gPars::general.output_folder + gPars::general.output_l_relax;
  ensure_file(fname);
  DataVector temp = DataVector(eq_data.fields, out_l);
  temp.scaleY(1e6);
  temp.write(fname, "Estimation from Boltzmann equaion\n"
                    "//Field (kV/cm)\tRelaxation distance (um)");

  fname = gPars::general.output_folder + gPars::general.output_l_relax_N;
  ensure_file(fname);
  temp = DataVector(eq_data.fields, out_l_N);
  temp.scaleY(1e6);
  temp.write(fname, "Estimation from average N collisions = 1\n"
                    "//Field (kV/cm)\tRelaxation distance (um)");

  fname = gPars::general.output_folder + gPars::general.output_tau_relax;
  ensure_file(fname);
  temp = DataVector(eq_data.fields, out_tau);
  temp.scaleY(1e9);
  temp.write(fname, "Estimation from Boltzmann equaion\n"
                    "//Field (kV/cm)\tRelaxation time (ns)");

  fname = gPars::general.output_folder + gPars::general.output_tau_relax_N;
  ensure_file(fname);
  temp = DataVector(eq_data.fields, out_tau_N);
  temp.scaleY(1e9);
  temp.write(fname, "Estimation from average N collisions = 1\n"
                    "//Field (kV/cm)\tRelaxation time (ns)");

  fname = gPars::general.output_folder + gPars::general.output_max_E_gradient;
  ensure_file(fname);
  temp = DataVector(eq_data.fields, out_grad);
  temp.write(fname, "Estimation from Boltzmann equaion\n"
                    "//Field (kV/cm)\tMaximum field gradient (kV/cm/cm)");

  fname = gPars::general.output_folder + gPars::general.output_max_E_gradient_N;
  ensure_file(fname);
  temp = DataVector(eq_data.fields, out_grad_N);
  temp.write(fname, "Estimation from average N collisions = 1\n"
                    "//Field (kV/cm)\tMaximum field gradient (kV/cm/cm)");

  std::cin.clear();
  std::cin.ignore(std::cin.rdbuf()->in_avail());
  std::cin.get();

  for (FILE* pipe : plot_pipes)
    if (nullptr != pipe)
      pclose(pipe);
  return 0;
}

/*
A C++ source file for OpenMC to give neutron spectra produced via the 7Li(p,n)7Be reaction.

Author: Max Conroy (m.j.conroy@pgr.bham.ac.uk)
*/

#include <memory> // for unique_ptr
#include <fstream>
#include <cmath>
#include <tuple>
#include <mutex>
#include "openmc/source.h"
#include "openmc/particle.h"
#include "spline.h"

// EDIT THIS TO CHANGE WHERE THE CROSS SECTION DATA ARE STORED
std::string FilePath = "/home/ADF/mjc970/Documents/High Flux Neutron Source/SOURCE_FILES/Compiled Source Method/src/";

class CompiledSource : public openmc::Source
{

public:
  CompiledSource(double Ep) : Ep_{Ep}
  {
    // Load data in constructor so that it can be accessed by all future sample calls across all threads
    loadData();
    // Set maximum depth variable for all sampling at this energy to use
    setMaxDepth();
  }

  // Samples from an instance of this class.
  openmc::SourceSite sample(uint64_t *seed) const
  {
    openmc::SourceSite particle;
    // Set particle type to neutron
    particle.particle = openmc::ParticleType::neutron();

    // Uncomment/comment to optionally smear energy

    // // Define a normal distribution to sample from
    // double E_sigma = 0.0075; // 1.2 keV --> MeV
    // openmc::Normal E_dist(Ep_, E_sigma);
    // // Then sample from it
    // double EE = E_dist.sample(seed);

    double EE = Ep_; // Uncomment if smearing is to be turned off
    // Interpolate strength of source (per mC) for input proton energy
    particle.wgt = linInterp(yield_E, yield_N, EE);
    // Sample particle starting position in x-y plane (10 cm disk)
    openmc::PowerLaw r_dist(0, 5, 1);
    double pos_r = r_dist.sample(seed).first;
    double pos_ang = 2.0 * M_PI * openmc::prn(seed);
    // Sample the position of iteraction in the x-y plane
    double x_pos = pos_r * std::cos(pos_ang); // r_dist.sample(seed);
    double y_pos = pos_r * std::sin(pos_ang); // r_dist.sample(seed);
    // Set placeholder z positions and angle. These will be changed later by reference
    double z_pos = 0.0;
    particle.u = {0.0, 0.0, 1.0};
    // Required variables
    double Eint = 0.0;
    double thetaCM = 1.0;
    double thetaLab = 1.0;
    double phiLab = 1.0;
    double energy = 1.0;
    bool excited = false;
    // Use a while loop to avoid returning erroneous samples
    while (true)
    {
      // First, find the energy of the proton at the point of interaction in the lithium (need to return the interaction depth for OpenMC)
      Eint = InteractionE(EE, seed, &z_pos, &excited);
      // Then, sample the neutron emission angle in the CM frame
      thetaCM = InteractionTheta(Eint, seed, &excited);
      // Since thetaCM returns 0 when it is invalid, skip and sample again
      if (thetaCM == 0)
      {
        continue;
      }
      thetaLab = ConvertThetaCMtoLab(Eint, thetaCM, &excited);
      if (thetaLab == 0)
      {
        continue;
      }
      // And calculate the final neutron energy in the lab frame based on this angle
      energy = NeutronEnergy(Eint, thetaCM, seed, &phiLab, &excited);
      // Since NeutronEnergy returns 0 when it is invalid, skip and sample again
      if (energy == 0)
      {
        continue;
      }
      break;
    }
    // Set neutron energy
    particle.E = energy * 1E6;
    // Set neutron direction
    particle.u = {std::sin(thetaLab) * std::cos(phiLab), std::sin(thetaLab) * std::sin(phiLab), std::cos(thetaLab)};
    // Set mutex to prevent cout from being accessed by multiple processes simultaneously
    // static std::mutex lock;
    // std::lock_guard guard{ lock };
    // std::cout << thetaLab*180/M_PI << std::endl; // Output cos(theta) for logging

    // Set neutron position
    // y: -42 cm for beam y offset from model origin
    // z: Need to convert z_pos from mm to cm
    // z: And -4.51 offset to account for model position
    particle.r = {x_pos, y_pos, z_pos / 10 - 4.51};
    // particle.r = {0., -42., 0.}; // Alternative to use point source starting location
    return particle;
  }

private:
  // Create required variables
  const double Ep_;
  // Cross section data
  std::vector<double> XS_t_x;
  std::vector<double> XS_t_y;
  std::vector<double> XS_0_x;
  std::vector<double> XS_0_y;
  std::vector<double> XS_1_x;
  std::vector<double> XS_1_y;
  // Angular distribution data
  std::vector<double> XS_ang0_E;
  std::vector<double> XS_ang0_A0;
  std::vector<double> XS_ang0_A1;
  std::vector<double> XS_ang0_A2;
  std::vector<double> XS_ang1_E;
  std::vector<double> XS_ang1_A0;
  std::vector<double> XS_ang1_A1;
  std::vector<double> XS_ang1_A2;
  // Splines
  tk::spline A0;
  tk::spline A1;
  tk::spline A2;
  // Yield data
  std::vector<double> yield_E;
  std::vector<double> yield_N;
  // Define required constants
  const double mp = 1.007276; // 1.007825;
  const double mn = 1.008665;
  const double mLi = 7.014358; // 7.016004;
  const double mBe = 7.013809; // 7.016929;
  const double amu = 931.4941;
  const double Q = -1.64424;
  const double Ethresh = 1.880356; // 1.8803; // MeV
  const int max_samples = 1000;
  double max_depth = 0.2; // mm

  void loadData()
  {
    std::cout << " Reading data for compiled source..." << std::endl;
    // Load total XS_t data
    // std::ifstream file0("/home/ADF/mjc970/Documents/High Flux Neutron Source/Python Module/hfadnef/hfadnef/src/XS_t_endf.txt");
    std::ifstream file0(FilePath + "XS_t_ideal.txt"); // XS_t_ideal.txt" l&p_xs
    if (!file0.is_open())
    {
      // Handle file opening error
      throw std::runtime_error("Failed to open file: XS_t_ideal.txt");
    }

    double value1, value2, value3, value4;
    while (file0 >> value1 >> value2)
    {
      // std::cout << "XS_t" << Ep_ << std::endl;
      XS_t_x.push_back(value1);
      XS_t_y.push_back(value2);
    }

    file0.close();

    // Load total XS_0 data
    std::ifstream file1(FilePath + "XS_0.txt");
    if (!file1.is_open())
    {
      // Handle file opening error
      throw std::runtime_error("Failed to open file: XS_0.txt");
    }
    while (file1 >> value1 >> value2)
    {
      // std::cout << "XS_0" << Ep_ << std::endl;
      XS_0_x.push_back(value1);
      XS_0_y.push_back(value2);
    }
    file1.close();

    // Load total XS_1 data
    std::ifstream file2(FilePath + "XS_1.txt");
    if (!file2.is_open())
    {
      // Handle file opening error
      throw std::runtime_error("Failed to open file: XS_1.txt");
    }
    while (file2 >> value1 >> value2)
    {
      // std::cout << "XS_1" << Ep_ << std::endl;
      XS_1_x.push_back(value1);
      XS_1_y.push_back(value2);
    }
    file2.close();

    // Load angular distribution data for ground state
    std::ifstream file3(FilePath + "XS_ang0.txt");
    if (!file3.is_open())
    {
      // Handle file opening error
      throw std::runtime_error("Failed to open file: XS_ang0.txt");
    }
    while (file3 >> value1 >> value2 >> value3 >> value4)
    {
      // std::cout << "1: " << value1 << "2: " << value2 << "3: " << value3 << "4: " << value4 << std::endl;
      XS_ang0_E.push_back(value1);
      XS_ang0_A0.push_back(value2);
      XS_ang0_A1.push_back(value3);
      XS_ang0_A2.push_back(value4);
    }

    file3.close();
    // Create splines for sampling
    A0.set_points(XS_ang0_E, XS_ang0_A0);
    A1.set_points(XS_ang0_E, XS_ang0_A1);
    A2.set_points(XS_ang0_E, XS_ang0_A2);

    // Load angular distribution data for excited state
    std::ifstream file4(FilePath + "XS_ang1.txt");
    if (!file4.is_open())
    {
      // Handle file opening error
      throw std::runtime_error("Failed to open file: XS_ang1.txt");
    }
    while (file4 >> value1 >> value2 >> value3 >> value4)
    {
      XS_ang1_E.push_back(value1);
      XS_ang1_A0.push_back(value2);
      XS_ang1_A1.push_back(value3);
      XS_ang1_A2.push_back(value4);
    }
    file4.close();

    // Load yield data
    std::ifstream file5(FilePath + "yield.txt");
    if (!file5.is_open())
    {
      // Handle file opening error
      throw std::runtime_error("Failed to open file: yield.txt");
    }
    while (file5 >> value1 >> value2)
    {
      yield_E.push_back(value1);
      yield_N.push_back(value2);
    }
    file5.close();
  };

  void setMaxDepth()
  {
    // std::cout << "Setting maximum depth" << std::endl;
    // Coefficients for linear fit to negative inverse stopping power
    const double m = -0.00512339;
    const double c = -0.00328388;
    // Integral of negative inverse stopping power used to calculate maximum range
    max_depth = 10 * ((m / 2) * (Ethresh * Ethresh - Ep_ * Ep_) + c * (Ethresh - Ep_)) + 0.001; // Add 1 micrometer extra range ensure threshold is covered
  }

  double linInterp(const std::vector<double> &x_values, const std::vector<double> &y_values, double x) const
  {
    // Check whether vectors are the same size (can probably remove for performance)
    if (x_values.size() != y_values.size() || x_values.size() < 2)
    {
      throw std::invalid_argument("Invalid input vectors");
    }

    // Check if x is within the range of x values

    if (x < x_values.front() || x > x_values.back())
    {
      // std::cout << x << " " << x_values.front()<< " " << x_values.back() << std::endl;
      throw std::out_of_range("Input x value is outside the range of x values");
    }
    // Find the interval where x lies
    size_t i = 0;
    while (x > x_values[i + 1])
    {
      ++i;
    }
    // Linear interpolation formula
    double m = (y_values[i + 1] - y_values[i]) / (x_values[i + 1] - x_values[i]);
    double c = y_values[i] - m * x_values[i];
    double y = m * x + c;
    return y;
  }

  double InteractionE(double Ep, uint64_t *seed, double *z_pos, bool *excited) const
  {

    // Integer number of attempts to rejection sample
    int i = 0;
    // Maximum cross section from XS_t for rejection sampling
    int max_XS = 1000;
    // Rejection sample until found or otherwise
    while (true)
    {
      // Choose a random depth in the lithium for the reaction to occur
      double trial_z = max_depth * openmc::prn(seed);
      // Set neutron z position to trial_z
      *z_pos = trial_z;
      // Calculate the energy of the proton at this depth in the lithium
      double energy = Energyloss(Ep, trial_z);
      // Check whether an interaction can occur at this proton energy
      if (energy <= Ethresh)
      {
        // If not, then use 'continue' to skip to the next iteration of the while loop
        continue;
      }
      double XS = 0.;
      // If the proton energy is low, then use theoretical calculations instead of tabulated data
      if (energy <= 1.925)
      {
        // Constants and equation from Lee and Zhou
        double C0 = 6;      // Consistent with Gibbons and Macklin XS data
        double A = 164.913; // mb MeV/sr
        double x = C0 * std::sqrt(1 - Ethresh / energy);
        XS = 4. * M_PI * (A * x) / (energy * (1 + x) * (1 + x));
      }
      // Otherwise, sample the cross section from tablated data
      else
      {
        XS = linInterp(XS_t_x, XS_t_y, energy);
      }
      // Now check whether we could have populated an excited state in berilium-8
      // This ratio describes the probability of being in an excited state
      // Sample from this probability (ie. random number above or below)
      if (energy > 2.37)
      {
        double ratio = linInterp(XS_1_x, XS_1_y, energy) / linInterp(XS_0_x, XS_0_y, energy);
        if (openmc::prn(seed) < ratio)
        {
          *excited = true; // Set excited variable by pointer reference
        }
      }
      else
      {
        *excited = false;
      }
      // Use calculated/tabulated XS to rejection sample
      // If we know what the cross section is at this energy, we can generate a random XS between 0 and 1000
      // If this is greater than the calculated XS, accept, if lower then reject
      if (XS > max_XS)
      {
        std::cout << "Max XS is exceeded!!" << std::endl;
      }
      if (XS > openmc::prn(seed) * max_XS)
      {

        return energy;
      }
      i++;
      if (i > max_samples)
      {
        std::cout << "Rejection sampling not working, taken more than 1000 samples" << std::endl;
        return 0;
      }
    }
    return 0;
  };

  double Energyloss(double Ep, double z) const
  {
    // Compute maximum range of proton with energy Ep in lithium (functional fit from SRIM)
    double range = 0.0237 * Ep * Ep + 0.0369 * Ep - 0.0102; // mm
    if (z > range)
    {
      return 0; // Fully stopped before Z
    }
    // Find range of proton at the point of interaction
    range -= z;
    // Calculate the energy of the proton at this point, based on its remaining range (functional fit from SRIM)
    double Eout = -7.98 * range * range + 10.0761 * range + 0.6029;
    return Eout;
  }

  double InteractionTheta(double Ep, uint64_t *seed, bool *excited) const
  {

    // We need to rejection sample from the differential cross section data. This is energy dependent, and requires
    // parameters A_1, A_2, A_3.
    // We first define an array of these values
    double A[3] = {1, 0, 0};
    // Define a maximum cross section for rejection sampling
    double XSmax = 3;
    // Integer number of attempts to rejection sample
    int i = 0;
    // Below where we have measured data we assume (fairly) that neutron emission is isotropic, so A = [1,0,0]
    // Above this point, we must sample from experimental data
    if (Ep >= 1.95)
    {
      // If an excited state of Be is produced, we take the A_i parameters from the excited distributions
      // Since we only have two data points for these, we linearly interpolate
      if (*excited == true)
      {
        // static std::mutex lock;
        // std::lock_guard guard{ lock };
        // std::cout << "Sampled angle" << std::endl; // Output cos(theta) for logging
        if (Ep >= 2.5)
        {
          // static std::mutex lock;
          // std::lock_guard guard{ lock };
          // std::cout << "linterp" << std::endl; // Output cos(theta) for logging
          A[0] = linInterp(XS_ang1_E, XS_ang1_A0, Ep);
          A[1] = linInterp(XS_ang1_E, XS_ang1_A1, Ep);
          A[2] = linInterp(XS_ang1_E, XS_ang1_A2, Ep);
        }
        // Where we have no data, assume isotropic
        else
        {
          // static std::mutex lock;
          // std::lock_guard guard{ lock };
          // std::cout << "isotropic" << std::endl; // Output cos(theta) for logging
          A[0] = 1.;
          A[1] = 0.;
          A[2] = 0.;
        }
      }
      // Otherwise, we must have produced the Be in the ground state
      // Here, we sample from the cubic splines that we have previously found
      else
      {
        A[0] = A0(Ep);
        A[1] = A1(Ep);
        A[2] = A2(Ep);
        // Alternative code to linearly interpolate instead
        // A[0] = linInterp(XS_ang0_E, XS_ang0_A0, Ep);
        // A[1] = linInterp(XS_ang0_E, XS_ang0_A1, Ep);
        // A[2] = linInterp(XS_ang0_E, XS_ang0_A2, Ep);
      }
    }
    else
    {
      A[0] = 1.;
      A[1] = 0.;
      A[2] = 0.;
    }
    // Now, perform rejection sampling]
    while (true)
    {
      // Sample a random emission angle
      double random_theta = 2 * M_PI * openmc::prn(seed);
      // Calculate XS based on Legendre polynomial expansion given by Liskien
      double XS = std::sin(random_theta) * (A[0] * 1 + A[1] * std::cos(random_theta) + A[2] * 0.5 * (3 * std::cos(random_theta) * std::cos(random_theta) - 1));
      // Sample an XS for rejection sampling
      double XS_sample = XSmax * openmc::prn(seed);
      // The following case means an error has occured
      if (XS > XSmax)
      {
        std::cout << "Sampling error for angle" << std::endl;
      }
      if (XS_sample < XS)
      {
        return random_theta;
      }
      i++;
      // If sampling takes too long, reject it
      if (i > max_samples)
      {
        return 0;
      }
    }
  }

  double ConvertThetaCMtoLab(double Ep, double thetaCM, bool *excited) const
  {

    // First, check if Be is excited. If so, this energy needs to be accounted for in the kinematics
    double Ex = 0.0;
    if (*excited == true)
    {
      Ex = 0.429;
    }
    // Now, complete kinematics
    // We can compute the ratio of the velocities of the CM and of the neutron (vCM/vn)
    // double gamma = std::sqrt(((mp*mn)/(mBe*mLi))*(Ep/(Ep+(Q-Ex)*(1.+mp/mLi)))); //  Calculate gamma for the conversion to lab COM
    double gamma = std::sqrt((mp * mn * Ep * (mn + mBe)) / (mBe * (mp + mLi) * (mp + mLi) * ((Q - Ex) + Ep * (mLi / (mp + mLi))))); // Correct gamma without assumption (mLi + mp = mBe + mn)
    // And use this in the angle conversion formula
    double theta = std::atan2((std::sin(thetaCM)), (std::cos(thetaCM) + gamma));
    // std::cout << "Calculated theta lab: " << theta << std::endl; // Output cos(theta) for logging
    return theta;
  }

  double NeutronEnergy(double Ep, double thetaCM, uint64_t *seed, double *phi, bool *excited) const
  {

    // First, check if Be is excited. If so, this energy needs to be accounted for in the kinematics
    double Ex = 0.0;
    if (*excited == true)
    {
      // std::cout << "Exciting..." << std::endl;
      Ex = 0.429;
    }
    // Compute centre of mass energy
    double ECM = mLi * Ep / (mLi + mp) + Q - Ex;
    // std::cout << "ECM: " << ECM << std::endl;
    if (ECM < 0)
    {
      // std::cout << "Below threshold! Something is wrong!" << std::endl;
      return 0;
    }
    // Calculate neutron total momentum
    double p_nT = std::sqrt(2.0 * ECM * mn * (mBe / (mn + mBe)));
    // Sample a random phi vector
    double sample_phi = 2 * M_PI * openmc::prn(seed);
    // Use both theta and phi to compute vector components of neutron momentum
    double pn[3] = {p_nT * std::sin(thetaCM) * std::cos(sample_phi), p_nT * std::sin(thetaCM) * std::sin(sample_phi), p_nT * std::cos(thetaCM)};
    // Boost the z component into the lab frame
    pn[2] += (mp / (mp + mLi)) * std::sqrt(2.0 * mp * Ep);
    // Compute lab energy
    double En = 0;
    for (int i = 0; i < 3; i++)
    {
      En += pn[i] * pn[i];
    }
    En /= (2.0 * mn);
    *phi = sample_phi;
    return En;
  }
};

// Create a smart pointer to the CompiledSource object
extern "C" std::unique_ptr<CompiledSource> openmc_create_source(std::string parameter)
{
  double Ep = std::stod(parameter); // Convert input parameter from string to double
  return std::make_unique<CompiledSource>(Ep);
}

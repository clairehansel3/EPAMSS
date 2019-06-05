// Copyright (C) 2019 Claire Hansel
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#ifndef EPAMSS_INPUT_HXX
#define EPAMSS_INPUT_HXX

#include <boost/mpi/communicator.hpp>
#include <cstddef>
#include <string>

class Parameters {

public:

  Parameters(const char* input_path, boost::mpi::communicator& world);

  // specified in input
  double rho_ion_initial_si;
  double plasma_length_si;
  double beam_energy_initial_gev;
  double acceleration_gradient_gev_per_m;
  double bennett_radius_initial_si;
  double cross_section_radius_si;
  double unperturbed_plasma_density_si;
  double integration_tolerance;
  double vartheta_cutoff;
  std::size_t ion_atomic_number;
  std::size_t minimum_steps_per_betatron_period;
  std::size_t particles_target;
  std::size_t analysis_points_target;
  std::size_t spline_points;
  unsigned seed;
  unsigned max_order;
  unsigned max_integration_depth;
  std::string output_filename;
  std::string statistics_filename;
  std::string phase_space_filename;
  bool output_phase_space;
  bool modified_bennett;

  // computed from other parameters
  double plasma_frequency_si;
  double plasma_angular_wavenumber_si;
  double plasma_skin_depth_si;
  double rho_ion_initial;
  double plasma_length;
  double gamma_initial;
  double gamma_prime;
  double bennett_radius_initial;
  double cross_section_radius;
  double delta;
  double gamma_final;
  double bennett_radius_final;
  double rho_ion_final;
  double betatron_frequency_final;
  double betatron_period_final;
  double betatron_frequency_final_si;
  double betatron_period_final_si;
  double step_size;
  double step_size_si;
  double gamma_minimum_angle;
  double omega_off_axis;
  double omega_on_axis_initial;
  double max_scattering_r_div_a_initial;
  double sigma_r_initial; // NOTE: drop initial ALSO TURN ON SCATTERNG
  double sigma_r_prime_initial;
  std::size_t steps;
  std::size_t stride;
  std::size_t compute_processes;
  std::size_t actual_particles;
  std::size_t particles_per_process;
  std::size_t actual_analysis_points;

  void writeOutputFile(double seconds);

private:

  void computeDependentParameters(int processes);

};

#endif

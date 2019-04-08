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
  double maximum_ion_density;               // [m^-3]
  double maximum_electron_density;
  double plasma_length;                     // [m]
  double beam_energy;                       // [GeV]
  double bennett_radius;                    // [m]
  double interaction_radius;                // [m]
  double integration_tolerance;
  double vartheta_cutoff;
  double unperturbed_plasma_density;
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

  // computed from other parameters
  double ion_linear_density;
  double electron_linear_density;
  double gamma;
  double alpha;
  double lambda;
  double betatron_frequency;
  double betatron_period;
  double step_size;
  double cross_section;
  double minimum_angle;
  double max_scattering_r_div_a;
  double percent_with_scattering;
  double omega_on_axis;
  double sigma_dist;
  double sigma;
  std::size_t steps;
  std::size_t stride;
  std::size_t actual_particles;
  std::size_t particles_per_process;
  std::size_t compute_processes;
  std::size_t actual_analysis_points;

  void writeOutputFile(double seconds);

private:

  void computeDependentParameters(int processes);

};

#endif

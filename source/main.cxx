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

#include "input.hxx"
#include "random.hxx"
#include "scattering.hxx"
#include "solver.hxx"
#include "statistics.hxx"
#include <algorithm>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <ctime>
#include <memory>
#include <stdexcept>

static void runMainProcess(boost::mpi::communicator& world, Parameters& p)
{
  // allocate statistics
  auto statistics = std::make_unique<Statistics[]>(2 * p.actual_analysis_points);
  {
    // no scattering
    auto temp_statistics = std::make_unique<Statistics[]>(p.actual_analysis_points);
    for (std::size_t i = 0; i != p.compute_processes; ++i) {
      world.recv(i + 1, 1, reinterpret_cast<char*>(temp_statistics.get()),
        sizeof(Statistics) * p.actual_analysis_points);
      std::cout << "received (ns) statistics for process " << i + 1 << std::endl;
      for (std::size_t j = 0; j != p.actual_analysis_points; ++j)
        statistics[j] = Statistics(statistics[j], temp_statistics[j]);
    }
    // scattering
    for (std::size_t i = 0; i != p.compute_processes; ++i) {
      world.recv(i + 1, 1, reinterpret_cast<char*>(temp_statistics.get()),
        sizeof(Statistics) * p.actual_analysis_points);
    std::cout << "received (s) statistics for process " << i + 1 << std::endl;
      for (std::size_t j = 0; j != p.actual_analysis_points; ++j)
        statistics[p.actual_analysis_points + j] = Statistics(
          statistics[p.actual_analysis_points + j], temp_statistics[j]);
    }
  }
  // write to file
  std::cout << "writing to file" << std::endl;
  std::ofstream statistics_file;
  statistics_file.exceptions(std::ofstream::failbit | std::ofstream::badbit);
  statistics_file.open(p.statistics_filename, std::ofstream::binary);
  for (std::size_t i = 0; i != 2 * p.actual_analysis_points; ++i)
    statistics[i].writeToFile(statistics_file);
}

static void runComputeProcess(boost::mpi::communicator& world, Parameters& p)
{
  // initialize scattering
  Scattering scattering{p.max_order, p.max_integration_depth,
    p.integration_tolerance, p.spline_points, p.vartheta_cutoff,
    p.omega_off_axis, p.gamma_minimum_angle, world.rank() == 1};

  // initialize beam
  auto beam = std::make_unique<Particle[]>(p.particles_per_process);
  initializeBeam(beam.get(), p.particles_per_process, p.bennett_radius_initial,
    p.gamma_initial, p.sigma_r_initial, p.sigma_r_prime_initial,
    p.modified_bennett);

  // open phase space file
  std::ofstream phase_space_file;
  std::ofstream* phase_space_file_ptr = nullptr;
  if (p.output_phase_space) {
    phase_space_file.exceptions(std::ofstream::failbit | std::ofstream::badbit);
    phase_space_file.open(p.phase_space_filename + "_" +
      std::to_string(world.rank()), std::ofstream::binary);
    phase_space_file_ptr = &phase_space_file;
  }

  // no scattering
  {
    auto beam_copy = std::make_unique<Particle[]>(p.particles_per_process);
    std::copy_n(beam.get(), p.particles_per_process, beam_copy.get());
    auto statistics = std::make_unique<Statistics[]>(p.actual_analysis_points);
    solve(beam_copy.get(), statistics.get(), scattering, phase_space_file_ptr,
      p.particles_per_process, p.steps, p.stride, p.ion_atomic_number,
      p.step_size, p.bennett_radius_initial, p.rho_ion_initial, p.gamma_initial,
      p.gamma_prime, p.delta, false, world.rank() == 1);
    world.send(0, 1, reinterpret_cast<char*>(statistics.get()),
       sizeof(Statistics) * p.actual_analysis_points);
  }

  // scattering
  {
    auto statistics = std::make_unique<Statistics[]>(p.actual_analysis_points);
    solve(beam.get(), statistics.get(), scattering, phase_space_file_ptr,
      p.particles_per_process, p.steps, p.stride, p.ion_atomic_number,
      p.step_size, p.bennett_radius_initial, p.rho_ion_initial, p.gamma_initial,
      p.gamma_prime, p.delta, true, world.rank() == 1);
    world.send(0, 1, reinterpret_cast<char*>(statistics.get()),
       sizeof(Statistics) * p.actual_analysis_points);
  }
}

int main(int argc, char* argv[])
{
  // Start timer
  std::time_t begin;
  std::time(&begin);

  // Initialize MPI
  boost::mpi::environment env;
  boost::mpi::communicator world;

  // read input file
  if (argc != 2)
    throw std::runtime_error("expected a single argument containing the path to"
      "the input file");
  Parameters parameters{argv[1], world};

  // seed random number generator
  std::cout << "process: " << world.rank() << " seed: " << parameters.seed <<
    std::endl;
  seedRandom(world.rank(), parameters.seed);

  // run program
  if (world.rank() == 0) {
    runMainProcess(world, parameters);
    std::time_t end;
    std::time(&end);
    double diff = std::difftime(end, begin);
    parameters.writeOutputFile(diff);
    return EXIT_SUCCESS;
  }
  else {
    runComputeProcess(world, parameters);
    return EXIT_SUCCESS;
  }
}

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

#include "solver.hxx"
#include "random.hxx"
#include "scattering.hxx"
#include "statistics.hxx"
#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <iostream>

void initializeBeam(Particle* beam, std::size_t particles, double sigma,
  double bennett_radius, double sigma_dist)
{
  for (std::size_t particle = 0; particle != particles; ++particle) {
    double r;
    while (true) {
      auto test_r = bennett_radius / std::sqrt((1 / randomUniform()) - 1);
      if (randomUniform() < std::exp(-test_r * test_r / (2 * sigma * sigma))) {
        r = test_r;
        break;
      }
    }
    double theta = 2 * boost::math::constants::pi<double>() * randomUniform();
    double r_prime = randomNormal(0, sigma_dist);
    double r_theta_prime = randomNormal(0, sigma_dist);
    beam[particle].x = r * std::cos(theta);
    beam[particle].y = r * std::sin(theta);
    beam[particle].vx = r_prime * std::cos(theta) - r_theta_prime * std::sin(
        theta);
    beam[particle].vy = r_prime * std::sin(theta) + r_theta_prime * std::cos(
        theta);
  }
}

void solve(Particle* beam, Statistics* statistics, Scattering& scattering,
  std::ofstream* phase_space_file, std::size_t particles, std::size_t steps,
  std::size_t stride, double bennett_radius, double step_size, double alpha,
  double lambda, double maximum_ion_density, double cross_section,
  double minimum_angle, bool enable_scattering, bool print_progress)
{
  double a2 = bennett_radius * bennett_radius;

  int percent = -1;

  for (std::size_t step = 0; step != steps + 1; ++step) {

    // Print progress
    if (print_progress) {
      int new_percent = (100 * step) / steps;
      if (percent != new_percent) {
        percent = new_percent;
        std::cout << percent << '%' << std::endl;
      }
    }

    // If the following condition is true, this step is an analysis point and
    // information about the beam is recorded.
    if (step % stride == 0) {
      if (phase_space_file) {
        // Write particles to phase space file.
        phase_space_file->write(reinterpret_cast<char*>(beam),
          sizeof(Particle) * particles);
      }
      // Compute statistics
      statistics[step / stride] = Statistics(beam, particles);
    }

    // track each particle
    for (std::size_t particle = 0; particle != particles; ++particle) {

      // get old positions and velocities
      double x_old  = beam[particle].x;
      double vx_old = beam[particle].vx;
      double y_old  = beam[particle].y;
      double vy_old = beam[particle].vy;
      double r2_old = x_old * x_old + y_old * y_old;

      // compute force
      double value = -alpha * (lambda + 1 / (1 + r2_old / a2));
      double fx_old = x_old * value;
      double fy_old = y_old * value;

      // update positions
      double x_new = x_old + step_size * (vx_old + 0.5 * fx_old * step_size);
      double y_new = y_old + step_size * (vy_old + 0.5 * fy_old * step_size);
      double r2_new = x_new * x_new + y_new * y_new;

      // compute new forces
      value = -alpha * (lambda + 1 / (1 + r2_new / a2));
      double fx_new = x_new * value;
      double fy_new = y_new * value;

      // compute new velocities
      double vx_new = vx_old + 0.5 * (fx_old + fx_new) * step_size;
      double vy_new = vy_old + 0.5 * (fy_old + fy_new) * step_size;

      // compute scattering
      if (enable_scattering)
        scattering.scatter(x_new, y_new, vx_new, vy_new,
          bennett_radius, maximum_ion_density, cross_section, step_size,
          minimum_angle, lambda);

      // write result to beam
      beam[particle].x = x_new;
      beam[particle].vx = vx_new;
      beam[particle].y = y_new;
      beam[particle].vy = vy_new;
    }
  }
}

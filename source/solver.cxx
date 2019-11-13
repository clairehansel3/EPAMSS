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

void initializeBeam(Particle* beam, std::size_t particles,
  double bennett_radius_initial, double gamma_initial, double gamma_prime,
  double sigma_r, double sigma_r_prime_initial, bool modified_bennett)
{
  double sqrt_gamma = std::sqrt(gamma_initial);

  for (std::size_t particle = 0; particle != particles; ++particle) {

    #ifdef EPAMSS_TEST_SCATTERING

    beam[particle].x = 0;
    beam[particle].y = 0;
    beam[particle].vx = 0;
    beam[particle].vy = 0;

    #else

    // sample r using rejection sampling
    double r;
    while (true) {
      auto test_r = bennett_radius_initial / std::sqrt((1 / randomUniform()) - 1);
      if (!modified_bennett || randomUniform() < std::exp(-test_r * test_r / (2 * sigma_r * sigma_r))) {
        r = test_r;
        break;
      }
    }

    // sample remaining coordinates
    double theta = 2 * boost::math::constants::pi<double>() * randomUniform();
    double r_prime = randomNormal(0, sigma_r_prime_initial);
    double r_theta_prime = randomNormal(0, sigma_r_prime_initial);

    // convert to phase space
    double x = r * std::cos(theta);
    double y = r * std::sin(theta);
    double vx = (r_prime * std::cos(theta) - r_theta_prime * std::sin(theta));
    double vy = (r_prime * std::sin(theta) + r_theta_prime * std::cos(theta));
    beam[particle].x = sqrt_gamma * x;
    beam[particle].vx = vx * sqrt_gamma + 0.5 * x * gamma_prime / sqrt_gamma;
    beam[particle].y = sqrt_gamma * y;
    beam[particle].vy = vy * sqrt_gamma + 0.5 * y * gamma_prime / sqrt_gamma;

    #endif

  }
}


void solve(Particle* beam, Statistics* statistics, Scattering& scattering,
  std::ofstream* phase_space_file, std::size_t particles, std::size_t steps,
  std::size_t stride, std::size_t ion_atomic_number, double step_size,
  double bennett_radius_initial, double rho_ion_div_n0, double gamma_initial,
  double gamma_prime, double delta, double drive_amplitude,
  double drive_angular_frequency, bool scattering_enabled, bool print_progress)
{

  int percent = -1;

  for (std::size_t step = 0; step != steps + 1; ++step) {

    // print progress
    if (print_progress) {
      int new_percent = (100 * step) / steps;
      if (percent != new_percent) {
        percent = new_percent;
        std::cout << percent << '%' << std::endl;
      }
    }

    // compute z dependent parameters
    double gamma_over_gamma_initial = 1 + (gamma_prime * step_size * step / gamma_initial);
    double gamma = gamma_initial * gamma_over_gamma_initial;
    double bennett_radius = bennett_radius_initial * std::pow(gamma_over_gamma_initial, -0.25);
    double gamma_over_gamma_initial_next = 1 + (gamma_prime * step_size * (step + 1) / gamma_initial);
    double gamma_next = gamma_initial * gamma_over_gamma_initial_next;
    double bennett_radius_next = bennett_radius_initial * std::pow(gamma_over_gamma_initial_next, -0.25);
    double drive_multiplier = 1 + drive_amplitude * std::sin(drive_angular_frequency * step_size * step);
    double drive_multiplier_next = 1 + drive_amplitude * std::sin(drive_angular_frequency * step_size * (step + 1));

    // if the step is an analysis step, record data
    if (step % stride == 0) {
      if (phase_space_file) {
        //write particles to phase space file.
        phase_space_file->write(reinterpret_cast<char*>(beam),
          sizeof(Particle) * particles);
      }
      // compute statistics
      statistics[step / stride] = Statistics(beam, particles, gamma, gamma_prime);
    }

    for (std::size_t particle = 0; particle != particles; ++particle) {

      // get old phase space coordinates
      double x_old  = beam[particle].x;
      double vx_old = beam[particle].vx;
      double y_old  = beam[particle].y;
      double vy_old = beam[particle].vy;

      #ifdef EPAMSS_TEST_SCATTERING

      double x_new = x_old + step_size * vx_old;
      double y_new = y_old + step_size * vy_old;
      double vx_new = vx_old;
      double vy_new = vy_old;

      #else

      // compute force
      double r2_old = x_old * x_old + y_old * y_old;
      double value_old = (-0.5 * ion_atomic_number / gamma) * drive_multiplier * (delta + rho_ion_div_n0 / (1 + (r2_old / (gamma * bennett_radius * bennett_radius))));
      value_old += 0.25 * gamma_prime * gamma_prime / (gamma * gamma);
      double fx_old = value_old * x_old;
      double fy_old = value_old * y_old;


      // compute new positions
      double x_new = x_old + step_size * (vx_old + 0.5 * step_size * fx_old);
      double y_new = y_old + step_size * (vy_old + 0.5 * step_size * fy_old);

      // compute new force
      double r2_new = x_new * x_new + y_new * y_new;
      double value_new = (-0.5 * ion_atomic_number / gamma_next) * drive_multiplier_next * (delta + rho_ion_div_n0 / (1 + (r2_new / (gamma_next * bennett_radius_next * bennett_radius_next))));
      value_new += 0.25 * gamma_prime * gamma_prime / (gamma_next * gamma_next);
      double fx_new = value_new * x_new;
      double fy_new = value_new * y_new;

      // compute new momenta
      double vx_new = vx_old + 0.5 * step_size * (fx_old + fx_new);
      double vy_new = vy_old + 0.5 * step_size * (fy_old + fy_new);

      #endif

      // compute scattering
      if (scattering_enabled)
        scattering.scatter(x_new, vx_new, y_new, vy_new, gamma_next, gamma_prime, bennett_radius_next, rho_ion_div_n0, delta, drive_multiplier_next);

      // write result to beam
      beam[particle].x = x_new;
      beam[particle].vx = vx_new;
      beam[particle].y = y_new;
      beam[particle].vy = vy_new;
    }
  }
}

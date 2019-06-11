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
  double bennett_radius_initial, double gamma_initial, double sigma_r,
  double sigma_r_prime_initial, bool modified_bennett)
{
  for (std::size_t particle = 0; particle != particles; ++particle) {

    #ifdef EPAMSS_TEST_SCATTERING

    beam[particle].x = 0;
    beam[particle].y = 0;
    beam[particle].px = 0;
    beam[particle].py = 0;

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
    beam[particle].x = r * std::cos(theta);
    beam[particle].y = r * std::sin(theta);
    beam[particle].px = gamma_initial * (r_prime * std::cos(theta) - r_theta_prime * std::sin(theta));
    beam[particle].py = gamma_initial * (r_prime * std::sin(theta) + r_theta_prime * std::cos(theta));

    #endif

  }
}


void solve(Particle* beam, Statistics* statistics, Scattering& scattering,
  std::ofstream* phase_space_file, std::size_t particles, std::size_t steps,
  std::size_t stride, std::size_t ion_atomic_number, double step_size,
  double bennett_radius_initial, double rho_ion_initial, double gamma_initial,
  double gamma_prime, double delta, bool scattering_enabled,
  bool print_progress)
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

    // if the step is an analysis step, record data
    if (step % stride == 0) {
      if (phase_space_file) {
        //write particles to phase space file.
        phase_space_file->write(reinterpret_cast<char*>(beam),
          sizeof(Particle) * particles);
      }
      // compute statistics
      statistics[step / stride] = Statistics(beam, particles);
    }

    // compute z dependent parameters
    double gamma_over_gamma_initial = 1 + (gamma_prime * step_size * step / gamma_initial);
    double gamma = gamma_initial * gamma_over_gamma_initial;
    double bennett_radius = bennett_radius_initial * std::pow(gamma_over_gamma_initial, -0.25);
    double rho_ion = rho_ion_initial * std::sqrt(gamma_over_gamma_initial);
    double gamma_over_gamma_initial_next = 1 + (gamma_prime * step_size * (step + 1) / gamma_initial);
    double gamma_next = gamma_initial * gamma_over_gamma_initial_next;
    double bennett_radius_next = bennett_radius_initial * std::pow(gamma_over_gamma_initial_next, -0.25);
    double rho_ion_next = rho_ion_initial * std::sqrt(gamma_over_gamma_initial_next);

    for (std::size_t particle = 0; particle != particles; ++particle) {

      // get old phase space coordinates
      double x_old  = beam[particle].x;
      double px_old = beam[particle].px;
      double y_old  = beam[particle].y;
      double py_old = beam[particle].py;

      #ifdef EPAMSS_TEST_SCATTERING

      double x_new = x_old + step_size * (px_old + 0.5) / gamma;
      double y_new = y_old + step_size * (py_old + 0.5) / gamma;
      double px_new = px_old;
      double py_new = py_old;

      #else

      // compute force
      double r2_old = x_old * x_old + y_old * y_old;
      double value_old = -0.5 * ion_atomic_number * (delta + rho_ion / (1 + (r2_old / (bennett_radius * bennett_radius))));
      double fx_old = value_old * x_old;
      double fy_old = value_old * y_old;


      // compute new positions
      double x_new = x_old + step_size * (px_old + 0.5 * step_size * fx_old) / gamma;
      double y_new = y_old + step_size * (py_old + 0.5 * step_size * fy_old) / gamma;

      // compute new force
      double r2_new = x_new * x_new + y_new * y_new;
      double value_new = -0.5 * ion_atomic_number * (delta + rho_ion_next / (1 + (r2_new / (bennett_radius_next * bennett_radius_next))));
      double fx_new = value_new * x_new;
      double fy_new = value_new * y_new;

      // compute new momenta
      double px_new = px_old + 0.5 * step_size * (fx_old + fx_new);
      double py_new = py_old + 0.5 * step_size * (fy_old + fy_new);

      #endif

      // compute scattering
      if (scattering_enabled)
        scattering.scatter(x_new, px_new, y_new, py_new, gamma_next, bennett_radius_next, rho_ion_next, delta);

      // write result to beam
      beam[particle].x = x_new;
      beam[particle].px = px_new;
      beam[particle].y = y_new;
      beam[particle].py = py_new;
    }
  }
}

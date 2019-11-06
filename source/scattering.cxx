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

#include "scattering.hxx"
#include "LambertW/LambertW.h"
#include "random.hxx"
#include "solver.hxx"
#include <boost/math/constants/constants.hpp>
#include <cmath>

Scattering::Scattering(unsigned max_order, unsigned max_integration_depth,
           double integration_tolerance, unsigned long spline_points,
           double vartheta_cutoff, double omega_off_axis,
           double gamma_minimum_angle, bool print_info)
: m_moliere{max_order, max_integration_depth, integration_tolerance,
            spline_points, vartheta_cutoff, print_info},
  m_omega_off_axis{omega_off_axis},
  m_gamma_minimum_angle{gamma_minimum_angle}
{}

void Scattering::scatter(double x, double& vx, double y, double& vy,
  double gamma, double gamma_prime, double bennett_radius, double rho_ion_div_n0,
  double delta, double drive_multiplier)
{
  #ifdef EPAMSS_TEST_SCATTERING
  double omega = m_omega_off_axis;
  #else
  double r2 = x * x + y * y;
  double density = drive_multiplier * (delta + rho_ion_div_n0 * std::pow(1 + r2 / (gamma * bennett_radius * bennett_radius), -2));
  double omega = m_omega_off_axis * density;
  #endif
  double minimum_angle = m_gamma_minimum_angle / gamma;
  if (omega > 25) {
    double constant = -std::exp(-2 * (1 -
      boost::math::constants::euler<double>()));
    double b = -utl::LambertW<-1>(constant / omega);
    double vartheta = m_moliere.sample(b);
    double theta = minimum_angle * std::sqrt(omega * b) * vartheta;
    double beta = 2 * boost::math::constants::pi<double>() * randomUniform();
    double tan_phi_x = std::tan(theta) * std::cos(beta);
    double tan_phi_y = std::tan(theta) * std::sin(beta);
    double vx_actual = (vx - 0.5 * x * gamma_prime / gamma) / std::sqrt(gamma);
    double vy_actual = (vy - 0.5 * y * gamma_prime / gamma) / std::sqrt(gamma);
    double delta_vx_actual = tan_phi_x * (1 + vx_actual * vx_actual) / (1 - vx_actual * tan_phi_x);
    double delta_vy_actual = tan_phi_y * (1 + vy_actual * vy_actual) / (1 - vy_actual * tan_phi_y);
    vx += std::sqrt(gamma) * delta_vx_actual;
    vy += std::sqrt(gamma) * delta_vy_actual;
  }
}

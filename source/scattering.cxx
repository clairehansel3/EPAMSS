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
           double vartheta_cutoff, bool print_info)
: m_moliere{max_order, max_integration_depth, integration_tolerance,
            spline_points, vartheta_cutoff, print_info}
{}

void Scattering::scatter(double x, double y, double& vx, double& vy,
  double bennett_radius, double maximum_ion_density, double cross_section,
  double step_size, double minimum_angle)
{
  double r2 = x * x + y * y;
  double temp = 1 + (r2 / (bennett_radius * bennett_radius));
  double density = maximum_ion_density / (temp * temp);
  double omega = density * cross_section * step_size;
  if (omega > 25) {
    double constant = -std::exp(-2 * (1 -
      boost::math::constants::euler<double>()));
    double B =  -utl::LambertW<-1>(constant / omega);
    double vartheta = m_moliere.sample(B);
    double theta = minimum_angle * std::sqrt(omega * B) * vartheta;
    double beta = 2 * boost::math::constants::pi<double>() * randomUniform();
    double tan_phi_x = std::tan(theta) * std::cos(beta);
    double tan_phi_y = std::tan(theta) * std::sin(beta);
    double delta_vx = tan_phi_x * (1 + vx * vx) / (1 - vx * tan_phi_x);
    double delta_vy = tan_phi_y * (1 + vy * vy) / (1 - vy * tan_phi_y);
    vx += delta_vx;
    vy += delta_vy;
  }
}

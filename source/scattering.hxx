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

#ifndef EPAMSS_SCATTERING_HXX
#define EPAMSS_SCATTERING_HXX

#include "moliere.hxx"

struct Particle;

class Scattering {

public:

  Scattering(unsigned max_order, unsigned max_integration_depth,
             double integration_tolerance, unsigned long spline_points,
             double vartheta_cutoff, double omega_off_axis,
             double gamma_minimum_angle, bool print_info);

  void scatter(double x, double& vx, double y, double& vy, double gamma,
    double gamma_prime, double bennett_radius, double rho_ion_div_n0, double delta,
    double drive_multiplier);

private:

  Moliere m_moliere;

  double m_omega_off_axis, m_gamma_minimum_angle;

};

#endif

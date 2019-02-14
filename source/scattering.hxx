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
             double vartheta_cutoff, bool print_info);

 void scatter(double x, double y, double& vx, double& vy, double bennett_radius,
   double maximum_ion_density, double cross_section, double step_size,
   double minimum_angle);

private:

  Moliere m_moliere;

};

#endif

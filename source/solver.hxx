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

#ifndef EPAMSS_SOLVER_HXX
#define EPAMSS_SOLVER_HXX

#include <cstddef>
#include <fstream>

struct Particle {
  double x, vx, y, vy;
};

class Scattering;
class Statistics;

void initializeBeam(Particle* beam, std::size_t particles,
  double bennett_radius_initial, double gamma_initial, double gamma_prime,
  double sigma_r, double sigma_r_prime_initial, bool modified_bennett);

void solve(Particle* beam, Statistics* statistics, Scattering& scattering,
  std::ofstream* phase_space_file, std::size_t particles, std::size_t steps,
  std::size_t stride, std::size_t ion_atomic_number, double step_size,
  double bennett_radius_initial, double rho_ion_div_n0_initial, double gamma_initial,
  double gamma_prime, double delta, bool scattering_enabled,
  bool print_progress);

#endif

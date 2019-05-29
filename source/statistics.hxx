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

#ifndef EPAMSS_STATISTICS_HXX
#define EPAMSS_STATISTICS_HXX

#include <array>
#include <cstddef>
#include <fstream>

struct Particle;

class Statistics {
// A statistics object contains statistics about a beam. It contains the number
// of particles in the beam, the means of x, px, y, py, and the 4x4 sigma matrix
// corresponding to the variables x, px, y, and py.

public:

  Statistics();
  // Initializes a statistics object corresponding to an empty beam.

  Statistics(Particle* beam, std::size_t particles);
  // Initializes a statistics object corresponding to the specified beam.

  Statistics(Statistics a, Statistics b);
  // From the statistics of two beams, Initializes a statistics object
  // corresponding to the combined beam.

  void writeToFile(std::ofstream& file);
  // Writes statistics to file as binary. The following values are converted to
  // doubles and written in order to the file: mean(x), mean(px), mean(y),
  // mean(py), cov(x, x), cov(x, px), cov(x, y), cov(x, py), cov(px, px),
  // cov(px, y), cov(px, py), cov(y, y), cov(y, py), cov(py, py)

private:

  std::array<double, 4> m_means;
  // The means of x, px, y, py.

  std::array<double, 10> m_covariance_matrix;
  // The entries of the 4x4 covariance matrix for the variables x, px, y, py
  // Since the covariance matrix is symmetric, there are 10 instead of 16
  // entries.

  std::size_t m_particles;
  // The numer of particles in the beam.

};

#endif

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

#include "statistics.hxx"
#include "solver.hxx"

Statistics::Statistics()
: m_means{0.0, 0.0, 0.0, 0.0},
  m_covariance_matrix{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  m_particles{0}
{}

Statistics::Statistics(Particle* beam, std::size_t particles)
: m_means{0.0, 0.0, 0.0, 0.0},
  m_covariance_matrix{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  m_particles{particles}
{
  // compute means
  for (std::size_t particle = 0; particle != particles; ++particle) {
    m_means[0] += beam[particle].x;
    m_means[1] += beam[particle].px;
    m_means[2] += beam[particle].y;
    m_means[3] += beam[particle].py;
  }
  m_means[0] /= particles;
  m_means[1] /= particles;
  m_means[2] /= particles;
  m_means[3] /= particles;

  // compute sigma matrix
  for (std::size_t particle = 0; particle != particles; ++particle) {
    double delta_x  = beam[particle].x  - m_means[0];
    double delta_px = beam[particle].px - m_means[1];
    double delta_y  = beam[particle].y  - m_means[2];
    double delta_py = beam[particle].py - m_means[3];
    m_covariance_matrix[0] += delta_x * delta_x;
    m_covariance_matrix[1] += delta_x * delta_px;
    m_covariance_matrix[2] += delta_x * delta_y;
    m_covariance_matrix[3] += delta_x * delta_py;
    m_covariance_matrix[4] += delta_px * delta_px;
    m_covariance_matrix[5] += delta_px * delta_y;
    m_covariance_matrix[6] += delta_px * delta_py;
    m_covariance_matrix[7] += delta_y * delta_y;
    m_covariance_matrix[8] += delta_y * delta_py;
    m_covariance_matrix[9] += delta_py * delta_py;
  }
  m_covariance_matrix[0] /= particles;
  m_covariance_matrix[1] /= particles;
  m_covariance_matrix[2] /= particles;
  m_covariance_matrix[3] /= particles;
  m_covariance_matrix[4] /= particles;
  m_covariance_matrix[5] /= particles;
  m_covariance_matrix[6] /= particles;
  m_covariance_matrix[7] /= particles;
  m_covariance_matrix[8] /= particles;
  m_covariance_matrix[9] /= particles;
}

Statistics::Statistics(Statistics a, Statistics b)
: m_particles{a.m_particles + b.m_particles}
{
  // combine means
  for (int i = 0; i != 4; ++i) {
    m_means[i] = (a.m_particles * a.m_means[i] + b.m_particles * b.m_means[i])
      / m_particles;
  }

  // combine covariance matrices
  std::array<double, 4> mean_diffs;
  for (int i = 0; i != 4; ++i) {
    mean_diffs[i] = a.m_means[i] - b.m_means[i];
  }
  double value = static_cast<double>(a.m_particles) *
    static_cast<double>(b.m_particles) / static_cast<double>(m_particles);
  m_covariance_matrix[0] = (a.m_particles * a.m_covariance_matrix[0] +
    b.m_particles * b.m_covariance_matrix[0] + mean_diffs[0] * mean_diffs[0] *
    value) / m_particles;
  m_covariance_matrix[1] = (a.m_particles * a.m_covariance_matrix[1] +
    b.m_particles * b.m_covariance_matrix[1] + mean_diffs[0] * mean_diffs[1] *
    value) / m_particles;
  m_covariance_matrix[2] = (a.m_particles * a.m_covariance_matrix[2] +
    b.m_particles * b.m_covariance_matrix[2] + mean_diffs[0] * mean_diffs[2] *
    value) / m_particles;
  m_covariance_matrix[3] = (a.m_particles * a.m_covariance_matrix[3] +
    b.m_particles * b.m_covariance_matrix[3] + mean_diffs[0] * mean_diffs[3] *
    value) / m_particles;
  m_covariance_matrix[4] = (a.m_particles * a.m_covariance_matrix[4] +
    b.m_particles * b.m_covariance_matrix[4] + mean_diffs[1] * mean_diffs[1] *
    value) / m_particles;
  m_covariance_matrix[5] = (a.m_particles * a.m_covariance_matrix[5] +
    b.m_particles * b.m_covariance_matrix[5] + mean_diffs[1] * mean_diffs[2] *
    value) / m_particles;
  m_covariance_matrix[6] = (a.m_particles * a.m_covariance_matrix[6] +
    b.m_particles * b.m_covariance_matrix[6] + mean_diffs[1] * mean_diffs[3] *
    value) / m_particles;
  m_covariance_matrix[7] = (a.m_particles * a.m_covariance_matrix[7] +
    b.m_particles * b.m_covariance_matrix[7] + mean_diffs[2] * mean_diffs[2] *
    value) / m_particles;
  m_covariance_matrix[8] = (a.m_particles * a.m_covariance_matrix[8] +
    b.m_particles * b.m_covariance_matrix[8] + mean_diffs[2] * mean_diffs[3] *
    value) / m_particles;
  m_covariance_matrix[9] = (a.m_particles * a.m_covariance_matrix[9] +
    b.m_particles * b.m_covariance_matrix[9] + mean_diffs[3] * mean_diffs[3] *
    value) / m_particles;
}

void Statistics::writeToFile(std::ofstream& file)
{
  std::array<double, 14> values = {
    m_means[0], m_means[1], m_means[2], m_means[3], m_covariance_matrix[0],
    m_covariance_matrix[1], m_covariance_matrix[2], m_covariance_matrix[3],
    m_covariance_matrix[4], m_covariance_matrix[5], m_covariance_matrix[6],
    m_covariance_matrix[7], m_covariance_matrix[8], m_covariance_matrix[9]
  };
  file.write(reinterpret_cast<char*>(values.data()), sizeof(values));
}

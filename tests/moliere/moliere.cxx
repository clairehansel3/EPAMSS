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

#include "moliere.hxx"
#include "random.hxx"
#include <ctime>
#include <fstream>
#include <iostream>
#include <stdexcept>

int main(void)
{
  seedRandom(0, std::time(nullptr));
  double integration_tolerance = 1E-10;
  double B = 10;
  double vartheta_cutoff = 10;
  unsigned long spline_points = 1000;
  unsigned max_order = 3;
  unsigned max_integration_depth = 15;
  double approx_pdf_points = 10000;
  std::size_t samples = 1000000;
  Moliere moliere{max_order, max_integration_depth, integration_tolerance,
    spline_points, vartheta_cutoff, true};
  try {
    std::cout << "computing pdf" << std::endl;
    std::ofstream file1;
    file1.exceptions(file1.failbit | file1.badbit);
    file1.open("pdf");
    double step = vartheta_cutoff / approx_pdf_points;
    for (double vartheta = 0; vartheta < vartheta_cutoff; vartheta += step)
      file1 << vartheta << " " << moliere.pdf(vartheta, B, 0) << " " << \
        moliere.pdf(vartheta, B, 1) << " " << moliere.pdf(vartheta, B, 2) << \
        " " << moliere.pdf(vartheta, B, 3) << std::endl;
    file1.close();
    std::cout << "generating samples" << std::endl;
    std::ofstream file2;
    file2.exceptions(file2.failbit | file2.badbit);
    file2.open("samples");
    for (std::size_t i = 0; i != samples; ++i)
      file2 << moliere.sample(B, 0) << " " << moliere.sample(B, 1) << " " << \
        moliere.sample(B, 2) << " " << moliere.sample(B, 3) << std::endl;
    file2.close();
  } catch (const std::ifstream::failure&) {
    throw std::runtime_error("unable to write output file");
  }
}

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

#ifndef EPAMSS_MOLIERE_HXX
#define EPAMSS_MOLIERE_HXX

#include <boost/math/interpolators/cubic_b_spline.hpp>
#include <functional>
#include <memory>
#include <optional>
#include <vector>

class Moliere {

public:

  Moliere(unsigned max_order, unsigned max_integration_depth,
          double integration_tolerance, unsigned long spline_points,
          double vartheta_cutoff, bool print_info);

  double pdf(double theta, double B);

  double sample(double B);

  double pdf(double theta, double B, unsigned max_order);

  double sample(double B, unsigned max_order);

private:

  class Spline {

  public:

    Spline(std::function<double(double)> function, double cutoff,
           unsigned long spline_points);

    double operator()(double);

  private:

    std::unique_ptr<double[]> m_table;

    std::optional<boost::math::cubic_b_spline<double>> m_spline;

  };

  unsigned m_max_order;

  double m_vartheta_cutoff;

  std::vector<Spline> m_splines;

};

#endif

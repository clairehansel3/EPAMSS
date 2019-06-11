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
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <memory>

static double compute_f_0(double theta)
{
  return 2 * theta * std::exp(-theta * theta);
}

static double compute_f_n(double theta, unsigned order,
                          unsigned max_integration_depth,
                          double integration_tolerance)
{
  assert(order >= 1);
  auto function = [theta, order](double eta){
    double x = eta * eta / 4;
    return std::exp(-x) * std::pow(x * std::log(x), order) *
      boost::math::cyl_bessel_j(0, eta * theta) * eta;
  };
  double error;
  double value = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(
    function, 0, std::numeric_limits<double>::infinity(), max_integration_depth,
    integration_tolerance, &error
  );
  return theta * value / boost::math::factorial<double>(order);
}

Moliere::Moliere(unsigned max_order, unsigned max_integration_depth,
                 double integration_tolerance, unsigned long spline_points,
                 double vartheta_cutoff, bool print_info)
: m_max_order{max_order},
  m_vartheta_cutoff{vartheta_cutoff}
{
  for (unsigned order = 1; order <= max_order; ++order) {
    if (print_info)
      std::cout << "constructing spline for f_" << order << std::endl;
    auto f_n = [=](double theta){
      return compute_f_n(theta, order, max_integration_depth,
                         integration_tolerance);
    };
    m_splines.emplace_back(f_n, vartheta_cutoff, spline_points);
  }
  if (max_order && print_info)
    std::cout << "done constructing splines" << std::endl;
}

double Moliere::pdf(double theta, double B)
{
  double value = compute_f_0(theta);
  if (theta < m_vartheta_cutoff)
    for (unsigned order = 1; order <= m_max_order; ++order) {
      value += m_splines[order - 1](theta)
               * std::pow(B, -static_cast<double>(order));
    }
  return value;
}

double Moliere::pdf(double theta, double B, unsigned max_order)
{
  double value = compute_f_0(theta);
  if (theta < m_vartheta_cutoff)
    for (unsigned order = 1; order <= max_order; ++order) {
      value += m_splines[order - 1](theta)
               * std::pow(B, -static_cast<double>(order));
    }
  return value;
}

double Moliere::sample(double B)
{
  if (m_max_order == 0) {
    double vartheta = std::sqrt(-std::log(1 - randomUniform()));
    assert(std::isfinite(vartheta));
    return vartheta;
  }
  while (true) {
    double vartheta = m_vartheta_cutoff * randomUniform();
    if (randomUniform() <= pdf(vartheta, B)) {
      return vartheta;
    }
  }
}


double Moliere::sample(double B, unsigned max_order)
{
  if (max_order == 0) {
    double vartheta = std::sqrt(-std::log(1 - randomUniform()));
    assert(std::isfinite(vartheta));
    return vartheta;
  }
  while (true) {
    double vartheta = m_vartheta_cutoff * randomUniform();
    if (randomUniform() <= pdf(vartheta, B, max_order)) {
      return vartheta;
    }
  }
}

Moliere::Spline::Spline(std::function<double(double)> function, double cutoff,
                        unsigned long spline_points)
: m_table{std::make_unique<double[]>(spline_points)}
{
  assert(cutoff > 0);
  double step_size = cutoff / (spline_points - 1);
  for (unsigned long i = 0; i != spline_points; ++i)
    m_table[i] = function(i * step_size);
  m_spline = boost::math::cubic_b_spline<double>(
    m_table.get(), spline_points, 0.0, step_size
  );
}

double Moliere::Spline::operator()(double x)
{
  return (*m_spline)(x);
}

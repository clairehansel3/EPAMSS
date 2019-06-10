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

#include "input.hxx"
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <ctime>
#include <fstream>
#include <map>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>

static std::optional<std::string> readFile(const char* path)
// Reads the file located at 'path' into a string. If there were any errors
// opening or reading the file, an empty optional is returned. Otherwise, an
// optional containing the string is returned.
{
  assert(path);
  try {
    std::ifstream file;
    file.exceptions(file.failbit | file.badbit);
    file.open(path);
    file.seekg(0, file.end);
    auto file_size = boost::numeric_cast<std::size_t>(std::streamoff(file.tellg()));
    file.seekg(0, file.beg);
    std::string result;
    result.resize(file_size);
    file.read(&result[0], file_size);
    return result;
  } catch (const std::ifstream::failure&) {
    return std::nullopt;
  }
}

static std::map<std::string, std::string> inputToDict(std::string input)
// Reads a string containing a number of lines of the form <key> = <value> into
// a dictionary mapping keys to values. The exact grammar is:
// inputfile      := { [ key-value-pair ] , "\n" } , [ key-value-pair ] , EOF ;
// key-value-pair := word , "=" , word ;
// word           := { space } , { character } , { space } ;
// space          := " " | "\t" | "\v" | "\f" | "\r" ;
// character      := ( byte - space - "\n" - "=") ;
{
  std::map<std::string, std::string> dict;
  auto iterator = input.begin();

  // iterate through the lines
  while (true) {

    // find start of key
    while (true) {
      if (iterator == input.end())
        return dict;
      if (*iterator == '=')
        throw std::runtime_error("bad input file: expected key before =");
      if (!std::isspace(*iterator))
        break;
      ++iterator;
    }
    auto key_begin = iterator;
    ++iterator;

    // find end of key
    while (true) {
      if (iterator == input.end() || *iterator == '\n')
        throw std::runtime_error("bad input file: expected = after key");
      if (std::isspace(*iterator) || *iterator == '=')
        break;
      ++iterator;
    }
    auto key_end = iterator;
    while (*iterator != '=') {
      if (iterator == input.end() || *iterator == '\n' || !std::isspace(*iterator))
        throw std::runtime_error("bad input file: expected = after key");
      ++iterator;
    }
    ++iterator;

    // find start of value
    while (true) {
      if (iterator == input.end() || *iterator == '\n')
        throw std::runtime_error("bad input file: expected value after =");
      if (*iterator == '=')
        throw std::runtime_error("bad input file: more than one = on line");
      if (!std::isspace(*iterator))
        break;
      ++iterator;
    }
    auto value_begin = iterator;

    // find end of value
    while (true) {
      if (iterator == input.end() || std::isspace(*iterator))
        break;
      if (*iterator == '=')
        throw std::runtime_error("bad input file: = in value");
      ++iterator;
    }
    auto value_end = iterator;

    // add value to dictionary
    std::string key{key_begin, key_end};
    std::string value{value_begin, value_end};
    if (dict.find(key) != dict.end())
      throw std::runtime_error("bad input file: duplicate keys");
    dict[key] = value;

    // advance to end of line
    while (true) {
      if (iterator == input.end())
        break;
      if (*iterator == '\n') {
        ++iterator;
        break;
      }
      if (!std::isspace(*iterator))
        throw std::runtime_error("bad input file: line must end after value");
      ++iterator;
    }

  }
}

static void get_from_dict(std::map<std::string, std::string>& dict,
  std::string name, std::string& value)
{
  auto value_iter = dict.find(name);
  if (value_iter == dict.end())
    throw std::runtime_error(std::string("parameter ") + name +
      " not in input file");
  value = value_iter->second;
  dict.erase(value_iter);
}

static void get_from_dict(std::map<std::string, std::string>& dict,
  std::string name, double& value)
{
  std::string value_str;
  get_from_dict(dict, std::move(name), value_str);
  try {
    value = std::stod(value_str);
  } catch (const std::logic_error&) {
    throw std::runtime_error(std::string("unable to convert value of parameter "
      ) + name + " to real number");
  }
}

template <typename T>
static std::enable_if_t<std::is_unsigned_v<T>, void> get_from_dict(
  std::map<std::string, std::string>& dict, std::string name, T& value)
{
  std::string value_str;
  get_from_dict(dict, std::move(name), value_str);
  try {
    value = boost::numeric_cast<T>(std::stoull(value_str));
  } catch (const std::logic_error&) {
    throw std::runtime_error(std::string("unable to convert value of parameter "
    ) + name + " to (unsigned) integer");
  } catch (boost::numeric::bad_numeric_cast&) {
    throw std::runtime_error(std::string("value of parameter ") + name +
      "is too large for c++ type of corresponding parameter");
  }
}

static void get_from_dict(std::map<std::string, std::string>& dict,
  std::string name, bool& value)
{
  std::string value_str;
  get_from_dict(dict, std::move(name), value_str);
  if (value_str == "true" || value_str == "True")
    value = true;
  else if (value_str == "false" || value_str == "False")
    value = false;
  else
    throw std::runtime_error(std::string("value of parameter ") + name +
      " must be either true or false");
}

#define EPAMSS_READ_PARAMETER(dict, param) get_from_dict(dict, #param, param)


Parameters::Parameters(const char* input_path, boost::mpi::communicator& world)
{
  auto input = readFile(input_path);
  if (!input)
    throw std::runtime_error("unable to open input");
  auto dict = inputToDict(*input);

  EPAMSS_READ_PARAMETER(dict, rho_ion_initial_si);
  EPAMSS_READ_PARAMETER(dict, plasma_length_si);
  EPAMSS_READ_PARAMETER(dict, beam_energy_initial_gev);
  EPAMSS_READ_PARAMETER(dict, acceleration_gradient_gev_per_m);
  EPAMSS_READ_PARAMETER(dict, bennett_radius_initial_si);
  EPAMSS_READ_PARAMETER(dict, cross_section_radius_si);
  EPAMSS_READ_PARAMETER(dict, unperturbed_plasma_density_si);
  EPAMSS_READ_PARAMETER(dict, integration_tolerance);
  EPAMSS_READ_PARAMETER(dict, vartheta_cutoff);
  EPAMSS_READ_PARAMETER(dict, ion_atomic_number);
  EPAMSS_READ_PARAMETER(dict, minimum_steps_per_betatron_period);
  EPAMSS_READ_PARAMETER(dict, particles_target);
  EPAMSS_READ_PARAMETER(dict, analysis_points_target);
  EPAMSS_READ_PARAMETER(dict, spline_points);

  try {
    std::string seed_str;
    get_from_dict(dict, "seed", seed_str);
    seed = boost::numeric_cast<unsigned>(std::stoull(seed_str));
  } catch (std::runtime_error&) {
    if (world.rank() == 0) {
      seed = std::time(nullptr);
      for (int i = 1; i != world.size(); ++i) {
        world.send(i, 0, &seed, 1);
      }
    } else {
      world.recv(0, 0, &seed, 1);
    }
  }
  EPAMSS_READ_PARAMETER(dict, max_order);
  EPAMSS_READ_PARAMETER(dict, max_integration_depth);
  EPAMSS_READ_PARAMETER(dict, output_filename);
  EPAMSS_READ_PARAMETER(dict, statistics_filename);
  EPAMSS_READ_PARAMETER(dict, phase_space_filename);
  EPAMSS_READ_PARAMETER(dict, output_phase_space);
  EPAMSS_READ_PARAMETER(dict, modified_bennett);

  if (world.rank() == 0) {
    for (auto pair : dict) {
      std::cerr << "WARNING: extraneous parameter '" << pair.first <<
        "' defined in input file with value '" << pair.second << "'" <<
        std::endl;
    }
  }

  computeDependentParameters(world.size());
}

#define EPAMSS_WRITE(file, param, units) file << #param " = " << param << " " units "\n"
#define EPAMSS_WRITE_BOOL(file, param) file << #param " = " << (param ? "true\n" : "false\n")

void Parameters::writeOutputFile(double seconds)
{
  try {
    std::ofstream file;
    file.exceptions(file.failbit | file.badbit);
    file.open(output_filename);

    EPAMSS_WRITE(file, rho_ion_initial_si, "m^-3");
    EPAMSS_WRITE(file, plasma_length_si, "m");
    EPAMSS_WRITE(file, beam_energy_initial_gev, "GeV");
    EPAMSS_WRITE(file, acceleration_gradient_gev_per_m, "GeV/m");
    EPAMSS_WRITE(file, bennett_radius_initial_si, "m");
    EPAMSS_WRITE(file, cross_section_radius_si, "m");
    EPAMSS_WRITE(file, unperturbed_plasma_density_si, "m^-3");
    EPAMSS_WRITE(file, integration_tolerance, "");
    EPAMSS_WRITE(file, vartheta_cutoff, "");
    EPAMSS_WRITE(file, ion_atomic_number, "");
    EPAMSS_WRITE(file, minimum_steps_per_betatron_period, "");
    EPAMSS_WRITE(file, particles_target, "");
    EPAMSS_WRITE(file, analysis_points_target, "");
    EPAMSS_WRITE(file, spline_points, "");
    EPAMSS_WRITE(file, seed, "");
    EPAMSS_WRITE(file, max_order, "");
    EPAMSS_WRITE(file, max_integration_depth, "");
    EPAMSS_WRITE(file, output_filename, "");
    EPAMSS_WRITE(file, statistics_filename, "");
    EPAMSS_WRITE(file, phase_space_filename, "");
    EPAMSS_WRITE_BOOL(file, output_phase_space);
    EPAMSS_WRITE_BOOL(file, modified_bennett);

    EPAMSS_WRITE(file, plasma_frequency_si, "s^-1");
    EPAMSS_WRITE(file, plasma_angular_wavenumber_si, "m^-1");
    EPAMSS_WRITE(file, plasma_skin_depth_si, "m");
    EPAMSS_WRITE(file, rho_ion_initial, "");
    EPAMSS_WRITE(file, plasma_length, "");
    EPAMSS_WRITE(file, gamma_initial, "");
    EPAMSS_WRITE(file, gamma_prime, "");
    EPAMSS_WRITE(file, bennett_radius_initial, "");
    EPAMSS_WRITE(file, cross_section_radius, "");
    EPAMSS_WRITE(file, delta, "");
    EPAMSS_WRITE(file, gamma_final, "");
    EPAMSS_WRITE(file, bennett_radius_final, "");
    EPAMSS_WRITE(file, rho_ion_final, "");
    EPAMSS_WRITE(file, betatron_frequency_final, "");
    EPAMSS_WRITE(file, betatron_period_final, "");
    EPAMSS_WRITE(file, betatron_frequency_final_si, "m^-1");
    EPAMSS_WRITE(file, betatron_period_final_si, "m");
    EPAMSS_WRITE(file, step_size, "");
    EPAMSS_WRITE(file, step_size_si, "m");
    EPAMSS_WRITE(file, gamma_minimum_angle, "");
    EPAMSS_WRITE(file, omega_off_axis, "");
    EPAMSS_WRITE(file, omega_on_axis_initial, "");
    EPAMSS_WRITE(file, max_scattering_r_div_a_initial, "");
    EPAMSS_WRITE(file, sigma_r, "");
    EPAMSS_WRITE(file, sigma_r_prime_initial, "");
    EPAMSS_WRITE(file, steps, "");
    EPAMSS_WRITE(file, stride, "");
    EPAMSS_WRITE(file, compute_processes, "");
    EPAMSS_WRITE(file, actual_particles, "");
    EPAMSS_WRITE(file, particles_per_process, "");
    EPAMSS_WRITE(file, actual_analysis_points, "");

    file << "seconds_elapsed = " << seconds << '\n';
    file << "minutes_elapsed = " << (seconds / 60) << '\n';
    file << "hours_elapsed = " << (seconds / (60 * 60)) << '\n';
    file << "approx_core_hours_elapsed = " << ((compute_processes * seconds) /
      (60 * 60)) << '\n';
    file.close();
  } catch (const std::ifstream::failure&) {
    throw std::runtime_error("unable to write output file");
  }
}

void Parameters::computeDependentParameters(int processes)
{
  constexpr double elementary_charge = 1.60217662E-19;
  constexpr double vacuum_permittivity = 8.8541878128E-12;
  constexpr double electron_mass = 9.1093837015E-31;
  constexpr double c_light = 299792458.0;
  constexpr double electron_rest_energy_mev = 0.5109989461;
  constexpr double classical_electron_radius = 2.8179403227E-15;
  plasma_frequency_si = elementary_charge * std::sqrt(unperturbed_plasma_density_si / (vacuum_permittivity * electron_mass));
  plasma_angular_wavenumber_si = plasma_frequency_si / c_light;
  plasma_skin_depth_si = 1 / plasma_angular_wavenumber_si;
  rho_ion_initial = rho_ion_initial_si / unperturbed_plasma_density_si;
  plasma_length = plasma_length_si * plasma_angular_wavenumber_si;
  gamma_initial = beam_energy_initial_gev * 1000 / electron_rest_energy_mev;
  gamma_prime = plasma_skin_depth_si * acceleration_gradient_gev_per_m * 1000 / electron_rest_energy_mev;
  assert(gamma_prime >= 0);
  bennett_radius_initial = bennett_radius_initial_si * plasma_angular_wavenumber_si;
  cross_section_radius = cross_section_radius_si * plasma_angular_wavenumber_si;
  delta = modified_bennett ? 1.0 : 0.0;
  double gamma_final_over_gamma_initial = 1 + (gamma_prime * plasma_length / gamma_initial);
  gamma_final = gamma_initial * gamma_final_over_gamma_initial;
  bennett_radius_final = bennett_radius_initial * std::pow(gamma_final_over_gamma_initial, -0.25);
  rho_ion_final = rho_ion_initial * std::sqrt(gamma_final_over_gamma_initial);
  betatron_frequency_final = std::sqrt((ion_atomic_number * (rho_ion_final + delta) / (2 * gamma_final)) - (gamma_prime * gamma_prime / (4 * gamma_final * gamma_final)));
  betatron_period_final = 2 * boost::math::constants::pi<double>() / betatron_frequency_final;
  betatron_frequency_final_si = plasma_angular_wavenumber_si * betatron_frequency_final;
  betatron_period_final_si = plasma_skin_depth_si * betatron_period_final;
  steps = static_cast<std::size_t>(std::ceil(
    minimum_steps_per_betatron_period * plasma_length * betatron_frequency_final
    / (2 * boost::math::constants::pi<double>())
  ));
  stride = steps / analysis_points_target;
  step_size = plasma_length / steps;
  step_size_si = step_size * plasma_skin_depth_si;
  gamma_minimum_angle = 2 * ion_atomic_number * classical_electron_radius / cross_section_radius_si;
  omega_off_axis = boost::math::constants::pi<double>() * cross_section_radius_si * cross_section_radius_si * step_size_si * unperturbed_plasma_density_si;
  omega_on_axis_initial = omega_off_axis * (rho_ion_initial + delta);
  if (omega_on_axis_initial < 25) {
    max_scattering_r_div_a_initial = 0;
  }
  else if (delta * omega_off_axis > 25) {
    max_scattering_r_div_a_initial = std::numeric_limits<double>::infinity();
  }
  else {
    max_scattering_r_div_a_initial = std::sqrt(std::sqrt(rho_ion_initial * omega_off_axis / (25 - delta * omega_off_axis)) - 1);
  }
  sigma_r = 0.5 * bennett_radius_initial * std::sqrt(rho_ion_initial);
  sigma_r_prime_initial = std::sqrt(boost::math::constants::pi<double>() * classical_electron_radius * bennett_radius_initial_si * bennett_radius_initial_si * ion_atomic_number * rho_ion_initial_si / (2 * gamma_initial));
  compute_processes = static_cast<std::size_t>(processes) - 1;
  actual_particles = particles_target + compute_processes - particles_target % compute_processes;
  particles_per_process = actual_particles / compute_processes;
  assert(actual_particles % compute_processes == 0);
  actual_analysis_points = (steps / stride) + 1;
}

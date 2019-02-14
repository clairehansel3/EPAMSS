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


Parameters::Parameters(const char* input_path, int processes)
{
  auto input = readFile(input_path);
  if (!input)
    throw std::runtime_error("unable to open input");
  auto dict = inputToDict(*input);

  EPAMSS_READ_PARAMETER(dict, maximum_ion_density);
  EPAMSS_READ_PARAMETER(dict, plasma_length);
  EPAMSS_READ_PARAMETER(dict, beam_energy);
  EPAMSS_READ_PARAMETER(dict, electron_linear_density);
  EPAMSS_READ_PARAMETER(dict, bennett_radius);
  EPAMSS_READ_PARAMETER(dict, interaction_radius);
  EPAMSS_READ_PARAMETER(dict, integration_tolerance);
  EPAMSS_READ_PARAMETER(dict, vartheta_cutoff);
  EPAMSS_READ_PARAMETER(dict, ion_atomic_number);
  EPAMSS_READ_PARAMETER(dict, minimum_steps_per_betatron_period);
  EPAMSS_READ_PARAMETER(dict, particles);
  EPAMSS_READ_PARAMETER(dict, analysis_points_target);
  EPAMSS_READ_PARAMETER(dict, spline_points);
  try {
    std::string seed_str;
    get_from_dict(dict, "seed", seed_str);
    seed = boost::numeric_cast<unsigned>(std::stoull(seed_str));
  } catch (std::runtime_error&) {
    seed = std::time(NULL);
  }
  EPAMSS_READ_PARAMETER(dict, max_order);
  EPAMSS_READ_PARAMETER(dict, max_integration_depth);
  EPAMSS_READ_PARAMETER(dict, output_filename);
  EPAMSS_READ_PARAMETER(dict, statistics_filename);
  EPAMSS_READ_PARAMETER(dict, phase_space_filename);
  EPAMSS_READ_PARAMETER(dict, output_phase_space);

  if (dict.size() != 0) {
    std::string error_msg("unknown extra parameters defined in input file: ");
    for (auto pair : dict) {
      error_msg += pair.first;
      error_msg += ' ';
    }
    throw std::runtime_error(error_msg);
  }

  computeDependentParameters(processes);
}

#define EPAMSS_WRITE(file, param) file << #param " = " << param << '\n'

void Parameters::writeOutputFile(double seconds)
{
  try {
    std::ofstream file;
    file.exceptions(file.failbit | file.badbit);
    file.open(output_filename);
    EPAMSS_WRITE(file, maximum_ion_density);
    EPAMSS_WRITE(file, plasma_length);
    EPAMSS_WRITE(file, beam_energy);
    EPAMSS_WRITE(file, electron_linear_density);
    EPAMSS_WRITE(file, bennett_radius);
    EPAMSS_WRITE(file, interaction_radius);
    EPAMSS_WRITE(file, integration_tolerance);
    EPAMSS_WRITE(file, vartheta_cutoff);
    EPAMSS_WRITE(file, ion_atomic_number);
    EPAMSS_WRITE(file, minimum_steps_per_betatron_period);
    EPAMSS_WRITE(file, particles);
    EPAMSS_WRITE(file, analysis_points_target);
    EPAMSS_WRITE(file, spline_points);
    EPAMSS_WRITE(file, seed);
    EPAMSS_WRITE(file, max_order);
    EPAMSS_WRITE(file, max_integration_depth);
    EPAMSS_WRITE(file, output_filename);
    EPAMSS_WRITE(file, statistics_filename);
    EPAMSS_WRITE(file, phase_space_filename);
    file << "output_phase_space = " << (output_phase_space ? "true" : "false")
      << '\n';
    EPAMSS_WRITE(file, ion_linear_density);
    EPAMSS_WRITE(file, gamma);
    EPAMSS_WRITE(file, alpha);
    EPAMSS_WRITE(file, step_size);
    EPAMSS_WRITE(file, cross_section);
    EPAMSS_WRITE(file, minimum_angle);
    EPAMSS_WRITE(file, betatron_frequency);
    EPAMSS_WRITE(file, betatron_period);
    EPAMSS_WRITE(file, max_scattering_r_div_a);
    EPAMSS_WRITE(file, percent_with_scattering);
    EPAMSS_WRITE(file, omega_on_axis);
    EPAMSS_WRITE(file, steps);
    EPAMSS_WRITE(file, stride);
    EPAMSS_WRITE(file, actual_particles);
    EPAMSS_WRITE(file, particles_per_process);
    EPAMSS_WRITE(file, compute_processes);
    EPAMSS_WRITE(file, analysis_points);
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
  constexpr double electron_rest_energy_mev = 0.5109989461;
  constexpr double classical_electron_radius = 2.8179403227E-15;
  ion_linear_density = boost::math::constants::pi<double>() * bennett_radius *
    bennett_radius * maximum_ion_density;
  gamma = beam_energy * 1000 / electron_rest_energy_mev;
  alpha = 2 * classical_electron_radius * (ion_atomic_number *
    ion_linear_density - electron_linear_density / (gamma * gamma)) / gamma;
  steps = static_cast<std::size_t>(std::ceil(
    minimum_steps_per_betatron_period * plasma_length * std::sqrt(alpha) /
    (2 * boost::math::constants::pi<double>() * bennett_radius)
  ));
  stride = steps / analysis_points_target;
  step_size = plasma_length / steps;
  cross_section = boost::math::constants::pi<double>() *
    interaction_radius * interaction_radius;
  minimum_angle = 2 * ion_atomic_number * classical_electron_radius / (gamma *
    interaction_radius);
  compute_processes = static_cast<std::size_t>(processes) - 1;
  actual_particles = particles + compute_processes - particles % compute_processes;
  particles_per_process = actual_particles / compute_processes;
  assert(actual_particles % compute_processes == 0);
  betatron_frequency = std::sqrt(alpha) / bennett_radius;
  betatron_period = 2 * boost::math::constants::pi<double>() /
    betatron_frequency;
  omega_on_axis = maximum_ion_density * interaction_radius * step_size;
  max_scattering_r_div_a = std::sqrt((std::sqrt(omega_on_axis) / 5) - 1);
  percent_with_scattering = 1 - 5 / std::sqrt(omega_on_axis);
  analysis_points = (steps / stride) + 1;
}

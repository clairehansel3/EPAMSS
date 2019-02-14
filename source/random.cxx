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

#include "random.hxx"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <ctime>

static boost::random::mt19937 twister;

void seedRandom(unsigned thread_rank)
{
  twister.seed(std::time(0));
  boost::random::uniform_int_distribution<unsigned> dist;
  auto seed = dist(twister) + thread_rank;
  twister.seed(seed);
}

double randomUniform()
{
  boost::random::uniform_real_distribution<double> dist;
  return dist(twister);
}


double randomNormal(double mean, double std)
{
  boost::random::normal_distribution dist{mean, std};
  return dist(twister);
}

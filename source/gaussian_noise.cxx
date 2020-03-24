#include "gaussian_noise.hxx"
#include "random.hxx"

std::unique_ptr<double[]> get_density(double density, std::size_t n, double factor)
{
  auto densities = std::unique_ptr<double[]>(new double[n]);
  for (std::size_t i = 0; i != n; ++i) {
    densities[i] = density * (1 + randomNormal(0.0, factor));
  }
  return densities;
}

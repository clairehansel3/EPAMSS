#ifndef EPAMSS_GAUSSIAN_NOISE_HXX
#define EPAMSS_GAUSSIAN_NOISE_HXX

#include <cstddef>
#include <memory>

std::unique_ptr<double[]> get_density(double density, std::size_t n, double factor);

#endif

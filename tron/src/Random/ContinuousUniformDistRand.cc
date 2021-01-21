#include "Random/ContinuousUniformDistRand.hh"

Tron::ContinuousUniformDistRand::ContinuousUniformDistRand() {
  std::random_device seed;
  fMersenneTwister         = std::mt19937(seed());
  fUniformRealDistribution = std::uniform_real_distribution<>(0.0, 1.0);
}

Tron::ContinuousUniformDistRand::ContinuousUniformDistRand(unsigned int seed) {
  fMersenneTwister         = std::mt19937(seed);
  fUniformRealDistribution = std::uniform_real_distribution<>(0.0, 1.0);
}

void Tron::ContinuousUniformDistRand::SetParams(double xmin, double xmax) {
  fUniformRealDistribution = std::uniform_real_distribution<>(xmin, xmax);
}

Tron::ContinuousUniformDistRand::ResultType Tron::ContinuousUniformDistRand::operator()() {
  return fUniformRealDistribution(fMersenneTwister);
}

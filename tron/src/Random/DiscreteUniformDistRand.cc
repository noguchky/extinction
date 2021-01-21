#include "Random/DiscreteUniformDistRand.hh"

Tron::DiscreteUniformDistRand::DiscreteUniformDistRand() {
  std::random_device seed;
  fMersenneTwister        = std::mt19937(seed());
  fUniformIntDistribution = std::uniform_int_distribution<>(0, 1);
}

Tron::DiscreteUniformDistRand::DiscreteUniformDistRand(unsigned int seed) {
  fMersenneTwister        = std::mt19937(seed);
  fUniformIntDistribution = std::uniform_int_distribution<>(0, 1);
}

void Tron::DiscreteUniformDistRand::SetParams(int xmin, int xmax) {
  fUniformIntDistribution = std::uniform_int_distribution<>(xmin, xmax);
}

Tron::DiscreteUniformDistRand::ResultType Tron::DiscreteUniformDistRand::operator()() {
  return fUniformIntDistribution(fMersenneTwister);
}

#include "Random/NormalDistRand.hh"
#include "Math.hh"

Tron::NormalDistRand::NormalDistRand() {
  std::random_device seed;
  fMersenneTwister    = std::mt19937(seed());
  fNormalDistribution = std::normal_distribution<>(0.0, 1.0);
}

Tron::NormalDistRand::NormalDistRand(unsigned int seed) {
  fMersenneTwister    = std::mt19937(seed);
  fNormalDistribution = std::normal_distribution<>(0.0, 1.0);
}

void Tron::NormalDistRand::SetParams(double mean, double sigma) {
  fNormalDistribution = std::normal_distribution<>(mean, sigma);
}
  
Tron::NormalDistRand::ResultType Tron::NormalDistRand::operator()() {
  return fNormalDistribution(fMersenneTwister);
}

Tron::NormalDistRand0::NormalDistRand0() {
  std::random_device seed;
  fContinuousUniformDistRand = ContinuousUniformDistRand(seed());
  fMean  = 0.0;
  fSigma = 1.0;
}

Tron::NormalDistRand0::NormalDistRand0(unsigned int seed) {
  fContinuousUniformDistRand = ContinuousUniformDistRand(seed);
  fMean  = 0.0;
  fSigma = 1.0;
}

void Tron::NormalDistRand0::SetParams(double mean, double sigma) {
  fMean  = mean;
  fSigma = sigma;
}

Tron::NormalDistRand0::ResultType Tron::NormalDistRand0::operator()() {
  const double x = fContinuousUniformDistRand();
  const double y = fContinuousUniformDistRand();
  const double value = std::sqrt(-2.0 * std::log(x)) * std::sin(twopi * y); // Box Muller Method
  return fMean + fSigma * value;
}

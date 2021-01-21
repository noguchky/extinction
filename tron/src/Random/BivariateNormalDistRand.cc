#include <random>
#include "TMath.h"
#include "Random/BivariateNormalDistRand.hh"

Tron::BivariateNormalDistRand::BivariateNormalDistRand() {
  std::random_device seed;
  fNormalDistRand = NormalDistRand(seed());
  fMeanX       = 0.0;
  fMeanY       = 0.0;
  fSigmaX      = 1.0;
  fSigmaY      = 1.0;
  fCorrelation = 0.0;
}

Tron::BivariateNormalDistRand::BivariateNormalDistRand(unsigned int seed) {
  fNormalDistRand = NormalDistRand(seed);
  fMeanX       = 0.0;
  fMeanY       = 0.0;
  fSigmaX      = 1.0;
  fSigmaY      = 1.0;
  fCorrelation = 0.0;
}

void Tron::BivariateNormalDistRand::SetParams(double meanX, double meanY, double sigmaX, double sigmaY, double correlation) {
  fMeanX       = meanX;
  fMeanY       = meanY;
  fSigmaX      = sigmaX;
  fSigmaY      = sigmaY;
  fCorrelation = correlation;
}
  
Tron::BivariateNormalDistRand::ResultType Tron::BivariateNormalDistRand::operator()() {
  ResultType result;

  if (fSigmaX) {
    // Create X
    const double meanX  = fMeanX;
    const double sigmaX = fSigmaX;
    result.X = meanX + sigmaX * fNormalDistRand();

    // Create Y from X (Y Depends on X)
    const double meanY  = fMeanY + fCorrelation * fSigmaY / fSigmaX * (result.X - fMeanX);
    const double sigmaY = TMath::Sqrt(1.0 - fCorrelation * fCorrelation) * fSigmaY;
    result.Y = meanY + sigmaY * fNormalDistRand();
  } else {
    // Create X
    result.X = fMeanX;
    
    // Create Y
    result.Y = fMeanY + fSigmaY * fNormalDistRand();
  }
  
  return result;
}

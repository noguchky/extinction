#include <random>
#include "TMath.h"
#include "Random/BivariateUniformDistRand.hh"
#include "Ellipse.hh"

Tron::BivariateUniformDistRand::BivariateUniformDistRand() {
  std::random_device seed;
  fUniformDistRand = ContinuousUniformDistRand(seed());
  fMeanX       = 0.0;
  fMeanY       = 0.0;
  fSigmaX      = 1.0;
  fSigmaY      = 1.0;
  fCorrelation = 0.0;
}

Tron::BivariateUniformDistRand::BivariateUniformDistRand(unsigned int seed) {
  fUniformDistRand = ContinuousUniformDistRand(seed);
  fMeanX       = 0.0;
  fMeanY       = 0.0;
  fSigmaX      = 1.0;
  fSigmaY      = 1.0;
  fCorrelation = 0.0;
}

void Tron::BivariateUniformDistRand::SetParams(double meanX, double meanY, double sigmaX, double sigmaY, double correlation) {
  fMeanX       = meanX;
  fMeanY       = meanY;
  fSigmaX      = sigmaX;
  fSigmaY      = sigmaY;
  fCorrelation = correlation;
}

Tron::BivariateUniformDistRand::ResultType Tron::BivariateUniformDistRand::operator()() {
  ResultType result;

  if (fSigmaX == 0 || fSigmaY == 0) {
    result.X = fMeanX + fSigmaX * (2.0 * fUniformDistRand() - 1.0);
    result.Y = fMeanY + fSigmaY * (2.0 * fUniformDistRand() - 1.0);
  } else {
    const double theta = fUniformDistRand() * TMath::TwoPi();
    const double rho   = TMath::Sqrt(fUniformDistRand());

    result.X = rho * TMath::Cos(theta);
    result.Y = rho * TMath::Sin(theta);

    result.X *= Ellipse::GetSemiAxisX(fSigmaX, fSigmaY, fCorrelation);
    result.Y *= Ellipse::GetSemiAxisY(fSigmaX, fSigmaY, fCorrelation);

    result.Rotate(Ellipse::GetRotationAngle(fSigmaX, fSigmaY, fCorrelation));

    result.X += fMeanX;
    result.Y += fMeanY;
  }

  return result;
}

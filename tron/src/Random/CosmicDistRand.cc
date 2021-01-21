#include <cmath>
#include <random>
#include "Random/CosmicDistRand.hh"
#include "Math.hh"

Tron::CosmicDistRand::CosmicDistRand(){
  std::random_device seed;
  fNormalDistRand = NormalDistRand(seed());
  fSourceRadius   = 1.0;
  fSourceDistance = 1.0;
}

Tron::CosmicDistRand::CosmicDistRand(unsigned int seed) {
  fNormalDistRand = NormalDistRand(seed);
  fSourceRadius   = 1.0;
  fSourceDistance = 1.0;
}

void Tron::CosmicDistRand::SetParams(double radius, double distance) {
  fSourceRadius   = radius;
  fSourceDistance = distance;
}

Tron::CosmicDistRand::ResultType Tron::CosmicDistRand::operator()() {
  ResultType result;

  const double sourceRho = fSourceRadius * std::sqrt(fNormalDistRand());
  const double sourcePhi = twopi * fNormalDistRand();
      
  result.Position.X = sourceRho * std::cos(sourcePhi);
  result.Position.Y = sourceRho * std::sin(sourcePhi);
  result.Position.Z = fSourceDistance;

  result.Direction.X =  0.0;
  result.Direction.Y =  0.0;
  result.Direction.Z = -1.0;

  const double cosmicZenithAngle  = std::acos(std::cbrt(fNormalDistRand()));
  const double cosmicAzimuthAngle = twopi * fNormalDistRand();

  result.Position.RotateY(cosmicZenithAngle);
  result.Position.RotateZ(cosmicAzimuthAngle);

  result.Direction.RotateY(cosmicZenithAngle);
  result.Direction.RotateZ(cosmicAzimuthAngle);

  return result;
}

#include "Math/SMatrix.h"
#include "Coordinate.hh"

Tron::Coordinate::Coordinate() {
  ResetOrigin();
  ResetAxes();
}

Tron::Coordinate::~Coordinate() {
}

TVector3 Tron::Coordinate::GetLocalPosition(const TVector3& globalPosition) const {
  const Double_t arrayA[] = {
    fXaxis.X(), fYaxis.X(), fZaxis.X(),
    fXaxis.Y(), fYaxis.Y(), fZaxis.Y(),
    fXaxis.Z(), fYaxis.Z(), fZaxis.Z(),
  };
  
  const ROOT::Math::SMatrix<Double_t, 3, 3> matrixA(arrayA, arrayA + sizeof(arrayA) / sizeof(*arrayA));

  const TVector3 vectorb = globalPosition - fOrigin;

  const Double_t arrayB[] = {
    vectorb.X(),
    vectorb.Y(),
    vectorb.Z(),
  };

  const ROOT::Math::SMatrix<Double_t, 3, 1> matrixB(arrayB, arrayB + sizeof(arrayB) / sizeof(*arrayB));
  
  int isFailed = false;
  const ROOT::Math::SMatrix<Double_t, 3, 1> vectorx = matrixA.Inverse(isFailed) * matrixB;
  
  return TVector3(vectorx.Array()[0], vectorx.Array()[1], vectorx.Array()[2]);
}

TVector3 Tron::Coordinate::GetGlobalPosition(const TVector3& localPosition) const {
  return fOrigin
    + localPosition.X() * fXaxis
    + localPosition.Y() * fYaxis
    + localPosition.Z() * fZaxis;
}

const TVector3& Tron::Coordinate::GetOrigin() const {
  return fOrigin;
}

void Tron::Coordinate::SetOrigin(const TVector3& origin) {
  fOrigin = origin;
}

void Tron::Coordinate::ShiftOrigin(const TVector3 dx) {
  fOrigin += dx;
}

void Tron::Coordinate::ResetOrigin() {
  fOrigin.SetXYZ(0.0, 0.0, 0.0);
}

const TVector3& Tron::Coordinate::GetXaxis() const {
  return fXaxis;
}

const TVector3& Tron::Coordinate::GetYaxis() const {
  return fYaxis;
}

const TVector3& Tron::Coordinate::GetZaxis() const {
  return fZaxis;
}

void Tron::Coordinate::SetXaxis(const TVector3& xaxis) {
  fXaxis = xaxis;
}

void Tron::Coordinate::SetYaxis(const TVector3& yaxis) {
  fYaxis = yaxis;
}

void Tron::Coordinate::SetZaxis(const TVector3& zaxis) {
  fZaxis = zaxis;
}

void Tron::Coordinate::RotateAxes(Double_t theta, const TVector3& axis) {
  fXaxis.Rotate(theta, axis);
  fYaxis.Rotate(theta, axis);
  fZaxis.Rotate(theta, axis);
}

void Tron::Coordinate::ResetAxes() {
  fXaxis.SetXYZ(1.0, 0.0, 0.0);
  fYaxis.SetXYZ(0.0, 1.0, 0.0);
  fZaxis.SetXYZ(0.0, 0.0, 1.0);
}

#include "TGraph.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "TMarker.h"
#include "Ellipse.hh"

Tron::Ellipse::Ellipse()
  : fSigmaX(1), fSigmaY(1), fCorrelation(0) {
}

void Tron::Ellipse::SetMuX(Double_t mu) {
  fMuX = mu;
}

void Tron::Ellipse::SetMuY(Double_t mu) {
  fMuY = mu;
}

void Tron::Ellipse::SetSigmaX(Double_t sigma) {
  fSigmaX = sigma;
}

void Tron::Ellipse::SetSigmaY(Double_t sigma) {
  fSigmaY = sigma;
}

void Tron::Ellipse::SetCorrelation(Double_t correlation) {
  fCorrelation = correlation;
}

TGraph* Tron::Ellipse::CreateGraph(const std::string& name) {
  TGraph* graph = new TGraph(fNpx);
  graph->SetName(name.data());
  graph->SetLineColor(kRed);
  graph->SetMarkerColor(kRed);
  graph->SetMarkerStyle(kDot);
  graph->SetLineWidth(2);
  graph->SetLineStyle(kSolid);

  const Double_t a   = GetSemiAxisX(fSigmaX, fSigmaY, fCorrelation);
  const Double_t b   = GetSemiAxisY(fSigmaX, fSigmaY, fCorrelation);
  const Double_t phi = GetRotationAngle(fSigmaX, fSigmaY, fCorrelation);
  for (Int_t i = 0; i <= fNpx; ++i) {
    const Double_t theta = TMath::TwoPi() * i / fNpx;
    const Double_t x     = GetX(theta, a, b, phi) + fMuX;
    const Double_t y     = GetY(theta, a, b, phi) + fMuY;
    graph->SetPoint(i, x, y);
  }

  return graph;
}

TF1* Tron::Ellipse::CreateFuncX(const std::string& name) {
  TF1* func = new TF1(name.data(),
                      "[0]"
                      " + [1] * TMath::Cos(x) * TMath::Cos([3])"
                      " + [2] * TMath::Sin(x) * TMath::Sin([3])",
                      0.0, TMath::TwoPi());
  func->SetParNames("MuX",
                    "SemiAxisX",
                    "SemiAxisY",
                    "RotationAngle");
  func->SetParameters(fMuX,
                      GetSemiAxisX    (fSigmaX, fSigmaY, fCorrelation),
                      GetSemiAxisY    (fSigmaX, fSigmaY, fCorrelation),
                      GetRotationAngle(fSigmaX, fSigmaY, fCorrelation));
  return func;
}

TF1* Tron::Ellipse::CreateFuncY(const std::string& name) {
  TF1* func = new TF1(name.data(),
                      "[0]"
                      " + [1] * TMath::Cos(x) * TMath::Sin([3])"
                      " + [2] * TMath::Sin(x) * TMath::Cos([3])",
                      0.0, TMath::TwoPi());
  func->SetParNames("MuY",
                    "SemiAxisX",
                    "SemiAxisY",
                    "RotationAngle");
  func->SetParameters(fMuY,
                      GetSemiAxisX    (fSigmaX, fSigmaY, fCorrelation),
                      GetSemiAxisY    (fSigmaX, fSigmaY, fCorrelation),
                      GetRotationAngle(fSigmaX, fSigmaY, fCorrelation));
  return func;
}

TF2* Tron::Ellipse::CreateFunc(const std::string& name) {
  const Double_t xmin = fMuX - fSigmaX * 3.0;
  const Double_t xmax = fMuX + fSigmaX * 3.0;
  const Double_t ymin = fMuY - fSigmaY * 3.0;
  const Double_t ymax = fMuY + fSigmaY * 3.0;
  TF2* func = new TF2(name.data(),
                      "((x-[0])*(x-[0])/[2]/[2] - 2*[4]*(x-[0])*(y-[1])/[2]/[3] + (y-[1])*(y-[1])/[3]/[3]) / (1-[4]*[4])",
                      xmin, xmax, ymin, ymax);
  func->SetParNames("MuX",
                    "MuY",
                    "SigmaX",
                    "SigmaY",
                    "Correlation");
  func->SetParameters(fMuX, fMuY, fSigmaX, fSigmaY, fCorrelation);
  return func;
}

TMarker* Tron::Ellipse::CreateCenterMarker() {
  TMarker* marker = new TMarker(fMuX, fMuY, kFullCircle);
  marker->SetMarkerColor(kRed);
  return marker;
}

Double_t Tron::Ellipse::GetSemiAxisX(Double_t sigmaX, Double_t sigmaY, Double_t correlation) {
  const Double_t sigmaX2 = sigmaX * sigmaX;
  const Double_t sigmaY2 = sigmaY * sigmaY;
  const Double_t rho2    = correlation * correlation;
  return TMath::Sqrt((sigmaX2 * sigmaX2 - rho2 * sigmaY2 * sigmaY2) / (sigmaX2 - rho2 * sigmaY2));
}

Double_t Tron::Ellipse::GetSemiAxisY(Double_t sigmaX, Double_t sigmaY, Double_t correlation) {
  const Double_t sigmaX2 = sigmaX * sigmaX;
  const Double_t sigmaY2 = sigmaY * sigmaY;
  const Double_t rho2    = correlation * correlation;
  return TMath::Sqrt((1.0 - rho2) * sigmaX2 * sigmaY2 / (sigmaX2 - rho2 * sigmaY2));
}

Double_t Tron::Ellipse::GetRotationAngle(Double_t sigmaX, Double_t sigmaY, Double_t correlation) {
  return TMath::ATan2(correlation * sigmaY, sigmaX);
}

Double_t Tron::Ellipse::GetX(Double_t theta, Double_t semiAxisX, Double_t semiAxisY, Double_t rotationAngle) {
  return
    semiAxisX * TMath::Cos(theta) * TMath::Cos(rotationAngle) +
    semiAxisY * TMath::Sin(theta) * TMath::Sin(rotationAngle);
}

Double_t Tron::Ellipse::GetY(Double_t theta, Double_t semiAxisX, Double_t semiAxisY, Double_t rotationAngle) {
  return
    semiAxisX * TMath::Cos(theta) * TMath::Sin(rotationAngle) +
    semiAxisY * TMath::Sin(theta) * TMath::Cos(rotationAngle);
}

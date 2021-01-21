#include "TMath.h"
#include "HoughTransformator.hh"

Double_t Tron::HoughLine::Y(Double_t* x, Double_t* p) {
  const Double_t
    _theta = p[0],
    _r     = p[1],
    _x     = x[0],
    _y     = (_r - _x * TMath::Cos(_theta)) / TMath::Sin(_theta);
  return _y;
}

Double_t Tron::HoughLine::R(Double_t* x, Double_t* p) {
  const Double_t
    _x     = p[0],
    _y     = p[1],
    _theta = x[0],
    _r     = _x * TMath::Cos(_theta) + _y * TMath::Sin(_theta);
  return _r;
}

Double_t Tron::HoughLine::Distance(Double_t theta, Double_t r, Double_t x, Double_t y) {
  return TMath::Abs(x * TMath::Cos(theta) + y * TMath::Sin(theta) - r);
}

Tron::HoughTransformator::HoughTransformator(const Char_t* name,
                                             Int_t nbinst, Double_t tlow, Double_t tup,
                                             Int_t nbinsr, Double_t rlow, Double_t rup)
  : fTheta(0.0), fR(0.0){
  //--- Create hist
  fTransResult = new TH2F(Form("%s_tr", name), "Transform Result;#theta [rad];r",
                          nbinst, tlow, tup, nbinsr, rlow, rup);

  //--- Create func
  fChiSquare   = new TF2 (Form("%s_cs", name), this, &HoughTransformator::ChiSquare,
                          tlow, tup, rlow, rup, 0, "HoughTransformator", "ChiSquare");
}

Tron::HoughTransformator::~HoughTransformator() {
  if (fTransResult) { 
    delete fTransResult;
    fTransResult = nullptr;
  }

  if (fChiSquare) {
    delete fChiSquare;
    fChiSquare = nullptr;
  }
}

Int_t Tron::HoughTransformator::GetN() const {
  return fW.size() - std::count(fW.begin(), fW.end(), 0.0);
}

void Tron::HoughTransformator::SetPoints(Int_t n, const Double_t* x, const Double_t* y, const Double_t* w) {
  fX.assign(x, x + n);
  fY.assign(y, y + n);
  fW.assign(w, w + n);
}

void Tron::HoughTransformator::SetPoints(Int_t n, const Double_t* x, const Double_t* y, Double_t w) {
  fX.assign(x, x + n);
  fY.assign(y, y + n);
  fW.assign(n, w);
}

const TH2* Tron::HoughTransformator::Transform(Int_t smoothTimes) {
  fTransResult->Reset();

  //--- Vote
  for (Int_t i = 0, n = fW.size(); i < n; ++i) {
    if (fW[i]) {
      Vote(fX[i], fY[i], fW[i]);
    }
  }

  //--- Smoothing
  if (smoothTimes > 0) {
    fTransResult->Smooth(smoothTimes);
  }

  //--- Get (theta, r)
  Int_t binX, binY, binZ;
  fTransResult->GetMaximumBin(binX, binY, binZ);
  fTheta = fTransResult->GetXaxis()->GetBinCenter(binX);
  fR     = fTransResult->GetYaxis()->GetBinCenter(binY);

  return fTransResult;
}

Int_t Tron::HoughTransformator::RemoveFarPoints(Double_t theta, Double_t r, Double_t threshold) {
  for (Int_t i = 0, n = fW.size(); i < n; ++i) {
    // Check Already Removed
    if (!fW[i]) {
      continue;
    }
    // Check Distance
    const Double_t distance = HoughLine::Distance(theta, r, fX[i], fY[i]);
    if (distance <= threshold) {
      continue;
    }
    // Remove Point
    fW[i] = 0.0;
  }
  return GetN();
}

const TF2* Tron::HoughTransformator::MinimizeChiSquare() {
  fChiSquare->GetMinimumXY(fTheta, fR);
  return fChiSquare;
}

TF1* Tron::HoughTransformator::CreateHoughLine(const Char_t* name, Double_t tmin, Double_t tmax) {
  TF1* func = new TF1(name, HoughLine::R, tmin, tmax, 2);
  func->SetParNames("x", "y");
  func->SetParameters(0.0, 0.0);
  return func;
}

TF1* Tron::HoughTransformator::CreateTrack(const Char_t* name, Double_t xmin, Double_t xmax) {
  TF1* func = new TF1(name, HoughLine::Y, xmin, xmax, 2);
  func->SetParNames("#theta", "r");
  func->SetParameters(fTheta, fR);
  return func;
}

void Tron::HoughTransformator::Vote(Double_t x, Double_t y, Double_t w) {
  const Int_t weight = 1.0 / w;
  const Int_t nbinst = fTransResult->GetXaxis()->GetNbins();
  Double_t p[2] = { x, y };

  //--- Get (theta, r) @ lower edge
  Double_t tlow    = fTransResult->GetXaxis()->GetBinLowEdge(0);
  Double_t rlow    = HoughLine::R(&tlow, p);
  Double_t rbinlow = fTransResult->GetYaxis()->FindBin(rlow);
  
  for (Int_t tbin = 0; tbin < nbinst; ++tbin) {
    //--- Get (theta, r) @ upper edge
    Double_t tup    = fTransResult->GetXaxis()->GetBinUpEdge(tbin);
    Double_t rup    = HoughLine::R(&tup, p);
    Double_t rbinup = fTransResult->GetYaxis()->FindBin(rup);

    //--- Fill points from (tlow, rlow) to (tup, rup) as linear
    const Double_t theta = (tlow + tup) * 0.5;
    for (Int_t rbin = rbinlow; rbin <= rbinup; rbin++) {
      const Double_t r = fTransResult->GetYaxis()->GetBinCenter(rbin);
      fTransResult->Fill(theta, r, weight);
    }

    //--- Update Parameters
    tlow    = tup;
    rlow    = rup;
    rbinlow = rbinup;
  }
}

Double_t Tron::HoughTransformator::ChiSquare(Double_t* x, Double_t*) {
  return EvalChiSquare(x[0], x[1]);
}

Double_t Tron::HoughTransformator::EvalChiSquare(Double_t theta, Double_t r) {
  Double_t sum = 0;
  for (Int_t i = 0, n = fW.size(); i < n; ++i) {
    if (fW[i]) {
      const Double_t distance = HoughLine::Distance(theta, r, fX[i], fY[i]) / fW[i];
      sum += distance * distance;
    }
  }
  return sum;
}

#include "TSF1.hh"
#include "Math.hh"

//==================================================
Tron::TSF1::TSF1(const Char_t* name, const TH1* hist, Double_t xmin, Double_t xmax)
  : TF1(name, this, &TSF1::Function, xmin, xmax, kNpar, "TSF1", "Function") {
  fSpline = new TSpline3(hist);

  //--- 
  SetParNames("amp", "mu", "calib", "base");
  SetParameters(1, 0, 1, 0);
}

Tron::TSF1::TSF1(const Char_t* name, const TSpline3* spline, Double_t xmin, Double_t xmax)
  : TF1(name, this, &TSF1::Function, xmin, xmax, kNpar, "TSF1", "Function") {
  fSpline = new TSpline3(*spline);

  //---
  SetParNames("amp", "mu", "calib", "base");
  SetParameters(1, 0, 1, 0);
}

Tron::TSF1::TSF1(const TSF1& sf1)
  : TF1() {
  fSpline = new TSpline3(*(sf1.GetSpline()));

  ((TF1&)sf1).Copy(*this);
}

Tron::TSF1::~TSF1() {
  delete fSpline; fSpline = nullptr;
}


//==================================================
const Int_t Tron::TSF1::kNpar = 4;

Double_t Tron::TSF1::Function(Double_t *x, Double_t *p) {

  //--- fitting parameters
  const Double_t* amp   = p;
  const Double_t* mu    = p + 1;
  const Double_t* calib = p + 2;
  const Double_t* base  = p + 3;

  //--- other parameters
  const Double_t xx = (*x - *mu) / *calib;
  const Double_t xxx = Math::Bind(xx, fSpline->GetXmin(), fSpline->GetXmax());

  return *amp * fSpline->Eval(xxx) / *calib + *base;
}

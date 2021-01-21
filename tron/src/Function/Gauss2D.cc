#include "TMath.h"
#include "Function/Gauss2D.hh"

const Int_t Tron::Gauss2D::kNpar = 6;

Int_t Tron::Gauss2D::Npar() const {
  return kNpar;
}

Double_t Tron::Gauss2D::Func(Double_t* x, Double_t* p) const {
  return Tron::Gauss2D::Eval(x[0], x[1], p[0], p[1], p[2], p[3], p[4], p[5]);
}

Double_t Tron::Gauss2D::Eval(Double_t x, Double_t y, Double_t height, Double_t muX, Double_t muY, Double_t sigmaX, Double_t sigmaY, Double_t correlation, Double_t, Double_t, Double_t, Double_t, Double_t) const {
  //--- Convert sigma
  sigmaX = TMath::Abs(sigmaX);
  sigmaY = TMath::Abs(sigmaY);

  //--- Convert x, y
  const Double_t X = (x - muX) / sigmaX;
  const Double_t Y = (y - muY) / sigmaY;

  return height * TMath::Exp(- (X * X - 2.0 * correlation * X * Y + Y * Y) / 2. / (1. - correlation * correlation));
}

TF2* Tron::Gauss2D::Create(const Char_t* name) const {
  TF2* func = new TF2(name, this, &Gauss2D::Func, 0.0, 1.0, 0.0, 1.0, kNpar);
  func->SetParNames("height", "#mu_{x}", "#mu_{y}", "#sigma_{x}", "#sigma_{y}", "#rho");
  return func;
}

TF2* Tron::Gauss2D::Fit(TH2* hist, Double_t xmin, Double_t xmax, Param_t* param) const {
  if (!hist) {
    return nullptr;
  }

  //--- Set range
  if (xmin < xmax) {
    hist->GetXaxis()->SetRangeUser(xmin, xmax);
  } else {
    xmin = hist->GetXaxis()->GetXmin();
    xmax = hist->GetXaxis()->GetXmax();
  }

  Double_t ymin = 0, ymax = 0;
  if (ymin < xmax) {
    hist->GetYaxis()->SetRangeUser(ymin, ymax);
  } else {
    ymin = hist->GetYaxis()->GetXmin();
    ymax = hist->GetYaxis()->GetXmax();
  }

  //--- Get fit function
  TF2* func = Create("Tron::Gauss2D");

  //--- Set initial parameters
  Param_t initParam;
  if (param) {
    initParam = *param;
  } else {
    Int_t binx, biny, binz;
    hist->GetMaximumBin(binx, biny, binz);
  
    initParam.Height      = hist->GetMaximum();
    initParam.MuX         = hist->GetXaxis()->GetBinCenter(binx);
    initParam.MuY         = hist->GetYaxis()->GetBinCenter(biny);
    initParam.SigmaX      = hist->GetRMS(1);
    initParam.SigmaY      = hist->GetRMS(2);
    initParam.Correlation = 0.0;
  }

  //--- Set limits
  func->SetParLimits (0,
                      initParam.Height *  0.50,
                      initParam.Height *  2.00);
  func->SetParLimits (1,
                      xmin,
                      xmax);
  func->SetParLimits (2,
                      ymin,
                      ymax);
  func->SetParLimits (3,
                      initParam.SigmaX *  0.05,
                      initParam.SigmaX * 10.00);
  func->SetParLimits (4,
                      initParam.SigmaY *  0.05,
                      initParam.SigmaY * 10.00);
  func->SetParLimits (5,
                      -0.00,
                      +0.00);
  
  //--- Set parameters
  func->SetParameters(initParam.Height,
                      initParam.MuX,
                      initParam.MuY,
                      initParam.SigmaX,
                      initParam.SigmaY,
                      initParam.Correlation);

  //--- Fitting
  hist->Fit(func, "Q", "goff");
  
  return func;
}

TF2* Tron::Gauss2D::Fit(TGraph2D* /*graph*/, Double_t /*xmin*/, Double_t /*xmax*/, Param_t* /*param*/) const {
  return nullptr;
}

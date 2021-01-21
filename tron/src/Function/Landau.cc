#include "TMath.h"
#include "Function/Landau.hh"
#include "PeakAnalyzer.hh"

const Int_t Tron::Landau::kNpar = 3;

Int_t Tron::Landau::Npar() const {
  return kNpar;
}

Double_t Tron::Landau::Func(Double_t* x, Double_t* p) const {
  return Tron::Landau::Eval(x[0], p[0], p[1], p[2]);
}

Double_t Tron::Landau::Eval(Double_t x, Double_t height, Double_t location, Double_t sigma, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t) const {
  //--- Convert sigma
  sigma = TMath::Abs(sigma);
  
  return height * TMath::Landau(x, location, sigma, kFALSE);
}

TF1* Tron::Landau::Create(const Char_t* name) const {
  TF1* func = new TF1(name, this, &Landau::Func, 0.0, 1.0, kNpar);
  func->SetParNames("height", "#mu", "#sigma");
  return func;
}

TF1* Tron::Landau::Fit(TH1* hist, Double_t xmin, Double_t xmax, Param_t* param) const {
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

  //--- Get fit function
  TF1* func = Create("Tron::Landau");

  //--- Set initial parameters
  Param_t initParam;
  if (param) {
    initParam = *param;
  } else {
    const Double_t invLmpv  = 1.0 / 4.01865;
    
    const Int_t    maxbin = hist->GetMaximumBin();
    const Double_t maxx   = hist->GetBinCenter(maxbin);
    const Double_t maxy   = hist->GetBinContent(maxbin);
    const Double_t fwhm   = PeakAnalyzer::GetFWHMFast(hist, xmin, xmax, 1);
  
    initParam.Height   = 1.0;
    initParam.Location = maxx;
    initParam.Sigma    = fwhm * invLmpv;
  
    func->SetParameters(initParam.Height, initParam.Location, initParam.Sigma);
    initParam.Height   = maxy / func->GetMaximum();
  }

  //--- Set limits
  func->SetParLimits (0,
                      initParam.Height *  0.50,
                      initParam.Height *  2.00);
  func->SetParLimits (1,
                      xmin,
                      xmax);
  func->SetParLimits (2,
                      initParam.Sigma  *  0.05,
                      initParam.Sigma  * 10.00);
  
  //--- Set parameters
  func->SetParameters(initParam.Height,
                      initParam.Location,
                      initParam.Sigma);

  //--- Fitting
  hist->Fit(func, "Q", "goff", xmin, xmax);
  
  return func;
}

TF1* Tron::Landau::Fit(TGraph* /*graph*/, Double_t /*xmin*/, Double_t /*xmax*/, Param_t* /*param*/) const {
  return nullptr;
}

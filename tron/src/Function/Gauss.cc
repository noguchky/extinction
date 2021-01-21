#include "TMath.h"
#include "Function/Gauss.hh"
#include "PeakAnalyzer.hh"

const Int_t Tron::Gauss::kNpar = 3;

Int_t Tron::Gauss::Npar() const {
  return kNpar;
}

Double_t Tron::Gauss::Func(Double_t* x, Double_t* p) const {
  return Tron::Gauss::Eval(x[0], p[0], p[1], p[2]);
}

Double_t Tron::Gauss::Eval(Double_t x, Double_t height, Double_t mu, Double_t sigma, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t) const {
  //--- Convert sigma
  sigma = TMath::Abs(sigma);
  
  return height * TMath::Gaus(x, mu, sigma, kFALSE);
}

TF1* Tron::Gauss::Create(const Char_t* name) const {
  TF1* func = new TF1(name, this, &Gauss::Func, 0.0, 1.0, kNpar);
  func->SetParNames("height", "#mu", "#sigma");
  return func;
}

TF1* Tron::Gauss::Fit(TH1* hist, Double_t xmin, Double_t xmax, Param_t* param) const {
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
  TF1* func = Create("Tron::Gauss");

  //--- Set initial parameters
  Param_t initParam;
  if (param) {
    initParam = *param;
  } else {
    const Double_t fwhm2rms = 0.5 / TMath::Sqrt(2.0 * TMath::Log(2.0));

    const Int_t    maxbin = hist->GetMaximumBin();
    const Double_t maxx   = hist->GetBinCenter(maxbin);
    const Double_t maxy   = hist->GetBinContent(maxbin);
    const Double_t fwhm   = PeakAnalyzer::GetFWHMFast(hist, xmin, xmax, 1);

    initParam.Height = maxy;
    initParam.Mu     = maxx;
    initParam.Sigma  = fwhm * fwhm2rms;
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
                      initParam.Mu,
                      initParam.Sigma);

  //--- Fitting
  hist->Fit(func, "Q", "goff", xmin, xmax);
  
  return func;
}

TF1* Tron::Gauss::Fit(TGraph* /*graph*/, Double_t /*xmin*/, Double_t /*xmax*/, Param_t* /*param*/) const {
  return nullptr;
}

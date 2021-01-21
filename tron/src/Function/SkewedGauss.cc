#include "TMath.h"
#include "Function/SkewedGauss.hh"
#include "PeakAnalyzer.hh"

const Int_t Tron::SkewedGauss::kNpar = 4;

Int_t Tron::SkewedGauss::Npar() const {
  return kNpar;
}

Double_t Tron::SkewedGauss::Func(Double_t* x, Double_t* p) const {
  return Tron::SkewedGauss::Eval(x[0], p[0], p[1], p[2], p[3]);
}

Double_t Tron::SkewedGauss::Eval(Double_t x, Double_t height, Double_t mu, Double_t sigma, Double_t slope, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t) const {
  //--- Convert sigma
  sigma = TMath::Abs(sigma);

  const Double_t invSqrt2 = 1.0 / TMath::Sqrt2();
  const Double_t xx = (x - mu) / slope;
  const Double_t ss = sigma / slope;
  
  return
    height
    * (0.5 / slope)
    * TMath::Exp(0.5 * ss * ss + xx)
    * (1.0 + TMath::Erf(-invSqrt2 * (xx - ss * ss) / ss));
}

TF1* Tron::SkewedGauss::Create(const Char_t* name) const {
  TF1* func = new TF1(name, this, &SkewedGauss::Func, 0.0, 1.0, kNpar);
  func->SetParNames("height", "#mu", "#sigma", "slope");
  return func;
}

TF1* Tron::SkewedGauss::Fit(TH1* hist, Double_t xmin, Double_t xmax, Param_t* param) const {
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
  TF1* func = Create("Tron::SkewedGauss");

  //--- Set initial parameters
  Param_t initParam;
  if (param) {
    initParam = *param;
  } else {
    const Double_t fwhm2rms = 0.5 / TMath::Sqrt(2.0 * TMath::Log(2.0));
    
    const Double_t hwhmRatio = PeakAnalyzer::GetHWHMFast(hist, xmax, xmin) / PeakAnalyzer::GetHWHMFast(hist, xmin, xmax);

    initParam.Height = hist->GetMaximum();
    initParam.Mu     = hist->GetBinCenter(hist->GetMaximumBin());
    initParam.Sigma  = PeakAnalyzer::GetFWHMFast(hist, xmin, xmax) * fwhm2rms;
    initParam.Slope  = (-2.84739 + 3.40816 * hwhmRatio) * initParam.Sigma;

    func->SetParameters(initParam.Height,
                        initParam.Mu,
                        initParam.Sigma,
                        initParam.Slope);
    initParam.Height = func->GetParameter(0) * hist->GetMaximum() / func->GetMaximum();
    initParam.Mu     = func->GetParameter(1) * hist->GetBinCenter(hist->GetMaximumBin()) / func->GetMaximumX();
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
  func->SetParLimits (3,
                      initParam.Slope  *  0.05,
                      initParam.Slope  * 10.00);
  
  //--- Set parameters
  func->SetParameters(initParam.Height,
                      initParam.Mu,
                      initParam.Sigma,
                      initParam.Slope);

  //--- Fitting
  hist->Fit(func, "Q", "goff", xmin, xmax);
  
  return func;
}

TF1* Tron::SkewedGauss::Fit(TGraph* /*graph*/, Double_t /*xmin*/, Double_t /*xmax*/, Param_t* /*param*/) const {
  return nullptr;
}

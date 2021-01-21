#include "TMath.h"
#include "Function/LandauGauss.hh"
#include "Integrator.hh"
#include "PeakAnalyzer.hh"

const Int_t Tron::LandauGauss::kNpar = 4;

Int_t Tron::LandauGauss::Npar() const {
  return kNpar;
}

Double_t Tron::LandauGauss::Func(Double_t* x, Double_t* p) const {
  return Tron::LandauGauss::Eval(x[0], p[0], p[1], p[2], p[3]);
}

Double_t Tron::LandauGauss::Eval(Double_t x, Double_t height, Double_t location, Double_t sigmaLandau, Double_t sigmaGauss, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t) const {
  //--- Convert sigma
  sigmaLandau = TMath::Abs(sigmaLandau);
  sigmaGauss  = TMath::Abs(sigmaGauss );

  //--- Convert height
  height = 0.3989422804014 * height / sigmaLandau / sigmaGauss;

  //--- Convert location
  Double_t locationLandau = location + 0.22278298 * sigmaLandau;
  Double_t locationGauss  = x;

  //--- Integral range
  struct {
    Double_t lower;
    Double_t upper;
  } rangeLandau, rangeGauss, *rangeInner, *rangeOuter;

  // rangeLandau.lower = locationLandau -  2.84464714 * sigmaLandau;
  // rangeLandau.upper = locationLandau + 25.587555   * sigmaLandau;
  // rangeGauss .lower = locationGauss  -  3.03485425 * sigmaGauss;
  // rangeGauss .upper = locationGauss  +  3.03485425 * sigmaGauss;
  rangeLandau.lower = locationLandau - 3.0 * sigmaLandau;
  rangeLandau.upper = locationLandau + 6.0 * sigmaLandau;
  rangeGauss .lower = locationGauss  - 3.0 * sigmaGauss;
  rangeGauss .upper = locationGauss  + 3.0 * sigmaGauss;

  rangeInner = sigmaLandau < sigmaGauss ? &rangeLandau : &rangeGauss;
  rangeOuter = sigmaLandau < sigmaGauss ? &rangeGauss  : &rangeLandau;

  //--- Integral width
  Double_t dxInner = 0.16 * TMath::Min(sigmaLandau, sigmaGauss);
  Double_t dxOuter = 0.16 * TMath::Max(sigmaLandau, sigmaGauss);

  //--- Integral function
  const auto func = [=](Double_t _xi) {
    return TMath::Landau(_xi, locationLandau, sigmaLandau, kFALSE) * TMath::Gaus(_xi, locationGauss, sigmaGauss, kFALSE);
  };

  //--- Integrate
  Double_t xi = 0.0, sum = 0.0;
  static const auto& integrator = Integrator::SimpsonMethod::Get();
  if (rangeOuter->lower < rangeInner->lower) {
    sum += integrator.Integrate(func, xi, rangeInner->lower, rangeOuter->lower, dxOuter);
  }
  {
    sum += integrator.Integrate(func, xi, rangeInner->lower, rangeInner->upper, dxInner);
  }
  if (xi < rangeOuter->upper) {
    sum += integrator.Integrate(func, xi, xi,                rangeOuter->upper, dxOuter);
  }

  return height * sum;
}

TF1* Tron::LandauGauss::Create(const Char_t* name) const {
  TF1* func = new TF1(name, this, &LandauGauss::Func, 0.0, 1.0, kNpar);
  func->SetParNames("height", "#mu", "#sigma_{Landau}", "#sigma_{Gauss}");
  return func;
}

TF1* Tron::LandauGauss::Fit(TH1* hist, Double_t xmin, Double_t xmax, Param_t* param) const {
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
  TF1* func = Create("Tron::LandauGauss");

  //--- Set initial parameters
  Param_t initParam;
  if (param) {
    initParam = *param;
  } else {
    // const Double_t invLmax  = 1.0 / 1.80655633634815849e-01;
    const Double_t invLmpv  = 1.0 / 4.01865;
    const Double_t fwhm2rms = 0.5 / TMath::Sqrt(2.0 * TMath::Log(2.0));

    const Int_t    maxbin = hist->GetMaximumBin();
    const Double_t maxx   = hist->GetBinCenter(maxbin);
    const Double_t maxy   = hist->GetBinContent(maxbin);
    const Double_t fwhm   = PeakAnalyzer::GetFWHMFast(hist, xmin, xmax, 1);

    initParam.Height      = 1.0;
    initParam.Location    = maxx;
    initParam.SigmaLandau = fwhm * invLmpv;
    initParam.SigmaGauss  = fwhm * fwhm2rms;

    func->SetParameters(initParam.Height, initParam.Location, initParam.SigmaLandau, initParam.SigmaGauss);
    initParam.Height      = maxy / func->GetMaximum();
  }

  //--- Set limits
  func->SetParLimits (0,
                      initParam.Height      *  0.50,
                      initParam.Height      *  2.00);
  func->SetParLimits (1,
                      xmin,
                      xmax);
  func->SetParLimits (2,
                      initParam.SigmaLandau *  0.05,
                      initParam.SigmaLandau * 10.00);
  func->SetParLimits (3,
                      initParam.SigmaGauss  *  0.05,
                      initParam.SigmaGauss  * 10.00);

  //--- Set parameters
  func->SetParameters(initParam.Height,
                      initParam.Location,
                      initParam.SigmaLandau,
                      initParam.SigmaGauss);

  //--- Fitting
  hist->Fit(func, "Q", "goff", xmin, xmax);

  return func;
}

TF1* Tron::LandauGauss::Fit(TGraph* /*graph*/, Double_t /*xmin*/, Double_t /*xmax*/, Param_t* /*param*/) const {
  return nullptr;
}

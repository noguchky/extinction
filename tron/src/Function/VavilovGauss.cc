#include "TMath.h"
#include "Function/VavilovGauss.hh"
#include "Integrator.hh"

const Int_t Tron::VavilovGauss::kNpar = 6;

Int_t Tron::VavilovGauss::Npar() const {
  return kNpar;
}

Double_t Tron::VavilovGauss::Func(Double_t* x, Double_t* p) const {
  return Tron::VavilovGauss::Eval(x[0], p[0], p[1], p[2], p[3], p[4], p[5]);
}

Double_t Tron::VavilovGauss::Eval(Double_t x, Double_t integral, Double_t location, Double_t kappa, Double_t beta2, Double_t sigmaVavilov, Double_t sigmaGauss, Double_t, Double_t, Double_t, Double_t, Double_t) const {
  // Convert sigma
  sigmaVavilov = TMath::Abs(sigmaVavilov);
  sigmaGauss   = TMath::Abs(sigmaGauss  );

  // Convert integral
  // integral = integral;

  // Convert location
  Double_t locationVavilov = location;
  Double_t locationGauss   = x;

  // Integral range
  struct {
    Double_t lower;
    Double_t upper;
  } rangeLandau, rangeGauss, *rangeInner, *rangeOuter;

  rangeLandau.lower = locationVavilov - 3.0 * sigmaVavilov;
  rangeLandau.upper = locationVavilov + 5.0 * sigmaVavilov;
  rangeGauss .lower = locationGauss   - 4.0 * sigmaGauss;
  rangeGauss .upper = locationGauss   + 4.0 * sigmaGauss;

  rangeInner = sigmaVavilov < sigmaGauss ? &rangeLandau : &rangeGauss;
  rangeOuter = sigmaVavilov < sigmaGauss ? &rangeGauss  : &rangeLandau;

  // Integral width
  Double_t dxInner = 0.16 * TMath::Min(sigmaVavilov, sigmaGauss);
  Double_t dxOuter = 0.16 * TMath::Max(sigmaVavilov, sigmaGauss);

  // Integral function
  const auto func = [=](Double_t _xi) {
    Double_t _xx = (_xi - locationVavilov) / sigmaVavilov;
    return TMath::Vavilov(_xx, kappa, beta2) * TMath::Gaus(_xi, locationGauss, sigmaGauss, kFALSE);
  };

  // Integrate
  Double_t xi = 0.0, sumOuter = 0.0, sumInner = 0.0;
  static const auto& integrator = Integrator::SimpsonMethod::Get();
  if (rangeOuter->lower < rangeInner->lower) {
    sumOuter += integrator.Integrate(func, xi, rangeInner->lower, rangeOuter->lower, dxOuter);
  }
  {
    sumInner += integrator.Integrate(func, xi, rangeInner->lower, rangeInner->upper, dxInner);
  }
  if (xi < rangeOuter->upper) {
    sumOuter += integrator.Integrate(func, xi, xi,                rangeOuter->upper, dxOuter);
  }

  return integral * (sumOuter * dxOuter + sumInner * dxInner);
}

TF1* Tron::VavilovGauss::Create(const Char_t* name) const {
  TF1* func = new TF1(name, this, &VavilovGauss::Func, 0.0, 1.0, kNpar);
  func->SetParNames("integral", "location", "#kappa", "#beta^{2}", "#sigma_{Vavilov}", "#sigma_{Gauss}");
  return func;
}

TF1* Tron::VavilovGauss::Fit(TH1* /*hist*/, Double_t /*xmin*/, Double_t /*xmax*/, Param_t* /*param*/) const {
  // TODO extends
  return nullptr;
}

TF1* Tron::VavilovGauss::Fit(TGraph* /*graph*/, Double_t /*xmin*/, Double_t /*xmax*/, Param_t* /*param*/) const {
  // TODO extends
  return nullptr;
}


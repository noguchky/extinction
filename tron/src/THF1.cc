#include "TMath.h"
#include "THF1.hh"

//==================================================
Tron::THF1::THF1(const Char_t* name, const TH1* hist, Double_t xmin, Double_t xmax, Bool_t norm)
  : TF1(name, this, &THF1::Function, xmin, xmax, kNpar, "THF1", "Function") {
  if (hist) {
    const Int_t nbinsx = hist->GetNbinsX();

    //--- 
    fBinLow = hist->FindFirstBinAbove(0);
    fBinUp  = hist->FindLastBinAbove (0);
    const Double_t integralInverted = (!norm) ? 1.: 1./hist->Integral(fBinLow, fBinUp, "W");

    //--- 
    fLowEdge  = std::vector<Double_t>(nbinsx + 2, 0.);
    fCenter   = std::vector<Double_t>(nbinsx + 1, 0.);
    fContents = std::vector<Double_t>(nbinsx + 1, 0.);
    for (Int_t bin = 1; bin <= nbinsx; ++bin) {
      fLowEdge [bin] = hist->GetBinLowEdge(bin);
      fCenter  [bin] = hist->GetBinCenter (bin);
      fContents[bin] = hist->GetBinContent(bin) * integralInverted;
    }
    fLowEdge[nbinsx + 1] = hist->GetBinLowEdge(nbinsx + 1);
  }

  //--- 
  SetParNames("amp", "mu", "calib", "base", "smear");
  SetParameters(1, 0, 1, 0, 0);
}

Tron::THF1::~THF1() {
}


//==================================================
const Int_t Tron::THF1::kNpar = 5;

// smeared Hist by Gaussian with sigma
Double_t Tron::THF1::Function(Double_t *x, Double_t *p) {

  //--- Function parameters
  const Double_t* amp   = p;
  const Double_t* mu    = p + 1;
  const Double_t* calib = p + 2;
  const Double_t* base  = p + 3;
  const Double_t* smear = p + 4;

  const Double_t xx = (*x - *mu) / *calib;

  if (*smear) {
    //--- Other parameters
    const Double_t factor = TMath::Sqrt(0.5) / *smear;

    //--- Leakage Function
    auto fLeakage = [&](Int_t bin) {
                      return TMath::Erf((xx - fLowEdge[bin]) * factor);
                    };

    //--- Summation
    Double_t sum = 0.0, leakageLower, leakageUpper;
    Int_t bin = fBinLow;
    leakageUpper = fLeakage(bin);
    for (; bin <= fBinUp; ++bin) {
      leakageLower = leakageUpper;
      leakageUpper = fLeakage(bin + 1);
      sum += fContents[bin] * (leakageLower - leakageUpper);
    }
    return 0.5 * *amp * sum / *calib + *base;
  } else {
    if (xx < fLowEdge.front()) {
      return *base;
    }
    for (Int_t bin = fBinLow; bin < fBinUp; ++bin) {
      if (xx < fLowEdge[bin + 1]) {
        return *amp * fContents[bin] + *base;
      }
    }
    return *base;
  }
}

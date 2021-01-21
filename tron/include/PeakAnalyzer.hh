#ifndef Tron_PeakAnalyzer_hh
#define Tron_PeakAnalyzer_hh

#include "TH1.h"
#include "TF1.h"
#include "Types.hh"

namespace Tron {

  namespace PeakAnalyzer {

    using Func_t = std::function<Double_t(Double_t)>;
    Double_t GetFWHM(const Func_t& func, Double_t from, Double_t to, Int_t sign = 1, Double_t base = 0.0, Double_t dx = 1e-3);
    Double_t GetHWHM(const Func_t& func, Double_t from, Double_t to, Int_t sign = 1, Double_t base = 0.0, Double_t dx = 1e-3);

    Int_t    FindMaximumBin(const TH1* h1, Double_t from, Double_t to);
    Int_t    FindMinimumBin(const TH1* h1, Double_t from, Double_t to);
    Double_t Find          (const TH1* h1, Double_t from, Double_t to, Double_t y, Double_t dx = 1e-3, Double_t dy = 1e-3);
    Double_t Find          (const TF1* f1, Double_t from, Double_t to, Double_t y, Double_t dx = 1e-3, Double_t dy = 1e-3);
    Double_t FindByLSM     (const TH1* h1, Double_t from, Double_t to, Double_t y);
    Double_t GetFWHM       (const TH1* h1, Double_t from, Double_t to, Int_t sign = 1, Double_t base = 0.0);
    Double_t GetFWHMFast   (const TH1* h1, Double_t from, Double_t to, Int_t sign = 1, Double_t base = 0.0);
    Double_t GetHWHMFast   (const TH1* h1, Double_t from, Double_t to, Int_t sign = 1, Double_t base = 0.0);

  }

}

#endif

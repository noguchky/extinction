#ifndef Tron_TF1Generator_hh
#define Tron_TF1Generator_hh

#include "TF2.h"
#include "TH2.h"
#include "TGraph2D.h"
#include "Singleton.hh"

namespace Tron {

  template <typename T>
  class TF2Generator {
    TRON_SINGLETONIZE(TF2Generator);
  public:
    using Param_t = T;

    virtual Int_t    Npar() const { return 0; }
    virtual Double_t Func(Double_t* /*x*/, Double_t* /*p*/) const { return 0.0; }
    virtual Double_t Eval(Double_t /*x*/, Double_t /*y*/, Double_t /*p0*/, Double_t /*p1*/, Double_t /*p2*/, Double_t /*p3*/, Double_t /*p4*/, Double_t /*p5*/, Double_t /*p6*/, Double_t /*p7*/, Double_t /*p8*/, Double_t /*p9*/, Double_t /*p10*/) const { return 0.0; }
    virtual TF2*     Create(const Char_t* /*name*/) const { return nullptr; }
    virtual TF2*     Fit(TH2* /*hist*/, Double_t /*xmin*/ = 0.0, Double_t /*xmax*/ = 0.0, T* /*param*/ = nullptr) const { return nullptr; }
    virtual TF2*     Fit(TGraph2D* /*graph*/, Double_t /*xmin*/ = 0.0, Double_t /*xmax*/ = 0.0, T* /*param*/ = nullptr) const { return nullptr; }
    
  };
  
}

#endif

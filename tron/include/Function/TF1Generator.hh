#ifndef Tron_TF1Generator_hh
#define Tron_TF1Generator_hh

#include "TF1.h"
#include "TH1.h"
#include "TGraph.h"
#include "Singleton.hh"

namespace Tron {

  template <typename T>
  class TF1Generator {
    TRON_SINGLETONIZE(TF1Generator);
  public:
    using Param_t = T;
    
    virtual Int_t    Npar() const { return 0; }
    virtual Double_t Func(Double_t* /*x*/, Double_t* /*p*/) const { return 0.0; }
    virtual Double_t Eval(Double_t /*x*/, Double_t /*p0*/, Double_t /*p1*/, Double_t /*p2*/, Double_t /*p3*/, Double_t /*p4*/, Double_t /*p5*/, Double_t /*p6*/, Double_t /*p7*/, Double_t /*p8*/, Double_t /*p9*/, Double_t /*p10*/) const { return 0.0; }
    virtual TF1*     Create(const Char_t* /*name*/) const { return nullptr; }
    virtual TF1*     Fit(TH1* /*hist*/, Double_t /*xmin*/ = 0.0, Double_t /*xmax*/ = 0.0, Param_t* /*param*/ = nullptr) const { return nullptr; }
    virtual TF1*     Fit(TGraph* /*graph*/, Double_t /*xmin*/ = 0.0, Double_t /*xmax*/ = 0.0, Param_t* /*param*/ = nullptr) const { return nullptr; }
    
  };
  
}

#endif

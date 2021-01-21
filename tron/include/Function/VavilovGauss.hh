#ifndef Tron_VavilovGauss_hh
#define Tron_VavilovGauss_hh

#include "TF1Generator.hh"

namespace Tron {

  struct VavilovGaussParam {
    Double_t Integral;
    Double_t Location;
    Double_t Kappa;
    Double_t Beta2;
    Double_t SigmaVavilov;
    Double_t SigmaGauss;
  };

  class VavilovGauss : public TF1Generator<VavilovGaussParam> {
    TRON_SINGLETONIZE(VavilovGauss);
  private:
    static const Int_t kNpar;

  public:
    virtual Int_t    Npar() const override;
    virtual Double_t Func(Double_t* x, Double_t* p) const override;
    virtual Double_t Eval(Double_t x, Double_t integral, Double_t location, Double_t kappa, Double_t beta2, Double_t sigmaVavilov, Double_t sigmaGauss, Double_t = 0.0, Double_t = 0.0, Double_t = 0.0, Double_t = 0.0, Double_t = 0.0) const override;
    virtual TF1*     Create(const Char_t* name) const override;
    virtual TF1*     Fit(TH1* hist, Double_t xmin = 0, Double_t xmax = 0, Param_t* param = nullptr) const override;
    virtual TF1*     Fit(TGraph* graph, Double_t xmin = 0, Double_t xmax = 0, Param_t* param = nullptr) const override;
    
  };
  
}

#endif

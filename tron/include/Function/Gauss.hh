#ifndef Tron_Gauss_hh
#define Tron_Gauss_hh

#include "TF1Generator.hh"

namespace Tron {

  struct GaussParam {
    Double_t Height;
    Double_t Mu;
    Double_t Sigma;
  };
  
  class Gauss : public TF1Generator<GaussParam> {
    TRON_SINGLETONIZE(Gauss);
  private:
    static const Int_t kNpar;

  public:
    virtual Int_t    Npar() const override;
    virtual Double_t Func(Double_t* x, Double_t* p) const override;
    virtual Double_t Eval(Double_t x, Double_t height, Double_t mu, Double_t sigma, Double_t = 0.0, Double_t = 0.0, Double_t = 0.0, Double_t = 0.0, Double_t = 0.0, Double_t = 0.0, Double_t = 0.0, Double_t = 0.0) const override;
    virtual TF1*     Create(const Char_t* name) const override;
    virtual TF1*     Fit(TH1* hist, Double_t xmin = 0.0, Double_t xmax = 0.0, Param_t* param = nullptr) const override;
    virtual TF1*     Fit(TGraph* graph, Double_t xmin = 0.0, Double_t xmax = 0.0, Param_t* param = nullptr) const override;
  };
  
}

#endif

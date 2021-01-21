#ifndef Tron_Gauss2D_hh
#define Tron_Gauss2D_hh

#include "TF2Generator.hh"

namespace Tron {

  struct Gauss2DParam {
    Double_t Height;
    Double_t MuX;
    Double_t MuY;
    Double_t SigmaX;
    Double_t SigmaY;
    Double_t Correlation;
  };

  class Gauss2D : public TF2Generator<Gauss2DParam> {
    TRON_SINGLETONIZE(Gauss2D);
  private:
    static const Int_t kNpar;
    
  public:
    virtual Int_t    Npar() const override;
    virtual Double_t Func(Double_t* x, Double_t* p) const override;
    virtual Double_t Eval(Double_t x, Double_t y, Double_t height, Double_t muX, Double_t muY, Double_t sigmaX, Double_t sigmaY, Double_t correlation, Double_t = 0, Double_t = 0, Double_t = 0, Double_t = 0, Double_t = 0) const override;
    virtual TF2*     Create(const Char_t* name) const override;
    virtual TF2*     Fit(TH2* hist, Double_t xmin = 0.0, Double_t xmax = 0.0, Param_t* param = nullptr) const override;
    virtual TF2*     Fit(TGraph2D* graph, Double_t xmin = 0.0, Double_t xmax = 0.0, Param_t* param = nullptr) const override;
    
  };
  
}

#endif

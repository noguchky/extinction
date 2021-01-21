#ifndef Tron_PeriodicSpline3_hh
#define Tron_PeriodicSpline3_hh

#include <vector>
#include "Types.hh"

namespace Tron {

  class PeriodicSpline3 {
  private:
    Double_t              fFrom;
    Double_t              fPeriod;
    std::vector<Double_t> fX;
    std::vector<Double_t> fY;
    std::vector<Double_t> fA;
    std::vector<Double_t> fB;
    std::vector<Double_t> fC;
    std::vector<Double_t> fD;

  public:
    PeriodicSpline3();
    PeriodicSpline3(Int_t n, Double_t* x, Double_t* y);

    inline Double_t GetXmin() const {
      return fFrom;
    }

    inline Double_t GetXmax() const {
      return fFrom + fPeriod;
    }

    inline Double_t GetPeriod() const {
      return fPeriod;
    }

    Double_t operator()(Double_t x) const;

  private:
    Double_t Eval(Double_t x, Int_t i) const;
    void     Build();

  };
  
}

#endif

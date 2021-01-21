#ifndef Tron_TSF1_hh
#define Tron_TSF1_hh

#include "TH1.h"
#include "TF1.h"
#include "TSpline.h"

namespace Tron {

  class TSF1 : public TF1 {
  private:
    TSpline3* fSpline;
  
    static const Int_t kNpar;
    Double_t Function(Double_t* x, Double_t* p);
  
  public:
    TSF1(const Char_t* name, const TH1*      hist,   Double_t xmin = 0, Double_t xmax = 1);
    TSF1(const Char_t* name, const TSpline3* spline, Double_t xmin = 0, Double_t xmax = 1);
    TSF1(const TSF1& sf1);
    virtual ~TSF1();

    inline TSpline3* GetSpline() const {
      return fSpline;
    }
  };

}

#endif

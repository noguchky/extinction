#ifndef Tron_THF1_hh
#define Tron_THF1_hh

#include <vector>
#include "TH1.h"
#include "TF1.h"

namespace Tron {

  class THF1 : public TF1 {
  private:
    Int_t                 fBinLow;
    Int_t                 fBinUp;
    std::vector<Double_t> fLowEdge;
    std::vector<Double_t> fCenter;
    std::vector<Double_t> fContents;
  
    static const Int_t kNpar;
    Double_t Function(Double_t* x, Double_t* p);

  public:
    THF1(const Char_t* name, const TH1* hist, Double_t xmin = 0, Double_t xmax = 1, Bool_t norm = kFALSE);
    virtual ~THF1();
  };

}

#endif

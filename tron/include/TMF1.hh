#ifndef Tron_TMF1_hh
#define Tron_TMF1_hh

#include <string>
#include <vector>
#include "TF1.h"

namespace Tron {

  class TMF1 : public TF1 {
  private:
    std::vector<TF1*>  fFuncs;
    std::vector<Int_t> fIparOffset;

  public:
    TMF1();
    TMF1(const Char_t* name, Int_t n, const TF1** funcs, Double_t xmin = 0, Double_t xmax = 1);
    TMF1(const Char_t* name, const std::initializer_list<const TF1*>& funcs, Double_t xmin = 0, Double_t xmax = 1);
    TMF1(const TMF1& mf1);
    virtual ~TMF1();

    inline  Int_t         GetIpar(Int_t ifun, Int_t ipar) const {
      return fIparOffset[ifun] + ipar;
    }

    inline  Int_t         GetNfun() const {
      return fFuncs.size();
    }

    inline  TF1*          GetFunc(Int_t ifun = 0) {
      return fFuncs[ifun];
    }
    
    inline  const TF1*    GetFunc(Int_t ifun = 0) const {
      return fFuncs[ifun];
    }
    
    virtual void          Draw(Option_t* option = "") override;

    virtual void          Update() override;

  private:
    Double_t              Func(Double_t* x, Double_t* p);
    void                  ExportParams();

  };

}

#endif

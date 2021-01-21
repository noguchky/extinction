#ifndef Tron_ProgressBar_hh
#define Tron_ProgressBar_hh

#include <iostream>
#include <cmath>
#include "Types.hh"

namespace Tron {

  namespace ProgressHelper {

    Bool_t IsLogTiming(Long64_t progress) {
      return progress <= 0L || (progress % (Long64_t)std::pow(10.0, std::floor(std::log10(progress))) == 0L);
    }

  }

  class ProgressBar {
  public:
    static const Double_t kInterval;
    static const Double_t kTotal;
    static const Int_t    kWidth;
    static const Char_t   kDone;
    static const Char_t   kNotYet;

  private:
    enum class ShowType {
      None,
      Real,
      Integer,
    };
    Double_t fPrevious;
    Double_t fInterval;
    Double_t fTotal;
    Int_t    fWidth;
    Char_t   fDone;
    Char_t   fNotYet;
    ShowType fFormat;

  public:
    ProgressBar() {
      fPrevious = 0.0;
      fInterval = kInterval;
      fTotal    = kTotal;
      fWidth    = kWidth;
      fDone     = kDone;
      fNotYet   = kNotYet;
      fFormat   = ShowType::None;
    }

    inline void     SetInterval(Double_t interval) {
      fInterval = interval;
    }

    inline void     SetTotal(Double_t total) {
      fTotal = total;
    }

    inline void     SetWidth(Int_t width) {
      fWidth = width;
    }

    inline void     SetStyle(Char_t done, Char_t notYet) {
      fDone   = done;
      fNotYet = notYet;
    }

    void            Show(Double_t progress);
    void            Show(Int_t    progress);
    void            Terminate();
  };

}

#endif

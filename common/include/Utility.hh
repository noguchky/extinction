#ifndef Extinction_Utility_hh
#define Extinction_Utility_hh

#include "Units.hh"
#include "TCanvas.h"
#include "TH1.h"
#include "TPaveStats.h"

namespace Extinction {

  namespace Utility {

    void ResizeStats(TH1* hist, Double_t factor = 60.0 * perCent) {
      if (gPad) {
        gPad->Update();
        if (TPaveStats* st = dynamic_cast<TPaveStats*>(hist->FindObject("stats"))) {
          const Double_t x1 = st->GetX1();
          const Double_t x2 = st->GetX2();
          const Double_t y1 = st->GetY1();
          const Double_t y2 = st->GetY2();
          st->SetY1(y2 - (y2 - y1) * factor);
          st->SetX1(x2 - (x2 - x1) * factor);
          gPad->Modified();
          gPad->Update();
        }
      }
    }
    
  }
  
}

#endif

#ifndef Tron_TViewFrame_hh
#define Tron_TViewFrame_hh

#include "TH2C.h"
#include "TGraph2D.h"

namespace Tron {

  class TViewFrame2D : public TH2C {
  public:
    TViewFrame2D()
      : TH2C() {
      this->SetDrawOption("AXIS");
      this->SetStats(0);
    }

    TViewFrame2D(const Char_t* name, const Char_t* title,
                 Double_t xmin, Double_t xmax,
                 Double_t ymin, Double_t ymax)
      : TH2C(name, title, 1, xmin, xmax, 1, ymin, ymax) {
      this->SetDrawOption("AXIS");
      this->SetStats(0);
    }

    TViewFrame2D(const TViewFrame2D& frame)
      : TH2C(frame) {
    }

    virtual ~TViewFrame2D() {
    }

  };

  class TViewFrame3D : public TGraph2D  {
  public:
    TViewFrame3D()
      : TGraph2D(2) {
    }

    TViewFrame3D(const Char_t* name, const Char_t* title,
                 Double_t xmin, Double_t xmax,
                 Double_t ymin, Double_t ymax,
                 Double_t zmin, Double_t zmax)
      : TGraph2D(name, title, 2,
                 new Double_t[2] { xmin, xmax },
                 new Double_t[2] { ymin, ymax },
                 new Double_t[2] { zmin, zmax }) {
      this->SetDrawOption("A");
    }

    TViewFrame3D(const TViewFrame3D& frame)
      : TGraph2D(frame) {
    }

    virtual ~TViewFrame3D() {
    }

    TViewFrame2D* Projection(const Char_t* name, Option_t* option = "xy");

  };

}

#endif

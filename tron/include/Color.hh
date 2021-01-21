#ifndef Tron_Color_hh
#define Tron_Color_hh

#include <vector>
#include "TColor.h"
#include "TMath.h"
#include "TVector3.h"
#include "Types.hh"
#include "Math.hh"
#include "Vector3.hh"

namespace Tron {

  class ColorBar {
  private:
    using DataColor      = Tron::Vector3<Double_t>;
    using InterfaceColor = Tron::Vector3<Int_t>;

  private:
    std::vector<DataColor> fColors;
    Double_t               fMinValue = 0.0;
    Double_t               fMaxValue = 1.0;

  public:
    ColorBar() = default;

    inline void SetColors(const std::initializer_list<InterfaceColor>& colors) {
      fColors = std::vector<DataColor>(colors.size());
      std::size_t i = 0;
      for (auto& color : colors) {
        fColors[i].X = Math::Bind(color.X, 0, 255) / 255.;
        fColors[i].Y = Math::Bind(color.Y, 0, 255) / 255.;
        fColors[i].Z = Math::Bind(color.Z, 0, 255) / 255.;
        ++i;
      }
    }

    inline void SetValueRange(Double_t min, Double_t max) {
      fMinValue = min;
      fMaxValue = max;
    }

    inline Int_t GetColor(Double_t value) {
      const Double_t binded = Math::Bind(value, fMinValue, fMaxValue);
      const Double_t diff   = fMaxValue - fMinValue;
      const Double_t nstep  = fColors.size() - 1.0;
      const Int_t    i      = (binded - fMinValue) / diff * nstep;
      const Double_t val1   =  i      * diff / nstep + fMinValue;
      const Double_t val2   = (i + 1) * diff / nstep + fMinValue;
      const auto color = Math::Interpolate(binded, val1, fColors[i], val2, fColors[i + 1]);

      return TColor::GetColor((Float_t)color[0], (Float_t)color[1], (Float_t)color[2]);
    }
  };

  class ColorCircle {
  public:
    enum class RGB { R, G, B };

  private:
    Double_t fMinFrac  = 0.0;
    Double_t fMaxFrac  = 1.0;
    Double_t fMinTheta = 0.0;
    Double_t fMaxTheta = 1.0;
    Double_t fMinValue = 0.0;
    Double_t fMaxValue = 1.0;

  public:
    ColorCircle() = default;

    inline void SetFracRange(Double_t min, Double_t max) {
      fMinFrac  = Math::Bind(min, 0.0, 1.0);
      fMaxFrac  = Math::Bind(max, 0.0, 1.0);
    }

    inline void SetThetaRange(Double_t min, Double_t max) {
      fMinTheta = Math::Bind(min, 0.0, 1.0);
      fMaxTheta = Math::Bind(max, 0.0, 1.0);
    }

    inline void SetValueRange(Double_t min, Double_t max) {
      fMinValue = min;
      fMaxValue = max;
    }

    inline Double_t GetFrac(RGB rgb, Double_t theta) {
      theta += 1.0 / 6.0 + (Double_t)rgb / 3.0;
      theta += - TMath::Floor(theta);
      theta *= 6.0;
      Double_t frac;
      if      (theta < 2.0) { frac = 1.0;         } // 1.0
      else if (theta < 3.0) { frac = 3.0 - theta; } // 1.0->0.0
      else if (theta < 5.0) { frac = 0.0;         } // 0.0
      else                  { frac = theta - 5.0; } // 0.0->1.0
      return frac * (fMaxFrac - fMinFrac) + fMinFrac;
    }

    inline Int_t GetColor(Double_t value) {
      const Double_t binded     = Math::Bind(value, fMinValue, fMaxValue);
      const Double_t normalized = (binded - fMinValue) / (fMaxValue - fMinValue);
      const Double_t theta      = normalized * (fMaxTheta - fMinTheta) + fMinTheta;
      const Float_t r = GetFrac(RGB::R, theta);
      const Float_t g = GetFrac(RGB::G, theta);
      const Float_t b = GetFrac(RGB::B, theta);
      return TColor::GetColor(r, g, b);
    }
  };

}  

#endif

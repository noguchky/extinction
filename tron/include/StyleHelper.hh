#ifndef Tron_StyleHelper_hh
#define Tron_StyleHelper_hh

#include <set>
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGaxis.h"

namespace Tron {

  namespace StyleHelper {

    enum class BaseStyle  { Analyze, Visual };

    enum class CanvasAlign { None, Left, Center, Right };

    struct OptStatStyle {
      enum {
            None        =          0,
            Default     =          1,
            Name        = 1000000001,
            Entries     =         10,
            Mean        =        100,
            MeanErr     =        200,
            StdDev      =       1000,
            StdDevErr   =       2000,
            Underflow   =      10000,
            Overflow    =     100000,
            Integral    =    1000000,
            IntegralW   =    2000000,
            Skewness    =   10000000,
            SkewnessErr =   20000000,
            Kurtosis    =  100000000,
            KurtosisErr =  200000000,
      };
    };

    struct OptFitStyle {
      enum {
            None           =     0,
            Default        =     1,
            Probability    = 10001,
            Chisquare      =    10,
            Errors         =   100,
            NonFixedParams =  1000,
            AllParams      =  2000,
      };
    };

    enum class TextAlignH { Right = 30, Center = 20, Left = 10 };
    enum class TextAlignV { Top = 3, Center = 2, Bottom = 1 };

    Int_t operator+(TextAlignH horizontal, TextAlignV vertical) {
      return (Int_t)horizontal + (Int_t)vertical;
    }
    Int_t operator+(TextAlignV vertical, TextAlignH horizontal) {
      return (Int_t)horizontal + (Int_t)vertical;
    }

    inline void SetBaseStyle(BaseStyle style) {
      switch (style) {
      case BaseStyle::Analyze:
        gStyle->SetOptStat(OptStatStyle::Default);
        gStyle->SetOptFit (OptFitStyle ::Default);
        break;
      case BaseStyle::Visual:
        gStyle->SetOptStat(OptStatStyle::None);
        gStyle->SetOptFit (OptFitStyle ::None);
        break;
      default:
        break;
      }
      gStyle->SetNdivisions(505, "x");
      gStyle->SetNdivisions(505, "y");
      gStyle->SetNdivisions(505, "z");
      gStyle->SetStripDecimals(false);
      TGaxis::SetMaxDigits(3);
    }

    inline void SetCanvasSquare(TVirtualPad* pad, Double_t margin, CanvasAlign align = CanvasAlign::Center) {
      const Double_t w = pad->GetWw();
      const Double_t h = pad->GetWh();
      const Double_t restH = (1.0 - 2.0 * margin) * h;
      Double_t leftMargin, rightMargin;
      switch (align) {
      case CanvasAlign::Center:
        leftMargin  = (w - restH) / w / 2.;
        rightMargin = leftMargin;
        break;
      case CanvasAlign::Left:
        leftMargin  = margin;
        rightMargin = (w - restH) / w - margin;
        break;
      case CanvasAlign::Right:
        leftMargin  = (w - restH) / w - margin;
        rightMargin = margin;
        break;
      default:
        leftMargin  = margin;
        rightMargin = margin;
        break;
      }
      pad->SetMargin(leftMargin, rightMargin, margin, margin);
    }

    inline void SetMarkerStyle(TAttMarker* marker, Color_t color, Style_t style, Size_t size) {
      marker->SetMarkerColor(color);
      marker->SetMarkerStyle(style);
      marker->SetMarkerSize (size );
    }

    inline void SetLineStyle(TAttLine* line, Color_t color, Style_t style, Width_t width) {
      line->SetLineColor(color);
      line->SetLineStyle(style);
      line->SetLineWidth(width);
    }

    inline void SetFillStyle(TAttFill* fill, Color_t color, Style_t style) {
      fill->SetFillColor(color);
      fill->SetFillStyle(style);
    }

  }

}

#endif

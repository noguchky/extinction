#ifndef Tron_MathUnits_hh
#define Tron_MathUnits_hh

#include "Types.hh"

namespace Tron {

  //
  //
  //
  static constexpr Double_t     pi      = 3.14159265358979323846;
  static constexpr Double_t  twopi      = 2 * pi;
  static constexpr Double_t halfpi      = pi / 2;
  static constexpr Double_t     pi2     = pi * pi;

  static constexpr Double_t goldenRatio = 1.6180339887498948482045868343656;

  //
  // Angle
  //
  static constexpr Double_t radian      = 1.;
  static constexpr Double_t milliradian = 1.e-3 * radian;
  static constexpr Double_t degree      = (pi / 180.0) * radian;

  static constexpr Double_t   steradian = 1.;

  // symbols
  static constexpr Double_t rad         = radian;
  static constexpr Double_t mrad        = milliradian;
  static constexpr Double_t sr          = steradian;
  static constexpr Double_t deg         = degree;

  //
  // Miscellaneous
  //
  static constexpr Double_t perCent          = 0.01;
  static constexpr Double_t perThousand      = 0.001;
  static constexpr Double_t perMillion       = 0.000001;

}

#endif

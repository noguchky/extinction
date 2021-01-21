#ifndef Extinction_Units_hh
#define Extinction_Units_hh

#include "Rtypes.h"

namespace Extinction {

  const Double_t  sec = 1.0;
  const Double_t msec = 1.0e-3 * sec;
  const Double_t usec = 1.0e-6 * sec;
  const Double_t nsec = 1.0e-9 * sec;

  const Double_t  Hz = 1.0 / sec;
  const Double_t kHz = 1.0e3 * Hz;
  const Double_t MHz = 1.0e6 * Hz;
  const Double_t GHz = 1.0e9 * Hz;

  const Double_t  m = 1.0;
  const Double_t cm = 1.0e-2 * m;
  const Double_t mm = 1.0e-3 * m;
  const Double_t um = 1.0e-6 * m;
  const Double_t nm = 1.0e-9 * m;

  const Double_t  m2 = 1.0 *  m *  m;
  const Double_t cm2 = 1.0 * cm * cm;
  const Double_t mm2 = 1.0 * mm * mm;
  const Double_t um2 = 1.0 * um * um;
  const Double_t nm2 = 1.0 * nm * nm;

  const Double_t  m3 = 1.0 *  m *  m *  m;
  const Double_t cm3 = 1.0 * cm * cm * cm;
  const Double_t mm3 = 1.0 * mm * mm * mm;
  const Double_t um3 = 1.0 * um * um * um;
  const Double_t nm3 = 1.0 * nm * nm * nm;

  const Double_t perCent = 0.01;

}

#endif

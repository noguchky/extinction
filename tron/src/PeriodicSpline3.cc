#include <iostream>
#include <cmath>
#include "PeriodicSpline3.hh"
#include "Math.hh"
#include "VMatrix.hh"

Tron::PeriodicSpline3::PeriodicSpline3() {      
}

Tron::PeriodicSpline3::PeriodicSpline3(Int_t n, Double_t* x, Double_t* y) {
  fX = std::vector<Double_t>(x, x + n);
  fY = std::vector<Double_t>(y, y + n);
  Build();
}

Double_t Tron::PeriodicSpline3::operator()(Double_t x) const {
  if (!fA.size()) {
    return 0.0;
  }
  const Double_t xx = Math::PeriodicBind(x, fFrom, fFrom + fPeriod);
  const Int_t    ii = Math::FindIndex(xx, fX.size(), fX.data());
  return Eval(xx, ii);
}

Double_t Tron::PeriodicSpline3::Eval(Double_t x, Int_t i) const {
  const Double_t dx = x - fX[i];
  return fA[i] * dx * dx * dx + fB[i] * dx * dx + fC[i] * dx + fD[i];
}

void Tron::PeriodicSpline3::Build() {
  fFrom   = fX.front();
  fPeriod = fX.back() - fFrom;

  const UInt_t n = fX.size() - 1;

  // Calc Array h, g
  std::vector<Double_t> h(n);
  std::vector<Double_t> g(n);
  for (UInt_t i = 0; i < n; ++i) {
    h[i] = fX[i + 1] - fX[i];
    g[i] = fY[i + 1] - fY[i];

    if (h[i] <= 0.0) {
      std::cout << "PeriodicSpline3::Build() [error] x array must be ascending order" << std::endl;
      return;
    }
  }

  // Calc Matrix M
  VMatrix<Double_t> m(n, n);
  for (UInt_t i = 0; i < n; ++i) {
    const UInt_t ip1 = (i + 1) % n;
    const UInt_t ip2 = (i + 2) % n;
    m(i, i  ) = h[i];
    m(i, ip1) = 2.0 * (h[i] + h[ip1]);
    m(i, ip2) = h[ip1];
  }

  // Calc Vector v
  VMatrix<Double_t> v(n, 1);
  for (UInt_t i = 0; i < n; ++i) {
    const UInt_t ip1 = (i + 1) % n;
    v[i] = (g[ip1] / h[ip1] - g[i] / h[i]) * 6.0;
  }

  // Calc Vector u
  VMatrix<Double_t> u = m.Inverse() * v;
            
  // Calc Array A, B, C, D
  fA = std::vector<Double_t>(n);
  fB = std::vector<Double_t>(n);
  fC = std::vector<Double_t>(n);
  fD = std::vector<Double_t>(n);
  for (UInt_t i = 0; i < n; ++i) {
    const UInt_t ip1 = (i + 1) % n;
    fA[i] = (u[ip1] - u[i]) / h[i] / 6.0;
    fB[i] = u[i] * 0.5;
    fC[i] = g[i] / h[i] - (u[ip1] + 2.0 * u[i]) * h[i] / 6.0;
    fD[i] = fY[i];
  }
}

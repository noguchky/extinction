#include <cmath>
#include "TGraph.h"
#include "Integrator.hh"

double Tron::Integrator::Integrator::Integrate(const Func_t&, double&, double, double, double) const {
  return 0.0;
}

double Tron::Integrator::Integrator::Integrate(const Func_t&, double, double, double) const {
  return 0.0;
}

double Tron::Integrator::TrapezoidalMethod::Integrate(const Func_t& func, double from, double to, double dx) const {
  double x = from;
  return Integrate(func, x, from, to, dx);
}

double Tron::Integrator::TrapezoidalMethod::Integrate(const Func_t& func, double& x, double from, double to, double dx) const {
  //--- Convert dx
  dx = std::abs(dx) * (to - from > 0.0 ? +1.0 : -1.0);

  //--- Summation
  double sum = 0.0;
  x = from;
  sum -= func(x) * 0.5;
  if (from < to) {
    while (x < to) {
      sum += func(x);
      x += dx;
    }
  } else {
    while (x > to) {
      sum += func(x);
      x += dx;
    }
  }
  sum += func(x) * 0.5;

  return sum * dx;
}

double Tron::Integrator::SimpsonMethod::Integrate(const Func_t& func, double from, double to, double dx) const {
  double x = from;
  return Integrate(func, x, from, to, dx);
}

double Tron::Integrator::SimpsonMethod::Integrate(const Func_t& func, double& x, double from, double to, double dx) const {
  //--- Convert dx
  dx = std::abs(dx) * (to - from > 0.0 ? +1.0 : -1.0);

  //--- Summation
  double sum = 0.0;
  x = from;
  sum -= func(x) * 0.5;
  if (from < to) {
    while (x < to) {
      sum += func(x);
      x += dx;
      sum += func(x) * 2.0;
      x += dx;
    }
  } else {
    while (x > to) {
      sum += func(x);
      x += dx;
      sum += func(x) * 2.0; 
      x += dx;
    }
  }
  sum += func(x) * 0.5;

  return sum * dx * 2.0 / 3.0;
}

double Tron::Integrator::PowerSimpsonMethod::Integrate(const Func_t& func, double from, double to, double nx) const {
  double x = from;
  return Integrate(func, x, from, to, nx);
}

double Tron::Integrator::PowerSimpsonMethod::Integrate(const Func_t& func, double& x, double from, double to, double nx) const {
  if (std::abs(from - to) <= std::numeric_limits<double>::epsilon()) {
    return 0.0;
  }

  const double sign = from * to;
  if (sign < 0.0) {
    return 0.0;
  }

  //--- Convert dx
  nx = std::ceil(nx);

  //--- Summation
  double  sum = 0.0;
  double dsum = 0.0;
  const double ax = std::pow(to / from, 1.0 / nx);
  const double term1 = (ax + 1.0) / ax;
  const double term2 = 4.0;
  const double term3 = 0.5 * (ax + 1.0) * from;

  sum += func(from);
  sum += term2 * func(term3); 
  double powa = 1.0;
  for (long long i = 1; i < nx; ++i) {
    powa *= ax;
    sum += powa * (term1 * func(powa * from ) +
                   term2 * func(powa * term3));
  }
  sum +=  powa / ax * func(powa * from);
  dsum = (to - powa * from) * func((to + powa * from) * 0.5);
  x = to;

  return sum * (ax - 1.0) * from / 6.0 + dsum;
}

double Tron::Integrator::LineGraphMethod::Integrate(const TGraph* graph, Double_t from, Double_t to) const {
  //--- Summation
  const int     n  = graph->GetN();
  const double* xs = graph->GetX();
  const double* ys = graph->GetY();

  double sum = 0.0;
  if (from < to) {
    int ip = 1;
    for (; ip < n; ++ip) {
      if (from < xs[ip]) {
        const double yy = (ys[ip] - ys[ip - 1]) / (xs[ip] - xs[ip - 1]) * (from - xs[ip]) + ys[ip];
        sum += (ys[ip] + yy) * (xs[ip] - from);
        break;
      }
    }
    for (++ip; ip < n; ++ip) {
      if (xs[ip] < to) {
        sum += (ys[ip] + ys[ip - 1]) * (xs[ip] - xs[ip - 1]);
      } else {
        const double yy = (ys[ip] - ys[ip - 1]) / (xs[ip] - xs[ip - 1]) * (to - xs[ip]) + ys[ip];
        sum += (yy + ys[ip - 1]) * (to - xs[ip - 1]);
        break;
      }
    }
  } else {
    int ip = n - 2;
    for (; 0 <= ip; --ip) {
      if (xs[ip] < from) {
        const double yy = (ys[ip + 1] - ys[ip]) / (xs[ip + 1] - xs[ip]) * (from - xs[ip]) + ys[ip];
        sum += (yy + ys[ip]) * (from - xs[ip]);
        break;
      }
    }
    for (--ip; 0 <= ip; --ip) {
      if (to < xs[ip]) {
        sum += (ys[ip + 1] + ys[ip]) * (xs[ip + 1] - xs[ip]);
      } else {
        const double yy = (ys[ip + 1] - ys[ip]) / (xs[ip + 1] - xs[ip]) * (to - xs[ip]) + ys[ip];
        sum += (ys[ip + 1] + yy) * (xs[ip + 1] - to);
        break;
      }
    }
  }

  return sum * 0.5;
}

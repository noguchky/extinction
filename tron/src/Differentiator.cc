#include "Differentiator.hh"

double Tron::Differentiator::Differentiator::Differentiate(const Func_t&, double, double) const {
  return 0.0;
}

double Tron::Differentiator::NewtonsDifference::Differentiate(const Func_t& func, double x, double dx) const {
  const double df = (func(x + dx) - func(x));
  return df / dx;
}

double Tron::Differentiator::SymmetricDifference::Differentiate(const Func_t& func, double x, double dx) const {
  const double df = (func(x + dx) - func(x - dx)) * 0.5;
  return df / dx;
}

double Tron::Differentiator::FivePointStencil::Differentiate(const Func_t& func, double x, double dx) const {
  const double df = (func(x - 2.0 * dx) + 8.0 * func(x - dx) + 8.0 * func(x + dx) + func(x + 2.0 * dx)) / 12.0;
  return df / dx;
}

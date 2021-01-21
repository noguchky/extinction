#include "RootFinder.hh"
#include "Math.hh"

double Tron::RootFinder::RootFinder::Find(const Func_t&, double, double, double, double) const {
  return 0.0;
}

double Tron::RootFinder::LinearSearchMethod::Find(const Func_t& func, double from, double to, double dx, double) const {
  //--- Convert dx, dy
  dx = std::abs(dx) * (to - from > 0.0 ? +1.0 : -1.0);

  //--- Get Range
  double solution = to + dx;
  double x1 = from, x2 = x1;
  double y1, y2 = func(x1);
  while (Math::Between(x1, from, to)) {
    x1 = x2, x2 = x1 + dx;
    y1 = y2, y2 = func(x2);
    //--- Crossing Check
    if (y1 * y2 <= 0) {
      solution = (x1 + x2) * 0.5;
      break;
    }
  }
  
  return solution;
}

double Tron::RootFinder::BisectionMethod::fRoughSearchFactor = 500.0;

double Tron::RootFinder::BisectionMethod::Find(const Func_t& func, double from, double to, double dx, double dy) const {
  //--- Convert dx, dy
  dx = std::abs(dx) * (to - from > 0.0 ? +1.0 : -1.0);
  dy = std::abs(dy);

  //--- Rough Find Root
  static const auto& finder = LinearSearchMethod::Get();
  double root = finder.Find(func, from, to, fRoughSearchFactor * dx, dy);
  
  //--- Get Range
  double from2 = Math::Bind(root - fRoughSearchFactor * dx, from, to);
  double to2   = Math::Bind(root + fRoughSearchFactor * dx, from, to);
  
  const double xmin = std::min(from2, to2);
  const double xmax = std::max(from2, to2);

  double solution = to + dx;
  double xm, ym;
  double x1 = xmin, y1 = func(x1);
  double x2 = xmax, y2 = func(x2);
  if (y1 * y2 <= 0) {
    while (std::abs(x1 - x2) > dx || std::abs(y1 - y2) > dy) {
      xm = (x1 + x2) * 0.5;
      ym = func(xm);

      if (y1 * ym <= 0) {
        x2 = xm; y2 = ym;
      } else if (y2 * ym <= 0) {
        x1 = xm; y1 = ym;
      } else {
        break;
      }
    }
    
    //--- Solution Check
    if (std::abs(x1 - x2) <= dx && std::abs(y1 - y2) <= dy) {
      solution = (x1 + x2) * 0.5;
    }
  }
  
  return solution;
}
  
// double Tron::RootFinder::SecantMethod::Find(const Func_t& func, double from, double to, double dx, double dy) const {
//   return 0.0;
// }
  
// double Tron::RootFinder::NewtonMethod::Find(const Func_t& func, double from, double to, double dx, double dy) const {
//   return 0.0;  
// }

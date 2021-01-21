#include "ExtremumSearcher.hh"
#include "Math.hh"

Tron::ExtremumSearcher::Result_t Tron::ExtremumSearcher::ExtremumSearcher::FindMinimum(const Func_t&, double, double, double, double) const {
  return Result_t(0.0, 0.0);
}

Tron::ExtremumSearcher::Result_t Tron::ExtremumSearcher::ExtremumSearcher::FindMaximum(const Func_t&, double, double, double, double) const {
  return Result_t(0.0, 0.0);
}

Tron::ExtremumSearcher::Result_t Tron::ExtremumSearcher::LinearSearchMethod::FindMinimum(const Func_t& func, double from, double to, double dx, double) const {
  double x, y, minx, miny;  

  minx = from, miny = func(minx);
  if (from < to) {
    for (x = from + dx; x < to; x += dx) {
      y = func(x);
      if (y < miny) {
        minx = x;
        miny = y;
      }
    }
  } else {
    for (x = from - dx; x < to; x -= dx) {
      y = func(x);
      if (y < miny) {
        minx = x;
        miny = y;
      }
    }
  }

  return Result_t(minx, miny);
}

Tron::ExtremumSearcher::Result_t Tron::ExtremumSearcher::LinearSearchMethod::FindMaximum(const Func_t& func, double from, double to, double dx, double) const {
  double x, y, maxx, maxy;  

  maxx = from, maxy = func(maxx);
  if (from < to) {
    for (x = from + dx; x < to; x += dx) {
      y = func(x);
      if (y > maxy) {
        maxx = x;
        maxy = y;
      }
    }
  } else {
    for (x = from - dx; x < to; x -= dx) {
      y = func(x);
      if (y > maxy) {
        maxx = x;
        maxy = y;
      }
    }
  }

  return Result_t(maxx, maxy);
}

double Tron::ExtremumSearcher::GoldenSectionSearchMethod::fRoughSearchFactor = 500.0;

Tron::ExtremumSearcher::Result_t Tron::ExtremumSearcher::GoldenSectionSearchMethod::FindMinimum(const Func_t& func, double from, double to, double dx, double dy) const {
  // Rough Search
  static const auto& searcher = LinearSearchMethod::Get();
  auto extremum = searcher.FindMinimum(func, from, to, fRoughSearchFactor * dx, fRoughSearchFactor * dy);

  double x1 = extremum.X - fRoughSearchFactor * dx;
  double x4 = extremum.X + fRoughSearchFactor * dx;;
  double x2 = Math::Average(x1, x4, goldenRatio, 1.0);
  double x3 = Math::Average(x1, x4, 1.0, goldenRatio);

  double y1 = func(x1);
  double y2 = func(x2);
  double y3 = func(x3);
  double y4 = func(x4);

  do {
    if (y2 < y3) {
      x4 = x3, x3 = x2, x2 = Math::Average(x1, x4, goldenRatio, 1.0);
      y4 = y3, y3 = y2, y2 = func(x2);
    } else {
      x1 = x2, x2 = x3, x3 = Math::Average(x1, x4, 1.0, goldenRatio);
      y1 = y2, y2 = y3, y3 = func(x3);
    }
  } while (std::abs(x1 - x4) > dx || std::abs(y1 - y4) > dy);

  return y2 < y2 ? Result_t(x2, y2) : Result_t(x3, y3);
}

Tron::ExtremumSearcher::Result_t Tron::ExtremumSearcher::GoldenSectionSearchMethod::FindMaximum(const Func_t& func, double from, double to, double dx, double dy) const {
  // Rough Search
  static const auto& searcher = LinearSearchMethod::Get();
  auto extremum = searcher.FindMaximum(func, from, to, fRoughSearchFactor * dx, fRoughSearchFactor * dy);

  double x1 = extremum.X - fRoughSearchFactor * dx;
  double x4 = extremum.X + fRoughSearchFactor * dx;
  double x2 = Math::Average(x1, x4, goldenRatio, 1.0);
  double x3 = Math::Average(x1, x4, 1.0, goldenRatio);

  double y1 = func(x1);
  double y2 = func(x2);
  double y3 = func(x3);
  double y4 = func(x4);

  do {
    if (y2 < y3) {
      x1 = x2, x2 = x3, x3 = Math::Average(x1, x4, 1.0, goldenRatio);
      y1 = y2, y2 = y3, y3 = func(x3);
    } else {
      x4 = x3, x3 = x2, x2 = Math::Average(x1, x4, goldenRatio, 1.0);
      y4 = y3, y3 = y2, y2 = func(x2);
    }
  } while (std::abs(x1 - x4) > dx || std::abs(y1 - y4) > dy);

  return y2 > y3 ? Result_t(x2, y2) : Result_t(x3, y3);
}

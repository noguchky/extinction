#ifndef Tron_LinearLeastSquareMethod_hh
#define Tron_LinearLeastSquareMethod_hh

#include "Types.hh"
#include "Param.hh"
#include "Math.hh"

namespace Tron {

  template <typename Number_t>
  struct LinearFittingParam_t {
    Param_t<Number_t> Intercept;
    Param_t<Number_t> Slope;

    Param_t<Number_t> Eval(Number_t x) {
      const Number_t value  = Intercept.Value + Slope.Value * x;
      const Number_t errorx = Intercept.Error;
      const Number_t errory = Slope.Error * x;
      const Number_t error  = Math::SqrtOfSumOfSquared(errorx, errory);
      return Param_t<Number_t>(value, error);
    }

    Param_t<Number_t> Find(Double_t y) {
      const Number_t value  = (y - Intercept.Value) / Slope.Value;
      const Number_t errorx = Intercept.Error / Slope.Value;
      const Number_t errory = Slope.Error * (y - Intercept.Value) / Slope.Value / Slope.Value;
      const Number_t error  = Math::SqrtOfSumOfSquared(errorx, errory);
      return Param_t<Number_t>(value, error);
    }

  };

  namespace LinearLeastSquareMethod {

    template <typename Number_t>
    LinearFittingParam_t<Number_t> Fit(Int_t n, Number_t* x, Number_t* y);

    template <typename Number_t>
    LinearFittingParam_t<Number_t> Fit(Int_t n, Number_t* x, Number_t* y, Number_t* yerr);

    // template <typename Number_t>
    // LinearFittingParam_t<Number_t> Fit(Int_t n, Number_t* x, Number_t* y, Number_t* xerr, Number_t* yerr);

  }

}

template <typename Number_t>
Tron::LinearFittingParam_t<Number_t> Tron::LinearLeastSquareMethod::Fit(Int_t n, Number_t *x, Number_t *y) {
  Number_t sumN = 0.0, sumX = 0.0, sumY = 0.0, sumXX = 0.0, sumXY = 0.0;
  for (Int_t i = 0; i < n; ++i, ++x, ++y) {
    sumN  += 1.0;
    sumX  += *x;
    sumY  += *y;
    sumXX += *x * *x;
    sumXY += *x * *y;
  }

  LinearFittingParam_t<Number_t> result;
  {
    const Number_t denominator = 1.0 / (sumN * sumXX - sumX * sumX);
    result.Intercept.Value = (sumN  * sumXY - sumX  * sumY) * denominator;
    result.Slope    .Value = (sumXX * sumY  - sumXY * sumX) * denominator;
    result.Intercept.Error = std::sqrt(sumN  * denominator);
    result.Slope    .Error = std::sqrt(sumXX * denominator);
  }

  return result;
}

template <typename Number_t>
Tron::LinearFittingParam_t<Number_t> Tron::LinearLeastSquareMethod::Fit(Int_t n, Number_t *x, Number_t *y, Number_t* yerr) {
  Number_t sumN = 0.0, sumX = 0.0, sumY = 0.0, sumXX = 0.0, sumXY = 0.0;
  for (Int_t i = 0; i < n; ++i, ++x, ++y, ++yerr) {
    const Number_t weight = 1.0 / *yerr / *yerr;
    sumN  += 1.0     * weight;
    sumX  += *x      * weight;
    sumY  += *y      * weight;
    sumXX += *x * *x * weight;
    sumXY += *x * *y * weight;
  }

  LinearFittingParam_t<Number_t> result;
  {
    const Number_t denominator = 1.0 / (sumN * sumXX - sumX * sumX);
    result.Intercept.Value = (sumN  * sumXY - sumX  * sumY) * denominator;
    result.Slope    .Value = (sumXX * sumY  - sumXY * sumX) * denominator;
    result.Intercept.Error = std::sqrt(sumN  * denominator);
    result.Slope    .Error = std::sqrt(sumXX * denominator);
  }

  return result;
}

#endif

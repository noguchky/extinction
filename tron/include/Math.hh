#ifndef Tron_Math_hh
#define Tron_Math_hh

#include <numeric>
#include <cmath>
#include <tuple>
#include <algorithm>
#include "Types.hh"
#include "MathUnits.hh"

namespace Tron {

  namespace Math {
    //==================================================
    template <typename Number_t>
    inline Number_t RevertEndian(Number_t x) {
      union {
        Number_t value;
        UChar_t  bytes[sizeof(Number_t)];
      } buff;

      buff.value = x;
      std::reverse(std::begin(buff.bytes), std::end(buff.bytes));
      return buff.value;
    }

    //==================================================
    template <typename Number_t>
    inline Number_t Bind(Number_t x, Number_t xmin, Number_t xmax) {
      return xmin < xmax ? (x < xmin ? xmin : (x > xmax ? xmax : x)) : (x < xmax ? xmax : (x > xmin ? xmin : x));
    }

    template <typename Number_t>
    inline Number_t PeriodicBind(Number_t x, Number_t xmin, Number_t xmax) {
      x = std::fmod(x - xmin, xmax - xmin);
      if (x < 0.0) {
        x += xmax - xmin;
      }
      return x + xmin;
    }

    template <typename Number_t>
    inline Bool_t Between(Number_t x, Number_t xmin, Number_t xmax) {
      return xmin < xmax ? (xmin <= x && x <= xmax) : (xmax <= x && x <= xmin);
    }

    //==================================================
    template <typename Number_t>
    inline Number_t Greatest(Number_t first) {
      return first;
    }

    template <typename Number_t, typename... Numbers_t>
    inline Number_t Greatest(Number_t first, Numbers_t... rest) {
      return std::max(first, Greatest(rest...));
    }

    template <typename Number_t>
    inline Number_t Least(Number_t first) {
      return first;
    }

    template <typename Number_t, typename... Numbers_t>
    inline Number_t Least(Number_t first, Numbers_t... rest) {
      return std::min(first, Least(rest...));
    }

    //==================================================
    template <typename Number_t>
    inline std::size_t MaxIndex(std::size_t n, const Number_t* array) {
      auto iter = std::max_element(array, array + n);
      return std::distance(array, iter);
    }

    template <typename Number_t>
    inline std::size_t MinIndex(std::size_t n, const Number_t* array) {
      auto iter = std::min_element(array, array + n);
      return std::distance(array, iter);
    }

    template <typename Number_t>
    inline Int_t FindIndex(Number_t x, std::size_t n, const Number_t* sortedArray) {
      std::size_t imin = 0;
      std::size_t imax = n - 1;
      for (std::size_t limit = 0; limit < n; ++limit) {
        const Int_t imid = (imin + imax) * 0.5;
        if        (x < sortedArray[imid    ]     ) {
          imax = imid;
        } else if (    sortedArray[imid + 1] <= x) {
          imin = imid + 1;
        } else { // sortedArray[imid] <= x < sortedArray[imid + 1]
          return imid;
        }
        if         (imax == 0    ) {
          return -1;
        } else  if (imin == n - 1) {
          return n - 1;
        }
      }
      return -1;
    }

    //==================================================
    template <typename Number_t>
    Number_t Sum(std::size_t n, const Number_t* array) {
      return std::accumulate(array, array + n, Number_t());
    }

    template <typename Number_t>
    inline Number_t SumOfSquared(std::size_t n, Number_t* array) {
      return std::accumulate(array, array + n, Number_t(), [] (Number_t sum, Number_t val) { return sum += val * val; });
    }

    template <typename Number_t>
    inline Number_t SumOfSquared(Number_t first) {
      return first * first;
    }

    template <typename Number_t, typename... Numbers_t>
    inline Number_t SumOfSquared(Number_t first, Numbers_t... rest) {
      return SumOfSquared(first) + SumOfSquared(rest...);
    }

    template <typename Number_t>
    inline Number_t SqrtOfSumOfSquared(std::size_t n, Number_t* array) {
      return std::sqrt(std::accumulate(array, array + n, Number_t(), [](Number_t sum, Number_t val) { return sum + val * val; }));
    }

    template <typename Number_t, typename... Numbers_t>
    inline Number_t SqrtOfSumOfSquared(Number_t first, Numbers_t... rest) {
      return std::sqrt(SumOfSquared<Double_t>(first, rest...));
    }

    //==================================================
    template <typename Value_t>
    inline Value_t Average(Value_t x1, Value_t x2, Double_t r1 = 1.0, Double_t r2 = 1.0) {
      return (x1 * r1 + x2 * r2) / (r1 + r2);
    }

    template <typename Value_t>
    Value_t Average(std::size_t n, const Value_t* x) {
      return std::accumulate(x, x + n, Value_t()) / n;
    }

    template <typename Value_t>
    Value_t Average(std::size_t n, const Value_t* x, const Double_t* w) {
      Value_t   sumx = Value_t();
      Double_t  sumw = 0.0;
      Double_t* wIt = (Double_t*)w;
      for (auto xIt = x, xEnd = x + n;
           xIt != xEnd;
           sumx += *xIt++ * *wIt, sumw += *wIt++);
      return sumx / sumw;
    }

    template <typename Value_t>
    inline Value_t Interpolate(Double_t x, Double_t x1, Value_t y1, Double_t x2, Value_t y2) {
      return Average(y1, y2, x2 - x, x - x1);
    }

    template <typename Number_t>
    void MovingAverage(std::size_t l, std::size_t n, Number_t* from, Number_t* to) {
      const std::size_t m = 2 * l + 1;

      for (auto tIt = to, fIt = from, fEnd = from + l;
           fIt != fEnd;
           *tIt++ = *fIt++);

      for (auto tIt = to + l, fIt = from, fEnd = from + n - 2 * l;
           fIt != fEnd;
           *tIt++ = Sum(m, fIt++) / m);

      for (auto tIt = to + n - 1 - l, fIt = from + n - 1 - l, fEnd = from + n;
           fIt != fEnd;
           *tIt++ = *fIt++);
    }

    //==================================================
    template <typename Number_t>
    Double_t Variance(std::size_t n, const Number_t* x, const Number_t* pMuX = nullptr) {
      const Double_t muX = pMuX ? *pMuX : Average(n, x);

      Double_t sigmaXX = 0.0;
      for (std::size_t i = 0; i < n; ++i) {
        sigmaXX += (x[i] - muX) * (x[i] - muX);
      }
      sigmaXX /= n;

      return sigmaXX;
    }

    template <typename Number_t>
    std::tuple<Double_t, Double_t, Double_t> Variance(std::size_t n,
                                                      const Number_t* x,
                                                      const Number_t* y,
                                                      const Number_t* pMuX = nullptr,
                                                      const Number_t* pMuY = nullptr) {
      const Double_t muX = (pMuX && pMuY) ? *pMuX : Average(n, x);
      const Double_t muY = (pMuX && pMuY) ? *pMuY : Average(n, y);

      Double_t sigmaXX = 0.0, sigmaYY = 0.0, sigmaXY = 0.0;
      for (std::size_t i = 0; i < n; ++i) {
        sigmaXX += (x[i] - muX) * (x[i] - muX);
        sigmaXY += (x[i] - muX) * (y[i] - muY);
        sigmaYY += (y[i] - muY) * (y[i] - muY);
      }
      sigmaXX /= n;
      sigmaXY /= n;
      sigmaYY /= n;

      return { sigmaXX, sigmaXY, sigmaYY };
    }

    //==================================================
    inline Bool_t IsSquared(Long_t n) {
      switch (n % 256) {
      case   0:  case   1:  case   4:  case   9:  case  16:  case  17:  case  25:  case  33:  case  36:  case  41:
      case  49:  case  57:  case  64:  case  65:  case  68:  case  73:  case  81:  case  89:  case  97:  case 100:
      case 105:  case 113:  case 121:  case 129:  case 132:  case 137:  case 144:  case 145:  case 153:  case 161:
      case 164:  case 169:  case 177:  case 185:  case 193:  case 196:  case 201:  case 209:  case 217:  case 225:
      case 228:  case 233:  case 241:  case 249:
        break;
      default:
        return false;
      }

      switch (n % 9) {
      case 0:  case 1:  case 4:  case 7:
        break;
      default:
        return false;
      }

      switch (n % 5) {
      case 0:  case 1: case 4:
        break;
      default:
        return false;
      }

      switch (n % 7) {
      case 0:  case 1: case 2: case 4:
        break;
      default:
        return false;
      }

      switch (n % 13) {
      case  0:  case  1:  case  3:  case  4:  case  9:  case 10:  case 12:
        break;
      default:
        return false;
      }

      switch (n % 17) {
      case  0:  case  1:  case  2:  case  4:  case  8:  case  9:  case 13:  case 15:  case 16:
        break;
      default:
        return false;
      }

      const Long_t sqrt = std::sqrt((Double_t)n);
      return sqrt * sqrt == n;
    }

  }

}

#endif

#ifndef Tron_ExtremumSearcher_hh
#define Tron_ExtremumSearcher_hh

#include <functional>
#include "Types.hh"
#include "Singleton.hh"
#include "Vector2.hh"

namespace Tron {

  namespace ExtremumSearcher {
    using Func_t   = std::function<Double_t(Double_t)>;
    using Result_t = Vector2<Double_t>;

    class ExtremumSearcher {
      TRON_SINGLETONIZE(ExtremumSearcher);
    public:
      virtual Result_t FindMinimum(const Func_t& func, Double_t from, Double_t to, Double_t dx = 1.0e-3, Double_t dy = 1.0e-3) const;
      virtual Result_t FindMaximum(const Func_t& func, Double_t from, Double_t to, Double_t dx = 1.0e-3, Double_t dy = 1.0e-3) const;
    };

    class LinearSearchMethod : public ExtremumSearcher {
      TRON_SINGLETONIZE(LinearSearchMethod);
    public:
      virtual Result_t FindMinimum(const Func_t& func, Double_t from, Double_t to, Double_t dx = 1.0e-3, Double_t dy = 1.0e-3) const override;
      virtual Result_t FindMaximum(const Func_t& func, Double_t from, Double_t to, Double_t dx = 1.0e-3, Double_t dy = 1.0e-3) const override;
    };
  
    class GoldenSectionSearchMethod : public ExtremumSearcher {
      TRON_SINGLETONIZE(GoldenSectionSearchMethod);
    public:
      virtual Result_t FindMinimum(const Func_t& func, Double_t from, Double_t to, Double_t dx = 1.0e-3, Double_t dy = 1.0e-3) const override;
      virtual Result_t FindMaximum(const Func_t& func, Double_t from, Double_t to, Double_t dx = 1.0e-3, Double_t dy = 1.0e-3) const override;
    private:
      static Double_t fRoughSearchFactor;
    public:
      inline static void   SetRoughSearchFactor(Double_t factor) {
        fRoughSearchFactor = factor;
      }
      inline static Double_t GetRoughSearchFactor() {
        return fRoughSearchFactor;
      }
    };

  }
    
}

#endif

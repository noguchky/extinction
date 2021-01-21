#ifndef Tron_RootFinder_hh
#define Tron_RootFinder_hh

#include <functional>
#include "Types.hh"
#include "Singleton.hh"

namespace Tron {

  namespace RootFinder {

    using Func_t = std::function<Double_t(Double_t)>;

    class RootFinder {
      TRON_SINGLETONIZE(RootFinder);
    public:
      virtual Double_t Find(const Func_t& func, Double_t from, Double_t to, Double_t dx = 1e-3, Double_t dy = 1e-3) const;
    };

    class LinearSearchMethod : public RootFinder {
      TRON_SINGLETONIZE(LinearSearchMethod);
    public:
      virtual Double_t Find(const Func_t& func, Double_t from, Double_t to, Double_t dx = 1e-3, Double_t dy = 1e-3) const override;
    };
  
    class BisectionMethod : public RootFinder {
      TRON_SINGLETONIZE(BisectionMethod);
    public:
      virtual Double_t Find(const Func_t& func, Double_t from, Double_t to, Double_t dx = 1e-3, Double_t dy = 1e-3) const override;
    private:
      static Double_t fRoughSearchFactor;
    public:
      static void     SetRoughSearchFactor(Double_t factor) {
        fRoughSearchFactor = factor;
      }
      static Double_t GetRoughSearchFactor() {
        return fRoughSearchFactor;
      }
    };
  
    // class SecantMethod : public RootFinder {
    //   TRON_SINGLETONIZE(SecantMethod);
    // public:
    //   virtual Double_t Find(const Func_t& func, Double_t from, Double_t to, Double_t dx = 1e-3, Double_t dy = 1e-3) const override;
    // };
  
    // class NewtonMethod : public RootFinder {
    //   TRON_SINGLETONIZE(NewtonMethod);
    // public:
    //   virtual Double_t Find(const Func_t& func, Double_t from, Double_t to, Double_t dx = 1e-3, Double_t dy = 1e-3) const override;
    // };

  }
  
}

#endif

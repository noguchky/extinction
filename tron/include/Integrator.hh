#ifndef Tron_Integrator_hh
#define Tron_Integrator_hh

#include <functional>
#include "Types.hh"
#include "Singleton.hh"

class TGraph;

namespace Tron {

  namespace Integrator {

    using Func_t = std::function<Double_t(Double_t)>;

    class Integrator {
      TRON_SINGLETONIZE(Integrator);
    public:
      virtual Double_t Integrate(const Func_t& func, Double_t from, Double_t to, Double_t dx) const;
      virtual Double_t Integrate(const Func_t& func, Double_t& x, Double_t from, Double_t to, Double_t dx) const;
    };

    class TrapezoidalMethod : public Integrator {
      TRON_SINGLETONIZE(TrapezoidalMethod);
    public:
      virtual Double_t Integrate(const Func_t& func, Double_t from, Double_t to, Double_t dx) const override;
      virtual Double_t Integrate(const Func_t& func, Double_t& x, Double_t from, Double_t to, Double_t dx) const override;
    };

    class SimpsonMethod : public Integrator {
      TRON_SINGLETONIZE(SimpsonMethod);
    public:
      virtual Double_t Integrate(const Func_t& func, Double_t from, Double_t to, Double_t dx) const override;
      virtual Double_t Integrate(const Func_t& func, Double_t& x, Double_t from, Double_t to, Double_t dx) const override;
    };

    class PowerSimpsonMethod : public Integrator {
      TRON_SINGLETONIZE(PowerSimpsonMethod);
    public:
      virtual Double_t Integrate(const Func_t& func, Double_t from, Double_t to, Double_t nx) const override;
      virtual Double_t Integrate(const Func_t& func, Double_t& x, Double_t from, Double_t to, Double_t nx) const override;
    };

    class LineGraphMethod {
      TRON_SINGLETONIZE(LineGraphMethod);
    public:
      // Graph must be sorted
      virtual Double_t Integrate(const TGraph* graph, Double_t from, Double_t to) const;
    };

  }
  
}

#endif
